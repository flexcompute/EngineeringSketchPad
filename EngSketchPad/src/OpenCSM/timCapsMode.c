/*
 ************************************************************************
 *                                                                      *
 * timCapsMode -- Tool Integration Module for CAPS within ESP           *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2024  John F. Dannenhoffer, III (Syracuse University)
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
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "OpenCSM.h"
#include "caps.h"
#include "tim.h"
#include "emp.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#ifdef WIN32
   #define LONG long long
   #define SLEEP(msec)  Sleep(msec)
   #define SLASH '\\'
#else
   #define LONG long
   #include <unistd.h>
   #define SLEEP(msec)  usleep(1000*msec)
   #define SLASH '/'
#endif

/* global variables */
static  int          outLevel=1;
#define MAX_BRCH_REV 100

typedef struct {
    char    *projName;                  // name of project
    char    *curPhase;                  // name of current  Phase
    char    *parPhase;                  // name of parent's Phase
    int     branch;                     // current Branch
    int     revision;                   // current Revision
    capsObj projObj;                    // CAPS Project object
} capsMode_T;

/* declarations for routines defined below */
static void addToResponse(char text[], char response[], int *max_resp_len, int *response_len);
static int  getToken(char text[], int nskip, char sep, char token[]);
static int  makeCsmForCaps(capsMode_T *capsMode, char filename[]);

/* CAPS routines not prototyped in caps.h */
extern int  caps_statFile(const char *path);   /* EGADS_OUTSIDE(1)=directory, EGADS_SUCCESS(0)=file, EGADS_NOTFOUND(-1)=neither */
extern int  caps_mkDir(const char *path);
extern int  caps_cpDir(const char *src, const char *dst);
extern int  caps_rmDir(const char *path);
extern int  caps_cpFile(const char *src, const char *dst);
extern int  caps_rmFile(const char *path);
extern int  caps_modifiedDesPmtrs(capsObj pobj, int *nignores, int *ignores[]);


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        void  *data)                    /* (in)  project[:branch[.revision]][*] */
{
    int    status=0;                    /* (out) return status */

    int            i, j, clearModl, nerror, npmtr, nbrch, nbody, ans;
    int            buildTo, builtTo, nchange, *changes, nhist, nline;
    short          datetime[6];
    char           *tempFile=NULL, *tempFile2=NULL;
    char           *token1=NULL, *token2=NULL, *buffer=NULL;
    char           *myPhase, *myProcess, *myProcID, *myUserID, **lines;
    FILE           *parent_fp, *files_fp;
    modl_T         *MODL;
    CAPSLONG       mySeqNum;
    capsErrs       *errors;
    capsOwn        *hists;
    capsObj        intentObj;

    capsMode_T     *capsMode;

    ROUTINE(timLoad(capsMode));

    /* --------------------------------------------------------------- */

    SPLINT_CHECK_FOR_NULL(data);

#ifdef DEBUG
    printf("timLoad(data=%s)\n", (char *)data);
#endif

    outLevel = ocsmSetOutLevel(-1);

#define MAX_BUFFER_LEN  10*MAX_FILENAME_LEN

    MALLOC(tempFile,  char, MAX_FILENAME_LEN);
    MALLOC(tempFile2, char, MAX_FILENAME_LEN);
    MALLOC(buffer,    char, MAX_BUFFER_LEN);

    GetToken(data, 0, '#', &token1);
    GetToken(data, 1, '#', &token2);

    SPLINT_CHECK_FOR_NULL(token1);
    SPLINT_CHECK_FOR_NULL(token2);

    if (ESP == NULL) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: cannot run without serveESP\n");
        goto cleanup;
    }

    status = ocsmInfo(ESP->MODL, &nbrch, &npmtr, &nbody);
    if (status != SUCCESS) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: ocsmInfo -> status=%d", status);
        goto cleanup;
    }

    /* create the capsMode_T structure */
    if (ESP->nudata >= MAX_TIM_NESTING) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: cannot nest more than %d TIMs\n", MAX_TIM_NESTING);
        goto cleanup;
    }

    ESP->nudata++;
    MALLOC(ESP->udata[ESP->nudata-1], capsMode_T, 1);

    strcpy(ESP->timName[ESP->nudata-1], "capsMode");

    capsMode = (capsMode_T *) (ESP->udata[ESP->nudata-1]);

    capsMode->projName = NULL;
    capsMode->curPhase = NULL;
    capsMode->parPhase = NULL;

    MALLOC(capsMode->projName, char, MAX_FILENAME_LEN);
    MALLOC(capsMode->curPhase, char, MAX_FILENAME_LEN);
    MALLOC(capsMode->parPhase, char, MAX_FILENAME_LEN);

    /* initialize it */
    capsMode->projName[0] = '\0';
    capsMode->curPhase[0] = '\0';
    capsMode->parPhase[0] = '\0';
    capsMode->branch      = 0;
    capsMode->revision    = 0;
    capsMode->projObj     = NULL;

    /* parse the user-specified Project:Branch.Revision */
    clearModl = 0;

    /* pull out Project name */
    i = 0;
    j = 0;
    while (i < strlen(token1)) {
        if ((token1[i] >= 'A' && token1[i] <= 'Z') ||
            (token1[i] >= 'a' && token1[i] <= 'z') ||
            (token1[i] >= '0' && token1[i] <= '9')   ) {
            capsMode->projName[j  ] = token1[i];
            capsMode->projName[j+1] = '\0';
            i++;
            j++;
        } else if (token1[i] == ':') {
            i++;
            break;
        } else if (token1[i] == '*') {
            clearModl = 1;
            i++;
            break;
        } else {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: bad character (%c) in \"%s\" while extracting Project name",
                     token1[i], token1);
            goto cleanup;
        }
    }
    if (strlen(capsMode->projName) == 0) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: Project name cannot be blank");
        goto cleanup;
    }

    /* pull out Branch if it exists */
    while (i < strlen(token1)) {
        if (token1[i] >= '0' && token1[i] <= '9') {
            capsMode->branch = 10 * capsMode->branch + token1[i] - '0';
            i++;
        } else if (token1[i] == '.') {
            i++;
            break;
        } else if (token1[i] == '*') {
            clearModl = 1;
            i++;
            break;
        } else {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: bad character (%c) in \"%s\" while extracting Branch",
                     token1[i], token1);
            goto cleanup;
        }
    }

    /* pull out Revision if it exists */
    while (i < strlen(token1)) {
        if (token1[i] >= '0' && token1[i] <= '9') {
            capsMode->revision = 10 * capsMode->revision + token1[i] - '0';
            i++;
        } else if (token1[i] == '.') {
            i++;
            break;
        } else if (token1[i] == '*') {
            clearModl = 1;
            i++;
            break;
        } else {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: bad character (%c) in \"%s\" while extracting Revision",
                     token1[i], token1);
            goto cleanup;
        }
    }

    /* clear the MODL if the user appended a "*" to the Project:Branch */
    if (clearModl == 1) {
        token2[0] = '\0';

        if (ESP->MODL != NULL) {
            status = ocsmFree(ESP->MODL);
            CHECK_STATUS(ocsmFree);
        }
        status = ocsmLoad("", (void **)(&ESP->MODL));
        CHECK_STATUS(ocsmLoad);

        tim_bcst("capsMode", "timDraw|");

        status = ocsmInfo(ESP->MODL, &nbrch, &npmtr, &nbody);
        CHECK_STATUS(ocsmInfo);
    }

    /* steal the lock if required by user */
    if (strcmp(token2, "<stealLock>") == 0) {
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d%ccapsLock", capsMode->projName, SLASH, capsMode->branch, capsMode->revision, SLASH);

        (void) caps_rmFile(tempFile);
    }

    /* determine if CAPS Project exists */
    status = caps_statFile(capsMode->projName);

    /* directory does not exist, so starting a new Project */
    if (status != EGADS_OUTSIDE) {

        /* check that no starting Phase was given */
       if (capsMode->branch != 0 || capsMode->revision != 0) {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: cannot specify a starting Phase for a new Project");
            goto cleanup;

        /* check that a .csm file was given */
        } else if (nbrch == 0 && npmtr <= 1) {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: there must be a .csm file for a new Project");
            goto cleanup;

            /* set up for new project */
        } else {
            strcpy(capsMode->parPhase, "0.0");
            strcpy(capsMode->curPhase, "1.1");
        }

    /* directory for project appears to exist */
    } else {

        /* make sure there is a Phase 1.1 */
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c1.1", capsMode->projName, SLASH);
        status = caps_statFile(tempFile);
        if (status != EGADS_OUTSIDE) {        // not a directory
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: Project \"%s\" does not contain a Phase \"1.1\"",
                     capsMode->projName);
            goto cleanup;
        }

        /* a specific starting Phase was given */
        if (capsMode->branch != 0 && capsMode->revision != 0) {

            /* check if Phase exists */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, capsMode->branch, capsMode->revision);
            status = caps_statFile(tempFile);
            if (status != EGADS_OUTSIDE) {    // not a directory
                snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: Project \"%s\" does not contain Phase \"%d.%d\"",
                         capsMode->projName, capsMode->branch, capsMode->revision);
                goto cleanup;
            }

            /* check if we can create a new Phase with the same Branch */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, capsMode->branch, capsMode->revision+1);
            status = caps_statFile(tempFile);

            /* next Phase on same Branch is possible */
            if (status != EGADS_OUTSIDE) {
                snprintf(capsMode->parPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, capsMode->revision  );
                snprintf(capsMode->curPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, capsMode->revision+1);

                /* we need to start a new Branch */
            } else {
                snprintf(capsMode->parPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, capsMode->revision);
                snprintf(capsMode->curPhase, MAX_FILENAME_LEN, "%d.%d", 0,                0                 );

                for (i = capsMode->branch+1; i < MAX_BRCH_REV; i++) {
                    snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, i, 1);
                    status = caps_statFile(tempFile);
                    if (status != EGADS_OUTSIDE) {
                        snprintf(capsMode->parPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, capsMode->revision);
                        snprintf(capsMode->curPhase, MAX_FILENAME_LEN, "%d.%d", i,                1                 );
                        break;
                    }
                }

                if (strcmp(capsMode->curPhase, "0.0") == 0) {
                    snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: could not find available Branch for Project \"%s\"",
                             capsMode->projName);
                    goto cleanup;
                }
            }

        /* a specific starting Phase was not given, so create next Revision in this Branch */
        } else {
            if (capsMode->branch == 0) {
                capsMode->branch = 1;
            }

            /* check if the Branch exists */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, capsMode->branch, 1);
            status = caps_statFile(tempFile);
            if (status != EGADS_OUTSIDE) {
                snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: Project \"%s\" does not contain Phase \"%d.1\"",
                          capsMode->projName, capsMode->branch);
                goto cleanup;
            }

            /* find first available Revision in this Branch */
            for (i = 2; i < MAX_BRCH_REV; i++) {
                snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, capsMode->branch, i);
                status = caps_statFile(tempFile);
                if (status != EGADS_OUTSIDE) {
                    snprintf(capsMode->parPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, i-1);
                    snprintf(capsMode->curPhase, MAX_FILENAME_LEN, "%d.%d", capsMode->branch, i  );
                    break;
                }
            }
        }
    }

#ifdef DEBUG
    printf("capsMode->projName=%s=\n", capsMode->projName);
    printf("capsMode->curPhase=%s=\n", capsMode->curPhase);
    printf("capsMode->parPhase=%s=\n", capsMode->parPhase);
    printf("CapsMode->branch  =%d\n",  capsMode->branch  );
    printf("capsMode->revision=%d\n",  capsMode->revision);
#endif

    /* if capsMode->parPhase is 0.0, make a new Project */
    if (strcmp(capsMode->parPhase, "0.0") == 0) {
        status = caps_mkDir(capsMode->projName);
        if (status != SUCCESS) {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: could not create directory \"%s\"",
                     capsMode->projName);
            goto cleanup;
        }

        /* make the Phase directory.  then make copies of the (modified) .csm file and all used
           .udc files in the Phase directory */
        status = makeCsmForCaps(capsMode, token2);
        CHECK_STATUS(makeCsmForCaps);

        /* open a new Project after setting up capsCSMFiles */
        for (j = strlen(token2)-1; j > 0; j--) {
            if (token2[j] == SLASH) {
                j++;
                break;
            }
        }
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%c%s", capsMode->projName, SLASH, capsMode->curPhase, SLASH, &(token2[j]));
        SPRINT3(1, "\n--> enter caps_open(%s:%s) new (from \"%s\")", capsMode->projName, capsMode->curPhase, tempFile);
        status = caps_open(capsMode->projName, capsMode->curPhase, 5, NULL, outLevel, &capsMode->projObj, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_open);

        /* write the parent file */
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%cparent.txt", capsMode->projName, SLASH, capsMode->curPhase, SLASH);
        parent_fp = fopen(tempFile, "w");
        if (parent_fp != NULL) {
            fprintf(parent_fp, "0.0\n");
            fclose(parent_fp);
        }

        /* save the pointer to CAPS and the MODL */
        ESP->MODLorig = ESP->MODL;
        ESP->CAPS     = (capsObj)capsMode->projObj;
        SPLINT_CHECK_FOR_NULL(ESP->CAPS);
        MODL = ESP->MODL = ((capsProblem*) (capsMode->projObj->blind))->modl;

        /* tell the browsers about the current Project and Phase */
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|%s|%s||", capsMode->projName, capsMode->curPhase);
        tim_bcst("capsMode", buffer);

    /* check the state of the capsMode->parPhase */
    } else {
        status = caps_phaseState(capsMode->projName, capsMode->parPhase, &ans);
        if (status != SUCCESS) {
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: caps_phaseState(%s, %s) -> status=%d",
                     capsMode->projName, capsMode->parPhase, status);
            goto cleanup;
        }

        /* capsMode->parPhase is locked by someone else */
        if (ans%2 == 1) {
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsLock", capsMode->projName, SLASH, capsMode->parPhase, SLASH);
            files_fp = fopen(tempFile, "r");
            if (files_fp != NULL) {
                FREE(  token1);
                MALLOC(token1, char, MAX_EXPR_LEN);
                if (fgets(token1, MAX_EXPR_LEN, files_fp) == NULL) {
                    strcpy(token1, "<unknown>");
                }
                fclose(files_fp);
            }

            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: Phase \"%s\" for Project \"%s\" is locked by:\n%s\nDo you want to steal the lock?",
                     capsMode->parPhase, capsMode->projName, token1);
            goto cleanup;

        /* capsMode->parPhase is open, so start CAPS in continuation mode */
        } else if (ans < 2) {
            SPRINT2(0, "WARNING:: Phase \"%s\" for Project \"%s\" will use continuation mode", capsMode->parPhase, capsMode->projName);

            /* move the pyscript.py file to CaPsTeMpFiLe.py (becuse we will be writing a new pyscript.py file) */
            snprintf(tempFile,  MAX_FILENAME_LEN, "%s%c%s%cpyscript.py", capsMode->projName, SLASH, capsMode->parPhase, SLASH);
            snprintf(tempFile2, MAX_FILENAME_LEN, "CaPsTeMpFiLe.py");
            (void) caps_rmFile(tempFile2);
            (void) caps_cpFile(tempFile, tempFile2);
            (void) caps_rmFile(tempFile);

            /* open CAPS in continuation mode */
            SPRINT2(1, "\n--> enter caps_open(%s:%s) continuation", capsMode->projName, capsMode->curPhase);
            status = caps_open(capsMode->projName, capsMode->parPhase, 4, capsMode->curPhase,  outLevel, &capsMode->projObj, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_open);

            SPLINT_CHECK_FOR_NULL(capsMode->projObj);

            strcpy(capsMode->curPhase, capsMode->parPhase);

            /* get the current intent phrase (with debug mode on so that it is not journalled) */
            status = caps_debug(capsMode->projObj);
            assert(status == 1);

            status = caps_childByName(capsMode->projObj, VALUE, PARAMETER, "__intent__",
                                      &intentObj, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_childByName);

            status = caps_getHistory(intentObj, &nhist, &hists);
            CHECK_STATUS(caps_getHistory);

            status = caps_ownerInfo(capsMode->projObj, hists[0], &myPhase, &myProcess, &myProcID, &myUserID,
                                    &nline, &lines, datetime, &mySeqNum);
            CHECK_STATUS(caps_ownerInfo);

            status = caps_debug(capsMode->projObj);
            assert(status == 0);

            /* save the pointer to CAPS and the MODL */
            ESP->MODLorig = ESP->MODL;
            ESP->CAPS     = (capsObj)capsMode->projObj;
            SPLINT_CHECK_FOR_NULL(ESP->CAPS);
            MODL = ESP->MODL = ((capsProblem*) (capsMode->projObj->blind))->modl;

            /* update the thread using the context */
            status = EG_updateThread(MODL->context);
            CHECK_STATUS(EG_updateThread);

            /* tell the browsers about the current Project and Phase */
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|%s|%s|<haveIntent>|", capsMode->projName, capsMode->curPhase);
            tim_bcst("capsMode", buffer);

            snprintf(buffer, MAX_BUFFER_LEN, "postMessage|Phase \"%s:%s\" being continued", capsMode->projName, capsMode->curPhase);
            tim_bcst("capsMode", buffer);

            /* execute CaPsTeMpFiLe.py */
            status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
            CHECK_STATUS(tim_load);

            status = tim_mesg("pyscript", "execute|");
            CHECK_STATUS(tim_mesg);

            status = tim_quit("pyscript");
            CHECK_STATUS(tim_quit);

            /* remove CaPsTeMpFiLe.py */
            (void) caps_rmFile(tempFile2);

        /* start a new Phase after reloading a new .csm file */
        } else if (strlen(token2) > 0 && strcmp(token2, "undefined") != 0) {

            /* make the Phase directory.  then make copies of the (modified) .csm file and all used
               .udc files in the Phase directory */
            status = makeCsmForCaps(capsMode, token2);
            CHECK_STATUS(makeCsmForCaps);

            SPRINT4(1, "\n--> enter caps_open(%s:%s) from Phase %s (from \"%s\")", capsMode->projName, capsMode->curPhase, capsMode->parPhase, token2);
            status = caps_open(capsMode->projName, capsMode->curPhase, 5, capsMode->parPhase, outLevel, &capsMode->projObj, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_open);

            /* write the parent file */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%cparent.txt", capsMode->projName, SLASH, capsMode->curPhase, SLASH);
            parent_fp = fopen(tempFile, "w");
            if (parent_fp != NULL) {
                fprintf(parent_fp, "%s\n", capsMode->parPhase);
                fclose(parent_fp);
            }

            /* remove pyscript.py file (copied by CAPS from previous Phase) */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%cpyscript.py", capsMode->projName, SLASH, capsMode->curPhase, SLASH);
            (void) caps_rmFile(tempFile);

            /* save the pointer to CAPS and the MODL */
            ESP->MODLorig = ESP->MODL;
            ESP->CAPS     = (capsObj)capsMode->projObj;
            SPLINT_CHECK_FOR_NULL(ESP->CAPS);
            MODL = ESP->MODL = ((capsProblem*) (capsMode->projObj->blind))->modl;

            /* tell the browsers about the current Project and Phase */
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|%s|%s||", capsMode->projName, capsMode->curPhase);
            tim_bcst("capsMode", buffer);

            /* post a message about the load */
            status = caps_modifiedDesPmtrs(capsMode->projObj, &nchange, &changes);
            CHECK_STATUS(caps_modifiedDesPmtrs);

            snprintf(buffer, MAX_BUFFER_LEN, "postMessage|Phase \"%s:%s\" being loaded", capsMode->projName, capsMode->curPhase);
            for (i = 0; i < nchange; i++) {
                snprintf(tempFile, MAX_FILENAME_LEN, "\n   CAPS overrides .csm value for %s", MODL->pmtr[changes[i]].name);
                strncat(buffer, tempFile, MAX_BUFFER_LEN);
            }
            tim_bcst("capsMode", buffer);

        /* start a new Phase with the old .csm file */
        } else {

            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%ccapsCSMLoad", capsMode->projName, SLASH, capsMode->parPhase, SLASH, SLASH);
            parent_fp = fopen(tempFile, "r");
            if (parent_fp != NULL) {
                fscanf(parent_fp, "%s", tempFile2);
                fclose(parent_fp);
            }
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%c%s", capsMode->projName, SLASH, capsMode->parPhase, SLASH, SLASH, tempFile2);

            /* make the Phase directory.  then make copies of the (modified) .csm file and all used
               .udc files in the Phase directory */
            status = makeCsmForCaps(capsMode, tempFile);
            CHECK_STATUS(makeCsmForCaps);

            SPRINT3(1, "\n--> enter caps_open(%s:%s) from Phase %s", capsMode->projName, capsMode->curPhase, capsMode->parPhase);
            status = caps_open(capsMode->projName, capsMode->curPhase, 5, capsMode->parPhase,  outLevel, &capsMode->projObj, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_open);

            /* write the parent file */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%cparent.txt", capsMode->projName, SLASH, capsMode->curPhase, SLASH);
            parent_fp = fopen(tempFile, "w");
            if (parent_fp != NULL) {
                fprintf(parent_fp, "%s\n", capsMode->parPhase);
                fclose(parent_fp);
            }

            /* remove pyscript.py file (copied by CAPS from previous Phase) */
            snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%cpyscript.py", capsMode->projName, SLASH, capsMode->curPhase, SLASH);
            (void) caps_rmFile(tempFile);

            /* save the pointer to CAPS and the MODL */
            ESP->MODLorig = ESP->MODL;
            ESP->CAPS     = (capsObj)capsMode->projObj;
            SPLINT_CHECK_FOR_NULL(ESP->CAPS);
            MODL = ESP->MODL = ((capsProblem*) (capsMode->projObj->blind))->modl;

            /* tell the browsers about the current Project and Phase */
            snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|%s|%s||", capsMode->projName, capsMode->curPhase);
            tim_bcst("capsMode", buffer);

            /* post a message about the load */
            status = caps_modifiedDesPmtrs(capsMode->projObj, &nchange, &changes);
            CHECK_STATUS(caps_modifiedDesPmtrs);

            snprintf(buffer, MAX_BUFFER_LEN, "postMessage|Phase \"%s:%s\" being loaded", capsMode->projName, capsMode->curPhase);
            for (i = 0; i < nchange; i++) {
                snprintf(tempFile, MAX_FILENAME_LEN, "\n   CAPS overrides .csm value for %s", MODL->pmtr[changes[i]].name);
                strncat(buffer, tempFile, MAX_BUFFER_LEN);
            }
            tim_bcst("capsMode", buffer);
        }
    }

    SPLINT_CHECK_FOR_NULL(MODL);

    /* inform ESP of the .csm and .udc files */
    snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%cfilenames.txt",
             capsMode->projName, SLASH, capsMode->curPhase, SLASH, SLASH);

    files_fp = fopen(tempFile, "r");
    if (files_fp == NULL) {
        SPRINT1(0, "ERROR:: \"%s\" could not be opened for reading", tempFile);
        status = OCSM_FILE_NOT_FOUND;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(files_fp);

    fscanf(files_fp, "%s", buffer);
    fclose(files_fp);

    tim_bcst("capsMode", buffer);

    /* build the Bodys in the MODL */
    buildTo = 0;
    nbody   = 0;
    status = ocsmBuild(MODL, buildTo, &builtTo, &nbody, NULL);
    if (status < SUCCESS) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: could not buld  MODL (status=%s)",
                 ocsmGetText(status));
        goto cleanup;
    }

    status = ocsmTessellate(MODL, 0);
    if (status < SUCCESS) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: could not tessellate MODL (status=%s)",
                 ocsmGetText(status));
        goto cleanup;
    }

    /* output statistics about the MODL */
    status = ocsmInfo(MODL, &nbrch, &npmtr, &nbody);
    if (status < SUCCESS) {
        snprintf(buffer, MAX_BUFFER_LEN, "timLoad|capsMode|ERROR:: could not get info out of MODL (status=%s)",
                 ocsmGetText(status));
        goto cleanup;
    }

    SPRINT4(1, "--> caps_open(%s) -> nbrch=%d, npmtr=%d, nbody=%d", token2, nbrch, npmtr, nbody);

    tim_bcst("capsMode", "timDraw|");

    /* hold the UI when executing */
    status = 1;

cleanup:
    if (buffer != NULL && strncmp(buffer, "timLoad|capsMode|ERROR::", 24) == 0) {
        tim_bcst("capsMode", buffer);

        (void) tim_quit("capsMode");

        if (ESP != NULL) {
            FREE(ESP->udata[ESP->nudata-1]);
            ESP->timName[   ESP->nudata-1][0] = '\0';
            ESP->nudata--;
        }

        status = SUCCESS;
    }

    FREE(token1);
    FREE(token2);

    FREE(tempFile );
    FREE(tempFile2);
    FREE(buffer   );

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timMesg - get command, process, and return response               */
/*                                                                     */
/***********************************************************************/

int
timMesg(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        char  command[])                /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int            *tempInt, nerror, irev, count, nbrch, npmtr, nbody, ipmtr, builtTo;
    int            type1, type2, ncol1, ncol2, nrow1, nrow2, irow, icol, nchange;
    int            nline, nanal, ianal, nbound, ibound, nhist, ihist, nhist2, ihist2;
    int            i, ibrch, nrow, ncol, max_resp_len=0, response_len=0, itype;
    int            ncval, icval, nUsed, *brchUsed=NULL, *revUsed=NULL;
    int            nvset, ivset, nGpts, nDpts;
    const int      *partial;
    short          datetime[6];
    bool           *tempBool;
    double         *tempDbl, value1, value2, dot1, dot2;
    const double   *reals;
    char           *arg1=NULL, *arg2=NULL, *arg3=NULL, *arg4=NULL;
    char           *response=NULL, *filename=NULL, *newPhase=NULL, *tempName=NULL, *phaseName=NULL, *modlName=NULL;
    char           *name=NULL, *entry=NULL, *tempFile=NULL, *tempFile2=NULL;
    char           *myName, *myPhase, *myPhase2, *myProcess, *myProcID, *myUserID, **lines, *tempStr;
    char           *pEnd, begend[4], *name3, prefix, *vname, *aname;
    char           name1[MAX_NAME_LEN], name2[MAX_NAME_LEN];
    const char     *units;
    const void     *data;
    FILE           *parent_fp, *temp_fp;

    capsObj        valObj, intentObj, capsAnal, link, parent, capsBound, tempProject;
    capsObj        cValue, capsVset, tempObj;
    enum capsoType otype;
    enum capssType stype;
    enum capsvType vtype;
    CAPSLONG       mySeqNum;
    capsErrs       *errors;
    capsOwn        *hists, *hists2;
    capsOwn        last;

    capsMode_T *capsMode = NULL;

    ROUTINE(timMesg(capsMode));

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("timMesg(command=%s)\n", command);
#endif

    for (i = 0; i < ESP->nudata; i++) {
        if (strcmp(ESP->timName[i], "capsMode") == 0) {
            capsMode = (capsMode_T *)(ESP->udata[i]);
        }
    }
    SPLINT_CHECK_FOR_NULL(capsMode);

    MALLOC(arg1, char, MAX_EXPR_LEN);
    MALLOC(arg2, char, MAX_EXPR_LEN);
    MALLOC(arg3, char, MAX_EXPR_LEN);
    MALLOC(arg4, char, MAX_EXPR_LEN);

    MALLOC(name, char, MAX_EXPR_LEN);
    MALLOC(entry, char, MAX_STR_LEN);

    MALLOC(filename,  char, MAX_FILENAME_LEN);
    MALLOC(newPhase,  char, MAX_FILENAME_LEN);
    MALLOC(tempName,  char, MAX_FILENAME_LEN);
    MALLOC(phaseName, char, MAX_FILENAME_LEN);
    MALLOC(modlName,  char, MAX_FILENAME_LEN);
    MALLOC(tempFile,  char, MAX_FILENAME_LEN);
    MALLOC(tempFile2, char, MAX_FILENAME_LEN);

    max_resp_len = 4096;
    MALLOC(response, char, max_resp_len+1);
    response[0] = '\0';

    /* "unlock|" */
    if (strncmp(command, "unlock|", 7) == 0) {

        /* removethe lock */
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsLock.py", capsMode->projName, SLASH, capsMode->parPhase, SLASH);
        (void) caps_rmFile(tempFile);

        /* move the pyscript.py file to CaPsTeMpFiLe.py (becuse we will be writing a new pyscript.py file) */
        snprintf(tempFile,  MAX_FILENAME_LEN, "%s%c%s%cpyscript.py", capsMode->projName, SLASH, capsMode->parPhase, SLASH);
        snprintf(tempFile2, MAX_FILENAME_LEN, "CaPsTeMpFiLe.py");
        (void) caps_rmFile(tempFile2);
        (void) caps_cpFile(tempFile, tempFile2);
        (void) caps_rmFile(tempFile);

        /* open CAPS in continuation mode */
        SPRINT2(1, "\n--> enter caps_open(%s:%s) continuation", capsMode->projName, capsMode->curPhase);
        status = caps_open(capsMode->projName, capsMode->parPhase, 4, capsMode->curPhase,  outLevel, &capsMode->projObj, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_open);

        SPLINT_CHECK_FOR_NULL(capsMode->projObj);

        strcpy(capsMode->curPhase, capsMode->parPhase);

        /* get the current intent phrase (with debug mode on so that it is not journalled) */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        status = caps_childByName(capsMode->projObj, VALUE, PARAMETER, "__intent__",
                                  &intentObj, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_childByName);

        status = caps_getHistory(intentObj, &nhist, &hists);
        CHECK_STATUS(caps_getHistory);

        status = caps_ownerInfo(capsMode->projObj, hists[0], &myPhase, &myProcess, &myProcID, &myUserID,
                                &nline, &lines, datetime, &mySeqNum);
        CHECK_STATUS(caps_ownerInfo);

        status = caps_debug(capsMode->projObj);
        assert(status == 0);

        /* create temporary pyscript */
        temp_fp = fopen("CaPsTeMpFiLe.py", "w");
        if (temp_fp != NULL) {
            fprintf(temp_fp, "# autogenerated from: %s\n\n", command);

            fprintf(temp_fp, "import pyCAPS\n");
            fprintf(temp_fp, "from   pyOCSM import esp\n\n");

            fprintf(temp_fp, "myProblem = pyCAPS.Problem(problemName = \"foo\", capsFile = \"foo\", outLevel = 1);\n\n");

            fprintf(temp_fp, "myProblem.intentPhrase([\"%s\"])\n\n", lines[0]);

            fprintf(temp_fp, "if \"__intent__\" in myProblem.parameter:\n");
            fprintf(temp_fp, "   myProblem.parameter[\"__intent__\"].value = myProblem.parameter[\"__intent__\"].value + 1\n");
            fprintf(temp_fp, "else:\n");
            fprintf(temp_fp, "   myProblem.parameter.create(\"__intent__\", 1)\n");

            fclose(temp_fp);
        }

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile("CaPsTeMpFiLe.py");

        /* save the pointer to CAPS and the MODL */
        ESP->MODLorig = ESP->MODL;
        ESP->CAPS     = (capsObj)capsMode->projObj;
        SPLINT_CHECK_FOR_NULL(ESP->CAPS);
        ESP->MODL = ((capsProblem*) (capsMode->projObj->blind))->modl;

        /* tell the browsers about the current Project and Phase */
        snprintf(response, MAX_BUFFER_LEN, "timLoad|capsMode|%s|%s||", capsMode->projName, capsMode->curPhase);
        tim_bcst("capsMode", response);

        snprintf(response, MAX_BUFFER_LEN, "postMessage|Phase \"%s:%s\" being continued", capsMode->projName, capsMode->curPhase);
        tim_bcst("capsMode", response);

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile(tempFile2);

    /* "getCvals|" */
    } else if (strncmp(command, "getCvals|", 9) == 0) {

        /* set debug mode so that caps_ calls are not journalled */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        /* build the response in JSON format */
        status = caps_size(ESP->CAPS, VALUE, PARAMETER, &ncval, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_size);

        if (ncval > 0) {
            snprintf(response, max_resp_len, "getCvals|[");
        } else {
            snprintf(response, max_resp_len, "getCvals|");
        }
        response_len = STRLEN(response);

        for (icval = 1; icval <= ncval; icval++) {
            status = caps_childByIndex(ESP->CAPS, VALUE, PARAMETER, icval, &cValue);
            CHECK_STATUS(caps_childByIndex);

            status = caps_info(cValue, &name3, &otype, &stype, &link, &parent, &last);
            CHECK_STATUS(caps_info);

            /* do not show dunder (double-under) Cvals */
            if (strncmp(name3, "__", 2) == 0) continue;

            status = caps_getValue(cValue, &vtype, &nrow, &ncol, &data, &partial, &units, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_getValue);

            reals = data;

            snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"nrow\":%d,\"ncol\":%d,\"value\":[",
                     name3, nrow, ncol);
            addToResponse(entry, response, &max_resp_len, &response_len);

            for (i = 0; i < nrow*ncol-1; i++) {
                snprintf(entry, MAX_STR_LEN, "%lg,", reals[i]);
                addToResponse(entry, response, &max_resp_len, &response_len);
            }
            if (icval < ncval) {
                snprintf(entry, MAX_STR_LEN, "%lg]},", reals[nrow*ncol-1]);
            } else {
                snprintf(entry, MAX_STR_LEN, "%lg]}]", reals[nrow*ncol-1]);
            }
            addToResponse(entry, response, &max_resp_len, &response_len);
        }

        /* unset caps debug mode */
        status = caps_debug(capsMode->projObj);
        assert(status == 0);

    /* "newCval|name|nrow|ncol|units|value1|..." */
    } else if (strncmp(command, "newCval|", 8) == 0) {

        /* extract arguments */
        nrow = 0;
        ncol = 0;

        if (getToken(command,  1, '|', name  ) == 0) name[0] = '\0';
        if (getToken(command,  2, '|', arg1  )     ) nrow = strtol(arg1, &pEnd, 10);
        if (getToken(command,  3, '|', arg2  )     ) ncol = strtol(arg2, &pEnd, 10);
        if (getToken(command,  4, '|', arg3  ) == 0) arg3[0] = '\0';

        /* create temporary pyscript */
        temp_fp = fopen("CaPsTeMpFiLe.py", "w");
        if (temp_fp != NULL) {
            fprintf(temp_fp, "# autogenerated from: %s\n\n", command);

            fprintf(temp_fp, "import pyCAPS\n");
            fprintf(temp_fp, "from   pyOCSM import esp\n\n");

            fprintf(temp_fp, "myProblem = pyCAPS.Problem(problemName = \"foo\", capsFile = \"foo\", outLevel = 1);\n\n");

            fprintf(temp_fp, "myProblem.parameter.create(\"%s\", ", name);

            if (nrow == 1 && ncol == 1) {
                getToken(command, 5, '|', arg4);
                fprintf(temp_fp, "float(%s))\n", arg4);
            } else if (nrow == 1 || ncol == 1) {
                fprintf(temp_fp, "(");
                for (i = 0; i < nrow*ncol; i++) {
                    getToken(command, i+5, '|', arg4);
                    if (i < nrow*ncol-1) {
                        fprintf(temp_fp, "float(%s), ", arg4);
                    } else {
                        fprintf(temp_fp, "float(%s)))\n", arg4);
                    }
                }
            } else {
                fprintf(temp_fp, "((");
                i = 0;
                for (irow = 0; irow < nrow; irow++) {
                    for (icol = 0; icol < ncol; icol++) {
                        getToken(command, i+5, '|', arg4);
                        if (icol < ncol-1) {
                            fprintf(temp_fp, "float(%s), ", arg4);
                        } else {
                            fprintf(temp_fp, "float(%s))", arg4);
                        }
                        i++;
                    }
                    if (irow < nrow-1) {
                        fprintf(temp_fp, ", (");
                    } else {
                        fprintf(temp_fp, "))\n");
                    }
                }
            }

            fclose(temp_fp);
        }

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile("CaPsTeMpFiLe.py");

    /* "setCval|name|nrow|ncol|value1|... " */
    } else if (strncmp(command, "setCval|", 8) == 0) {

        /* extract arguments */
        nrow  = 0;
        ncol  = 0;

        getToken(command, 1, '|', name);
        if (getToken(command, 2, '|', arg2)) nrow  = strtol(arg2, &pEnd, 10);
        if (getToken(command, 3, '|', arg3)) ncol  = strtol(arg3, &pEnd, 10);
        getToken(command, 4, '|', arg4);

        /* create temporary pyscript */
        temp_fp = fopen("CaPsTeMpFiLe.py", "w");
        if (temp_fp != NULL) {
            fprintf(temp_fp, "# autogenerated from: %s\n\n", command);

            fprintf(temp_fp, "import pyCAPS\n");
            fprintf(temp_fp, "from   pyOCSM import esp\n\n");

            fprintf(temp_fp, "myProblem = pyCAPS.Problem(problemName = \"foo\", capsFile = \"foo\", outLevel = 1);\n\n");

            fprintf(temp_fp, "myProblem.parameter[\"%s\"].value = ", name);

            if (nrow == 1 && ncol == 1) {
                getToken(command, 5, '|', arg4);
                fprintf(temp_fp, "float(%s)\n", arg4);
            } else if (nrow == 1 || ncol == 1) {
                fprintf(temp_fp, "(");
                for (i = 0; i < nrow*ncol; i++) {
                    getToken(command, i+5, '|', arg4);
                    if (i < nrow*ncol-1) {
                        fprintf(temp_fp, "float(%s), ", arg4);
                    } else {
                        fprintf(temp_fp, "float(%s))\n", arg4);
                    }
                }
            } else {
                fprintf(temp_fp, "((");
                i = 0;
                for (irow = 0; irow < nrow; irow++) {
                    for (icol = 0; icol < ncol; icol++) {
                        getToken(command, i+5, '|', arg4);
                        if (icol < ncol-1) {
                            fprintf(temp_fp, "float(%s), ", arg4);
                        } else {
                            fprintf(temp_fp, "float(%s))", arg4);
                        }
                        i++;
                    }
                    if (irow < nrow-1) {
                        fprintf(temp_fp, ", (");
                    } else {
                        fprintf(temp_fp, ")\n");
                    }
                }
            }

            fclose(temp_fp);
        }

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile("CaPsTeMpFiLe.py");

    /* "setPmtr|name|nrow|ncol|value1|... " */
    } else if (strncmp(command, "setPmtr|", 8) == 0) {

        /* extract arguments */
        nrow  = 0;
        ncol  = 0;

        getToken(command, 1, '|', name);
        if (getToken(command, 2, '|', arg2)) nrow  = strtol(arg2, &pEnd, 10);
        if (getToken(command, 3, '|', arg3)) ncol  = strtol(arg3, &pEnd, 10);
        getToken(command, 4, '|', arg4);

        /* create temporary pyscript */
        temp_fp = fopen("CaPsTeMpFiLe.py", "w");
        if (temp_fp != NULL) {
            fprintf(temp_fp, "# autogenerated from: %s\n\n", command);

            fprintf(temp_fp, "import pyCAPS\n");
            fprintf(temp_fp, "from   pyOCSM import esp\n\n");

            fprintf(temp_fp, "myProblem = pyCAPS.Problem(problemName = \"foo\", capsFile = \"foo\", outLevel = 1);\n\n");

            fprintf(temp_fp, "myProblem.geometry.despmtr[\"%s\"].value = ", name);

            if (nrow == 1 && ncol == 1) {
                getToken(command, 5, '|', arg4);
                fprintf(temp_fp, "float(%s)\n", arg4);
            } else if (nrow == 1 || ncol == 1) {
                fprintf(temp_fp, "(");
                for (i = 0; i < nrow*ncol; i++) {
                    getToken(command, i+5, '|', arg4);
                    if (i < nrow*ncol-1) {
                        fprintf(temp_fp, "float(%s), ", arg4);
                    } else {
                        fprintf(temp_fp, "float(%s))\n", arg4);
                    }
                }
            } else {
                fprintf(temp_fp, "((");
                i = 0;
                for (irow = 0; irow < nrow; irow++) {
                    for (icol = 0; icol < ncol; icol++) {
                        getToken(command, i+5, '|', arg4);
                        if (icol < ncol-1) {
                            fprintf(temp_fp, "float(%s), ", arg4);
                        } else {
                            fprintf(temp_fp, "float(%s))", arg4);
                        }
                        i++;
                    }
                    if (irow < nrow-1) {
                        fprintf(temp_fp, ", (");
                    } else {
                        fprintf(temp_fp, ")\n");
                    }
                }
            }

            fclose(temp_fp);
        }

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile("CaPsTeMpFiLe.py");

    /* "intent|message|" */
    } else if (strncmp(command, "intent|", 7) == 0) {
        getToken(command, 1, '|', arg2);

        /* create temporary pyscript */
        temp_fp = fopen("CaPsTeMpFiLe.py", "w");
        if (temp_fp != NULL) {
            fprintf(temp_fp, "# autogenerated from: %s\n\n", command);

            fprintf(temp_fp, "import pyCAPS\n");
            fprintf(temp_fp, "from   pyOCSM import esp\n\n");

            fprintf(temp_fp, "myProblem = pyCAPS.Problem(problemName = \"foo\", capsFile = \"foo\", outLevel = 1);\n\n");

            fprintf(temp_fp, "myProblem.intentPhrase([\"%s\"])\n\n", arg2);

            fprintf(temp_fp, "if \"__intent__\" in myProblem.parameter:\n");
            fprintf(temp_fp, "   myProblem.parameter[\"__intent__\"].value = myProblem.parameter[\"__intent__\"].value + 1\n");
            fprintf(temp_fp, "else:\n");
            fprintf(temp_fp, "   myProblem.parameter.create(\"__intent__\", 1)\n");

            fclose(temp_fp);
        }

        /* execute CaPsTeMpFiLe.py */
        status = tim_load("pyscript", ESP, "CaPsTeMpFiLe.py");
        CHECK_STATUS(tim_load);

        status = tim_mesg("pyscript", "execute|");
        CHECK_STATUS(tim_mesg);

        status = tim_quit("pyscript");
        CHECK_STATUS(tim_quit);

        /* remove CaPsTeMpFiLe.py */
        (void) caps_rmFile("CaPsTeMpFiLe.py");

    /* "commit" */
    } else if (strncmp(command, "commit|", 7) == 0) {

        /* copy the DESPMTRs and CFGPMTRs from ESP->MODL back into ESP->MODLorig */
        status = ocsmInfo(ESP->MODLorig, &nbrch, &npmtr, &nbody);
        CHECK_STATUS(ocsmInfo);

        nchange = 0;
        for (ipmtr = 1; ipmtr <= npmtr; ipmtr++) {
            status = ocsmGetPmtr(ESP->MODLorig, ipmtr, &type1, &nrow1, &ncol1, name1);
            CHECK_STATUS(ocsmGetPmtr);

            if (type1 != OCSM_DESPMTR && type1 != OCSM_CFGPMTR) continue;

            status = ocsmGetPmtr(ESP->MODL,     ipmtr, &type2, &nrow2, &ncol2, name2);
            CHECK_STATUS(ocsmGetPmtr);

            if (type1 == type2 && nrow1 == nrow2 && ncol1 == ncol2 && strcmp(name1, name2) == 0) {
                for (irow = 1; irow <= nrow1; irow++) {
                    for (icol = 1; icol <= ncol1; icol++) {
                        status = ocsmGetValu(ESP->MODLorig,  ipmtr, irow, icol, &value1, &dot1);
                        CHECK_STATUS(ocsmGetValu);

                        status = ocsmGetValu(ESP->MODL,      ipmtr, irow, icol, &value2, &dot2);
                        CHECK_STATUS(ocsmGetValu);

                        if (value1 != value2) {
                            status = ocsmSetValuD(ESP->MODLorig, ipmtr, irow, icol, value2);
                            CHECK_STATUS(ocsmSetValuD);

                            nchange++;
                        }
                    }
                }
            }
        }

        if (nchange > 0) {
            nbody = 0;
            status = ocsmBuild(ESP->MODLorig, 0, &builtTo, &nbody, NULL);
            CHECK_STATUS(ocsmBuild);

            tim_bcst("capsMode", "timDraw|");
        }

        /* close the current Phase */
        status = caps_close(capsMode->projObj, 1, NULL);
        CHECK_STATUS(caps_close);

        /* we no longer have a CAPS Problem */
        ESP->CAPS = NULL;
        ESP->MODL = ESP->MODLorig;

        snprintf(response, max_resp_len, "timQuit|capsMode|%d|", nchange);
        response_len = STRLEN(response);

    /* "quit|" */
    } else if (strncmp(command, "quit|", 5) == 0) {

        /* close the current Phase (with remove flag on) */
        status = caps_close(capsMode->projObj, -1, NULL);
        CHECK_STATUS(caps_close);

        /* if this is phase 1.1, remove the Project directory too */
        if (strcmp(capsMode->curPhase, "1.1") == 0) {
            (void) caps_rmDir(capsMode->projName);
        }

        /* we no longer have a CAPS Problem */
        ESP->CAPS = NULL;
        ESP->MODL = ESP->MODLorig;

        snprintf(response, max_resp_len, "timQuit|capsMode|");
        response_len = STRLEN(response);

    /* "suspend|" */
    } else if (strncmp(command, "suspend|", 8) == 0) {

        /* copy the DESPMTRs and CFGPMTRs from ESP->MODL back into ESP->MODLorig */
        status = ocsmInfo(ESP->MODLorig, &nbrch, &npmtr, &nbody);
        CHECK_STATUS(ocsmInfo);

        nchange = 0;
        for (ipmtr = 1; ipmtr <= npmtr; ipmtr++) {
            status = ocsmGetPmtr(ESP->MODLorig, ipmtr, &type1, &nrow1, &ncol1, name1);
            CHECK_STATUS(ocsmGetPmtr);

            if (type1 != OCSM_DESPMTR && type1 != OCSM_CFGPMTR) continue;

            status = ocsmGetPmtr(ESP->MODL,     ipmtr, &type2, &nrow2, &ncol2, name2);
            CHECK_STATUS(ocsmGetPmtr);

            if (type1 == type2 && nrow1 == nrow2 && ncol1 == ncol2 && strcmp(name1, name2) == 0) {
                for (irow = 1; irow <= nrow1; irow++) {
                    for (icol = 1; icol <= ncol1; icol++) {
                        status = ocsmGetValu(ESP->MODLorig,  ipmtr, irow, icol, &value1, &dot1);
                        CHECK_STATUS(ocsmGetValu);

                        status = ocsmGetValu(ESP->MODL,      ipmtr, irow, icol, &value2, &dot2);
                        CHECK_STATUS(ocsmGetValu);

                        if (value1 != value2) {
                            status = ocsmSetValuD(ESP->MODLorig, ipmtr, irow, icol, value2);
                            CHECK_STATUS(ocsmSetValuD);

                            nchange++;
                        }
                    }
                }
            }
        }

        if (nchange > 0) {
            nbody = 0;
            status = ocsmBuild(ESP->MODLorig, 0, &builtTo, &nbody, NULL);
            CHECK_STATUS(ocsmBuild);

            tim_bcst("capsMode", "timDraw|");
        }

        /* close the current Phase (leaving it open) */
        status = caps_close(capsMode->projObj, 0, NULL);
        CHECK_STATUS(caps_close);

        /* we no longer have a CAPS Problem */
        ESP->CAPS = NULL;
        ESP->MODL = ESP->MODLorig;

        snprintf(response, max_resp_len, "timQuit|capsMode|");
        response_len = STRLEN(response);

    /* "listPhases|" */
    } else if (strncmp(command, "listPhases|", 11) == 0) {

        /* set debug mode so that caps_ calls are not journalled */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        /* make a list of all Branches/Revisions used in current Phase */
        MALLOC(brchUsed, int, MAX_BRCH_REV);
        MALLOC(revUsed,  int, MAX_BRCH_REV);

        nUsed = 0;
        ibrch = 0;
        irev  = 0;
        if (getToken(capsMode->curPhase, 0, '.', arg1)) ibrch = strtol(arg1, &pEnd, 10);
        if (getToken(capsMode->curPhase, 1, '.', arg2)) irev  = strtol(arg2, &pEnd, 10);

        while (ibrch > 0 && irev > 0) {
            if (nUsed < MAX_BRCH_REV-1) {
                brchUsed[nUsed] = ibrch;
                revUsed[ nUsed] = irev;
                nUsed++;
            }

            snprintf(tempName, MAX_FILENAME_LEN, "%s%c%d.%d%cparent.txt", capsMode->projName, SLASH, ibrch, irev, SLASH);
            parent_fp = fopen(tempName, "r");
            if (parent_fp != NULL) {
                fscanf(parent_fp, "%s", response);
                fclose(parent_fp);
            }

            if (getToken(response, 0, '.', arg1)) ibrch = strtol(arg1, &pEnd, 10);
            if (getToken(response, 1, '.', arg2)) irev  = strtol(arg2, &pEnd, 10);
        }

        /* header */
        snprintf(response, max_resp_len, "caps|listPhases|List of Phases for Project \"%s\"\n", capsMode->projName);
        response_len = STRLEN(response);

        addToResponse("  Phase    Parent   Model                Intent phrase\n", response, &max_resp_len, &response_len);
        addToResponse("  -------- -------- -------------------- -------------\n", response, &max_resp_len, &response_len);

        /* look through all possible Phases */
        for (ibrch = 1; ibrch < MAX_BRCH_REV; ibrch++) {
            for (irev = 1; irev < MAX_BRCH_REV; irev++) {
                snprintf(tempName, MAX_FILENAME_LEN, "%s%c%d.%d", capsMode->projName, SLASH, ibrch, irev);
                status = caps_statFile(tempName);
                if (status != EGADS_OUTSIDE) break;

                snprintf(tempName, MAX_FILENAME_LEN, "%s%c%d.%d%cparent.txt", capsMode->projName, SLASH, ibrch, irev, SLASH);
                parent_fp = fopen(tempName, "r");
                if (parent_fp != NULL) {
                    fscanf(parent_fp, "%s", capsMode->parPhase);
                    fclose(parent_fp);
                }

                snprintf(phaseName, MAX_FILENAME_LEN, "%d.%d", ibrch, irev);

                if (strcmp(phaseName, capsMode->curPhase) == 0) {
                    tempProject = capsMode->projObj;
                } else {
                    status = caps_open(capsMode->projName, phaseName, 7, NULL, 0, &tempProject, &nerror, &errors);
                    (void) caps_printErrors(0, nerror, errors);
                    CHECK_STATUS(caps_open);
                }

                status = caps_childByName(tempProject, VALUE, PARAMETER, "__intent__",
                                          &intentObj, &nerror, &errors);
                (void) caps_printErrors(0, nerror, errors);
                CHECK_STATUS(caps_childByName);

                status = caps_getHistory(intentObj, &nhist, &hists);
                CHECK_STATUS(caps_getHistory);

                status = caps_ownerInfo(tempProject, hists[nhist-1], &myPhase, &myProcess, &myProcID, &myUserID,
                                        &nline, &lines, datetime, &mySeqNum);
                CHECK_STATUS(caps_ownerInfo);

                /* get the .csm filename for the Phase */
                snprintf(tempName, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%ccapsCSMLoad",
                         capsMode->projName, SLASH, myPhase, SLASH, SLASH);
                parent_fp = fopen(tempName, "r");
                if (parent_fp != NULL) {
                    fscanf(parent_fp, "%s", modlName);
                    fclose(parent_fp);
                }

                /* current Phase */
                prefix = ' ';
                if (strcmp(phaseName, capsMode->curPhase) == 0) {
                    prefix = '*';

                /* other Phases */
                } else {
                    for (i = 0; i < nUsed; i++) {
                        snprintf(tempName, MAX_FILENAME_LEN, "%d.%d", brchUsed[i], revUsed[i]);
                        if (strcmp(phaseName, tempName) == 0) {
                            prefix = '-';
                            break;
                        }
                    }
                }

                if (ibrch == 1 && irev == 1) {
                    snprintf(tempName, MAX_FILENAME_LEN, "%c %-8s %-8s %-20s %s\n",
                             prefix, phaseName, "", modlName, lines[nline-1]);
                } else {
                    snprintf(tempName, MAX_FILENAME_LEN, "%c %-8s %-8s %-20s %s\n",
                             prefix, phaseName, capsMode->parPhase, modlName, lines[nline-1]);
                }
                addToResponse(tempName, response, &max_resp_len, &response_len);
            }
            if (status != EGADS_OUTSIDE && irev == 1) break;
        }

        FREE(brchUsed);
        FREE(revUsed );

        /* unset caps debug mode */
        status = caps_debug(capsMode->projObj);
        assert(status == 0);

        status = SUCCESS;

    /* "listAnalyses" */
    } else if (strncmp(command, "listAnalyses|", 13) == 0) {

        /* set debug mode so that caps_ calls are not journalled */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        /* header */
        snprintf(response, max_resp_len, "caps|listAnalyses|List of Analysis objects for Project \"%s\"\n", capsMode->projName);
        response_len = STRLEN(response);

        status = caps_size(capsMode->projObj, ANALYSIS, NONE, &nanal, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_size);

        /* loop through all Analyses */
        count = 0;
        for (ianal = 1; ianal <= nanal; ianal++) {
            status = caps_childByIndex(capsMode->projObj, ANALYSIS, NONE, ianal, &capsAnal);
            CHECK_STATUS(caps_childByIndex);

            status = caps_info(capsAnal, &myName, &otype, &stype, &link, &parent, &last);
            CHECK_STATUS(caps_info);

            snprintf(tempName, MAX_FILENAME_LEN, "  Analysis %2d: %s\n", ianal, myName);
            addToResponse(tempName, response, &max_resp_len, &response_len);

            count++;

            /* add info about the heritage of the Analysis */
            strcpy(phaseName, capsMode->curPhase);

            strcpy(tempFile2, "");
            while (strcmp(phaseName, "1.1") != 0) {
                snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%c%s",
                         capsMode->projName, SLASH, phaseName, SLASH, myName);
                status = caps_statFile(tempFile);
                if (status == EGADS_OUTSIDE) {
                    if (strlen(tempFile2) > 0) {
                        snprintf(tempName, MAX_FILENAME_LEN, "    updated in phase %s\n", tempFile2);
                        addToResponse(tempName, response, &max_resp_len, &response_len);
                    }
                    strcpy(tempFile2, phaseName);
                }

                snprintf(tempName, MAX_FILENAME_LEN, "%s%c%s%cparent.txt", capsMode->projName, SLASH, phaseName, SLASH);
                parent_fp = fopen(tempName, "r");
                if (parent_fp != NULL) {
                    fscanf(parent_fp, "%s", phaseName);
                    fclose(parent_fp);
                }
            }

            snprintf(tempName, MAX_FILENAME_LEN, "    created in phase %s\n", tempFile2);
            addToResponse(tempName, response, &max_resp_len, &response_len);
        }

        if (count == 0) {
            addToResponse("  <none>\n", response, &max_resp_len, &response_len);
        }

        /* unset caps debug mode */
        status = caps_debug(capsMode->projObj);
        assert(status == 0);

    /* "listBounds|" */
    } else if (strncmp(command, "listBounds|", 11) == 0) {

        /* set debug mode so that caps_ calls are not journalled */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        /* header */
        snprintf(response, max_resp_len, "caps|listBounds|List of Bound objects for Project \"%s\"\n", capsMode->projName);
        response_len = STRLEN(response);

        status = caps_size(capsMode->projObj, BOUND, NONE, &nbound, &nerror, &errors);
        (void) caps_printErrors(0, nerror, errors);
        CHECK_STATUS(caps_size);

        /* loop through all Bounds */
        count = 0;
        for (ibound = 1; ibound <= nbound; ibound++) {
            status = caps_childByIndex(capsMode->projObj, BOUND, NONE, ibound, &capsBound);
            CHECK_STATUS(caps_childByIndex);

            status = caps_info(capsBound, &myName, &otype, &stype, &link, &parent, &last);
            CHECK_STATUS(caps_info);

            snprintf(tempName, MAX_FILENAME_LEN, "  Bound %2d: %s\n", ibound, myName);
            addToResponse(tempName, response, &max_resp_len, &response_len);

            /* list the associated Analyses */
            status = caps_size(capsBound, VERTEXSET, CONNECTED, &nvset, &nerror, &errors);
            (void) caps_printErrors(0, nerror, errors);
            CHECK_STATUS(caps_size);

            for (ivset = 1; ivset <= nvset; ivset++) {
                status = caps_childByIndex(capsBound, VERTEXSET, CONNECTED, ivset, &capsVset);
                CHECK_STATUS(caps_childByIndex);

                status = caps_info(capsVset, &vname, &otype, &stype, &link, &parent, &last);
                CHECK_STATUS(caps_info);

                status = caps_vertexSetInfo(capsVset, &nGpts, &nDpts, &tempObj, &capsAnal);
                CHECK_STATUS(caps_vertexSetInfo);

                status = caps_info(capsAnal, &aname, &otype, &stype, &link, &parent, &last);
                CHECK_STATUS(caps_info);

                snprintf(tempName, MAX_FILENAME_LEN, "    associated with vset=%s, anal=%s\n",
                         vname, aname);
                addToResponse(tempName, response, &max_resp_len, &response_len);
            }

            count++;
        }

        if (count == 0) {
            addToResponse("  <none>\n", response, &max_resp_len, &response_len);
        }

        /* unset caps debug mode */
        status = caps_debug(capsMode->projObj);
        assert(status == 0);

    /* "listHistory|name|" */
    } else if (strncmp(command, "listHistory|", 12) == 0) {
        getToken(command, 1, '|', arg2);

        /* set debug mode so that caps_ calls are not journalled */
        status = caps_debug(capsMode->projObj);
        assert(status == 1);

        /* header */
        snprintf(response, max_resp_len, "caps|listHistory|History of \"%s\"\n", arg2);
        response_len = STRLEN(response);

        if (strlen(arg2) == 0) {
            addToResponse("  <name not given>\n", response, &max_resp_len, &response_len);
        } else {
            /* try a Cval, and if not found then try a DESPMTR or OUTPMTR */
            itype = 0;
            status = caps_childByName(    capsMode->projObj, VALUE, PARAMETER,   arg2, &valObj, &nerror, &errors);
            if (status != SUCCESS) {
                itype = 1;
                status = caps_childByName(capsMode->projObj, VALUE, GEOMETRYIN,  arg2, &valObj, &nerror, &errors);
            }
            if (status != SUCCESS) {
                itype = 2;
                status = caps_childByName(capsMode->projObj, VALUE, GEOMETRYOUT, arg2, &valObj, &nerror, &errors);
            }
            if (status != SUCCESS) {
                addToResponse("  <not found>\n", response, &max_resp_len, &response_len);
                itype = -1;
            }

            /* we found "valObj", which is the required value */
            if (itype >= 0) {
                (void) caps_printErrors(0, nerror, errors);

                /* get the history associated with __intent__ (since it contains the intent phrase
                   associated with each Phase).  this will be needed as part of ListHistory's output */
                status = caps_childByName(capsMode->projObj, VALUE, PARAMETER, "__intent__", &intentObj, &nerror, &errors);
                (void) caps_printErrors(0, nerror, errors);
                CHECK_STATUS(caps_childByName);

                status = caps_getHistory(intentObj, &nhist2, &hists2);
                CHECK_STATUS(caps_getHistory);

                /* now find the history entries associated with the user-specified value (valObj) */
                status = caps_getHistory(valObj, &nhist, &hists);
                CHECK_STATUS(caps_getHistory);

                for (ihist = 0; ihist < nhist; ihist++) {

                    /* find the phase (myPhase) associated with this history entry */
                    status = caps_ownerInfo(capsMode->projObj, hists[ihist], &myPhase, &myProcess, &myProcID, &myUserID,
                                            &nline, &lines, datetime, &mySeqNum);
                    CHECK_STATUS(caps_ownerInfo);

                    if (strcmp(lines[0], "New Phase -- reload CSM") == 0 ||
                        strcmp(lines[0], "Initial Phase")           == 0   ) {
                        continue;
                    } else {
                        strcpy(begend, "end");
                    }

                    /* get the last intent phrase that matches myPhase */
                    for (ihist2 = nhist2-1; ihist2 >= 0;  ihist2--) {
                        status =  caps_ownerInfo(capsMode->projObj, hists2[ihist2], &myPhase2, &myProcess, &myProcID, &myUserID,
                                                 &nline, &lines, datetime, &mySeqNum);
                        CHECK_STATUS(caps_ownerIfo);
                        if (strcmp(myPhase2, myPhase) == 0) break;
                    }

                    /* get the intent phrase and value associated with myPhase2 */
                    status = caps_open(capsMode->projName, myPhase, 7, NULL, 0, &tempProject, &nerror, &errors);
                    if (status == SUCCESS) {
                        if        (itype == 0) {
                            status = caps_childByName(tempProject, VALUE, PARAMETER,   arg2, &valObj, &nerror, &errors);
                        } else if (itype == 1) {
                            status = caps_childByName(tempProject, VALUE, GEOMETRYIN,  arg2, &valObj, &nerror, &errors);
                        } else if (itype == 2) {
                            status = caps_childByName(tempProject, VALUE, GEOMETRYOUT, arg2, &valObj, &nerror, &errors);
                        }
                        CHECK_STATUS(caps_childByName);
                    } else {
                        tempProject = NULL;

                        if        (itype == 0) {
                            status = caps_childByName(capsMode->projObj, VALUE, PARAMETER,   arg2, &valObj, &nerror, &errors);
                        } else if (itype == 1) {
                            status = caps_childByName(capsMode->projObj, VALUE, GEOMETRYIN,  arg2, &valObj, &nerror, &errors);
                        } else if (itype == 2) {
                            status = caps_childByName(capsMode->projObj, VALUE, GEOMETRYOUT, arg2, &valObj, &nerror, &errors);
                        }
                        CHECK_STATUS(caps_childByName);
                    }

                    status = caps_getValue(valObj, &vtype, &nrow, &ncol, &data, &partial, &units, &nerror, &errors);
                    (void) caps_printErrors(0, nerror, errors);
                    CHECK_STATUS(caps_getValue);

                    /* print the result[0] for myPhase */
                    if        (vtype == Boolean) {
                        tempBool = (bool *)data;
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=%-15d  size=[%d*%d] at %s of Phase %s: %s\n",
                                 (int)(tempBool[0]), nrow, ncol, begend, myPhase, lines[0]);
                    } else if (vtype == Integer) {
                        tempInt = (int *)data;
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=%-15d  size=[%d*%d] at %s of Phase %s: %s\n",
                                 tempInt[0], nrow, ncol, begend, myPhase, lines[0]);
                    } else if (vtype == Double) {
                        tempDbl = (double *)data;
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=%15.6f  size=[%d*%d] at %s of Phase %s: %s\n",
                                 tempDbl[0], nrow, ncol, begend, myPhase, lines[0]);
                    } else if (vtype == String) {
                        tempStr = (char *)data;
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=%-15s at %s of Phase %s: %s\n",
                                 tempStr, begend, myPhase, lines[0]);
                    } else if (vtype == Tuple) {
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=<Tuple> at %s of Phase %s: %s\n",
                                 begend, myPhase, lines[0]);
                    } else if (vtype == Pointer) {
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=<Pointer> at %s of Phase %s: %s\n",
                                 begend, myPhase, lines[0]);
                    } else if (vtype == DoubleDeriv) {
                        tempDbl = (double *)data;
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=%15.6f  size=[%d*%d] at %s of Phase %s: %s\n",
                                 tempDbl[0], nrow, ncol, begend, myPhase, lines[0]);
                    } else if (vtype == PointerMesh) {
                        snprintf(tempName, MAX_FILENAME_LEN, "  value=<PointerMesh> at %s of Phase %s: %s\n",
                                 begend, myPhase, lines[0]);
                    }
                    addToResponse(tempName, response, &max_resp_len, &response_len);

                    if (tempProject != NULL) {
                        status = caps_close(tempProject, -1, NULL);
                        CHECK_STATUS(caps_close);
                    }
                }
            }
        }

        /* unset caps debug mode */
        status = caps_debug(capsMode->projObj);
        assert(status == 0);

    /* "attachList|" */
    } else if (strncmp(command, "attachList|", 11) == 0) {

    /* "attachFile|filename|" */
    } else if (strncmp(command, "attachFile|", 11) == 0) {
        getToken(command, 1, '|', arg2);

        /* check to see if the filename has one of the supported filetypes */
        if (strcmp(&arg2[strlen(arg2)-4], ".txt" ) != 0 &&
            strcmp(&arg2[strlen(arg2)-4], ".png" ) != 0 &&
            strcmp(&arg2[strlen(arg2)-4], ".pdf" ) != 0 &&
            strcmp(&arg2[strlen(arg2)-5], ".html") != 0   ) {
            snprintf(response, MAX_BUFFER_LEN, "timMesg|capsMode|ERROR:: filetype is not .txt, .png, .pdf, or .html\n");
            goto cleanup;
        }

        /* check to see if the file exists */
        status = caps_statFile(arg2);
        if (status != EGADS_SUCCESS) {
            snprintf(response, MAX_BUFFER_LEN, "timMesg|capsMode|ERROR:: file \"%s\" does not exist", arg2);
            goto cleanup;
        }

        /* check to see if "attachments" directory exists */
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%cattachments", capsMode->projName, SLASH);
        status = caps_statFile(tempFile);

        /* if a file name "attachments" exists, we have an error */
        if (status == EGADS_SUCCESS) {
            snprintf(response, MAX_BUFFER_LEN, "timMesg|capsMode|ERROR:: cannot have a file named attachments\n");
            goto cleanup;

        /* if a file or directory does not exist, make the directory now */
        } else if (status == EGADS_NOTFOUND) {
            status = caps_mkDir(tempFile);
            CHECK_STATUS(caps_mkDir);
        }

        /* find name of attachment without the directory */
        strncpy(tempFile2, arg2, MAX_FILENAME_LEN);
        for (i = strlen(arg2); i >= 0; i--) {
            if (arg2[i] == SLASH) {
                strncpy(tempFile2, &arg2[i+1], MAX_FILENAME_LEN);
                break;
            }
        }

        /* if the file already exists in the attachments directory, delete it (since caps_cpFile cannot overwrite) */
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%cattachments%c%s", capsMode->projName, SLASH, SLASH, tempFile2);
        status = caps_statFile(tempFile);
        if (status == EGADS_SUCCESS) {
            (void) caps_rmFile(tempFile);
        } else if (status == EGADS_OUTSIDE) {
            (void) caps_rmDir(tempFile);
        }

        /* make a copy of the file into the attachments directory */
        (void) caps_cpFile(arg2, tempFile);

    /* "attachOpen|filename|" */
    } else if (strncmp(command, "attachOpen|", 11) == 0) {

    }

    /* send the response to the browsers */
    if (strlen(response) > 0) {
        tim_bcst("capsMode", response);
    }

cleanup:
    FREE(response);

    FREE(arg1);
    FREE(arg2);
    FREE(arg3);
    FREE(arg4);

    FREE(name );
    FREE(entry);

    FREE(filename );
    FREE(newPhase );
    FREE(tempName );
    FREE(phaseName);
    FREE(modlName );
    FREE(tempFile );
    FREE(tempFile2);

    FREE(brchUsed);
    FREE(revUsed );

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i;

    capsMode_T *capsMode;

    ROUTINE(timSave(capsMode));

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("timSave()\n");
#endif

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "capsMode") != 0) {
        printf("WARNING:: TIM on top of stack is not \"capsMode\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    } else {
        capsMode = (capsMode_T *)(ESP->udata[ESP->nudata-1]);
    }

    if (capsMode == NULL) {
        goto cleanup;
    }

    /* free up the capsMode_T structure */
    FREE(capsMode->projName);
    FREE(capsMode->curPhase);
    FREE(capsMode->parPhase);

    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    tim_bcst("capsMode", "timSave|capsMode|");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(esp_T *ESP,                     /* (in)  pointer to ESP structure */
/*@unused@*/int   unload)               /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i;

    capsMode_T *capsMode;

    ROUTINE(timQuit(capsMode));

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("timQuit()\n");
#endif

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "capsMode") != 0) {
        printf("WARNING:: TIM on top of stack is not \"capsMode\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    } else {
        capsMode = (capsMode_T *)(ESP->udata[ESP->nudata-1]);
    }

    if (capsMode == NULL) {
        goto cleanup;
    }

    /* free up the capsMode_T structure */
    FREE(capsMode->projName);
    FREE(capsMode->curPhase);
    FREE(capsMode->parPhase);

    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    tim_bcst("capsMode", "timQuit|capsMode|");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   addToResponse - add text to response (with buffer length protection) */
/*                                                                     */
/***********************************************************************/

static void
addToResponse(char   text[],            /* (in)  text to add */
              char   response[],        /* (both) current response */
              int    *max_resp_len,     /* (both) length of response buffer */
              int    *response_len)     /* (both) length of response */
{
    int    status=SUCCESS, text_len;

    void   *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(addToResponse);

    /* --------------------------------------------------------------- */

    text_len = STRLEN(text);

    while (*response_len+text_len > *max_resp_len-2) {
        *max_resp_len += 4096;
        SPRINT1(2, "increasing max_resp_len=%d", *max_resp_len);

        RALLOC(response, char, *max_resp_len+1);
    }

    strcpy(&(response[*response_len]), text);

    *response_len += text_len;

cleanup:
    if (status != SUCCESS) {
        SPRINT0(0, "ERROR:: max_resp_len could not be increased");
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   getToken - get a token from a string                              */
/*                                                                     */
/***********************************************************************/

static int
getToken(char   text[],                 /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   sep,                    /* (in)  separator character */
         char   token[])                /* (out) token */
{
    int    status = SUCCESS;

    int    lentok, i, j, count, iskip;
    char   *newText=NULL;

    ROUTINE(getToken);

    /* --------------------------------------------------------------- */

    token[0] = '\0';
    lentok   = 0;

    MALLOC(newText, char, strlen(text)+2);

    /* convert tabs to spaces, remove leading white space, and
       compress other white space */
    j = 0;
    for (i = 0; i < STRLEN(text); i++) {

        /* convert tabs and newlines */
        if (text[i] == '\t' || text[i] == '\n') {
            newText[j++] = ' ';
        } else {
            newText[j++] = text[i];
        }

        /* remove leading white space */
        if (j == 1 && newText[0] == ' ') {
            j--;
        }

        /* compress white space */
        if (j > 1 && newText[j-2] == ' ' && newText[j-1] == ' ') {
            j--;
        }
    }
    newText[j] = '\0';

    if (strlen(newText) == 0) goto cleanup;

    /* count the number of separators */
    count = 0;
    for (i = 0; i < STRLEN(newText); i++) {
        if (newText[i] == sep) {
            count++;
        }
    }

    if (count < nskip) {
        goto cleanup;
    } else if (count == nskip && newText[strlen(newText)-1] == sep) {
        goto cleanup;
    }

    /* skip over nskip tokens */
    i = 0;
    for (iskip = 0; iskip < nskip; iskip++) {
        while (newText[i] != sep) {
            i++;
        }
        i++;
    }

    /* if token we are looking for is empty, return 0 */
    if (newText[i] == sep) goto cleanup;

    /* extract the token we are looking for */
    while (newText[i] != sep && newText[i] != '\0') {
        token[lentok++] = newText[i++];
        token[lentok  ] = '\0';

        if (lentok >= MAX_EXPR_LEN-1) {
            SPRINT0(0, "ERROR:: token exceeds MAX_EXPR_LEN");
            break;
        }
    }

    status =  STRLEN(token);

cleanup:

    FREE(newText);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   makeCsmForCaps - make .csm and .udc files for caps mode           */
/*                                                                     */
/***********************************************************************/

static int
makeCsmForCaps(capsMode_T *capsMode,      /* (in)  pointer to CAPS structure */
               char    filename[])      /* (in)  name of .csm file */
{
    int       status = SUCCESS;

    int       i, j, hasPrefix;
    char      *tempFile=NULL, *prefix=NULL, *dirname=NULL, *tempFilelist=NULL;
    char      *tok1=NULL, *tok2=NULL, *buf1=NULL, *buf2=NULL;
    FILE      *fp_src=NULL, *fp_tgt=NULL;
    modl_T    *tempMODL;

    ROUTINE(makeCsmForCaps);

    /* --------------------------------------------------------------- */

    MALLOC(tempFile, char, MAX_FILENAME_LEN);
    MALLOC(prefix,   char, MAX_FILENAME_LEN);
    MALLOC(dirname,  char, MAX_FILENAME_LEN);
    MALLOC(tok1,     char, MAX_EXPR_LEN    );
    MALLOC(tok2,     char, MAX_EXPR_LEN    );
    MALLOC(buf1,     char, MAX_LINE_LEN    );
    MALLOC(buf2,     char, MAX_LINE_LEN    );

    /* do a dummy load of the .csm file to make sure all the parts
       are there and that the file loads correctly */
    status = ocsmLoad(filename, (void**)(&tempMODL));
    CHECK_STATUS(ocsmLoad);

    /* get a list of the files used by capsMode->csmFile */
    status =  ocsmGetFilelist(tempMODL, &tempFilelist);
    CHECK_STATUS(ocsmGetFilelist);

    SPLINT_CHECK_FOR_NULL(tempFilelist);

    /* get rid of the temporary MODL */
    status = ocsmFree(tempMODL);
    CHECK_STATUS(ocsmFree);

    /* make sure the Phase subdirectory does not exist */
    snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s",
             capsMode->projName, SLASH, capsMode->curPhase);
    status = caps_statFile(tempFile);
    if        (status == EGADS_OUTSIDE) {
        SPRINT1(0, "ERROR:: directory \"%s\" already exists", tempFile);
        goto cleanup;
    } else if (status == EGADS_SUCCESS) {
        SPRINT1(0, "ERROR:: file \"%s\" already exists", tempFile);
        goto cleanup;
    }

    /* make the Phase subdirectory */
    status = caps_mkDir(tempFile);
    CHECK_STATUS(caps_mkDir);

    /* make the capsCSMFiles subdirectory */
    snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles",
             capsMode->projName, SLASH, capsMode->curPhase, SLASH);
    status = caps_mkDir(tempFile);
    CHECK_STATUS(caps_mkDir);

    /* copy all .cpc, .csm, and .udc files, modifing the UDPARG and UDPRIM
       statements to change the primname to the form: $/primname (in same directory
       as .csm file) */
    prefix[0] = '\0';

    for (i = 0; i < 100; i++) {
        getToken(tempFilelist, i, '|', tok1);

        if (strlen(tok1) == 0) break;

        /* if this is the first (which is presumed to be the .csm or .cpc file)
           remember te part up to the last slash */
        if (i == 0) {
            strncpy(prefix, tok1, MAX_FILENAME_LEN-1);

            for (j = STRLEN(prefix)-1; j >= 0; j--) {
                if (prefix[j] == SLASH) {
                    prefix[j+1] = '\0';
                    break;
                }
            }
        }

        /* pull out the filename (beyond the prefix) */
        if (strncmp(prefix, tok1, strlen(prefix)) == 0) {
            hasPrefix = 1;
            j = strlen(prefix);

        /* pull out the filename (without the directory) */
        } else {
            hasPrefix = 0;
            for (j = strlen(tok1)-1; j > 0; j--) {
                if (tok1[j] == SLASH) {
                    j++;
                    break;
                }
            }
        }

        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%c%s",
                 capsMode->projName, SLASH, capsMode->curPhase, SLASH, SLASH, &(tok1[j]));

        /* if tempFile includes a slash, make the directory(s) if they do not exist already */
        if (hasPrefix == 1) {
            for (j = strlen(capsMode->projName)+strlen(capsMode->curPhase)+4; j < strlen(tempFile); j++) {
                if (tempFile[j] == SLASH) {
                    strcpy(dirname, tempFile);
                    dirname[j] = '\0';

                    status = caps_statFile(dirname);
                    if (status == EGADS_NOTFOUND) {             // does not exist
                        status = caps_mkDir(dirname);
                        if (status != EGADS_SUCCESS) {
                            goto cleanup;
                        }
                    } else if (status == EGADS_OUTSIDE) {       // directory
                    } else if (status == EGADS_SUCCESS) {       // file
                        SPRINT1(0, "ERROR:: \"%s\" cannot be a directory since file already exists", dirname);
                        status = OCSM_FILE_NOT_FOUND;
                        goto cleanup;
                    }
                }
            }
        }

        fp_src = fopen(tok1,     "r");
        fp_tgt = fopen(tempFile, "w");

        if        (fp_src == NULL) {
            SPRINT1(0, "ERROR:: \"%s\" could not be opened for reading", tok1);
            status = OCSM_FILE_NOT_FOUND;
            goto cleanup;
        } else if (fp_tgt == NULL) {
            SPRINT1(0, "ERROR:: \"%s\" could not be opened for writing", tempFile);
            status = OCSM_FILE_NOT_FOUND;
            goto cleanup;
        }

        SPLINT_CHECK_FOR_NULL(fp_src);
        SPLINT_CHECK_FOR_NULL(fp_tgt);

        while (fgets(buf1, MAX_LINE_LEN, fp_src) != NULL) {
            getToken(buf1, 0, ' ', tok1);
            getToken(buf1, 1, ' ', tok2);

            if (strlen(tok1) != 6) {
                strcpy(buf2, buf1);
            } else if ((strcmp(tok1, "udparg") == 0 || strcmp(tok1, "UDPARG") == 0 ||
                        strcmp(tok1, "udprim") == 0 || strcmp(tok1, "UDPRIM") == 0   ) &&
                       (strncmp(tok2, "$$/", 3) == 0                                 )   ) {
                strcpy(buf2, buf1);
            } else if ((strcmp(tok1, "udparg") == 0 || strcmp(tok1, "UDPARG") == 0 ||
                        strcmp(tok1, "udprim") == 0 || strcmp(tok1, "UDPRIM") == 0   ) &&
                       (strncmp(tok2, "$/", 2) == 0                                  )   ) {
                strcpy(buf2, buf1);
            } else if ((strcmp(tok1, "udparg") == 0 || strcmp(tok1, "UDPARG") == 0 ||
                        strcmp(tok1, "udprim") == 0 || strcmp(tok1, "UDPRIM") == 0   ) &&
                       (tok2[0] == '$' || tok2[0] == '/'                             )   ) {
                strcpy(buf2, tok1);

                for (j = strlen(tok2)-1; j >= 0; j--) {
                    if (tok2[j] == '/') {
                        j++;
                        break;
                    }
                }

                snprintf(tok1, MAX_EXPR_LEN, " $/%s",
                         &(tok2[j]));
                strcat(buf2, tok1);

                for (j = 2; j < 100; j++) {
                    getToken(buf1, j, ' ', tok1);
                    if (strlen(tok1) == 0) break;

                    strcat(buf2, " ");
                    strcat(buf2, tok1);
                }
                strcat(buf2, "    # <modified>\n");
            } else {
                strcpy(buf2, buf1);
            }

            fputs(buf2, fp_tgt);
        }

        fclose(fp_src);   fp_src = NULL;
        fclose(fp_tgt);   fp_tgt = NULL;

        /* append this filename to the filenames.txt file */
        strncpy(tok1, tempFile, MAX_EXPR_LEN);
        snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%cfilenames.txt",
                 capsMode->projName, SLASH, capsMode->curPhase, SLASH, SLASH);

        if (i == 0) {
            fp_tgt = fopen(tempFile, "w");
        } else {
            fp_tgt = fopen(tempFile, "a");
        }
        if (fp_tgt == NULL) {
            SPRINT1(0, "ERROR:: \"%s\" could not be opened for writing", tempFile);
            status = OCSM_FILE_NOT_FOUND;
            goto cleanup;
        }

        SPLINT_CHECK_FOR_NULL(fp_tgt);

        if (i == 0) {
            fprintf(fp_tgt, "getFilenames|%s|", tok1);
        } else {
            fprintf(fp_tgt, "%s|", tok1);
        }

        fclose( fp_tgt);   fp_tgt = NULL;
    }

    /* write the file that tells CAPS the name of the .csm file */
    snprintf(tempFile, MAX_FILENAME_LEN, "%s%c%s%ccapsCSMFiles%ccapsCSMLoad",
             capsMode->projName, SLASH, capsMode->curPhase, SLASH, SLASH);

    fp_tgt = fopen(tempFile, "w");
    if (fp_tgt == NULL) {
        SPRINT1(0, "ERROR \"%s\" could not be opened for writing", tempFile);
        status = OCSM_FILE_NOT_FOUND;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(fp_tgt);

    getToken(tempFilelist, 0, '|', tok1);
    for (j = strlen(tok1)-1; j > 0; j--) {
        if (tok1[j] == SLASH) {
            j++;
            break;
        }
    }

    fprintf(fp_tgt, "%s\n", &(tok1[j]));

    fclose(fp_tgt);   fp_tgt = NULL;

cleanup:
    if (fp_src != NULL) fclose(fp_src);
    if (fp_tgt != NULL) fclose(fp_tgt);

    EG_free(tempFilelist);

    FREE(tempFile);
    FREE(prefix  );
    FREE(dirname );
    FREE(tok1    );
    FREE(tok2    );
    FREE(buf1    );
    FREE(buf2    );

    return status;
}
