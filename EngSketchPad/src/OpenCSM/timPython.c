/*
 ************************************************************************
 *                                                                      *
 * timPython -- Tool Integration Module for embedded Python             *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

// for blastoff
#include <unistd.h>

#include "egads.h"
#include "common.h"
#include "OpenCSM.h"
#include "tim.h"
#include "emp.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define REDIRECT_STDOUT_STDERR 1

/* macros */
static void *realloc_temp=NULL;              /* used by RALLOC macro */


/***********************************************************************/
/*                                                                     */
/* global variable holding the active MODL                             */
/*                                                                     */
/***********************************************************************/

static void      *oldMODL=NULL;        /* pointer to old active MODL */
static void      *newMODL=NULL;        /* pointer to new active MODL */

void executePython(void *esp);


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        void  *data)                    /* (in)  script name */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int      len;
    char     *filename = (char *) data;

    ROUTINE(timLoad(python));

    /* --------------------------------------------------------------- */

    /* remember the filename */
    len = strlen(filename);

    MALLOC(ESP->udata, char, len+1);

    strcpy(ESP->udata, filename);

cleanup:
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

    ROUTINE(timSave(python));

    /* --------------------------------------------------------------- */

    /* free up the filename */
    FREE(ESP->udata);

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timQuit(python));

    /* --------------------------------------------------------------- */

    /* free up the filename */
    FREE(ESP->udata);

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timMesg - get command, process, and return response               */
/*                                                                     */
/***********************************************************************/

int
timMesg(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        char  command[],                /* (in)  command */
        int   *max_resp_len,            /* (in)  length of response */
        char  *response[])              /* (out) response */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    char     templine[1024];
    FILE     *fp;
    void     *temp;

    ROUTINE(timMesg(python));

    /* --------------------------------------------------------------- */

    /* "execute" */
    if        (strncmp(command, "execute|", 8) == 0) {
        EMP_ThreadCreate(executePython, (void *)ESP);

        snprintf(*response, *max_resp_len, "execute");

    /* "stdout" */
    } else if  (strncmp(command, "stdout|", 7) == 0) {

        snprintf(*response, *max_resp_len, "stdout|");

        fp = fopen("stdout.txt", "r");
        if (fp != NULL) {
            while (1) {
                temp = fgets(templine, 1023, fp);
                if (temp == NULL) {
                    fclose(fp);
                    remove("stdout.txt");
                    break;
                }

                if (strlen(*response)+strlen(templine) > *max_resp_len-2) {
                    (*max_resp_len) += 4096;
                    RALLOC(*response, char, *max_resp_len);
                }

                strcat(*response, templine);
            }
        }

    /* "stderr" */
    } else if (strncmp(command, "stderr|", 7) == 0) {

        snprintf(*response, *max_resp_len, "stderr|");

        fp = fopen("stderr.txt", "r");
        if (fp != NULL) {
            while (1) {
                temp = fgets(templine, 1023, fp);
                if (temp == NULL) {
                    fclose(fp);
                    remove("stderr.txt");
                    break;
                }

                if (strlen(*response)+strlen(templine) > *max_resp_len-2) {
                    (*max_resp_len) += 4096;
                    RALLOC(*response, char, *max_resp_len);
                }

                strcat(*response, templine);
            }
        }

    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timGetModl - get the active MODL                                  */
/*                                                                     */
/***********************************************************************/

int
timGetModl(void **myModl)               /* (out) pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(GetModl);

    /* --------------------------------------------------------------- */

    *myModl = oldMODL;

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSetModl - set the active MODL                                  */
/*                                                                     */
/***********************************************************************/

int
timSetModl(void *myModl)                /* (in)  pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(SetModl);

    /* --------------------------------------------------------------- */

    newMODL = myModl;

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timViewModl - view the active MODL in serveESP                    */
/*                                                                     */
/***********************************************************************/

int
timViewModl(void *myModl)               /* (in)  pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int i;

    ROUTINE(ViewModl);

    /* --------------------------------------------------------------- */

    newMODL = myModl;

    printf("Counting down...\n");
    for (i = 10; i > 0; i--) {
        printf("     %d\n", i);
        usleep(1000000);
    }
    printf("     ...Blastoff!!\n");

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   executePython - execute python is a separate thread               */
/*                                                                     */
/***********************************************************************/

void executePython(void *esp)
{
    int status = SUCCESS;

    int      saved_stdout, saved_stderr, oldLevel, outLevel=1;
    FILE     *fp, *fp_stdout, *fp_stderr;
    PyConfig config;

    esp_T *ESP      = (esp_T *)esp;
    char  *filename = (char *)ESP->udata;

    ROUTINE(executePython);

    /* --------------------------------------------------------------- */

    oldMODL = ESP->MODL;
    newMODL = ESP->MODL;

    /* get teh current outLevel */
    oldLevel = ocsmSetOutLevel(1);
    (void)     ocsmSetOutLevel(oldLevel);
    outLevel = oldLevel;

    /* update the thread using the context */
    if (ESP->MODL->context != NULL) {
        status = EG_updateThread(ESP->MODL->context);
        CHECK_STATUS(EG_updateThread);
    }

    /* redirect stdout & stderr */
    if (REDIRECT_STDOUT_STDERR == 1) {
        saved_stdout = dup(fileno(stdout));
        saved_stderr = dup(fileno(stderr));

        fp_stdout = freopen("stdout.txt", "w", stdout);
        fp_stderr = freopen("stderr.txt", "w", stderr);
    }

    /* make sure we have a python file */
    if (strstr(filename, ".py") == NULL) {
        fprintf(stderr, "\"%s\" does not end with \".py\"\n", filename);
        status = -2;
        goto cleanup;
    }

    /* open the file containing the python script */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open \"%s\"\n", filename);
        status = -3;
        goto cleanup;
    }

    /* initialize the python interpreter */
    PyConfig_InitPythonConfig(&config);
    config.buffered_stdio       = 0;
#ifdef WIN32
    config.legacy_windows_stdio = 1;
#endif
    Py_InitializeFromConfig(&config);

    /* run from the file */
    PyRun_SimpleFile(fp, filename);

    fclose(fp);

    /* if an error occurred, print out the traceback */
    if (PyErr_Occurred()) {
        PyErr_Print();
    }

    /* undo all initializations made by Py_Initialize */
    if (Py_FinalizeEx() < 0) {
        status = -3;
        goto cleanup;
    }

    /* if the MODL has been changed, delete the oldMODL and
       inform serveESP about the newMODL */
    if (oldMODL != newMODL) {
        status = ocsmFree(oldMODL);
        CHECK_STATUS(ocsmFree);

        ESP->MODL = newMODL;

        oldMODL = newMODL;
        newMODL = NULL;
    }

    /* put stderr & stdout back */
    if (REDIRECT_STDOUT_STDERR == 1) {
        fflush(fp_stdout);
        fflush(fp_stderr);

        dup2(saved_stdout, fileno(stdout));
        dup2(saved_stderr, fileno(stderr));

//$$$        printf("^^^^^ start of stdout ^^^^^\n"); fflush(stdout);
//$$$        system("cat   stdout.txt");
//$$$        printf("vvvvv end   of stdout vvvvv\n");
//$$$
//$$$        printf("^^^^^ start of stderr ^^^^^\n"); fflush(stdout);
//$$$        system("cat   stderr.txt");
//$$$        printf("vvvvv end   of stderr vvvvv\n");
    }

    /* tell serveESP taht we are done */
    SPRINT1(2, "\n<<< browserToServer(text=%s)", "timMesg|python|executeDone|");

    wv_broadcastText("timMesg|python|executeDone|");

cleanup:
    if (status < SUCCESS) {
        SPRINT1(0, "ERROR:: status=%d in executePython", status);
    }
}
