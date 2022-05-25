/*
 ************************************************************************
 *                                                                      *
 * timPyscript -- Tool Integration Module for embedded Python           *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#if        PY_MAJOR_VERSION < 3
   #error  "requires at least python version 3.8"
#elif      PY_MINOR_VERSION < 8
   #error  "requires at least python version 3.8"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "egads.h"
#include "common.h"
#include "OpenCSM.h"
#include "tim.h"
#include "emp.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define REDIRECT_STDOUT_STDERR 1
#define SHOW_STDOUT_STDERR     1

#ifdef WIN32
   #define LONG long long
   #define SLEEP(msec)  Sleep(msec)
#else
   #define LONG long
   #include <unistd.h>
   #define SLEEP(msec)  usleep(1000*msec)
#endif

/* macros */
static void *realloc_temp=NULL;              /* used by RALLOC macro */

void executePyscript(void *esp);
static void tee(void *arg);

static int    initialized   = 0;
static int    killTee       = 0;
static int    pythonRunning = 0;

#define       MMODLS  10
static int    nmodls = 0;
static modl_T *modls[MMODLS];


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        void  *data)                    /* (in)  script name */
{
    int    status=0;                    /* (out) return status */

    int      len, i, count, nbuffer=1024;
    char     *env, *buffer=NULL, templine[MAX_LINE_LEN], *filename=(char*)data;
    void     *temp;
    FILE     *fp;

    ROUTINE(timLoad(pyscript));

    /* --------------------------------------------------------------- */

    /* remember the filename */
    len = strlen(filename);

    MALLOC(ESP->udata, char, len+1);

    strcpy(ESP->udata, filename);

    /* make sure we have a pyscript file */
    if (strcmp(&filename[len-3], ".py") != 0) {
        MALLOC(buffer, char, nbuffer+1);
        snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: \"%s\" does not end with \".py\"", filename);
        tim_bcst("pyscript", buffer);
        FREE(buffer);

        FREE(ESP->udata);
        goto cleanup;
    }

    /* open the file containing the pyscript script */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        MALLOC(buffer, char, nbuffer+1);
        snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: Could not open \"%s\"", filename);
        tim_bcst("pyscript", buffer);
        FREE(buffer);

        FREE(ESP->udata);
        goto cleanup;
    }
    fclose(fp);

    /* set Python Home if not set */
    if (getenv("PYTHONHOME") == NULL) {
        env = getenv("PYTHONINC");
        if (env == NULL) {
            MALLOC(buffer, char, nbuffer+1);
            snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: neither PYTHONHOME nor PYTHONINC are set");
            tim_bcst("pyscript", buffer);
            FREE(buffer);

            FREE(ESP->udata);
            goto cleanup;
        }

        len   = strlen(env);
        count = 0;
        for (i = len-1; i > 0; i--) {
#ifdef WIN32
            if (env[i] == '\\') break;
#else
            if (env[i] == '/') {
                count++;
                if (count == 2) break;
            }
#endif
        }

        if (i <= 0) {
            MALLOC(buffer, char, nbuffer+1);
            snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: PYTHONINC (%s) does not contain a path", env);
            tim_bcst("pyscript", buffer);
            FREE(buffer);

            FREE(ESP->udata);
            goto cleanup;
        }
    }

    /* send the pyscript file over to the browser */
    MALLOC(buffer, char, nbuffer);
    snprintf(buffer, nbuffer, "timLoad|pyscript|%s|", filename);

    fp = fopen(filename, "r");
    while (!feof(fp)) {
        temp = fgets(templine, MAX_LINE_LEN, fp);
        if (temp == NULL) break;

        while (strlen(buffer)+strlen(templine) > nbuffer-1) {
            nbuffer += 1024;
            RALLOC(buffer, char, nbuffer);
        }

        strncat(buffer, templine, nbuffer);
    }
    fclose(fp);

    tim_bcst("pyscript", buffer);
    FREE(buffer);

    /* do not hold the UI when executing */
    status = 0;

cleanup:
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

    int      response_len, ichar;
    char     templine[1024], *filename=NULL;
    char     *response=NULL;
    static FILE  *fp=NULL;;
    void     *temp;

    ROUTINE(timMesg(pyscript));

    /* --------------------------------------------------------------- */

    /* "fileBeg|filename|..." */
    if        (strncmp(command, "fileBeg|", 8) == 0) {

        /* extract argument */
        GetToken(command, 1, '|', &filename);

        fp = fopen(filename, "w");
        SPLINT_CHECK_FOR_NULL(fp);

        ichar = 9 + strlen(filename);
        while (command[ichar] != '\0') {
            fprintf(fp, "%c", command[ichar++]);
        }

    /* "fileMid|..." */
    } else if (strncmp(command, "fileMid|", 8) == 0) {
        SPLINT_CHECK_FOR_NULL(fp);

        ichar = 8;
        while (command[ichar] != '\0') {
            fprintf(fp, "%c", command[ichar++]);
        }

    /* "fileEnd|..." */
    } else if (strncmp(command, "fileEnd|", 8) == 0) {
        SPLINT_CHECK_FOR_NULL(fp);

        fclose(fp);
        fp = NULL;

        tim_bcst("pyscript", "timMesg|pyscript|fileEnd|");

    /* "execute" */
    } else if (strncmp(command, "execute|", 8) == 0) {
        executePyscript(ESP);

        tim_bcst("pyscript", "timMesg|pyscript|execute");

    /* "stderr" */
    } else if (strncmp(command, "stderr|", 7) == 0) {

        response_len = 1000;
        MALLOC(response, char, response_len);

        strcpy(response, "timMesg|pyscript|stderr|");

        fp = fopen("stderr.txt", "r");
        if (fp != NULL) {
            while (1) {
                temp = fgets(templine, 1023, fp);
                if (temp == NULL) {
                    fclose(fp);
                    fp = NULL;
                    remove("stderr.txt");
                    break;
                }

                if (strlen(response)+strlen(templine) > response_len-2) {
                    response_len += 4096;
                    RALLOC(response, char, response_len);
                }

                strcat(response, templine);
            }
        }

        tim_bcst("pyscript", response);
    }

cleanup:
    FREE(filename);
    FREE(response);

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

    ROUTINE(timSave(pyscript));

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
timQuit(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        int   unload)                   /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timQuit(pyscript));

    /* --------------------------------------------------------------- */

    /* free up the filename */
    FREE(ESP->udata);

    /* if last call, finalize pyscript */
    if (unload == 1) {
        /* if we are holding, release the hold now */
        (void) tim_lift("pyscript");

        /* if python is still running, we cannot call Py_Finalize since
           it will cause a seg fault.  so, we will have to live with
           the resultig memory leak */
        if (pythonRunning == 0) {
            if (Py_FinalizeEx() < 0) {
                status = -3;
                goto cleanup;
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
timGetModl(void  **myModl,              /* (out) pointer to active MODL */
           esp_T *ESP)                  /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timGetModl);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        printf("WARNING:: not running via serveESP\n");
    } else {
        *myModl = ESP->MODL;
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSetModl - set the active MODL                                  */
/*                                                                     */
/***********************************************************************/

int
timSetModl(void  *myModl,               /* (in)  pointer to active MODL */
           esp_T *ESP)                  /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    imodl, jmodl;

    ROUTINE(timSetModl);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        printf("WARNING:: not running via serveESP\n");
    } else if (ESP->MODL != myModl) {

        /* see if already in list */
        jmodl = -1;
        for (imodl = 0; imodl < nmodls; imodl++) {
            if (myModl == modls[imodl]) {
                jmodl = imodl;
                break;
            }
        }

        /* if not in the list, add it now */
        if (jmodl < 0) {
            if (nmodls+1 < MMODLS) {
                modls[nmodls] = myModl;
                nmodls++;
            } else {
                printf("WARNING:: maximum modls exceeded, so this modl will not be cleaned up\n");
            }
        }

        /* update the ESP structure */
        ESP->MODL = myModl;
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timGetCaps - get the active CAPS                                  */
/*                                                                     */
/***********************************************************************/

int
timGetCaps(void  **myCaps,              /* (out) pointer to active CAPS */
           esp_T *ESP)                  /* (in)  pointer to CAPS structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timGetCaps);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        printf("WARNING:: not running via serveESP\n");
    } else {
        *myCaps = ESP->CAPS;
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSetCaps - set the active CAPS                                  */
/*                                                                     */
/***********************************************************************/

int
timSetCaps(void  *myCaps,               /* (in)  pointer to active CAPS */
           esp_T *ESP)                  /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timSetCaps);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        printf("WARNING:: not running via serveESP\n");
    } else {
        ESP->CAPS = myCaps;
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   executePyscript - execute python is a separate thread             */
/*                                                                     */
/***********************************************************************/

void executePyscript(void *esp)
{
    int status = SUCCESS;

    int         saved_stdout, saved_stderr, outLevel;
    int         i, len, count;
    char        *buffer, *env;
    FILE        *fp, *fp_stdout, *fp_stderr;
    wchar_t     *pHome=NULL;
    PyPreConfig preConfig;
    PyConfig    config;
    void        *thread;


    esp_T *ESP      = (esp_T *)esp;
    char  *filename = (char *)ESP->udata;

    ROUTINE(executePyscript);

    /* --------------------------------------------------------------- */

    /* set flag that we are running */
    pythonRunning = 1;

    /* initialize the list of MODLs that we may have to clean up */
    modls[0] = ESP->MODL;
    nmodls = 1;

    /* get the current outLevel */
    outLevel = ocsmSetOutLevel(-1);

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

    /* set Python Home if not set */
    if (getenv("PYTHONHOME") == NULL) {
        PyPreConfig_InitPythonConfig(&preConfig);
        Py_PreInitialize(&preConfig);

        env = getenv("PYTHONINC");
        if (env == NULL) {
            fprintf(stderr, "Fatal error: PYTHONINC not set\n");
            status = -3;
            goto cleanup;
        }

        len = strlen(env);
        count = 0;
#ifdef WIN32
        for (i = len-1; i > 0; i--) {
            if (env[i] == '\\') break;
        }
#else
        for (i = len-1; i > 0; i--) {
            if (env[i] == '/') {
                count++;
                if (count == 2) break;
            }
        }
#endif

        len    = i+1;
        buffer = (char *) malloc(len*sizeof(char));
        if (buffer == NULL) {
            fprintf(stderr, "Cannot allocate %d bytes\n", i);
            status = -3;
            goto cleanup;
        }

        for (i = 0; i < len-1; i++) {
            buffer[i] = env[i];
        }
        buffer[len-1] = 0;

        pHome = Py_DecodeLocale(buffer, NULL);
        free(buffer);

        Py_SetPythonHome(pHome);
    }

    /* initialize the python interpreter */
    if (initialized == 0) {
        initialized = 1;

        PyConfig_InitPythonConfig(&config);
        config.buffered_stdio       = 0;
#ifdef WIN32
        config.legacy_windows_stdio = 1;
#endif
        Py_InitializeFromConfig(&config);
    }

    /* create a thread to broadcast stdout back to ESP */
    thread = EMP_ThreadCreate(tee, "stdout.txt");
    if (thread == NULL) exit(0);

    /* run from the file */
    fp = fopen(filename, "r");

    PyRun_SimpleFile(fp, filename);

    fclose(fp);

    /* clean up all variables */
    PyRun_SimpleString("for JfD3key in dir():\n    if JfD3key[0:2] != \"__\":\n        del globals()[JfD3key]\ndel JfD3key\n");

    /* wait a second and then tell the tee thread to kill itself */
    SLEEP(1000);

    killTee = 1;

    /* if an error occurred, print out the traceback */
    if (PyErr_Occurred()) {
        PyErr_Print();
    }

    /* undo all initializations made by Py_Initialize */
    if (pHome != NULL) {
        PyMem_RawFree(pHome);
    }

    /* this cannot be done here because of a bug in numpy */
//$$$    if (Py_FinalizeEx() < 0) {
//$$$        status = -3;
//$$$        goto cleanup;
//$$$    }

    /* delete all modls except the last one saved */
    for (i = 0; i < nmodls; i++) {
        if (modls[i] != ESP->MODL && modls[i] != NULL) {
            status = ocsmFree(modls[i]);
            if (status < EGADS_SUCCESS) {
                printf("ERROR:: ocsmFree failed for i=%d\n", i);
            }
            modls[i] = NULL;
        }
    }

    for (i = nmodls-1; i >= 0; i--) {
        if (modls[i] != NULL) break;
        nmodls--;
    }

    /* put stderr & stdout back */
    if (REDIRECT_STDOUT_STDERR == 1) {
        fflush(fp_stdout);
        fflush(fp_stderr);

        dup2(saved_stdout, fileno(stdout));
        dup2(saved_stderr, fileno(stderr));

        if (SHOW_STDOUT_STDERR) {
            printf("^^^^^ start of stdout ^^^^^\n"); fflush(stdout);
#ifndef WIN32
            system("cat   stdout.txt");
#else
            system("type  stdout.txt");
#endif
            printf("vvvvv end   of stdout vvvvv\n");

            printf("^^^^^ start of stderr ^^^^^\n"); fflush(stdout);
#ifndef WIN32
            system("cat   stderr.txt");
#else
            system("type  stderr.txt");
#endif
            printf("vvvvv end   of stderr vvvvv\n");
        }
    }

cleanup:
    /* unset flag that we are running */
    pythonRunning = 0;

    if (status < SUCCESS) {
        SPRINT1(0, "ERROR:: status=%d in executePyscript", status);
    }
}


/***********************************************************************/
/*                                                                     */
/*   getline - version for Windoze                                     */
/*                                                                     */
/***********************************************************************/

#ifdef WIN32
static int
getline(char   **linep,
        size_t *linecapp,
        FILE   *stream)
{
    int  i=0, count=0;
    char *buffer;

    /* --------------------------------------------------------------- */

    buffer = (char *) realloc(*linep, 1026*sizeof(char));
    if (buffer == NULL) return -1;

    *linep    = buffer;
    *linecapp = 1024;
    while (count < 1024 && !ferror(stream)) {
        i = fread(&buffer[count], sizeof(char), 1, stream);
        if (i == 0) break;
        count++;
        if ((buffer[count-1] == 10) ||
            (buffer[count-1] == 13) ||
            (buffer[count-1] ==  0)   ) break;
    }

    buffer[count] = 0;
    return count == 0 ? -1 : count;
}
#endif


/***********************************************************************/
/*                                                                     */
/*   tee - process to broadcast stdout back to ESP                     */
/*                                                                     */
/***********************************************************************/

static void
tee(void *arg)
{
    FILE   *fp;
    LONG   fpos, strt, len=0;
    char   *name = (char *) arg;
    size_t linecap=0;
    char   *newline, *line=NULL;


    /* --------------------------------------------------------------- */

    /* infinite loop */
    while (1) {

        /* open stdout.txt */
        fp = fopen(name, "r");
        if (fp == NULL) {
            printf(" Cannot open %s for read!\n", name);
            return;
        }

        /* find the character position at the end of the file */
        strt = 0;
#ifdef WIN32
        _fseeki64(fp, strt, SEEK_END);
        fpos = _ftelli64(fp);
#else
        fseek(fp, strt, SEEK_END);
        fpos =  ftell(fp);
#endif

        /* if this is not the same as last time, process this */
        if (fpos != len) {

            /* find the length of the file */
#ifdef WIN32
            _fseeki64(fp, len, SEEK_SET);
#else
            fseek(fp, len, SEEK_SET);
#endif

            /* process each line (with an added \n) */
            while (getline(&line, &linecap, fp) > 0) {
                newline = malloc((linecap+10)*sizeof(char));
                if (newline == NULL) {
                    printf("Trouble with malloc\n");
                    return;
                }

                sprintf(newline, "%s\n", line);
                wv_broadcastText(newline);
                free(newline);
            }
        }

        /* remember how the big the file was for next time and close it */
        len = fpos;
        fclose(fp);

        /* wait a short while */
        SLEEP(250);

        /* if the main process wants this process killed, do it now */
        if (killTee == 1) {
            free(line);

            killTee = 0;
            return;
        }
    }
}
