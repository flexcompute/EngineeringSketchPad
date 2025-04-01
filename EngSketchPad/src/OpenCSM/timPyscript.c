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

#include "OpenCSM.h"
#include "tim.h"
#include "emp.h"
#include "caps.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define REDIRECT_STDOUT_STDERR 1
#define SHOW_STDOUT_STDERR     1

#ifdef WIN32
   #define LONG long long
   #define SLEEP(msec)     Sleep(msec)
   #define FSEEK(fp,st,op) _fseeki64(fp, st, op)
   #define FTELL(fp)       _ftelli64(fp)
   #define SLASH '\\'

   #define getcwd          _getcwd
   #define PATH_MAX        _MAX_PATH
#else
   #include <unistd.h>
   #include <limits.h>
   #define LONG long
   #define SLEEP(msec)     usleep(1000*msec)
   #define FSEEK(fp,st,op) fseek(fp, st, op)
   #define FTELL(fp)       ftell(fp)
   #define SLASH '/'
#endif

/* macros */
void executePyscript(void *esp);
static void tee(void *arg);

/* python requires a globbaly-defined main thread state */
static PyThreadState *mainThreadState = NULL;

static int    killTee       = 0;
static void   *thread       = NULL;

#define       MMODLS  10
static int    nmodls = 0;
static modl_T *modls[MMODLS];

static char   stdout_name[PATH_MAX+12];
static char   stderr_name[PATH_MAX+12];

/* globals to keep track of the lines executed */
static int    curLine = 0;                   /* last    line executed */
static int    maxLine = 0;                   /* maximum line executed */
static char   pyModule[MAX_FILENAME_LEN];    /* name of module being executed */

/* trace functions to keep track of line numbers */
#if PY_MINOR_VERSION > 8
    static int
    traceFunc(PyObject *obj, PyFrameObject *frame, int event, PyObject *arg)
    {
        const char   *module;
        PyCodeObject *code;

        if (event != PyTrace_LINE) return 0;

        code   = PyFrame_GetCode(frame);
        module = PyUnicode_AsUTF8(code->co_filename);
        if (strcmp(module, pyModule) == 0) {
            curLine = PyFrame_GetLineNumber(frame);
            if (curLine > maxLine) maxLine = curLine;
            if (REDIRECT_STDOUT_STDERR == 0) {
                printf("---------executing line %d\n", curLine);
            }
        }

        return 0;
    }
#else
    #include "frameobject.h"
    static int
    traceFunc(PyObject *obj, struct _frame *frame, int event, PyObject *arg)
    {
        const char *module;

        if (event != PyTrace_LINE) return 0;
        module = PyUnicode_AsUTF8(frame->f_code->co_filename);
        if (strcmp(module,pyModule) == 0) {
            curLine = frame->f_lineno;
            if (curLine > maxLine) maxLine = curLine;
            if (REDIRECT_STDOUT_STDERR == 0) {
                printf("---------executing line %d\n", curLine);
            }
        }

        return 0;
    }
#endif


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
    char     cwd[PATH_MAX];
    void     *temp;
    void     *realloc_temp = NULL;            /* used by RALLOC macro */
    FILE     *fp;

    ROUTINE(timLoad(pyscript));

    /* --------------------------------------------------------------- */

    if (ESP->nudata >= MAX_TIM_NESTING) {
        printf("ERROR:: cannot nest more than %d TIMs\n", MAX_TIM_NESTING);
        exit(0);
    }

    /* set up the names of the stdout and stderr files */
    (void) getcwd(cwd, PATH_MAX);

    snprintf(stdout_name, PATH_MAX+12, "%s%cstdout.txt", cwd, SLASH);
    snprintf(stderr_name, PATH_MAX+12, "%s%cstderr.txt", cwd, SLASH);
//$$$    (void) getcwd(stderr_name, PATH_MAX);     // suspiciously looks like a typo

    /* remember the filename */
    len = strlen(filename);

    ESP->nudata++;
    MALLOC(ESP->udata[ESP->nudata-1], char, len+1);

    strcpy(ESP->timName[ESP->nudata-1], "pyscript");

    strcpy(ESP->udata[ESP->nudata-1], filename);

//$$$    /* make sure we have a pyscript file */
//$$$    if (strcmp(&filename[len-3], ".py") != 0) {
//$$$        MALLOC(buffer, char, nbuffer+1);
//$$$        snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: \"%s\" does not end with \".py\"", filename);
//$$$        tim_bcst("pyscript", buffer);
//$$$        FREE(buffer);
//$$$
//$$$        FREE(ESP->udata[ESP->nudata-1]);
//$$$        ESP->timName[   ESP->nudata-1][0] = '\0';
//$$$        ESP->nudata--;
//$$$        goto cleanup;
//$$$    }

//$$$    /* try open the file containing the pyscript script */
//$$$    fp = fopen(filename, "r");
//$$$    if (fp == NULL) {
//$$$        MALLOC(buffer, char, nbuffer+1);
//$$$        snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: Could not open \"%s\"", filename);
//$$$        tim_bcst("pyscript", buffer);
//$$$        FREE(buffer);
//$$$
//$$$        FREE(ESP->udata[ESP->nudata-1]);
//$$$        ESP->timName[   ESP->nudata-1][0] = '\0';
//$$$        ESP->nudata--;
//$$$        goto cleanup;
//$$$    }
//$$$    fclose(fp);

    /* set Python Home if not set */
    if (getenv("PYTHONHOME") == NULL) {
        env = getenv("PYTHONINC");
        if (env == NULL) {
            MALLOC(buffer, char, nbuffer+1);
            snprintf(buffer, nbuffer, "timLoad|pyscript|ERROR:: neither PYTHONHOME nor PYTHONINC are set");
            tim_bcst("pyscript", buffer);
            FREE(buffer);

            FREE(ESP->udata[ESP->nudata-1]);
            ESP->timName[   ESP->nudata-1][0] = '\0';
            ESP->nudata--;
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

            FREE(ESP->udata[ESP->nudata-1]);
            ESP->timName[   ESP->nudata-1][0] = '\0';
            ESP->nudata--;
            goto cleanup;
        }
    }

    /* do not run if CaPsTeMpFiLe.py */
    if (strcmp(filename, "CaPsTeMpFiLe.py") == 0) {
        status = 0;
        goto cleanup;
    }

    /* wait a second to make sure that the browser has had
       a chance to change its mode */
    SLEEP(1000);

    /* send the pyscript file over to the browser */
    MALLOC(buffer, char, nbuffer);
    snprintf(buffer, nbuffer, "timLoad|pyscript|%s|", filename);

    fp = fopen(filename, "r");
    if (fp != NULL) {
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
    }

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
    void     *realloc_temp = NULL;            /* used by RALLOC macro */

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

    /* "execute|" */
    } else if (strncmp(command, "execute|", 8) == 0) {
        executePyscript(ESP);

        if (strcmp((char *)(ESP->udata[ESP->nudata-1]), "CaPsTeMpFiLe.py") != 0) {
        tim_bcst("pyscript", "timMesg|pyscript|execute|");
        }

    /* "stderr|" */
    } else if (strncmp(command, "stderr|", 7) == 0) {

        response_len = 1000;
        MALLOC(response, char, response_len);

        strcpy(response, "timMesg|pyscript|stderr|");

        fp = fopen(stderr_name, "r");
        if (fp != NULL) {
            while (1) {
                temp = fgets(templine, 1023, fp);
                if (temp == NULL) {
                    fclose(fp);
                    fp = NULL;
                    remove(stderr_name);
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

    /* "lineNums" */
    } else if (strncmp(command, "lineNums|", 9) == 0) {

        response_len = 1000;
        MALLOC(response, char, response_len);

        snprintf(response, response_len, "timMesg|pyscript|lineNums|%d|%d|", curLine, maxLine);

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

    int    i;

    ROUTINE(timSave(pyscript));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "pyscript") != 0) {
        printf("WARNING:: TIM on top of stack is not \"pyscript\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    }

    /* free up the filename */
    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

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
        int   unload)                   /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i;

    ROUTINE(timQuit(pyscript));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "pyscript") != 0) {
        printf("WARNING:: TIM on top of stack is not \"pyscript\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    }

    /* if we entered pyscript from somewhere other than capsMode,
       restore the original MODL */
    if (ESP->nudata == 1 || strcmp(ESP->timName[ESP->nudata-2], "capsMode") != 0) {
        ESP->MODL = ESP->MODLorig;
        ESP->CAPS = NULL;

        tim_bcst("pyscript", "returnMessage|build|0|");

        /* also clear the udp cache */
        status = ocsmFree(NULL);
    }

    /* update the display to the way it was before pyscript was executed if
       not in capsMode */
    tim_bcst("pyscript", "returnMessage|getPmtrs|");
    tim_bcst("pyscript", "returnMessage|getBrchs|");

    /* free up the filename */
    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    /* if last call, finalize pyscript */
    if (unload == 1) {
        /* if we are holding, release the hold now */
        (void) tim_lift("pyscript");
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
        ESP->MODL     = myModl;
        ESP->MODLorig = myModl;
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

//$$$    int         builtTo, nbody=0;
    int         oclass, mtype;
    ego         topRef, prev, next;
    capsObject  *myObject;
    capsProblem *myProblem;

    ROUTINE(timSetCaps);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        printf("WARNING:: not running via serveESP\n");
    } else {
        ESP->CAPS = myCaps;
        myObject  = (capsObject *) myCaps;
        myProblem = (capsProblem *) myObject->blind;
        ESP->MODL = myProblem->modl;

        status = EG_getInfo((ego)(ESP->MODL), &oclass, &mtype, &topRef, &prev, &next);
        if (status == EGADS_SUCCESS && oclass == MODEL) {
            printf("WARNING:: starting from a .egads file is not supported\n");
            goto cleanup;
        } else {
            status = EGADS_SUCCESS;
        }

        /* there is no need to build the configuration since CAPS will
           build it in a lazy manner whenever it is needed */

        /* if there are no Bodys but Branches, build and tessellate now */
//$$$        if (ESP->MODL->nbrch > 0 && ESP->MODL->nbody == 0) {
//$$$            printf("this is ocsmBuild(4)\n");
//$$$            status = ocsmBuild(ESP->MODL, 0, &builtTo, &nbody, NULL);
//$$$            CHECK_STATUS(ocsmBuild);
//$$$
//$$$            status = ocsmTessellate(ESP->MODL, 0);
//$$$            CHECK_STATUS(ocsmTessellate);
//$$$        }
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timBegPython - initialize python (in the main thread)             */
/*                                                                     */
/***********************************************************************/

int
timBegPython()
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int         i, len, start, count;
    char        *env, *buffer;
    wchar_t     *pHome=NULL;
    PyPreConfig preConfig;
    PyConfig    config;

    ROUTINE(timBegPython);

    /* --------------------------------------------------------------- */

    /* initialize the preconfiguration */
    PyPreConfig_InitPythonConfig(&preConfig);
    Py_PreInitialize(&preConfig);

    /* initialize the python interpreter */
    PyConfig_InitPythonConfig(&config);
    config.buffered_stdio          = 0;
    config.install_signal_handlers = 0;
#ifdef WIN32
    config.legacy_windows_stdio    = 1;
#endif

    /* set Python Home if not set */
    if (getenv("PYTHONHOME") == NULL) {
        env = getenv("PYTHONINC");
        if (env == NULL) {
            fprintf(stderr, "Fatal error: PYTHONINC not set\n");
            status = -3;
            goto cleanup;
        }

        len = strlen(env);
        count = 0;
        start = 0;
#ifdef WIN32
        if (env[0] == '"') start = 1;
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

        for (i = start; i < len-1; i++) {
            buffer[i-start] = env[i];
        }
        buffer[i-start] = 0;

        fprintf(stdout, "config.home being set to \"%s\"\n", buffer);

        pHome = Py_DecodeLocale(buffer, NULL);
        free(buffer);

        config.home = pHome;
//$$$        Py_SetPythonHome(pHome);
    }

    Py_InitializeFromConfig(&config);

    /* import numpy (which MUST be done from the main thread) */
    PyRun_SimpleString("try:\n import numpy\n from pyEGADS import egads\n from pyOCSM import ocsm\n from pyOCSM import esp\nexcept ImportError:\n pass\n");

    /* store the main thread state and release the GIL */
    mainThreadState = PyEval_SaveThread();

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timEndPython - finalize python (in the main thread)               */
/*                                                                     */
/***********************************************************************/

int
timEndPython()
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    ROUTINE(timEndPython);

    /* --------------------------------------------------------------- */

    /* restore the main thread */
    PyEval_RestoreThread(mainThreadState);

    /* stop python */
    if (Py_FinalizeEx() < 0) {
        status = -333;
    }
    printf("--> Py_FinalizeEx -> status=%d\n", status);

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   executePyscript - execute python in a separate thread             */
/*                                                                     */
/***********************************************************************/

void executePyscript(void *esp)
{
    int status = SUCCESS;

    int         saved_stdout, saved_stderr, i;
    char        templine[MAX_LINE_LEN];
    char        cmd[PATH_MAX+20];
    FILE        *fp, *fp2, *fp_stdout, *fp_stderr;
    wchar_t     *pHome=NULL;
    void        *temp;

    esp_T *ESP      = (esp_T *)esp;
    char  *filename = NULL;

    typedef struct {
        char    *projName;                  // name of project
        char    *curPhase;                  // name of current  Phase
        char    *parPhase;                  // name of parent's Phase
        int     branch;                     // current Branch
        int     revision;                   // current Revision
        void    *projObj;                   // CAPS Project object
    } capsMode_T;

    PyThreadState *myThreadState = NULL;

    ROUTINE(executePyscript);

    /* --------------------------------------------------------------- */

    /* create a new thread state based on the main interpreter */
    myThreadState = PyThreadState_New(mainThreadState->interp);

    /* set the thread state for this thread and acquire the GIL */
    PyEval_RestoreThread(myThreadState);

    /* get the name of the user's pyscript */
    for (i = 0; i < ESP->nudata; i++) {
        if (strcmp(ESP->timName[i], "pyscript") == 0) {
            filename = (char *)(ESP->udata[i]);
            break;
        }
    }
    SPLINT_CHECK_FOR_NULL(filename);

    /* if in CAPS mode, add contents of filename to projName/curPhase/pyscript.py */
    capsMode_T *capsMode = NULL;
    for (i = 0; i < ESP->nudata; i++) {
        if (strcmp(ESP->timName[i], "capsMode") == 0) {
            capsMode = (capsMode_T *)(ESP->udata[i]);
            break;
        }
    }

    if (capsMode != NULL) {
        snprintf(templine, MAX_LINE_LEN, "%s%c%s%cpyscript.py", capsMode->projName, SLASH, capsMode->curPhase, SLASH);

        /* do not write if source and destination are the same (happens during
           execution of pyscript during continuation mode) */
        if (strcmp(templine, filename) != 0) {
            fp  = fopen(filename, "r");
            fp2 = fopen(templine, "a");

            /* if templine is not at the beginning, add the lines that tell
               it to clear memory before appending the current script */
            if (FTELL(fp2) != 0) {
                fprintf(fp2, "\n# --------------------------------------------------\n");
                fprintf(fp2,   "for JfD3key in dir():\n");
                fprintf(fp2,   "    if JfD3key[0:2] != \"__\":\n");
                fprintf(fp2,   "        del globals()[JfD3key]\n");
                fprintf(fp2,   "del JfD3key\n");
                fprintf(fp2,   "import gc; gc.collect()\n");
                fprintf(fp2,   "# --------------------------------------------------\n\n");
            }

            /* write banner for current script */
            fprintf(fp2, "\n# ==================================================\n");
            fprintf(fp2,   "# executing file \"%s\"\n", filename);
            fprintf(fp2,   "# ==================================================\n\n");

            if (strcmp(filename, "CaPsTeMpFiLe.py") != 0) {
                /* put big try/except around script */
                fprintf(fp2, "try:\n");

                /* append the current script */
                while (!feof(fp)) {
                    temp  = fgets(templine, MAX_LINE_LEN, fp);
                    if (temp == NULL) break;

                    fprintf(fp2, "    %s", templine);
                }

                fprintf(fp2, "\nexcept:\n");
                fprintf(fp2, "    pass\n");
            } else {
                /* append the current script */
                while (!feof(fp)) {
                    temp  = fgets(templine, MAX_LINE_LEN, fp);
                    if (temp == NULL) break;

                    fprintf(fp2, "%s", templine);
                }
            }

            fclose(fp );
            fclose(fp2);
        }
    }

    /* initialize the list of MODLs that we may have to clean up */
    modls[0] = ESP->MODL;
    nmodls = 1;

    /* update the thread using the context */
    if (ESP->MODL->context != NULL) {
        status = EG_updateThread(ESP->MODL->context);
        CHECK_STATUS(EG_updateThread);
    }

    /* redirect stdout & stderr */
    if (REDIRECT_STDOUT_STDERR == 1) {
        saved_stdout = dup(fileno(stdout));
        saved_stderr = dup(fileno(stderr));

        fp_stdout = freopen(stdout_name, "w", stdout);
        fp_stderr = freopen(stderr_name, "w", stderr);
    }

    /* set an auxiliary pointer so that caps_open will know that we are
       in timPyscript and so it will use ESP->CAPS instead of the given
       project/filename */
    (void) ocsmSetAuxPtr(ESP->CAPS);

    /* create a thread to broadcast stdout back to ESP */
    thread = EMP_ThreadCreate(tee, stdout_name);
    if (thread == NULL) exit(0);

//$$$    /* import esp into python */
//$$$    PyRun_SimpleString("from pyOCSM import esp\n");

    /* set up to trace lines in the file as they are executed */
    strncpy(pyModule, filename, MAX_FILENAME_LEN);
    curLine = 0;
    maxLine = 0;

    PyEval_SetTrace(traceFunc, NULL);

    /* inform python that we are running from pyscript */
    PyRun_SimpleString("from pyOCSM import ocsm\nfrom pyOCSM import esp\nocsm.PyScRiPt = 1");

    /* run from the file */
    fp = fopen(filename, "r");
    PyRun_SimpleFile(fp, filename);
    fclose(fp);

    /* update the display */
    PyRun_SimpleString("esp.UpdateESP()\n");

    /* clean up all variables */
    PyRun_SimpleString("for JfD3key in dir():\n    if JfD3key[0:2] != \"__\":\n        del globals()[JfD3key]\ndel JfD3key\nimport gc\ngc.collect()\n");

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

    /* delete all modls except the last one saved */
    for (i = 0; i < nmodls; i++) {
        if (modls[i] != ESP->MODL && modls[i] != ESP->MODLorig && modls[i] != NULL) {
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
            snprintf(cmd, PATH_MAX+20, "cat  %s", stdout_name);
            system(cmd);
#else
            snprintf(cmd, PATH_MAX+20, "type %s", stdout_name);
            system(cmd);
#endif
            printf("vvvvv end   of stdout vvvvv\n");

            printf("^^^^^ start of stderr ^^^^^\n"); fflush(stdout);
#ifndef WIN32
            snprintf(cmd, PATH_MAX+20, "cat  %s", stderr_name);
            system(cmd);
#else
            snprintf(cmd, PATH_MAX+20, "type %s", stderr_name);
            system(cmd);
#endif
            printf("vvvvv end   of stderr vvvvv\n");
        }
    }

    /* cleanup and delete the thread state structure */
    PyThreadState_Clear(myThreadState);
    PyThreadState_DeleteCurrent();

    /* if this was not run from CAPSmode, restore from the original MODL pointer */
    if (ESP->nudata == 1) {
        ESP->MODL = ESP->MODLorig;
    }

cleanup:
    if (status < SUCCESS) {
        printf("ERROR:: status=%d in executePyscript\n", status);
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
        FSEEK(fp, strt, SEEK_END);
        fpos =  FTELL(fp);

        /* if this is not the same as last time, process this */
        if (fpos != len) {

            /* find the length of the file */
            FSEEK(fp, len, SEEK_SET);

            /* process each line (with an added \n) */
            while (getline(&line, &linecap, fp) > 0) {
                newline = malloc((linecap+10)*sizeof(char));
                if (newline == NULL) {
                    printf("Trouble with malloc\n");
                    return;
                }

                snprintf(newline, linecap+9, "%s\n", line);
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

            if (thread != NULL) {
                EMP_ThreadDestroy(thread);
                thread = NULL;
            }

            killTee = 0;
            return;
        }
    }
}
