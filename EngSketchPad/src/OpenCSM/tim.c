/*
 ************************************************************************
 *                                                                      *
 * tim.c -- functions for the Tool Integration Modules                  *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
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
#include <assert.h>

#ifdef WIN32
   #define DLL        HINSTANCE
   #define WIN32_LEAN_AND_MEAN
   #include <windows.h>
   #include <winsock2.h>
   #define snprintf     _snprintf
   #define strcasecmp   stricmp
   #define SLEEP(msec)  Sleep(msec)
#else
   #define DLL void *
   #include <dlfcn.h>
   #include <dirent.h>
   #include <limits.h>
   #include <unistd.h>
   #define SLEEP(msec)  usleep(1000*msec)
#endif

#define  SHOW_STATES  0

#include "egads.h"
#include "OpenCSM.h"
#include "common.h"
#include "tim.h"
#include "emp.h"

#include <unistd.h>

#include "wsserver.h"

#define MAX_TIMS 32

typedef int (*timDLLfunc) (void);

static char      *timName[  MAX_TIMS];
static DLL        timDLL[   MAX_TIMS];
static int        timState [MAX_TIMS];
static void      *timEsp[   MAX_TIMS];
static void      *timThread[MAX_TIMS];
static void      *timMutex[ MAX_TIMS];
static int        timUnset[ MAX_TIMS];
static int        timHold[  MAX_TIMS];

// func pointers  name in DLL       name in OpenCSM     description
// -------------  ---------------   ------------------  -----------------------------------------------

// tim_Load[]     timLoad()         tim_load()          open a tim
typedef int (*TIM_Load) (esp_T *ESP, /*@null@*/void *data);
static        TIM_Load  tim_Load[MAX_TIMS];

// tim_Mesg[]     timMesg()         tim_mesg()          get command, process, and return response
typedef int (*TIM_Mesg) (esp_T *ESP, char command[]);
static        TIM_Mesg  tim_Mesg[MAX_TIMS];

// tim_Save[]     timSave()         tim_save()          save tim data and close tim instance
typedef int (*TIM_Save) (esp_T *ESP);
static        TIM_Save  tim_Save[MAX_TIMS];

// tim_Quit[]     timQuit()         tim_quit()          close tim instance without saving
typedef int (*TIM_Quit) (esp_T *ESP, int unload);
static        TIM_Quit  tim_Quit[MAX_TIMS];

//                                  tim_free()          frees up memory associated with TIMs

static int tim_nTim = 0;

/* structure to communicate between threads */
typedef struct {
    TIM_Mesg routine;
    char     *command;
    void     *data;
    int      *state;
} emp_T;

/* prototypes */
static /*@null@*/DLL timDLopen(const char *name);
static void          timDLclose(/*@unused@*/ /*@only@*/ DLL dll);
static timDLLfunc    timDLget(DLL dll, const char *symname, /*@null@*/ const char *name);
static int           timDLoaded(const char *name);
static int           timDYNload(const char *name);

static void          timExec(void *);

/*@-mustfreefresh@*/


/***********************************************************************/
/*                                                                     */
/*   tim_load - load a TIM                                             */
/*                                                                     */
/***********************************************************************/

int
tim_load(char       timName[],          /* (in)  name of TIM */
         esp_T      *ESP,               /* (in)  pointer to ESP structure */
/*@null@*/void      *data)              /* (in)  user-provided data */
{
    int i, status;

    /* --------------------------------------------------------------- */

//$$$    printf("tim_load(%s)\n", timName);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, try to load it now */
    if (i == -1) {
        i = timDYNload(timName);
        if (i < 0) {
            printf("ERROR:: tim_load(%s) could not be dynamically loaded\n", timName);
            status = i;
            goto cleanup;
        }

        timState[i] = TIM_INACTIVE;
        if (SHOW_STATES) printf("in tim_load: setting TIM_INACTIVE\n");
    }

    /* if not currently "inactive", return an error */
    if (timState[i] != TIM_INACTIVE) {
        printf("ERROR:: tim_load(%s) is not inactive, state=%d\n", timName, timState[i]);
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    timState[i] = TIM_LOADING;
    if (SHOW_STATES) printf("in tim_load: setting TIM_LOADING\n");

    /* save the data and call the load routine in the TIM */
    timEsp[i] = ESP;

    status = tim_Load[i](ESP, data);
    if (status < EGADS_SUCCESS) {
        printf("ERROR:: tim_load(%s) returned status=%d\n", timName, status);
        goto cleanup;
    }

    /* remember if we should run this tim in the same or a new thread */
    timHold[i] = status;

    /* mark the TIM as being ready */
    timState[i] = TIM_READY;
    if (SHOW_STATES) printf("in tim_load: setting TIM_READY\n");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_getEsp - get the ESP structure (needed by esp.py)             */
/*                                                                     */
/***********************************************************************/

/*@null@*/void *
tim_getEsp(char timName[])
{
    int  i;

    /* --------------------------------------------------------------- */

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("WARNING:: tim_getEsp(%s) is not currently loaded\n", timName);
        return NULL;
    }

//cleanup:
    return timEsp[i];
}


/***********************************************************************/
/*                                                                     */
/*   tim_mesg - process a message in the TIM (main execution method)   */
/*                                                                     */
/***********************************************************************/

int
tim_mesg(char   timName[],              /* (in)  name of TIM */
         char   command[])              /* (in)  command to process */
{
    int i, status=EGADS_SUCCESS;

    emp_T  *myEmp=NULL;

    ROUTINE(tim_mesg);

    /* --------------------------------------------------------------- */

//$$$    printf("tim_mesg(%s, %s)\n", timName, command);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_mesg(%s) is not currently loaded\n", timName);
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    /* wait here if a previous spawned process is still running */
    if (timThread[i] != NULL) {
        EMP_ThreadWait(timThread[i]);
    }

    /* make sure TIM is currently ready */
    if (timState[i] != TIM_READY) {
        printf("ERROR:: tim_mesg(%s) is not in ready, state=%d\n", timName, timState[i]);
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    /* destroy a thread that may be left over from previous call */
    if (timThread[i] != NULL) {
        EMP_ThreadDestroy(timThread[i]);
        timThread[i] = NULL;
    }

    timState[i] = TIM_EXECUTING;
    if (SHOW_STATES) printf("in tim_mesg: setting TIM_EXECUTING\n");

    /* if we are holding until released, call the message routine directly */
    if (timHold[i] == 1) {
        status = tim_Mesg[i](timEsp[i], command);
        if (status < EGADS_SUCCESS) {
            printf("ERROR:: tim_mesg(%s, %s) returned status=%d\n", timName, command, status);
            goto cleanup;
        }

        /* mark the TIM as being ready */
        timState[i] = TIM_READY;
        if (SHOW_STATES) printf("in tim_mesg: setting TIM_READY\n");

    /* otherwise we need to create a thread in which the message routine will run */
    } else {

        MALLOC(myEmp, emp_T, 1);

        myEmp->routine = tim_Mesg[i];
        myEmp->command = NULL;
        myEmp->data    = timEsp[i];
        /*@ignore@*/
        myEmp->state   = &(timState[i]);
        /*@end@*/

        MALLOC(myEmp->command, char, STRLEN(command)+1);
        strcpy(myEmp->command, command);

        timThread[i] = EMP_ThreadCreate(timExec, (void*)(myEmp));
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_hold - hold the TIM until the lock is released                */
/*                                                                     */
/***********************************************************************/

int
tim_hold(char   timName[],              /* (in)  name of TIM to hold */
         char   overlay[])              /* (in)  name of overlay to lift hold */
{
    int i, status=EGADS_SUCCESS;

    char  message[MAX_EXPR_LEN];
    esp_T *ESP=NULL;

    ROUTINE(tim_mesg);

    /* --------------------------------------------------------------- */

//$$$    printf("tim_hold(%s, %s)\n", timName, overlay);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_hold(%s) is not currently loaded\n", timName);
        goto cleanup;
    }

    /* send a message to the browsers to tell them which hold to lift */
    snprintf(message, MAX_EXPR_LEN, "overlayBeg|%s|%s|", timName, overlay);
    tim_bcst(timName, message);

    /* now wait until the browswer lifts the hold */
    EMP_LockSet(timMutex[i]);

    /* this thread should not set the lock, so it is immediately unlocked */
    EMP_LockRelease(timMutex[i]);

    /* wait for the lock to be set again */
    while (EMP_LockTest(timMutex[i]) == 0) {
        SLEEP(100);
    }

    /* reset the EGADS thread */
        /* update the thread using the context */
    ESP = (esp_T *)(timEsp[i]);
    if (ESP->MODL->context != NULL) {
        status = EG_updateThread(ESP->MODL->context);
        CHECK_STATUS(EG_updateThread);
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_lift - ask tim_lock to lift the hold                          */
/*                                                                     */
/***********************************************************************/

int
tim_lift(char   timName[])            /* (in)  name of TIM */
{
    int i, status=EGADS_SUCCESS;

    ROUTINE(tim_lift);

    /* --------------------------------------------------------------- */

//$$$    printf("tim_lift(%s)\n", timName);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_lift(%s) is not currently loaded\n", timName);
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    timUnset[i] = 1;

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_bcst - broadcast a message to all browsers                    */
/*                                                                     */
/***********************************************************************/

int
tim_bcst(char   timName[],              /* (in)  name of TIM */
         char   text[])                 /* (in)  text to broadcast */
{

    int    i, status=EGADS_SUCCESS;
    esp_T  *ESP=NULL;

    /* --------------------------------------------------------------- */

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_bcst(%s) is not currently loaded\n", timName);
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    ESP = (esp_T *) timEsp[i];

    /* only broadcast if wv is currently running */
    if (ESP->cntxt == NULL) {
        goto cleanup;
    }

    /* if we have a MODL and its outLevel >= 2, write tracing info to STDOUT */
    if (ESP->MODL != NULL) {
        if (ocsmSetOutLevel(-1) >= 2) {
            printf("\n<<< server2browser: %s\n", text);
        }
    }

    /* broadcast the message */
    wv_broadcastText(text);

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_save - save data associated with TIM and close it             */
/*                                                                     */
/***********************************************************************/

int
tim_save(char   timName[])              /* (in)  name of TIM */
{
    int    i, status;

    /* --------------------------------------------------------------- */

//$$$    printf("tim_save(%s)\n", timName);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_save(%s) is not currently loaded\n", timName);
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    /* wait here if a previous spawned process is still running */
    if (timThread[i] != NULL) {
        EMP_ThreadWait(timThread[i]);
    }

    /* make sure TIM is currently ready */
    if (timState[i] != TIM_READY) {
        printf("ERROR:: tim_save(%s) is not in ready, state=%d\n", timName, timState[i]);
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    timState[i] = TIM_CLOSING;
    if (SHOW_STATES) printf("in tim_save: setting TIM_CLOSING\n");

    /* call the save routine in the TIM */
    status = tim_Save[i](timEsp[i]);
    if (status < EGADS_SUCCESS) {
        printf("ERROR:: tim_save(%s) returned status=%d\n", timName, status);
        goto cleanup;
    }

    /* destroy the thread associated with this TIM */
    if (timThread[i] != NULL) {
        EMP_ThreadDestroy(timThread[i]);
        timThread[i] = NULL;
    }

    /* mark the TIM as being inactive */
    timState[i] = TIM_INACTIVE;
    if (SHOW_STATES) printf("in tim_save: setting TIM_INACTIVE\n");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_quit - remove data associated with TIM and close it           */
/*                                                                     */
/***********************************************************************/

int
tim_quit(char   timName[])              /* (in)  name of TIM */
{
    int   i, status;

    /* --------------------------------------------------------------- */

//$$$    printf("tim_quit(%s)\n", timName);

    /* determine if dynamic library is already loaded */
    i = timDLoaded(timName);

    /* if not loaded, return an error */
    if (i == -1) {
        printf("ERROR:: tim_quit(%s) is not currently loaded\n", timName);
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    /* wait here if a previous spawned process is still running */
    if (timThread[i] != NULL) {
        EMP_ThreadWait(timThread[i]);
    }

    /* make sure TIM is currently ready */
    if (timState[i] != TIM_READY) {
        printf("ERROR:: tim_quit(%s) is not in ready, state=%d\n", timName, timState[i]);
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    timState[i] = TIM_CLOSING;
    if (SHOW_STATES) printf("in tim_quit: setting TIM_CLOSING\n");

    /* call the quit routine in the TIM */
    status = tim_Quit[i](timEsp[i], 0);
    if (status < EGADS_SUCCESS) {
        printf("ERROR:: tim_quit(%s) returned status=%d\n", timName, status);
        goto cleanup;
    }

    /* destroy the thread associated with this TIM and reset te parent pointer */
    if (timThread[i] != NULL) {
        EMP_ThreadDestroy(timThread[i]);
        timThread[i] = NULL;
    }

    /* mark the TIM as being inactive */
    timState[i] = TIM_INACTIVE;
    if (SHOW_STATES) printf("in tim_quit: setting TIM_INACTIVE\n");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   tim_lock - set lock on all TIMs                                   */
/*                                                                     */
/***********************************************************************/

void
tim_lock()
{
    int i; //, ival;

    ROUTINE(tim_lock);

    /* --------------------------------------------------------------- */

    /* loop through the TIMs that are currently loaded */
    for (i = 0; i < tim_nTim; i++) {

        /* if the flag asks that the lock be released, do so now */
        if (timUnset[i] == 1) {
            EMP_LockRelease(timMutex[i]);

        /* if not set, set it now */
        } else if (EMP_LockTest(timMutex[i]) == 0) {
            EMP_LockSet(timMutex[i]);
        }

        timUnset[i] = 0;
    }

//cleanup:
    return;
}


/***********************************************************************/
/*                                                                     */
/*   tim_free - free all data associated with TIMs (called at end)     */
/*                                                                     */
/***********************************************************************/

void
tim_free()
{
    int i, status;

    /* --------------------------------------------------------------- */

//$$$    printf("tim_free()\n");

    /* loop through the TIMs that are currently loaded, and quit them all */
    for (i = 0; i < tim_nTim; i++) {

        /* free up the data associated with the TIM */
        if (timEsp[i] != NULL) {
            status = tim_Quit[i](timEsp[i], 1);
            if (status < EGADS_SUCCESS) {
                printf("ERROR:: tim_free(%s) quit returned status=%d\n", timName[i], status);
                goto cleanup;
            }
        }

        /* free up any threads */
        if (timThread[i] != NULL) {
            EMP_ThreadDestroy(timThread[i]);
            timThread[i] = NULL;
        }

        /* free up any mutexs */
        if (timMutex[i] != NULL) {
            EMP_LockDestroy(timMutex[i]);
            timMutex[i] = NULL;
        }

        /* remove the TIM name from the dispatch table */
        if (timName[i] != NULL) {
            free(timName[i]);
        }

        /* close the dynamically loaded routines */
        timDLclose(timDLL[i]);
    }

    /* we have no more active TIMs */
    tim_nTim = 0;

cleanup:
    return;
}


/***********************************************************************/
/*                                                                     */
/*   Utility Functions                                                 */
/*       timDLopen  - open a TIM                                       */
/*       timDLclose - close a TIM                                      */
/*       timDLget   - get function associated with a TIM               */
/*       timDLoaded - determine if a TIM is loaded                     */
/*       timDYNload - dynamically load a TIM                           */
/*       GetToken   - extract a token from a string                    */
/*                                                                     */
/***********************************************************************/

static /*@null@*/ DLL
timDLopen(const char *name)
{
    int             i, len;
    DLL             dll;
    char            *full, *env;
#ifdef WIN32
    WIN32_FIND_DATA ffd;
    char            dir[MAX_PATH];
    size_t          length_of_arg;
    HANDLE hFind =  INVALID_HANDLE_VALUE;
#else
    /*@-unrecog@*/
    char            dir[PATH_MAX];
    /*@+unrecog@*/
    struct dirent   *de;
    DIR             *dr;
#endif

    /* --------------------------------------------------------------- */

    env = getenv("ESP_ROOT");
    if (env == NULL) {
        printf("WARNING:: Could not find $ESP_ROOT\n");
        return NULL;
    }

    if (name == NULL) {
        printf("WARNING:: Dynamic Loader invoked with NULL name\n");
        return NULL;
    }
    len  = strlen(name);
    full = (char *) malloc((len+5)*sizeof(char));
    if (full == NULL) {
        printf("WARNING:: Dynamic Loader MALLOC Error\n");
        return NULL;
    }

    for (i = 0; i < len; i++) {
        full[i] = name[i];
    }
    full[len  ] = '.';
#ifdef WIN32
    full[len+1] = 'D';
    full[len+2] = 'L';
    full[len+3] = 'L';
    full[len+4] =  0;
    snprintf(dir, MAX_PATH, "%s\\lib\\*", env);
    hFind = FindFirstFile(dir, &ffd);
    if (INVALID_HANDLE_VALUE == hFind) {
        printf("WARNING:: Dynamic Loader could not open %s\n", dir);
        free(full);
        return NULL;
    }
    i = 0;
    do {
        if (strcasecmp(full, ffd.cFileName) == 0) i++;
    } while (FindNextFile(hFind, &ffd) != 0);
    FindClose(hFind);
    if (i > 1) {
        printf("WARNING:: Dynamic Loader more than 1 file: %s\n", full);
        free(full);
        return NULL;
    }

    if (i == 1) {
        hFind = FindFirstFile(dir, &ffd);
        do {
            if (strcasecmp(full, ffd.cFileName) == 0) break;
        } while (FindNextFile(hFind, &ffd) != 0);
        dll = LoadLibrary(ffd.cFileName);
        FindClose(hFind);
    } else {
        dll = LoadLibrary(full);
    }

#else

    full[len+1] = 's';
    full[len+2] = 'o';
    full[len+3] =  0;
    snprintf(dir, PATH_MAX, "%s/lib", env);
    dr = opendir(dir);
    if (dr == NULL) {
        printf("WARNING:: Dynamic Loader could not open %s\n", dir);
        free(full);
        return NULL;
    }
    i = 0;
    while ((de = readdir(dr)) != NULL)
        /*@-unrecog@*/
        if (strcasecmp(full, de->d_name) == 0) i++;
        /*@+unrecog@*/
    closedir(dr);
    if (i > 1) {
        printf("WARNING:: Dynamic Loader more than 1 file: %s\n", full);
        free(full);
        return NULL;
    }

    if (i == 1) {
        dr = opendir(dir);
        while ((de = readdir(dr)) != NULL)
            if (strcasecmp(full, de->d_name) == 0) break;
        /*@-nullderef@*//*@-unrecog@*/
        dll = dlopen(de->d_name, RTLD_NOW | RTLD_NODELETE | RTLD_GLOBAL);
        /*@+nullderef@*//*@+unrecog@*/
        closedir(dr);
    } else {
        /*@-nullderef@*//*@-unrecog@*/
        dll = dlopen(full, RTLD_NOW | RTLD_NODELETE | RTLD_GLOBAL);
        /*@+nullderef@*//*@+unrecog@*/
    }
#endif

    if (!dll) {
        printf("WARNING:: Dynamic Loader for %s not found\n", full);
#ifndef WIN32
        printf("              %s\n", dlerror());
#endif
        free(full);
        return NULL;
    }
    free(full);

    return dll;
}

/* ------------------------------------------------------------------------- */

static void
timDLclose(/*@unused@*/ /*@only@*/ DLL dll)
{
#ifdef WIN32
    FreeLibrary(dll);
#else
    dlclose(dll);
#endif
}

/* ------------------------------------------------------------------------- */

static timDLLfunc
timDLget(DLL dll,
         const char *symname,
         /*@null@*/ const char *name)
{
    timDLLfunc data;

    /* --------------------------------------------------------------- */

#ifdef WIN32
    data = (timDLLfunc) GetProcAddress(dll, symname);
#else
    /*@-castfcnptr@*/
    data = (timDLLfunc) dlsym(dll, symname);
    /*@+castfcnptr@*/
#endif
    if ((data == NULL) && (name != NULL))
        printf("WARNING:: Couldn't get symbol %s in %s\n", symname, name);
    return data;
}

/* ------------------------------------------------------------------------- */

static int
timDLoaded(const char *name)
{
    int i;

    /* --------------------------------------------------------------- */

    for (i = 0; i < tim_nTim; i++) {
        if (strcasecmp(name, timName[i]) == 0) {
            return i;
        }
    }

    return -1;
}

/* ------------------------------------------------------------------------- */

static int
timDYNload(const char *name)
{
    int i, len, ret;
    DLL dll;

    /* --------------------------------------------------------------- */

    if (tim_nTim >= MAX_TIMS) {
        printf("WARNING:: Number of Primitives %d > %d\n", tim_nTim, MAX_TIMS);
        return EGADS_INDEXERR;
    }
    dll = timDLopen(name);
    if (dll == NULL) return EGADS_NULLOBJ;

    ret = tim_nTim;
    tim_Load[ret] = (TIM_Load) timDLget(dll, "timLoad", name);
    tim_Mesg[ret] = (TIM_Mesg) timDLget(dll, "timMesg", name);
    tim_Save[ret] = (TIM_Save) timDLget(dll, "timSave", name);
    tim_Quit[ret] = (TIM_Quit) timDLget(dll, "timQuit", name);
    if ((tim_Load[ret] == NULL) ||
        (tim_Mesg[ret] == NULL) ||
        (tim_Save[ret] == NULL) ||
        (tim_Quit[ret] == NULL)   ) {
        timDLclose(dll);
        return EGADS_EMPTY;
    }

    len  = strlen(name) + 1;
    timName[ ret] = (char *) malloc(len*sizeof(char));
    if (timName[ret] == NULL) {
        printf("WARNING:: timName could not be malloced\n");
        timDLclose(dll);
        return EGADS_MALLOC;
    }
    for (i = 0; i < len; i++) {
        timName[ret][i] = name[i];
    }

    timDLL[   ret] = dll;
    timState[ ret] = TIM_INACTIVE;
    timEsp[   ret] = NULL;
    timThread[ret] = NULL;
    timMutex[ ret] = EMP_LockCreate();
    if (timMutex[ret] == NULL) {
        printf("WARNING:: timMutex could not be created\n");
        /*@-dependenttrans@*/
        timDLclose(dll);
        /*@+dependenttrans@*/
        free(timName[ret]);
        return EGADS_MALLOC;
    }
    timUnset[ ret] = 0;
    timHold[  ret] = 0;
    tim_nTim++;

    if (SHOW_STATES) printf("in timDYNload: setting TIM_INACTIVE\n");

    return ret;
}


/***********************************************************************/
/*                                                                     */
/*   GetToken - pull a token out of a string                           */
/*                                                                     */
/***********************************************************************/

int
GetToken(char   text[],                 /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   sep,                    /* (in)  separator character */
         char   *token[])               /* (out) token (freeable) */
{
    int    status = SUCCESS;

    int    lentok, i, count, iskip;

    ROUTINE(GetToken);

    /* --------------------------------------------------------------- */

    MALLOC(*token, char, strlen(text));

    (*token)[0] = '\0';
    lentok   = 0;

    /* convert tabs to spaces */
    for (i = 0; i < STRLEN(text); i++) {
        if (text[i] == '\t') {
            text[i] = ' ';
        }
    }

    /* count the number of separators */
    count = 0;
    for (i = 0; i < STRLEN(text); i++) {
        if (text[i] == sep) {
            count++;
        }
    }

    if (count < nskip) return 0;

    /* skip over nskip tokens */
    i = 0;
    for (iskip = 0; iskip < nskip; iskip++) {
        while (text[i] != sep) {
            i++;
        }
        i++;
    }

    /* if token we are looking for is empty, return 0 */
    if (text[i] == sep) {
        (*token)[0] = '0';
        (*token)[1] = '\0';
    }

    /* extract the token we are looking for */
    while (text[i] != sep && text[i] != '\0') {
        (*token)[lentok++] = text[i++];
        (*token)[lentok  ] = '\0';
    }

cleanup:
    if (status == SUCCESS) {
        return STRLEN(*token);
    } else {
        return status;
    }
}


/***********************************************************************/
/*                                                                     */
/*   timExec - routine that actually executes a message                */
/*                                                                     */
/***********************************************************************/

static void
timExec(void *s)
{
    emp_T  *myEmp = (emp_T *) s;

    ROUTINE(timExec);

    /* --------------------------------------------------------------- */

    /* execute timMesg */
    myEmp->routine((esp_T *)(myEmp->data), myEmp->command);

    /* we are now redy */
    *(myEmp->state) = TIM_READY;
    if (SHOW_STATES) printf("in timExec: setting TIM_READY\n");

    /* free up the emp_T structure (that was allocated in tim_mesg) */
    FREE(myEmp->command);
    FREE(myEmp);
}
