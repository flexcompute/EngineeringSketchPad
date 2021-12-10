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
 * Copyright (C) 2010/2021  John F. Dannenhoffer, III (Syracuse University)
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

#ifdef WIN32
   #define DLL        HINSTANCE
   #define WIN32_LEAN_AND_MEAN
   #include <windows.h>
   #include <winsock2.h>
   #define snprintf   _snprintf
   #define strcasecmp stricmp
#else
   #define DLL void *
   #include <dlfcn.h>
   #include <dirent.h>
   #include <limits.h>
#endif

#include "egads.h"
#include "OpenCSM.h"
#include "common.h"
#include "tim.h"

#include "wsserver.h"

#define MAX_TIMS 32

typedef int (*timDLLfunc) (void);

static char      *timName[MAX_TIMS];
static DLL        timDLL[ MAX_TIMS];
static void      *timData[MAX_TIMS];

// func pointers  name in DLL       name in OpenCSM     description
// -------------  ---------------   ------------------  -----------------------------------------------

// tim_Load[]     timLoad()         tim_load()          open a tim
typedef int (*TIM_Load) (esp_T *ESP, /*@null@*/void *data);
static        TIM_Load  tim_Load[MAX_TIMS];

// tim_Save[]     timSave()         tim_save()          save tim data and close tim instance
typedef int (*TIM_Save) (esp_T *ESP);
static        TIM_Save  tim_Save[MAX_TIMS];

// tim_Quit[]     timQuit()         tim_quit()          close tim instance without saving
typedef int (*TIM_Quit) (esp_T *ESP);
static        TIM_Quit  tim_Quit[MAX_TIMS];

// tim_Mesg[]     timMesg()         tim_mesg()          get command, process, and return response
typedef int (*TIM_Mesg) (esp_T *ESP, char command[], int *len_response, char *response[]);
static        TIM_Mesg  tim_Mesg[MAX_TIMS];

//                                  tim_free()          frees up memory associated with TIMs

static int tim_nTim = 0;


/***********************************************************************/
/*                                                                     */
/*   Utility Functions                                                 */
/*       timDLopen  - open a TIM                                       */
/*       timDLclose - close a TIM                                      */
/*       timDLget   - get function associated with a TIM               */
/*       timDLoaded - determine if a TIM is loaded                     */
/*       timDYNload - dynamically load a TIM                           */
/*                                                                     */
/***********************************************************************/

static /*@null@*/ DLL timDLopen(const char *name)
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
        printf(" Information: Could not find $ESP_ROOT\n");
        return NULL;
    }

    if (name == NULL) {
        printf(" Information: Dynamic Loader invoked with NULL name!\n");
        return NULL;
    }
    len  = strlen(name);
    full = (char *) malloc((len+5)*sizeof(char));
    if (full == NULL) {
        printf(" Information: Dynamic Loader MALLOC Error!\n");
        return NULL;
    }
    for (i = 0; i < len; i++) full[i] = name[i];
    full[len  ] = '.';
#ifdef WIN32
    full[len+1] = 'D';
    full[len+2] = 'L';
    full[len+3] = 'L';
    full[len+4] =  0;
    snprintf(dir, MAX_PATH, "%s\\lib\\*", env);
    hFind = FindFirstFile(dir, &ffd);
    if (INVALID_HANDLE_VALUE == hFind) {
        printf(" Information: Dynamic Loader could not open %s\n", dir);
        free(full);
        return NULL;
    }
    i = 0;
    do {
        if (strcasecmp(full, ffd.cFileName) == 0) i++;
    } while (FindNextFile(hFind, &ffd) != 0);
    FindClose(hFind);
    if (i > 1) {
        printf(" Information: Dynamic Loader more than 1 file: %s!\n", full);
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
        printf(" Information: Dynamic Loader could not open %s\n", dir);
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
        printf(" Information: Dynamic Loader more than 1 file: %s!\n", full);
        free(full);
        return NULL;
    }

    if (i == 1) {
        dr = opendir(dir);
        while ((de = readdir(dr)) != NULL)
            if (strcasecmp(full, de->d_name) == 0) break;
/*@-nullderef@*/
        dll = dlopen(de->d_name, RTLD_NOW | RTLD_NODELETE);
/*@+nullderef@*/
        closedir(dr);
    } else {
        dll = dlopen(full, RTLD_NOW | RTLD_NODELETE);
    }
#endif

    if (!dll) {
        printf(" Information: Dynamic Loader for %s not found\n", full);
#ifndef WIN32
/*@-mustfreefresh@*/
        printf("              %s\n", dlerror());
/*@+mustfreefresh@*/
#endif
        free(full);
/*@-mustfreefresh@*/
        return NULL;
/*@+mustfreefresh@*/
    }
    free(full);

/*@-mustfreefresh@*/
    return dll;
/*@+mustfreefresh@*/
}

/* ------------------------------------------------------------------------- */

static void timDLclose(/*@unused@*/ /*@only@*/ DLL dll)
{
#ifdef WIN32
    FreeLibrary(dll);
#else
    dlclose(dll);
#endif
}

/* ------------------------------------------------------------------------- */

static timDLLfunc timDLget(DLL dll,   const char *symname,
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
        printf(" Information: Couldn't get symbol %s in %s\n", symname, name);
    return data;
}

/* ------------------------------------------------------------------------- */

static int timDLoaded(const char *name)
{
    int i;

    /* --------------------------------------------------------------- */

    for (i = 0; i < tim_nTim; i++)
        if (strcasecmp(name, timName[i]) == 0) return i;

    return -1;
}

/* ------------------------------------------------------------------------- */

static int timDYNload(const char *name)
{
    int i, len, ret;
    DLL dll;

    /* --------------------------------------------------------------- */

    if (tim_nTim >= MAX_TIMS) {
        printf(" Information: Number of Primitives > %d!\n", MAX_TIMS);
        return EGADS_INDEXERR;
    }
    dll = timDLopen(name);
    if (dll == NULL) return EGADS_NULLOBJ;

    ret = tim_nTim;
    tim_Load[ret] = (TIM_Load) timDLget(dll, "timLoad", name);
    tim_Save[ret] = (TIM_Save) timDLget(dll, "timSave", name);
    tim_Quit[ret] = (TIM_Quit) timDLget(dll, "timQuit", name);
    tim_Mesg[ret] = (TIM_Mesg) timDLget(dll, "timMesg", name);
    if ((tim_Load[ret] == NULL) ||
        (tim_Save[ret] == NULL) ||
        (tim_Quit[ret] == NULL) ||
        (tim_Mesg[ret] == NULL)   ) {
        timDLclose(dll);
        return EGADS_EMPTY;
    }

    len  = strlen(name) + 1;
    timName[ret] = (char *) malloc(len*sizeof(char));
    timData[ret] = NULL;
    if (timName[ret] == NULL) {
        timDLclose(dll);
        return EGADS_MALLOC;
    }
    for (i = 0; i < len; i++) timName[ret][i] = name[i];
    timDLL[ret] = dll;
    tim_nTim++;

    return ret;
}


/***********************************************************************/
/*                                                                     */
/*   Exposed (to serveESP) Functions                                   */
/*       tim_load   - dynamically load a TIM                           */
/*       tim_save   - leave a TIM and save  its data                   */
/*       tim_quit   - leave a TIM and loose its data                   */
/*       tim_mesg   - send a message to a TIM                          */
/*       tim_data   - get header data associated with a TIM            */
/*       tim_free   - free up data associated with a TIM               */
/*                                                                     */
/***********************************************************************/

int
tim_load(char       timName[],
         esp_T      *ESP,
/*@null@*/void      *data)
{
    int i, status;

    /* --------------------------------------------------------------- */

    i = timDLoaded(timName);

    if (i == -1) {
        i = timDYNload(timName);
        if (i < 0) {
            status = i;
            goto cleanup;
        }
    }

    timData[i] = ESP;

    status = tim_Load[i](ESP, data);

cleanup:
    return status;
}

/* ------------------------------------------------------------------------- */

int
tim_save(char   timName[])
{
    int i, status;

    /* --------------------------------------------------------------- */

    i = timDLoaded(timName);
    if (i == -1) {
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    status = tim_Save[i](timData[i]);

cleanup:
    return status;
}

/* ------------------------------------------------------------------------- */

int
tim_quit(char   timName[])
{
    int i, status;

    /* --------------------------------------------------------------- */

    i = timDLoaded(timName);
    if (i == -1) {
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    status = tim_Quit[i](timData[i]);

cleanup:
    return status;
}

/* ------------------------------------------------------------------------- */

int
tim_mesg(char   timName[],
         char   command[],
         int    *len_response,
         char   *response[])
{
    int i, status;

    /* --------------------------------------------------------------- */

    i = timDLoaded(timName);
    if (i == -1) {
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    status = tim_Mesg[i](timData[i], command, len_response, response);

cleanup:
    return status;
}

/* ------------------------------------------------------------------------- */

int
tim_data(char   timName[],
         esp_T  **ESP)
{
    int i, status;

    /* --------------------------------------------------------------- */

    i = timDLoaded(timName);
    if (i == -1) {
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    *ESP = (esp_T *)(timData[i]);

    status = EGADS_SUCCESS;

cleanup:
    return status;

}

/* ------------------------------------------------------------------------- */

void
tim_free()
{
    int i;

    /* --------------------------------------------------------------- */

    if (tim_nTim == 0) return;

    for (i = 0; i < tim_nTim; i++) {
        if (timData[i] != NULL) {
            tim_Quit[i](timData[i]);
        }
        if (timName[i] != NULL) {
            free(timName[i]);
        }
        timDLclose(timDLL[i]);
    }

    tim_nTim = 0;
}


/***********************************************************************/
/*                                                                     */
/*   Other support Functions                                           */
/*       GetToken   - extract a token from a string                    */
/*                                                                     */
/***********************************************************************/

int
GetToken(char   text[],                 /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   sep,                    /* (in)  separator character */
         char   token[])                /* (out) token */
{
    int    lentok, i, count, iskip;

    /* --------------------------------------------------------------- */

    token[0] = '\0';
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

    if (count < nskip+1) return 0;

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
        token[0] = '0';
        token[1] = '\0';
    }

    /* extract the token we are looking for */
    while (text[i] != sep) {
        token[lentok++] = text[i++];
        token[lentok  ] = '\0';

        if (lentok >= MAX_EXPR_LEN-1) {
            break;
        }
    }

    return STRLEN(token);
}
