/*
 ************************************************************************
 *                                                                      *
 * OpenCSM/EGADS -- Extended User Defined Primitive/Function Interface  *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2024  John F. Dannenhoffer, III (Syracuse University)
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
   #include <windows.h>
   #define snprintf   _snprintf
   #define strcasecmp stricmp
#else
   #define DLL void *
   #include <dlfcn.h>
   #include <dirent.h>
   #include <limits.h>
#endif

#include "egads.h"
#include "udp.h"

typedef int (*udpDLLfunc) (void);
static char *udpName[MAXPRIM];
static DLL   udpDLL[ MAXPRIM];

// func pointers  name in DLL       name in OpenCSM     description
// -------------  ---------------   ------------------  -----------------------------------------------

// ----- defined in udpUtilities.c -----

// udpInit[]      udpInitialize()   udp_initialize()    initialize and get info about list of arguments
typedef int (*udpI) (int *nArgs, char **namex[], int *typex[], int *idefault[], double *ddefault[], udp_T *Udps[]);
static        udpI  udpInit[ MAXPRIM];

// udpNumB[]      udpNumBodys()     udp_numBodys()      number of Bodys expected in first arg to udpExecute
typedef int (*udpN) ();
static        udpN  udpNumB[ MAXPRIM];

// udpBodyL[]     udpBodyList()     udp_bodyList()      list of Bodys input to a UDF
typedef int (*udpB) (ego, /*@null@*/const int**, int numudp, udp_T Udps[]);
static        udpB  udpBodyL[MAXPRIM];

// udpReset[]     udpReset()        udp_clrArguments()  reset the argument to their defaults
typedef int (*udpR) (/*@null@*/int *NumUdp, /*@null@*/udp_T Udps[]);
static        udpR  udpReset[MAXPRIM];

// udpFree[]      udpFree()        udp_free()           free memory associated with UDP/UDFs
typedef int (*udpF) (int numudp, udp_T Udps[]);
static        udpF  udpFree[MAXPRIM];

// udpClean[]     udpClean()        udp_clean()         clean the udp cache
typedef int (*udpC) (int *NumUdp, udp_T Udps[]);
static        udpC  udpClean[MAXPRIM];

// udpSet[]       udpSetArgument()  udp_setArgument()   sets value of an argument
typedef int (*udpS) (char name[], void *value, int nrow, int ncol, /*@null@*/char message[], udp_T Udps[]);
static        udpS  udpSet[MAXPRIM];

// udpGet[]       udpGet()          udp_getOutput()     return an output parameter
typedef int (*udpG) (ego ebody, char name[], int *nrow, int *ncol, void **val, void **dot, char message[], int numudp, udp_T Udps[]);
static        udpG  udpGet[MAXPRIM];

// udpMesh[]      udpMesh()         udp_getMesh()       return mesh associated with primitive
typedef int (*udpM) (ego ebody, int iMesh, int *imax, int *jmax, int *kmax, double *mesh[], int *NumUdp, udp_T Udps[]);
static        udpM  udpGrid[MAXPRIM];

// udpVel[]       udpVel()          udp_setVelocity()   set velocity of an argument
typedef int (*udpV) (ego ebody, char name[], double dot[], int ndot, int numudp, udp_T Udps[]);
static        udpV  udpVel[MAXPRIM];

// udpPost[]      udpPost()         udp_postSens()      reset ndotchg flag */
typedef int (*udpO) (ego ebody, int numudp, udp_T Udps[]);
static        udpO  udpPost[MAXPRIM];

// ----- defined in udpXXX -----

// udpExec[]      udpExecute()      udp_executePrim()   execute the primitive
typedef int (*udpE) (ego context, ego *ebody, int *nMesh, char *string[], int *NumUdp, udp_T *Udps[]);
static        udpE  udpExec[MAXPRIM];

// udpSens[]      udpSensitivity    udp_sensitivity()   return sensitivity derivatives for the "real" argument
typedef int (*udpP) (ego ebody, int npnt, int entType, int entIndex, /*@null@*/double uvs[], double vels[], int *NumUdp, udp_T Udps[]);
static        udpP  udpSens[MAXPRIM];

static int udp_nPrim = 0;

/* ************************* Utility Functions ***************************** */

static /*@null@*/ DLL udpDLopen(const char *name)
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


static void udpDLclose(/*@unused@*/ /*@only@*/ DLL dll)
{
#ifdef WIN32
    FreeLibrary(dll);
#else
    dlclose(dll);
#endif
}


static udpDLLfunc udpDLget(DLL dll,   const char *symname,
                        /*@null@*/ const char *name)
{
    udpDLLfunc data;

#ifdef WIN32
    data = (udpDLLfunc) GetProcAddress(dll, symname);
#else
/*@-castfcnptr@*/
    data = (udpDLLfunc) dlsym(dll, symname);
/*@+castfcnptr@*/
#endif
    if ((data == NULL) && (name != NULL))
        printf(" Information: Couldn't get symbol %s in %s\n", symname, name);
    return data;
}


static int udpDLoaded(const char *name)
{
    int i;

    for (i = 0; i < udp_nPrim; i++)
        if (strcasecmp(name, udpName[i]) == 0) return i;

    return -1;
}


static int udpDYNload(const char *name)
{
    int i, len, ret;
    DLL dll;

    if (udp_nPrim >= MAXPRIM) {
        printf(" Information: Number of Primitives > %d!\n", MAXPRIM);
        return EGADS_INDEXERR;
    }
    dll = udpDLopen(name);
    if (dll == NULL) return EGADS_NULLOBJ;

    ret = udp_nPrim;
    udpInit[ ret] = (udpI) udpDLget(dll, "udpInitialize",  name);
    udpNumB[ ret] = (udpN) udpDLget(dll, "udpNumBodys",    name);
    udpBodyL[ret] = (udpB) udpDLget(dll, "udpBodyList",    name);
    udpReset[ret] = (udpR) udpDLget(dll, "udpReset",       name);
    udpFree[ ret] = (udpF) udpDLget(dll, "udpFree",        name);
    udpClean[ret] = (udpC) udpDLget(dll, "udpClean",       name);
    udpSet[  ret] = (udpS) udpDLget(dll, "udpSet",         name);
    udpExec[ ret] = (udpE) udpDLget(dll, "UdpExecute",     name);
    udpGet[  ret] = (udpG) udpDLget(dll, "udpGet",         NULL);
    udpGrid[ ret] = (udpM) udpDLget(dll, "UdpMesh",        NULL);
    udpVel[  ret] = (udpV) udpDLget(dll, "udpVel",         name);
    udpPost[ ret] = (udpO) udpDLget(dll, "udpPost",        name);
    udpSens[ ret] = (udpP) udpDLget(dll, "UdpSensitivity", name);
    if ((udpInit[ ret] == NULL) ||
        (udpNumB[ ret] == NULL) ||
        (udpBodyL[ret] == NULL) ||
        (udpReset[ret] == NULL) ||
        (udpFree[ ret] == NULL) ||
        (udpClean[ret] == NULL) ||
        (udpSet[  ret] == NULL) ||
        (udpExec[ ret] == NULL) ||
        (udpVel[  ret] == NULL) ||
        (udpPost[ ret] == NULL) ||
        (udpSens[ ret] == NULL)   ) {
        udpDLclose(dll);
        return EGADS_EMPTY;
    }

    len  = strlen(name) + 1;
    udpName[ret] = (char *) malloc(len*sizeof(char));
    if (udpName[ret] == NULL) {
        udpDLclose(dll);
        return EGADS_MALLOC;
    }
    for (i = 0; i < len; i++) udpName[ret][i] = name[i];
    udpDLL[ret] = dll;
    udp_nPrim++;

    return ret;
}


/* ************************* Exposed Functions ***************************** */


int udp_initialize(const char primName[],
                   modl_T     *MODL,
                   int        *nArgs,
                   char       **name[],
                   int        *type[],
                   int        *idefault[],
                   double     *ddefault[])
{
    int i;

    if (UDP_TRACE) printf("udp_initialize(primName=%s)\n", primName);

    i = udpDLoaded(primName);

    if (i == -1) {
        i = udpDYNload(primName);
        if (i < 0) return i;
    }

    return udpInit[i](nArgs, name, type, idefault, ddefault, &(MODL->Udps[i]));
}


int udp_numBodys(const char primName[])
{
    int i;

    if (UDP_TRACE) printf("udp_numBodys(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

/*@-nullpass@*/
    return udpNumB[i]();
/*@+nullpass@*/
}


int udp_bodyList(const char primName[],
                 modl_T     *MODL,
                 ego        body,
                 const int  **bodyList)
{
    int i;

    if (UDP_TRACE) printf("udp_bodyList(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

/*@-nullpass@*/
    return udpBodyL[i](body, bodyList, MODL->NumUdp[i], MODL->Udps[i]);
/*@+nullpass@*/
}


int udp_clrArguments(const char primName[],
                     modl_T     *MODL)
{
    int i;

    if (UDP_TRACE) printf("udp_clrArguments(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpReset[i](&(MODL->NumUdp[i]), MODL->Udps[i]);
}


int udp_clean(const char primName[],
              modl_T     *MODL)
{
    int i;

    if (UDP_TRACE) printf("udp_clean(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpClean[i](&(MODL->NumUdp[i]), MODL->Udps[i]);
}


int udp_setArgument(const char primName[],
                    modl_T     *MODL,
                    char       name[],
                    void       *value,
                    int        nrow,
                    int        ncol,
          /*@null@*/char       message[])
{
    int i;

    if (UDP_TRACE) printf("udp_setArgument(primName=%s, name=%s, nrow=%d, ncol=%d)\n", primName, name, nrow, ncol);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpSet[i](name, value, nrow, ncol, message, MODL->Udps[i]);
}


int udp_free(modl_T     *MODL)
{
    int i, status;

    if (UDP_TRACE) printf("udp_free()\n");

    for (i = 0; i < udp_nPrim; i++) {
        status = udpFree[i](MODL->NumUdp[i], MODL->Udps[i]);
        if (status != EGADS_SUCCESS) return status;
    }

    return EGADS_SUCCESS;
}


int udp_executePrim(const char primName[],
                    modl_T     *MODL,
                    ego        context,
                    ego        *body,
                    int        *nMesh,
                    char       *string[])
{
    int i;

    if (UDP_TRACE) printf("udp_executePrim(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpExec[i](context, body, nMesh, string, &(MODL->NumUdp[i]), &(MODL->Udps[i]));
}


int udp_getOutput(const char primName[],
                  modl_T     *MODL,
                  ego        body,
                  char       name[],
                  int        *nrow,
                  int        *ncol,
                  void       *val[],
                  void       *dot[],
                  char       message[])
{
    int i;

    if (UDP_TRACE) printf("udp_getOutput(primName=%s, body=%llx)\n", primName, (long long)body);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;
    if (udpGet[i] == NULL) return EGADS_EMPTY;

    return udpGet[i](body, name, nrow, ncol, val, dot, message, MODL->NumUdp[i], MODL->Udps[i]);
}


int udp_getMesh(const char primName[],
                modl_T     *MODL,
                ego        body,
                int        imesh,
                int        *imax,
                int        *jmax,
                int        *kmax,
                double     *mesh[])
{
    int i;

    if (UDP_TRACE) printf("udp_getMesh(primName=%s, body=%llx\n", primName, (long long)body);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;
    if (udpGrid[i] == NULL) return EGADS_EMPTY;

    return udpGrid[i](body, imesh, imax, jmax, kmax, mesh, &(MODL->NumUdp[i]), MODL->Udps[i]);
}


int udp_setVelocity(const char primName[],
                    modl_T     *MODL,
                    ego        body,
                    char       name[],
                    double     value[],
                    int        nvalue)
{
    int i;

    if (UDP_TRACE) printf("udp_setVelocity(primName=%s, body=%llx, name=%s, nvalue=%d)\n", primName, (long long)body, name, nvalue);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpVel[i](body, name, value, nvalue, MODL->NumUdp[i], MODL->Udps[i]);
}


int udp_sensitivity(const char primName[],
                    modl_T     *MODL,
                    ego        body,
                    int        npts,
                    int        entType,
                    int        entIndex,
                    double     uvs[],
                    double     vels[])
{
    int i;

    if (UDP_TRACE) printf("udp_sensitivity(primName=%s, body=%llx, npts=%d, entType=%d, entIndex=%d)\n", primName, (long long)body, npts, entType, entIndex);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpSens[i](body, npts, entType, entIndex, uvs, vels, &(MODL->NumUdp[i]), MODL->Udps[i]);
}


int udp_postSens(const char primName[],
                 modl_T     *MODL,
                 ego        body)
{
    int i;

    if (UDP_TRACE) printf("udp_postSens(primName=%s)\n", primName);

    i = udpDLoaded(primName);
    if (i == -1) return EGADS_NOTFOUND;

    return udpPost[i](body, MODL->NumUdp[i], MODL->Udps[i]);
}


void udp_cleanupAll()
{
    int i;

    if (UDP_TRACE) printf("udp_cleanupAll()\n");

    if (udp_nPrim == 0) return;

    for (i = 0; i < udp_nPrim; i++) {
        free(udpName[i]);
        udpDLclose(udpDLL[i]);
    }

    udp_nPrim = 0;
}
