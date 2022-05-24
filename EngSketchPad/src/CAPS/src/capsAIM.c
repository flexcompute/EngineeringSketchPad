/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Dynamic Subsystem
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "capsTypes.h"
#include "capsErrors.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp

#else

#include <dlfcn.h>
#include <dirent.h>
#include <limits.h>
#endif

extern void caps_freeDiscr(capsDiscr *discr);


typedef int (*DLLFunc) (void);


/* ************************* Utility Functions ***************************** */

static /*@null@*/ DLL aimDLopen(const char *name)
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
/*@-nullderef@*/
    closedir(dr);
  } else {
    dll = dlopen(full, RTLD_NOW | RTLD_NODELETE);
  }
#endif
  
  if (!dll) {
    printf(" Information: Dynamic Loader Error for %s\n", full);
#ifndef WIN32
/*@-mustfreefresh@*/
    printf("              %s\n", dlerror());
/*@+mustfreefresh@*/
#endif
    free(full);
    return NULL;
  }
  free(full);
  
  return dll;
}


static void aimDLclose(/*@only@*/ DLL dll)
{
#ifdef WIN32
  FreeLibrary(dll);
#else
  dlclose(dll);
#endif
}


static DLLFunc aimDLget(DLL dll, const char *symname)
{
  DLLFunc data;
  
#ifdef WIN32
  data = (DLLFunc) GetProcAddress(dll, symname);
#else
/*@-castfcnptr@*/
  data = (DLLFunc) dlsym(dll, symname);
/*@+castfcnptr@*/
#endif
  
  return data;
}


static int aimDLoaded(aimContext cntxt, const char *name)
{
  int i;
  
  for (i = 0; i < cntxt.aim_nAnal; i++)
    if (strcasecmp(name, cntxt.aimName[i]) == 0) return i;
  
  return -1;
}


static int aimDYNload(aimContext *cntxt, const char *name)
{
  int i, len, ret;
  DLL dll;
  
  if (cntxt->aim_nAnal >= MAXANAL) {
    printf(" Information: Number of AIMs > %d!\n", MAXANAL);
    return EGADS_INDEXERR;
  }
  dll = aimDLopen(name);
  if (dll == NULL) return EGADS_NULLOBJ;
  
  ret                     = cntxt->aim_nAnal;
  cntxt->aimInit[ret]     = (aimI)  aimDLget(dll, "aimInitialize"    );
  cntxt->aimDiscr[ret]    = (aimD)  aimDLget(dll, "aimDiscr"         );
  cntxt->aimFreeD[ret]    = (aimF)  aimDLget(dll, "aimFreeDiscrPtr"  );
  cntxt->aimLoc[ret]      = (aimL)  aimDLget(dll, "aimLocateElement" );
  cntxt->aimInput[ret]    = (aimIn) aimDLget(dll, "aimInputs"        );
  cntxt->aimPAnal[ret]    = (aimA)  aimDLget(dll, "aimPreAnalysis"   );
  cntxt->aimExec[ret]     = (aimEx) aimDLget(dll, "aimExecute"       );
#ifdef ASYNCEXEC
  cntxt->aimCheck[ret]    = (aimEx) aimDLget(dll, "aimCheck"         );
#endif
  cntxt->aimPost[ret]     = (aimPo) aimDLget(dll, "aimPostAnalysis"  );
  cntxt->aimOutput[ret]   = (aimO)  aimDLget(dll, "aimOutputs"       );
  cntxt->aimCalc[ret]     = (aimC)  aimDLget(dll, "aimCalcOutput"    );
  cntxt->aimXfer[ret]     = (aimT)  aimDLget(dll, "aimTransfer"      );
  cntxt->aimIntrp[ret]    = (aimP)  aimDLget(dll, "aimInterpolation" );
  cntxt->aimIntrpBar[ret] = (aimP)  aimDLget(dll, "aimInterpolateBar");
  cntxt->aimIntgr[ret]    = (aimG)  aimDLget(dll, "aimIntegration"   );
  cntxt->aimIntgrBar[ret] = (aimG)  aimDLget(dll, "aimIntegrateBar"  );
  cntxt->aimBdoor[ret]    = (aimBd) aimDLget(dll, "aimBackdoor"      );
  cntxt->aimClean[ret]    = (aimCU) aimDLget(dll, "aimCleanup"       );
  if ((cntxt->aimInit[ret]   == NULL) || (cntxt->aimClean[ret] == NULL) ||
      (cntxt->aimInput[ret]  == NULL) || (cntxt->aimPAnal[ret] == NULL) ||
      (cntxt->aimOutput[ret] == NULL) || (cntxt->aimCalc[ret]  == NULL) ||
      (cntxt->aimPost[ret]   == NULL)) {
    aimDLclose(dll);
    if (cntxt->aimInit[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimInitialize' in %s\n", name);
    if (cntxt->aimInput[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimInputs' in %s\n", name);
    if (cntxt->aimOutput[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimOutputs' in %s\n", name);
    if (cntxt->aimPAnal[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimPreAnalysis' in %s\n", name);
    if (cntxt->aimPost[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimPostAnalysis' in %s\n", name);
    if (cntxt->aimCalc[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimCalcOutput' in %s\n", name);
    if (cntxt->aimClean[ret]   == NULL)
      printf(" Error: Missing non-optional symbol 'aimCleanup' in %s\n", name);
    return EGADS_EMPTY;
  }
  
  len  = strlen(name) + 1;
  cntxt->aimName[ret] = (char *) malloc(len*sizeof(char));
  if (cntxt->aimName[ret] == NULL) {
    aimDLclose(dll);
    return EGADS_MALLOC;
  }
  for (i = 0; i < len; i++) cntxt->aimName[ret][i] = name[i];
  cntxt->aimDLL[ret]    = dll;
  cntxt->aim_nInst[ret] = 0;
  cntxt->aim_nAnal++;
  
  return ret;
}


/* ************************* Exposed Functions ***************************** */

int
aim_Initialize(aimContext *cntxt,
               const char *analysisName,
               int        *qeFlag,      /* query/execute flag */
    /*@null@*/ const char *unitSys,     /* Unit System requested */
    /*@null@*/ void       *aimStruc,    /* the AIM context */
               int        *major,       /* the returned major version */
               int        *minor,       /* the returned minor version */
               int        *nIn,         /* returned number of inputs */
               int        *nOut,        /* returned number of outputs */
               int        *nField,      /* returned number of DataSet fields */
               char       ***fnames,    /* returned pointer to field strings */
               int        **franks,     /* returned pointer to field ranks */
               int        **fInOut,     /* returned pointer to field in/out */
               void       **instStore)  /* returned instance storage */
{
  int i, ninst, stat, qFlag;
  
  i = aimDLoaded(*cntxt, analysisName);
  
  if (i == -1) {
    i = aimDYNload(cntxt, analysisName);
    if (i < 0) return i;
  }
  qFlag   = *qeFlag;
  *qeFlag = 1;
  if (cntxt->aimExec[i] == NULL) *qeFlag = 0;

  if (qFlag != 0) {
    ninst = -1;
  } else {
    ninst = cntxt->aim_nInst[i];
  }
  
  stat = cntxt->aimInit[i](ninst, unitSys, aimStruc, instStore,
                           major, minor, nIn, nOut,
                           nField, fnames, franks, fInOut);
  if ((qFlag != 0) || (stat != CAPS_SUCCESS)) return stat;

  cntxt->aim_nInst[i]++;
  return ninst;
}


int
aim_Index(aimContext cntxt, const char *analysisName)
{
  return aimDLoaded(cntxt, analysisName);
}


int
aim_Discr(aimContext cntxt,
          const char *analysisName,
          char       *bname,            /* Bound name */
          capsDiscr  *discr)            /* the structure to fill */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                 == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimDiscr[i] == NULL) {
    printf("aimDiscr not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimDiscr[i](bname, discr);
}


int
aim_FreeDiscr(aimContext cntxt,
              const char *analysisName,
              capsDiscr  *discr)        /* the structure to free up */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  if (cntxt.aimDiscr[i] != NULL)
    if (discr->ptrm != NULL) {
      cntxt.aimFreeD[i](discr->ptrm);
      discr->ptrm = NULL;
    }
  caps_freeDiscr(discr);
  
  return CAPS_SUCCESS;
}


int
aim_LocateElement(aimContext cntxt,
                  const char *analysisName,
                  capsDiscr  *discr,    /* the input discrete structure */
                  double     *params,   /* the input global parametric space */
                  double     *param,    /* the requested position */
                  int        *bIndex,   /* the returned body index */
                  int        *eIndex,   /* the returned element index */
                  double     *bary)     /* the barycentric coordinates */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i               == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimLoc[i] == NULL) {
    printf("aimLocateElement not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimLoc[i](discr, params, param, bIndex, eIndex, bary);
}


int
aim_LocateElIndex(aimContext cntxt,
                  int        i,
                  capsDiscr  *discr,    /* the input discrete structure */
                  double     *params,   /* the input global parametric space */
                  double     *param,    /* the requested position */
                  int        *bIndex,   /* the returned body index */
                  int        *eIndex,   /* the returned element index */
                  double     *bary)     /* the barycentric coordinates */
{
  if ((i < 0) || (i >= cntxt.aim_nAnal)) return EGADS_RANGERR;
  if (cntxt.aimLoc[i] == NULL) {
    printf("aimLocateElement not implemented in AIM %s\n", cntxt.aimName[i]);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimLoc[i](discr, params, param, bIndex, eIndex, bary);
}


int
aim_Inputs(aimContext cntxt,
           const char *analysisName,
/*@null@*/ void       *instStore,       /* instance storage */
           void       *aimStruc,        /* the AIM context */
           int        index,            /* the input index [1-nIn] */
           char       **ainame,         /* pointer to the returned name */
           capsValue  *defaultVal)      /* pointer to default value (filled) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  return cntxt.aimInput[i](instStore, aimStruc, index, ainame, defaultVal);
}


int
aim_PreAnalysis(aimContext cntxt,
                const char *analysisName,
     /*@null@*/ void       *instStore,  /* instance storage */
                void       *aimStruc,   /* the AIM context */
     /*@null@*/ capsValue  *inputs)     /* complete suite of analysis inputs */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  return cntxt.aimPAnal[i](instStore, aimStruc, inputs);
}


int
aim_Execute(aimContext cntxt,
            const char *analysisName,
 /*@null@*/ void       *instStore,  /* instance storage */
            void       *aimStruc,   /* the AIM context */
            int        *state)      /* the returned state of the execution */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  if (cntxt.aimExec[i] == NULL) {
    printf("aimExecute not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimExec[i](instStore, aimStruc, state);
}


#ifdef ASYNCEXEC
int
aim_Check(aimContext cntxt,
          const char *analysisName,
/*@null@*/void       *instStore,    /* instance storage */
          void       *aimStruc,     /* the AIM context */
          int        *state)        /* the returned state of the execution */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  if (cntxt.aimCheck[i] == NULL) {
    printf("aimCheck not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimCheck[i](instStore, aimStruc, state);
}
#endif


int
aim_Outputs(aimContext cntxt,
            const char *analysisName,
 /*@null@*/ void       *instStore,      /* instance storage */
            void       *aimStruc,       /* the AIM context */
            int        index,           /* the output index [1-nOut] */
            char       **aoname,        /* pointer to the returned name */
            capsValue  *formVal)        /* pointer to form/units (filled) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  return cntxt.aimOutput[i](instStore, aimStruc, index, aoname, formVal);
}


int
aim_PostAnalysis(aimContext cntxt,
                 const char *analysisName,
      /*@null@*/ void       *instStore,  /* instance storage */
                 void       *aimStruc,   /* the AIM context */
                 int        restart,     /* 0 - normal, 1 - restart */
      /*@null@*/ capsValue  *inputs)     /* complete suite of analysis inputs */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  return cntxt.aimPost[i](instStore, aimStruc, restart, inputs);
}


int
aim_CalcOutput(aimContext cntxt,
               const char *analysisName,
    /*@null@*/ void       *instStore,   /* instance storage */
               void       *aimStruc,    /* the AIM context */
               int        index,        /* the output index [1-nOut] */
               capsValue  *value)       /* pointer to value struct to fill */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return CAPS_NOTFOUND;
  
  return cntxt.aimCalc[i](instStore, aimStruc, index, value);
}


int
aim_Transfer(aimContext cntxt,
             const char *analysisName,
             capsDiscr  *discr,         /* the input discrete structure */
             char       *name,          /* the field name */
             int        npts,           /* number of points to be filled */
             int        rank,           /* the rank of the data */
             double     *data,          /* pointer to the memory to be filled */
             char       **units)        /* the units string returned */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimXfer[i] == NULL) {
    printf("aimTransfer not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimXfer[i](discr, name, npts, rank, data, units);
}


int
aim_Interpolation(aimContext cntxt,
                  const char *analysisName,
                  capsDiscr  *discr,    /* the input discrete structure */
                  const char *name,     /* the dataset name */
                  int        bIndex,   /* the input discretized body (1-bias) */
                  int        eIndex,    /* the input element (1-bias) */
                  double     *bary,     /* the barycentric coordinates */
                  int        rank,      /* the rank of the data */
                  double     *data,     /* global discrete support */
                  double     *result)   /* the result (rank in length) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                 == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimIntrp[i] == NULL) {
    printf("aimInterpolation not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntrp[i](discr, name, bIndex, eIndex, bary, rank, data, result);
}


int
aim_InterpolIndex(aimContext cntxt,
                  int        i,
                  capsDiscr  *discr,    /* the input discrete structure */
                  const char *name,     /* the dataset name */
                  int        bIndex,    /* the input discretized body (1-bias) */
                  int        eIndex,    /* the input element (1-bias) */
                  double     *bary,     /* the barycentric coordinates */
                  int        rank,      /* the rank of the data */
                  double     *data,     /* global discrete support */
                  double     *result)   /* the result (rank in length) */
{
  if ((i < 0) || (i >= cntxt.aim_nAnal)) return EGADS_RANGERR;
  if (cntxt.aimIntrp[i] == NULL) {
    printf("aimInterpolation not implemented in AIM %s\n", cntxt.aimName[i]);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntrp[i](discr, name, bIndex, eIndex, bary, rank, data, result);
}


int
aim_InterpolateBar(aimContext cntxt,
                   const char *analysisName,
                   capsDiscr  *discr,   /* the input discrete structure */
                   const char *name,    /* the dataset name */
                   int        bIndex,   /* the input discretized body (1-bias) */
                   int        eIndex,   /* the input element (1-bias) */
                   double     *bary,    /* the barycentric coordinates */
                   int        rank,     /* the rank of the data */
                   double     *r_bar,   /* input d(objective)/d(result) */
                   double     *d_bar)   /* returned d(objective)/d(data) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                    == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimIntrpBar[i] == NULL) {
    printf("aimInterpolateBar not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntrpBar[i](discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int
aim_InterpolIndBar(aimContext cntxt,
                   int        i,
                   capsDiscr  *discr,   /* the input discrete structure */
                   const char *name,    /* the dataset name */
                   int        bIndex,   /* the input discretized body (1-bias) */
                   int        eIndex,   /* the input element (1-bias) */
                   double     *bary,    /* the barycentric coordinates */
                   int        rank,     /* the rank of the data */
                   double     *r_bar,   /* input d(objective)/d(result) */
                   double     *d_bar)   /* returned d(objective)/d(data) */
{
  if ((i < 0) || (i >= cntxt.aim_nAnal)) return EGADS_RANGERR;
  if (cntxt.aimIntrpBar[i] == NULL) {
    printf("aimInterpolateBar not implemented in AIM %s\n", cntxt.aimName[i]);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntrpBar[i](discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int
aim_Integration(aimContext cntxt,
                const char *analysisName,
                capsDiscr  *discr,      /* the input discrete structure */
                const char  *name,      /* the dataset name */
                int        bIndex,      /* the input discretized body (1-bias) */
                int        eIndex,      /* the input element (1-bias) */
                int        rank,        /* the rank of the data */
     /*@null@*/ double     *data,       /* global discrete support */
                double     *result)     /* the result (rank in length) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                 == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimIntgr[i] == NULL) {
    printf("aimIntegration not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntgr[i](discr, name, bIndex, eIndex, rank, data, result);
}


int
aim_IntegrIndex(aimContext cntxt,
                int        i,
                capsDiscr  *discr,      /* the input discrete structure */
                const char *name,       /* the dataset name */
                int        bIndex,      /* the input discretized body (1-bias) */
                int        eIndex,      /* the input element (1-bias) */
                int        rank,        /* the rank of the data */
     /*@null@*/ double     *data,       /* global discrete support */
                double     *result)     /* the result (rank in length) */
{
  if ((i < 0) || (i >= cntxt.aim_nAnal)) return EGADS_RANGERR;
  if (cntxt.aimIntgr[i] == NULL) {
    printf("aimIntegration not implemented in AIM %s\n", cntxt.aimName[i]);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntgr[i](discr, name, bIndex, eIndex, rank, data, result);
}


int
aim_IntegrateBar(aimContext cntxt,
                 const char *analysisName,
                 capsDiscr  *discr,     /* the input discrete structure */
                 const char *name,      /* the dataset name */
                 int        bIndex,     /* the input discretized body (1-bias) */
                 int        eIndex,     /* the input element (1-bias) */
                 int        rank,       /* the rank of the data */
                 double     *r_bar,     /* input d(objective)/d(result) */
                 double     *d_bar)     /* returned d(objective)/d(data) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                    == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimIntgrBar[i] == NULL) {
    printf("aimIntegrateBar not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntgrBar[i](discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}


int
aim_IntegrIndBar(aimContext cntxt,
                 int        i,
                 capsDiscr  *discr,     /* the input discrete structure */
                 const char *name,      /* the dataset name */
                 int        bIndex,     /* the input discretized body (1-bias) */
                 int        eIndex,     /* the input element (1-bias) */
                 int        rank,       /* the rank of the data */
                 double     *r_bar,     /* input d(objective)/d(result) */
                 double     *d_bar)     /* returned d(objective)/d(data) */
{
  if ((i < 0) || (i >= cntxt.aim_nAnal)) return EGADS_RANGERR;
  if (cntxt.aimIntgrBar[i] == NULL) {
    printf("aimIntegrateBar not implemented in AIM %s\n", cntxt.aimName[i]);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimIntgrBar[i](discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}


int
aim_Backdoor(aimContext cntxt,
             const char *analysisName,
  /*@null@*/ void       *instStore,     /* instance storage */
             void       *aimStruc,      /* the AIM context */
             const char *JSONin,        /* the input(s) */
             char       **JSONout)      /* the output(s) */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i                 == -1)   return CAPS_NOTFOUND;
  if (cntxt.aimBdoor[i] == NULL) {
    printf("aimBackdoor not implemented in AIM %s\n", analysisName);
    return CAPS_NOTIMPLEMENT;
  }
  
  return cntxt.aimBdoor[i](instStore, aimStruc, JSONin, JSONout);
}


void
aim_cleanup(aimContext cntxt,
            const char *analysisName,
            void       *instStore)      /* instance data */
{
  int i;
  
  i = aimDLoaded(cntxt, analysisName);
  if (i == -1) return;
  cntxt.aimClean[i](instStore);
}


void aim_cleanupAll(aimContext *cntxt)
{
  int i;
  
  if (cntxt->aim_nAnal == 0) return;

  for (i = 0; i < cntxt->aim_nAnal; i++) {
    free(cntxt->aimName[i]);
    aimDLclose(cntxt->aimDLL[i]);
  }
  
  cntxt->aim_nAnal = 0;
}
