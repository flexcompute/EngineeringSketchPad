/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Utility Functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef WIN32
/* needs Advapi32.lib & Ws2_32.lib */
#include <Windows.h>
#include <direct.h>
#define strcasecmp stricmp
#define snprintf   _snprintf
#define access     _access
#define F_OK       0
#define SLASH     '\\'
#else
#include <unistd.h>
#include <dirent.h>
#include <dlfcn.h>
//#include <limits.h>
#include <sys/stat.h>
#define SLASH     '/'
#endif


#include "capsTypes.h"
#include "udunits.h"
#include "capsAIM.h"
#include "aimUtil.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"


#define EBUFSIZE        4096
#define UNIT_BUFFER_MAX  257

/*@-incondefs@*/
extern void ut_free(/*@only@*/ ut_unit* const unit);
extern void cv_free(/*@only@*/ cv_converter* const converter);
/*@+incondefs@*/


#ifdef WIN32
static int aim_flipSlash(const char *src, char *back)
{
  int len, i;

  len = strlen(src);
  if (len >= PATH_MAX) return EGADS_INDEXERR;

  for (i = 0; i <= len; i++)
    if (src[i] == '/') {
      back[i] = '\\';
    } else {
      back[i] = src[i];
    }

  return EGADS_SUCCESS;
}
#endif


void
aim_capsRev(int *major, int *minor)
{
  *major = CAPSMAJOR;
  *minor = CAPSMINOR;
}


int
aim_getRootPath(void *aimStruc, const char **fullPath)
{
  aimInfo     *aInfo;
  capsProblem *problem;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem   = aInfo->problem;
  *fullPath = problem->root;

  return CAPS_SUCCESS;
}


int
aim_file(void *aimStruc, const char *file, char *aimFile)
{
  int          status;
  aimInfo      *aInfo;
  capsAnalysis *analysis;
  const char   *filename = file;
#ifdef WIN32
  char         back[PATH_MAX];

  status = aim_flipSlash(file, back);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimStruc, "File path exceeds max length!");
    return CAPS_DIRERR;
  }
  filename = back;
#endif

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  analysis = (capsAnalysis *) aInfo->analysis;

  status = snprintf(aimFile, PATH_MAX, "%s%c%s",
                    analysis->fullPath, SLASH, filename);
  if (status >= PATH_MAX) {
    AIM_ERROR(aimStruc, "File path exceeds max length!");
    return CAPS_DIRERR;
  }

  return CAPS_SUCCESS;
}


static int
aim_fileSP(void *aimStruc, const char *file, char *aimFile)
{
  int          status;
  aimInfo      *aInfo;
  capsAnalysis *analysis;
  const char   *filename = file;
#ifdef WIN32
  char         back[PATH_MAX];

  status = aim_flipSlash(file, back);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimStruc, "File path exceeds max length!");
    return CAPS_DIRERR;
  }
  filename = back;
#endif

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  analysis = (capsAnalysis *) aInfo->analysis;

  status = snprintf(aimFile, PATH_MAX, "\"%s%c\"%s",
                    analysis->fullPath, SLASH, filename);
  if (status >= PATH_MAX) {
    AIM_ERROR(aimStruc, "File path exceeds max length!");
    return CAPS_DIRERR;
  }

  return CAPS_SUCCESS;
}


int
aim_isDir(void *aimStruc, const char *path)
{
  int   status;
  char  aimDir[PATH_MAX];
#ifdef WIN32
  DWORD dwAttr;
#else
  DIR   *dr = NULL;
#endif

  status = aim_file(aimStruc, path, aimDir);
  if (status != CAPS_SUCCESS) return status;

  /* are we a directory */
#ifdef WIN32
  dwAttr = GetFileAttributes(aimDir);
  if ((dwAttr != (DWORD) -1) && ((dwAttr & FILE_ATTRIBUTE_DIRECTORY) != 0))
    return CAPS_SUCCESS;
#else
/*@-mustfreefresh@*/
  dr = opendir(aimDir);
  if (dr != NULL) {
    closedir(dr);
    return CAPS_SUCCESS;
  }
/*@+mustfreefresh@*/
#endif

  return CAPS_NOTFOUND;
}


int
aim_mkDir(void *aimStruc, const char *path)
{
  int  status;
  char aimDir[PATH_MAX];

  status = aim_file(aimStruc, path, aimDir);
  if (status != CAPS_SUCCESS) return status;

#ifdef WIN32
  status = _mkdir(aimDir);
#else
  status =  mkdir(aimDir, S_IRWXU);
#endif

  if (status != 0) {
/*@-unrecog@*/
    if (errno != EEXIST) {
/*@+unrecog@*/
      AIM_ERROR(aimStruc, "Unable to make: %s", aimDir);
      return CAPS_DIRERR;
    }
  }

  return CAPS_SUCCESS;
}


int
aim_rmDir(void *aimStruc, const char *path)
{
  int  status;
  size_t i, len;
  char cmd[PATH_MAX+14];
  char aimDir[PATH_MAX];
  int  wild = (int) false;

  len = strlen(path);
  for (i = 0; i < len; i++) {
    if (path[i] == '*' || path[i] == '?') wild = (int) true;
  }
  if (wild == (int) false) {
    status = aim_isDir(aimStruc, path);
    if (status != CAPS_SUCCESS) return status;
  }
  status = aim_fileSP(aimStruc, path, aimDir);
  if (status != CAPS_SUCCESS) return status;

#ifdef WIN32
  snprintf(cmd, PATH_MAX+14, "rmdir /Q /S %s", aimDir);
  fflush(NULL);
#else
  if (status != CAPS_SUCCESS) return status;
  snprintf(cmd, PATH_MAX+14, "rm -rf %s", aimDir);
#endif
  status = system(cmd);
  if ((status == -1) || (status == 127)) return CAPS_DIRERR;

  return CAPS_SUCCESS;
}


int
aim_isFile(void *aimStruc, const char *file)
{
  int  status;
  char aimFile[PATH_MAX];

  status = aim_file(aimStruc, file, aimFile);
  if (status != CAPS_SUCCESS) return status;

  if (access(aimFile, F_OK) == 0) return CAPS_SUCCESS;

  return CAPS_NOTFOUND;
}


int
aim_rmFile(void *aimStruc, const char *file)
{
  int  status;
  char aimFile[PATH_MAX];

  status = aim_file(aimStruc, file, aimFile);
  if (status != CAPS_SUCCESS) return status;

  remove(aimFile);

  return CAPS_SUCCESS;
}


int aim_cpFile(void *aimStruc, const char *src, const char *dst)
{
  int  status;
  char cmd[2*PATH_MAX+14], aimDst[PATH_MAX];
#ifdef WIN32
  char sback[PATH_MAX], dback[PATH_MAX];
#endif

  if (strlen(src) > PATH_MAX) {
    AIM_ERROR(aimStruc, "File src path exceeds max length!");
    return CAPS_IOERR;
  }
  if (strlen(dst) > PATH_MAX) {
    AIM_ERROR(aimStruc, "File dst path exceeds max length!");
    return CAPS_IOERR;
  }

  status = aim_file(aimStruc, dst, aimDst);
  if (status != CAPS_SUCCESS) return status;

#ifdef WIN32
  status = aim_flipSlash(src, sback);
  if (status != EGADS_SUCCESS)  return status;
  status = aim_flipSlash(aimDst, dback);
  if (status != EGADS_SUCCESS)  return status;
  snprintf(cmd, 2*PATH_MAX+14, "copy /Y \"%s\" \"%s\"", sback, dback);
  fflush(NULL);
#else
  snprintf(cmd, 2*PATH_MAX+14, "cp '%s' '%s'", src, aimDst);
#endif
  status = system(cmd);
  if (status != 0) {
    AIM_ERROR(aimStruc, "Could not execute: %s", cmd);
    return CAPS_IOERR;
  }

  return CAPS_SUCCESS;
}


int aim_relPath(void *aimStruc, const char *src,
                /*@null@*/ const char *dst, char *relPath)
{
  int         i, j, len, status;
  char        aimDst[PATH_MAX], *ptr;
  aimInfo     *aInfo;
  capsProblem *problem;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem = aInfo->problem;

  if (strlen(src) > PATH_MAX) {
    AIM_ERROR(aimStruc, "File src path exceeds max length!");
    return CAPS_IOERR;
  }
  if (relPath == NULL) {
    AIM_ERROR(aimStruc, "NULL relPath!");
    return CAPS_IOERR;
  }
  relPath[0] = '\0';
  relPath[1] = '\0';
  relPath[2] = '\0';

  /* get the destination file */
  if (dst != NULL) {
    if (strlen(dst) > PATH_MAX) {
      AIM_ERROR(aimStruc, "File dst path exceeds max length!");
      return CAPS_IOERR;
    }
    if (strlen(dst) > 0) {
      status = aim_file(aimStruc, dst, aimDst);
    } else {
      status = aim_file(aimStruc, ".", aimDst);
    }
  } else {
    status   = aim_file(aimStruc, ".", aimDst);
  }
  AIM_STATUS(aimStruc, status);

  len = strlen(problem->root);
  if ((len >= strlen(src)) || (len >= strlen(aimDst))) {
    AIM_ERROR(aimStruc, "File not in rootPath!");
    return CAPS_IOERR;
  }
  for (i = 0; i < len; i++)
    if ((problem->root[i] != src[i]) || (problem->root[i] != aimDst[i])) {
      AIM_ERROR(aimStruc, "Path mismatch!");
      return CAPS_IOERR;
    }
  len++; /* skip the slash */

  j = 0;
  ptr = aimDst+len;
  while((ptr = strchr(ptr, SLASH)) != NULL) {
    relPath[j++] = '.';
    relPath[j++] = '.';
    relPath[j++] = SLASH;
    ptr++;
  }

  for (i = len; i <= strlen(src); i++, j++) relPath[j] = src[i];

cleanup:
  return status;
}


int aim_symLink(void *aimStruc, const char *src, /*@null@*/ const char *dst)
{
#ifdef WIN32
  if (dst == NULL) return aim_cpFile(aimStruc, src, ".");
  return      aim_cpFile(aimStruc, src, dst);
#else
  int         i, j, k, m, len, status;
  char        cmd[2*PATH_MAX+14], aimDst[PATH_MAX], relSrc[PATH_MAX];
  aimInfo     *aInfo;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;

  if (strlen(src) > PATH_MAX) {
    AIM_ERROR(aimStruc, "File src path exceeds max length!");
    return CAPS_IOERR;
  }
  if (access(src, F_OK) != 0) {
    AIM_ERROR(aimStruc, "%s Not a File!", src);
    return CAPS_IOERR;
  }

  /* convert the absolute src path to a relative path */
  status = aim_relPath(aimStruc, src, dst, relSrc);
  AIM_STATUS(aimStruc, status);

  /* get the destination file */
  if (dst != NULL) {
    if (strlen(dst) > PATH_MAX) {
      AIM_ERROR(aimStruc, "File dst path exceeds max length!");
      return CAPS_IOERR;
    }
    if (strlen(dst) > 0) {
      status = aim_file(aimStruc, dst, aimDst);
    } else {
      status = aim_file(aimStruc, ".", aimDst);
    }
  } else {
    status   = aim_file(aimStruc, ".", aimDst);
  }
  AIM_STATUS(aimStruc, status);

  /* remove any old links */
  len = strlen(aimDst);
  if ((aimDst[len-1] == '.') && (aimDst[len-2] == '/')) {
    for (i = strlen(relSrc); i > 0; i--)
      if (relSrc[i] == '/') break;
    m = len-2;
    if ((i == 0) && (relSrc[0] != '/')) m = len-1;
    j = strlen(relSrc);
    for (k = i; k < j; k++) aimDst[m+k-i] = relSrc[k];
    aimDst[m+j-i] = 0;
  }
  unlink(aimDst);

  snprintf(cmd, 2*PATH_MAX+14, "ln -s '%s' '%s'", relSrc, aimDst);
  status = system(cmd);
  if (status != 0) {
    AIM_ERROR(aimStruc, "Could not execute: %s", cmd);
    return CAPS_IOERR;
  }

cleanup:
  return status;
#endif
}


/*@null@*/ /*@out@*/ /*@only@*/ FILE *
aim_fopen(void *aimStruc, const char *path, const char *mode)
{
  int          status;
  char         fullPath[PATH_MAX];
  aimInfo      *aInfo;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return NULL;
  if (aInfo->magicnumber != CAPSMAGIC) return NULL;
  analysis = (capsAnalysis *) aInfo->analysis;

  status   = snprintf(fullPath, PATH_MAX, "%s%c%s",
                      analysis->fullPath, SLASH, path);
  if (status >= PATH_MAX) return NULL;

/*@-dependenttrans@*/
  return fopen(fullPath, mode);
/*@+dependenttrans@*/
}


int
aim_system(void *aimStruc, /*@null@*/ const char *rpath, const char *command)
{
  size_t       status;
  char         *fullcommand;
  size_t       len;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  len = 9 + strlen(problem->root) + 1 + strlen(analysis->path) +
        4 + strlen(command) + 3;
  if (rpath == NULL) {
    fullcommand = EG_alloc(len*sizeof(char));
    if (fullcommand == NULL) return EGADS_MALLOC;
#ifdef WIN32
    status = snprintf(fullcommand, len, "%c: && cd \"%s\\%s\" && %s",
                      problem->root[0], &problem->root[2], analysis->path,
                      command);
#else
    status = snprintf(fullcommand, len, "cd '%s/%s' && %s",
                      problem->root, analysis->path, command);
#endif
  } else {
    len += strlen(rpath) + 1;
    fullcommand = EG_alloc(len*sizeof(char));
    if (fullcommand == NULL) return EGADS_MALLOC;
#ifdef WIN32
    status = snprintf(fullcommand, len, "%c: && cd \"%s\\%s\\%s\" && %s",
                      problem->root[0], &problem->root[2], analysis->path,
                      rpath, command);
#else
    status = snprintf(fullcommand, len, "cd '%s/%s/%s' && %s",
                      problem->root, analysis->path, rpath, command);
#endif
  }
  if (status >= len) {
    EG_free(fullcommand);
    return CAPS_BADVALUE;
  }

  status = system(fullcommand);
  EG_free(fullcommand);
  return status == 0 ? CAPS_SUCCESS : CAPS_EXECERR;
}


int
aim_getUnitSys(void *aimStruc, char **unitSys)
{
  aimInfo      *aInfo;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  analysis = (capsAnalysis *) aInfo->analysis;

  *unitSys = analysis->unitSys;
  return CAPS_SUCCESS;
}


int
aim_getBodies(void *aimStruc, const char **intents, int *nBody, ego **bodies)
{
  aimInfo      *aInfo;
  capsAnalysis *analysis;

  *nBody  = 0;
  *bodies = NULL;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  analysis = (capsAnalysis *) aInfo->analysis;

  *intents = analysis->intents;
  *nBody   = analysis->nBody;
  *bodies  = analysis->bodies;
  return CAPS_SUCCESS;
}


int
aim_newGeometry(void *aimStruc)
{
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                               return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)             return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  if (problem->geometry.sNum < analysis->pre.sNum) return CAPS_CLEAN;
  return CAPS_SUCCESS;
}


int
aim_newAnalysisIn(void *aimStruc, int index)
{
  aimInfo      *aInfo;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                               return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)             return CAPS_BADOBJECT;
  analysis = (capsAnalysis *) aInfo->analysis;
  if (index <= 0 || index > analysis->nAnalysisIn) return CAPS_RANGEERR;

  if (analysis->analysisIn[index-1]->last.sNum < analysis->pre.sNum)
    return CAPS_CLEAN;
  return CAPS_SUCCESS;
}


int
aim_numInstance(void *aimStruc)
{
  int          i;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  for (i = 0; i < problem->aimFPTR.aim_nAnal; i++)
    if (strcasecmp(analysis->loadName, problem->aimFPTR.aimName[i]) == 0) break;
  if (i == problem->aimFPTR.aim_nAnal) return CAPS_BADINIT;

  return problem->aimFPTR.aim_nInst[i];
}


int
aim_getInstance(void *aimStruc)
{
  aimInfo *aInfo;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;

  return aInfo->instance;
}


int
aim_convert(void *aimStruc, const int count,
            /*@null@*/ const char  *inUnits, double  *inValue,
            /*@null@*/ const char *outUnits, double *outValue)
{
  int          i;
  aimInfo      *aInfo;
  capsProblem  *problem;
  ut_unit      *utunit1, *utunit2;
  cv_converter *converter;

  if (inValue  == NULL) return CAPS_NULLVALUE;
  if (outValue == NULL) return CAPS_NULLVALUE;

  /* check for simple equality */
  if ((inUnits == NULL) && (outUnits == NULL)) {
    for (i = 0; i < count; i++) outValue[i] = inValue[i];
    return CAPS_SUCCESS;
  }
  if (inUnits  == NULL) return CAPS_UNITERR;
  if (outUnits == NULL) return CAPS_UNITERR;
  if (strcmp(outUnits, inUnits) == 0) {
    for (i = 0; i < count; i++) outValue[i] = inValue[i];
    return CAPS_SUCCESS;
  }

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem = aInfo->problem;

  utunit1 = ut_parse((ut_system *) problem->utsystem,  inUnits, UT_ASCII);
  utunit2 = ut_parse((ut_system *) problem->utsystem, outUnits, UT_ASCII);
  converter = ut_get_converter(utunit1, utunit2);
  if (converter == NULL) {
    ut_free(utunit1);
    ut_free(utunit2);
    AIM_ERROR(aimStruc, "Cannot convert units '%s' to '%s'", inUnits, outUnits);
    return CAPS_UNITERR;
  }

  (void)cv_convert_doubles(converter, inValue, count, outValue);
  cv_free(converter);
  ut_free(utunit1);
  ut_free(utunit2);

  if (ut_get_status() != UT_SUCCESS) {
    AIM_ERROR(aimStruc, "Cannot convert units '%s' to '%s'", inUnits, outUnits);
    return CAPS_UNITERR;
  }

  return CAPS_SUCCESS;
}


int
aim_unitMultiply(void *aimStruc, const char *inUnits1, const char *inUnits2,
                 char **outUnits)
{
  int         status;
  aimInfo     *aInfo;
  capsProblem *problem;
  ut_unit     *utunit1, *utunit2, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnits1 == NULL) ||
      (inUnits2 == NULL) || (outUnits == NULL)) return CAPS_NULLNAME;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                            return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)          return CAPS_BADOBJECT;
  problem = aInfo->problem;

  *outUnits = NULL;

  utunit1 = ut_parse((ut_system *) problem->utsystem, inUnits1, UT_ASCII);
  utunit2 = ut_parse((ut_system *) problem->utsystem, inUnits2, UT_ASCII);
  utunit  = ut_multiply(utunit1, utunit2);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit2);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit2);
  ut_free(utunit);

  if (status < UT_SUCCESS) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
aim_unitDivide(void *aimStruc, const char *inUnits1, const char *inUnits2,
               char **outUnits)
{
  int         status;
  aimInfo     *aInfo;
  capsProblem *problem;
  ut_unit     *utunit1, *utunit2, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnits1 == NULL) ||
      (inUnits2 == NULL) || (outUnits == NULL)) return CAPS_NULLNAME;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                            return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)          return CAPS_BADOBJECT;
  problem = aInfo->problem;

  *outUnits = NULL;

  utunit1 = ut_parse((ut_system *) problem->utsystem, inUnits1, UT_ASCII);
  utunit2 = ut_parse((ut_system *) problem->utsystem, inUnits2, UT_ASCII);
  utunit  = ut_divide(utunit1, utunit2);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit2);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit2);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
aim_unitInvert(void *aimStruc, const char *inUnit, char **outUnits)
{
  int         status;
  aimInfo     *aInfo;
  capsProblem *problem;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLNAME;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                          return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)        return CAPS_BADOBJECT;
  problem = aInfo->problem;

  *outUnits = NULL;

  utunit1 = ut_parse((ut_system *) problem->utsystem, inUnit, UT_ASCII);
  utunit  = ut_invert(utunit1);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
aim_unitRaise(void *aimStruc, const char *inUnit, const int power,
              char **outUnits)
{
  int         status;
  aimInfo     *aInfo;
  capsProblem *problem;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLNAME;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                          return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)        return CAPS_BADOBJECT;
  problem = aInfo->problem;

  *outUnits = NULL;

  utunit1 = ut_parse((ut_system *) problem->utsystem, inUnit, UT_ASCII);
  utunit  = ut_raise(utunit1, power);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
aim_unitOffset(void *aimStruc, const char *inUnit, const double offset,
               char **outUnits)
{
  int         status;
  aimInfo     *aInfo;
  capsProblem *problem;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLNAME;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                          return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)        return CAPS_BADOBJECT;
  problem = aInfo->problem;

  *outUnits = NULL;

  utunit1 = ut_parse((ut_system *) problem->utsystem, inUnit, UT_ASCII);
  utunit  = ut_offset(utunit1, offset);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
aim_getIndex(void *aimStruc, /*@null@*/ const char *name,
             enum capssType subtype)
{
  int          i, nobj;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsObject   **objs;

  aInfo = (aimInfo *) aimStruc;
  if ((subtype != GEOMETRYIN) && (subtype != GEOMETRYOUT) &&
      (subtype != ANALYSISIN) && (subtype != ANALYSISOUT) &&
      (subtype != ANALYSISDYNO))       return CAPS_BADTYPE;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;
  if (subtype == GEOMETRYIN) {
    nobj = problem->nGeomIn;
    objs = problem->geomIn;
  } else if (subtype == GEOMETRYOUT) {
    nobj = problem->nGeomOut;
    objs = problem->geomOut;
  } else if (subtype == ANALYSISIN)  {
    nobj = analysis->nAnalysisIn;
    objs = analysis->analysisIn;
  } else if (subtype == ANALYSISOUT) {
    nobj = analysis->nAnalysisOut;
    objs = analysis->analysisOut;
  } else {
    nobj = analysis->nAnalysisDynO;
    objs = analysis->analysisDynO;
  }
  if (name == NULL) return nobj;

  for (i = 0; i < nobj; i++) {
    if (objs[i]->name == NULL) {
      printf(" aimUtil Info: object %d with NULL name!\n", i+1);
      continue;
    }
    if (strcmp(objs[i]->name, name) == 0) return i+1;
  }

  return CAPS_NOTFOUND;
}


int
aim_getValue(void *aimStruc, int index, enum capssType subtype,
             capsValue **value)
{
  int          nobj;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsObject   **objs;
  const char *objName = "";

  aInfo = (aimInfo *) aimStruc;
  if ((subtype != GEOMETRYIN) && (subtype != GEOMETRYOUT) &&
      (subtype != ANALYSISIN) && (subtype != ANALYSISOUT) &&
      (subtype != ANALYSISDYNO))       return CAPS_BADTYPE;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (index <= 0)                      return CAPS_BADINDEX;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;
  if (subtype == GEOMETRYIN) {
    nobj = problem->nGeomIn;
    objs = problem->geomIn;
    objName = "GEOMETRYIN";
  } else if (subtype == GEOMETRYOUT) {
    nobj = problem->nGeomOut;
    objs = problem->geomOut;
    objName = "GEOMETRYOUT";
  } else if (subtype == ANALYSISIN)  {
    nobj = analysis->nAnalysisIn;
    objs = analysis->analysisIn;
    objName = "ANALYSISIN";
  } else if (subtype == ANALYSISOUT) {
    nobj = analysis->nAnalysisOut;
    objs = analysis->analysisOut;
    objName = "ANALYSISOUT";
  } else {
    nobj = analysis->nAnalysisDynO;
    objs = analysis->analysisDynO;
    objName = "ANALYSISDYNO";
  }
  if (index > nobj) {
    AIM_ERROR(aimStruc, "%s Index (%d > %d) out-of-range!", objName, index, nobj);
    return CAPS_BADINDEX;
  }

  *value = (capsValue *) objs[index-1]->blind;
  return CAPS_SUCCESS;
}


int
aim_initValue(capsValue *value)
{
  if (value == NULL) return CAPS_NULLVALUE;

  value->length          = value->nrow = value->ncol = 1;
  value->type            = Integer;
  value->dim             = value->pIndex = value->index = 0;
  value->lfixed          = value->sfixed = Fixed;
  value->nullVal         = NotAllowed;
  value->units           = NULL;
  value->meshWriter      = NULL;
  value->link            = NULL;
  value->vals.reals      = NULL;
  value->limits.dlims[0] = value->limits.dlims[1] = 0.0;
  value->linkMethod      = Copy;
  value->gInType         = 0;
  value->partial         = NULL;
  value->nderiv          = 0;
  value->derivs          = NULL;

  return CAPS_SUCCESS;
}


int
aim_copyValue(capsValue *value, capsValue *copy)
{
  int i, j, len;

  if (value == NULL) return CAPS_NULLVALUE;
  if (copy  == NULL) return CAPS_NULLVALUE;

  aim_initValue(copy);

  copy->length          = value->length;
  copy->nrow            = value->nrow;
  copy->ncol            = value->ncol;
  copy->type            = value->type;
  copy->dim             = value->dim;
  copy->pIndex          = value->pIndex;
  copy->index           = value->index;
  copy->lfixed          = value->lfixed;
  copy->sfixed          = value->sfixed;
  copy->nullVal         = value->nullVal;
  copy->units           = EG_strdup(value->units);
  copy->meshWriter      = value->meshWriter;
  copy->link            = value->link;

  if (copy->type == Double || copy->type == DoubleDeriv) {
    if (value->length > 1) {
      copy->vals.reals = (double*)EG_alloc(value->length*sizeof(double));
      if (copy->vals.reals == NULL) return EGADS_MALLOC;
      for (i = 0; i < value->length; i++)
        copy->vals.reals[i] = value->vals.reals[i];
    } else {
      copy->vals.real = value->vals.real;
    }
  } else if (copy->type == Integer) {
    if (value->length > 1) {
      copy->vals.integers = (int*)EG_alloc(value->length*sizeof(int));
      if (copy->vals.integers == NULL) return EGADS_MALLOC;
      for (i = 0; i < value->length; i++)
        copy->vals.integers[i] = value->vals.integers[i];
    } else {
      copy->vals.integer = value->vals.integer;
    }
  } else if (copy->type == String) {
    if (value->length > 1) {
      len = 0;
      for (i = 0; i < value->length; i++)
        len += strlen(&value->vals.string[len])+1;

      copy->vals.string = (char*)EG_alloc(len*sizeof(char));
      if (copy->vals.string == NULL) return EGADS_MALLOC;
      for (i = 0; i < len; i++)
        copy->vals.string[i] = value->vals.string[i];
    } else {
      copy->vals.string = EG_strdup(value->vals.string);
    }
  } else if (copy->type == Tuple) {

    copy->vals.tuple = (capsTuple*)EG_alloc(value->length*sizeof(capsTuple));
    if (copy->vals.tuple == NULL) return EGADS_MALLOC;
    for (i = 0; i < value->length; i++) {
      copy->vals.tuple[i].name = EG_strdup(value->vals.tuple[i].name);
      copy->vals.tuple[i].value = EG_strdup(value->vals.tuple[i].value);
    }
  } else {
    return CAPS_NOTIMPLEMENT;
  }

  copy->limits.dlims[0] = value->limits.dlims[0];
  copy->limits.dlims[1] = value->limits.dlims[1];
  copy->linkMethod      = value->linkMethod;
  copy->gInType         = value->gInType;
  if (value->partial != NULL) {
    copy->partial = (int*)EG_alloc(value->length*sizeof(int));
    if (copy->partial == NULL) return EGADS_MALLOC;

    for (i = 0; i < value->length; i++)
      copy->partial[i] = value->partial[i];
  }
  copy->nderiv = value->nderiv;
  if (value->derivs != NULL) {

    copy->derivs = (capsDeriv*)EG_alloc(value->nderiv*sizeof(capsDeriv));
    if (copy->derivs == NULL) return EGADS_MALLOC;

    for (i = 0; i < value->length; i++) {
      copy->derivs[i].len_wrt = value->derivs[i].len_wrt;
      copy->derivs[i].name = EG_strdup(copy->derivs[i].name);

      len = value->length*value->derivs[i].len_wrt;
      copy->derivs[i].deriv = (double*)EG_alloc(len*sizeof(double));
      for (j = 0; j < len; j++) {
        copy->derivs[i].deriv[j] = value->derivs[i].deriv[j];
      }
    }
  }

  return CAPS_SUCCESS;
}


void
aim_freeValue(capsValue *value)
{
  int i;

  if (value == NULL) return;

  EG_free(value->units);

  if (value->length > 1) {
    if (value->type == Double || value->type == DoubleDeriv) {
      EG_free(value->vals.reals);
    } else if (value->type == Integer) {
      EG_free(value->vals.integers);
    } else if (value->type == String) {
      EG_free(value->vals.string);
    } else if (value->type == Tuple) {
      for (i = 0; i < value->length; i++) {
        EG_free(value->vals.tuple[i].name);
        EG_free(value->vals.tuple[i].value);
      }
      EG_free(value->vals.tuple);
    }
  }
  EG_free(value->partial);
  if (value->derivs != NULL) {
    for (i = 0; i < value->nderiv; i++) {
      EG_free(value->derivs[i].name);
      EG_free(value->derivs[i].deriv);
    }
    EG_free(value->derivs);
  }

  aim_initValue(value);
}


int
aim_makeDynamicOutput(void *aimStruc, const char *dynObjName, capsValue *value)
{
  int          i;
  aimInfo      *aInfo;
  capsAnalysis *analysis;
  capsProblem  *problem;
  capsValue    *val;
  capsObject   *obj, **tmp, *pobject;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (dynObjName == NULL)              return CAPS_NULLNAME;
  if (value      == NULL)              return CAPS_NULLVALUE;
  if (aInfo->inPost != 1)              return CAPS_STATEERR;
  problem  = aInfo->problem;
  pobject  = problem->mySelf;
  analysis = (capsAnalysis *) aInfo->analysis;

  for (i = 0; i < analysis->nAnalysisDynO; i++) {
    obj = analysis->analysisDynO[i];
    if (obj       == NULL) continue;
    if (obj->name == NULL) continue;
    if (strcmp(obj->name, dynObjName) == 0) return CAPS_BADNAME;
  }
  for (i = 0; i < problem->nAnalysis; i++) {
    if (problem->analysis[i] == NULL) continue;
    if (problem->analysis[i]->blind == analysis) break;
  }
  if (i == problem->nAnalysis) return CAPS_NOTFOUND;

  /* make the object and fill it in */

  val = (capsValue  *) EG_alloc(sizeof(capsValue));
  obj = (capsObject *) EG_alloc(sizeof(capsObject));
  if ((obj == NULL) || (val == NULL)) {
    if (obj != NULL) EG_free(obj);
    if (val != NULL) EG_free(val);
    return EGADS_MALLOC;
  }

  if (analysis->analysisDynO == NULL) {
    analysis->analysisDynO = (capsObject  **) EG_alloc(sizeof(capsObject *));
    if (analysis->analysisDynO == NULL) {
      EG_free(obj);
      EG_free(val);
      return EGADS_MALLOC;
    }
    analysis->nAnalysisDynO = 0;
  } else {
    tmp = (capsObject  **) EG_reall( analysis->analysisDynO,
                                    (analysis->nAnalysisDynO+1)*
                                    sizeof(capsObject *));
    if (tmp == NULL) {
      EG_free(obj);
      EG_free(val);
      return EGADS_MALLOC;
    }
    analysis->analysisDynO = tmp;
  }
  /* copy contents */
  (*val) = (*value);
  val->index = analysis->nAnalysisDynO+1;

  obj->magicnumber = CAPSMAGIC;
  obj->type        = VALUE;
  obj->subtype     = ANALYSISDYNO;
  obj->delMark     = 0;
  obj->name        = EG_strdup(dynObjName);
  obj->attrs       = NULL;
  obj->blind       = val;
  obj->flist       = NULL;
  obj->parent      = problem->analysis[i];
  obj->nHistory    = 0;
  obj->history     = NULL;
  obj->last.index  = -1;
  obj->last.pname  = EG_strdup(pobject->last.pname);
  obj->last.pID    = EG_strdup(pobject->last.pID);
  obj->last.user   = EG_strdup(pobject->last.user);
  obj->last.sNum   = pobject->last.sNum;
  for (i = 0; i < 6; i++)
    obj->last.datetime[i] = pobject->last.datetime[i];

  analysis->analysisDynO[analysis->nAnalysisDynO] = obj;
  analysis->nAnalysisDynO++;

  /* clear out the input value structure */
  (void) aim_initValue(value);
  return CAPS_SUCCESS;
}


int
aim_getName(void *aimStruc, int index, enum capssType subtype,
            const char **name)
{
  int          nobj;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsObject   **objs;

  aInfo = (aimInfo *) aimStruc;
  if ((subtype != GEOMETRYIN) && (subtype != GEOMETRYOUT) &&
      (subtype != ANALYSISIN) && (subtype != ANALYSISOUT) &&
      (subtype != ANALYSISDYNO))       return CAPS_BADTYPE;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (index <= 0)                      return CAPS_BADINDEX;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;
  if (subtype == GEOMETRYIN) {
    nobj = problem->nGeomIn;
    objs = problem->geomIn;
  } else if (subtype == GEOMETRYOUT) {
    nobj = problem->nGeomOut;
    objs = problem->geomOut;
  } else if (subtype == ANALYSISIN)  {
    nobj = analysis->nAnalysisIn;
    objs = analysis->analysisIn;
  } else if (subtype == ANALYSISOUT) {
    nobj = analysis->nAnalysisOut;
    objs = analysis->analysisOut;
  } else {
    nobj = analysis->nAnalysisDynO;
    objs = analysis->analysisDynO;
  }
  if (index > nobj)                    return CAPS_BADINDEX;

  *name = objs[index-1]->name;
  return CAPS_SUCCESS;
}


int
aim_getGeomInType(void *aimStruc, int index)
{
  aimInfo     *aInfo;
  capsObject  *vobj;
  capsValue   *value;
  capsProblem *problem;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (index <= 0)                      return CAPS_BADINDEX;
  problem = aInfo->problem;
  if (index > problem->nGeomIn)        return CAPS_BADINDEX;
  vobj = problem->geomIn[index-1];
  if (vobj->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (vobj->type        != VALUE)      return CAPS_BADTYPE;
  if (vobj->blind       == NULL)       return CAPS_NULLBLIND;
  value = (capsValue *) vobj->blind;

  return value->gInType;
}


int
aim_newTess(void *aimStruc, ego tess)
{
  int          stat, i, state, npts;
  aimInfo      *aInfo;
  capsAnalysis *analysis;
  ego          body, tbody;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (tess == NULL)                    return EGADS_NULLOBJ;
  if (tess->magicnumber  != MAGIC)     return EGADS_NOTOBJ;
  analysis = (capsAnalysis *) aInfo->analysis;

  /* allocate tessellation storage if it does not exist */
  if (analysis->tess == NULL) {
    analysis->tess = (ego *) EG_alloc(analysis->nBody*sizeof(ego));
    if (analysis->tess == NULL) return EGADS_MALLOC;
    analysis->nTess = analysis->nBody;
    for (i = 0; i < analysis->nBody; i++) analysis->tess[i] = NULL;
  }

  /* associate a tessellation object */
  stat = EG_statusTessBody(tess, &body, &state, &npts);
  if (stat <  EGADS_SUCCESS) return stat;
  if (stat == EGADS_OUTSIDE) return CAPS_SOURCEERR;

  for (i = 0; i < analysis->nBody; i++) {
    if (body != analysis->bodies[i]) continue;
    if (analysis->tess[i] != NULL) {
      EG_deleteObject(analysis->tess[i]);
    }
    analysis->tess[i] = tess;
    return CAPS_SUCCESS;
  }

  /* look for AIM created bodies */
  for (i = analysis->nBody-1; i < analysis->nTess; i++) {
    stat = EG_statusTessBody(analysis->tess[i], &tbody, &state, &npts);
    if (stat <  EGADS_SUCCESS) return stat;
    if (stat == EGADS_OUTSIDE) return CAPS_SOURCEERR;
    if (tbody != body) continue;
    if (analysis->tess[i] != NULL) {
      EG_deleteObject(analysis->tess[i]);
    }
    analysis->tess[i] = tess;
    return CAPS_SUCCESS;
  }

  /* not in the Body list -- extend the list of tessellations */
  analysis->tess = (ego *) EG_reall(analysis->tess, (analysis->nTess+1)*sizeof(ego));
  if (analysis->tess == NULL) return EGADS_MALLOC;
  analysis->tess[analysis->nTess] = tess;
  analysis->nTess++;

  return CAPS_SUCCESS;
}


int
aim_getDiscr(void *aimStruc, const char *bname, capsDiscr **discr)
{
  int           i, j;
  aimInfo       *aInfo;
  capsProblem   *problem;
  capsBound     *bound;
  capsAnalysis  *analysis, *analy;
  capsVertexSet *vertexset;

  *discr = NULL;
  aInfo  = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  /* find bound with same name */
  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]       == NULL) continue;
    if (problem->bounds[i]->name == NULL) continue;
    if (strcmp(bname, problem->bounds[i]->name) != 0) continue;
    /* find our analysis in the Bound */
    bound = (capsBound *) problem->bounds[i]->blind;
    if (bound == NULL) return CAPS_NULLOBJ;
    for (j = 0; j < bound->nVertexSet; j++) {
      if (bound->vertexSet[j]              == NULL)      continue;
      if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
      if (bound->vertexSet[j]->type        != VERTEXSET) continue;
      if (bound->vertexSet[j]->blind       == NULL)      continue;
      vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
      if (vertexset           == NULL) continue;
      if (vertexset->analysis == NULL) continue;
      analy = (capsAnalysis *) vertexset->analysis->blind;
      if (analy != analysis) continue;
      *discr = vertexset->discr;
      return CAPS_SUCCESS;
    }
  }

  return CAPS_NOTFOUND;
}


int
aim_getDiscrState(void *aimStruc, const char *bname)
{
  int         i;
  aimInfo     *aInfo;
  capsProblem *problem;
  capsObject  *bobject;
  capsBound   *bound;

  aInfo  = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem = aInfo->problem;

  /* find bound with same name */
  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]       == NULL) continue;
    if (problem->bounds[i]->name == NULL) continue;
    if (strcmp(bname, problem->bounds[i]->name) != 0) continue;
    bobject = problem->bounds[i];
    bound   = (capsBound *) bobject->blind;
    if (bound == NULL) return CAPS_NULLOBJ;
    if (bobject->last.sNum < problem->geometry.sNum) return CAPS_DIRTY;
    return CAPS_SUCCESS;
  }

  return CAPS_NOTFOUND;
}


int
aim_getDataSet(capsDiscr *discr, const char *dname, enum capsdMethod *method,
               int *npts, int *rank, double **data, char **units)
{
  int           i, j, k;
  aimInfo       *aInfo;
  capsObject    *link, *dobject;
  capsProblem   *problem;
  capsBound     *bound;
  capsAnalysis  *analysis, *analy;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;

  *method = Interpolate;
  *npts   = *rank = 0;
  *data   = NULL;
  if (discr == NULL)                   return CAPS_NULLOBJ;
  aInfo   = discr->aInfo;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  /* find dataset in discr */
  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i] == NULL) continue;
    bound = (capsBound *) problem->bounds[i]->blind;
    if (bound == NULL) continue;
    for (j = 0; j < bound->nVertexSet; j++) {
      if (bound->vertexSet[j]              == NULL)      continue;
      if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
      if (bound->vertexSet[j]->type        != VERTEXSET) continue;
      if (bound->vertexSet[j]->blind       == NULL)      continue;
      vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
      if (vertexset           == NULL) continue;
      if (vertexset->analysis == NULL) continue;
      analy = (capsAnalysis *) vertexset->analysis->blind;
      if (analy            != analysis) continue;
      if (vertexset->discr != discr)    continue;
      for (k = 0; k < vertexset->nDataSets; k++) {
        if (vertexset->dataSets[k]       == NULL) continue;
        if (vertexset->dataSets[k]->name == NULL) continue;
        if (strcmp(dname, vertexset->dataSets[k]->name) != 0) continue;
        dobject = vertexset->dataSets[k];
        dataset = (capsDataSet *) dobject->blind;
        if (dataset == NULL) return CAPS_NULLOBJ;
        /* are we dirty? */
        if (dataset->ftype == FieldIn) {
          /* check linked source DataSets */
          link = dataset->link;
          if (link == NULL                  ) return CAPS_SOURCEERR;
          if (link->magicnumber != CAPSMAGIC) return CAPS_SOURCEERR;
          if (link->type        != DATASET  ) return CAPS_SOURCEERR;
          if ((link->last.sNum == 0) &&
              (dataset->startup != NULL)) {
            /* special startup data */
            *method = dataset->linkMethod;
            *rank   = dataset->rank;
            *npts   = 1;
            *data   = dataset->startup;
            return CAPS_SUCCESS;
          }
          if ((link->last.sNum > dobject->last.sNum) ||
              (dataset->data == NULL)) return CAPS_DIRTY;
        }
        *method = dataset->linkMethod;
        *rank   = dataset->rank;
        *npts   = dataset->npts;
        *data   = dataset->data;
        *units  = dataset->units;
        return CAPS_SUCCESS;
      }
    }
  }

  return CAPS_NOTFOUND;
}


int
aim_getBounds(void *aimStruc, int *nTname, char ***tnames)
{
  int           i, j, k;
  char          **names;
  aimInfo       *aInfo;
  capsProblem   *problem;
  capsBound     *bound;
  capsAnalysis  *analysis, *analy;
  capsVertexSet *vertexset;

  *tnames = NULL;
  *nTname = 0;
  aInfo   = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  if (problem->nBound == 0) return CAPS_SUCCESS;
  names = (char **) EG_alloc(problem->nBound*sizeof(char *));
  if (names == NULL) return EGADS_MALLOC;

  /* look at all bounds and find those with this analysis */
  for (k = i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]       == NULL) continue;
    if (problem->bounds[i]->name == NULL) continue;
    /* find our analysis in the Bound */
    bound = (capsBound *) problem->bounds[i]->blind;
    if (bound == NULL) continue;
    for (j = 0; j < bound->nVertexSet; j++) {
      if (bound->vertexSet[j]              == NULL)      continue;
      if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
      if (bound->vertexSet[j]->type        != VERTEXSET) continue;
      if (bound->vertexSet[j]->blind       == NULL)      continue;
      vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
      if (vertexset           == NULL) continue;
      if (vertexset->analysis == NULL) continue;
      analy = (capsAnalysis *) vertexset->analysis->blind;
      if (analy != analysis) continue;
      names[k] = problem->bounds[i]->name;
      k++;
      break;
    }
  }
  if (k == 0) {
    EG_free(names);
    return CAPS_SUCCESS;
  }

  *nTname = k;
  *tnames = names;
  return CAPS_SUCCESS;
}


static int
aim_fillAttrs(capsObject *object, int *n, char ***names, capsValue **values)
{
  int       i, j, num, len, atype, slen;
  char      **aName;
  capsValue *value;

  num    = object->attrs->nattrs;
  if (num <= 0) return CAPS_SUCCESS;
  aName  = (char **) EG_alloc(num*sizeof(char *));
  if (aName == NULL) return EGADS_MALLOC;
  value = (capsValue *) EG_alloc(num*sizeof(capsValue));
  if (value == NULL) {
    EG_free(aName);
    return EGADS_MALLOC;
  }
  for (i = 0; i < num; i++) {
    len   = object->attrs->attrs[i].length;
    atype = object->attrs->attrs[i].type;
    value[i].length     = 1;
    value[i].type       = Integer;
    value[i].nrow       = value[i].ncol   = value[i].dim = 0;
    value[i].index      = value[i].pIndex = 0;
    value[i].lfixed     = value[i].sfixed = Fixed;
    value[i].nullVal    = NotAllowed;
    value[i].meshWriter = NULL;
    value[i].units      = NULL;
    value[i].link       = NULL;
    aName[i]            = EG_strdup(object->attrs->attrs[i].name);
    if (aName[i] == NULL) goto bail;
    value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
    value[i].linkMethod      = Copy;
    value[i].gInType         = 0;
    value[i].partial         = NULL;
    value[i].nderiv          = 0;
    value[i].derivs          = NULL;
    if (atype == ATTRINT) {
      if (len == 1) {
        value[i].vals.integer = object->attrs->attrs[i].vals.integer;
      } else {
        value[i].vals.integers = (int *) EG_alloc(len*sizeof(int));
        if (value[i].vals.integers == NULL) goto bail;
        for (j = 0; j < len; j++)
          value[i].vals.integers[j] = object->attrs->attrs[i].vals.integers[j];
      }
    } else if (atype == ATTRREAL) {
      value[i].type = Double;
      if (len == 1) {
        value[i].vals.real = object->attrs->attrs[i].vals.real;
      } else {
        value[i].vals.reals = (double *) EG_alloc(len*sizeof(double));
        if (value[i].vals.reals == NULL) goto bail;
        for (j = 0; j < len; j++)
          value[i].vals.reals[j] = object->attrs->attrs[i].vals.reals[j];
      }
    } else {
      value[i].type = String;
      for (slen = j = 0; j < value[i].length; j++)
        slen += strlen(object->attrs->attrs[i].vals.string + slen) + 1;
      value[i].vals.string = (char*) EG_alloc(slen*sizeof(char));
      if (value[i].vals.string == NULL) goto bail;
      for (j = 0; j < slen; j++)
        value[i].vals.string[j] = object->attrs->attrs[i].vals.string[j];
    }
    value[i].length = len;
  }

  *n      = num;
  *names  = aName;
  *values = value;
  return CAPS_SUCCESS;

bail:
  for (j = 0; j < i; j++) {
    EG_free(aName[j]);
    if (value[j].type == Integer) {
      if (value[j].length == 1) continue;
      EG_free(value[j].vals.integers);
    } else if (value[j].type == Double) {
      if (value[j].length == 1) continue;
      EG_free(value[j].vals.reals);
    } else {
      EG_free(value[j].vals.string);
    }
  }
/*@-dependenttrans@*/
  EG_free(aName);
  EG_free(value);
/*@+dependenttrans@*/
  return EGADS_MALLOC;
}


int
aim_valueAttrs(void *aimStruc, int index, enum capssType stype, int *nValue,
               char ***names, capsValue **values)
{
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsObject   *object = NULL;

  *nValue = 0;
  *names  = NULL;
  *values = NULL;
  aInfo   = (aimInfo *) aimStruc;
  if (aInfo == NULL)                    return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (index < 1)                        return CAPS_BADINDEX;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  if (stype == GEOMETRYIN) {
    if (index > problem->nGeomIn)        return CAPS_BADINDEX;
    object = problem->geomIn[index-1];
  } else if (stype == GEOMETRYOUT) {
    if (index > problem->nGeomOut)       return CAPS_BADINDEX;
    object = problem->geomOut[index-1];
  } else if (stype == ANALYSISIN) {
    if (index > analysis->nAnalysisIn)   return CAPS_BADINDEX;
    object = analysis->analysisIn[index-1];
  } else if (stype == ANALYSISOUT) {
    if (index > analysis->nAnalysisOut)  return CAPS_BADINDEX;
    object = analysis->analysisOut[index-1];
  } else if (stype == ANALYSISDYNO) {
    if (index > analysis->nAnalysisDynO) return CAPS_BADINDEX;
    object = analysis->analysisDynO[index-1];
  } else {
    return CAPS_BADTYPE;
  }
  if (object == NULL)                   return CAPS_NULLOBJ;
  if (object->attrs == NULL)            return CAPS_SUCCESS;

  return aim_fillAttrs(object, nValue, names, values);
}


int
aim_analysisAttrs(void *aimStruc, int *nValue, char ***names,
                  capsValue **values)
{
  int          i;
  aimInfo      *aInfo;
  capsProblem  *problem;
  capsAnalysis *analysis;

  *nValue  = 0;
  *names   = NULL;
  *values  = NULL;
  aInfo    = (aimInfo *) aimStruc;
  if (aInfo == NULL)                    return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  for (i = 0; i < problem->nAnalysis; i++) {
    if (problem->analysis[i] == NULL) continue;
    if (problem->analysis[i]->blind == analysis) break;
  }
  if (i == problem->nAnalysis) return CAPS_NOTFOUND;

  return aim_fillAttrs(problem->analysis[i], nValue, names, values);
}


void
aim_freeAttrs(int nValue, char **names, capsValue *values)
{
  int i;

  for (i = 0; i < nValue; i++) {
    EG_free(names[i]);
    if (values[i].type == Integer) {
      if (values[i].length == 1) continue;
      EG_free(values[i].vals.integers);
    } else if (values[i].type == Double) {
      if (values[i].length == 1) continue;
      EG_free(values[i].vals.reals);
    } else {
      EG_free(values[i].vals.string);
    }
  }
  EG_free(names);
  EG_free(values);
}


int
aim_setSensitivity(void *aimStruc, const char *GIname, int irow, int icol)
{
  int         stat, i, j, k, m, n, ipmtr, nbrch, npmtr, nrow, ncol, type;
  int         buildTo, builtTo, nbody, i_wrt, len_wrt, outLevel;
  char        name[MAX_NAME_LEN];
  double      *reals;
  modl_T      *MODL;
  aimInfo     *aInfo;
  capsProblem *problem;
  capsValue   *value;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (GIname == NULL)                  return CAPS_NULLNAME;
  if ((irow < 1) || (icol < 1))        return CAPS_BADINDEX;
  problem  = aInfo->problem;
  MODL     = (modl_T *) problem->modl;
  if (MODL == NULL)                    return CAPS_NOTPARMTRIC;

  /* find the OpenCSM Parameter */
  stat = ocsmInfo(problem->modl, &nbrch, &npmtr, &nbody);
  if (stat != SUCCESS) return stat;
  for (ipmtr = i = 0; i < npmtr; i++) {
    stat = ocsmGetPmtr(problem->modl, i+1, &type, &nrow, &ncol, name);
    if (stat != SUCCESS) continue;
    if (type != OCSM_DESPMTR) continue;
    if (strcmp(name, GIname) != 0) continue;
    ipmtr = i+1;
    break;
  }
  if (ipmtr == 0)  return CAPS_NOSENSITVTY;
  if (irow > nrow) return CAPS_BADINDEX;
  if (icol > ncol) return CAPS_BADINDEX;

  aInfo->pIndex = 0;
  aInfo->irow   = 0;
  aInfo->icol   = 0;

  /* clear all then set */
  stat = ocsmSetDtime(problem->modl, 0);
  if (stat != SUCCESS) return stat;
  stat = ocsmSetVelD(problem->modl, 0,     0,    0,    0.0);
  if (stat != SUCCESS) return stat;
  stat = ocsmSetVelD(problem->modl, ipmtr, irow, icol, 1.0);
  if (stat != SUCCESS) return stat;
  buildTo = 0;
  nbody   = 0;
  outLevel = ocsmSetOutLevel(0);
  printf(" CAPS Info: Building sensitivity information for: %s[%d,%d]\n", name, irow, icol);
  stat    = ocsmBuild(problem->modl, buildTo, &builtTo, &nbody, NULL);
  fflush(stdout);
  if (stat != SUCCESS) return stat;
  ocsmSetOutLevel(outLevel);

  /* fill in GeometryOut Dot values -- same as caps_geomOutSensit */
  for (i = 0; i < problem->nRegGIN; i++) {
    if (problem->geomIn[problem->regGIN[i].index-1] == NULL) continue;
    value = (capsValue *) problem->geomIn[problem->regGIN[i].index-1]->blind;
    if (value                   == NULL)  continue;
    if (value->pIndex           != ipmtr) continue;
    if (problem->regGIN[i].irow != irow && problem->regGIN[i].irow > 0)  continue;
    if (problem->regGIN[i].icol != icol && problem->regGIN[i].icol > 0)  continue;

    if ((problem->regGIN[i].irow == 0) &&
        (problem->regGIN[i].icol == 0))
      i_wrt = value->ncol*(irow-1) + (icol-1);
    else
      i_wrt = 0;

    for (j = 0; j < problem->nGeomOut; j++) {
      if (problem->geomOut[j] == NULL) continue;
      value = (capsValue *) problem->geomOut[j]->blind;
      if (value                  == NULL) continue;
      if (value->derivs          == NULL) continue;

      len_wrt = value->derivs[i].len_wrt;
      if (value->derivs[i].deriv == NULL)
        value->derivs[i].deriv = (double *) EG_alloc(value->length*len_wrt*sizeof(double));
      if (value->derivs[i].deriv == NULL) { return EGADS_MALLOC; }

      reals = value->vals.reals;
      if (value->length == 1) reals = &value->vals.real;
      for (n = k = 0; k < value->nrow; k++)
        for (m = 0; m < value->ncol; m++, n++)
          ocsmGetValu(problem->modl, value->pIndex, k+1, m+1,
                      &reals[n], &value->derivs[i].deriv[len_wrt*n + i_wrt]);
    }
  }

  aInfo->pIndex = ipmtr;
  aInfo->irow   = irow;
  aInfo->icol   = icol;
  return CAPS_SUCCESS;
}


int
aim_getSensitivity(void *aimStruc, ego tess, int ttype, int index, int *npts,
                   double **dxyz)
{
  int          i, stat, ibody, npt, ntri, type, outLevel;
  double       *dsen;
  const int    *ptype, *pindex, *tris, *tric;
  const double *xyzs, *parms;
  ego          body, oldtess;
  egTessel     *btess;
  modl_T       *MODL;
  aimInfo      *aInfo;
  capsProblem  *problem;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aInfo->pIndex      == 0)         return CAPS_STATEERR;
  if ((ttype < -2) || (ttype > 2))     return CAPS_RANGEERR;
  if (index < 1)                       return CAPS_BADINDEX;
  problem  = aInfo->problem;
  MODL     = (modl_T *) problem->modl;
  if (MODL == NULL)                    return CAPS_NOTPARMTRIC;
  if (tess == NULL)                    return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)      return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION)    return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                   return EGADS_NOTFOUND;
  body = btess->src;
  if (body == NULL)                    return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC)      return EGADS_NOTOBJ;
  if (body->oclass != BODY)            return EGADS_NOTBODY;

  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    if (MODL->body[ibody].ebody   == body) break;
  }
  if (ibody > MODL->nbody) return CAPS_NOTFOUND;

  /* set the OCSM tessellation to the one given & save the old */
  oldtess = MODL->body[ibody].etess;
  MODL->body[ibody].etess = tess;

  /* get the length */
  if (ttype == 0) {
    type = OCSM_NODE;
    npt  = 1;
  } else if (abs(ttype) == 1) {
    type = OCSM_EDGE;
    stat = EG_getTessEdge(MODL->body[ibody].etess, index, &npt, &xyzs, &parms);
    if (stat != EGADS_SUCCESS) {
      MODL->body[ibody].etess = oldtess;
      return stat;
    }
  } else {
    type = OCSM_FACE;
    stat = EG_getTessFace(MODL->body[ibody].etess, index, &npt, &xyzs, &parms,
                          &ptype, &pindex, &ntri, &tris, &tric);
    if (stat != EGADS_SUCCESS) {
      MODL->body[ibody].etess = oldtess;
      return stat;
    }
  }
  if (npt <= 0) {
    MODL->body[ibody].etess = oldtess;
    return CAPS_NULLVALUE;
  }

  /* get the memory to store the sensitivity */
  dsen = (double *) EG_alloc(3*npt*sizeof(double));
  if (dsen == NULL) {
    MODL->body[ibody].etess = oldtess;
    return EGADS_MALLOC;
  }

  /* get and return the requested sensitivity */
  if (ttype <= 0) {
    outLevel = ocsmSetOutLevel(0);
    stat = ocsmGetVel(problem->modl, ibody, type, index, npt, NULL, dsen);
    ocsmSetOutLevel(outLevel);
    if (stat != SUCCESS) {
      EG_free(dsen);
      MODL->body[ibody].etess = oldtess;
      return stat;
    }
  } else {
    outLevel = ocsmSetOutLevel(0);
    stat = ocsmGetTessVel(problem->modl, ibody, type, index, &xyzs);
    ocsmSetOutLevel(outLevel);
    if (stat != SUCCESS) {
      EG_free(dsen);
      MODL->body[ibody].etess = oldtess;
      return stat;
    }
    for (i = 0; i < 3*npt; i++) dsen[i] = xyzs[i];
  }
  /* reset OCSMs tessellation */
  MODL->body[ibody].etess = oldtess;
  if (MODL->dtime != 0)
    printf(" CAPS Info: Sensitivity finite differenced\n");

  *npts = npt;
  *dxyz = dsen;
  return CAPS_SUCCESS;
}


int
aim_tessSensitivity(void *aimStruc, const char *GIname, int irow, int icol,
                    ego tess, int *npts, double **dxyz)
{
  int          i, j, ipmtr, ibody, stat, nbrch, npmtr, nbody, nrow, ncol, type;
  int          npt, nface, np, ntris, global, outLevel;
  const int    *ptype, *pindex, *tris, *tric;
  char         name[MAX_NAME_LEN];
  double       *dsen;
  const double *xyzs, *xyz, *uv;
  ego          body, oldtess;
  egTessel     *btess;
  modl_T       *MODL;
  aimInfo      *aInfo;
  capsProblem  *problem;

  *npts = 0;
  *dxyz = NULL;
  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (GIname == NULL)                  return CAPS_NULLNAME;
  if ((irow < 1) || (icol < 1))        return CAPS_BADINDEX;
  problem  = aInfo->problem;
  MODL     = (modl_T *) problem->modl;
  if (MODL == NULL)                    return CAPS_NOTPARMTRIC;
  if (tess == NULL)                    return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)      return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION)    return EGADS_NOTTESS;
  stat = EG_statusTessBody(tess, &body, &i, &npt);
  if (stat == EGADS_OUTSIDE)           return EGADS_TESSTATE;
  if (body == NULL)                    return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC)      return EGADS_NOTOBJ;
  if (body->oclass != BODY)            return EGADS_NOTBODY;
  stat = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
  if (stat != EGADS_SUCCESS)           return stat;

  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    if (MODL->body[ibody].ebody   == body) break;
  }
  if (ibody > MODL->nbody) return CAPS_NOTFOUND;

  /* find the OpenCSM Parameter */
  stat = ocsmInfo(problem->modl, &nbrch, &npmtr, &nbody);
  if (stat != SUCCESS) return stat;
  for (ipmtr = i = 0; i < npmtr; i++) {
    stat = ocsmGetPmtr(problem->modl, i+1, &type, &nrow, &ncol, name);
    if (stat != SUCCESS) continue;
    if (type != OCSM_DESPMTR) continue;
    if (strcmp(name, GIname) != 0) continue;
    ipmtr = i+1;
    break;
  }
  if (ipmtr == 0)  return CAPS_NOSENSITVTY;
  if (irow > nrow) return CAPS_BADINDEX;
  if (icol > ncol) return CAPS_BADINDEX;

  /* set the Parameter if not already set */
  if ((aInfo->pIndex != ipmtr) || (aInfo->irow != irow) ||
      (aInfo->icol   != icol)) {
    stat = aim_setSensitivity(aimStruc, GIname, irow, icol);
    if (stat != CAPS_SUCCESS) return stat;
  }

  dsen = (double *) EG_alloc(3*npt*sizeof(double));
  if (dsen == NULL) return EGADS_MALLOC;
  for (i = 0; i < 3*npt; i++) dsen[i] = 0.0;

  /* return the requested sensitivity */
  oldtess = MODL->body[ibody].etess;
  MODL->body[ibody].etess = tess;

  btess = (egTessel *) tess->blind;
  if (btess->nFace == 0) {
    for (i = 1; i <= btess->nEdge; i++) {
      stat = EG_getTessEdge(tess, i, &np, &xyz, &uv);
      if (stat != EGADS_SUCCESS) {
        EG_free(dsen);
        return stat;
      }
      outLevel = ocsmSetOutLevel(0);
      stat = ocsmGetTessVel(problem->modl, ibody, OCSM_EDGE, i, &xyzs);
      ocsmSetOutLevel(outLevel);
      if (stat != SUCCESS) {
        EG_free(dsen);
        return stat;
      }
      for (j = 1; j <= np; j++) {
        stat = EG_localToGlobal(tess, -i, j, &global);
        if (stat != EGADS_SUCCESS) {
          EG_free(dsen);
          return stat;
        }
        dsen[3*global-3] = xyzs[3*j-3];
        dsen[3*global-2] = xyzs[3*j-2];
        dsen[3*global-1] = xyzs[3*j-1];
      }
    }
  } else {
    for (i = 1; i <= nface; i++) {
      stat = EG_getTessFace(tess, i, &np, &xyz, &uv, &ptype, &pindex, &ntris,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        EG_free(dsen);
        return stat;
      }
      outLevel = ocsmSetOutLevel(0);
      stat = ocsmGetTessVel(problem->modl, ibody, OCSM_FACE, i, &xyzs);
      ocsmSetOutLevel(outLevel);
      if (stat != SUCCESS) {
        EG_free(dsen);
        return stat;
      }
      for (j = 1; j <= np; j++) {
        stat = EG_localToGlobal(tess, i, j, &global);
        if (stat != EGADS_SUCCESS) {
          EG_free(dsen);
          return stat;
        }
        dsen[3*global-3] = xyzs[3*j-3];
        dsen[3*global-2] = xyzs[3*j-2];
        dsen[3*global-1] = xyzs[3*j-1];
      }
    }
  }
  MODL->body[ibody].etess = oldtess;
  if (MODL->dtime != 0)
    printf(" CAPS Info: Sensitivity finite difference used for: %s[%d,%d]\n", name, irow, icol);

  *npts = npt;
  *dxyz = dsen;
  return CAPS_SUCCESS;
}


int
aim_isNodeBody(ego body, double *xyz)
{
  int    status, oclass, mtype, nchild, *psens = NULL;
  double data[4];
  ego    *pchldrn, ref, loop, edge;

  /* Determine of the body is a "node body", i.e. a degenerate wire body */
  status = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchild, &pchldrn,
                          &psens);
  if (status != EGADS_SUCCESS) return status;

  /* should be a body */
  if (oclass != BODY) return EGADS_NOTBODY;

  /* done if it's not a wire body */
  if (mtype != WIREBODY) return EGADS_OUTSIDE;

  /* there should be only one child, and it is a loop */
  if (nchild != 1) return EGADS_OUTSIDE;

  loop = pchldrn[0];

  /* get the edges */
  status = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nchild, &pchldrn,
                          &psens);
  if (status != CAPS_SUCCESS) return status;

  /* there should be only one child, and it is an edge */
  if (nchild != 1) return EGADS_OUTSIDE;

  edge = pchldrn[0];

  /* get the node(s) of the edge */
  status = EG_getTopology(edge, &ref, &oclass, &mtype, data, &nchild, &pchldrn,
                          &psens);
  if (status != CAPS_SUCCESS) return status;

  /* the edge should be degenerate */
  if (mtype != DEGENERATE) return EGADS_OUTSIDE;

  /* something is really wrong if the edge does not have one node */
  if (nchild != 1) return EGADS_GEOMERR;

  /* retrieve the coordinates of the node */
  status = EG_getTopology(pchldrn[0], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens);
  if (status != EGADS_SUCCESS) return status;

  return EGADS_SUCCESS;
}


static void
aim_addErrorLine(void *aimStruc, enum capseType etype, const char *line)
{
  int       index, len;
  char      **ltmp;
  capsError *etmp;
  aimInfo   *aInfo;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return;
  if (aInfo->magicnumber != CAPSMAGIC) return;

  if (etype == CONTINUATION) {
    if (aInfo->errs.nError == 0) {
      printf(" CAPS Internal: Continuation without a start!\n");
      return;
    }
    index = aInfo->errs.nError - 1;
    len   = aInfo->errs.errors[index].nLines + 1;
    ltmp  = (char **) EG_reall(aInfo->errs.errors[index].lines,
                               len*sizeof(char *));
    if (ltmp == NULL) return;
    ltmp[len-1] = EG_strdup(line);
    aInfo->errs.errors[index].lines  = ltmp;
    aInfo->errs.errors[index].nLines = len;
    return;
  }

  index = aInfo->errs.nError;
  if (index == 0) {
    aInfo->errs.errors = (capsError *) EG_alloc(sizeof(capsError));
    if (aInfo->errs.errors == NULL) return;
  } else {
    etmp = (capsError *) EG_reall(aInfo->errs.errors,
                                  (index+1)*sizeof(capsError));
    if (etmp == NULL) return;
    aInfo->errs.errors = etmp;
  }
  aInfo->errs.errors[index].errObj = NULL;
  aInfo->errs.errors[index].eType  = etype;
  aInfo->errs.errors[index].index  = 0;
  aInfo->errs.errors[index].nLines = 1;
  aInfo->errs.errors[index].lines  = (char **) EG_alloc(sizeof(char *));
  if (aInfo->errs.errors[index].lines != NULL)
    aInfo->errs.errors[index].lines[0] = EG_strdup(line);

  aInfo->errs.nError = index+1;
}


void
aim_status(void *aimInfo, const int status, const char *file,
           const int line, const char *func, const int narg, ...)
{
  va_list args;
  char    buffer[EBUFSIZE] = {'\0'}, buffer2[EBUFSIZE] = {'\0'};
  const char *format = NULL;

  snprintf(buffer2, EBUFSIZE, "%s:%d in %s(): Error status = %d",
           file, line, func, status);
  aim_addErrorLine(aimInfo, CSTAT, buffer2);

  if (narg > 0) {
    va_start(args, narg);
    format = va_arg(args, const char *);
    vsnprintf(buffer, EBUFSIZE, format, args);
    va_end(args);

    snprintf(buffer2, EBUFSIZE, "%s", buffer);
    aim_addErrorLine(aimInfo, CONTINUATION, buffer2);
  }
}


void
aim_message(void *aimStruc, enum capseType etype, int index, const char *file,
            const int line, const char *func, const char *format, ...)
{
  int     i;
  aimInfo *aInfo;
  capsAnalysis *analysis;

  va_list args;
  char    buffer[EBUFSIZE] = {'\0'}, buffer2[EBUFSIZE] = {'\0'};
  const char  *eType[4] = {"Info",
                           "Warning",
                           "Error",
                           "Possible Developer Error"};

  // clang-analyzer
  if (etype == CONTINUATION) return;

  snprintf(buffer2, EBUFSIZE, "%s:%d in %s():", file, line, func);
  aim_addErrorLine(aimStruc, etype, buffer2);

  if (index > 0) {
    aInfo = (aimInfo *) aimStruc;
    if (aInfo == NULL)                   return;
    if (aInfo->magicnumber != CAPSMAGIC) return;

    analysis = (capsAnalysis *) aInfo->analysis;
    if (index <= analysis->nAnalysisIn) {
      snprintf(buffer, EBUFSIZE, "ANALYSISIN: %s",
               analysis->analysisIn[index-1]->name);
      aim_addErrorLine(aimStruc, CONTINUATION, buffer);
    }
  }

  va_start(args, format);
  vsnprintf(buffer, EBUFSIZE, format, args);
  va_end(args);

  snprintf(buffer2, EBUFSIZE, "%s: %s", eType[etype], buffer);
  aim_addErrorLine(aimStruc, CONTINUATION, buffer2);

  if (index > 0) {
    i = aInfo->errs.nError-1;
    aInfo->errs.errors[i].index = index;
  }
}


void
aim_addLine(void *aimInfo, const char *format, ...)
{
  va_list args;
  char    buffer[EBUFSIZE] = {'\0'};

  va_start(args, format);
  vsnprintf(buffer, EBUFSIZE, format, args);
  va_end(args);

  aim_addErrorLine(aimInfo, CONTINUATION, buffer);
}


void
aim_setIndexError(void *aimStruc, int index)
{
  int     i;
  aimInfo *aInfo;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return;
  if (aInfo->magicnumber != CAPSMAGIC) return;
  if (aInfo->errs.nError == 0)         return;

  i = aInfo->errs.nError - 1;
  aInfo->errs.errors[i].index = index;
}


void
aim_removeError(void *aimStruc)
{
  int     i, j, k;
  aimInfo *aInfo;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                   return;
  if (aInfo->magicnumber != CAPSMAGIC) return;

  for (j = i = 0; i < aInfo->errs.nError; i++)
    if (aInfo->errs.errors[i].eType == CERROR ||
        aInfo->errs.errors[i].eType == CSTAT) {
      for (k = 0; k < aInfo->errs.errors[i].nLines; k++)
        EG_free(aInfo->errs.errors[i].lines[k]);
      EG_free(aInfo->errs.errors[i].lines);
    } else {
      aInfo->errs.errors[j] = aInfo->errs.errors[i];
      j++;
    }

  if (j == 0) {
    EG_free(aInfo->errs.errors);
    aInfo->errs.errors = NULL;
  }
  aInfo->errs.nError = j;
}
