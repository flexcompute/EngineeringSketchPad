/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Base Object Functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef WIN32
/* needs Advapi32.lib & Ws2_32.lib */
#include <Windows.h>
#define getpid     _getpid
#define snprintf   _snprintf
#define strcasecmp  stricmp
#define access     _access
#define PATH_MAX   _MAX_PATH
#define F_OK       0
#else
#include <unistd.h>
#include <dirent.h>
#include <dlfcn.h>
#include <limits.h>
#include <sys/stat.h>
#endif

#include "udunits.h"

#include "capsTypes.h"
#include "capsFunIDs.h"
#include "capsAIM.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"


#define STRING(a)       #a
#define STR(a)          STRING(a)

static char  *CAPSprop[2] = {STR(CAPSPROP),
                   "\nCAPSprop: Copyright 2014-2022 MIT. All Rights Reserved."};

/*@-incondefs@*/
extern void ut_free(/*@only@*/ ut_unit* const unit);
/*@+incondefs@*/

extern /*@null@*/ /*@only@*/ char *EG_strdup(/*@null@*/ const char *str);

extern void caps_jrnlWrite(int funID, capsProblem *problem, capsObject *obj,
                           int status, int nargs, capsJrnl *args, CAPSLONG sNum0,
                           CAPSLONG sNum);
extern int  caps_jrnlEnd(capsProblem *problem);
extern int  caps_jrnlRead(int funID, capsProblem *problem, capsObject *obj,
                          int nargs, capsJrnl *args, CAPSLONG *sNum, int *stat);
extern int  caps_filter(capsProblem *problem, capsAnalysis *analysis);
extern int  caps_Aprx1DFree(/*@only@*/ capsAprx1D *approx);
extern int  caps_Aprx2DFree(/*@only@*/ capsAprx2D *approx);
extern int  caps_build(capsObject *pobject, int *nErr, capsErrs **errors);
extern int  caps_dumpGeomVals(capsProblem *problem, int flag);
extern void caps_freeFList(capsObject *obj);
extern int  caps_analysisInfX(const capsObj aobject, char **apath, char **unSys,
                              int *major, int *minor, char **intents,
                              int *nField, char ***fnames, int **ranks,
                              int **fInOut, int *execution, int *status);
extern int  caps_execX(capsObject *aobject, int *nErr, capsErrs **errors);
extern int  caps_writeValueObj(capsProblem *problem, capsObject *valobj);
extern int  caps_writeBound(const capsObject *bobject);
extern int  caps_writeAnalysisObj(capsProblem *problem, capsObject *aobject);

/* used internally */
int caps_freeError(/*@only@*/ capsErrs *errs);



#ifdef WIN32
static int caps_flipSlash(const char *src, char *back)
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


static int caps_statFileX(const char *path)
{
  int   stat;
#ifdef WIN32
  DWORD dwAttr;
#else
  DIR   *dr;
#endif

/*@-unrecog@*/
  stat = access(path, F_OK);
/*@+unrecog@*/
  if (stat == -1) return EGADS_NOTFOUND;

  /* are we a directory */
#ifdef WIN32
  dwAttr = GetFileAttributes(path);
  if ((dwAttr != (DWORD) -1) && ((dwAttr & FILE_ATTRIBUTE_DIRECTORY) != 0))
    return EGADS_OUTSIDE;
#else
  dr = opendir(path);
  if (dr != NULL) {
    closedir(dr);
    return EGADS_OUTSIDE;
  }
#endif

  return EGADS_SUCCESS;
}


int caps_statFile(const char *path)
{
#ifdef WIN32
  int  stat;
  char back[PATH_MAX];

  stat = caps_flipSlash(path, back);
  if (stat != EGADS_SUCCESS) return stat;
  return caps_statFileX(back);
#else
  return caps_statFileX(path);
#endif
}


int caps_rmFile(const char *path)
{
#ifdef WIN32
  int  stat;
  char back[PATH_MAX];

  stat = caps_flipSlash(path, back);
  if (stat != EGADS_SUCCESS) return stat;
  if (_unlink(back) == -1) return EGADS_NOTFOUND;
#else
  if ( unlink(path) == -1) return EGADS_NOTFOUND;
#endif

  return EGADS_SUCCESS;
}


int caps_rmDir(const char *path)
{
  int  stat;
/*@-unrecog@*/
  char cmd[PATH_MAX+15];
/*@+unrecog@*/
#ifdef WIN32
  char back[PATH_MAX];

  stat = caps_flipSlash(path, back);
  if (stat != EGADS_SUCCESS) return stat;
  stat = caps_statFileX(back);
  if (stat != EGADS_OUTSIDE) return EGADS_NOTFOUND;
  snprintf(cmd, PATH_MAX+15, "rmdir /Q /S \"%s\"", path);
  fflush(NULL);
#else
  stat = caps_statFileX(path);
  if (stat != EGADS_OUTSIDE) return EGADS_NOTFOUND;
  snprintf(cmd, PATH_MAX+15, "rm -rf '%s'", path);
#endif
  stat = system(cmd);
  if ((stat == -1) || (stat == 127)) return EGADS_REFERCE;
  
  return EGADS_SUCCESS;
}


void caps_rmWild(const char *path, const char *wild)
{
/*@-unrecog@*/
  char cmd[PATH_MAX+16];
/*@+unrecog@*/
#ifdef WIN32
  int  stat;
  char back[PATH_MAX];

  stat = caps_flipSlash(path, back);
  if (stat != EGADS_SUCCESS) {
    printf(" CAPS Error: caps_flipSlash = %d (caps_rmWild)!\n", stat);
    return;
  }
  snprintf(cmd, PATH_MAX+16, "del /Q \"%s\"\\%s 2>NUL", path, wild);
  fflush(NULL);
#else
  snprintf(cmd, PATH_MAX+16,  "rm -f '%s'/%s", path, wild);
#endif
  (void) system(cmd);
}


int caps_mkDir(const char *path)
{
  int  stat;
#ifdef WIN32
  char back[PATH_MAX];

  stat = caps_flipSlash(path, back);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(back);
  if (stat != EGADS_NOTFOUND) return EGADS_EXISTS;
  stat = mkdir(back);
#else
  stat = caps_statFileX(path);
  if (stat != EGADS_NOTFOUND) return EGADS_EXISTS;
  stat = mkdir(path, 0777);
#endif
  if (stat == -1) return EGADS_NOTFOUND;

  return EGADS_SUCCESS;
}


int caps_cpFile(const char *src, const char *dst)
{
  int  stat;
  char cmd[2*PATH_MAX+14];
#ifdef WIN32
  char sback[PATH_MAX], dback[PATH_MAX];

  stat = caps_flipSlash(src, sback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(sback);
  if (stat != EGADS_SUCCESS)  return EGADS_NOTFOUND;
  stat = caps_flipSlash(dst, dback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(dback);
  if (stat != EGADS_NOTFOUND) return EGADS_NOTFOUND;
  snprintf(cmd, 2*PATH_MAX+14, "copy /Y \"%s\" \"%s\"", sback, dback);
  fflush(NULL);
#else
  stat = caps_statFileX(src);
  if (stat != EGADS_SUCCESS)  return EGADS_NOTFOUND;
  stat = caps_statFileX(dst);
  if (stat != EGADS_NOTFOUND) return EGADS_NOTFOUND;
  snprintf(cmd, 2*PATH_MAX+14, "cp '%s' '%s'", src, dst);
#endif
  stat = system(cmd);
  if (stat != 0) return EGADS_REFERCE;
  
  return EGADS_SUCCESS;
}


int caps_cpDir(const char *src, const char *dst)
{
  int  stat;
  char cmd[2*PATH_MAX+24];
#ifdef WIN32
  char sback[PATH_MAX], dback[PATH_MAX];

  stat = caps_flipSlash(src, sback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(sback);
  if (stat != EGADS_OUTSIDE)  return EGADS_NOTFOUND;
  stat = caps_flipSlash(dst, dback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(dback);
  if (stat != EGADS_NOTFOUND) return EGADS_NOTFOUND;
  snprintf(cmd, 2*PATH_MAX+24, "xcopy \"%s\" \"%s\" /I /E /Q /Y", sback, dback);
  fflush(NULL);
#else
  stat = caps_statFileX(src);
  if (stat != EGADS_OUTSIDE)  return EGADS_NOTFOUND;
  stat = caps_statFileX(dst);
  if (stat != EGADS_NOTFOUND) return EGADS_NOTFOUND;
  snprintf(cmd, 2*PATH_MAX+24, "cp -R -p '%s' '%s'", src, dst);
#endif
  stat = system(cmd);
  if (stat != 0) return EGADS_REFERCE;
  
  return EGADS_SUCCESS;
}


int caps_rename(const char *src, const char *dst)
{
  int  stat;
#ifdef WIN32
  char sback[PATH_MAX], dback[PATH_MAX];
  
  stat = caps_flipSlash(src, sback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(sback);
  if (stat == EGADS_NOTFOUND) return stat;
  stat = caps_flipSlash(dst, dback);
  if (stat != EGADS_SUCCESS)  return stat;
  stat = caps_statFileX(dback);
  if (stat != EGADS_NOTFOUND) {
    if (stat == EGADS_SUCCESS) {
      caps_rmFile(dback);
    } else {
      caps_rmDir(dback);
    }
  }
  stat = rename(sback, dback);
#else
  stat = rename(src,   dst);
#endif
  if (stat != 0) return CAPS_DIRERR;
  
  return EGADS_SUCCESS;
}


void caps_getStaticStrings(char ***signature, char **pID, char **user)
{
#ifdef WIN32
  int     len;
  WSADATA wsaData;
#endif
  char    name[257], ID[301];

  *signature = CAPSprop;
  
#ifdef WIN32
  WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif
  gethostname(name, 257);
#ifdef WIN32
  WSACleanup();
#endif
  snprintf(ID, 301, "%s:%d", name, getpid());
  *pID = EG_strdup(ID);

#ifdef WIN32
  len = 256;
  GetUserName(name, &len);
#else
/*@-mustfreefresh@*/
  snprintf(name, 257, "%s", getlogin());
/*@+mustfreefresh@*/
#endif
  *user = EG_strdup(name);
}


void caps_fillDateTime(short *datetime)
{
#ifdef WIN32
  SYSTEMTIME sysTime;

  GetLocalTime(&sysTime);
  datetime[0] = sysTime.wYear;
  datetime[1] = sysTime.wMonth;
  datetime[2] = sysTime.wDay;
  datetime[3] = sysTime.wHour;
  datetime[4] = sysTime.wMinute;
  datetime[5] = sysTime.wSecond;

#else
  time_t     tloc;
  struct tm  tim, *dum;
  
  tloc = time(NULL);
/*@-unrecog@*/
  dum  = localtime_r(&tloc, &tim);
/*@+unrecog@*/
  if (dum == NULL) {
    datetime[0] = 1900;
    datetime[1] = 0;
    datetime[2] = 0;
    datetime[3] = 0;
    datetime[4] = 0;
    datetime[5] = 0;
  } else {
    datetime[0] = tim.tm_year + 1900;
    datetime[1] = tim.tm_mon  + 1;
    datetime[2] = tim.tm_mday;
    datetime[3] = tim.tm_hour;
    datetime[4] = tim.tm_min;
    datetime[5] = tim.tm_sec;
  }
#endif
}


void caps_fillLengthUnits(capsProblem *problem, ego body, char **lunits)
{
  int          status, atype, alen;
  const int    *aints;
  const double *areals;
  const char   *astr;
  ut_unit      *utunit;

  *lunits = NULL;
  status  = EG_attributeRet(body, "capsLength", &atype, &alen, &aints, &areals,
                            &astr);
  if ((status != EGADS_SUCCESS) && (status != EGADS_NOTFOUND)) {
    printf(" CAPS Warning: EG_attributeRet = %d in fillLengthUnits!\n", status);
    return;
  }
  if (status == EGADS_NOTFOUND) return;
  if (atype != ATTRSTRING) {
    printf(" CAPS Warning: capsLength w/ incorrect type in fillLengthUnits!\n");
    return;
  }
  
  utunit = ut_parse((ut_system *) problem->utsystem, astr, UT_ASCII);
  if (utunit == NULL) {
    printf(" CAPS Warning: capsLength %s is not a valid unit!\n", astr);
    return;
  }
/*
  printf(" Body with length units = %s\n", astr);
 */
  ut_free(utunit);
  *lunits = EG_strdup(astr);
}


void caps_geomOutUnits(char *name, /*@null@*/ char *lunits, char **units)
{
  int         i, len;
  static char *names[21] = {"@xmin", "@xmax", "@ymin", "@ymax", "@zmin",
                            "@zmax", "@length", "@area", "@volume", "@xcg",
                            "@ycg", "@zcg", "@Ixx", "@Ixy", "@Ixz", "@Iyx",
                            "@Iyy", "@Iyz", "@Izx", "@Izy", "@Izz" };
  static int  power[21]  = {1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 4, 4, 4, 4, 4,
                            4, 4, 4, 4};

  *units = NULL;
  if (lunits == NULL) return;
  
  for (i = 0; i < 21; i++)
    if (strcmp(name, names[i]) == 0) {
      if (power[i] == 1) {
        *units = EG_strdup(lunits);
      } else {
        len    = strlen(lunits) + 3;
        *units = (char *) EG_alloc(len*sizeof(char));
        if (*units == NULL) return;
        snprintf(*units, len, "%s^%d", lunits, power[i]);
      }
      return;
    }
}


int caps_makeTuple(int n, capsTuple **tuple)
{
  int       i;
  capsTuple *tmp;
  
  *tuple = NULL;
  if (n < 1) return CAPS_RANGEERR;
  
  tmp = (capsTuple *) EG_alloc(n*sizeof(capsTuple));
  if (tmp == NULL) return EGADS_MALLOC;
  
  for (i = 0; i < n; i++) tmp[i].name = tmp[i].value = NULL;
  *tuple = tmp;
  
  return CAPS_SUCCESS;
}


void caps_freeTuple(int n, /*@only@*/ capsTuple *tuple)
{
  int i;
  
  if (tuple == NULL) return;
  for (i = 0; i < n; i++) {
    if (tuple[i].name  != NULL) EG_free(tuple[i].name);
    if (tuple[i].value != NULL) EG_free(tuple[i].value);
  }
  EG_free(tuple);
}


void caps_freeOwner(capsOwn *own)
{
  own->index = -1;
  if (own->pname != NULL) {
    EG_free(own->pname);
    own->pname = NULL;
  }
  if (own->pID != NULL) {
    EG_free(own->pID);
    own->pID = NULL;
  }
  if (own->user != NULL) {
    EG_free(own->user);
    own->user = NULL;
  }
}


void caps_freeHistory(capsObject *obj)
{
  int i;
  
  if ((obj->nHistory <= 0) || (obj->history == NULL)) return;
  
  for (i = 0; i < obj->nHistory; i++) {
    EG_free(obj->history[i].pname);
    EG_free(obj->history[i].pID);
    EG_free(obj->history[i].user);
  }
  EG_free(obj->history);
  obj->history  = NULL;
  obj->nHistory = 0;
}


void caps_freeAttrs(egAttrs **attrx)
{
  int     i;
  egAttrs *attrs;
  
  attrs = *attrx;
  if (attrs == NULL) return;
  
  *attrx = NULL;
  /* remove any attributes */
  for (i = 0; i < attrs->nseqs; i++) {
    EG_free(attrs->seqs[i].root);
    EG_free(attrs->seqs[i].attrSeq);
  }
  if (attrs->seqs != NULL) EG_free(attrs->seqs);
  for (i = 0; i < attrs->nattrs; i++) {
    if (attrs->attrs[i].name != NULL) EG_free(attrs->attrs[i].name);
    if (attrs->attrs[i].type == ATTRINT) {
      if (attrs->attrs[i].length > 1) EG_free(attrs->attrs[i].vals.integers);
    } else if (attrs->attrs[i].type == ATTRREAL) {
      if (attrs->attrs[i].length > 1) EG_free(attrs->attrs[i].vals.reals);
    } else {
      EG_free(attrs->attrs[i].vals.string);
    }
  }
  EG_free(attrs->attrs);
  EG_free(attrs);
}


void caps_freeValueObjects(int vflag, int nObjs, capsObject **objects)
{
  int       i, j;
  capsValue *value, *vArray;

  if (objects == NULL) return;
  vArray = (capsValue *) objects[0]->blind;

  for (i = 0; i < nObjs; i++) {
    /* remove any journal read pointers */
    caps_freeFList(objects[i]);
    /* clean up any allocated data arrays */
    value = (capsValue *) objects[i]->blind;
    if (value != NULL) {
      objects[i]->blind = NULL;
      if ((value->type == Boolean) || (value->type == Integer)) {
        if (value->length > 1) EG_free(value->vals.integers);
      } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
        if (value->length > 1) EG_free(value->vals.reals);
      } else if (value->type == String) {
        EG_free(value->vals.string);
      } else if (value->type == Tuple) {
        caps_freeTuple(value->length, value->vals.tuple);
      } else {
        /* pointer type -- nothing should be done here... */
      }
      if (value->units      != NULL) EG_free(value->units);
      if (value->meshWriter != NULL) EG_free(value->meshWriter);
      if (value->partial    != NULL) EG_free(value->partial);
      if (value->derivs     != NULL) {
        for (j = 0; j < value->nderiv; j++) {
          if (value->derivs[j].name  != NULL) EG_free(value->derivs[j].name);
          if (value->derivs[j].deriv != NULL) EG_free(value->derivs[j].deriv);
        }
        EG_free(value->derivs);
      }
      if (vflag == 1) EG_free(value);
    }
    
    /* cleanup and invalidate the object */
    caps_freeHistory(objects[i]);
    caps_freeAttrs(&objects[i]->attrs);
    caps_freeOwner(&objects[i]->last);
    objects[i]->magicnumber = 0;
    EG_free(objects[i]->name);
    objects[i]->name = NULL;
    EG_free(objects[i]);
  }
  
  if (vflag == 0) EG_free(vArray);
  EG_free(objects);
}


void caps_freeEleType(capsEleType *eletype)
{
  eletype->nref  = 0;
  eletype->ndata = 0;
  eletype->ntri  = 0;
  eletype->nseg  = 0;
  eletype->nmat  = 0;

  EG_free(eletype->gst);
  eletype->gst = NULL;
  EG_free(eletype->dst);
  eletype->dst = NULL;
  EG_free(eletype->matst);
  eletype->matst = NULL;
  EG_free(eletype->tris);
  eletype->tris = NULL;
  EG_free(eletype->segs);
  eletype->segs = NULL;
}


void caps_initDiscr(capsDiscr *discr)
{
  discr->dim       = 0;
  discr->instStore = NULL;
  discr->aInfo     = NULL;
  discr->nPoints   = 0;
  discr->nVerts    = 0;
  discr->verts     = NULL;
  discr->celem     = NULL;
  discr->nDtris    = 0;
  discr->dtris     = NULL;
  discr->nDsegs    = 0;
  discr->dsegs     = NULL;
  discr->nTypes    = 0;
  discr->types     = NULL;
  discr->nBodys    = 0;
  discr->bodys     = NULL;
  discr->tessGlobal= NULL;
  discr->ptrm      = NULL;
}


void caps_freeDiscr(capsDiscr *discr)
{
  int           i;
  capsBodyDiscr *discBody = NULL;

  /* free up this capsDiscr */

  EG_free(discr->verts);
  discr->verts = NULL;
  EG_free(discr->celem);
  discr->celem = NULL;
  EG_free(discr->dtris);
  discr->dtris = NULL;
  EG_free(discr->dsegs);
  discr->dsegs = NULL;

  discr->nPoints = 0;
  discr->nVerts  = 0;
  discr->nDtris  = 0;
  discr->nDsegs  = 0;

  /* free Types */
  if (discr->types != NULL) {
    for (i = 0; i < discr->nTypes; i++)
      caps_freeEleType(discr->types + i);
    EG_free(discr->types);
    discr->types = NULL;
  }
  discr->nTypes = 0;

  EG_free(discr->tessGlobal);
  discr->tessGlobal = NULL;

  /* free Body discretizations */
  for (i = 0; i < discr->nBodys; i++) {
    discBody = &discr->bodys[i];
    EG_free(discBody->elems);
    EG_free(discBody->gIndices);
    EG_free(discBody->dIndices);
    EG_free(discBody->poly);
  }
  EG_free(discr->bodys);
  discr->bodys  = NULL;
  discr->nBodys = 0;

  /* aim must free discr->ptrm */
  if (discr->ptrm != NULL)
    printf(" CAPS Warning: discr->ptrm is not NULL (caps_freeDiscr)!\n");
}


void caps_freeAnalysis(int flag, /*@only@*/ capsAnalysis *analysis)
{
  int         i, j, state, npts;
  const char  *eType[4] = {"Info",
                           "Warning",
                           "Error",
                           "Possible Developer Error"};
  ego body;
  capsProblem *problem;

  if (analysis == NULL) return;
  
  for (i = 0; i < analysis->info.wCntxt.aimWriterNum; i++) {
    EG_free(analysis->info.wCntxt.aimWriterName[i]);
#ifdef WIN32
    FreeLibrary(analysis->info.wCntxt.aimWriterDLL[i]);
#else
    dlclose(analysis->info.wCntxt.aimWriterDLL[i]);
#endif
  }
  problem = analysis->info.problem;
  if (analysis->instStore != NULL)
    aim_cleanup(problem->aimFPTR, analysis->loadName, analysis->instStore);
  for (i = 0; i < analysis->nField; i++) EG_free(analysis->fields[i]);
  EG_free(analysis->fields);
  EG_free(analysis->ranks);
  EG_free(analysis->fInOut);
  if (analysis->intents  != NULL) EG_free(analysis->intents);
  if (analysis->loadName != NULL) EG_free(analysis->loadName);
  if (analysis->unitSys  != NULL) EG_free(analysis->unitSys);
  if (analysis->fullPath != NULL) EG_free(analysis->fullPath);
  if (analysis->path     != NULL) EG_free(analysis->path);
  if (analysis->bodies   != NULL) EG_free(analysis->bodies);

  if (analysis->tess != NULL) {
    for (j = 0; j < analysis->nTess; j++)
      if (analysis->tess[j] != NULL) {
        /* delete the body in the tessellation if not on the OpenCSM stack */
        body = NULL;
        if (j >= analysis->nBody) {
          (void) EG_statusTessBody(analysis->tess[j], &body, &state, &npts);
        }
        EG_deleteObject(analysis->tess[j]);
        analysis->tess[j] = NULL;
        if (body != NULL) EG_deleteObject(body);
      }
    EG_free(analysis->tess);
    analysis->tess   = NULL;
    analysis->nTess  = 0;
  }

  if (analysis->info.errs.errors != NULL) {
    printf(" Note: Lost AIM Communication ->\n");
    for (i = 0; i < analysis->info.errs.nError; i++) {
      for (j = 0; j < analysis->info.errs.errors[i].nLines; j++) {
        if (j == 0) {
          printf("   %s: %s\n", eType[analysis->info.errs.errors[i].eType],
                 analysis->info.errs.errors[i].lines[j]);
        } else {
          printf("            %s\n", analysis->info.errs.errors[i].lines[j]);
        }
        EG_free(analysis->info.errs.errors[i].lines[j]);
      }
      EG_free(analysis->info.errs.errors[i].lines);
    }
    EG_free(analysis->info.errs.errors);
    analysis->info.errs.nError = 0;
    analysis->info.errs.errors = NULL;
  }
  if (flag == 1) return;

  if (analysis->analysisIn   != NULL)
    caps_freeValueObjects(0, analysis->nAnalysisIn,   analysis->analysisIn);
  
  if (analysis->analysisOut  != NULL)
    caps_freeValueObjects(0, analysis->nAnalysisOut,  analysis->analysisOut);
  
  if (analysis->analysisDynO != NULL)
    caps_freeValueObjects(1, analysis->nAnalysisDynO, analysis->analysisDynO);
  
  caps_freeOwner(&analysis->pre);
  EG_free(analysis);
}


int caps_makeObject(capsObject **objs)
{
  capsObject *objects;

  *objs   = NULL;
  objects = (capsObject *) EG_alloc(sizeof(capsObject));
  if (objects == NULL) return EGADS_MALLOC;

  objects->magicnumber = CAPSMAGIC;
  objects->type        = UNUSED;
  objects->subtype     = NONE;
  objects->delMark     = 0;
  objects->name        = NULL;
  objects->attrs       = NULL;
  objects->blind       = NULL;
  objects->flist       = NULL;
  objects->parent      = NULL;
  objects->nHistory    = 0;
  objects->history     = NULL;
  objects->last.index  = -1;
  objects->last.pname  = NULL;
  objects->last.pID    = NULL;
  objects->last.user   = NULL;
  objects->last.sNum   = 0;
  caps_fillDateTime(objects->last.datetime);

  *objs = objects;
  return CAPS_SUCCESS;
}


int
caps_makeVal(enum capsvType type, int len, const void *data, capsValue **val)
{
  char             *chars;
  int              *ints, j, slen;
  double           *reals;
  enum capsBoolean *bools;
  capsValue        *value;
  capsTuple        *tuple;
  
  *val  = NULL;
  if (data  == NULL) return CAPS_NULLVALUE;
  value = (capsValue *) EG_alloc(sizeof(capsValue));
  if (value == NULL) return EGADS_MALLOC;
  value->length     = len;
  value->type       = type;
  value->nrow       = value->ncol   = value->dim = 0;
  value->index      = value->pIndex = 0;
  value->lfixed     = value->sfixed = Fixed;
  value->nullVal    = NotAllowed;
  value->units      = NULL;
  value->meshWriter = NULL;
  value->link       = NULL;
  value->limits.dlims[0] = value->limits.dlims[1] = 0.0;
  value->linkMethod      = Copy;
  value->gInType         = 0;
  value->partial         = NULL;
  value->nderiv          = 0;
  value->derivs          = NULL;

  if (data == NULL) {
    value->nullVal = IsNull;
    if (type == Boolean) {
      if (value->length <= 1) {
        value->vals.integer = false;
      } else {
        value->vals.integers = (int *) EG_alloc(value->length*sizeof(int));
        if (value->vals.integers == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.integers[j] = false;
      }
    } else if (type == Integer) {
      if (value->length <= 1) {
        value->vals.integer = 0;
      } else {
        value->vals.integers = (int *) EG_alloc(value->length*sizeof(int));
        if (value->vals.integers == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.integers[j] = 0;
      }
    } else if ((type == Double) || (type == DoubleDeriv)) {
      if (value->length <= 1) {
        value->vals.real = 0.0;
      } else {
        value->vals.reals = (double *) EG_alloc(value->length*sizeof(double));
        if (value->vals.reals == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.reals[j] = 0.0;
      }
    } else if (type == String) {
      value->vals.string = (char *) EG_alloc(2*value->length*sizeof(char));
      if (value->vals.string == NULL) {
        EG_free(value);
        return EGADS_MALLOC;
      }
      for (j = 0; j < 2*value->length; j++) value->vals.string[j] = 0;
    } else if (type == Tuple) {
      value->vals.tuple = NULL;
      if (len > 0) {
        j = caps_makeTuple(len, &value->vals.tuple);
        if ((j != CAPS_SUCCESS) || (value->vals.tuple == NULL)) {
          if (value->vals.tuple == NULL) j = CAPS_NULLVALUE;
          EG_free(value);
          return j;
        }
      }
    } else {
      value->vals.AIMptr = NULL;
    }
    
  } else {
    
    if (type == Boolean) {
      bools = (enum capsBoolean *) data;
      if (value->length == 1) {
        value->vals.integer = bools[0];
      } else {
        value->vals.integers = (int *) EG_alloc(value->length*sizeof(int));
        if (value->vals.integers == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.integers[j] = bools[j];
      }
    } else if (type == Integer) {
      ints = (int *) data;
      if (value->length == 1) {
        value->vals.integer = ints[0];
      } else {
        value->vals.integers = (int *) EG_alloc(value->length*sizeof(int));
        if (value->vals.integers == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.integers[j] = ints[j];
      }
    } else if ((type == Double) || (type == DoubleDeriv)) {
      reals = (double *) data;
      if (value->length == 1) {
        value->vals.real = reals[0];
      } else {
        value->vals.reals = (double *) EG_alloc(value->length*sizeof(double));
        if (value->vals.reals == NULL) {
          EG_free(value);
          return EGADS_MALLOC;
        }
        for (j = 0; j < value->length; j++) value->vals.reals[j] = reals[j];
      }
    } else if (type == String) {
      chars = (char *) data;
      for (slen = j = 0; j < value->length; j++)
        slen += strlen(chars + slen)+1;
      value->vals.string = (char *) EG_alloc(slen*sizeof(char));
      if (value->vals.string == NULL) {
        EG_free(value);
        return EGADS_MALLOC;
      }
      for (j = 0; j < slen; j++) value->vals.string[j] = chars[j];
    } else if (type == Tuple) {
      value->vals.tuple = NULL;
      if (len > 0) {
        j = caps_makeTuple(len, &value->vals.tuple);
        if ((j != CAPS_SUCCESS) || (value->vals.tuple == NULL)) {
          if (value->vals.tuple == NULL) j = CAPS_NULLVALUE;
          EG_free(value);
          return j;
        }
        tuple = (capsTuple *) data;
        for (j = 0; j < len; j++) {
          value->vals.tuple[j].name  = EG_strdup(tuple[j].name);
          value->vals.tuple[j].value = EG_strdup(tuple[j].value);
          if ((tuple[j].name  != NULL) &&
              (value->vals.tuple[j].name  == NULL)) {
            EG_free(value);
            return EGADS_MALLOC;
          }
          if ((tuple[j].value != NULL) &&
              (value->vals.tuple[j].value == NULL)) {
            EG_free(value);
            return EGADS_MALLOC;
          }
        }
      }
    } else {
      value->vals.AIMptr = (void *) data;
    }
  }
  if (value->length > 1) value->dim = Vector;
    
  *val = value;
  return CAPS_SUCCESS;
}


int
caps_findProblem(const capsObject *object, int funID, capsObject **pobject)
{
  capsObject  *pobj;
  capsProblem *problem;
  
  *pobject = NULL;
  if (object == NULL) return CAPS_NULLOBJ;

  pobj = (capsObject *) object;
  do {
    if (pobj->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (pobj->type        == PROBLEM) {
      if (pobj->blind == NULL) return CAPS_NULLBLIND;
      if (funID != 9999) {
        problem  = (capsProblem *) pobj->blind;
        problem->funID = funID;
      }
      *pobject = pobj;
      return CAPS_SUCCESS;
    }
    pobj = pobj->parent;
  } while (pobj != NULL);

  return CAPS_NOTPROBLEM;
}


void
caps_makeSimpleErr(/*@null@*/ capsObject *object, enum capseType type,
                   const char *line1, /*@null@*/ const char *line2,
                   /*@null@*/ const char *line3, capsErrs **errs)
{
  int       i, index = 0;
  capsError *tmp;
  capsErrs  *error;
  
  if (*errs == NULL) {
    error = (capsErrs *) EG_alloc(sizeof(capsErrs));
    if (error == NULL) return;
    error->nError = 1;
    error->errors = (capsError *) EG_alloc(sizeof(capsError));
    if (error->errors == NULL) {
      EG_free(error);
      return;
    }
  } else {
    error = *errs;
    i     = error->nError;
    tmp   = (capsError *) EG_reall(error->errors, (i+1)*sizeof(capsError));
    if (tmp == NULL) return;
    error->errors = tmp;
    error->nError = i+1;
    index         = i;
  }
  
  i = 0;
  if (line1 != NULL) i++;
  if (line2 != NULL) i++;
  if (line3 != NULL) i++;
  error->errors[index].errObj = object;
  error->errors[index].index  = 0;
  error->errors[index].eType  = type;
  error->errors[index].nLines = i;
  error->errors[index].lines  = NULL;
  if (i != 0) {
    error->errors[index].lines = (char **) EG_alloc(i*sizeof(char *));
    if (error->errors[index].lines == NULL) return;
    i = 0;
    if (line1 != NULL) {
      error->errors[index].lines[i] = EG_strdup(line1);
      i++;
    }
    if (line2 != NULL) {
      error->errors[index].lines[i] = EG_strdup(line2);
      i++;
    }
    if (line3 != NULL) error->errors[index].lines[i] = EG_strdup(line3);
  }
  
  *errs = error;
}


int
caps_addHistory(capsObject *object, capsProblem *problem)
{
  int        i = 0;
  char       *pname, *pID, *user;
  capsOwn    *tmp;
  capsObject *pobj;

  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  pobj = problem->mySelf;

  caps_fillDateTime(object->last.datetime);
  
  /* get intent Phrase -- if none bail */
  if (problem->iPhrase < 0) return CAPS_SUCCESS;
  
  /* is it the same as the last -- if so just update the DateTime & sNum */
  if ((object->history != NULL) && (object->nHistory > 0) &&
      (problem->nPhrase > 0)) {
    if  (problem->iPhrase != object->history[object->nHistory-1].index) i = 1;
    if ((object->history[object->nHistory-1].pname != NULL) &&
        (object->last.pname != NULL))
      if (strcmp(object->history[object->nHistory-1].pname,
                 object->last.pname) != 0) i = 1;
    if ((object->history[object->nHistory-1].pID != NULL) &&
        (object->last.pID != NULL))
      if (strcmp(object->history[object->nHistory-1].pID,
                 object->last.pID) != 0) i = 1;
    if ((object->history[object->nHistory-1].user != NULL) &&
        (object->last.user != NULL))
      if (strcmp(object->history[object->nHistory-1].user,
                 object->last.user) != 0) i = 1;
    
    if (i == 0) {
      for (i = 0; i < 6; i++)
        object->history[object->nHistory-1].datetime[i] = object->last.datetime[i];
      object->history[object->nHistory-1].sNum          = object->last.sNum;
      return CAPS_SUCCESS;
    }
  }

  /* make the room */
  if (object->history == NULL) {
    object->nHistory = 0;
    object->history  = (capsOwn *) EG_alloc(sizeof(capsOwn));
    if (object->history == NULL) return EGADS_MALLOC;
  } else {
    tmp = (capsOwn *) EG_reall( object->history,
                               (object->nHistory+1)*sizeof(capsOwn));
    if (tmp == NULL) return EGADS_MALLOC;
    object->history = tmp;
  }
  
  pname = object->last.pname;
  pID   = object->last.pID;
  user  = object->last.user;
  if (pname == NULL) {
    pname = pobj->last.pname;
    pID   = pobj->last.pID;
    user  = pobj->last.user;
  }

  object->history[object->nHistory]       = object->last;
  object->history[object->nHistory].index = problem->iPhrase;
  object->history[object->nHistory].pname = EG_strdup(pname);
  object->history[object->nHistory].pID   = EG_strdup(pID);
  object->history[object->nHistory].user  = EG_strdup(user);
  object->nHistory += 1;

  return CAPS_SUCCESS;
}


int
caps_freeBound(capsObject *object)
{
  int           i, j, status;
  char          filename[PATH_MAX], temp[PATH_MAX];
  capsObject    *pobject;
  capsProblem   *problem;
  capsBound     *bound;
  capsAnalysis  *analysis;
  capsVertexSet *vertexSet;
  capsDataSet   *dataSet;
  FILE          *fp;
  
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->type        != BOUND)     return CAPS_BADTYPE;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;

  status = caps_findProblem(object, 9999, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  caps_freeFList(object);

  bound = (capsBound *) object->blind;
  if (bound->curve   != NULL) caps_Aprx1DFree(bound->curve);
  if (bound->surface != NULL) caps_Aprx2DFree(bound->surface);
  if (bound->lunits  != NULL) EG_free(bound->lunits);
  for (i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;

    vertexSet = (capsVertexSet *) bound->vertexSet[i]->blind;
    if ((vertexSet->analysis != NULL) && (vertexSet->discr != NULL))
      if (vertexSet->analysis->blind != NULL) {
        analysis = (capsAnalysis *) vertexSet->analysis->blind;
        aim_FreeDiscr(problem->aimFPTR, analysis->loadName, vertexSet->discr);
        EG_free(vertexSet->discr);
      }
    for (j = 0; j < vertexSet->nDataSets; j++) {
      if (vertexSet->dataSets[j]->magicnumber != CAPSMAGIC) continue;
      if (vertexSet->dataSets[j]->blind       == NULL)      continue;

      dataSet = (capsDataSet *) vertexSet->dataSets[j]->blind;
      if (dataSet->data    != NULL) EG_free(dataSet->data);
      if (dataSet->units   != NULL) EG_free(dataSet->units);
      if (dataSet->startup != NULL) EG_free(dataSet->startup);
      EG_free(dataSet);
      
      caps_freeHistory(vertexSet->dataSets[j]);
      caps_freeAttrs(&vertexSet->dataSets[j]->attrs);
      caps_freeOwner(&vertexSet->dataSets[j]->last);
      vertexSet->dataSets[j]->magicnumber = 0;
      EG_free(vertexSet->dataSets[j]->name);
      vertexSet->dataSets[j]->name = NULL;
      EG_free(vertexSet->dataSets[j]);
    }
    EG_free(vertexSet->dataSets);
    vertexSet->dataSets = NULL;
    EG_free(vertexSet);

    caps_freeHistory(bound->vertexSet[i]);
    caps_freeAttrs(&bound->vertexSet[i]->attrs);
    caps_freeOwner(&bound->vertexSet[i]->last);
    bound->vertexSet[i]->magicnumber = 0;
    EG_free(bound->vertexSet[i]->name);
    bound->vertexSet[i]->name = NULL;
    EG_free(bound->vertexSet[i]);
  }
  EG_free(bound->vertexSet);
  bound->vertexSet = NULL;
  EG_free(bound);

  /* remove the bound from the list of bounds in the problem */
  for (i = j = 0; i < problem->nBound; i++) {
    if (problem->bounds[i] == object) continue;
    problem->bounds[j++] = problem->bounds[i];
  }
  problem->nBound--;
  
  if (problem->funID != CAPS_CLOSE) {
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\bound.txt", problem->root);
    snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx",  problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/bound.txt",   problem->root);
    snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",    problem->root);
#endif
    fp = fopen(temp, "w");
    if (fp == NULL) {
      printf(" CAPS Warning: Cannot open %s (caps_freeBound)\n", filename);
    } else {
      fprintf(fp, "%d %d\n", problem->nBound, problem->mBound);
      if (problem->bounds != NULL)
        for (i = 0; i < problem->nBound; i++) {
          bound = (capsBound *) problem->bounds[i]->blind;
          j = 0;
          if (bound != NULL) j = bound->index;
          fprintf(fp, "%d %s\n", j, problem->bounds[i]->name);
        }
      fclose(fp);
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS)
        printf(" CAPS Warning: Cannot rename %s!\n", filename);
    }
  }
  
  /* cleanup and invalidate the object */
  
  caps_freeHistory(object);
  caps_freeAttrs(&object->attrs);
  caps_freeOwner(&object->last);
  object->magicnumber = 0;
  EG_free(object->name);
  object->name = NULL;
  EG_free(object);
  
  return CAPS_SUCCESS;
}


int
caps_writeSerialNum(const capsProblem *problem)
{
  char    filename[PATH_MAX];
  size_t  n;
  FILE    *fp;
#ifdef WIN32
  __int64 fpos;

  snprintf(filename, PATH_MAX, "%s\\capsRestart\\Problem",  problem->root);
#else
  long    fpos;
  
  snprintf(filename, PATH_MAX, "%s/capsRestart/Problem",    problem->root);
#endif

  fp = fopen(filename, "ab");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s (caps_writeSerialNum)!\n", filename);
    return CAPS_DIRERR;
  }
  fpos = 0;
#ifdef WIN32
  _fseeki64(fp, fpos, SEEK_SET);
#else
  fseek(fp,     fpos, SEEK_SET);
#endif

  n = fwrite(&problem->sNum, sizeof(CAPSLONG), 1, fp);
  fclose(fp);
  if (n != 1) return CAPS_IOERR;
  
  return CAPS_SUCCESS;
}


/************************* CAPS exposed functions ****************************/


void
caps_revision(int *major, int *minor)
{
  *major = CAPSMAJOR;
  *minor = CAPSMINOR;
}


int
caps_info(capsObject *object, char **name, enum capsoType *type,
          enum capssType *subtype, capsObject **link, capsObject **parent,
          capsOwn *last)
{
  int         i, stat, ret;
  CAPSLONG    sNum;
  capsValue   *value;
  capsProblem *problem;
  capsObject  *pobj;
  capsJrnl    args[2];
  
  *name       = NULL;
  *type       = UNUSED;
  *subtype    = NONE;
  *link       = *parent     = NULL;
  last->index = -1;
  last->user  = last->pname = last->pID = NULL;
  for (i = 0; i < 6; i++) last->datetime[i] = 0;

  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  stat     = caps_findProblem(object, CAPS_INFO, &pobj);
  if (stat != CAPS_SUCCESS) return stat;
  problem  = (capsProblem *) pobj->blind;
  
  *type    = object->type;
  *subtype = object->subtype;
  *name    = object->name;
  *parent  = object->parent;
  
  args[0].type = jObject;
  args[1].type = jOwn;
  if (problem->dbFlag == 0) {
    stat       = caps_jrnlRead(CAPS_INFO, problem, object, 2, args, &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *link = args[0].members.obj;
        *last = args[1].members.own;
      }
      return ret;
    }
  }
  
  if (object->type == VALUE) {
    value  = (capsValue *) object->blind;
    *link  = value->link;
  }
  if (object->last.pname == NULL) {
    if (stat == CAPS_SUCCESS) {
      last->pname = pobj->last.pname;
      last->user  = pobj->last.user;
      last->pID   = pobj->last.pID;
    }
  } else {
    last->pname = object->last.pname;
    last->user  = object->last.user;
    last->pID   = object->last.pID;
  }
  last->index = object->last.index;
  last->sNum  = object->last.sNum;
  for (i = 0; i < 6; i++) last->datetime[i] = object->last.datetime[i];
  if (problem->dbFlag == 1) return CAPS_SUCCESS;
  
  args[0].members.obj = *link;
/*@-nullret@*/
  args[1].members.own = *last;
/*@+nullret@*/
  caps_jrnlWrite(CAPS_INFO, problem, object, CAPS_SUCCESS, 2, args,
                 problem->sNum, problem->sNum);
  
  return CAPS_SUCCESS;
}


static int
caps_sizeX(capsObject *object, enum capsoType type, enum capssType stype,
           int *size, int *nErr, capsErrs **errors)
{
  int          i, status, major, minor, nField, exec, dirty;
  int          *ranks, *fInOut;
  char         *intents, *apath, *unitSys, **fnames;
  egAttrs       *attrs;
  capsObject    *pobject;
  capsProblem   *problem;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexSet;

  *nErr   = 0;
  *errors = NULL;
  *size   = 0;
  
  if (object->type == PROBLEM) {
  
    problem = (capsProblem *) object->blind;
    if (type == BODIES) {
      /* make sure geometry is up-to-date */
      if (problem->dbFlag == 0) {
        status = caps_build(object, nErr, errors);
        if ((status != CAPS_SUCCESS) && (status != CAPS_CLEAN)) return status;
      }
      *size = problem->nBodies;
    } else if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    } else if (type == VALUE) {
      if (stype == GEOMETRYIN)  *size = problem->nGeomIn;
      if (stype == GEOMETRYOUT) *size = problem->nGeomOut;
      if (stype == PARAMETER)   *size = problem->nParam;
      if (stype == USER)        *size = problem->nUser;
    } else if (type == ANALYSIS) {
      *size = problem->nAnalysis;
    } else if (type == BOUND) {
      *size = problem->nBound;
    }
    
  } else if (object->type == VALUE) {

    if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    }
    
  } else if (object->type == ANALYSIS) {

    status = caps_findProblem(object, CAPS_SIZE, &pobject);
    if (status != CAPS_SUCCESS) return status;
    analysis = (capsAnalysis *) object->blind;
    if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    } else if (type == VALUE) {
      problem = (capsProblem *) pobject->blind;
      if (problem == NULL) return CAPS_NULLBLIND;
      if (stype == ANALYSISIN)   *size = analysis->nAnalysisIn;
      if (stype == ANALYSISOUT)  *size = analysis->nAnalysisOut;
      if (stype == ANALYSISDYNO) {
        if (problem->dbFlag == 0) {
          /* check to see if analysis is CLEAN, and auto-execute if possible */
          status = caps_analysisInfX(object, &apath, &unitSys, &major, &minor,
                                     &intents, &nField, &fnames, &ranks,
                                     &fInOut, &exec, &dirty);
          if (status != CAPS_SUCCESS) return status;
          if (dirty > 0) {
            /* auto execute if available */
            if ((exec == 2) && (dirty < 5)) {
              status = caps_execX(object, nErr, errors);
              if (status != CAPS_SUCCESS) return status;
            } else {
              return CAPS_DIRTY;
            }
          }
        }
        *size = analysis->nAnalysisDynO;
      }
    } else if (type == BODIES) {

      problem = (capsProblem *) pobject->blind;
      if (problem == NULL) return CAPS_NULLBLIND;
      if (problem->dbFlag == 1) {
        *size = analysis->nBody;
        return CAPS_SUCCESS;
      }
      
      /* make sure geometry is up-to-date */
      status = caps_build(pobject, nErr, errors);
      if ((status != CAPS_SUCCESS) && (status != CAPS_CLEAN)) return status;

      if ((problem->nBodies > 0) && (problem->bodies != NULL)) {
        if (analysis->bodies == NULL) {
          status = caps_filter(problem, analysis);
          if (status != CAPS_SUCCESS) return status;
        }
        *size = analysis->nBody;
      }
    }
    
  } else if (object->type == BOUND) {
    
    if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    } else if (type == VERTEXSET) {
      bound = (capsBound *) object->blind;
      for (*size = i = 0; i < bound->nVertexSet; i++)
        if (bound->vertexSet[i]->subtype == stype) *size += 1;
    }
    
  } else if (object->type == VERTEXSET) {
    
    if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    } else if (type == DATASET) {
      vertexSet = (capsVertexSet *) object->blind;
      *size = vertexSet->nDataSets;
    }
    
  } else if (object->type == DATASET) {

    if (type == ATTRIBUTES) {
      attrs = object->attrs;
      if (attrs != NULL) *size = attrs->nattrs;
    }

  }

  return CAPS_SUCCESS;
}


int
caps_size(capsObject *object, enum capsoType type, enum capssType stype,
          int *size, int *nErr, capsErrs **errors)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[3];

  *nErr   = 0;
  *errors = NULL;
  *size   = 0;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  status = caps_findProblem(object, CAPS_SIZE, &pobject);
  if (status != CAPS_SUCCESS)           return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)
    return caps_sizeX(object, type, stype, size, nErr, errors);
  
  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jErr;
  status       = caps_jrnlRead(CAPS_SIZE, problem, object, 3, args, &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    *size   = args[0].members.integer;
    *nErr   = args[1].members.integer;
    *errors = args[2].members.errs;
    return ret;
  }
  
  sNum = problem->sNum;
  ret  = caps_sizeX(object, type, stype, size, nErr, errors);
  if (ret == CAPS_SUCCESS) {
    args[0].members.integer = *size;
    args[1].members.integer = *nErr;
    args[2].members.errs    = *errors;
  }
  caps_jrnlWrite(CAPS_SIZE, problem, object, ret, 3, args, sNum, problem->sNum);
  
  return ret;
}


static int
caps_childByIndX(capsObject *object, capsProblem *problem, enum capsoType type,
                 enum capssType stype, int index, capsObject **child)
{
  int           i, j;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexSet;
  
  *child = NULL;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  if (index               <= 0)         return CAPS_RANGEERR;
  
  if (object->type == PROBLEM) {

    if (type == VALUE) {
      if (stype == GEOMETRYIN) {
        if (index > problem->nGeomIn)   return CAPS_RANGEERR;
        *child = problem->geomIn[index-1];
      }
      if (stype == GEOMETRYOUT) {
        if (index > problem->nGeomOut)  return CAPS_RANGEERR;
        *child = problem->geomOut[index-1];
      }
      if (stype == PARAMETER) {
        if (index > problem->nParam)    return CAPS_RANGEERR;
        *child = problem->params[index-1];
      }
      if (stype == USER) {
        if (index > problem->nUser)     return CAPS_RANGEERR;
        *child = problem->users[index-1];
      }
    } else if (type == ANALYSIS) {
      if (index > problem->nAnalysis)   return CAPS_RANGEERR;
      *child = problem->analysis[index-1];
    } else if (type == BOUND) {
      if (index  > problem->nBound)     return CAPS_RANGEERR;
      *child = problem->bounds[index-1];
    }
    
  } else if (object->type == ANALYSIS) {
    
    analysis = (capsAnalysis *) object->blind;
    if (type == VALUE) {
      if (stype == ANALYSISIN) {
        if (index > analysis->nAnalysisIn)   return CAPS_RANGEERR;
        *child = analysis->analysisIn[index-1];
      }
      if (stype == ANALYSISOUT) {
        if (index > analysis->nAnalysisOut)  return CAPS_RANGEERR;
        *child = analysis->analysisOut[index-1];
      }
      if (stype == ANALYSISDYNO) {
        if (index > analysis->nAnalysisDynO) return CAPS_RANGEERR;
        *child = analysis->analysisDynO[index-1];
      }
    }
    
  } else if (object->type == BOUND) {
    
    bound = (capsBound *) object->blind;
    for (j = i = 0; i < bound->nVertexSet; i++) {
      if (bound->vertexSet[i]->subtype == stype) j++;
      if (j != index) continue;
      *child = bound->vertexSet[i];
      break;
    }
    
  } else if (object->type == VERTEXSET) {
    
    vertexSet = (capsVertexSet *) object->blind;
    if (type == DATASET) {
      if (index > vertexSet->nDataSets) return CAPS_RANGEERR;
      *child = vertexSet->dataSets[index-1];
    }
    
  }
  
  if (*child == NULL) return CAPS_NOTFOUND;
  return CAPS_SUCCESS;
}


int
caps_childByIndex(capsObject *object, enum capsoType type, enum capssType stype,
                  int index, capsObject **child)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[1];
  
  *child = NULL;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  if (index               <= 0)         return CAPS_RANGEERR;

  status = caps_findProblem(object, CAPS_CHILDBYINDEX, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)
    return caps_childByIndX(object, problem, type, stype, index, child);
  
  args[0].type = jObject;
  status       = caps_jrnlRead(CAPS_CHILDBYINDEX, problem, object, 1, args,
                               &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    *child = args[0].members.obj;
    return ret;
  }
  
  ret = caps_childByIndX(object, problem, type, stype, index, child);
  args[0].members.obj = *child;
  caps_jrnlWrite(CAPS_CHILDBYINDEX, problem, object, ret, 1, args,
                 problem->sNum, problem->sNum);
  
  return ret;
}


static int
caps_findByName(const char *name, int len, capsObject **objs,
                capsObject **child)
{
  int i;
  
  if (objs == NULL) return CAPS_NOTFOUND;
  
  for (i = 0; i < len; i++) {
    if (objs[i]       == NULL) continue;
    if (objs[i]->name == NULL) continue;
    if (strcmp(objs[i]->name, name) == 0) {
      *child = objs[i];
      return CAPS_SUCCESS;
    }
  }
  
  return CAPS_NOTFOUND;
}


int
caps_childByName(capsObject *object, enum capsoType type, enum capssType stype,
                 const char *name, capsObject **child,
                 int *nErr, capsErrs **errors)
{
  int           stat, ret, major, minor, nField, exec, dirty;
  int           *ranks, *fInOut;
  char          *intents, *apath, *unitSys, **fnames;
  CAPSLONG      sNum;
  capsObject    *pobject;
  capsProblem   *problem;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexSet;
  capsJrnl      args[1];
  
  if (child               == NULL)      return CAPS_NULLOBJ;
  if (nErr                == NULL)      return CAPS_NULLOBJ;
  if (errors              == NULL)      return CAPS_NULLOBJ;
  *child = NULL;
  *nErr = 0;
  *errors = NULL;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  if (name                == NULL)      return CAPS_NULLNAME;
  
  stat = caps_findProblem(object, CAPS_CHILDBYNAME, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;
  
  args[0].type = jObject;
  if (problem->dbFlag == 0) {
    stat       = caps_jrnlRead(CAPS_CHILDBYNAME, problem, object, 1, args,
                               &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *child = args[0].members.obj;
      return ret;
    }
  }
  
  if (object->type == PROBLEM) {
    
    if (type == VALUE) {
      if (stype == GEOMETRYIN) {
        ret = caps_findByName(name, problem->nGeomIn,
                                    problem->geomIn,  child);
        goto retrn;
      }
      if (stype == GEOMETRYOUT) {
        ret = caps_findByName(name, problem->nGeomOut,
                                    problem->geomOut, child);
        goto retrn;
      }
      if (stype == PARAMETER) {
        ret = caps_findByName(name, problem->nParam,
                                    problem->params,  child);
        goto retrn;
      }
    } else if (type == ANALYSIS) {
      ret = caps_findByName(name, problem->nAnalysis,
                                  problem->analysis,  child);
      goto retrn;
    } else if (type == BOUND) {
      ret = caps_findByName(name, problem->nBound,
                                  problem->bounds,    child);
      goto retrn;
    }
    
  } else if (object->type == ANALYSIS) {
    
    analysis = (capsAnalysis *) object->blind;
    if (type == VALUE) {
      if (stype == ANALYSISIN) {
        ret = caps_findByName(name, analysis->nAnalysisIn,
                                    analysis->analysisIn,   child);
        goto retrn;
      }
      if (stype == ANALYSISOUT) {
        ret = caps_findByName(name, analysis->nAnalysisOut,
                                    analysis->analysisOut,  child);
        goto retrn;
      }
      if (stype == ANALYSISDYNO) {

        if (problem->dbFlag == 0) {
          /* check to see if analysis is CLEAN, and auto-execute if possible */
          stat = caps_analysisInfX(object, &apath, &unitSys, &major, &minor,
                                   &intents, &nField, &fnames, &ranks, &fInOut,
                                   &exec, &dirty);
          if (stat != CAPS_SUCCESS) return stat;
          if (dirty > 0) {
            /* auto execute if available */
            if ((exec == 2) && (dirty < 5)) {
              stat = caps_execX(object, nErr, errors);
              if (stat != CAPS_SUCCESS) return stat;
            } else {
              return CAPS_DIRTY;
            }
          }
        }

        ret = caps_findByName(name, analysis->nAnalysisDynO,
                                    analysis->analysisDynO, child);
        goto retrn;
      }
    }
    
  } else if (object->type == BOUND) {
    
    bound = (capsBound *) object->blind;
    if ((type == VERTEXSET) && (stype == CONNECTED)) {
      ret = caps_findByName(name, bound->nVertexSet,
                                  bound->vertexSet, child);
      goto retrn;
    }
    
  } else if (object->type == VERTEXSET) {
    
    vertexSet = (capsVertexSet *) object->blind;
    if (type == DATASET) {
      ret = caps_findByName(name, vertexSet->nDataSets,
                                  vertexSet->dataSets, child);
      goto retrn;
    }
    
  }
  ret = CAPS_NOTFOUND;
  if (problem->dbFlag == 1) return ret;

retrn:
  args[0].members.obj = *child;
  caps_jrnlWrite(CAPS_CHILDBYNAME, problem, object, ret, 1, args, problem->sNum,
                 problem->sNum);
  
  return ret;
}


int
caps_bodyByIndex(capsObject *object, int index, ego *body, char **lunits)
{
  int          i, status, ret, oclass, mtype, *senses;
  CAPSLONG     sNum;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsJrnl     args[2];
  ego          ref, *bodies;
  
  *body   = NULL;
  *lunits = NULL;
  if  (object              == NULL)      return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if ((object->type        != PROBLEM) &&
      (object->type        != ANALYSIS)) return CAPS_BADTYPE;
  if  (object->blind       == NULL)      return CAPS_NULLBLIND;

  status = caps_findProblem(object, CAPS_BODYBYINDEX, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  
  args[0].type = jEgos;
  args[1].type = jString;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_BODYBYINDEX, problem, object, 2, args,
                                 &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      status = EG_getTopology(args[0].members.model, &ref, &oclass, &mtype, NULL,
                              &i, &bodies, &senses);
      if (status != EGADS_SUCCESS) {
        printf(" CAPS Warning: EG_getTopology = %d (caps_bodyByIndex)\n", status);
      } else {
        *body = bodies[0];
      }
      *lunits = args[1].members.string;
      return ret;
    }
  }
  
  status = CAPS_RANGEERR;
  if (object->type == PROBLEM) {
    
    if (index               <= 0) goto done;
    if (index > problem->nBodies) goto done;
    *body   = problem->bodies[index-1];
    *lunits = problem->lunits[index-1];
    
  } else {
    
    if (index             <= 0) goto done;
    analysis = (capsAnalysis *) object->blind;
    if ((problem->nBodies > 0) && (problem->bodies != NULL)) {
      if (analysis->bodies == NULL) {
        status = caps_filter(problem, analysis);
        if (status != CAPS_SUCCESS) goto done;
        status = CAPS_RANGEERR;
      }
      if (index > analysis->nBody) goto done;
      *body = analysis->bodies[index-1];
      for (i = 0; i < problem->nBodies; i++)
        if (*body == problem->bodies[i]) {
          *lunits = problem->lunits[i];
          break;
        }
    }
    
  }
  
  status = CAPS_SUCCESS;
  if (problem->dbFlag == 1) return status;
  
done:
  args[0].members.model  = *body;
  args[1].members.string = *lunits;
  caps_jrnlWrite(CAPS_BODYBYINDEX, problem, object, status, 2, args,
                 problem->sNum, problem->sNum);
  
  return status;
}


int
caps_ownerInfo(const capsObject *pobject, const capsOwn owner, char **phase,
               char **pname, char **pID, char **userID, int *nLines,
               char ***lines, short *datetime, CAPSLONG *sNum)
{
  int         i;
  capsProblem *problem;
  
  *phase  = NULL;
  *pname  = NULL;
  *pID    = NULL;
  *userID = NULL;
  *nLines = 0;
  *lines  = 0;
  *sNum   = 0;
  for (i = 0; i < 6; i++) datetime[i] = 0;
  
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if ((owner.index >= 0) && (owner.index < problem->nPhrase)) {
    *phase  = problem->phrases[owner.index].phase;
    *nLines = problem->phrases[owner.index].nLines;
    *lines  = problem->phrases[owner.index].lines;
  }
  *pname  = owner.pname;
  *pID    = owner.pID;
  *userID = owner.user;
  *sNum   = owner.sNum;
  for (i = 0; i < 6; i++) datetime[i] = owner.datetime[i];
  
  return CAPS_SUCCESS;
}


int
caps_delete(capsObject *object)
{
  int         ret, stat;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl     args[1];    /* not really used */

  if  (object              == NULL)       return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if ((object->type        != ANALYSIS) && (object->type != VALUE) &&
      (object->type        != BOUND))     return CAPS_BADTYPE;
  if ((object->type        == VALUE) &&
      (object->subtype     != PARAMETER)) return CAPS_BADTYPE;
  if  (object->blind       == NULL)       return CAPS_NULLBLIND;
  if  (object->parent      == NULL)       return CAPS_NULLOBJ;
  pobject = (capsObject *)   object->parent;
  if  (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)          return CAPS_READONLYERR;

  args[0].type = jString;
  stat         = caps_jrnlRead(CAPS_DELETE, problem, object, 0, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL)    return ret;

  object->delMark = 1;
  /* write the object */
  if (object->type == ANALYSIS) {
    ret = caps_writeAnalysisObj(problem, object);
  } else if (object->type == BOUND) {
    ret = caps_writeBound(object);
  } else {
    ret = caps_writeValueObj(problem, object);
  }
  
  caps_jrnlWrite(CAPS_DELETE, problem, object, ret, 0, args, problem->sNum,
                 problem->sNum);

  return ret;
}


int
caps_errorInfo(capsErrs *errs, int eIndex, capsObject **errObj, int *eType,
               int *nLines, char ***lines)
{
  *errObj = NULL;
  *nLines = 0;
  *lines  = NULL;
  if ((eIndex < 1) || (eIndex > errs->nError)) return CAPS_BADINDEX;
  
  *errObj = errs->errors[eIndex-1].errObj;
  *eType  = errs->errors[eIndex-1].eType;
  *nLines = errs->errors[eIndex-1].nLines;
  *lines  = errs->errors[eIndex-1].lines;
  return CAPS_SUCCESS;
}


int
caps_freeError(/*@only@*/ capsErrs *errs)
{
  int i, j;
  
  if (errs == NULL) return CAPS_SUCCESS;
  for (i = 0; i < errs->nError; i++) {
    for (j = 0; j < errs->errors[i].nLines; j++) {
      EG_free(errs->errors[i].lines[j]);
    }

    EG_free(errs->errors[i].lines);
  }
  
  EG_free(errs->errors);
  errs->nError = 0;
  errs->errors = NULL;
  EG_free(errs);

  return CAPS_SUCCESS;
}


int
caps_printErrors(/*@null@*/ FILE *fp, int nErr, capsErrs *errors)
{
  int         i, j, stat, eType, nLines;
  char        **lines;
  capsObj     obj;
  static char *type[5] = {"Cont:   ", "Info:   ", "Warning:", "Error:  ",
                          "Status: "};

  if (errors == NULL) return CAPS_SUCCESS;

  for (i = 1; i <= nErr; i++) {
    stat = caps_errorInfo(errors, i, &obj, &eType, &nLines, &lines);
    if (stat != CAPS_SUCCESS) {
      if (fp == NULL) {
        printf(" printErrors: %d/%d caps_errorInfo = %d\n", i, nErr, stat);
      } else {
        fprintf(fp, " printErrors: %d/%d caps_errorInfo = %d\n", i, nErr, stat);
      }
      caps_freeError(errors);
      return stat;
    }
    for (j = 0; j < nLines; j++)
      if (fp == NULL) {
        if (j == 0) {
          printf(" CAPS %s ", type[eType+1]);
        } else {
          printf("               ");
        }
        printf("%s\n", lines[j]);
      } else {
        if (j == 0) {
          fprintf(fp, " CAPS %s ", type[eType+1]);
        } else {
          fprintf(fp, "               ");
        }
        fprintf(fp, "%s\n", lines[j]);
      }

  }
  
  caps_freeError(errors);
  return CAPS_SUCCESS;
}


int
caps_writeParameters(const capsObject *pobject, char *filename)
{
  int         status;
  capsProblem *problem;
  
  if (pobject              == NULL)       return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)    return CAPS_BADTYPE;
  if (pobject->subtype     != PARAMETRIC) return CAPS_BADTYPE;
  if (pobject->blind       == NULL)       return CAPS_NULLBLIND;
  if (filename             == NULL)       return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_WRITEPARAMETERS;
  
  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }
  
  return ocsmSaveDespmtrs(problem->modl, filename);
}


int
caps_readParameters(const capsObject *pobject, char *filename)
{
  int         i, j, k, n, type, nrow, ncol, fill, status;
  char        name[MAX_NAME_LEN];
  double      real, dot, *reals;
  capsObject  *object;
  capsProblem *problem;
  capsValue   *value;
  
  if (pobject              == NULL)       return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)    return CAPS_BADTYPE;
  if (pobject->subtype     != PARAMETRIC) return CAPS_BADTYPE;
  if (pobject->blind       == NULL)       return CAPS_NULLBLIND;
  if (filename             == NULL)       return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)               return CAPS_READONLYERR;
  problem->funID = CAPS_READPARAMETERS;
  
  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }
  
  status  = ocsmUpdateDespmtrs(problem->modl, filename);
  if (status < SUCCESS) {
    printf(" CAPS Error: ocsmSaveDespmtrs = %d (caps_readParameters)!\n",
           status);
    return status;
  }
  
  /* need to reload all GeomIn Values */

  problem->sNum += 1;
  for (i = 0; i < problem->nGeomIn; i++) {
    object = problem->geomIn[i];
    if (object              == NULL)      continue;
    if (object->magicnumber != CAPSMAGIC) continue;
    if (object->type        != VALUE)     continue;
    if (object->blind       == NULL)      continue;
    value = (capsValue *) object->blind;
    if (value               == NULL)      continue;
    if (value->type         != Double)    continue;
    status = ocsmGetPmtr(problem->modl, value->pIndex, &type,
                         &nrow, &ncol, name);
    if (status              != SUCCESS)   continue;
    fill = 0;

    /* has the shape changed? */
    if ((nrow != value->nrow) || (ncol != value->ncol)) {
      reals = NULL;
      if (nrow*ncol != 1) {
        reals = (double *) EG_alloc(nrow*ncol*sizeof(double));
        if (reals == NULL) {
          printf(" CAPS Warning: %s resize %d %d Malloc(caps_readParameters)\n",
                 object->name, nrow, ncol);
          continue;
        }
      }
      if (value->length != 1) EG_free(value->vals.reals);
      value->length = nrow*ncol;
      value->nrow   = nrow;
      value->ncol   = ncol;
      if (value->length != 1) value->vals.reals = reals;
      fill = 1;
    }
    
    /* check if values changed */
    if (fill == 0) {
      reals = value->vals.reals;
      if (value->length == 1) reals = &value->vals.real;
      for (n = k = 0; k < nrow; k++) {
        for (j = 0; j < ncol; j++, n++) {
          status = ocsmGetValu(problem->modl, value->pIndex, k+1, j+1,
                               &real, &dot);
          if (status != SUCCESS) {
            printf(" CAPS Warning: %s GetValu[%d,%d] = %d (caps_readParameters)\n",
                   object->name, k+1, j+1, status);
            continue;
          }
          if (real != reals[n]) {
            fill = 1;
            break;
          }
        }
        if (fill == 1) break;
      }
    }
    
    /* refill the values */
    if (fill == 0) continue;
    reals = value->vals.reals;
    if (value->length == 1) reals = &value->vals.real;
    for (n = k = 0; k < nrow; k++)
      for (j = 0; j < ncol; j++, n++)
        ocsmGetValu(problem->modl, value->pIndex, k+1, j+1, &reals[n], &dot);
    
    caps_freeOwner(&object->last);
    object->last.sNum = problem->sNum;
    caps_fillDateTime(object->last.datetime);
  }
  status = caps_dumpGeomVals(problem, 1);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_dumpGeomVals = %d (caps_readParameters)\n",
           status);
  status = caps_writeSerialNum(problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeSerialNum = %d (caps_readParameters)\n",
           status);
  
  return status;
}


static int
caps_writeGeometrX(capsObject *object, int flag, const char *filename,
                   int *nErr, capsErrs **errors)
{
  int          i, j, nBody, nTess=0, len, idot, stat;
  ego          context, model, *bodies, *tess, *newBodies;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  FILE         *fp;
  
  *nErr    = 0;
  *errors  = NULL;

  stat = caps_findProblem(object, CAPS_WRITEGEOMETRY, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  
  /* find the extension */
  len = strlen(filename);
  for (idot = len-1; idot > 0; idot--)
    if (filename[idot] == '.') break;
  if (idot == 0) return CAPS_BADNAME;
/*@-unrecog@*/
  if ((strcasecmp(&filename[idot],".iges")  != 0) &&
      (strcasecmp(&filename[idot],".igs")   != 0) &&
      (strcasecmp(&filename[idot],".step")  != 0) &&
      (strcasecmp(&filename[idot],".stp")   != 0) &&
      (strcasecmp(&filename[idot],".brep")  != 0) &&
      (strcasecmp(&filename[idot],".egads") != 0)) return CAPS_BADNAME;
/*@+unrecog@*/
  if (strcasecmp(&filename[idot],".egads") == 0)
    if ((flag != 0) && (flag != 1)) return CAPS_RANGEERR;
  
  /* make sure geometry is up-to-date */
  stat = caps_build(pobject, nErr, errors);
  if ((stat != CAPS_SUCCESS) && (stat != CAPS_CLEAN)) return stat;

  if (object->type == PROBLEM) {
    problem  = (capsProblem *) object->blind;
    nBody    = problem->nBodies;
    bodies   = problem->bodies;
    tess     = NULL;
  } else {
    analysis = (capsAnalysis *) object->blind;
    problem  = (capsProblem *) pobject->blind;
    if (analysis->bodies == NULL) {
      stat = caps_filter(problem, analysis);
      if (stat != CAPS_SUCCESS) return stat;
    }
    nBody    = analysis->nBody;
    bodies   = analysis->bodies;
    tess     = analysis->tess;
  }
  context    = problem->context;
  
  if ((nBody <= 0) || (bodies == NULL)) return CAPS_NOBODIES;

  /* remove existing file (if any) */
  fp = fopen(filename, "r");
  if (fp != NULL) {
    fclose(fp);
    caps_rmFile(filename);
  }
  
  if (nBody == 1 && (tess == NULL || flag == 0)) {
    stat =  EG_saveModel(bodies[0], filename);
    if (stat != EGADS_SUCCESS) {
      printf(" CAPS Error: EG_saveModel = %d (caps_writeGeometry)!\n", stat);
      return stat;
    }
  } else {
    newBodies = (ego *) EG_alloc(2*nBody*sizeof(ego));
    if (newBodies == NULL) {
      printf(" CAPS Error: MALLOC on %d Bodies (caps_writeGeometry)!\n", nBody);
      return EGADS_MALLOC;
    }
    for (i = 0; i < nBody; i++) {
      stat = EG_copyObject(bodies[i], NULL, &newBodies[i]);
      if (stat != EGADS_SUCCESS) {
        printf(" CAPS Error: EG_copyObject %d/%d = %d (caps_writeGeometry)!\n",
               i+1, nBody, stat);
        for (j = 0; j < i; j++) EG_deleteObject(newBodies[j]);
        EG_free(newBodies);
        return stat;
      }
    }

    nTess = 0;
    if ((strcasecmp(&filename[idot],".egads") == 0) &&
        (tess != NULL) &&
        (flag == 1)) {
      for (i = nBody; i < 2*nBody; i++) {
        if (tess[i-nBody] == NULL) continue;
        stat = EG_copyObject(tess[i-nBody], newBodies[i-nBody],
                             &newBodies[nBody+nTess++]);
        if (stat != EGADS_SUCCESS) {
          printf(" CAPS Error: EG_copyObject %d/%d = %d (caps_writeGeometry)!\n",
                 i+1, 2*nBody, stat);
          for (j = 0; j < i; j++) EG_deleteObject(newBodies[j]);
          EG_free(newBodies);
          return stat;
        }
      }
    }

    /* make a Model */
    stat = EG_makeTopology(context, NULL, MODEL, nBody+nTess, NULL, nBody,
                           newBodies, NULL, &model);
    EG_free(newBodies);
    if (stat != EGADS_SUCCESS) {
      printf(" CAPS Error: EG_makeTopology %d = %d (caps_writeGeometry)!\n",
             nBody, stat);
      return stat;
    }
    stat = EG_saveModel(model, filename);
    EG_deleteObject(model);
    if (stat != EGADS_SUCCESS) {
      printf(" CAPS Error: EG_saveModel = %d (caps_writeGeometry)!\n", stat);
      return stat;
    }
  }

  return CAPS_SUCCESS;
}


int
caps_writeGeometry(capsObject *object, int flag, const char *filename,
                   int *nErr, capsErrs **errors)
{
  int         stat, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[2];
  
  *nErr    = 0;
  *errors  = NULL;
  if  (object              == NULL)      return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if ((object->type        != PROBLEM) &&
      (object->type        != ANALYSIS)) return CAPS_BADTYPE;
  if  (object->blind       == NULL)      return CAPS_NULLBLIND;
  if  (filename            == NULL)      return CAPS_NULLNAME;

  stat = caps_findProblem(object, CAPS_WRITEGEOMETRY, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1) return CAPS_READONLYERR;
  
  args[0].type = jInteger;
  args[1].type = jErr;
  stat         = caps_jrnlRead(CAPS_WRITEGEOMETRY, problem, object, 2, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *nErr   = args[0].members.integer;
    *errors = args[1].members.errs;
    return ret;
  }
  
  sNum = problem->sNum;
  ret  = caps_writeGeometrX(object, flag, filename, nErr, errors);
  args[0].members.integer = *nErr;
  args[1].members.errs    = *errors;
  caps_jrnlWrite(CAPS_WRITEGEOMETRY, problem, object, ret, 2, args, sNum,
                 problem->sNum);
  
  return ret;
}


int
caps_getHistory(capsObject *object, int *nHist, capsOwn **history)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[1];

  *nHist   = 0;
  *history = NULL;
  if (object              == NULL)         return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (object->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(object, CAPS_GETHISTORY, &pobject);
  if (status              != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jOwns;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_GETHISTORY, problem, object, 1, args,
                               &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *nHist   = args[0].num;
        *history = args[0].members.owns;
      }
      return ret;
    }
  }

  *nHist   = object->nHistory;
  *history = object->history;
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  args[0].num          = *nHist;
  args[0].members.owns = *history;
  caps_jrnlWrite(CAPS_GETHISTORY, problem, object, CAPS_SUCCESS, 1, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}
