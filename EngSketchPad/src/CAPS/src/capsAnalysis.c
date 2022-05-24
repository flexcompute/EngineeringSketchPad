/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Analysis Object Functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
#define snprintf _snprintf
#define strtok_r strtok_s
#define PATH_MAX _MAX_PATH
#else
#include <unistd.h>
#include <limits.h>
#endif

#include "capsBase.h"
#include "capsAIM.h"
#include "prm.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"

#include "egadsTris.h"

#define NOTFILLED        -1
#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


typedef struct {
  ego geom;                     /* geometry object */
  int gIndex;                   /* geometry object index */
  int n;                        /* number of objects */
  ego *objs;                    /* objects */
  int *indices;                 /* markers */
} bodyObjs;


extern /*@null@*/ /*@only@*/
       char *EG_strdup(/*@null@*/ const char *str);
extern void  EG_makeConnect(int k1, int k2, int *tri, int *kedge, int *ntable,
                            connect *etable, int face);
extern int   caps_statFile(const char *path);
extern int   caps_mkDir(const char *path);
extern int   caps_rmDir(const char *path);
extern int   caps_rename(const char *src, const char *dst);
extern int   caps_isNameOK(const char *name);
extern int   caps_dumpAnalysis(capsProblem *problem, capsObject *aobject);
extern int   caps_writeBound(const capsObject *bobject);
extern int   caps_writeVertexSet(capsObject *vobject);
extern int   caps_writeDataSet(capsObject *dobject);
extern void  caps_jrnlWrite(int funID, capsProblem *problem, capsObject *obj,
                            int status, int nargs, capsJrnl *args, CAPSLONG sNm0,
                            CAPSLONG sNum);
extern int   caps_jrnlEnd(capsProblem *problem);
extern int   caps_jrnlRead(int funID, capsProblem *problem, capsObject *obj,
                           int nargs, capsJrnl *args, CAPSLONG *sNum, int *stat);

extern int   caps_dupValues(capsValue *val1, capsValue *val2);
extern int   caps_transferValueX(capsObject *source, enum capstMethod method,
                                 capsObject *trgt, int *nErr, capsErrs **errs);
extern int   caps_snDataSets(const capsObject *aobject, int flag, CAPSLONG *sn);
extern int   caps_fillCoeff2D(int nrank, int nux, int nvx, double *fit,
                              double *coeff, double *r);
extern int   caps_Aprx1DFree(/*@only@*/ capsAprx1D *approx);
extern int   caps_Aprx2DFree(/*@only@*/ capsAprx2D *approx);
extern int   caps_invInterpolate1D(capsAprx1D *interp, double *sv, double *t);
extern int   caps_invInterpolate2D(capsAprx2D *interp, double *sv, double *uv);
extern int   caps_build(capsObject *pobject, int *nErr, capsErrs **errors);
extern int   caps_freeError(/*@only@*/ capsErrs *errs);
extern int   caps_checkValueObj(capsObject *object);
extern int   caps_getDataX(capsObject *dobject, int *npts, int *rank,
                           double **data, char **units, int *nErr,
                           capsErrs **errors);



void
caps_concatErrs(/*@null@*/ capsErrs *errs, capsErrs **errors)
{
  int       i, n;
  capsErrs  *oldErrs;
  capsError *newErr;

  if (errs    == NULL) return;
  if (*errors == NULL) {
    *errors = errs;
    return;
  }

  /* merge the Errors */
  oldErrs = *errors;
  n       = errs->nError + oldErrs->nError;
  newErr  = (capsError *) EG_reall(errs->errors, n*sizeof(capsError));
  if (newErr == NULL) {
    printf(" CAPS Internal: MALLOC on %d Errors (caps_concatErrs)!\n", n);
    caps_freeError(errs);
    return;
  }
  for (i = 0; i < oldErrs->nError; i++)
    newErr[i+errs->nError] = oldErrs->errors[i];
  EG_free(oldErrs->errors);
  oldErrs->errors = newErr;
  oldErrs->nError = n;
  EG_free(errs);
}


void
caps_getAIMerrs(capsAnalysis *analysis, int *nErr, capsErrs **errors)
{
  capsErrs *errs;

  *nErr   = 0;
  *errors = NULL;
  if (analysis->info.errs.nError == 0) return;

  errs    = (capsErrs *) EG_alloc(sizeof(capsErrs));
  if (errs == NULL) {
    printf(" CAPS Warning: Allocation Error (caps_getAIMerrs)\n");
    return;
  }
  *nErr        = analysis->info.errs.nError;
  errs->nError = *nErr;
  errs->errors = analysis->info.errs.errors;
  *errors      = errs;
  analysis->info.errs.nError = 0;
  analysis->info.errs.errors = NULL;
}


int
caps_queryAnalysis(capsObject *pobject, const char *aname, int *nIn, int *nOut,
                   int *execute)
{
  int          nField, major, minor, *ranks, *fInOut;
  char         **fields;
  void         *instStore = NULL;
  capsProblem  *problem;

  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aname                == NULL)      return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_QUERYANALYSIS;

  major = CAPSMAJOR;
  minor = CAPSMINOR;

  /* try to load the AIM and get the info */
  *execute = 1;
  return aim_Initialize(&problem->aimFPTR, aname, execute, NULL, NULL, &major,
                        &minor, nIn, nOut, &nField, &fields, &ranks, &fInOut,
                        &instStore);
}


int
caps_getInput(capsObject *pobject, const char *aname, int index, char **ainame,
              capsValue *defaults)
{
  int         stat;
  capsProblem *problem;

  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aname                == NULL)      return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_GETINPUT;

  defaults->length          = defaults->nrow = defaults->ncol = 1;
  defaults->type            = Integer;
  defaults->dim             = defaults->pIndex = defaults->index = 0;
  defaults->lfixed          = defaults->sfixed = Fixed;
  defaults->nullVal         = NotAllowed;
  defaults->units           = NULL;
  defaults->meshWriter      = NULL;
  defaults->link            = NULL;
  defaults->vals.integer    = 0;
  defaults->limits.dlims[0] = defaults->limits.dlims[1] = 0.0;
  defaults->linkMethod      = Copy;
  defaults->gInType         = 0;
  defaults->partial         = NULL;
  defaults->nderiv          = 0;
  defaults->derivs          = NULL;

  stat = aim_Inputs(problem->aimFPTR, aname, NULL, NULL, index, ainame,
                    defaults);
  if (stat == CAPS_SUCCESS) defaults->length = defaults->ncol*defaults->nrow;
  return stat;
}


int
caps_getOutput(capsObject *pobject, const char *aname, int index, char **aoname,
               capsValue *form)
{
  int         stat;
  capsProblem *problem;

  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aname                == NULL)      return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_GETOUTPUT;

  form->length          = form->nrow = form->ncol = 1;
  form->type            = Integer;
  form->dim             = form->pIndex = form->index = 0;
  form->lfixed          = form->sfixed = Fixed;
  form->nullVal         = NotAllowed;
  form->units           = NULL;
  form->meshWriter      = NULL;
  form->link            = NULL;
  form->vals.integer    = 0;
  form->limits.dlims[0] = form->limits.dlims[1] = 0.0;
  form->linkMethod      = Copy;
  form->gInType         = 0;
  form->partial         = NULL;
  form->nderiv          = 0;
  form->derivs          = NULL;

  stat = aim_Outputs(problem->aimFPTR, aname, NULL, NULL, index, aoname, form);
  if (stat == CAPS_SUCCESS) form->length = form->ncol*form->nrow;
  return stat;
}


int
caps_AIMbackdoor(capsObject *aobject, const char *JSONin, char **JSONout)
{
  int          stat, ret;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  CAPSLONG     sNum;
  capsJrnl     args[1];

  *JSONout = NULL;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;
  pobject = (capsObject *)  aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  problem  = (capsProblem *)  pobject->blind;

  args[0].type   = jString;
  if (problem->dbFlag == 0) {
    stat           = caps_jrnlRead(CAPS_AIMBACKDOOR, problem, aobject, 1, args, &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *JSONout = EG_strdup(args[0].members.string);
      return ret;
    }
  }

  ret = aim_Backdoor(problem->aimFPTR, analysis->loadName, analysis->instStore,
                     &analysis->info, JSONin, JSONout);
  if (problem->dbFlag == 1) return ret;
  args[0].members.string = *JSONout;
  caps_jrnlWrite(CAPS_AIMBACKDOOR, problem, aobject, ret, 1, args,
                 problem->sNum, problem->sNum);

  return ret;
}


int
caps_system(capsObject *aobject, /*@null@*/ const char *rpath,
            const char *command)
{
  size_t       status, len;
  int          stat, ret;
  char         *fullcommand;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  CAPSLONG     sNum;
  capsJrnl     args[1];    /* not really used */

  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;
  pobject = (capsObject *)  aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  problem  = (capsProblem *)  pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  args[0].type = jString;
  stat         = caps_jrnlRead(CAPS_SYSTEM, problem, aobject, 0, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL)    return ret;

  sNum = problem->sNum;
  args[0].members.string = NULL;
  if (analysis->pre.sNum < aobject->last.sNum) {
    ret = CAPS_CLEAN;
    caps_jrnlWrite(CAPS_SYSTEM, problem, aobject, ret, 0, args, sNum,
                   problem->sNum);
    return ret;
  }
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
  ret = CAPS_SUCCESS;
  if (status != 0) {
    printf( "Failed to execute! (caps_system): %s\n",fullcommand);
    ret = CAPS_EXECERR;
  } else {
    problem->sNum += 1;
  }
  EG_free(fullcommand);
  caps_jrnlWrite(CAPS_SYSTEM, problem, aobject, ret, 0, args, sNum,
                 problem->sNum);

  return ret;
}


static int
caps_makeAnalysiX(capsObject *pobject, const char *aname,
                  /*@null@*/ const char *name, /*@null@*/ const char *unitSys,
                  /*@null@*/ const char *intents, int *exec,
                  capsObject **aobject, int *nErr, capsErrs **errors)
{
  int          i, j, status, nIn, nOut, *ranks, *fInOut, eFlag, nField = 0;
  int          major, minor, len = 0;
  char         dirName[PATH_MAX], filename[PATH_MAX], temp[PATH_MAX];
  char         **fields = NULL, *apath = NULL;
  void         *instStore;
  capsObject   *object, **tmp;
  capsProblem  *problem;
  capsAnalysis *analy, *analysis = NULL;
  capsValue    *value;
  FILE         *fp;

  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_MAKEANALYSIS;

  /* is the name unique? */
  if (name != NULL) {
    status = caps_isNameOK(name);
    if (status != CAPS_SUCCESS) return status;
    apath = EG_strdup(name);
    if (apath == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nAnalysis; i++) {
      if (problem->analysis[i]       == NULL) continue;
      if (problem->analysis[i]->name == NULL) continue;
      if (strcmp(apath, problem->analysis[i]->name) == 0) {
        status = CAPS_BADNAME;
        goto cleanup;
      }
    }
  } else {
    len   = strlen(aname) + 5;
    apath = (char *) EG_alloc(len*sizeof(char));
    if (apath == NULL) return EGADS_MALLOC;
    for (j = 0; j < 1000; j++) {
      snprintf(apath, len, "%s%d", aname, j);
      for (i = 0; i < problem->nAnalysis; i++) {
        if (problem->analysis[i]       == NULL) continue;
        if (problem->analysis[i]->name == NULL) continue;
        if (strcmp(apath, problem->analysis[i]->name) == 0) break;
      }
      if (i == problem->nAnalysis) break;
    }
    if (j == 1000) {
      status = CAPS_BADNAME;
      goto cleanup;
    }
  }

#ifdef WIN32
  snprintf(dirName, PATH_MAX, "%s\\%s", problem->root, apath);
#else
  snprintf(dirName, PATH_MAX, "%s/%s",  problem->root, apath);
#endif
  status = caps_statFile(dirName);
  if (status != EGADS_NOTFOUND) {
    snprintf(temp, PATH_MAX, "%s already exists in the phase (makeAnalysis)!",
             apath);
    caps_makeSimpleErr(pobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    status = CAPS_DIRERR;
    goto cleanup;
  }

  /* initialize the analysis structure */
  analysis = (capsAnalysis *) EG_alloc(sizeof(capsAnalysis));
  if (analysis == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  analysis->loadName         = EG_strdup(aname);
  analysis->fullPath         = EG_strdup(dirName);
  analysis->path             = apath;
  analysis->unitSys          = EG_strdup(unitSys);
  analysis->major            = CAPSMAJOR;
  analysis->minor            = CAPSMINOR;
  analysis->instStore        = NULL;
  analysis->autoexec         = 0;
  analysis->eFlag            = 0;
  analysis->reload           = 0;
  analysis->intents          = EG_strdup(intents);
  analysis->nField           = 0;
  analysis->fields           = NULL;
  analysis->ranks            = NULL;
  analysis->fInOut           = NULL;
  analysis->nAnalysisIn      = 0;
  analysis->analysisIn       = NULL;
  analysis->nAnalysisOut     = 0;
  analysis->analysisOut      = NULL;
  analysis->nAnalysisDynO    = 0;
  analysis->analysisDynO     = NULL;
  analysis->nBody            = 0;
  analysis->bodies           = NULL;
  analysis->nTess            = 0;
  analysis->tess             = NULL;
  analysis->pre.index        = -1;
  analysis->pre.pname        = NULL;
  analysis->pre.pID          = NULL;
  analysis->pre.user         = NULL;
  analysis->pre.sNum         = 0;
  analysis->info.magicnumber         = CAPSMAGIC;
  analysis->info.instance            = -1;
  analysis->info.inPost              = 0;
  analysis->info.problem             = problem;
  analysis->info.analysis            = analysis;
  analysis->info.pIndex              = 0;
  analysis->info.irow                = 0;
  analysis->info.icol                = 0;
  analysis->info.errs.nError         = 0;
  analysis->info.errs.errors         = NULL;
  analysis->info.wCntxt.aimWriterNum = 0;
  for (i = 0; i < 6; i++) analysis->pre.datetime[i] = 0;

  apath  = NULL;

  /* try to load the AIM */
  eFlag     = 0;
  nField    = 0;
  fields    = NULL;
  ranks     = NULL;
  fInOut    = NULL;
  major     = CAPSMAJOR;
  minor     = CAPSMINOR;
  instStore = NULL;
  status    = aim_Initialize(&problem->aimFPTR, aname, &eFlag, unitSys,
                             &analysis->info, &major, &minor, &nIn, &nOut,
                             &nField, &fields, &ranks, &fInOut, &instStore);
  if (status <  CAPS_SUCCESS) goto cleanup;
  if (nIn    <= 0) {
    status = CAPS_BADINIT;
    goto cleanup;
  }

  analysis->major        = major;
  analysis->minor        = minor;
  analysis->instStore    = instStore;
  analysis->eFlag        = eFlag;
  analysis->nField       = nField;
  analysis->fields       = fields;
  analysis->ranks        = ranks;
  analysis->fInOut       = fInOut;
  analysis->nAnalysisIn  = nIn;
  analysis->nAnalysisOut = nOut;
  if ((*exec == 1) && (eFlag == 1)) analysis->autoexec = 1;
  analysis->info.instance = status;

  *exec                  = eFlag;
  instStore              = NULL;
  fields                 = NULL;
  ranks                  = NULL;
  fInOut                 = NULL;

  /* allocate the objects for input */
  if (nIn != 0) {
    analysis->analysisIn = (capsObject **) EG_alloc(nIn*sizeof(capsObject *));
    if (analysis->analysisIn == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nIn; i++) analysis->analysisIn[i] = NULL;
    value = (capsValue *) EG_alloc(nIn*sizeof(capsValue));
    if (value == NULL) {
      EG_free(analysis->analysisIn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    problem->sNum += 1;
    for (i = 0; i < nIn; i++) {
      value[i].length          = value[i].nrow = value[i].ncol = 1;
      value[i].type            = Integer;
      value[i].dim             = value[i].pIndex = value[i].index = 0;
      value[i].lfixed          = value[i].sfixed = Fixed;
      value[i].nullVal         = NotAllowed;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv          = 0;
      value[i].derivs          = NULL;

      status = caps_makeObject(&object);
      if (status != CAPS_SUCCESS) {
        EG_free(value);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      if (i == 0) object->blind = value;
      object->parent    = NULL;
      object->name      = NULL;
      object->type      = VALUE;
      object->subtype   = ANALYSISIN;
      object->last.sNum = problem->sNum;
/*@-immediatetrans@*/
      object->blind     = &value[i];
/*@+immediatetrans@*/
      analysis->analysisIn[i] = object;
    }

    for (i = 0; i < nIn; i++) {
      status = aim_Inputs(problem->aimFPTR, aname, analysis->instStore,
                          &analysis->info, i+1, &analysis->analysisIn[i]->name,
                          &value[i]);
      if (status != CAPS_SUCCESS) goto cleanup;
      value[i].index = i+1;

      status = caps_checkValueObj(analysis->analysisIn[i]);
      if (status != CAPS_SUCCESS) {
        snprintf(temp, PATH_MAX, "Check Analysis '%s' Inputs Value Object %d!",
                 aname, status);
        caps_makeSimpleErr(pobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        goto cleanup;
      }
    }
  }

  /* allocate the objects for output */
  if (nOut != 0) {
    analysis->analysisOut = (capsObject **) EG_alloc(nOut*sizeof(capsObject *));
    if (analysis->analysisOut == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nOut; i++) analysis->analysisOut[i] = NULL;

    value = (capsValue *) EG_alloc(nOut*sizeof(capsValue));
    if (value == NULL) { status = EGADS_MALLOC; goto cleanup; }

    problem->sNum += 1;
    for (i = 0; i < nOut; i++) {
      value[i].length          = value[i].nrow = value[i].ncol = 1;
      value[i].type            = Integer;
      value[i].dim             = value[i].pIndex = value[i].index = 0;
      value[i].lfixed          = value[i].sfixed = Fixed;
      value[i].nullVal         = NotAllowed;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv          = 0;
      value[i].derivs          = NULL;

      status = caps_makeObject(&object);
      if (status != CAPS_SUCCESS) {
        EG_free(value);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      if (i == 0) object->blind = value;
      object->parent    = NULL;
      object->name      = NULL;
      object->type      = VALUE;
      object->subtype   = ANALYSISOUT;
      object->last.sNum = problem->sNum;
/*@-immediatetrans@*/
      object->blind     = &value[i];
/*@+immediatetrans@*/
      analysis->analysisOut[i] = object;
    }

    for (i = 0; i < nOut; i++) {
      status = aim_Outputs(problem->aimFPTR, aname, analysis->instStore,
                           &analysis->info, i+1,
                           &analysis->analysisOut[i]->name, &value[i]);
      if (status != CAPS_SUCCESS) goto cleanup;
      value[i].index = i+1;

      status = caps_checkValueObj(analysis->analysisOut[i]);
      if (status != CAPS_SUCCESS) {
        snprintf(temp, PATH_MAX, "Check Analysis '%s' Outputs Value Object %d!",
                 aname, status);
        caps_makeSimpleErr(pobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        goto cleanup;
      }
    }
  }

  /* get a place in the problem to store the data away */
  if (problem->analysis == NULL) {
    problem->analysis = (capsObject **) EG_alloc(sizeof(capsObject *));
    if (problem->analysis == NULL) { status = EGADS_MALLOC; goto cleanup; }
  } else {
    tmp = (capsObject **) EG_reall( problem->analysis,
                                   (problem->nAnalysis+1)*sizeof(capsObject *));
    if (tmp == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    problem->analysis = tmp;
  }

  /* make the directory */
  status = caps_mkDir(dirName);
  if (status != EGADS_SUCCESS) {
    snprintf(temp, PATH_MAX, "Cannot make %s (caps_makeAnalysis)", dirName);
    caps_makeSimpleErr(pobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    status = CAPS_DIRERR;
    goto cleanup;
  }

  /* get the analysis object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) goto cleanup;

  /* leave sNum 0 to flag we are unexecuted */
  object->parent = pobject;
  object->name   = EG_strdup(analysis->path);
  if (object->name == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  object->type   = ANALYSIS;
  object->blind  = analysis;
/*@-nullderef@*/
/*@-kepttrans@*/
  for (i = 0; i < nIn;  i++) analysis->analysisIn[i]->parent  = object;
  for (i = 0; i < nOut; i++) analysis->analysisOut[i]->parent = object;
/*@+kepttrans@*/
/*@+nullderef@*/
  *aobject = object;
  analysis = NULL;

  problem->analysis[problem->nAnalysis] = object;
  problem->nAnalysis += 1;

  /* setup for restarts */
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\analy.txt", problem->root);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx",  problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/analy.txt",   problem->root);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",    problem->root);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Warning: Cannot open %s (caps_makeAnalysis)\n", filename);
  } else {
    fprintf(fp, "%d\n", problem->nAnalysis);
    if (problem->analysis != NULL)
      for (i = 0; i < problem->nAnalysis; i++) {
        nIn   = nOut = 0;
        analy = (capsAnalysis *) problem->analysis[i]->blind;
        if (analy != NULL) {
          nIn  = analy->nAnalysisIn;
          nOut = analy->nAnalysisOut;
        }
        fprintf(fp, "%d %d %s\n", nIn, nOut, problem->analysis[i]->name);
      }
    fclose(fp);
    status = caps_rename(temp, filename);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: Cannot rename %s!\n", filename);
  }
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\AN-%s",
           problem->root, object->name);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/AN-%s",
           problem->root, object->name);
#endif
  status = caps_mkDir(filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Warning: Cant make dir %s (caps_makeAnalysis)\n", filename);
    status = CAPS_SUCCESS;
    goto cleanup;
  } else {
    status = caps_dumpAnalysis(problem, object);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Warning: caps_dumpAnalysis = %d (caps_makeAnalysis)\n",
             status);
      status = CAPS_SUCCESS;
      goto cleanup;
    }
  }

cleanup:
  EG_free(apath);
  if (fields != NULL) {
    for (i = 0; i < nField; i++) EG_free(fields[i]);
    EG_free(fields);
  }

  if (status != CAPS_SUCCESS) {
    if (analysis != NULL) caps_freeAnalysis(0, analysis);
  } else {
    caps_getAIMerrs((capsAnalysis *) object->blind, nErr, errors);
  }

  return status;
}


int
caps_makeAnalysis(capsObject *pobject, const char *aname,
                  /*@null@*/ const char *name, /*@null@*/ const char *unitSys,
                  /*@null@*/ const char *intents, int *exec,
                  capsObject **aobject, int *nErr, capsErrs **errors)
{
  int         stat, ret;
  CAPSLONG    sNum;
  capsProblem *problem;
  capsJrnl    args[4];

  *aobject = NULL;
  *nErr    = 0;
  *errors  = NULL;
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  if (aname                == NULL)      return CAPS_NULLNAME;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  args[0].type = jObject;
  args[1].type = jInteger;
  args[2].type = jInteger;
  args[3].type = jErr;
  stat         = caps_jrnlRead(CAPS_MAKEANALYSIS, problem, *aobject, 4, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *aobject = args[0].members.obj;
    *exec    = args[1].members.integer;
    *nErr    = args[2].members.integer;
    *errors  = args[3].members.errs;
    return ret;
  }

  sNum = problem->sNum;
  ret  = caps_makeAnalysiX(pobject, aname, name, unitSys, intents, exec,
                           aobject, nErr, errors);
  args[0].members.obj     = *aobject;
  args[1].members.integer = *exec;
  args[2].members.integer = *nErr;
  args[3].members.errs    = *errors;
  caps_jrnlWrite(CAPS_MAKEANALYSIS, problem, *aobject, ret, 4, args, sNum,
                 problem->sNum);

  return ret;
}


static int
caps_dupAnalysiX(capsObject *from, const char *name, capsObject **aobject)
{
  int          i, j, status, nIn, nOut, *ranks, *fInOut, eFlag, nField = 0;
  int          major, minor, index, len = 0;
  char         dirName[PATH_MAX], filename[PATH_MAX], temp[PATH_MAX];
  char         **fields = NULL, *apath = NULL;
  void         *instStore;
  capsObject   *pobject, *object, **tmp;
  capsAnalysis *froma, *analy, *analysis = NULL;
  capsProblem  *problem;
  capsValue    *value = NULL, *src;
  FILE         *fp;

  *aobject = NULL;
  froma    = (capsAnalysis *) from->blind;
  pobject  = (capsObject *)   from->parent;
  problem  = (capsProblem *)  pobject->blind;

  /* is the name unique? */
  if (name != NULL) {
    status = caps_isNameOK(name);
    if (status != CAPS_SUCCESS)       return status;
    apath = EG_strdup(name);
    if (apath == NULL) return EGADS_MALLOC;
  } else {
    index = aim_Index(problem->aimFPTR, froma->loadName);
    j     = problem->aimFPTR.aim_nInst[index];
    len   = strlen(froma->loadName) + 8;
    apath = (char *) EG_alloc(len*sizeof(char));
    if (apath == NULL) return EGADS_MALLOC;
    snprintf(apath, len, "%s%d", froma->loadName, j);
  }
  for (i = 0; i < problem->nAnalysis; i++) {
    if (problem->analysis[i]       == NULL) continue;
    if (problem->analysis[i]->name == NULL) continue;
    if (strcmp(apath, problem->analysis[i]->name) == 0) {
      status = CAPS_BADNAME;
      goto cleanup;
    }
  }
#ifdef WIN32
  snprintf(dirName, PATH_MAX, "%s\\%s", problem->root, apath);
#else
  snprintf(dirName, PATH_MAX, "%s/%s",  problem->root, apath);
#endif
  status = caps_statFile(dirName);
  if (status != EGADS_NOTFOUND) {
    printf(" CAPS Error: %s already exists in the phase!\n", apath);
    status = CAPS_DIRERR;
    goto cleanup;
  }

  /* initialize the analysis structure */
  analysis = (capsAnalysis *) EG_alloc(sizeof(capsAnalysis));
  if (analysis == NULL) { status = EGADS_MALLOC; goto cleanup; }

  analysis->loadName         = EG_strdup(froma->loadName);
  analysis->fullPath         = EG_strdup(dirName);
  analysis->path             = apath;
  analysis->unitSys          = EG_strdup(froma->unitSys);
  analysis->major            = CAPSMAJOR;
  analysis->minor            = CAPSMINOR;
  analysis->instStore        = NULL;
  analysis->autoexec         = froma->autoexec;
  analysis->eFlag            = froma->eFlag;
  analysis->reload           = 0;
  analysis->intents          = EG_strdup(froma->intents);
  analysis->nField           = 0;
  analysis->fields           = NULL;
  analysis->ranks            = NULL;
  analysis->fInOut           = NULL;
  analysis->nAnalysisIn      = 0;
  analysis->analysisIn       = NULL;
  analysis->nAnalysisOut     = 0;
  analysis->analysisOut      = NULL;
  analysis->nAnalysisDynO    = 0;
  analysis->analysisDynO     = NULL;
  analysis->nBody            = 0;
  analysis->bodies           = NULL;
  analysis->nTess            = 0;
  analysis->tess             = NULL;
  analysis->pre.index        = -1;
  analysis->pre.pname        = NULL;
  analysis->pre.pID          = NULL;
  analysis->pre.user         = NULL;
  analysis->pre.sNum         = 0;
  analysis->info.magicnumber         = CAPSMAGIC;
  analysis->info.inPost              = 0;
  analysis->info.problem             = problem;
  analysis->info.analysis            = analysis;
  analysis->info.pIndex              = 0;
  analysis->info.irow                = 0;
  analysis->info.icol                = 0;
  analysis->info.errs.nError         = 0;
  analysis->info.errs.errors         = NULL;
  analysis->info.wCntxt.aimWriterNum = 0;
  for (i = 0; i < 6; i++) analysis->pre.datetime[i] = 0;

  apath  = NULL;

  /* get a new instance AIM */
  eFlag     = 0;
  nField    = 0;
  fields    = NULL;
  ranks     = NULL;
  fInOut    = NULL;
  major     = CAPSMAJOR;
  minor     = CAPSMINOR;
  instStore = NULL;
  status    = aim_Initialize(&problem->aimFPTR, froma->loadName, &eFlag,
                             froma->unitSys, &analysis->info,
                             &major, &minor, &nIn, &nOut,
                             &nField, &fields, &ranks, &fInOut, &instStore);
  if (status <  CAPS_SUCCESS) goto cleanup;
  if (nIn    <= 0) {
    status = CAPS_BADINIT;
    goto cleanup;
  }

  analysis->major         = major;
  analysis->minor         = minor;
  analysis->instStore     = instStore;
  analysis->eFlag         = eFlag;
  analysis->nField        = nField;
  analysis->fields        = fields;
  analysis->ranks         = ranks;
  analysis->fInOut        = fInOut;
  analysis->nAnalysisIn   = nIn;
  analysis->nAnalysisOut  = nOut;
  analysis->info.instance = status;
  if ((analysis->autoexec == 1) && (eFlag == 0)) analysis->autoexec = 0;

  instStore = NULL;
  fields    = NULL;
  ranks     = NULL;
  fInOut    = NULL;

  /* allocate the objects for input */
  if (nIn != 0) {
    analysis->analysisIn = (capsObject **) EG_alloc(nIn*sizeof(capsObject *));
    if (analysis->analysisIn == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }

    for (i = 0; i < nIn; i++) analysis->analysisIn[i] = NULL;
    value = (capsValue *) EG_alloc(nIn*sizeof(capsValue));
    if (value == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }

    problem->sNum += 1;
    for (i = 0; i < nIn; i++) {
      value[i].length          = value[i].nrow = value[i].ncol = 1;
      value[i].type            = Integer;
      value[i].dim             = value[i].pIndex = value[i].index = 0;
      value[i].lfixed          = value[i].sfixed = Fixed;
      value[i].nullVal         = NotAllowed;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv          = 0;
      value[i].derivs          = NULL;

      status = caps_makeObject(&object);
      if (status != CAPS_SUCCESS) {
        EG_free(value);
        goto cleanup;
      }
      if (i == 0) object->blind = value;
      object->parent    = NULL;
      object->name      = EG_strdup(froma->analysisIn[i]->name);
      if (object->name == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
      object->type      = VALUE;
      object->subtype   = ANALYSISIN;
      object->last.sNum = problem->sNum;
/*@-immediatetrans@*/
      object->blind     = &value[i];
/*@+immediatetrans@*/
      analysis->analysisIn[i] = object;
    }

    src = (capsValue *) froma->analysisIn[0]->blind;
    if ((nIn != 0) && (src == NULL)) {
      status = CAPS_NULLBLIND;
      goto cleanup;
    }
    for (i = 0; i < nIn; i++) {
      status = caps_dupValues(&src[i], &value[i]);
      if (status != CAPS_SUCCESS) goto cleanup;
    }
  }

  /* allocate the objects for output */
  if (nOut != 0) {
    analysis->analysisOut = (capsObject **) EG_alloc(nOut*sizeof(capsObject *));
    if (analysis->analysisOut == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nOut; i++) analysis->analysisOut[i] = NULL;

    value = (capsValue *) EG_alloc(nOut*sizeof(capsValue));
    if (value == NULL) { status = EGADS_MALLOC; goto cleanup; }

    problem->sNum += 1;
    for (i = 0; i < nOut; i++) {
      value[i].length          = value[i].nrow = value[i].ncol = 1;
      value[i].type            = Integer;
      value[i].dim             = value[i].pIndex = value[i].index = 0;
      value[i].lfixed          = value[i].sfixed = Fixed;
      value[i].nullVal         = NotAllowed;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv          = 0;
      value[i].derivs          = NULL;

      status = caps_makeObject(&object);
      if (status != CAPS_SUCCESS) {
        EG_free(value);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      if (i == 0) object->blind = value;
      object->parent    = NULL;
      object->name      = EG_strdup(froma->analysisOut[i]->name);
      if (object->name == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
      object->type      = VALUE;
      object->subtype   = ANALYSISOUT;
      object->last.sNum = problem->sNum;
/*@-immediatetrans@*/
      object->blind     = &value[i];
/*@+immediatetrans@*/
      analysis->analysisOut[i] = object;
    }

    src = (capsValue *) froma->analysisOut[0]->blind;
    if ((nOut != 0) && (src == NULL)) {
      status = CAPS_NULLBLIND;
      goto cleanup;
    }
    for (i = 0; i < nOut; i++) {
      status = caps_dupValues(&src[i], &value[i]);
      if (status != CAPS_SUCCESS) goto cleanup;
    }
  }

  /* get a place in the problem to store the data away */
  if (problem->analysis == NULL) {
    problem->analysis = (capsObject **) EG_alloc(sizeof(capsObject *));
    if (problem->analysis == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
  } else {
    tmp = (capsObject **) EG_reall( problem->analysis,
                                   (problem->nAnalysis+1)*sizeof(capsObject *));
    if (tmp == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    problem->analysis = tmp;
  }

  /* make the directory */
  status = caps_mkDir(dirName);
  if (status != EGADS_SUCCESS) {
    printf(" CAPS Error: Cannot make %s (caps_dupAnalysis)\n", dirName);
    status = CAPS_DIRERR;
    goto cleanup;
  }

  /* get the analysis object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) goto cleanup;

  /* leave sNum 0 to flag we are unexecuted */
  object->parent = pobject;
  object->name   = EG_strdup(analysis->path);
  if (object->name == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  object->type   = ANALYSIS;
  object->blind  = analysis;
/*@-nullderef@*/
/*@-kepttrans@*/
  for (i = 0; i < nIn;  i++) analysis->analysisIn[i]->parent  = object;
  for (i = 0; i < nOut; i++) analysis->analysisOut[i]->parent = object;
/*@+kepttrans@*/
/*@+nullderef@*/

  *aobject = object;
  analysis = NULL;

  problem->analysis[problem->nAnalysis] = object;
  problem->nAnalysis += 1;

  /* setup for restarts */
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\analy.txt", problem->root);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx",  problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/analy.txt",   problem->root);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",    problem->root);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Warning: Cannot open %s (caps_dupAnalysis)\n", filename);
  } else {
    fprintf(fp, "%d\n", problem->nAnalysis);
    if (problem->analysis != NULL)
      for (i = 0; i < problem->nAnalysis; i++) {
        nIn   = nOut = 0;
        analy = (capsAnalysis *) problem->analysis[i]->blind;
        if (analy != NULL) {
          nIn  = analy->nAnalysisIn;
          nOut = analy->nAnalysisOut;
        }
        fprintf(fp, "%d %d %s\n", nIn, nOut, problem->analysis[i]->name);
      }
    fclose(fp);
    status = caps_rename(temp, filename);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: Cannot rename %s!\n", filename);
  }
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\AN-%s",
           problem->root, object->name);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/AN-%s",
           problem->root, object->name);
#endif
  status = caps_mkDir(filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Warning: Cant make dir %s (caps_dupAnalysis)\n", filename);
    status = CAPS_SUCCESS;
    goto cleanup;
  } else {
    status = caps_dumpAnalysis(problem, object);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Warning: caps_dumpAnalysis = %d (caps_dupAnalysis)\n",
             status);
      status = CAPS_SUCCESS;
      goto cleanup;
    }
  }

cleanup:
  EG_free(apath);
  if (fields != NULL) {
    for (i = 0; i < nField; i++) EG_free(fields[i]);
    EG_free(fields);
  }

  if (status != CAPS_SUCCESS)
    if (analysis != NULL) caps_freeAnalysis(0, analysis);

  return status;
}


int
caps_dupAnalysis(capsObject *from, const char *name, capsObject **aobject)
{
  int         stat, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[1];

  if (aobject           == NULL)      return CAPS_NULLOBJ;
  *aobject = NULL;
  if (from              == NULL)      return CAPS_NULLOBJ;
  if (from->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (from->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (from->blind       == NULL)      return CAPS_NULLBLIND;
  if (name              == NULL)      return CAPS_NULLNAME;
  if (from->parent      == NULL)      return CAPS_NULLOBJ;
  pobject = (capsObject *)  from->parent;
  if (pobject->blind    == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)      return CAPS_READONLYERR;

  args[0].type = jObject;
  stat         = caps_jrnlRead(CAPS_DUPANALYSIS, problem, *aobject, 1, args,
                              &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *aobject = args[0].members.obj;
    return ret;
  }

  sNum = problem->sNum;
  ret  = caps_dupAnalysiX(from, name, aobject);
  args[0].members.obj = *aobject;
  caps_jrnlWrite(CAPS_DUPANALYSIS, problem, *aobject, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_resetAnalysis(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int          status;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;

  if (nErr                 == NULL)       return CAPS_NULLVALUE;
  if (errors               == NULL)       return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if (aobject              == NULL)       return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)   return CAPS_BADTYPE;
  if (aobject->blind       == NULL)       return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  status   = caps_findProblem(aobject, CAPS_RESETANALYSIS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)          return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  status = caps_rmDir(analysis->fullPath);
  if (status != EGADS_SUCCESS) {
    caps_makeSimpleErr(pobject, CERROR,
                       "Cannot remove directory (caps_resetAnalysis):",
                       analysis->fullPath, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return status;
  }
  status = caps_mkDir(analysis->fullPath);
  if (status != EGADS_SUCCESS) {
    caps_makeSimpleErr(pobject, CERROR,
                       "Cannot make directory (caps_resetAnalysis):",
                       analysis->fullPath, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return status;
  }

  return CAPS_SUCCESS;
}


int
caps_analysisInfX(capsObject *aobject, char **apath, char **unitSys, int *major,
                  int *minor, char **intents, int *nField, char ***fnames,
                  int **ranks, int **fInOut, int *execute, int *status)
{
  int          i, stat, ret, gstatus;
  int          majorX, minorX, nFieldX, execX, dirtyX;
  int          *ranksX, *fInOutX;
  char         *intentsX, *apathX, *unitSysX, **fnamesX;
  CAPSLONG     sn;
  capsAnalysis *analysis;
  capsObject   *pobject, *source, *object, *last;
  capsValue    *value;
  capsProblem  *problem;

  *nField  = *status = *major = *minor = 0;
  *apath   = NULL;
  *unitSys = NULL;
  *fnames  = NULL;
  *ranks   = NULL;
  *fInOut  = NULL;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;
  pobject  = (capsObject *)   aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  *apath   = analysis->fullPath;
  *unitSys = analysis->unitSys;
  *execute = analysis->eFlag;
  *major   = analysis->major;
  *minor   = analysis->minor;
  *intents = analysis->intents;
  *nField  = analysis->nField;
  *fnames  = analysis->fields;
  *ranks   = analysis->ranks;
  *fInOut  = analysis->fInOut;
  if ((analysis->autoexec == 1) && (analysis->eFlag == 1)) *execute = 2;

  /* are we "geometry" clean? */
  gstatus = 0;
  if (pobject->subtype == PARAMETRIC) {
    /* check for dirty geometry inputs */
    for (i = 0; i < problem->nGeomIn; i++) {
      source = object = problem->geomIn[i];
      do {
        if (source->magicnumber != CAPSMAGIC) {
          ret = CAPS_BADOBJECT;
          goto done;
        }
        if (source->type        != VALUE) {
          ret = CAPS_BADTYPE;
          goto done;
        }
        if (source->blind       == NULL) {
          ret = CAPS_NULLBLIND;
          goto done;
        }
        value = (capsValue *) source->blind;
        if (value->link == object) {
          ret = CAPS_CIRCULARLINK;
          goto done;
        }
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      if (last->last.sNum > problem->geometry.sNum) {
        gstatus = 1;
        break;
      }
    }
  }

  /* are we "analysis" clean? */
  if (analysis->pre.sNum == 0) {
    *status = 1;
  } else {
    /* check for dirty inputs */
    for (i = 0; i < analysis->nAnalysisIn; i++) {
      source = object = analysis->analysisIn[i];
      do {
        if (source->magicnumber != CAPSMAGIC) {
          ret = CAPS_BADOBJECT;
          goto done;
        }
        if (source->type        != VALUE) {
          ret = CAPS_BADTYPE;
          goto done;
        }
        if (source->blind       == NULL) {
          ret = CAPS_NULLBLIND;
          goto done;
        }
        value = (capsValue *) source->blind;
        if (value->link == object) {
          ret = CAPS_CIRCULARLINK;
          goto done;
        }
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      if (last->last.sNum > analysis->pre.sNum) {
        *status = 1;
        break;
      }
      /* if the value is linked, check the state of the linked parent */
      if ((last != object) && (last->parent != NULL) &&
          (last->parent->type == ANALYSIS)) {
        stat = caps_analysisInfX(last->parent, &apathX, &unitSysX, &majorX,
                                 &minorX, &intentsX, &nFieldX, &fnamesX,
                                 &ranksX, &fInOutX, &execX, &dirtyX);
        if (stat != CAPS_SUCCESS) return stat;
        if (dirtyX > 0) {
          *status = 1;
          break;
        }
      }
    }
    if (*status == 0) {
      stat = caps_snDataSets(aobject, 0, &sn);
      if (stat == CAPS_SUCCESS)
        if (sn > analysis->pre.sNum) *status = 1;
    }
  }
  *status += gstatus*2;

  /* is geometry new? */
  if (*status == 0)
    if (problem->geometry.sNum > analysis->pre.sNum) *status = 4;

  /* is post required? */
  if (*status == 0)
    if (analysis->pre.sNum > aobject->last.sNum) {
      *status = 5;
      if (analysis->eFlag == 0) *status = 6;
    }

  ret = CAPS_SUCCESS;

done:
  return ret;
}


int
caps_analysisInfo(capsObject *aobject, char **apath, char **unitSys, int *major,
                  int *minor, char **intents, int *nField, char ***fnames,
                  int **ranks, int **fInOut, int *execute, int *status)
{
  int          stat, ret;
  CAPSLONG     sNum;
  capsAnalysis *analysis;
  capsObject   *pobject;
  capsProblem  *problem;
  capsJrnl     args[1];

  *nField  = *status = *major = *minor = 0;
  *apath   = NULL;
  *unitSys = NULL;
  *fnames  = NULL;
  *ranks   = NULL;
  *fInOut  = NULL;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;
  pobject  = (capsObject *)   aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  *apath   = analysis->fullPath;
  *unitSys = analysis->unitSys;
  *execute = analysis->eFlag;
  *major   = analysis->major;
  *minor   = analysis->minor;
  *intents = analysis->intents;
  *nField  = analysis->nField;
  *fnames  = analysis->fields;
  *ranks   = analysis->ranks;
  *fInOut  = analysis->fInOut;
  if ((analysis->autoexec == 1) && (analysis->eFlag == 1)) *execute = 2;

  args[0].type = jInteger;
  if (problem->dbFlag == 0) {
    stat = caps_jrnlRead(CAPS_ANALYSISINFO, problem, aobject, 1, args, &sNum,
                         &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *status = args[0].members.integer;
      return ret;
    }
  }

  sNum = problem->sNum;
  ret  = caps_analysisInfX(aobject, apath, unitSys, major, minor, intents,
                           nField, fnames, ranks, fInOut, execute, status);
  if (problem->dbFlag == 1) return ret;
  args[0].members.integer = *status;
  caps_jrnlWrite(CAPS_ANALYSISINFO, problem, aobject, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_circularAutoExecs(capsObject *asrc, /*@null@*/ capsObject *aobject)
{
  int           i, j, k, stat;
  int           found;
  capsAnalysis  *analysis;
  capsObject    *object, *source, *last, *pobject, *aObj=NULL;
  capsValue     *value;
  capsBound     *bound;
  capsVertexSet *vertexset, *vset_link;
  capsDataSet   *dataset;
  capsProblem   *problem;

  if (asrc->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (asrc->blind       == NULL)      return CAPS_NULLBLIND;

  /* startup */
  if (asrc->type == VALUE) asrc = asrc->parent;
  else if (asrc->type == DATASET) {

    if (asrc->blind == NULL)       return CAPS_NULLBLIND;
    dataset = (capsDataSet *) asrc->blind;
    if (dataset->ftype != FieldIn) return CAPS_BADTYPE;

    source = dataset->link->parent;
    if (source              == NULL)      return CAPS_NULLOBJ;
    if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (source->type        != VERTEXSET) return CAPS_BADTYPE;
    if (source->blind       == NULL)      return CAPS_NULLBLIND;
    vertexset = source->blind;
    if (vertexset->analysis == NULL)      return CAPS_NULLOBJ;
    asrc = vertexset->analysis;
  }
  if (aobject == NULL) aObj = asrc;
  else                 aObj = aobject;

  /* extract analysis from aobject */
  if (aObj->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aObj->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aObj->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aObj->blind;
  pobject  = aObj->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  /* check if we cannot auto-execute */
  if ((analysis->autoexec != 1) || (analysis->eFlag != 1))
    return CAPS_SUCCESS;

  /* asrc came back on it self, circular link! */
  if (asrc == aobject) return CAPS_CIRCULARLINK;

  /* we can auto execute, do we have any links?
   * first check for AnalysisIn links, then FieldIn links
   */

  /* check AnalysisIn links */
  for (i = 0; i < analysis->nAnalysisIn; i++) {
    source = object = analysis->analysisIn[i];
    do {
      if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (source->type        != VALUE)     return CAPS_BADTYPE;
      if (source->blind       == NULL)      return CAPS_NULLBLIND;
      value = (capsValue *) source->blind;
      if (value->link         == object)    return CAPS_CIRCULARLINK;
      last   = source;
      source = value->link;
    } while (value->link != NULL);
    if (last != object) {
      source = last->parent;
      if (source->type == ANALYSIS) {
        stat = caps_circularAutoExecs(asrc, source);
        if (stat != CAPS_SUCCESS) return stat;
      }
    }
  }

  /* check FieldIn links */
  for (i = 0; i < problem->nBound; i++) {

    bound = (capsBound*) problem->bounds[i]->blind;
    for (found = j = 0; j < bound->nVertexSet; j++) {
      vertexset = (capsVertexSet*) bound->vertexSet[j]->blind;
      if (vertexset->analysis == aObj) {
        found = 1;
        break;
      }
    }
    if (found == 0) continue;

    for (j = 0; j < bound->nVertexSet; j++) {
      vertexset = (capsVertexSet*) bound->vertexSet[j]->blind;
      if (vertexset->analysis != aObj) continue;

      for (k = 0; k < vertexset->nDataSets; k++) {

        source = vertexset->dataSets[k];
        if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
        if (source->type        != DATASET)   return CAPS_BADTYPE;
        if (source->blind       == NULL)      return CAPS_NULLBLIND;
        dataset = (capsDataSet *) source->blind;
        if (dataset->link       == source)    return CAPS_CIRCULARLINK;

        if (dataset->ftype != FieldIn) continue;
        if (dataset->link == NULL) continue;

        source = dataset->link->parent;
        if (source              == NULL)      return CAPS_NULLOBJ;
        if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
        if (source->type        != VERTEXSET) return CAPS_BADTYPE;
        if (source->blind       == NULL)      return CAPS_NULLBLIND;
        vset_link = source->blind;
        if (vset_link->analysis == NULL)      return CAPS_NULLOBJ;
        stat = caps_circularAutoExecs(asrc, vset_link->analysis);
        if (stat != CAPS_SUCCESS) return stat;
      }
    }
  }

  return CAPS_SUCCESS;
}


/* maps gIndices indexes to body global tessellation indexes */
static int
caps_makeTessGlobal(capsDiscr *discr)
{
  int           i, k, iface, vert, local, global, bIndex, typ, len, stat;
  int           state, nGlobal;
  ego           body;
  capsBodyDiscr *discBody = NULL;

  if (discr->nPoints == 0) return CAPS_SUCCESS;
  if (discr->nBodys  == 0) return CAPS_SUCCESS;

  discr->bodys[0].globalOffset = 0;

  stat = EGADS_MALLOC;
  discr->tessGlobal = (int *) EG_alloc(2*discr->nPoints*sizeof(int));
  if (discr->tessGlobal == NULL) goto cleanup;

  for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
    discBody = &discr->bodys[bIndex-1];

    /* fill in tessGlobal */
    for (i = 0; i < discBody->nElems; i++) {
      typ = discBody->elems[i].tIndex;
      len = discr->types[typ-1].nref;
      for (k = 0; k < len; k++) {
        iface  = discBody->elems[i].eIndex;
        vert   = discBody->elems[i].gIndices[2*k  ]-1;
        local  = discBody->elems[i].gIndices[2*k+1];
        stat = EG_localToGlobal(discBody->tess, iface, local, &global);
        if (stat != EGADS_SUCCESS) goto cleanup;
        discr->tessGlobal[2*vert  ] = bIndex;
        discr->tessGlobal[2*vert+1] = global;
      }
    }

    /* fill in globalOffset across bodies */
    if (bIndex < discr->nBodys) {
      stat = EG_statusTessBody(discBody->tess, &body, &state, &nGlobal);
      if (stat != EGADS_SUCCESS) goto cleanup;
      discr->bodys[bIndex].globalOffset = discBody->globalOffset + nGlobal;
    }
  }

  stat = CAPS_SUCCESS;

cleanup:
  return stat;
}


int
caps_checkDiscr(capsDiscr *discr, int l, char *line)
{
  int           i, j, len, typ, stat, bIndex, nEdge, nFace, last, state;
  int           nGlobal, npts, ntris;
  ego           body;
  capsBodyDiscr *discBody;
  const int     *ptype, *pindex, *tris, *tric;
  const double  *xyz, *prms;

  if (discr->types == NULL) {
    snprintf(line, l, "caps_checkDiscr: types is NULL!\n");
    return CAPS_NULLBLIND;
  }

  /* do the element types */
  for (i = 0; i < discr->nTypes; i++) {
    if (discr->types[i].gst == NULL) {
      snprintf(line, l, "caps_checkDiscr: types[%d].gst = NULL!", i+1);
      return CAPS_NULLVALUE;
    }
    if ((discr->types[i].dst == NULL) && (discr->types[i].ndata != 0)) {
      snprintf(line, l, "caps_checkDiscr: types[%d].dst = NULL!", i+1);
      return CAPS_NULLVALUE;
    }
    if ((discr->types[i].matst == NULL) && (discr->types[i].nmat != 0)) {
      snprintf(line, l, "caps_checkDiscr: types[%d].matst = NULL!", i+1);
      return CAPS_NULLVALUE;
    }
    if (discr->dim == 2)
      if (discr->types[i].tris == NULL) {
        snprintf(line, l, "caps_checkDiscr: types[%d].tris = NULL!", i+1);
        return CAPS_NULLVALUE;
      }
    for (j = 0; j < 3*discr->types[i].ntri; j++)
      if ((discr->types[i].tris[j] < 1) ||
          (discr->types[i].tris[j] > discr->types[i].nref)) {
        snprintf(line, l,
                 "caps_checkDiscr: types[%d].tris[%d] = %d out of range [1-%d] ",
                 i+1, j+1, discr->types[i].tris[j], discr->types[i].nref);
        return CAPS_BADINDEX;
      }
  }

  /* check vert element indices */
  if (discr->nVerts != 0) {
    if (discr->verts == NULL) {
      snprintf(line, l, "caps_checkDiscr: nVert = %d but verts = NULL!",
               discr->nVerts);
      return CAPS_NULLVALUE;
    }
    if (discr->celem == NULL) {
      snprintf(line, l, "caps_checkDiscr: nVert = %d but celem = NULL!",
               discr->nVerts);
      return CAPS_NULLVALUE;
    }
    for (i = 0; i < discr->nVerts; i++) {
      bIndex = discr->celem[2*i];
      if ((bIndex < 1) || (bIndex > discr->nBodys)) {
        snprintf(line, l,
                 "caps_checkDiscr:celem[2*%d] = %d out of range [1-%d] ",
                 i, discr->celem[i], discr->nBodys);
        return CAPS_BADINDEX;
      }
      discBody = &discr->bodys[discr->celem[2*i]-1];

      if ((discr->celem[2*i] < 1) || (discr->celem[2*i+1] > discBody->nElems)) {
        snprintf(line, l,
                 "caps_checkDiscr: celem[2*%d+1] = %d out of range [1-%d] ",
                 i, discr->celem[i], discBody->nElems);
        return CAPS_BADINDEX;
      }
    }
  }

  /* look at body discretizations */
  if (discr->bodys == NULL) {
    snprintf(line, l, "caps_checkDiscr: bodys is NULL!\n");
    return CAPS_NULLBLIND;
  }

  for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
    discBody = &discr->bodys[bIndex-1];

    if (discBody->tess == NULL) {
      snprintf(line, l, "caps_checkDiscr: body %d etess is NULL!\n", bIndex);
      return CAPS_NULLBLIND;
    }

    if (discBody->elems == NULL) {
      snprintf(line, l, "caps_checkDiscr: body %d elems is NULL!\n", bIndex);
      return CAPS_NULLBLIND;
    }

    stat = EG_statusTessBody(discBody->tess, &body, &state, &nGlobal);
    if (stat < EGADS_SUCCESS) {
      snprintf(line, l, "caps_checkDiscr: statusTessBody = %d for body %d!\n",
               stat, bIndex);
      return stat;
    }
    if (state == 0) {
      snprintf(line, l, "caps_checkDiscr: Tessellation is Open for body %d!\n",
               bIndex);
      return EGADS_TESSTATE;
    }

    /* look at the data associated with each body */
    stat = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    if (stat != EGADS_SUCCESS) {
      snprintf(line, l, "caps_checkDiscr: getBodyTopos (Face) = %d for %d!\n",
               stat, bIndex);
      return stat;
    }
    stat = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
    if (stat != EGADS_SUCCESS) {
      snprintf(line, l, "caps_checkDiscr: getBodyTopos (Edge) = %d for %d!\n",
               stat, bIndex);
      return stat;
    }

    /* do the elements */
    last = 0;
    for (i = 1; i <= discBody->nElems; i++) {
      if ((discBody->elems[i-1].tIndex < 1) ||
          (discBody->elems[i-1].tIndex > discr->nTypes)) {
        snprintf(line, l,
                 "caps_checkDiscr: elems[%d].tIndex = %d out of range [1-%d] ",
                 i, discBody->elems[i].tIndex, discr->nTypes);
        return CAPS_BADINDEX;
      }
      if (discr->dim == 1) {
        if ((discBody->elems[i-1].eIndex < 1) ||
            (discBody->elems[i-1].eIndex > nEdge)) {
          snprintf(line, l,
                   "caps_checkDiscr: elems[%d].eIndex = %d out of range [1-%d] ",
                   i, discBody->elems[i-1].eIndex, nEdge);
          return CAPS_BADINDEX;
        }
        if (discBody->elems[i-1].eIndex != last) {
/*@-nullpass@*/
          stat = EG_getTessEdge(discBody->tess, discBody->elems[i-1].eIndex,
                                &npts, &xyz, &prms);
/*@+nullpass@*/
          if (stat != EGADS_SUCCESS) {
            snprintf(line, l, "caps_checkDiscr: getTessEdge %d = %d for %d!\n",
                     discBody->elems[i-1].eIndex, stat, bIndex);
            return stat;
          }
          ntris = npts-1;
          last  = discBody->elems[i-1].eIndex;
        }
      } else if (discr->dim == 2) {
        if ((discBody->elems[i-1].eIndex < 1) ||
            (discBody->elems[i-1].eIndex > nFace)) {
          snprintf(line, l,
                   "caps_checkDiscr: elems[%d].eIndex = %d out of range [1-%d] ",
                   i, discBody->elems[i-1].eIndex, nFace);
          return CAPS_BADINDEX;
        }
        if (discBody->elems[i-1].eIndex != last) {
/*@-nullpass@*/
          stat = EG_getTessFace(discBody->tess,
                                discBody->elems[i-1].eIndex, &npts, &xyz, &prms,
                                &ptype, &pindex, &ntris, &tris, &tric);
/*@+nullpass@*/
          if (stat != EGADS_SUCCESS) {
            snprintf(line, l, "caps_checkDiscr: getTessFace %d = %d for %d!\n",
                     discBody->elems[i-1].eIndex, stat, bIndex);
            return stat;
          }
          last = discBody->elems[i-1].eIndex;
        }
      }
      typ = discBody->elems[i-1].tIndex;
      len = discr->types[typ-1].nref;
      for (j = 0; j < len; j++) {
        if ((discBody->elems[i-1].gIndices[2*j] < 1) ||
            (discBody->elems[i-1].gIndices[2*j] > discr->nPoints)) {
          snprintf(line, l,
                   "caps_checkDiscr: elems[%d].gIndices[%d]p = %d out of range [1-%d] ",
                   i, j, discBody->elems[i-1].gIndices[2*j], discr->nPoints);
          return CAPS_BADINDEX;
        }
#ifndef __clang_analyzer__
        if ((discBody->elems[i-1].gIndices[2*j+1] < 1) ||
            (discBody->elems[i-1].gIndices[2*j+1] > npts)) {
          snprintf(line, l,
                   "caps_checkDiscr: elems[%d].gIndices[%d]t = %d out of range [1-%d] ",
                   i, j+1, discBody->elems[i-1].gIndices[2*j+1], npts);
          return CAPS_BADINDEX;
        }
#endif
      }
      if (discr->verts != NULL) {
        len = discr->types[typ-1].ndata;
        if ((len != 0) && (discBody->elems[i-1].dIndices == NULL)) {
          snprintf(line, l,
                   "caps_checkDiscr: elems[%d].dIndices[%d] == NULL ",
                   i, j+1);
          return CAPS_NULLVALUE;
        }
        for (j = 0; j < len; j++)
          if ((discBody->elems[i-1].dIndices[j] < 1) ||
              (discBody->elems[i-1].dIndices[j] > discr->nVerts)) {
            snprintf(line, l,
                     "caps_checkDiscr: elems[%d].dIndices[%d] = %d out of range [1-%d] ",
                     i, j, discBody->elems[i-1].dIndices[j], discr->nPoints);
            return CAPS_BADINDEX;
          }
      }
      len = discr->types[typ-1].ntri;
      if (len > 2) {
        if (discBody->elems[i-1].eTris.poly == NULL) {
          snprintf(line, l, "caps_checkDiscr: elems[%d].poly = NULL!", i);
          return CAPS_NULLVALUE;
        }
        for (j = 0; j < len; j++) {
#ifndef __clang_analyzer__
          if ((discBody->elems[i-1].eTris.poly[j] < 1) ||
              (discBody->elems[i-1].eTris.poly[j] > ntris)) {
            snprintf(line, l,
                     "caps_checkDiscr: elems[%d].eTris[%d] = %d out of range [1-%d] ",
                     i, j+1, discBody->elems[i-1].eTris.poly[j], ntris);
            return CAPS_BADINDEX;
          }
#endif
        }

      } else {
        for (j = 0; j < len; j++) {
#ifndef __clang_analyzer__
          if ((discBody->elems[i-1].eTris.tq[j] < 1) ||
              (discBody->elems[i-1].eTris.tq[j] > ntris)) {
            snprintf(line, l,
                     "caps_checkDiscr: elems[%d].eTris[%d] = %d out of range [1-%d] ",
                     i, j+1, discBody->elems[i].eTris.tq[j], ntris);
            return CAPS_BADINDEX;
          }
#endif
        }
      }
    }
  }

  /* create the global tessellation index and check */
  if (discr->tessGlobal == NULL) {
    stat = caps_makeTessGlobal(discr);
    if (stat < EGADS_SUCCESS || discr->tessGlobal == NULL) {
      snprintf(line, l, "caps_checkDiscr: caps_makeTessGlobal stat = %d!\n",
               stat);
      return stat;
    }
  }

  /* look at body indices */
/*@-nullderef@*/
  for (i = 0; i < discr->nPoints; i++)
    if ((discr->tessGlobal[2*i] < 1) || (discr->tessGlobal[2*i] > discr->nBodys)) {
      snprintf(line, l, "caps_checkDiscr: body mapping %d = %d [1,%d]!\n",
               i+1, discr->tessGlobal[2*i], discr->nBodys);
      return CAPS_BADINDEX;
    }
/*@+nullderef@*/

  for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
    discBody = &discr->bodys[bIndex-1];

    stat = EG_statusTessBody(discBody->tess, &body, &state, &nGlobal);
    if (stat < EGADS_SUCCESS) {
      snprintf(line, l, "caps_checkDiscr: statusTessBody = %d for body %d!\n",
               stat, bIndex);
      return stat;
    }
    if (state == 0) {
      snprintf(line, l, "caps_checkDiscr: Tessellation is Open for body %d!\n",
               bIndex);
      return EGADS_TESSTATE;
    }

    /* check point tessGlobal is withing the expected range for this body */
/*@-nullderef@*/
    for (i = 0; i < discr->nPoints; i++) {
      if (discr->tessGlobal[2*i] != bIndex) continue;
      if ((discr->tessGlobal[2*i+1] < 1) || (discr->tessGlobal[2*i+1] > nGlobal)) {
        snprintf(line, l, "caps_checkDiscr: tessGlobal[2*%d+1] = %d [1,%d] for %d!\n",
                 i, discr->tessGlobal[2*i+1], nGlobal, bIndex);
        return CAPS_BADINDEX;
      }
    }
/*@+nullderef@*/

    if (bIndex == 1) {
      if (discBody->globalOffset != 0) {
        snprintf(line, l, "caps_checkDiscr: body %d globalOffset != 0!\n",
                 bIndex);
        return CAPS_BADINDEX;
      }
    }
    if (bIndex < discr->nBodys) {
      if (discr->bodys[bIndex].globalOffset != discBody->globalOffset + nGlobal) {
        snprintf(line, l, "caps_checkDiscr: body %d globalOffset != %d!\n",
                 bIndex+1, discBody->globalOffset + nGlobal);
        return CAPS_BADINDEX;
      }
    }
  }

  /* check data triangulation */
  if ((discr->nDtris != 0) && (discr->dtris != NULL) &&
      (discr->nVerts != 0))
    for (i = 0; i < discr->nDtris; i++)
      for (j = 0; j < 3; j++)
        if ((discr->dtris[3*i+j] < 1) ||
            (discr->dtris[3*i+j] > discr->nVerts)) {
          snprintf(line, l,
                   "caps_checkDiscr: dtris[%d %d] = %d out of range [1-%d] ",
                   i+1, j+1, discr->dtris[3*i+j], discr->nVerts);
          return CAPS_BADINDEX;
        }

  return CAPS_SUCCESS;
}


static void
caps_fillBuiltIn(capsObject *bobject, capsDiscr *discr, capsObject *dobject,
                 CAPSLONG sNum)
{
  int           i, j, k, m, stat, npts, rank, bIndex, pt, pi, typ, len;
  int           last, ntris, nptx, index, global;
  double        *values, xyz[3], st[2];
  capsObject    *pobject, *vobject;
  capsProblem   *problem;
  capsBound     *bound;
  aimInfo       *aInfo;
  capsAnalysis  *analysis;
  capsDataSet   *dataset, *ds;
  capsVertexSet *vertexset;
  capsBodyDiscr *discBody;
  const int     *ptype, *pindex, *tris, *tric;
  const double  *xyzs, *prms;

  aInfo     = (aimInfo *)       discr->aInfo;
  analysis  = (capsAnalysis *)  aInfo->analysis;
  bound     = (capsBound *)     bobject->blind;
  dataset   = (capsDataSet *)   dobject->blind;
  vobject   =                   dobject->parent;
  vertexset = (capsVertexSet *) vobject->blind;
  pobject   =                   bobject->parent;
  problem   = (capsProblem *)   pobject->blind;
  if ((strcmp(dobject->name, "xyzd")   == 0) && (discr->verts == NULL)) return;
  if ((strcmp(dobject->name, "paramd") == 0) && (discr->verts == NULL)) return;

  rank = dataset->rank;
  if (strcmp(dobject->name, "xyz") == 0) {
    npts = discr->nPoints;
  } else if (strcmp(dobject->name, "xyzd") == 0) {
    npts = discr->nVerts;
  } else if (strcmp(dobject->name, "param") == 0) {
    if (bound->state == MultipleError) return;
    npts = discr->nPoints;
  } else if (strcmp(dobject->name, "paramd") == 0) {
    if (bound->state == MultipleError) return;
    npts = discr->nVerts;
  } else {
    printf(" CAPS Internal: Unknown BuiltIn DataSet = %s\n", dobject->name);
    return;
  }
  if (npts == 0) return;

  values = (double *) EG_alloc(npts*rank*sizeof(double));
  if (values == NULL) {
    printf(" CAPS Internal: Malloc on %d %d  Dataset = %s\n", npts, rank,
           dobject->name);
    return;
  }

  if (strcmp(dobject->name, "xyz") == 0) {

    for (i = 0; i < npts; i++) {
      bIndex = discr->tessGlobal[2*i  ];
      global = discr->tessGlobal[2*i+1];
      discBody = &discr->bodys[bIndex-1];
      stat = EG_getGlobal(discBody->tess, global, &pt, &pi, &values[3*i]);
      if (stat != EGADS_SUCCESS)
        printf(" CAPS Internal: %d EG_getGlobal %d for %s = %d\n",
               bIndex, i+1, dobject->name, stat);
    }
    if (bound->lunits != NULL) {
      if (dataset->units != NULL) EG_free(dataset->units);
      dataset->units = EG_strdup(bound->lunits);
    }

  } else if (strcmp(dobject->name, "xyzd") == 0) {

    for (i = 0; i < npts; i++) {
      values[3*i  ] = discr->verts[3*i  ];
      values[3*i+1] = discr->verts[3*i+1];
      values[3*i+2] = discr->verts[3*i+2];
    }
    if (bound->lunits != NULL) {
      if (dataset->units != NULL) EG_free(dataset->units);
      dataset->units = EG_strdup(bound->lunits);
    }

  } else if (strcmp(dobject->name, "param") == 0) {

    if (bound->dim == 2) {
      for (i = 0; i < npts; i++) {
        values[2*i  ] = 0.0;
        values[2*i+1] = 0.0;
      }
      if (bound->state != Multiple) {
        for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
          discBody = &discr->bodys[bIndex-1];
          for (last = j = 0; j < discBody->nElems; j++) {
            if (discBody->elems[j].eIndex != last) {
              stat = EG_getTessFace(discBody->tess, discBody->elems[j].eIndex,
                                    &nptx, &xyzs, &prms,
                                    &ptype, &pindex, &ntris, &tris, &tric);
              if (stat != EGADS_SUCCESS) {
                printf(" CAPS Internal: getTessFace %d = %d for %d\n",
                       discBody->elems[j].eIndex, stat, bIndex);
                continue;
              }
              last = discBody->elems[j].eIndex;
            }
#ifndef __clang_analyzer__
            typ = discBody->elems[j].tIndex;
            len = discr->types[typ-1].nref;
            for (k = 0; k < len; k++) {
              i  = discBody->elems[j].gIndices[2*k  ]-1;
              pt = discBody->elems[j].gIndices[2*k+1]-1;
              values[2*i  ] = prms[2*pt  ];
              values[2*i+1] = prms[2*pt+1];
            }
#endif
          }
        }
      } else {
        for (i = 0; i < npts; i++) {
          bIndex = discr->tessGlobal[2*i  ];
          global = discr->tessGlobal[2*i+1];
          discBody = &discr->bodys[bIndex-1];
          stat = EG_getGlobal(discBody->tess, global, &pt, &pi, xyz);
          if (stat != EGADS_SUCCESS) {
            printf(" CAPS Internal: %d EG_getGlobal %d for %s = %d\n",
                   bIndex, i+1, dobject->name, stat);
            continue;
          }
          stat = caps_invInterpolate2D(bound->surface, xyz, &values[2*i]);
          if (stat != EGADS_SUCCESS)
            printf(" CAPS Internal: caps_invInterpolate2D %d for %s = %d\n",
                   i+1, dobject->name, stat);
        }
      }

    } else {

      if (bound->state == Single) {
        for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
          discBody = &discr->bodys[bIndex-1];
          for (last = j = 0; j < discBody->nElems; j++) {
            if (discBody->elems[j].eIndex != last) {
              stat = EG_getTessEdge(discBody->tess, discBody->elems[j].eIndex,
                                    &nptx, &xyzs, &prms);
              if (stat != EGADS_SUCCESS) {
                printf(" CAPS Internal: getTessEdge %d = %d for %d\n",
                       discBody->elems[j].eIndex, stat, bIndex);
                continue;
              }
              last = discBody->elems[j].eIndex;
            }
#ifndef __clang_analyzer__
            typ = discBody->elems[j].tIndex;
            len = discr->types[typ-1].nref;
            for (k = 0; k < len; k++) {
              i  = discBody->elems[j].gIndices[2*k  ]-1;
              pt = discBody->elems[j].gIndices[2*k+1]-1;
              values[i] = prms[pt];
            }
#endif
          }
        }
      } else {
        for (i = 0; i < npts; i++) {
          bIndex = discr->tessGlobal[2*i  ];
          global = discr->tessGlobal[2*i+1];
          discBody = &discr->bodys[bIndex-1];
          stat = EG_getGlobal(discBody->tess, global, &pt, &pi, xyz);
          if (stat != EGADS_SUCCESS) {
            printf(" CAPS Internal: %d EG_getGlobal %d for %s = %d\n",
                   bIndex, i+1, dobject->name, stat);
            continue;
          }
          stat = caps_invInterpolate1D(bound->curve, xyz, &values[i]);
          if (stat != EGADS_SUCCESS)
            printf(" CAPS Internal: caps_invInterpolate1D %d for %s = %d\n",
                   i+1, dobject->name, stat);
        }
      }
    }

  } else {

    for (i = 0; i < npts; i++) values[i] = 0.0;
    index = aim_Index(problem->aimFPTR, analysis->loadName);
    if (index < 0)
      printf(" CAPS Internal: aim_Index %s = %d\n", dobject->name, index);

    /* assume that we have the dataset "param" */
    ds = NULL;
    if (vertexset->dataSets[2] == NULL) {
      printf(" CAPS Internal: params obj == NULL for %s\n", dobject->name);
      index = -1;
    } else {
      ds = (capsDataSet *) vertexset->dataSets[2]->blind;
      if (ds == NULL) {
        printf(" CAPS Internal: params ds == NULL for %s\n", dobject->name);
        index = -1;
      }
    }

    if (bound->dim == 2) {
      for (i = 0; i < npts; i++) {
        values[2*i  ] = 0.0;
        values[2*i+1] = 0.0;
      }
      if (bound->state == Single) {
        if (index >= 0) {
          for (i = 0; i < npts; i++) {
            bIndex = discr->celem[2*i  ];
            k      = discr->celem[2*i+1] - 1;
            discBody = &discr->bodys[bIndex-1];
            m   = discBody->elems[k].tIndex - 1;
            len = discr->types[m].ndata;
            for (j = 0; j < len; j++)
              if (discBody->elems[k].dIndices[j] == i+1) {
                st[0] = discr->types[m].dst[2*j  ];
                st[1] = discr->types[m].dst[2*j+1];
                break;
              }
            if ((j == len) || (ds == NULL)) {
              printf(" CAPS Internal: data ref %d for %s not found!\n",
                     i+1, dobject->name);
              continue;
            }
            stat  = aim_InterpolIndex(problem->aimFPTR, index, discr, "paramd",
                                      discr->celem[2*i], discr->celem[2*i+1],
                                      st, 2, ds->data, &values[2*i]);
            if (stat != CAPS_SUCCESS)
              printf(" CAPS Internal: aim_InterpolIndex %d for %s = %d\n",
                     i+1, dobject->name, stat);
          }
        }
      } else {
        for (i = 0; i < npts; i++) {
          xyz[0] = discr->verts[3*i  ];
          xyz[1] = discr->verts[3*i+1];
          xyz[2] = discr->verts[3*i+2];
          stat   = caps_invInterpolate2D(bound->surface, xyz, &values[2*i]);
          if (stat != EGADS_SUCCESS)
            printf(" CAPS Internal: caps_invInterpolate2D %d for %s = %d\n",
                   i+1, dobject->name, stat);
        }
      }

    } else {
      for (i = 0; i < npts; i++) values[i] = 0.0;
      if (bound->state == Single) {
        if (index >= 0) {
          for (i = 0; i < discr->nVerts; i++) {
            bIndex = discr->celem[2*i  ];
            k      = discr->celem[2*i+1] - 1;
            discBody = &discr->bodys[bIndex-1];
            m   = discBody->elems[k].tIndex - 1;
            len = discr->types[m].ndata;
            for (j = 0; j < len; j++)
              if (discBody->elems[k].dIndices[j] == i+1) {
                st[0] = discr->types[m].dst[j];
                break;
              }
            if ((j == len) || (ds == NULL)) {
              printf(" CAPS Internal: data ref %d for %s not found!\n",
                     i+1, dobject->name);
              continue;
            }
            stat  = aim_InterpolIndex(problem->aimFPTR, index, discr, "paramd",
                                      discr->celem[2*i], discr->celem[2*i+1],
                                      st, 1, ds->data, &values[i]);
            if (stat != CAPS_SUCCESS)
              printf(" CAPS Internal: aim_InterpolIndex %d for %s = %d\n",
                     i+1, dobject->name, stat);
          }
        }
      } else {
        for (i = 0; i < npts; i++) {
          xyz[0] = discr->verts[3*i  ];
          xyz[1] = discr->verts[3*i+1];
          xyz[2] = discr->verts[3*i+2];
          stat = caps_invInterpolate1D(bound->curve, xyz, &values[i]);
          if (stat != EGADS_SUCCESS)
            printf(" CAPS Internal: caps_invInterpolate1D %d for %s = %d\n",
                   i+1, dobject->name, stat);
        }
      }
    }

  }

  dataset->data = values;
  dataset->npts = npts;
  caps_freeOwner(&dobject->last);
  dobject->last.sNum = sNum;
  stat = caps_addHistory(dobject, problem);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_fillBuiltIn)\n", stat);
  stat = caps_writeDataSet(dobject);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeDataSet = %d (caps_fillBuiltIn)\n", stat);
}


static void
caps_fillSensit(capsProblem *problem, capsDiscr *discr, capsDataSet *dataset)
{
  int           i, j, k, n, bIndex, index, ibody, stat, state, nGlobal, outLevel;
  int           ii, ni, oclass, mtype, nEdge, nFace, len, ntri, *bins, *senses;
  double        limits[4], *vel;
  const int     *ptype, *pindex, *tris, *tric;
  const double  *dxyz, *xyzs, *uvs;
  ego           topRef, prev, next, tess, oldtess, body, eobject, *objs;
  modl_T        *MODL;
  capsBodyDiscr *discBody;

  MODL = (modl_T *) problem->modl;
  for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
    discBody = &discr->bodys[bIndex-1];

    stat = EG_statusTessBody(discBody->tess, &body, &state, &nGlobal);
    if (stat != EGADS_SUCCESS) {
      printf(" caps_fillSensit abort: EG_statusTessBody = %d for body %d!\n",
             stat, bIndex);
      return;
    }

    stat = EG_getInfo(body, &oclass, &mtype, &topRef, &prev, &next);
    if (stat != EGADS_SUCCESS) {
      printf(" caps_fillSensit abort: getInfo = %d for body %d!\n",
             stat, bIndex);
      return;
    }
    nEdge = nFace = 0;
    if (mtype == WIREBODY) {
      stat = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
      if (stat != EGADS_SUCCESS) {
        printf(" caps_fillSensit abort: getBodyTopos (Edge) = %d for %d!\n",
               stat, bIndex);
        return;
      }
    } else {
      stat = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
      if (stat != EGADS_SUCCESS) {
        printf(" caps_fillSensit abort: getBodyTopos (Face) = %d for %d!\n",
               stat, bIndex);
        return;
      }
    }
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
      if (MODL->body[ibody].onstack != 1) continue;
      if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
      if (MODL->body[ibody].ebody   == body) break;
    }
    if (ibody > MODL->nbody) {
      printf(" caps_fillSensit abort: Body Not Found in OpenCSM stack!\n");
      return;
    }
    oldtess = MODL->body[ibody].etess;
    tess    = discBody->tess;
    if (tess == NULL) {
      printf(" caps_fillSensit abort: Body Tess %d Not Found!\n",
             bIndex);
      return;
    }
    MODL->body[ibody].etess = tess;

    bins = (int *) EG_alloc((nEdge+nFace)*sizeof(int));
    if (bins == NULL) {
      printf(" caps_fillSensit abort: %d allocating = %d ints!\n",
             ibody, nEdge+nFace);
      MODL->body[ibody].etess = oldtess;
      return;
    }

    for (ii = 0; ii < discBody->nElems; ii++) {
      j = discBody->elems[ii].eIndex - 1;
      bins[j]++;
    }

    if (nFace == 0) {

      for (index = 1; index <= nEdge; index++) {
        if (bins[index-1] == 0) continue;
        stat = EG_getTessEdge(tess, index, &len, &xyzs, &uvs);
        if (stat != SUCCESS) {
          printf(" caps_fillSensit EG_getTessFace Edge = %d for %d!\n",
                 stat, index);
          continue;
        }
        if (len == 0) continue;
        /* get the EGADS status of the Edge */
        stat = EG_objectBodyTopo(body, EDGE, index, &eobject);
        if (stat != EGADS_SUCCESS) {
          printf(" caps_fillSensit EG_objectBodyTopo Edge = %d for %d!\n",
                 stat, index);
          continue;
        }
        stat = EG_getTopology(eobject, &prev, &oclass, &mtype, limits, &n,
                              &objs, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" caps_fillSensit EG_getTopology Edge = %d for %d!\n",
                 stat, index);
          continue;
        }
        if (mtype == DEGENERATE) {
          printf(" caps_fillSensit DEGENERATE Edge for %d!\n", index);
          continue;
        }
        vel = (double *) EG_alloc(3*len*sizeof(double));
        if (vel == NULL) {
          printf(" caps_fillSensit Edge malloc using %d doubles for %d!\n",
                 len, index);
          continue;
        }
        if (dataset->ftype == TessSens) {
          outLevel = ocsmSetOutLevel(0);
          stat = ocsmGetTessVel(problem->modl, ibody, OCSM_EDGE, index, &dxyz);
          ocsmSetOutLevel(outLevel);
          if (stat != SUCCESS) {
            printf(" caps_fillSensit ocsmGetTessVel Edge = %d for %d!\n",
                   stat, index);
            EG_free(vel);
            continue;
          }
          for (k = 0; k < 3*len; k++) vel[k] = dxyz[k];
          /* readjust the Nodes */
          ii = EG_indexBodyTopo(body, objs[0]);
          if (ii <= EGADS_SUCCESS) {
            printf(" caps_fillSensit EG_indexBodyTopo Edge 0 = %d for %d!\n",
                   ii, index);
            EG_free(vel);
            continue;
          }
          outLevel = ocsmSetOutLevel(0);
          stat = ocsmGetTessVel(problem->modl, ibody, OCSM_NODE, ii, &dxyz);
          ocsmSetOutLevel(outLevel);
          if (stat != SUCCESS) {
            printf(" caps_fillSensit ocsmGetTessVel Node 0 = %d for %d!\n",
                   stat, index);
            EG_free(vel);
            continue;
          }
          vel[0] = dxyz[0];
          vel[1] = dxyz[1];
          vel[2] = dxyz[2];
          if (n > 1) {
            ii = EG_indexBodyTopo(body, objs[1]);
            if (ii <= EGADS_SUCCESS) {
              printf(" caps_fillSensit EG_indexBodyTopo Edge 1 = %d for %d!\n",
                     ii, index);
              EG_free(vel);
              continue;
            }
            outLevel = ocsmSetOutLevel(0);
            stat = ocsmGetTessVel(problem->modl, ibody, OCSM_NODE, ii, &dxyz);
            ocsmSetOutLevel(outLevel);
            if (stat != SUCCESS) {
              printf(" caps_fillSensit ocsmGetTessVel Node 1 = %d for %d!\n",
                     stat, index);
              EG_free(vel);
              continue;
            }
          }
          vel[3*len-3] = dxyz[0];
          vel[3*len-2] = dxyz[1];
          vel[3*len-1] = dxyz[2];
        } else {
          outLevel = ocsmSetOutLevel(0);
          stat = ocsmGetVel(problem->modl, ibody, OCSM_EDGE, index, len, NULL,
                            vel);
          ocsmSetOutLevel(outLevel);
          if (stat != SUCCESS) {
            printf(" caps_fillSensit ocsmGetVel Edge = %d for %d!\n",
                   stat, index);
            EG_free(vel);
            continue;
          }
        }
        for (ii = 0; ii < discBody->nElems; ii++) {
          if (discBody->elems[ii].eIndex != index) continue;
          ni = discr->types[discBody->elems[ii].tIndex-1].nref;
          for (k = 0; k < ni; k++) {
            i = discBody->elems[ii].gIndices[2*k  ]-1;
            j = discBody->elems[ii].gIndices[2*k+1]-1;
            dataset->data[3*i  ] = vel[3*j  ];
            dataset->data[3*i+1] = vel[3*j+1];
            dataset->data[3*i+2] = vel[3*j+2];
          }
        }
        EG_free(vel);
      }

    } else {

      for (index = 1; index <= nFace; index++) {
        if (bins[index-1] == 0) continue;
        stat = EG_getTessFace(tess, index, &len, &xyzs, &uvs, &ptype,
                              &pindex, &ntri, &tris, &tric);
        if (stat != SUCCESS) {
          printf(" caps_fillSensit EG_getTessFace Face = %d for %d!\n",
                 stat, index);
          continue;
        }
        if (len == 0) continue;
        vel = (double *) EG_alloc(3*len*sizeof(double));
        if (vel == NULL) {
          printf(" caps_fillSensit Face malloc using %d doubles for %d!\n",
                 len, index);
          continue;
        }
        if (dataset->ftype == TessSens) {
          outLevel = ocsmSetOutLevel(0);
          stat = ocsmGetTessVel(problem->modl, ibody, OCSM_FACE, index, &dxyz);
          ocsmSetOutLevel(outLevel);
          if (stat != SUCCESS) {
            printf(" caps_fillSensit ocsmGetTessVel Face = %d for %d!\n",
                   stat, index);
            EG_free(vel);
            continue;
          }
          for (k = 0; k < 3*len; k++) vel[k] = dxyz[k];
          /* readjust the non-interior velocities */
          for (state = ni = k = 0; k < len; k++)
            if (ptype[k] == 0) {
              outLevel = ocsmSetOutLevel(0);
              stat = ocsmGetTessVel(problem->modl, ibody, OCSM_NODE, pindex[k],
                                    &dxyz);
              ocsmSetOutLevel(outLevel);
              if (stat != SUCCESS) {
                printf(" caps_fillSensit ocsmGetTessVel Node = %d for %d - %d!\n",
                       stat, pindex[k], index);
                state++;
                break;
              }
              ni         = 0;
              vel[3*k  ] = dxyz[0];
              vel[3*k+1] = dxyz[1];
              vel[3*k+2] = dxyz[2];
            } else if (ptype[k] > 0) {
              if (pindex[k] != ni) {
                outLevel = ocsmSetOutLevel(0);
                stat = ocsmGetTessVel(problem->modl, ibody, OCSM_EDGE, pindex[k],
                                      &dxyz);
                ocsmSetOutLevel(outLevel);
                if (stat != SUCCESS) {
                  printf(" caps_fillSensit ocsmGetTessVel Edge = %d for %d - %d!\n",
                         stat, pindex[k], index);
                  state++;
                  break;
                }
                ni = pindex[k];
              }
              vel[3*k  ] = dxyz[3*ptype[k]-3];
              vel[3*k+1] = dxyz[3*ptype[k]-2];
              vel[3*k+2] = dxyz[3*ptype[k]-1];
            }
          if (state != 0) {
            EG_free(vel);
            continue;
          }
        } else {
          outLevel = ocsmSetOutLevel(0);
          stat = ocsmGetVel(problem->modl, ibody, OCSM_FACE, index, len, NULL,
                            vel);
          ocsmSetOutLevel(outLevel);
          if (stat != SUCCESS) {
            printf(" caps_fillSensit ocsmGetVel Face = %d for %d!\n",
                   stat, index);
            EG_free(vel);
            continue;
          }
        }
        for (ii = 0; ii < discBody->nElems; ii++) {
          if (discBody->elems[ii].eIndex != index) continue;
          ni = discr->types[discBody->elems[ii].tIndex-1].nref;
          for (k = 0; k < ni; k++) {
            i = discBody->elems[ii].gIndices[2*k  ]-1;
            j = discBody->elems[ii].gIndices[2*k+1]-1;
            dataset->data[3*i  ] = vel[3*j  ];
            dataset->data[3*i+1] = vel[3*j+1];
            dataset->data[3*i+2] = vel[3*j+2];
          }
        }
        EG_free(vel);
      }

    }
    EG_free(bins);

    MODL->body[ibody].etess = oldtess;
  }

}


static void
caps_freeBodyObjs(int nBodies, /*@only@*/ bodyObjs *bodies)
{
  int i;

  for (i = 0; i < nBodies; i++) {
    if (bodies[i].objs    != NULL) EG_free(bodies[i].objs);
    if (bodies[i].indices != NULL) EG_free(bodies[i].indices);
  }
  EG_free(bodies);
}


static double
caps_triangleArea3D(const double *xyz0, const double *xyz1, const double *xyz2)
{
  double x1[3], x2[3], n[3];

  x1[0] = xyz1[0] - xyz0[0];
  x2[0] = xyz2[0] - xyz0[0];
  x1[1] = xyz1[1] - xyz0[1];
  x2[1] = xyz2[1] - xyz0[1];
  x1[2] = xyz1[2] - xyz0[2];
  x2[2] = xyz2[2] - xyz0[2];
  CROSS(n, x1, x2);
  return 0.5*sqrt(DOT(n, n));
}


static int
caps_paramQuilt(capsBound *bound, int l, char *line)
{
  int           i, j, k, m, n, stat, npts, ntris, own, nu, nv, per, eType;
  int           ntrx, nptx, last, bIndex, pt, pi, *ppnts, *vtab, count, iVS;
  int           i0, i1, i2, global;
  double        coord[3], box[6], tol, rmserr, maxerr, dotmin, area, d;
  double        *grid, *r, *xyzs;
  capsDiscr     *quilt;
  capsAprx2D    *surface;
  capsVertexSet *vertexset;
  capsBodyDiscr *discBody;
  prmXYZ        *xyz;
  prmTri        *tris;
  prmUVF        *uvf;
  prmUV         *uv;
  connect       *etab;
  const int     *ptype, *pindex, *trix, *tric;
  const double  *xyzx, *prms;

  for (npts = ntris = i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL)      continue;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->type        != VERTEXSET) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    if (vertexset->analysis              == NULL)      continue;
    if (vertexset->discr                 == NULL)      continue;
    quilt  = vertexset->discr;

    for (bIndex = 1; bIndex <= quilt->nBodys; bIndex++) {
      discBody = &quilt->bodys[bIndex-1];
      for (j = 0; j < discBody->nElems; j++) {
        eType = discBody->elems[j].tIndex;
        if (quilt->types[eType-1].tris == NULL) {
          ntris++;
        } else {
          ntris += quilt->types[eType-1].ntri;
        }
      }
    }
    npts += quilt->nPoints;
  }

  if ((ntris == 0) || (npts == 0)) {
    snprintf(line, l,
             "caps_paramQuilt Error: nPoints = %d  nTris = %d", npts, ntris);
    return CAPS_NOTCONNECT;
  }

  uv    = (prmUV *)  EG_alloc(npts*sizeof(prmUV));
  if (uv   == NULL) {
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on = %d prmUV", npts);
    return EGADS_MALLOC;
  }
  uvf   = (prmUVF *) EG_alloc(ntris*sizeof(prmUVF));
  if (uvf  == NULL) {
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on = %d prmUVF", ntris);
    return EGADS_MALLOC;
  }
  xyz   = (prmXYZ *) EG_alloc(npts*sizeof(prmXYZ));
  if (xyz  == NULL) {
    EG_free(uvf);
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on = %d prmXYZ", npts);
    return EGADS_MALLOC;
  }
  xyzs  = (double *) xyz;
  tris  = (prmTri *) EG_alloc(ntris*sizeof(prmTri));
  if (tris  == NULL) {
    EG_free(xyz);
    EG_free(uvf);
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on = %d ints", 3*ntris);
    return EGADS_MALLOC;
  }

  /* find the best candidate VertexSet for fitting */
  iVS  = -1;
  area =  0.0;
  for (i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL)      continue;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->type        != VERTEXSET) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    if (vertexset->analysis              == NULL)      continue;
    if (vertexset->discr                 == NULL)      continue;
    quilt    = vertexset->discr;

    d        = 0.0;
    ntris    = 0;
    for (bIndex = 1; bIndex <= quilt->nBodys; bIndex++) {
      discBody = &quilt->bodys[bIndex-1];
      for (last = j = 0; j < discBody->nElems; j++) {
        eType = discBody->elems[j].tIndex;
        own   = discBody->elems[j].eIndex;
        if (own != last) {
          stat = EG_getTessFace(discBody->tess,
                                own, &nptx, &xyzx, &prms, &ptype, &pindex,
                                &ntrx, &trix, &tric);
          if (stat != EGADS_SUCCESS) {
            printf(" CAPS Internal: EG_getTessFace %d = %d for %d\n",
                   own, stat, bIndex);
            continue;
          }
/*        printf(" VS = %s\n   Face %d: ntris = %d  quilt eles = %d!\n",
                 bound->vertexSet[i]->parent->name, own, nptx, quilt->nElems); */
          last = own;
        }
        if (quilt->types[eType-1].tris == NULL) {
#ifndef __clang_analyzer__
          i0 = discBody->elems[j].gIndices[1] - 1;
          i1 = discBody->elems[j].gIndices[3] - 1;
          i2 = discBody->elems[j].gIndices[5] - 1;
          d += caps_triangleArea3D(&xyzx[3*i0], &xyzx[3*i1], &xyzx[3*i2]);
#endif
          ntris++;
        } else {
#ifndef __clang_analyzer__
          for (k = 0; k < quilt->types[eType-1].ntri; k++, ntris++) {
            n  = quilt->types[eType-1].tris [3*k  ] - 1;
            i0 = discBody->elems[j].gIndices[2*n+1] - 1;
            n  = quilt->types[eType-1].tris [3*k+1] - 1;
            i1 = discBody->elems[j].gIndices[2*n+1] - 1;
            n  = quilt->types[eType-1].tris [3*k+2] - 1;
            i2 = discBody->elems[j].gIndices[2*n+1] - 1;
            d += caps_triangleArea3D(&xyzx[3*i0], &xyzx[3*i1], &xyzx[3*i2]);
          }
#endif
        }
      }
    }
#ifdef DEBUG
    printf(" VertexSet %d: area = %lf  ntris = %d\n", i+1, d, ntris);
#endif
    if (d > area) {
      iVS  = i;
      area = d;
    }
  }
  if (iVS == -1) {
    EG_free(tris);
    EG_free(xyz);
    EG_free(uvf);
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: No VertexSet Selected!");
    return EGADS_NOTFOUND;
  }
#ifdef DEBUG
  printf(" selected VertexSet = %d\n", iVS+1);
#endif

/*
  for (count = npts = ntris = i = 0; i < bound->nVertexSet; i++) {  */
  count = npts = ntris = 0;
  for (i = iVS; i <= iVS; i++) {
    if (bound->vertexSet[i]              == NULL)      continue;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->type        != VERTEXSET) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    if (vertexset->analysis              == NULL)      continue;
    if (vertexset->discr                 == NULL)      continue;
    quilt    = vertexset->discr;
    prms     = NULL;
    for (bIndex = 1; bIndex <= quilt->nBodys; bIndex++) {
      discBody = &quilt->bodys[bIndex-1];
      for (last = j = 0; j < discBody->nElems; j++) {
        eType = discBody->elems[j].tIndex;
        own   = discBody->elems[j].eIndex;
        if (own != last) {
          count++;
/*        printf(" Analysis %d/%d Body %d/%d Face %d: count = %d, npts = %d\n",
                 i+1, bound->nVertexSet, bIndex, analysis->nBody, own, count,
                 npts);  */
          stat = EG_getTessFace(discBody->tess,
                                own, &nptx, &xyzx, &prms, &ptype, &pindex,
                                &ntrx, &trix, &tric);
          if (stat != EGADS_SUCCESS) {
            printf(" CAPS Internal: getTessFace %d = %d for %d\n",
                   own, stat, bIndex);
            continue;
          }
          last = own;
        }
        if (prms == NULL) continue;
        if (quilt->types[eType-1].tris == NULL) {
          m = discBody->elems[j].gIndices[1] - 1;
          uvf[ntris].u0          = prms[2*m  ];
          uvf[ntris].v0          = prms[2*m+1];
          tris[ntris].indices[0] = discBody->elems[j].gIndices[0] + npts;
          m = discBody->elems[j].gIndices[3] - 1;
          uvf[ntris].u1          = prms[2*m  ];
          uvf[ntris].v1          = prms[2*m+1];
          tris[ntris].indices[1] = discBody->elems[j].gIndices[2] + npts;
          m = discBody->elems[j].gIndices[5] - 1;
          uvf[ntris].u2          = prms[2*m  ];
          uvf[ntris].v2          = prms[2*m+1];
          tris[ntris].indices[2] = discBody->elems[j].gIndices[4] + npts;
          tris[ntris].own        = count;
          ntris++;
        } else {
          for (k = 0; k < quilt->types[eType-1].ntri; k++, ntris++) {
            n = quilt->types[eType-1].tris [3*k  ] - 1;
            m = discBody->elems[j].gIndices[2*n+1] - 1;
            uvf[ntris].u0          = prms[2*m  ];
            uvf[ntris].v0          = prms[2*m+1];
            tris[ntris].indices[0] = discBody->elems[j].gIndices[2*n] + npts;
            n = quilt->types[eType-1].tris [3*k+1] - 1;
            m = discBody->elems[j].gIndices[2*n+1] - 1;
            uvf[ntris].u1          = prms[2*m  ];
            uvf[ntris].v1          = prms[2*m+1];
            tris[ntris].indices[1] = discBody->elems[j].gIndices[2*n] + npts;
            n = quilt->types[eType-1].tris [3*k+2] - 1;
            m = discBody->elems[j].gIndices[2*n+1] - 1;
            uvf[ntris].u2          = prms[2*m  ];
            uvf[ntris].v2          = prms[2*m+1];
            tris[ntris].indices[2] = discBody->elems[j].gIndices[2*n] + npts;
            tris[ntris].own        = count;
          }
        }
      }
    }
    for (j = 0; j < quilt->nPoints; j++) {
      bIndex = quilt->tessGlobal[2*j  ];
      global = quilt->tessGlobal[2*j+1];
      discBody = &quilt->bodys[bIndex-1];
      stat = EG_getGlobal(discBody->tess, global, &pt, &pi, coord);
      if (stat != EGADS_SUCCESS) {
        printf(" CAPS Internal: %d EG_getGlobal %d = %d\n", bIndex, j+1, stat);
      } else {
        xyz[npts+j].x = coord[0];
        xyz[npts+j].y = coord[1];
        xyz[npts+j].z = coord[2];
      }
    }
    npts += quilt->nPoints;
  }

  /* make the neighbors */

  vtab = (int *) EG_alloc(npts*sizeof(int));
  if (vtab == NULL) {
    EG_free(tris);
    EG_free(xyz);
    EG_free(uvf);
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on table = %d ints", npts);
    return EGADS_MALLOC;
  }
  etab = (connect *) EG_alloc(ntris*3*sizeof(connect));
  if (etab == NULL) {
    EG_free(vtab);
    EG_free(tris);
    EG_free(xyz);
    EG_free(uvf);
    EG_free(uv);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on side table = %d connect",
             3*ntris);
    return EGADS_MALLOC;
  }
  n = NOTFILLED;
  for (j = 0; j < npts; j++) vtab[j] = NOTFILLED;
  for (i = 0; i < ntris;  i++) {
    tris[i].neigh[0] = tris[i].neigh[1] = tris[i].neigh[2] = i+1;
    EG_makeConnect( tris[i].indices[1], tris[i].indices[2],
                   &tris[i].neigh[0], &n, vtab, etab, 0);
    EG_makeConnect( tris[i].indices[0], tris[i].indices[2],
                   &tris[i].neigh[1], &n, vtab, etab, 0);
    EG_makeConnect( tris[i].indices[0], tris[i].indices[1],
                   &tris[i].neigh[2], &n, vtab, etab, 0);
  }
  /* find any unconnected triangle sides */
  for (j = 0; j <= n; j++) {
    if (etab[j].tri == NULL) continue;
/*  printf(" CAPS Info: Unconnected Side %d %d = %d\n",
           etab[j].node1+1, etab[j].node2+1, *etab[j].tri); */
    *etab[j].tri = 0;
  }
  EG_free(etab);
  EG_free(vtab);

  /* get tolerance */

  box[0] = box[3] = xyz[0].x;
  box[1] = box[4] = xyz[0].y;
  box[2] = box[5] = xyz[0].z;
  for (j = 1; j < npts; j++) {
    if (xyz[j].x < box[0]) box[0] = xyz[j].x;
    if (xyz[j].x > box[3]) box[3] = xyz[j].x;
    if (xyz[j].y < box[1]) box[1] = xyz[j].y;
    if (xyz[j].y > box[4]) box[4] = xyz[j].y;
    if (xyz[j].z < box[2]) box[2] = xyz[j].z;
    if (xyz[j].z > box[5]) box[5] = xyz[j].z;
  }
  n     = 1;
  grid  = NULL;
  ppnts = NULL;
  tol   = 1.e-7*sqrt((box[3]-box[0])*(box[3]-box[0]) +
                     (box[4]-box[1])*(box[4]-box[1]) +
                     (box[5]-box[2])*(box[5]-box[2]));

  /* reparameterize */

  stat  = prm_CreateUV(0, ntris, tris, uvf, npts, NULL, NULL, uv, xyz,
                       &per, &ppnts);
#ifdef DEBUG
  printf(" caps_paramQuilt: prm_CreateUV = %d  per = %d\n", stat, per);
#endif
  if (stat > 0) {
    n    = 2;
    stat = prm_SmoothUV(3, per, ppnts, ntris, tris, npts, 3, uv, xyzs);
#ifdef DEBUG
    printf(" caps_paraQuilt: prm_SmoothUV = %d\n", stat);
#endif
    if (stat == CAPS_SUCCESS) {
      n    = 3;
      stat = prm_NormalizeUV(0.05, per, npts, uv);
#ifdef DEBUG
      printf(" caps_paraQuilt: prm_NormalizeUV = %d\n", stat);
#endif
      if (stat == CAPS_SUCCESS) {
        n    = 4;
        nu   = 2*npts;
        nv   = 0;
        stat = prm_BestGrid(npts, 3, uv, xyzs, ntris, tris, tol, per, ppnts,
                            &nu, &nv, &grid, &rmserr, &maxerr, &dotmin);
        if (stat == PRM_TOLERANCEUNMET) {
          printf(" caps_paramQuilt: Tolerance not met: %lf (%lf)!\n",
                 maxerr, tol);
          stat = CAPS_SUCCESS;
        }
#ifdef DEBUG
        printf(" caps_paramQuilt: prm_BestGrid = %d  %d %d  %lf %lf (%lf)\n",
               stat, nu, nv, rmserr, maxerr, tol);
#endif
      }
    }
  }
  if (ppnts != NULL) EG_free(ppnts);
  EG_free(tris);
  EG_free(uvf);
  EG_free(xyz);
  EG_free(uv);
  if ((stat != CAPS_SUCCESS) || (grid == NULL)) {
    snprintf(line, l,
             "caps_paramQuilt: Create/Smooth/Normalize/BestGrid %d = %d!",
             n, stat);
    return stat;
  }

  /* make the surface approximation */

  surface = (capsAprx2D *) EG_alloc(sizeof(capsAprx2D));
  if (surface == NULL) {
    EG_free(grid);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on Surface!");
    return EGADS_MALLOC;
  }
  surface->nrank     = 3;
  surface->periodic  = per;
  surface->nus       = nu;
  surface->nvs       = nv;
  surface->urange[0] = 0.0;
  surface->urange[1] = nu-1;
  surface->vrange[0] = 0.0;
  surface->vrange[1] = nv-1;
  surface->num       = 0;
  surface->nvm       = 0;
  surface->uvmap     = NULL;
  surface->interp    = (double *) EG_alloc(3*4*nu*nv*sizeof(double));
  if (surface->interp == NULL) {
    EG_free(surface);
    EG_free(grid);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on Approx nu = %d, nv = %d",
             nu, nv);
    return EGADS_MALLOC;
  }
  n = nu;
  if (nv > nu) n = nv;
  r = (double *) EG_alloc(6*n*sizeof(double));
  if (r == NULL) {
    caps_Aprx2DFree(surface);
    EG_free(grid);
    snprintf(line, l,
             "caps_paramQuilt Error: Malloc on temp, size = %d", n);
    return EGADS_MALLOC;
  }
  stat = caps_fillCoeff2D(3, nu, nv, grid, surface->interp, r);
  EG_free(r);
  EG_free(grid);
  if (stat == 1) {
    caps_Aprx2DFree(surface);
    snprintf(line, l,
             "caps_paramQuilt Error: Failure in producing interpolant!");
    return CAPS_PARAMBNDERR;
  }
  bound->surface    = surface;
  bound->plimits[0] = bound->plimits[2] = 0.0;
  bound->plimits[1] = nu-1;
  bound->plimits[3] = nv-1;

  return CAPS_SUCCESS;
}


static int
caps_parameterize(capsProblem *problem, capsObject *bobject, int l, char *line)
{
  int           i, j, k, stat, last, top, bIndex, state, nGlobal;
  char          *units;
  double        plims[4];
  ego           entity, body;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsDiscr     *discr;
  capsBodyDiscr *discBody;
  bodyObjs      *bodies;

  bound = (capsBound *) bobject->blind;
  if ((bound->dim != 1) && (bound->dim != 2)) {
    snprintf(line, l, "caps_parameterize Error: Dimension = %d", bound->dim);
    return CAPS_BADINDEX;
  }
  bodies = (bodyObjs *) EG_alloc(problem->nBodies*sizeof(bodyObjs));
  if (bodies == NULL) {
    snprintf(line, l, "caps_parameterize Error: Allocating %d bodyObjs",
             problem->nBodies);
    return EGADS_MALLOC;
  }
  for (i = 0; i < problem->nBodies; i++) {
    bodies[i].n       = 0;
    bodies[i].gIndex  = 0;
    bodies[i].geom    = NULL;
    bodies[i].objs    = NULL;
    bodies[i].indices = NULL;
  }
  for (i = 0; i < problem->nBodies; i++) {
    if (bound->dim == 1) {
      stat = EG_getBodyTopos(problem->bodies[i], NULL, EDGE, &bodies[i].n,
                             &bodies[i].objs);
      if (stat != EGADS_SUCCESS) {
        snprintf(line, l,
                 "caps_parameterize Error: getBodyTopos EDGE for Body #%d",
                 i+1);
        caps_freeBodyObjs(problem->nBodies, bodies);
        return stat;
      }
    } else {
      stat = EG_getBodyTopos(problem->bodies[i], NULL, FACE, &bodies[i].n,
                             &bodies[i].objs);
      if (stat != EGADS_SUCCESS) {
        snprintf(line, l,
                 "caps_parameterize Error: getBodyTopos FACE for Body #%d",
                 i+1);
        caps_freeBodyObjs(problem->nBodies, bodies);
        return stat;
      }
    }
    if (bodies[i].n != 0) {
      bodies[i].indices = (int *) EG_alloc(bodies[i].n*sizeof(int));
      if (bodies[i].indices == NULL) {
        snprintf(line, l, "caps_parameterize Error: malloc %d ints for Body #%d",
                 bodies[i].n, i+1);
        caps_freeBodyObjs(problem->nBodies, bodies);
        return EGADS_MALLOC;
      }
      for (j = 0; j < bodies[i].n; j++) bodies[i].indices[j] = 0;
    }
  }

  for (i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL)      continue;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->type        != VERTEXSET) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    if (vertexset->analysis              == NULL)      continue;
    if (vertexset->discr                 == NULL)      continue;
    discr = vertexset->discr;
    if (discr->dim != bound->dim) {
      snprintf(line, l,
               "caps_parameterize Error: VertexSet %s Dimension = %d not %d",
               bound->vertexSet[i]->name, discr->dim, bound->dim);
      caps_freeBodyObjs(problem->nBodies, bodies);
      return CAPS_BADINDEX;
    }
    for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
      discBody = &discr->bodys[bIndex-1];

      stat = EG_statusTessBody(discBody->tess, &body, &state, &nGlobal);
      if (stat != EGADS_SUCCESS) {
        printf(" caps_parameterize: EG_statusTessBody = %d for body %d!\n",
               stat, bIndex);
        caps_freeBodyObjs(problem->nBodies, bodies);
        return stat;
      }

      for (k = 0; k < problem->nBodies; k++)
        if (problem->bodies[k] == body) {
          for (j = 0; j < discBody->nElems; j++) {
            bodies[k].indices[discBody->elems[j].eIndex-1]++;
            break;
          }
        }
    }
  }

  /* get the parameterization */
  entity            = NULL;
  bound->plimits[0] = bound->plimits[1] = bound->plimits[2] =
                      bound->plimits[3] = 0.0;

  /* do we have only a single Geometry entity -- in body? */
  for (last = i = 0; i < problem->nBodies; i++) {
    for (j = 0; j < bodies[i].n; j++) {
      if (bodies[i].indices[j] != 0) {
        if (bodies[i].gIndex == 0) {
          bodies[i].gIndex    = j+1;
          bodies[i].geom      = bodies[i].objs[j];
          plims[0] = plims[1] = plims[2] = plims[3] = 0.0;
          stat = EG_getRange(bodies[i].geom, plims, &top);
          if (stat != EGADS_SUCCESS) {
            snprintf(line, l,
                     "caps_parameterize Error: getRange for Body #%d %d",
                     i+1, bodies[i].gIndex);
            caps_freeBodyObjs(problem->nBodies, bodies);
            return stat;
          }
          if (last == 0) {
            bound->plimits[0] = plims[0];
            bound->plimits[1] = plims[1];
            bound->plimits[2] = plims[2];
            bound->plimits[3] = plims[3];
            entity            = bodies[i].geom;
            last++;
          } else {
            if (plims[0] < bound->plimits[0]) bound->plimits[0] = plims[0];
            if (plims[1] > bound->plimits[1]) bound->plimits[1] = plims[1];
            if (plims[2] < bound->plimits[2]) bound->plimits[2] = plims[2];
            if (plims[3] > bound->plimits[3]) bound->plimits[3] = plims[3];
          }
        } else {
          stat = EG_isSame(bodies[i].geom, bodies[i].objs[j]);
          if (stat < 0) {
            snprintf(line, l,
                     "caps_parameterize Error: isSame for Body #%d %d %d",
                     i+1, bodies[i].gIndex, j+1);
            caps_freeBodyObjs(problem->nBodies, bodies);
            return stat;
          }
          if (stat != EGADS_SUCCESS) {
            bodies[i].gIndex = -1;
            break;
          }
          plims[0] = plims[1] = plims[2] = plims[3] = 0.0;
          stat = EG_getRange(bodies[i].objs[j], plims, &top);
          if (stat != EGADS_SUCCESS) {
            snprintf(line, l,
                     "caps_parameterize Error: getRange for Body #%d %d",
                     i+1, j+1);
            caps_freeBodyObjs(problem->nBodies, bodies);
            return stat;
          }
          if (plims[0] < bound->plimits[0]) bound->plimits[0] = plims[0];
          if (plims[1] > bound->plimits[1]) bound->plimits[1] = plims[1];
          if (plims[2] < bound->plimits[2]) bound->plimits[2] = plims[2];
          if (plims[3] > bound->plimits[3]) bound->plimits[3] = plims[3];
        }
      }
    }
  }

  /* cross body */
  for (last = i = 0; i < problem->nBodies; i++) {
    if (bodies[i].gIndex ==  0) continue;
    if (bodies[i].gIndex == -1) {
      last = -1;
      break;
    }
    if (last == 0) {
      last = i+1;
    } else {
      stat = EG_isSame(bodies[last-1].geom, bodies[i].geom);
      if (stat < 0) {
        snprintf(line, l,
                 "caps_parameterize Error: isSame for Body #%d %d - #%d %d",
                 last, bodies[last-1].gIndex, i+1, bodies[i].gIndex);
        caps_freeBodyObjs(problem->nBodies, bodies);
        return stat;
      }
      if (stat != EGADS_SUCCESS) {
        last = -1;
        break;
      }
    }
  }
  if (last == 0) {
    bound->state = Empty;
    printf(" CAPS Info: No geometry for Bound -> %s!\n", bobject->name);
    caps_freeBodyObjs(problem->nBodies, bodies);
    return CAPS_SUCCESS;
  }

  /* single geometric entity */
  if (last != -1) {
    if (bound->lunits != NULL) EG_free(bound->lunits);
    bound->geom   = entity;
    bound->iBody  = last;
    bound->iEnt   = bodies[last-1].gIndex;
    bound->state  = Single;
    bound->lunits = EG_strdup(problem->lunits[last-1]);
    caps_freeBodyObjs(problem->nBodies, bodies);
    return CAPS_SUCCESS;
  }

  /* need to reparameterize */
  units = NULL;
  for (i = 0; i < problem->nBodies; i++) {
    if (bodies[i].gIndex == 0) continue;
    if (units == NULL) {
      units = problem->lunits[i];
      continue;
    }
    /* what do we do about units that don't match? */
    if (strcmp(problem->lunits[i], units) != 0) {
      printf(" CAPS Info: Units don't match for Bound -> %s -- %s %s\n",
             bobject->name, problem->lunits[i], units);
    }
  }
  if (bound->lunits != NULL) EG_free(bound->lunits);
  bound->geom   = NULL;
  bound->iBody  = 0;
  bound->iEnt   = 0;
  bound->state  = Multiple;
  bound->lunits = EG_strdup(units);
  caps_freeBodyObjs(problem->nBodies, bodies);

  if (bound->dim == 1) {

    stat = CAPS_SUCCESS;

  } else {
    stat = caps_paramQuilt(bound, l, line);
    if (stat == CAPS_SUCCESS)
      printf(" CAPS Info: Reparameterization Bound -> %s -- nu, nv = %d %d\n",
             bobject->name, bound->surface->nus, bound->surface->nvs);
  }
  if (stat != CAPS_SUCCESS) bound->state = MultipleError;

  return stat;
}


void
caps_geomOutSensit(capsProblem *problem, int ipmtr, int irow, int icol)
{
  int       i, j, k, m, n, len_wrt, i_wrt;
  double    *reals;
  capsValue *value;

  for (i = 0; i < problem->nRegGIN; i++) {
    if (problem->geomIn[problem->regGIN[i].index-1] == NULL) continue;
    value = (capsValue *) problem->geomIn[problem->regGIN[i].index-1]->blind;
    if (value                   == NULL)  continue;
    if (value->pIndex           != ipmtr) continue;
    if (problem->regGIN[i].irow != irow && problem->regGIN[i].irow > 0) continue;
    if (problem->regGIN[i].icol != icol && problem->regGIN[i].icol > 0) continue;

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
        value->derivs[i].deriv = (double *) EG_alloc(value->length*len_wrt*
                                                     sizeof(double));
      if (value->derivs[i].deriv == NULL) continue;

      reals = value->vals.reals;
      if (value->length == 1) reals = &value->vals.real;
      for (n = k = 0; k < value->nrow; k++)
        for (m = 0; m < value->ncol; m++, n++)
          ocsmGetValu(problem->modl, value->pIndex, k+1, m+1,
                      &reals[n], &value->derivs[i].deriv[len_wrt*n + i_wrt]);
    }
  }
}


int
caps_transferLinks(capsAnalysis *analysis, int *nErr, capsErrs **errors)
{
  int              i, status;
  char             temp[PATH_MAX];
  enum capstMethod method;
  capsValue        *value;
  capsObject       *source, *object, *last;

  /* fill in any values that have links */
  for (i = 0; i < analysis->nAnalysisIn; i++) {
    source = object = analysis->analysisIn[i];
    do {
      if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (source->type        != VALUE)     return CAPS_BADTYPE;
      if (source->blind       == NULL)      return CAPS_NULLBLIND;
      value = (capsValue *) source->blind;
      if (value->link == object)            return CAPS_CIRCULARLINK;
      last   = source;
      source = value->link;
    } while (value->link != NULL);
    if (last != object) {
      value  = (capsValue *) object->blind;
      source = value->link;
      method = value->linkMethod;
      status = caps_transferValueX(last, method, object, nErr, errors);
      value->link       = source;
      value->linkMethod = method;
      if (status != CAPS_SUCCESS) {
        snprintf(temp, PATH_MAX, "transferValues for %s from %s = %d",
                 object->name, source->name, status);
        caps_makeSimpleErr(object, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return status;
      }
      caps_freeOwner(&object->last);
      object->last       = last->last;
      object->last.pname = EG_strdup(last->last.pname);
      object->last.pID   = EG_strdup(last->last.pID);
      object->last.user  = EG_strdup(last->last.user);
    }
  }

  return CAPS_SUCCESS;
}


static int
caps_refillBound(capsProblem *problem, capsObject *bobject,
                 int *nErr, capsErrs **errors)
{
  int           i, j, k, m, n, irow, icol, status, len, outLevel;
  int           buildTo, builtTo, nbody, open;
  char          name[129], error[129], *str, **names, **ntmp;
  capsValue     *value;
  capsAnalysis  *anal;
  capsObject    *dobject;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;
  capsDiscr     *discr;
  capsErrs      *errs;
  modl_T        *MODL;

  /* invalidate/cleanup any geometry dependencies & remake the bound */
  bound = (capsBound *) bobject->blind;
  caps_Aprx1DFree(bound->curve);   bound->curve = NULL;
  caps_Aprx2DFree(bound->surface); bound->surface = NULL;

  for (j = 0; j < bound->nVertexSet; j++) {
    if (bound->vertexSet[j]              == NULL)      continue;
    if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[j]->type        != VERTEXSET) continue;
    if (bound->vertexSet[j]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;

    /* remove existing DataSets in the VertexSet */
    for (k = 0; k < vertexset->nDataSets; k++) {
      if (vertexset->dataSets[k]              == NULL)      continue;
      if (vertexset->dataSets[k]->magicnumber != CAPSMAGIC) continue;
      if (vertexset->dataSets[k]->type        != DATASET)   continue;
      if (vertexset->dataSets[k]->blind       == NULL)      continue;
      dobject = vertexset->dataSets[k];
      dataset = (capsDataSet *) dobject->blind;
      if ((dataset->ftype == User) &&
          (strcmp(dobject->name, "xyz") == 0)) continue;
      if (dataset->data != NULL) EG_free(dataset->data);
      dataset->npts = 0;
      dataset->data = NULL;
    }

    if (vertexset->analysis != NULL) {
      if (vertexset->analysis->blind != NULL) {
        anal = (capsAnalysis *) vertexset->analysis->blind;
        aim_FreeDiscr(problem->aimFPTR, anal->loadName, vertexset->discr);

        /* fill in any values that have links */
        status = caps_transferLinks(anal, nErr, errors);
        if (status != CAPS_SUCCESS) return status;

        status = aim_Discr(problem->aimFPTR, anal->loadName,
                           bobject->name, vertexset->discr);
        if (status != CAPS_SUCCESS) {
          aim_FreeDiscr(problem->aimFPTR, anal->loadName, vertexset->discr);
          snprintf(error, 129, "Bound = %s and Analysis = %s",
                   bobject->name, anal->loadName);
          caps_makeSimpleErr(bound->vertexSet[j], CERROR,
                             "caps_refillBound Error: aimDiscr fails!",
                             error, NULL, errors);
          if (*errors != NULL) {
            errs  = *errors;
            *nErr = errs->nError;
          }
          return status;
        }

        /* check the validity of the discretization just returned */
        status = caps_checkDiscr(vertexset->discr, 129, name);
        if (status != CAPS_SUCCESS) {
          snprintf(error, 129, "Bound = %s and Analysis = %s",
                   bobject->name, anal->loadName);
          caps_makeSimpleErr(bound->vertexSet[j], CERROR, name, error,
                             NULL, errors);
          if (*errors != NULL) {
            errs  = *errors;
            *nErr = errs->nError;
          }
          aim_FreeDiscr(problem->aimFPTR, anal->loadName, vertexset->discr);
          return status;
        }
        caps_freeOwner(&bound->vertexSet[j]->last);
        bound->vertexSet[j]->last.sNum = problem->sNum;
        status = caps_addHistory(bound->vertexSet[j], problem);
        if (status != CAPS_SUCCESS)
          printf(" CAPS Warning: caps_addHistory = %d (caps_refillBound)\n",
                 status);
        status = caps_writeVertexSet(bound->vertexSet[j]);
        if (status != CAPS_SUCCESS)
          printf(" CAPS Warning: caps_writeVertexSet = %d (caps_refillBound)\n",
                 status);
      }
    }
  }
  /* reparameterize the existing bounds (dim=1&2) for multiple entities */
  if (bound->dim != 3) {
    status = caps_parameterize(problem, bobject, 129, name);
    if (status != CAPS_SUCCESS) {
      snprintf(error, 129, "Bound = %s", bobject->name);
      caps_makeSimpleErr(bobject, CERROR,
                         "caps_refillBound: Bound Parameterization fails!",
                         error, NULL, errors);
      if (*errors != NULL) {
        errs  = *errors;
        *nErr = errs->nError;
      }
      return status;
    }
  }

  /* populate any built-in DataSet entries after the parameterization */
  for (j = 0; j < bound->nVertexSet; j++) {
    if (bound->vertexSet[j]              == NULL)      continue;
    if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[j]->type        != VERTEXSET) continue;
    if (bound->vertexSet[j]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
    for (k = 0; k < vertexset->nDataSets; k++) {
      if (vertexset->dataSets[k]                 == NULL)      continue;
      if (vertexset->dataSets[k]->magicnumber    != CAPSMAGIC) continue;
      if (vertexset->dataSets[k]->type           != DATASET)   continue;
      if (vertexset->dataSets[k]->blind          == NULL)      continue;
      dataset = (capsDataSet *) vertexset->dataSets[k]->blind;
      if (dataset->ftype                         != BuiltIn)   continue;
      caps_fillBuiltIn(bobject, vertexset->discr, vertexset->dataSets[k],
                       problem->sNum);
    }
  }

  caps_freeOwner(&bobject->last);
  bobject->last.sNum = problem->sNum;
  status = caps_addHistory(bobject, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_refillBound)\n", status);
  status = caps_writeBound(bobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeBound = %d (caps_refillBound)\n", status);

  /* populate any sensitivities in DataSets */
  names = NULL;
  for (n = j = 0; j < bound->nVertexSet; j++) {
    if (bound->vertexSet[j]              == NULL)      continue;
    if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[j]->type        != VERTEXSET) continue;
    if (bound->vertexSet[j]->blind       == NULL)      continue;
    vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
    if (vertexset->analysis              == NULL)      continue;
    if (vertexset->discr                 == NULL)      continue;
    discr = vertexset->discr;
    if (discr->nPoints                   == 0)         continue;
    if (discr->dim                       == 3)         continue;
    for (k = 0; k < vertexset->nDataSets; k++) {
      if (vertexset->dataSets[k]              == NULL)        continue;
      if (vertexset->dataSets[k]->magicnumber != CAPSMAGIC)   continue;
      if (vertexset->dataSets[k]->type        != DATASET)     continue;
      if (vertexset->dataSets[k]->blind       == NULL)        continue;
      dataset = (capsDataSet *) vertexset->dataSets[k]->blind;
      if ((dataset->ftype != GeomSens) &&
          (dataset->ftype != TessSens))                       continue;
      if (vertexset->dataSets[k]->last.sNum   >=
          problem->geometry.sNum)                             continue;
      if (names == NULL) {
        names = (char **) EG_alloc(sizeof(char *));
        if (names == NULL) continue;
        names[0] = vertexset->dataSets[k]->name;
        n++;
      } else {
        for (m = 0; m < n; m++)
          if (strcmp(names[m], vertexset->dataSets[k]->name) == 0) break;
        if (m == n) {
          ntmp = (char **) EG_reall(names, (n+1)*sizeof(char *));
          if (ntmp == NULL) continue;
          names    = ntmp;
          names[n] = vertexset->dataSets[k]->name;
          n++;
        }
      }
    }
  }
  if (names != NULL) {

    /* go through the list of sensitivity names */
    for (m = 0; m < n; m++) {
      irow = icol = 1;
      open = 0;
      str  = EG_strdup(names[m]);
      if (str == NULL) continue;
      len = strlen(str);
      for (j = 1; j < len; j++)
        if (str[j] == '[') {
          open = j;
          break;
        }
      if (open != 0) {
        str[open] = 0;
        for (i = open+1; i < len; i++)
          if (str[i] == ',') str[i] = ' ';
        sscanf(&str[open+1], "%d%d\n", &irow, &icol);
      }
      open = -1;
      for (i = 0; i < problem->nGeomIn; i++)
        if (strcmp(str,problem->geomIn[i]->name) == 0) {
          value = (capsValue *) problem->geomIn[i]->blind;
          if (value == NULL) continue;
          open  = value->pIndex;
          break;
        }
      EG_free(str);
      if (open == -1) continue;

      /* clear all then set */
      status = ocsmSetDtime(problem->modl, 0);
      if (status != SUCCESS) return status;
      status = ocsmSetVelD(problem->modl, 0,    0,    0,    0.0);
      if (status != SUCCESS) return status;
      status = ocsmSetVelD(problem->modl, open, irow, icol, 1.0);
      if (status != SUCCESS) return status;
      buildTo = 0;
      nbody   = 0;
      outLevel = ocsmSetOutLevel(0);
      printf(" CAPS Info: Building sensitivity information for: %s[%d,%d]\n",
             problem->geomIn[i]->name, irow, icol);
      status  = ocsmBuild(problem->modl, buildTo, &builtTo, &nbody, NULL);
      ocsmSetOutLevel(outLevel);
      printf(" CAPS Info: parameter %d %d %d sensitivity status = %d\n",
             open, irow, icol, status);
      fflush(stdout);
      if (status != SUCCESS) continue;
      caps_geomOutSensit(problem, open, irow, icol);

      for (j = 0; j < bound->nVertexSet; j++) {
        if (bound->vertexSet[j]              == NULL)      continue;
        if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
        if (bound->vertexSet[j]->type        != VERTEXSET) continue;
        if (bound->vertexSet[j]->blind       == NULL)      continue;
        vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
        if (vertexset->analysis              == NULL)      continue;
        if (vertexset->discr                 == NULL)      continue;
        discr = vertexset->discr;
        if (discr->nPoints                   == 0)         continue;
        if (discr->dim                       == 3)         continue;
        for (k = 0; k < vertexset->nDataSets; k++) {
          if (vertexset->dataSets[k]              == NULL)        continue;
          if (vertexset->dataSets[k]->magicnumber != CAPSMAGIC)   continue;
          if (vertexset->dataSets[k]->type        != DATASET)     continue;
          if (vertexset->dataSets[k]->blind       == NULL)        continue;
          dobject = vertexset->dataSets[k];
          dataset = (capsDataSet *) dobject->blind;
          if ((dataset->ftype != GeomSens) &&
              (dataset->ftype != TessSens))                       continue;
          if (vertexset->dataSets[k]->last.sNum   >=
              problem->geometry.sNum)                             continue;
          if (strcmp(names[m],dobject->name)      != 0)           continue;
          dataset->data = (double *)
            EG_alloc(3*discr->nPoints*sizeof(double));
          if (dataset->data                       == NULL)        continue;
          caps_fillSensit(problem, discr, dataset);
          dataset->npts = discr->nPoints;
          caps_freeOwner(&dobject->last);
          vertexset->dataSets[k]->last.sNum = problem->sNum;
          status = caps_addHistory(dobject, problem);
          if (status != CAPS_SUCCESS)
            printf(" CAPS Warning: caps_addHistory = %d (caps_refillBound)\n",
                   status);

          status = caps_writeDataSet(dobject);
          if (status != CAPS_SUCCESS)
            printf(" CAPS Warning: caps_writeDataSet = %d (caps_refillBound)\n",
                   status);
        }

        MODL = (modl_T *) problem->modl;
        if (MODL->dtime != 0)
          printf(" CAPS Info: Sensitivity finite difference used for: %s[%d,%d]\n",
                 problem->geomIn[i]->name, irow, icol);
      }
    }
    EG_free(names);
  }

  return CAPS_SUCCESS;
}


static int
caps_toktokcmp(const char *str1, const char *str2)
{
  char copy1[1025], copy2[1025], *word1, *word2, *last1, *last2, sep[3] = " ;";

  /* do a token/token compare against 2 strings */
  strncpy(copy1, str1, 1024);
  copy1[1024] = 0;
  word1 = strtok_r(copy1, sep, &last1);
  while (word1 != NULL) {
    strncpy(copy2, str2, 1024);
    copy2[1024] = 0;
    word2 = strtok_r(copy2, sep, &last2);
    while (word2 != NULL) {
      if (strcmp(word2, word1) == 0) return 0;
      word2 = strtok_r(NULL, sep, &last2);
    }
    word1 = strtok_r(NULL, sep, &last1);
  }

  return 1;
}


int
caps_filter(capsProblem *problem, capsAnalysis *analysis)
{
  int          i, nBody, alen, atype, status;
  ego          *bodies;
  const int    *aints;
  const double *areals;
  const char   *astr;

  bodies = (ego *) EG_alloc(problem->nBodies*sizeof(ego));
  if (bodies == NULL) return EGADS_MALLOC;
  for (i = 0; i < problem->nBodies; i++) bodies[i] = NULL;

  /* do any bodies have our "capsAIM" attribute? */
  for (nBody = i = 0; i < problem->nBodies; i++) {
    status = EG_attributeRet(problem->bodies[i], "capsAIM", &atype,
                             &alen, &aints, &areals, &astr);
    if (status != EGADS_SUCCESS) continue;
    if (atype  != ATTRSTRING)    continue;
    if (caps_toktokcmp(analysis->loadName, astr) != 0) continue;
    if (analysis->intents != NULL) {
      status = EG_attributeRet(problem->bodies[i], "capsIntent", &atype,
                               &alen, &aints, &areals, &astr);
      if (status != EGADS_SUCCESS) continue;
      if (atype  != ATTRSTRING)    continue;
      if (caps_toktokcmp(analysis->intents, astr) != 0) continue;
    }
    bodies[nBody] = problem->bodies[i];
    nBody++;
  }

  if (nBody == 0)
    if (analysis->intents == NULL) {
      printf(" caps_filter Warning: No bodies with capsAIM = %s!\n",
             analysis->loadName);
      /* print all "capsAIM" attribute */
      for (i = 0; i < problem->nBodies; i++) {
        status = EG_attributeRet(problem->bodies[i], "capsAIM", &atype,
                                 &alen, &aints, &areals, &astr);
        if (status != EGADS_SUCCESS) continue;
        if (atype  == ATTRSTRING)
          printf("   Body %d capsAIM = %s\n",i+1, astr);
        else
          printf("   Body %d capsAIM is not a string\n",i+1);
      }
    } else {
      printf(" caps_filter Warning: No bodies with capsAIM = %s and capsIntent = %s!\n",
             analysis->loadName, analysis->intents);
      /* print all "capsAIM" and "capsIntent" attribute */
      for (i = 0; i < problem->nBodies; i++) {
        status = EG_attributeRet(problem->bodies[i], "capsAIM", &atype,
                                 &alen, &aints, &areals, &astr);
        if (status == EGADS_SUCCESS) {
          if (atype  == ATTRSTRING)
            printf("   Body %d capsAIM = %s\n",i+1, astr);
          else
            printf("   Body %d capsAIM is not a string\n",i+1);
        }
        status = EG_attributeRet(problem->bodies[i], "capsIntent", &atype,
                                 &alen, &aints, &areals, &astr);
        if (status == EGADS_SUCCESS) {
          if (atype  == ATTRSTRING)
            printf("   Body %d capsIntent = %s\n",i+1, astr);
          else
            printf("   Body %d capsIntent is not a string\n",i+1);
        }
      }
    }
  analysis->bodies = bodies;
  analysis->nBody  = nBody;
  return CAPS_SUCCESS;
}


int
caps_getBodies(capsObject *aobject, int *nBody, ego **bodies,
               int *nErr, capsErrs **errors)
{
  int          i, stat, ret, oclass, mtype, *senses;
  CAPSLONG     sNum;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsJrnl     args[3];
  ego          model, ref, *egos = NULL;

  *nErr   = 0;
  *errors = NULL;
  *nBody  = 0;
  *bodies = NULL;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  stat = caps_findProblem(aobject, CAPS_GETBODIES, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;

  args[0].type   = jEgos;
  args[1].type   = jInteger;
  args[2].type   = jErr;
  if (problem->dbFlag == 0) {
    stat         = caps_jrnlRead(CAPS_GETBODIES, problem, aobject, 3, args,
                                 &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *nErr   = args[1].members.integer;
      *errors = args[2].members.errs;
      if (args[0].members.model != NULL) {
        i = EG_getTopology(args[0].members.model, &ref, &oclass, &mtype, NULL,
                           nBody, bodies, &senses);
        if (i != EGADS_SUCCESS)
          printf(" CAPS Warning: EG_getTopology = %d (caps_getBodies)\n", i);
      }
      return ret;
    }
  }

  analysis = (capsAnalysis *) aobject->blind;
  sNum     = problem->sNum;
  if ((analysis->bodies == NULL) && (problem->dbFlag == 0)) {
    /* make sure geometry is up-to-date */
    stat = caps_build(pobject, nErr, errors);
    if ((stat != CAPS_SUCCESS) && (stat != CAPS_CLEAN)) goto cleanup;
    stat = caps_filter(problem, analysis);
    if (stat != CAPS_SUCCESS) goto cleanup;
  }
  *nBody  = analysis->nBody;
  *bodies = analysis->bodies;
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  if ((*nBody == 0) || (*bodies == NULL)) {
    *nBody  = 0;
    *bodies = NULL;
    args[0].members.model = NULL;
  } else if (*nBody == 1) {
    args[0].members.model = *bodies[0];
  } else {
    egos = (ego *) EG_alloc(*nBody*sizeof(ego));
    if (egos == NULL) {
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < *nBody; i++) {
      stat = EG_copyObject((*bodies)[i], NULL, &egos[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }
    stat = EG_makeTopology(problem->context, NULL, MODEL, 0, NULL, *nBody,
                           egos, NULL, &model);
    if (stat != EGADS_SUCCESS) {
      printf(" CAPS Error: EG_makeTopology %d = %d (caps_getBodies)!\n",
             *nBody, stat);
      goto cleanup;
    }
    args[0].members.model = model;
  }
  stat = CAPS_SUCCESS;

cleanup:
  args[1].members.integer = *nErr;
  args[2].members.errs    = *errors;
  caps_jrnlWrite(CAPS_GETBODIES, problem, aobject, stat, 3, args, sNum,
                 problem->sNum);

  if (egos != NULL) EG_free(egos);

  return stat;
}


int
caps_getTessels(capsObject *aobject, int *nTessel, ego **tessels,
                int *nErr, capsErrs **errors)
{
  int          i, j, n, npts, stat, ret, oclass, mtype, *senses;
  CAPSLONG     sNum;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsJrnl     args[3];
  ego          model, ref, *bodies, *egos = NULL;

  *nErr    = 0;
  *errors  = NULL;
  *nTessel = 0;
  *tessels = NULL;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  stat = caps_findProblem(aobject, CAPS_GETTESSELS, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;

  args[0].type   = jEgos;
  args[1].type   = jInteger;
  args[2].type   = jErr;
  if (problem->dbFlag == 0) {
    stat         = caps_jrnlRead(CAPS_GETTESSELS, problem, aobject, 3, args,
                                 &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *nErr   = args[1].members.integer;
      *errors = args[2].members.errs;
      if (args[0].members.model != NULL) {
        i = EG_getTopology(args[0].members.model, &ref, &oclass, &mtype, NULL,
                           &n, &bodies, &senses);
        if (i != EGADS_SUCCESS)
          printf(" CAPS Warning: EG_getTopology = %d (caps_getTessels)\n", i);
        if (mtype > n) {
          *nTessel = mtype - n;
          *tessels = &bodies[n];
        }
      }
      return ret;
    }
  }

  analysis = (capsAnalysis *) aobject->blind;
  sNum     = problem->sNum;
  *nTessel = analysis->nTess;
  *tessels = analysis->tess;
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  if (analysis->nTess == 0) {
    args[0].members.model = NULL;
  } else {
    /* find all of the referenced bodies */
    egos = (ego *) EG_alloc(3*analysis->nTess*sizeof(ego));
    if (egos == NULL) {
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    bodies = &egos[2*analysis->nTess];
    for (n = i = 0; i < analysis->nTess; i++) {
      stat = EG_statusTessBody(analysis->tess[i], &ref, &j, &npts);
      if (stat != EGADS_SUCCESS) {
        printf(" CAPS Error: EG_statusTessBody %d = %d (caps_getTessels)!\n",
               i, stat);
        goto cleanup;
      }
      for (j = 0; j < n; j++)
        if (bodies[j] == ref) break;
      if (j == n) {
        bodies[n] = ref;
        n++;
      }
    }

    /* copy the body objects */
    for (i = 0; i < n; i++) {
      stat = EG_copyObject(bodies[i], NULL, &egos[i]);
      if (stat != EGADS_SUCCESS) {
        for (j = 0; j < i; j++) EG_deleteObject(egos[j]);
        goto cleanup;
      }
    }
    /* copy the tessellation objects*/
    for (i = 0; i < analysis->nTess; i++) {
      stat = EG_statusTessBody(analysis->tess[i], &ref, &j, &npts);
      if (stat != EGADS_SUCCESS) {
        printf(" CAPS Error: EG_statusTessBody %d = %d (caps_getTessels)!\n",
               i, stat);
        for (j = 0; j < n+i; j++) EG_deleteObject(egos[j]);
        goto cleanup;
      }
      for (j = 0; j < n; j++)
        if (bodies[j] == ref) break;
      if (j == n) {
        stat = EGADS_NOTFOUND;
        printf(" CAPS Error: Cannot find Body %d (caps_getTessels)!\n", i);
        for (j = 0; j < n+i; j++) EG_deleteObject(egos[j]);
        goto cleanup;
      }
      stat = EG_copyObject(analysis->tess[i], egos[j], &egos[n+i]);
      if (stat != EGADS_SUCCESS) {
        for (j = 0; j < n+i; j++) EG_deleteObject(egos[j]);
        goto cleanup;
      }
    }
    /* make the model */
    stat = EG_makeTopology(problem->context, NULL, MODEL, analysis->nTess+n,
                           NULL, n, egos, NULL, &model);
    if (stat != EGADS_SUCCESS) {
      printf(" CAPS Error: EG_makeTopology %d = %d (caps_getTessels)!\n",
             analysis->nTess+n, stat);
      for (j = 0; j < analysis->nTess+n; j++) EG_deleteObject(egos[j]);
      goto cleanup;
    }
    args[0].members.model = model;
  }
  stat = CAPS_SUCCESS;

cleanup:
  args[1].members.integer = *nErr;
  args[2].members.errs    = *errors;
  caps_jrnlWrite(CAPS_GETTESSELS, problem, aobject, stat, 3, args, sNum,
                 problem->sNum);

  if (egos != NULL) EG_free(egos);

  return stat;
}


static int
caps_preAnalysiX(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int           i, j, k, npts, rank, stat, status, gstatus;
  double        *data;
  char          *units;
  int           majorX, minorX, nFieldX, execX, dirtyX;
  int           *ranksX, *fInOutX;
  char          *intentsX, *apathX, *unitSysX, **fnamesX;
  CAPSLONG      sn;
  capsValue     *valIn, *value;
  capsProblem   *problem;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;
  capsObject    *pobject, *source, *object, *last;

  analysis = NULL;
  valIn    = NULL;
  if (aobject->type == PROBLEM) {
    pobject  = aobject;
    problem  = (capsProblem *)  pobject->blind;
  } else {
    analysis = (capsAnalysis *) aobject->blind;
    pobject  = (capsObject *)   aobject->parent;
    problem  = (capsProblem *)  pobject->blind;
    valIn    = (capsValue *)    analysis->analysisIn[0]->blind;
  }

  /* make sure geometry is up-to-date */
  gstatus = caps_build(pobject, nErr, errors);
  if (analysis == NULL) {
    if (gstatus == CAPS_CLEAN) return CAPS_CLEAN;
    if (problem->nBodies == 0)
      printf(" caps_preAnalysis Warning: No bodies generated!\n");
    return CAPS_SUCCESS;
  }
  if ((gstatus != CAPS_SUCCESS) && (gstatus != CAPS_CLEAN)) return gstatus;

  /* are we "analysis" clean? */
  status = 0;
  if (analysis != NULL) {
    if (analysis->pre.sNum == 0) {
      status = 1;
    } else {
      /* check for dirty analysis inputs */
      for (i = 0; i < analysis->nAnalysisIn; i++) {
        source = object = analysis->analysisIn[i];
        do {
          if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
          if (source->type        != VALUE)     return CAPS_BADTYPE;
          if (source->blind       == NULL)      return CAPS_NULLBLIND;
          value = (capsValue *) source->blind;
          if (value->link == object)            return CAPS_CIRCULARLINK;
          last   = source;
          source = value->link;
        } while (value->link != NULL);
        if (last->last.sNum > analysis->pre.sNum) {
          status = 1;
          break;
        }
        /* if the value is linked, check the state of the linked parent */
        if ((last != object) && (last->parent != NULL) &&
            (last->parent->type == ANALYSIS)) {
          stat = caps_analysisInfX(last->parent, &apathX, &unitSysX, &majorX,
                                   &minorX, &intentsX, &nFieldX, &fnamesX,
                                   &ranksX, &fInOutX, &execX, &dirtyX);
          if (stat != CAPS_SUCCESS) return stat;
          if (dirtyX > 0) {
            status = 1;
            break;
          }
        }
      }
      if (status == 0) {
        stat = caps_snDataSets(aobject, 1, &sn);
        if (stat == CAPS_SUCCESS)
          if (sn > analysis->pre.sNum) status = 1;
      }
    }
    if ((status == 0) && (gstatus == CAPS_CLEAN) &&
        (problem->geometry.sNum < analysis->pre.sNum)) return CAPS_CLEAN;

    /* fill in any values that have links */
    stat = caps_transferLinks(analysis, nErr, errors);
    if (stat != CAPS_SUCCESS) return stat;
  }

  if ((problem->nBodies <= 0) || (problem->bodies == NULL)) {
    printf(" caps_preAnalysis Warning: No bodies to %s!\n", analysis->loadName);
    if (analysis->bodies != NULL) EG_free(analysis->bodies);
    analysis->bodies = NULL;
    analysis->nBody  = 0;
  } else if (analysis->bodies == NULL) {
    stat = caps_filter(problem, analysis);
    if (stat != CAPS_SUCCESS) return stat;
  }

  /* do we need to fill FieldIn datasets? */
  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]              == NULL)      return CAPS_NULLOBJ;
    if (problem->bounds[i]->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (problem->bounds[i]->type        != BOUND)     return CAPS_BADTYPE;
    if (problem->bounds[i]->blind       == NULL)      return CAPS_NULLBLIND;
    bound = (capsBound *) problem->bounds[i]->blind;
    for (j = 0; j < bound->nVertexSet; j++) {
      if (bound->vertexSet[j]              == NULL)      return CAPS_NULLOBJ;
      if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (bound->vertexSet[j]->type        != VERTEXSET) return CAPS_BADTYPE;
      if (bound->vertexSet[j]->blind       == NULL)      return CAPS_NULLBLIND;
      vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
      if (vertexset->analysis != aobject) continue;
      for (k = 0; k < vertexset->nDataSets; k++) {
        if (vertexset->dataSets[k]              == NULL)      return CAPS_NULLOBJ;
        if (vertexset->dataSets[k]->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
        if (vertexset->dataSets[k]->type        != DATASET)   return CAPS_BADTYPE;
        if (vertexset->dataSets[k]->blind       == NULL)      return CAPS_NULLBLIND;
        dataset = (capsDataSet *) vertexset->dataSets[k]->blind;
        if (dataset->ftype != FieldIn) continue;
        stat = caps_getDataX(vertexset->dataSets[k], &npts, &rank, &data,
                             &units, nErr, errors);
        if (stat != CAPS_SUCCESS) return stat;
      }
    }
  }

  /* running preAnalysis means reload is no longer viable */
  analysis->reload = 0;

  /* do it! */
  stat = aim_PreAnalysis(problem->aimFPTR, analysis->loadName,
                         analysis->instStore, &analysis->info, valIn);
  caps_getAIMerrs(analysis, nErr, errors);
  if (stat == CAPS_SUCCESS) {
    caps_freeValueObjects(1, analysis->nAnalysisDynO, analysis->analysisDynO);
    analysis->analysisDynO  = NULL;
    analysis->nAnalysisDynO = 0;
    caps_freeOwner(&analysis->pre);
    problem->sNum     += 1;
    analysis->pre.sNum = problem->sNum;
    caps_fillDateTime(analysis->pre.datetime);
    stat = caps_dumpAnalysis(problem, aobject);
    if (stat != CAPS_SUCCESS)
      printf(" CAPS Error: caps_dumpAnalysis = %d (caps_preAnalysis)\n",
             stat);
  }

  return stat;
}


static int
caps_preAnalysiZ(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int          stat, ret;
  CAPSLONG     sNum;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsObject   *pobject;
  capsJrnl     args[2];

  *nErr    = 0;
  *errors  = NULL;
  analysis = NULL;
  if (aobject                          == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber             != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type                    == PROBLEM) {
    if (aobject->blind                 == NULL)      return CAPS_NULLBLIND;
    pobject  = aobject;
    problem  = (capsProblem *)  pobject->blind;
  } else {
    if (aobject->type                  != ANALYSIS)  return CAPS_BADTYPE;
    if (aobject->blind                 == NULL)      return CAPS_NULLBLIND;
    analysis = (capsAnalysis *) aobject->blind;
    if (aobject->parent                == NULL)      return CAPS_NULLOBJ;
    pobject  = (capsObject *)   aobject->parent;
    if (pobject->blind                 == NULL)      return CAPS_NULLBLIND;
    problem  = (capsProblem *)  pobject->blind;
    if (analysis->analysisIn[0]->blind == NULL)      return CAPS_NULLVALUE;
  }

  args[0].type = jInteger;
  args[1].type = jErr;
  stat         = caps_jrnlRead(CAPS_PREANALYSIS, problem, aobject, 2, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *nErr    = args[0].members.integer;
    *errors  = args[1].members.errs;
    return ret;
  }

  ret   = caps_preAnalysiX(aobject, nErr, errors);
  *nErr = 0;
  if (*errors != NULL) *nErr = (*errors)->nError;

  args[0].members.integer = *nErr;
  args[1].members.errs    = *errors;
  caps_jrnlWrite(CAPS_PREANALYSIS, problem, aobject, ret, 2, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_preAnalysis(capsObject *aobject, int *nErr, capsErrs **errors)
{
  char         temp[PATH_MAX];
  capsAnalysis *analysis;
  capsProblem  *problem;
  capsObject   *pobject;

  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  pobject  = (capsObject *) aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  analysis = (capsAnalysis *) aobject->blind;
  if ((analysis->autoexec  == 1) &&
      (analysis->eFlag     == 1)) {
    snprintf(temp, PATH_MAX, "Cannot call preAnalysis with AIM that auto-executes!");
    caps_makeSimpleErr(aobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_EXECERR;
  }

  return caps_preAnalysiZ(aobject, nErr, errors);
}


static int
caps_boundDependent(capsProblem *problem, capsObject *aobject,
                    capsObject  *oobject)
{
  int           i, j, jj, k, kk;
  capsBound     *bound;
  capsVertexSet *vs, *vso;
  capsDataSet   *ds, *dso;

  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]                 == NULL)      continue;
    if (problem->bounds[i]->magicnumber    != CAPSMAGIC) continue;
    if (problem->bounds[i]->type           != BOUND)     continue;
    if (problem->bounds[i]->blind          == NULL)      continue;
    bound = (capsBound *) problem->bounds[i]->blind;
    for (j = 0; j < bound->nVertexSet; j++) {
      if (bound->vertexSet[j]              == NULL)      continue;
      if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
      if (bound->vertexSet[j]->type        != VERTEXSET) continue;
      if (bound->vertexSet[j]->blind       == NULL)      continue;
      vs = (capsVertexSet *) bound->vertexSet[j]->blind;
      if (vs->analysis                     != aobject)   continue;
      for (k = 0; k < vs->nDataSets; k++) {
        if (vs->dataSets[k]                == NULL)      continue;
        if (vs->dataSets[k]->magicnumber   != CAPSMAGIC) continue;
        if (vs->dataSets[k]->type          != DATASET)   continue;
        if (vs->dataSets[k]->blind         == NULL)      continue;
        ds = (capsDataSet *) vs->dataSets[k]->blind;
        if (ds->ftype != FieldIn)                        continue;
        for (jj = 0; jj < bound->nVertexSet; jj++) {
          if (j == jj) continue;
          if (bound->vertexSet[jj]              == NULL)       continue;
          if (bound->vertexSet[jj]->magicnumber != CAPSMAGIC)  continue;
          if (bound->vertexSet[jj]->type        != VERTEXSET)  continue;
          if (bound->vertexSet[jj]->blind       == NULL)       continue;
          vso = (capsVertexSet *) bound->vertexSet[jj]->blind;
          if (vso->analysis                     != oobject)    continue;
          for (kk = 0; kk < vso->nDataSets; kk++) {
            if (vso->dataSets[kk]                == NULL)      continue;
            if (vso->dataSets[kk]->magicnumber   != CAPSMAGIC) continue;
            if (vso->dataSets[kk]->type          != DATASET)   continue;
            if (vso->dataSets[kk]->blind         == NULL)      continue;
            dso = (capsDataSet *) vso->dataSets[kk]->blind;
            if (dso->ftype                       != FieldOut)  continue;
            if (strcmp(vs->dataSets[k]->name,vso->dataSets[kk]->name) == 0)
              return CAPS_SUCCESS;
          }
        }
      }
    }
  }

  return CAPS_NOTFOUND;
}


int
caps_postAnalysiX(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int          i, major, minor, nField, *ranks, *fInOut, status, exec, dirty;
  char         *intents, *apath, *unitSys, **fnames;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsValue    *valIn = NULL;
  capsObject   *pobject;

  analysis = (capsAnalysis *) aobject->blind;
  pobject  = (capsObject *)   aobject->parent;
  problem  = (capsProblem *)  pobject->blind;

  /* check to see if we need to do post */
  status = caps_analysisInfX(aobject, &apath, &unitSys, &major, &minor, &intents,
                             &nField, &fnames, &ranks, &fInOut, &exec, &dirty);
  if (status != CAPS_SUCCESS) return status;
  if (dirty == 0) return CAPS_CLEAN;
  if (dirty <  5) return CAPS_DIRTY;

  /* pre was called with journaling, so links need to be re-populated */
  if (analysis->reload == 2) {
    /* fill in any values that have links */
    status = caps_transferLinks(analysis, nErr, errors);
    if (status != CAPS_SUCCESS) return status;
    analysis->reload = 0;
  }

  if (analysis->nAnalysisIn > 0)
    valIn = (capsValue *) analysis->analysisIn[0]->blind;

  /* call post in the AIM */
  analysis->info.inPost = 1;
  status = aim_PostAnalysis(problem->aimFPTR, analysis->loadName,
                            analysis->instStore, &analysis->info, 0, valIn);
  analysis->info.inPost = 0;
  caps_getAIMerrs(analysis, nErr, errors);
  if (status != CAPS_SUCCESS) {
    caps_freeValueObjects(1, analysis->nAnalysisDynO, analysis->analysisDynO);
    analysis->analysisDynO  = NULL;
    analysis->nAnalysisDynO = 0;
    return status;
  }

  /* set the time/date stamp */
  caps_freeOwner(&aobject->last);
  problem->sNum     += 1;
  aobject->last.sNum = problem->sNum;
  status = caps_addHistory(aobject, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_postAnalysis)\n", status);
  for (i = 0; i < analysis->nAnalysisDynO; i++) {
    analysis->analysisDynO[i]->last.sNum = problem->sNum;
    status = caps_addHistory(analysis->analysisDynO[i], problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_addHistory = %d (caps_postAnalysis)\n",
             status);
  }
  status = caps_dumpAnalysis(problem, aobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Error: caps_dumpAnalysis = %d (caps_postAnalysis)\n", status);

  return status;
}


static int
caps_postAnalysiZ(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int          stat, ret;
  CAPSLONG     sNum;
  capsProblem  *problem;
  capsObject   *pobject;
  capsJrnl     args[2];

  *nErr   = 0;
  *errors = NULL;
  if (aobject              == NULL)            return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC)       return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)        return CAPS_BADTYPE;
  if (aobject->blind       == NULL)            return CAPS_NULLBLIND;
  if (aobject->parent      == NULL)            return CAPS_NULLOBJ;
  pobject  = (capsObject *)   aobject->parent;
  if (pobject->blind       == NULL)            return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  args[0].type = jInteger;
  args[1].type = jErr;
  stat         = caps_jrnlRead(CAPS_POSTANALYSIS, problem, aobject, 2, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *nErr    = args[0].members.integer;
    *errors  = args[1].members.errs;
    return ret;
  }

  sNum = problem->sNum;
  ret  = caps_postAnalysiX(aobject, nErr, errors);
  args[0].members.integer = *nErr;
  args[1].members.errs    = *errors;
  caps_jrnlWrite(CAPS_POSTANALYSIS, problem, aobject, ret, 2, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_postAnalysis(capsObject *aobject, int *nErr, capsErrs **errors)
{
  char         temp[PATH_MAX];
  capsAnalysis *analysis;
  capsProblem  *problem;
  capsObject   *pobject;

  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  pobject  = (capsObject *) aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  analysis = (capsAnalysis *) aobject->blind;
  if ((analysis->autoexec == 1) &&
      (analysis->eFlag == 1)) {
    snprintf(temp, PATH_MAX, "Cannot call postAnalysis with AIM that auto-executes!");
    caps_makeSimpleErr(aobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_EXECERR;
  }

  return caps_postAnalysiZ(aobject, nErr, errors);
}


int
caps_buildBound(capsObject *bobject, int *nErr, capsErrs **errors)
{
  int           i, j, ia, it, stat, found = 0, dirty = 0;
  char          temp[PATH_MAX];
  capsObject    *pobject;
  capsProblem   *problem;
  capsBound     *bound;
  capsAnalysis  *analysis, *anal;
  capsVertexSet *vs;
  ego           tess;

  *nErr   = 0;
  *errors = NULL;
  if (bobject              == NULL)      return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)     return CAPS_BADTYPE;
  if (bobject->blind       == NULL)      return CAPS_NULLBLIND;
  bound   = (capsBound *) bobject->blind;
  if (bound->state         == Open)      return CAPS_STATEERR;
  pobject = bobject->parent;
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  if (bound->state == Open) {
    snprintf(temp, PATH_MAX, "Bound '%s' is Open (caps_buildBound)!",
             bobject->name);
    caps_makeSimpleErr(bobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_DIRTY;
  }

  stat = caps_build(pobject, nErr, errors);
  if ((stat != CAPS_CLEAN) && (stat != CAPS_SUCCESS)) return stat;

  if (bobject->last.sNum < problem->geometry.sNum) dirty = 1;

  for (i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL)      continue;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
    if (bound->vertexSet[i]->type        != VERTEXSET) continue;
    if (bound->vertexSet[i]->blind       == NULL)      continue;
    vs       = (capsVertexSet *) bound->vertexSet[i]->blind;
    analysis = (capsAnalysis *)  vs->analysis->blind;
    if (analysis                         == NULL)      continue;
    if (analysis->bodies == NULL) {
      stat = caps_filter(problem, analysis);
      if (stat != CAPS_SUCCESS) {
        snprintf(temp, PATH_MAX, "caps_filter = %d (caps_buildBound)!", stat);
        caps_makeSimpleErr(bobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return stat;
      }
      dirty = 1;
    }
    if (vs->discr == NULL) {
      dirty = 1;
      continue;
    }

    /* if the discr tessellation object is not found (or is EMPTY)
       then it has been replaced and the VertexSet needs to be updated */
    for (j = 0; j < vs->discr->nBodys; j++) {
      tess = vs->discr->bodys[j].tess;
      if (tess         == NULL ) { dirty = 1; continue; }
      if (tess->oclass == EMPTY) { dirty = 1; continue; }
      found = 0;
      for (ia = 0; ia < problem->nAnalysis && found == 0; ia++) {
        anal = (capsAnalysis *) problem->analysis[ia]->blind;
        for (it = 0; it < anal->nTess; it++) {
          if ( tess == anal->tess[it] ) {
            found = 1;
            break;
          }
        }
      }
      if (found == 0) dirty = 1;
    }
  }

  if (dirty == 0) return CAPS_SUCCESS;

  stat = caps_refillBound(problem, bobject, nErr, errors);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeBound(bobject);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeBound = %d (caps_buildBound)\n", stat);

  problem->sNum += 1;
  stat = caps_writeSerialNum(problem);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeSerialNum = %d (caps_buildBound)\n",
           stat);

  return stat;
}


static void
caps_orderAnalyses(int nAobj, /*@null@*/ capsObject **aobjs)
{
  int        i, hit;
  capsObject *aobj;

  if ((nAobj <= 1) || (aobjs == NULL)) return;

  do {
    for (hit = i = 0; i < nAobj-1; i++) {
      if (aobjs[i]->last.sNum > aobjs[i+1]->last.sNum) {
        aobj       = aobjs[i];
        aobjs[i]   = aobjs[i+1];
        aobjs[i+1] = aobj;
        hit++;
      }
    }
  } while (hit != 0);

}


static int
caps_dirtyAnalysiX(capsObject *object, capsProblem *problem,
                   int *nAobj, capsObject ***aobjs)
{
  int           i, j, execute, major, minor, nField, *ranks, *fInOut;
  int           stat, dirty;
  char          *intents, *apath, *unitSys, **fnames;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsAnalysis  *analysis;
  capsObject    *aobject, *source, *start, *last;
  capsValue     *value;

  *nAobj = 0;
  *aobjs = NULL;
  if  (object              == NULL)         return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if ((object->type        != PROBLEM) && (object->type != ANALYSIS) &&
      (object->type        != BOUND))       return CAPS_BADTYPE;
  if  (object->blind       == NULL)         return CAPS_NULLBLIND;

  if (object->type == PROBLEM) {
    if (problem->nAnalysis   == 0)         return CAPS_SUCCESS;

    for (i = 0; i < problem->nAnalysis; i++) {
      stat = caps_analysisInfX(problem->analysis[i], &apath, &unitSys, &major,
                               &minor, &intents, &nField, &fnames, &ranks,
                               &fInOut, &execute, &dirty);
      if (stat != CAPS_SUCCESS) {
        if (*aobjs != NULL) {
          EG_free(*aobjs);
          *nAobj = 0;
          *aobjs = NULL;
        }
        return stat;
      }
      if (dirty == 0) continue;
      if (*aobjs == NULL) {
        *aobjs = (capsObject **)
                 EG_alloc(problem->nAnalysis*sizeof(capsObject *));
        if (*aobjs == NULL) return EGADS_MALLOC;
      }
      (*aobjs)[*nAobj] = problem->analysis[i];
      *nAobj += 1;
    }

  } else if (object->type == BOUND) {

    /* for Bound Objects -- find dependent Analysis Objects */
    bound  = (capsBound *) object->blind;
    for (i = 0; i < bound->nVertexSet; i++) {
      if (bound->vertexSet[i]              == NULL)      continue;
      if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) continue;
      if (bound->vertexSet[i]->type        != VERTEXSET) continue;
      if (bound->vertexSet[i]->blind       == NULL)      continue;
      vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
      aobject   = vertexset->analysis;
      if (aobject                          == NULL)      continue;
      stat = caps_analysisInfX(aobject, &apath, &unitSys, &major, &minor,
                               &intents, &nField, &fnames, &ranks, &fInOut,
                               &execute, &dirty);
      if (stat != CAPS_SUCCESS) {
        if (*aobjs != NULL) {
          EG_free(*aobjs);
          *nAobj = 0;
          *aobjs = NULL;
        }
        return stat;
      }
      if (dirty == 0) continue;
      if (*aobjs == NULL) {
        *aobjs = (capsObject **)
                 EG_alloc(problem->nAnalysis*sizeof(capsObject *));
        if (*aobjs == NULL) return EGADS_MALLOC;
      }
      for (j = 0; j < *nAobj; j++)
        if (aobject == (*aobjs)[j]) break;
      if (j != *nAobj) continue;
      (*aobjs)[*nAobj] = aobject;
      *nAobj += 1;
    }

  } else {

    /* for Analysis Objects -- find dependent Analysis Objects in Bounds */
    for (i = 0; i < problem->nAnalysis; i++) {
      if (problem->analysis[i] == object) continue;
      stat = caps_boundDependent(problem, object, problem->analysis[i]);
      if (stat != CAPS_SUCCESS) continue;
      stat = caps_analysisInfX(problem->analysis[i], &apath, &unitSys, &major,
                               &minor, &intents, &nField, &fnames, &ranks,
                               &fInOut, &execute, &dirty);
      if (stat != CAPS_SUCCESS) {
        if (*aobjs != NULL) {
          EG_free(*aobjs);
          *nAobj = 0;
          *aobjs = NULL;
        }
        return stat;
      }
      if (dirty == 0) continue;
      if (*aobjs == NULL) {
        *aobjs = (capsObject **)
                 EG_alloc(problem->nAnalysis*sizeof(capsObject *));
        if (*aobjs == NULL) return EGADS_MALLOC;
      }
      (*aobjs)[*nAobj] = problem->analysis[i];
      *nAobj += 1;
    }
    /* now find dependent Analysis Objects in Links */
    analysis = (capsAnalysis *) object->blind;
    for (i = 0; i < analysis->nAnalysisIn; i++) {
      source = start = analysis->analysisIn[i];
      do {
        if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
        if (source->type        != VALUE)     return CAPS_BADTYPE;
        if (source->blind       == NULL)      return CAPS_NULLBLIND;
        value = (capsValue *) source->blind;
        if (value->link         == start)     return CAPS_CIRCULARLINK;
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      /* reject GeomOuts and save the Analysis Objects for AnalyOuts */
      source = last->parent;
      if ((source->type != ANALYSIS) || (source == object)) continue;
      stat = caps_analysisInfX(source, &apath, &unitSys, &major, &minor,
                               &intents, &nField, &fnames, &ranks, &fInOut,
                               &execute, &dirty);
      if (stat != CAPS_SUCCESS) {
        if (*aobjs != NULL) {
          EG_free(*aobjs);
          *nAobj = 0;
          *aobjs = NULL;
        }
        return stat;
      }
      if (dirty == 0) continue;
      if (*aobjs != NULL) {
        for (j = 0; j < *nAobj; j++)
          if ((*aobjs)[j] == source) break;
      } else {
        j = 0;
      }
      if (j == *nAobj) {
        if (*aobjs == NULL) {
          *aobjs = (capsObject **)
                   EG_alloc(problem->nAnalysis*sizeof(capsObject *));
          if (*aobjs == NULL) return EGADS_MALLOC;
        }
        (*aobjs)[*nAobj] = source;
        *nAobj += 1;
      }
    }
  }

  caps_orderAnalyses(*nAobj, *aobjs);

  return CAPS_SUCCESS;
}


int
caps_dirtyAnalysis(capsObject *object, int *nAobj, capsObject ***aobjs)
{
  int         i, ret, stat;
  CAPSLONG    sNum;
  capsProblem *problem;
  capsObject  *pobject, **objs;
  capsJrnl    args[1];

  *nAobj = 0;
  *aobjs = NULL;
  if  (object              == NULL)         return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if ((object->type        != PROBLEM) && (object->type != ANALYSIS) &&
      (object->type        != BOUND))       return CAPS_BADTYPE;
  if  (object->blind       == NULL)         return CAPS_NULLBLIND;
  stat  = caps_findProblem(object, CAPS_DIRTYANALYSIS, &pobject);
  if (stat                 != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jObjs;
  if (problem->dbFlag == 0) {
    stat       = caps_jrnlRead(CAPS_DIRTYANALYSIS, problem, object, 1, args,
                               &sNum, &ret);
    if (stat == CAPS_JOURNALERR) return stat;
    if (stat == CAPS_JOURNAL) {
      *nAobj = args[0].num;
      if (*nAobj != 0) {
        objs = (capsObject **) EG_alloc(*nAobj*sizeof(capsObject *));
        if (objs == NULL) return EGADS_MALLOC;
        for (i = 0; i < *nAobj; i++) objs[i] = args[0].members.objs[i];
        *aobjs = objs;
      }
      return ret;
    }
  }

  sNum = problem->sNum;
  ret  = caps_dirtyAnalysiX(object, problem, nAobj, aobjs);
  if (problem->dbFlag == 1) return ret;

  args[0].num          = *nAobj;
  args[0].members.objs = *aobjs;
  caps_jrnlWrite(CAPS_DIRTYANALYSIS, problem, object, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


static int
caps_runAnalysis(capsObject *aobject, int *stat, int *nErr, capsErrs **errors)
{
  int          status;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;

  if (aobject              == NULL)       return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)   return CAPS_BADTYPE;
  if (aobject->blind       == NULL)       return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  status   = caps_findProblem(aobject, CAPS_RUNANALYSIS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)          return CAPS_READONLYERR;

  /* ignore if restarting
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) return CAPS_SUCCESS; */

  /* call the AIM */
  status = aim_Execute(problem->aimFPTR, analysis->loadName,
                       analysis->instStore, &analysis->info, stat);
  caps_getAIMerrs(analysis, nErr, errors);
  if (status != CAPS_SUCCESS) return status;

  return CAPS_SUCCESS;
}


#ifdef ASYNCEXEC
int
caps_checkAnalysis(capsObject *aobject, int *phase, int *nErr, capsErrs **errors)
{
  int          status;
  capsObject   *pobject;
  capsProblem  *problem;
  capsAnalysis *analysis;

  if (nErr                 == NULL)       return CAPS_NULLVALUE;
  if (errors               == NULL)       return CAPS_NULLVALUE;
  *phase  = 0;
  *nErr   = 0;
  *errors = NULL;
  if (aobject              == NULL)       return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)   return CAPS_BADTYPE;
  if (aobject->blind       == NULL)       return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  status   = caps_findProblem(aobject, CAPS_CHECKANALYSIS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  /* call the AIM */
  status = aim_Check(problem->aimFPTR, analysis->loadName,
                     analysis->instStore, &analysis->info, phase);
  caps_getAIMerrs(analysis, nErr, errors);
  if (status != CAPS_SUCCESS) return status;

  return CAPS_SUCCESS;
}
#endif


int
caps_execute(capsObject *object, int *state, int *nErr, capsErrs **errors)
{
  int          stat;
  capsErrs     *errs = NULL;
  capsAnalysis *analysis;
  capsObject   *pobject;
  capsProblem  *problem;

  if (nErr                 == NULL)      return CAPS_NULLVALUE;
  if (errors               == NULL)      return CAPS_NULLVALUE;
  *state  = 0;
  *nErr   = 0;
  *errors = NULL;
  pobject = object;
  if  (object              == NULL)      return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if ((object->type        != PROBLEM) &&
      (object->type        != ANALYSIS)) return CAPS_BADTYPE;
  if  (object->type        == ANALYSIS) {
    if (object->blind      == NULL)      return CAPS_NULLBLIND;
    analysis = (capsAnalysis *) object->blind;
    if (analysis->eFlag    != 1)         return CAPS_EXECERR;
    pobject = (capsObject *) object->parent;
  }
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;
  if (problem->dbFlag      == 1)         return CAPS_READONLYERR;

  /* perform pre, exec, post -- use individual journaling */
  stat = caps_preAnalysiZ(object, nErr, errors);
  if (stat != CAPS_SUCCESS)    return stat;
  if (object->type == PROBLEM) return CAPS_SUCCESS;
  if (*nErr != 0) {
    errs    = *errors;
    *nErr   = 0;
    *errors = NULL;
  }

  stat = caps_runAnalysis(object, state, nErr, errors);
  caps_concatErrs(errs, errors);
  if (stat != CAPS_SUCCESS) {
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    return stat;
  }
  if (*nErr != 0) {
    errs    = *errors;
    *nErr   = 0;
    *errors = NULL;
  }

  stat = caps_postAnalysiZ(object, nErr, errors);
  caps_concatErrs(errs, errors);
  *nErr = 0;
  if (*errors != NULL) *nErr = (*errors)->nError;
  return stat;
}


int
caps_execX(capsObject *aobject, int *nErr, capsErrs **errors)
{
  int          stat, state;
  capsErrs     *errs = NULL;
  capsObject   *pobject;
  capsAnalysis *analysis;
  capsProblem  *problem;

  if (nErr                 == NULL)      return CAPS_NULLVALUE;
  if (errors               == NULL)      return CAPS_NULLVALUE;
  if (*nErr                != 0)         errs = *errors;

  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  if (analysis->autoexec   != 1)         return CAPS_EXECERR;
  pobject  = (capsObject *)   aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  /* perform pre, exec, post -- no journaling! */
  stat = caps_preAnalysiX(aobject, nErr, errors);
  caps_concatErrs(errs, errors);
  if (stat != CAPS_SUCCESS) {
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    return stat;
  }
  if (*nErr != 0) {
    errs    = *errors;
    *nErr   = 0;
    *errors = NULL;
  }

  /* call the AIM to execute */
  stat = aim_Execute(problem->aimFPTR, analysis->loadName,
                     analysis->instStore, &analysis->info, &state);
  caps_getAIMerrs(analysis, nErr, errors);
  caps_concatErrs(errs, errors);
  if (stat != CAPS_SUCCESS) {
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    return stat;
  }
  if (*nErr != 0) {
    errs    = *errors;
    *nErr   = 0;
    *errors = NULL;
  }

  stat = caps_postAnalysiX(aobject, nErr, errors);
  caps_concatErrs(errs, errors);
  if (stat != CAPS_SUCCESS) {
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    return stat;
  }

  return CAPS_SUCCESS;
}
