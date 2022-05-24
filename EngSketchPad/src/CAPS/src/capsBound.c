/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Bound, VertexSet & DataSet Object Functions
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
#define PATH_MAX _MAX_PATH
#else
#include <limits.h>
#endif

#include "capsBase.h"
#include "capsAIM.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"



typedef struct {
  int    bIndex;              /* the body index in quilt */
  int    eIndex;              /* the element index in quilt */
  double st[3];               /* position in element ref coords */
} capsTarget;

typedef struct {
  capsTarget source;          /* the positions in the source */
  capsTarget target;          /* the positions in the target */
} capsMatch;

typedef struct {
  char       *name;
  aimContext *aimFPTR;
  int        npts;
  int        irank;
  int        nrank;
  int        sindx;           /* source DLL index */
  int        tindx;           /* target DLL index */
  double     afact;           /* area penalty function weight */
  double     area_src;        /* area associated with source (output) */
  double     area_tgt;        /* area associated with target (output) */
  capsDiscr  *src;
  double     *prms_src;
  double     *data_src;
  capsDiscr  *tgt;
  double     *prms_tgt;
  /*@exposed@*/
  double     *data_tgt;
  int        nmat;            /* number of MatchPoints */
  capsMatch  *mat;            /* array  of MatchPoints */
} capsConFit;


extern /*@null@*/ /*@only@*/ char *EG_strdup(/*@null@*/ const char *str);

extern void caps_initDiscr(capsDiscr *discr);
extern int  caps_mkDir(const char *path);
extern int  caps_rename(const char *src, const char *dst);
extern int  caps_isNameOK(const char *name);
extern int  caps_dumpBound(capsObject *pobject, capsObject *bobject);
extern int  caps_writeVertexSet(capsObject *vobject);
extern int  caps_writeDataSet(capsObject *dobject);
extern void caps_jrnlWrite(int funID, capsProblem *problem, capsObject *obj,
                           int status, int nargs, capsJrnl *args, CAPSLONG sNm0,
                           CAPSLONG sNum);
extern int  caps_jrnlEnd(capsProblem *problem);
extern int  caps_jrnlRead(int funID, capsProblem *problem, capsObject *obj,
                          int nargs, capsJrnl *args, CAPSLONG *sNum, int *stat);
extern int  caps_buildBound(capsObject *bobject, int *nErr, capsErrs **errors);
extern int  caps_unitParse(/*@null@*/ const char *unit);
extern int  caps_getDataX(capsObject *dobj, int *npts, int *rank, double **data,
                          char **units, int *nErr, capsErrs **errors);
extern void caps_concatErrs(/*@null@*/ capsErrs *errs, capsErrs **errors);
extern void caps_getAIMerrs(capsAnalysis *analy, int *nErr, capsErrs **errors);
extern int  caps_analysisInfX(const capsObj aobject, char **apath, char **unSys,
                              int *major, int *minor, char **intents,
                              int *nField, char ***fnames, int **ranks,
                              int **fInOut, int *execution, int *status);
extern int  caps_postAnalysiX(capsObject *aobject, int *nErr, capsErrs **errors,
                              int flag);
extern int  caps_execX(capsObject *aobject, int *nErr, capsErrs **errors);
extern int  caps_circularAutoExecs(capsObject *asrc,
                                   /*@null@*/ capsObject *aobject);

extern int
  caps_conjGrad(int (*func)(int, double[], void *, double *,
                            /*@null@*/ double[]), /* (in)  objective function */
                void   *blind,                    /* (in)  blind pointer to struct */
                int    n,                         /* (in)  number of variables */
                double x[],                       /* (in)  initial set of variables */
                                                  /* (out) optimized variables */
                double ftol,                      /* (in)  convergence tolerance on
                                                           objective function */
                /*@null@*/
                FILE   *fp,                       /* (in)  FILE to write history */
                double *fopt);                    /* (out) objective function */


int
caps_integrateData(const capsObject *object, enum capstMethod method,
                   int *rank, double **data, char **units)
{
  int           i, j, bIndex, status;
  double        *result;
  capsProblem   *problem;
  capsDataSet   *dataset;
  capsObject    *vso;
  capsVertexSet *vertexset;
  capsAnalysis  *analysis;
  capsDiscr     *discr;

  *rank  = 0;
  *data  = NULL;
  *units = NULL;          /* what do I fill this in with? */

  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->type        != DATASET)   return CAPS_BADTYPE;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  dataset = (capsDataSet *) object->blind;
  vso     = object->parent;
  if (vso                 == NULL)      return CAPS_NULLOBJ;
  if (vso->magicnumber    != CAPSMAGIC) return CAPS_BADOBJECT;
  if (vso->type           != VERTEXSET) return CAPS_BADTYPE;
  if (vso->blind          == NULL)      return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vso->blind;
  if (vertexset->analysis == NULL)      return CAPS_NULLOBJ;
  analysis  = (capsAnalysis *) vertexset->analysis;
  discr     =                  vertexset->discr;
  if (discr               == NULL)      return CAPS_NULLVALUE;
  problem   = analysis->info.problem;
  if (problem             == NULL)      return CAPS_NULLOBJ;


  *rank = dataset->rank;
  i     = *rank*2;
  if (method == Average) i = *rank*3;
  result = (double *) EG_alloc(i*sizeof(double));
  if (result == NULL) return EGADS_MALLOC;

  for (i = 0; i < *rank; i++) result[i] = 0.0;
  if (method == Average) for (i = 0; i < *rank; i++) result[*rank*2+i] = 0.0;

  /* loop over all of the elements */
  for (bIndex = 1; bIndex <= discr->nBodys; bIndex++) {
    for (i = 0; i < discr->bodys[bIndex-1].nElems; i++) {
      status = aim_Integration(problem->aimFPTR, analysis->loadName, discr,
                               object->name, bIndex, i+1, *rank, dataset->data,
                               &result[*rank]);
      if (status != CAPS_SUCCESS) {
        printf(" caps_integrateData Warning: status = %d for %s/%s!\n",
               status, analysis->loadName, object->name);
        continue;
      }
      for (j = 0; j < *rank; j++) result[j] += result[*rank+j];
      if (method != Average) continue;
      status = aim_Integration(problem->aimFPTR, analysis->loadName, discr,
                               object->name, bIndex, i+1, *rank, NULL,
                               &result[*rank]);
      if (status != CAPS_SUCCESS) {
        printf(" caps_integrateData Warning: Status = %d for %s/%s!\n",
               status, analysis->loadName, object->name);
        continue;
      }
      for (j = 0; j < *rank; j++) result[*rank*2+j] += result[*rank+j];
    }
  }

  /* make weigted average, if called for */
  if (method == Average) {
    for (j = 0; j < *rank; j++) result[j] /= result[*rank*2+j];
    *units = EG_strdup(dataset->units);
  }

  *data = result;
  return CAPS_SUCCESS;
}


int
caps_boundInfo(capsObject *object, enum capsState *state, int *dim,
               double *plims)
{
  int         status, ret;
  double      *reals;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsBound   *bound;
  capsJrnl    args[3];

  *state = Empty;
  *dim   = 0;
  if (object              == NULL)         return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (object->type        != BOUND)        return CAPS_BADTYPE;
  if (object->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(object, CAPS_BOUNDINFO, &pobject);
  if (status              != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;

  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jPointer;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_BOUNDINFO, problem, object, 3, args, &sNum,
                               &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *state   = args[0].members.integer;
        *dim     = args[1].members.integer;
        reals    = args[2].members.pointer;
        if (*dim >= 1) {
          plims[0] = reals[0];
          plims[1] = reals[1];
        }
        if (*dim == 2) {
          plims[0] = reals[2];
          plims[1] = reals[3];
        }
      }
      return ret;
    }
  }

  bound  = (capsBound *) object->blind;
  *dim   = bound->dim;
  *state = bound->state;
  if (*dim >= 1) {
    plims[0] = bound->plimits[0];
    plims[1] = bound->plimits[1];
    args[2].length = 2*sizeof(double);
  }
  if (*dim == 2) {
    plims[2] = bound->plimits[2];
    plims[3] = bound->plimits[3];
    args[2].length = 4*sizeof(double);
  }
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  args[0].members.integer = *state;
  args[1].members.integer = *dim;
  args[2].members.pointer = plims;
  caps_jrnlWrite(CAPS_BOUNDINFO, problem, object, CAPS_SUCCESS, 3, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}


int
caps_makeBound(capsObject *pobject, int dim, const char *bname,
               capsObject **bobj)
{
  int         i, j, status, ret, *index;
  char        filename[PATH_MAX], temp[PATH_MAX];
  CAPSLONG    sNum;
  capsProblem *problem;
  capsBound   *bound, *bnd;
  capsObject  *object, **tmp;
  capsJrnl    args[1];
  FILE        *fp;

  *bobj = NULL;
  if (pobject == NULL)                   return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)          return CAPS_BADTYPE;
  if (pobject->blind == NULL)            return CAPS_NULLBLIND;
  if (bname == NULL)                     return CAPS_NULLNAME;
  if ((dim < 1) || (dim > 3))            return CAPS_RANGEERR;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)              return CAPS_READONLYERR;

  args[0].type = jObject;
  status       = caps_jrnlRead(CAPS_MAKEBOUND, problem, *bobj, 1, args, &sNum,
                               &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) *bobj = args[0].members.obj;
    return ret;
  }

  /* same name? */
  for (i = 0; i < problem->nBound; i++) {
    if (problem->bounds[i]       == NULL) continue;
    if (problem->bounds[i]->name == NULL) continue;
    if (strcmp(bname, problem->bounds[i]->name) == 0) {
      status = CAPS_BADNAME;
      goto bout;
    }
  }

  sNum  = problem->sNum;
  bound = (capsBound *) EG_alloc(sizeof(capsBound));
  if (bound == NULL) {
    status = EGADS_MALLOC;
    goto bout;
  }
  bound->dim        = dim;
  bound->state      = Open;
  bound->lunits     = NULL;
  bound->plimits[0] = bound->plimits[1] = bound->plimits[2] =
                      bound->plimits[3] = 0.0;
  bound->geom       = NULL;
  bound->iBody      = 0;
  bound->iEnt       = 0;
  bound->curve      = NULL;
  bound->surface    = NULL;
  bound->index      = problem->mBound+1;
  bound->nVertexSet = 0;
  bound->vertexSet  = NULL;
  if (bound->index > 9999) {
    index = (int *) EG_alloc(9998*sizeof(int));
    if (index == NULL) {
      status = EGADS_MALLOC;
      goto bout;
    }
    for (i = 0; i < 9998; i++) index[i] = 0;
    for (i = 0; i < problem->nBound; i++) {
      if (problem->bounds[i]        == NULL) continue;
      if (problem->bounds[i]->blind == NULL) continue;
      bnd = (capsBound *) problem->bounds[i]->blind;
      if ((bnd->index < 1) || (bnd->index > 9999)) continue;
      index[bnd->index-1]++;
    }
    for (i = 0; i < 9998; i++)
      if (index[i] == 0) {
        bound->index = i+1;
        break;
      }
    EG_free(index);
    if (i == 9998) {
      status = CAPS_BADINDEX;
      goto bout;
    }
  }

  /* make the object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(bound);
    goto bout;
  }

  if (problem->bounds == NULL) {
    problem->bounds = (capsObject  **) EG_alloc(sizeof(capsObject *));
    if (problem->bounds == NULL) {
      EG_free(object);
      EG_free(bound);
      status = EGADS_MALLOC;
      goto bout;
    }
  } else {
    tmp = (capsObject **) EG_reall( problem->bounds,
                                   (problem->nBound+1)*sizeof(capsObject *));
    if (tmp == NULL) {
      EG_free(object);
      EG_free(bound);
      status = EGADS_MALLOC;
      goto bout;
    }
    problem->bounds = tmp;
  }

  object->parent = pobject;
  object->name   = EG_strdup(bname);
  object->type   = BOUND;
  object->blind  = bound;

  problem->bounds[problem->nBound] = object;
  if (problem->mBound <= 9999) problem->mBound += 1;
  problem->nBound  += 1;
  problem->sNum    += 1;
  object->last.sNum = problem->sNum;
  status = caps_addHistory(object, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_makeBound)\n", status);

  /* setup for restarts */
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\bound.txt", problem->root);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx",  problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/bound.txt",   problem->root);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",    problem->root);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Warning: Cannot open %s (caps_makeBound)\n", filename);
  } else {
    fprintf(fp, "%d %d\n", problem->nBound, problem->mBound);
    if (problem->bounds != NULL)
      for (i = 0; i < problem->nBound; i++) {
        bnd = (capsBound *) problem->bounds[i]->blind;
        j   = 0;
        if (bnd != NULL) j = bnd->index;
        fprintf(fp, "%d %s\n", j, problem->bounds[i]->name);
      }
    fclose(fp);
    status = caps_rename(temp, filename);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: Cannot rename %s!\n", filename);
  }
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d",
           problem->root, bound->index);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d",
           problem->root, bound->index);
#endif
  status = caps_mkDir(filename);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: Cant make dir %s (caps_makeBound)\n", filename);

  *bobj  = object;
  status = CAPS_SUCCESS;

bout:
  args[0].members.obj = *bobj;
  caps_jrnlWrite(CAPS_MAKEBOUND, problem, *bobj, status, 1, args, sNum,
                 problem->sNum);

  return status;
}


int
caps_closeBound(capsObject *bobject)
{
  int           i, j, status;
  capsObject    *pobject;
  capsProblem   *problem;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;

  if (bobject              == NULL)         return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)        return CAPS_BADTYPE;
  if (bobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(bobject, CAPS_CLOSEBOUND, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)            return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR)   return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)               return CAPS_SUCCESS;
  }

  /* do we have any entries? */
  bound = (capsBound *) bobject->blind;
  if (bound->state         != Open)         return CAPS_STATEERR;
  if (bound->nVertexSet    == 0)            return CAPS_NOTFOUND;

  /* are the VertexSets OK? */
  for (i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL)      return CAPS_NULLOBJ;
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (bound->vertexSet[i]->type        != VERTEXSET) return CAPS_BADTYPE;
    if (bound->vertexSet[i]->blind       == NULL)      return CAPS_NULLBLIND;
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    for (j = 0; j < vertexset->nDataSets; j++) {
      if (vertexset->dataSets[j]              == NULL)      return CAPS_NULLOBJ;
      if (vertexset->dataSets[j]->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (vertexset->dataSets[j]->type        != DATASET)   return CAPS_BADTYPE;
      if (vertexset->dataSets[j]->blind       == NULL)      return CAPS_NULLBLIND;
    }
  }

  /* do all dependent DataSets have a link? */
  for (i = 0; i < bound->nVertexSet; i++) {
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    for (j = 0; j < vertexset->nDataSets; j++) {
      dataset = (capsDataSet *) vertexset->dataSets[j]->blind;
      if (dataset->ftype != FieldIn) continue;
      if (dataset->link == NULL) {
        printf(" caps_closeBound: No link for VertexSet %s, DataSet %s!\n",
               bound->vertexSet[i]->name, vertexset->dataSets[j]->name);
        return CAPS_SOURCEERR;
      }
    }
  }

  bound->state = Empty;
  status = caps_dumpBound(pobject, bobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_dumpBound = %d (caps_closeBound)!\n",
           status);

  return CAPS_SUCCESS;
}


static int
caps_makeDataSeX(capsObject *vobject, const char *dname, enum capsfType ftype,
                 int rank, capsObject **dobj, int *nErr, capsErrs **errors)
{
  int           i, j, status, len, type, nrow, ncol;
  int           open = 0, irow = 1, icol = 1;
  char          *pname = NULL, name[MAX_NAME_LEN], *hash, temp[PATH_MAX];
  capsObject    *bobject, *pobject, *aobject, *object, **tmp;
  capsVertexSet *vertexset;
  capsAnalysis  *analysis;
  capsValue     *value;
  capsProblem   *problem;
  capsDataSet   *dataset;

  *dobj    = NULL;
  *nErr    = 0;
  *errors  = NULL;
  if (dname                == NULL)      return CAPS_NULLNAME;
  if (vobject              == NULL)      return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET) return CAPS_BADTYPE;
  if (vobject->blind       == NULL)      return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vobject->blind;
  aobject   = vertexset->analysis;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis  = (capsAnalysis *) aobject->blind;
  bobject   = vobject->parent;
  if (bobject              == NULL)      return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)     return CAPS_BADTYPE;
  if (bobject->blind       == NULL)      return CAPS_NULLBLIND;
  pobject   = bobject->parent;
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;

  /* is that a legal name? */

  if (ftype == BuiltIn) {
    if ((strcmp(dname, "xyz")  != 0) && (strcmp(dname, "param")  != 0) &&
        (strcmp(dname, "xyzd") != 0) && (strcmp(dname, "paramd") != 0))
      return CAPS_BADDSETNAME;
  } else if ((ftype == GeomSens) || (ftype == TessSens)) {
    rank = 3;
    len = strlen(dname);
    for (j = 1; j < len; j++)
      if (dname[j] == '[') {
        open = j;
        break;
      }
    if (open != 0) {
      pname = (char *) malloc((len+1)*sizeof(char));
      if (pname == NULL) return EGADS_MALLOC;
      for (i = 0; i < len; i++) pname[i] = dname[i];
      pname[open] = 0;
      pname[len]  = 0;
      for (i = open+1; i < len; i++) {
        if (pname[i] == ' ') {
          EG_free(pname);
          return CAPS_BADNAME;
        }
        if (pname[i] == ',') pname[i] = ' ';
      }
      sscanf(&pname[open+1], "%d%d\n", &irow, &icol);
    } else {
      pname = (char *) dname;
    }
    for (i = 0; i < problem->nGeomIn; i++)
      if (strcmp(pname,problem->geomIn[i]->name) == 0) {
        value  = (capsValue *) problem->geomIn[i]->blind;
        status = ocsmGetPmtr(problem->modl, value->pIndex, &type, &nrow, &ncol,
                             name);
        if (status < SUCCESS) {
          snprintf(temp, PATH_MAX, "%s ocsmGetPmtr = %d (caps_makeDataSet)!",
                   dname, status);
          caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
          if (*errors != NULL) *nErr = (*errors)->nError;
          if (open != 0) EG_free(pname);
          return status;
        }
        if (type != OCSM_DESPMTR) {
          snprintf(temp, PATH_MAX, "%s is NOT a Design Parameter (caps_makeDataSet)!",
                   dname);
          caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
          if (*errors != NULL) *nErr = (*errors)->nError;
          if (open != 0) EG_free(pname);
          return CAPS_NOSENSITVTY;
        }
        if ((irow != 1) || (icol != 1)) {
          if (value == NULL) {
            if (open != 0) EG_free(pname);
            return CAPS_NULLVALUE;
          }
          if ((irow < 1) || (irow > value->nrow) ||
              (icol < 1) || (icol > value->ncol)) {
            if (open != 0) EG_free(pname);
            return CAPS_BADINDEX;
          }
        }
        break;
      }
    if (open != 0) EG_free(pname);
    if (i == problem->nGeomIn) {
      snprintf(temp, PATH_MAX, "%s NOT match GeometryInput (caps_makeDataSet)!",
               dname);
      caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_BADNAME;
    }
  } else {
    for (i = 0; i < problem->nGeomIn; i++)
      if (strcmp(dname,problem->geomIn[i]->name) == 0) {
        snprintf(temp, PATH_MAX, "%s matches GeometryInput (caps_makeDataSet)!",
                 dname);
        caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return CAPS_BADNAME;
      }

    if (ftype == FieldIn) {
      rank = -1;
      for (i = 0; i < analysis->nField; i++) {
        if (analysis->fInOut[i] == FieldOut) continue;
        if (strcmp(analysis->fields[i], dname) == 0) {
          rank = analysis->ranks[i];
          break;
        }
        /* check for a name with a numeric wild card */
        hash = strchr(analysis->fields[i],'#');
        if (hash != NULL) {
          strcpy(name, analysis->fields[i]);
          strcpy(name + (int)(hash-analysis->fields[i])/sizeof(char), "%d");
          status = sscanf(dname, name, &j);
          if (status == 1) {
            rank = analysis->ranks[i];
            break;
          }
        }
      }
      if (rank == -1) {
        snprintf(temp, PATH_MAX, "Analysis '%s' does not have a FieldIn '%s'!",
                 aobject->name, dname);
        caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return CAPS_BADNAME;
      }
    } else if (ftype == FieldOut) {
      for (i = 0; i < analysis->nField; i++) {
        if (analysis->fInOut[i] == FieldIn) continue;
        if (strcmp(analysis->fields[i], dname) == 0) {
          rank = analysis->ranks[i];
          break;
        }
        /* check for a name with a numeric wild card */
        hash = strchr(analysis->fields[i],'#');
        if (hash != NULL) {
          strcpy(name, analysis->fields[i]);
          strcpy(name + (int)(hash-analysis->fields[i])/sizeof(char), "%d");
          status = sscanf(dname, name, &j);
          if (status == 1) {
            rank = analysis->ranks[i];
            break;
          }
        }
      }
      if (rank == -1) {
        snprintf(temp, PATH_MAX, "Analysis '%s' does not have a FieldOut '%s'!",
                 aobject->name, dname);
        caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return CAPS_BADNAME;
      }
    } else if (ftype != User) {
      snprintf(temp, PATH_MAX, "Unknown Field Type ftype = %d (caps_makeDataSet)!",
               ftype);
      caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_BADTYPE;
    }
  }

  /* is this name unique? */
  for (i = 0; i < vertexset->nDataSets; i++)
    if (strcmp(dname,vertexset->dataSets[i]->name) == 0) {
      snprintf(temp, PATH_MAX, "%s is already registered (caps_makeDataSet)!",
               dname);
      caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      if (open != 0) EG_free(pname);
      return CAPS_BADNAME;
    }

  status = caps_isNameOK(dname);
  if (status != CAPS_SUCCESS) {
    snprintf(temp, PATH_MAX, "%s has illegal characters (caps_makeDataSet)!",
           dname);
    caps_makeSimpleErr(vobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return status;
  }

  /* fill in the dataset data */
  dataset = (capsDataSet *) EG_alloc(sizeof(capsDataSet));
  if (dataset == NULL) return EGADS_MALLOC;
  dataset->ftype      = ftype;
  dataset->npts       = 0;
  dataset->rank       = rank;
  dataset->data       = NULL;
  dataset->units      = NULL;
  dataset->startup    = NULL;
  dataset->linkMethod = Interpolate;
  dataset->link       = NULL;

  /* make the object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(dataset);
    return status;
  }

  if (vertexset->dataSets == NULL) {
    vertexset->dataSets = (capsObject  **) EG_alloc(sizeof(capsObject *));
    if (vertexset->dataSets == NULL) {
      EG_free(object);
      EG_free(dataset);
      return EGADS_MALLOC;
    }
  } else {
    tmp = (capsObject **) EG_reall( vertexset->dataSets,
                                   (vertexset->nDataSets+1)*sizeof(capsObject *));
    if (tmp == NULL) {
      EG_free(object);
      EG_free(dataset);
      return EGADS_MALLOC;
    }
    vertexset->dataSets = tmp;
  }

  object->parent = vobject;
  object->name   = EG_strdup(dname);
  object->type   = DATASET;
  object->blind  = dataset;

  vertexset->dataSets[vertexset->nDataSets] = object;
  vertexset->nDataSets += 1;
  object->last.sNum     = 0;
  status = caps_addHistory(object, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_makeDataSet)\n", status);

  *dobj = object;
  return CAPS_SUCCESS;
}


int
caps_makeDataSet(capsObject *vobject, const char *dname, enum capsfType ftype,
                 int rank, capsObject **dobj, int *nErr, capsErrs **errors)
{
  int           status, ret;
  CAPSLONG      sNum;
  capsObject    *bobject, *pobject, *aobject;
  capsVertexSet *vertexset;
  capsBound     *bound;
  capsProblem   *problem;
  capsJrnl      args[1];

  if (nErr                 == NULL)      return CAPS_NULLVALUE;
  if (errors               == NULL)      return CAPS_NULLVALUE;
  *dobj    = NULL;
  *nErr    = 0;
  *errors  = NULL;
  if (dname                == NULL)      return CAPS_NULLNAME;
  if (vobject              == NULL)      return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET) return CAPS_BADTYPE;
  if (vobject->blind       == NULL)      return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vobject->blind;
  aobject   = vertexset->analysis;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  bobject   = vobject->parent;
  if (bobject              == NULL)      return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)     return CAPS_BADTYPE;
  if (bobject->blind       == NULL)      return CAPS_NULLBLIND;
  pobject   = bobject->parent;
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADTYPE;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)              return CAPS_READONLYERR;

  args[0].type = jObject;
  status       = caps_jrnlRead(CAPS_MAKEDATASET, problem, *dobj, 1, args, &sNum,
                               &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) *dobj = args[0].members.obj;
    return ret;
  }

  ret   = CAPS_SUCCESS;
  sNum  = problem->sNum;
  bound = (capsBound *) bobject->blind;
  if (bound->state != Open) ret = CAPS_STATEERR;
  if (ret == CAPS_SUCCESS) {
    ret  = caps_makeDataSeX(vobject, dname, ftype, rank, dobj, nErr, errors);
    args[0].members.obj = *dobj;
  }
  caps_jrnlWrite(CAPS_MAKEDATASET, problem, *dobj, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_dataSetInfo(capsObject *dobject, enum capsfType *ftype, capsObject **link,
                 enum capsdMethod *dmeth)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsDataSet *dataset;
  capsJrnl    args[3];

  *link  = NULL;
  *ftype = BuiltIn;
  if  (dobject              == NULL)         return CAPS_NULLOBJ;
  if  (dobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if  (dobject->type        != DATASET)      return CAPS_BADTYPE;
  if  (dobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(dobject, CAPS_DATASETINFO, &pobject);
  if (status                != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;

  args[0].type = jInteger;
  args[1].type = jObject;
  args[2].type = jInteger;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_DATASETINFO, problem, dobject, 3, args,
                               &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *ftype = args[0].members.integer;
        *link  = args[1].members.obj;
        *dmeth = args[2].members.integer;
      }
      return ret;
    }
  }

  dataset = (capsDataSet *) dobject->blind;
  *ftype  = args[0].members.integer = dataset->ftype;
  *link   = args[1].members.obj     = dataset->link;
  *dmeth  = args[2].members.integer = dataset->linkMethod;
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  caps_jrnlWrite(CAPS_DATASETINFO, problem, dobject, CAPS_SUCCESS, 3, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}


int
caps_linkDataSet(capsObject *link, enum capsdMethod method,
                 capsObject *target, int *nErr, capsErrs **errors)
{
  int         status;
  char        temp[PATH_MAX];
  capsDataSet *tgtD, *srcD;
  capsObject  *tgtV, *srcV;
  capsObject  *pobject;
  capsProblem *problem;
  capsBound   *bound;

  if (nErr              == NULL)                return CAPS_NULLVALUE;
  if (errors            == NULL)                return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  status  = caps_findProblem(target, CAPS_LINKDATASET, &pobject);
  if (status != CAPS_SUCCESS)                   return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)                     return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR)       return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)                   return CAPS_SUCCESS;
  }

  /* look at link */
  if (link              == NULL)                return CAPS_NULLOBJ;
  if (link->magicnumber != CAPSMAGIC)           return CAPS_BADOBJECT;
  if (link->type        != DATASET)             return CAPS_BADTYPE;
  if (link->blind       == NULL)                return CAPS_NULLBLIND;
  srcD = (capsDataSet *) link->blind;
  if ((srcD->ftype != FieldOut) &&
      (srcD->ftype != User))                    return CAPS_BADTYPE;

  /* look at target */
  if (target              == NULL)              return CAPS_NULLOBJ;
  if (target->magicnumber != CAPSMAGIC)         return CAPS_BADOBJECT;
  if (target->type        != DATASET)           return CAPS_BADTYPE;
  if (target->blind       == NULL)              return CAPS_NULLBLIND;
  tgtD = (capsDataSet *) target->blind;
  if (tgtD->ftype != FieldIn)                   return CAPS_BADTYPE;

  /* check for compatability */
  if (srcD->rank != tgtD->rank)                 return CAPS_RANGEERR;

  if (link->parent                == NULL)      return CAPS_NULLOBJ;
  if (link->parent->magicnumber   != CAPSMAGIC) return CAPS_BADOBJECT;
  if (link->parent->type          != VERTEXSET) return CAPS_BADTYPE;
  srcV = link->parent;
  if (srcV->parent                == NULL)      return CAPS_NULLOBJ;
  if (srcV->parent->magicnumber   != CAPSMAGIC) return CAPS_BADOBJECT;
  if (srcV->parent->type          != BOUND)     return CAPS_BADTYPE;

  if (target->parent              == NULL)      return CAPS_NULLOBJ;
  if (target->parent->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (target->parent->type        != VERTEXSET) return CAPS_BADTYPE;
  tgtV = target->parent;
  if (tgtV->parent                == NULL)      return CAPS_NULLOBJ;
  if (tgtV->parent->magicnumber   != CAPSMAGIC) return CAPS_BADOBJECT;
  if (tgtV->parent->type          != BOUND)     return CAPS_BADTYPE;

  /* check that the bound object is the same */
  if (srcV->parent != tgtV->parent) {
    snprintf(temp, PATH_MAX, "link (%s) and target (%s) bound missmatch!\n",
             srcV->parent->name, tgtV->parent->name);
    caps_makeSimpleErr(tgtV, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_BADTYPE;
  }
  bound = (capsBound *) tgtV->parent->blind;
  if (bound->state                != Open)      return CAPS_STATEERR;

  /* set the link */
  tgtD->linkMethod  = method;
  tgtD->link        = link;

  /* look for circular links in auto execution */
  status = caps_circularAutoExecs(target, NULL);
  if (status != CAPS_SUCCESS) {
    tgtD->linkMethod  = Copy;
    tgtD->link        = NULL;
    return status;
  }

  caps_freeOwner(&target->last);
  problem->sNum    += 1;
  target->last.sNum = problem->sNum;
  status = caps_addHistory(target, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_linkDataSet)\n", status);

  status = caps_writeSerialNum(problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeSerialNum = %d (caps_linkDataSet)!\n",
           status);

  return CAPS_SUCCESS;
}


static int
caps_makeVertexSeX(capsObject *bobject, /*@null@*/ capsObject *aobject,
                   const char *name, capsObject **vobj,
                   int *nErr, capsErrs **errors, capsObject **ds)
{
  int           i, status;
  capsProblem   *problem;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsObject    *object, *pobject, **tmp;

  bound   = (capsBound *)   bobject->blind;
  pobject = bobject->parent;
  problem = (capsProblem *) pobject->blind;

  /* unique name? */
  for (i = 0; i < bound->nVertexSet; i++)
    if (strcmp(name, bound->vertexSet[i]->name) == 0) return CAPS_BADNAME;

  status = caps_isNameOK(name);
  if (status != CAPS_SUCCESS) return status;

  vertexset = (capsVertexSet *) EG_alloc(sizeof(capsVertexSet));
  if (vertexset == NULL) return EGADS_MALLOC;
  vertexset->analysis  = aobject;
  vertexset->discr     = NULL;
  vertexset->nDataSets = 0;
  vertexset->dataSets  = NULL;

  /* make the object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(vertexset);
    return status;
  }

  if (bound->vertexSet == NULL) {
    bound->vertexSet = (capsObject  **) EG_alloc(sizeof(capsObject *));
    if (bound->vertexSet == NULL) {
      EG_free(object);
      EG_free(vertexset);
      return EGADS_MALLOC;
    }
  } else {
    tmp = (capsObject **) EG_reall( bound->vertexSet,
                                   (bound->nVertexSet+1)*sizeof(capsObject *));
    if (tmp == NULL) {
      EG_free(object);
      EG_free(vertexset);
      return EGADS_MALLOC;
    }
    bound->vertexSet = tmp;
  }

  vertexset->discr = (capsDiscr *) EG_alloc(sizeof(capsDiscr));
  if (vertexset->discr == NULL) {
    EG_free(object);
    EG_free(vertexset);
    return EGADS_MALLOC;
  }
  caps_initDiscr(vertexset->discr);
  vertexset->discr->dim       = bound->dim;

/*@-kepttrans@*/
  object->parent  = bobject;
/*@+kepttrans@*/
  object->name    = EG_strdup(name);
  object->type    = VERTEXSET;
  object->subtype = UNCONNECTED;
  object->blind   = vertexset;

  status = caps_makeDataSeX(object, "xyz",  BuiltIn, 3, &ds[0], nErr, errors);
  if (status != CAPS_SUCCESS) {
    EG_free(vertexset->discr);
    EG_free(object->blind);
    EG_free(object);
    return EGADS_MALLOC;
  }
  if (aobject != NULL) {
    analysis = (capsAnalysis *) aobject->blind;
    if (analysis != NULL) {
      vertexset->discr->instStore =  analysis->instStore;
      vertexset->discr->aInfo     = &analysis->info;
    }
    object->subtype = CONNECTED;
    status = caps_makeDataSeX(object, "xyzd", BuiltIn, 3, &ds[1], nErr, errors);
    if (status != CAPS_SUCCESS) {
      EG_free(ds[0]->name);
      EG_free(ds[0]->blind);
      EG_free(vertexset->dataSets);
      EG_free(vertexset->discr);
      EG_free(object->blind);
      EG_free(object);
      return EGADS_MALLOC;
    }
  }
  if (bound->dim != 3) {
    status = caps_makeDataSeX(object, "param",  BuiltIn, bound->dim, &ds[2],
                              nErr, errors);
    if (status != CAPS_SUCCESS) {
      if (aobject != NULL) {
        EG_free(ds[1]->name);
        EG_free(ds[1]->blind);
      }
      EG_free(ds[0]->name);
      EG_free(ds[0]->blind);
      EG_free(vertexset->dataSets);
      EG_free(vertexset->discr);
      EG_free(object->blind);
      EG_free(object);
      return EGADS_MALLOC;
    }
    if (aobject != NULL) {
      status = caps_makeDataSeX(object, "paramd", BuiltIn, bound->dim, &ds[3],
                                nErr, errors);
      if (status != CAPS_SUCCESS) {
        EG_free(ds[2]->name);
        EG_free(ds[2]->blind);
        EG_free(ds[1]->name);
        EG_free(ds[1]->blind);
        EG_free(ds[0]->name);
        EG_free(ds[0]->blind);
        EG_free(vertexset->dataSets);
        EG_free(vertexset->discr);
        EG_free(object->blind);
        EG_free(object);
        return EGADS_MALLOC;
      }
    }
  }

  bound->vertexSet[bound->nVertexSet] = object;
  bound->nVertexSet += 1;
  object->last.sNum  = 0;
  status = caps_addHistory(object, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_makeVertexSet)\n",
           status);

  *vobj = object;
  return CAPS_SUCCESS;
}


int
caps_makeVertexSet(capsObject *bobject, /*@null@*/ capsObject *aobject,
                   /*@null@*/ const char *vname, capsObject **vobj,
                   int *nErr, capsErrs **errors)
{
  int         status, ret;
  const char  *name;
  CAPSLONG    sNum;
  capsBound   *bound;
  capsProblem *problem;
  capsObject  *pobject, *ds[4] = {NULL, NULL, NULL, NULL};
  capsJrnl    args[1];

  if (nErr                 == NULL)         return CAPS_NULLVALUE;
  if (errors               == NULL)         return CAPS_NULLVALUE;
  *vobj = NULL;
  *nErr = 0;
  *errors = NULL;
  if (bobject              == NULL)         return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)        return CAPS_BADTYPE;
  if (bobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(bobject, CAPS_MAKEVERTEXSET, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)            return CAPS_READONLYERR;

  name = vname;
  if (aobject != NULL) {
    /* connected vertex set */
    if (aobject              == NULL)      return CAPS_NULLOBJ;
    if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
    if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
    if (vname == NULL) name = aobject->name;
  }
  if (name                   == NULL)      return CAPS_NULLNAME;

  args[0].type = jObject;
  status       = caps_jrnlRead(CAPS_MAKEVERTEXSET, problem, *vobj, 1, args,
                               &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) *vobj = args[0].members.obj;
    return ret;
  }

  ret   = CAPS_SUCCESS;
  sNum  = problem->sNum;
  bound = (capsBound *) bobject->blind;
  if (bound->state != Open) ret = CAPS_STATEERR;
  if (ret == CAPS_SUCCESS) {
    ret = caps_makeVertexSeX(bobject, aobject, name, vobj, nErr, errors, ds);
    args[0].members.obj = *vobj;
  }
  caps_jrnlWrite(CAPS_MAKEVERTEXSET, problem, *vobj, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_vertexSetInfo(capsObject *vobject, int *nGpts, int *nDpts,
                   capsObject **bobj, capsObject **aobj)
{
  int           status, ret;
  CAPSLONG      sNum;
  capsProblem   *problem;
  capsObject    *pobject;
  capsVertexSet *vertexset;
  capsJrnl      args[2];

  *nGpts = *nDpts = 0;
  *bobj  = *aobj  = NULL;
  if (vobject              == NULL)         return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET)    return CAPS_BADTYPE;
  if (vobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(vobject, CAPS_VERTEXSETINFO, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem   = (capsProblem *)   pobject->blind;
  vertexset = (capsVertexSet *) vobject->blind;

  *bobj = vobject->parent;
  *aobj = vertexset->analysis;
  args[0].type = jInteger;
  args[1].type = jInteger;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_VERTEXSETINFO, problem, vobject, 2, args,
                               &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *nGpts = args[0].members.integer;
        *nDpts = args[1].members.integer;
      }
      return ret;
    }
  }

  if (vertexset->discr != NULL) {
    *nGpts = vertexset->discr->nPoints;
    *nDpts = vertexset->discr->nVerts;
  }
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  args[0].members.integer = *nGpts;
  args[1].members.integer = *nDpts;
  caps_jrnlWrite(CAPS_VERTEXSETINFO, problem, vobject, CAPS_SUCCESS, 2, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}


int
caps_fillUnVertexSet(capsObject *vobject, int npts, const double *xyzs)
{
  int           i, status;
  capsDataSet   *dataset;
  capsVertexSet *vertexset;
  capsObject    *bobject, *pobject, *dobject = NULL;
  capsProblem   *problem;

  if (vobject              == NULL)       return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET)  return CAPS_BADTYPE;
  if (vobject->blind       == NULL)       return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vobject->blind;
  if (vertexset->dataSets  == NULL)       return CAPS_BADMETHOD;
  dobject   = vertexset->dataSets[0];
  if (dobject->blind       == NULL)       return CAPS_NULLBLIND;
  dataset   = (capsDataSet *) dobject->blind;
  if (vertexset->analysis  != NULL)       return CAPS_NOTCONNECT;
  bobject   = vobject->parent;
  if (bobject              == NULL)       return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)      return CAPS_BADTYPE;
  pobject   = bobject->parent;
  if (pobject              == NULL)       return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)    return CAPS_BADTYPE;
  if (pobject->blind       == NULL)       return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)          return CAPS_READONLYERR;
  problem->funID = CAPS_FILLUNVERTEXSET;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  caps_freeOwner(&vobject->last);
  vobject->last.sNum  = 0;
  status = caps_addHistory(vobject, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_fillUnVertexSet)\n",
           status);

  /* clear this out */
  if (npts <= 0) {
    if (vertexset->discr == NULL) return CAPS_SUCCESS;
    if (vertexset->discr->verts != NULL) EG_free(vertexset->discr->verts);
    EG_free(vertexset->discr);
    vertexset->discr = NULL;
    if (dataset->data != NULL) EG_free(dataset->data);
    dataset->data = NULL;
    dataset->npts = 0;
    return CAPS_SUCCESS;
  }

  if (vertexset->discr->nVerts != npts) {
    vertexset->discr->nVerts = 0;
    if (vertexset->discr->verts != NULL) EG_free(vertexset->discr->verts);
    vertexset->discr->verts = (double *) EG_alloc(3*npts*sizeof(double));
    if (vertexset->discr->verts == NULL) return EGADS_MALLOC;
    vertexset->discr->nVerts = npts;
  }
  if (dataset->npts != npts) {
    dataset->npts = 0;
    if (dataset->data != NULL) EG_free(dataset->data);
    dataset->data = (double *) EG_alloc(3*npts*sizeof(double));
    if (dataset->data != NULL) dataset->npts = npts;
  }

  for (i = 0; i < npts; i++) {
    vertexset->discr->verts[3*i  ] = xyzs[3*i  ];
    vertexset->discr->verts[3*i+1] = xyzs[3*i+1];
    vertexset->discr->verts[3*i+2] = xyzs[3*i+2];
    if (dataset->data == NULL) continue;
    dataset->data[3*i  ] = xyzs[3*i  ];
    dataset->data[3*i+1] = xyzs[3*i+1];
    dataset->data[3*i+2] = xyzs[3*i+2];
  }
  problem->sNum     += 1;
  vobject->last.sNum = problem->sNum;
  status = caps_writeSerialNum(problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Error: caps_writeSerialNum = %d (caps_fillUnVertexSet)\n",
           status);
  status = caps_writeVertexSet(vobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Error: caps_writeVertexSet = %d (caps_fillUnVertexSet)\n",
           status);

  return CAPS_SUCCESS;
}


int
caps_initDataSet(capsObject *dobject, int rank, const double *startup,
                 int *nErr, capsErrs **errors)
{
  int         i, status;
  double      *data;
  capsObject  *pobject;
  capsProblem *problem;
  capsDataSet *dataset;

  if (nErr                 == NULL)         return CAPS_NULLVALUE;
  if (errors               == NULL)         return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if (dobject              == NULL)         return CAPS_NULLOBJ;
  if (dobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (dobject->type        != DATASET)      return CAPS_BADTYPE;
  if (dobject->blind       == NULL)         return CAPS_NULLBLIND;
  if (startup == NULL)                      return CAPS_NULLVALUE;
  dataset = (capsDataSet *) dobject->blind;
  if (dataset->startup     != NULL)         return CAPS_EXISTS;
  status = caps_findProblem(dobject, CAPS_INITDATASET, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)            return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR)   return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)               return CAPS_SUCCESS;
  }

  if (dataset->rank        != rank)         return CAPS_BADRANK;
  if (dataset->ftype       != FieldIn)      return CAPS_BADMETHOD;

  data = (double *) EG_alloc(rank*sizeof(double));
  if (data == NULL) return EGADS_MALLOC;

  for (i = 0; i < rank; i++) data[i] = startup[i];
  dataset->startup = data;

  return CAPS_SUCCESS;
}


int
caps_setData(capsObject *dobject, int nverts, int rank, const double *data,
             const char *units, int *nErr, capsErrs **errors)
{
  int         status;
  char        *uarray = NULL;
  double      *darray;
  capsObject  *vobject, *bobject, *pobject;
  capsDataSet *dataset;
  capsBound   *bound;
  capsProblem *problem;

  if (nErr                 == NULL)       return CAPS_NULLVALUE;
  if (errors               == NULL)       return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if ((nverts <= 0) || (data == NULL))    return CAPS_NULLVALUE;
  if (dobject              == NULL)       return CAPS_NULLOBJ;
  if (dobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (dobject->type        != DATASET)    return CAPS_BADTYPE;
  if (dobject->blind       == NULL)       return CAPS_NULLBLIND;
  dataset = (capsDataSet *) dobject->blind;
  if (dataset->rank        != rank)       return CAPS_BADRANK;
  if (dataset->ftype       != User)       return CAPS_BADTYPE;
  vobject = dobject->parent;
  if (vobject              == NULL)       return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET)  return CAPS_BADTYPE;
  if (vobject->blind       == NULL)       return CAPS_NULLBLIND;
  bobject = vobject->parent;
  if (bobject              == NULL)       return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)      return CAPS_BADTYPE;
  if (bobject->blind       == NULL)       return CAPS_NULLBLIND;
  bound   = (capsBound *) bobject->blind;
  if (bound->state         == Open)       return CAPS_STATEERR;
  pobject = bobject->parent;
  if (pobject              == NULL)       return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)    return CAPS_BADTYPE;
  if (pobject->blind       == NULL)       return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)          return CAPS_READONLYERR;
  problem->funID = CAPS_SETDATA;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  darray = (double *) EG_alloc(rank*nverts*sizeof(double));
  if (darray == NULL) return EGADS_MALLOC;
  if (units != NULL) {
    uarray = EG_strdup(units);
    if (uarray == NULL) {
      EG_free(darray);
      return EGADS_MALLOC;
    }
  }

  if (dataset->units != NULL) EG_free(dataset->units);
  if (dataset->data  != NULL) EG_free(dataset->data);
  dataset->npts  = nverts;
  dataset->rank  = rank;
  dataset->data  = darray;
  dataset->units = uarray;
  caps_freeOwner(&dobject->last);
  problem->sNum     += 1;
  dobject->last.sNum = problem->sNum;
  status = caps_addHistory(dobject, problem);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Warning: caps_addHistory = %d (caps_setData)\n", status);
    status = caps_writeSerialNum(problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Error: caps_writeSerialNum = %d (caps_setData)\n", status);
  }
  status = caps_writeDataSet(dobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Error: caps_writeDataSet = %d (caps_setData)\n", status);

  return CAPS_SUCCESS;
}


/*
 * obj_bar: compute objective function and gradient via reverse differentiation
 */
static int
obj_bar(int    n,                     /* (in)  number of design variables */
        double ftgt[],                /* (in)  design variables */
        void   *blind,                /* (in)  blind pointer to structure */
        double *obj,                  /* (out) objective function */
        double ftgt_bar[])            /* (out) gradient of objective function */
{
  int        status = CAPS_SUCCESS;   /* (out) return status */
  int        idat, bIndex, ibods, ielms, ibodt, ielmt, imat, irank, jrank;
  int        nrank, sindx, tindx;
  char       *name;
  double     f_src, f_tgt;
  double     area_src, area_tgt, area_tgt_bar, obj_bar1;
  double     *result, *result_bar = NULL, *data_bar = NULL;
  capsConFit *cfit = (capsConFit *) blind;

  irank  = cfit->irank;
  nrank  = cfit->nrank;
  sindx  = cfit->sindx;
  tindx  = cfit->tindx;
  name   = cfit->name;
  result = (double *) EG_alloc(nrank*sizeof(double));
  if (result == NULL) return EGADS_MALLOC;

  /* store ftgt into tgt structure */
  for (idat = 0; idat < n; idat++)
    cfit->data_tgt[nrank*idat+irank] = ftgt[idat];

  /* compute the area for srcD */
  area_src = 0.0;

  for (bIndex = 1; bIndex <= cfit->src->nBodys; bIndex++) {
    for (ielms = 0; ielms < cfit->src->bodys[bIndex-1].nElems; ielms++) {
      status = aim_IntegrIndex(*cfit->aimFPTR, sindx, cfit->src, name, bIndex,
                               ielms+1, nrank, cfit->data_src, result);
      if (status != CAPS_SUCCESS) goto cleanup;

      area_src += result[irank];
    }
  }
  cfit->area_src = area_src;

  /* compute the area for tgt */
  area_tgt = 0.0;

  for (bIndex = 1; bIndex <= cfit->tgt->nBodys; bIndex++) {
    for (ielmt = 0; ielmt < cfit->tgt->bodys[bIndex-1].nElems; ielmt++) {
      status = aim_IntegrIndex(*cfit->aimFPTR, tindx, cfit->tgt, name, bIndex,
                               ielmt+1, nrank, cfit->data_tgt, result);
      if (status != CAPS_SUCCESS) goto cleanup;

      area_tgt += result[irank];
    }
  }
  cfit->area_tgt = area_tgt;

  /* penalty function part of objective function */
  *obj = cfit->afact * pow(area_tgt-area_src, 2);

  /* minimize the difference between source and target at the Match points */
  for (imat = 0; imat < cfit->nmat; imat++) {
    ibods  = cfit->mat[imat].source.bIndex;
    ielms  = cfit->mat[imat].source.eIndex;
    ibodt  = cfit->mat[imat].target.bIndex;
    ielmt  = cfit->mat[imat].target.eIndex;
    if ((ielms == -1) || (ielmt == -1)) continue;
    status = aim_InterpolIndex(*cfit->aimFPTR, sindx, cfit->src, name, ibods,
                               ielms, cfit->mat[imat].source.st, nrank,
                               cfit->data_src, result);
    if (status != CAPS_SUCCESS) goto cleanup;
    f_src  = result[irank];

    status = aim_InterpolIndex(*cfit->aimFPTR, tindx, cfit->tgt, name, ibodt,
                               ielmt, cfit->mat[imat].target.st, nrank,
                               cfit->data_tgt, result);
    if (status != CAPS_SUCCESS) goto cleanup;
    f_tgt  = result[irank];

    *obj  += pow(f_tgt-f_src, 2);
  }

  /* if we do not need gradient, return now */
  if (ftgt_bar == NULL) goto cleanup;

  result_bar = (double *) EG_alloc(nrank*sizeof(double));
  if (result_bar == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  data_bar = (double *) EG_alloc(nrank*n*sizeof(double));
  if (data_bar == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  /* initialize the derivatives */
  obj_bar1  = 1.0;
  for (jrank = 0; jrank < nrank; jrank++) {
    for (idat = 0; idat < n; idat++)
      data_bar[nrank*idat+jrank] = 0.0;
    result_bar[jrank] = 0.0;
  }

  /* reverse: minimize the difference between the source and target
              at the Match points */
  for (imat = cfit->nmat-1; imat >= 0; imat--) {
    ibods  = cfit->mat[imat].source.bIndex;
    ielms  = cfit->mat[imat].source.eIndex;
    ibodt  = cfit->mat[imat].target.bIndex;
    ielmt  = cfit->mat[imat].target.eIndex;
    if ((ielms == -1) || (ielmt == -1)) continue;
    status = aim_InterpolIndex(*cfit->aimFPTR, sindx, cfit->src, name, ibods,
                               ielms, cfit->mat[imat].source.st, nrank,
                               cfit->data_src, result);
    if (status != CAPS_SUCCESS) goto cleanup;
    f_src  = result[irank];

    status = aim_InterpolIndex(*cfit->aimFPTR, tindx, cfit->tgt, name, ibodt,
                               ielmt, cfit->mat[imat].target.st, nrank,
                               cfit->data_tgt, result);
    if (status != CAPS_SUCCESS) goto cleanup;
    f_tgt  = result[irank];

    result_bar[irank] = (f_tgt - f_src) * 2 * obj_bar1;

    status = aim_InterpolIndBar(*cfit->aimFPTR, tindx, cfit->tgt, name, ibodt,
                                ielmt, cfit->mat[imat].target.st, nrank,
                                result_bar, data_bar);
    if (status != CAPS_SUCCESS) goto cleanup;
  }

  /* reverse: penalty function part of objective function */
  area_tgt_bar = cfit->afact * 2 * (area_tgt - area_src) * obj_bar1;

  result_bar[irank] = area_tgt_bar;

  /* reverse: compute the area for tgt */
  for (ibodt = cfit->tgt->nBodys-1; ibodt >= 0; ibodt--) {
    for (ielmt = cfit->tgt->bodys[ibodt].nElems-1; ielmt >= 0; ielmt--) {
      status = aim_IntegrIndBar(*cfit->aimFPTR, tindx, cfit->tgt, name, ibodt+1,
                                ielmt+1, nrank, result_bar, data_bar);
      if (status != CAPS_SUCCESS) goto cleanup;
    }
  }

  for (idat = 0; idat < n; idat++)
    ftgt_bar[idat] = data_bar[nrank*idat+irank];

cleanup:
  if (data_bar   != NULL) EG_free(data_bar);
  if (result_bar != NULL) EG_free(result_bar);
  EG_free(result);

  return status;
}


static int
caps_Match(capsConFit *fit, int dim)
{
  int       i, j, k, ib, n, t, npts, stat;
  double    pos[3], *st;
  capsMatch *match;

  /* count the positions */
  npts = 0;
  for (ib = 0; ib < fit->tgt->nBodys; ib++) {
    for (i = 0; i < fit->tgt->bodys[ib].nElems; i++) {
      t = fit->tgt->bodys[ib].elems[i].tIndex - 1;
      if (fit->tgt->types[t].nmat == 0) {
        npts += fit->tgt->types[t].nref;
      } else {
        npts += fit->tgt->types[t].nmat;
      }
    }
  }

  /* set the match locations */
  match = (capsMatch *) EG_alloc(npts*sizeof(capsMatch));
  if (match == NULL) return EGADS_MALLOC;
  for (i = 0; i < npts; i++) {
    match[i].source.bIndex = -1;
    match[i].source.eIndex = -1;
    match[i].source.st[0]  = 0.0;
    match[i].source.st[1]  = 0.0;
    match[i].source.st[2]  = 0.0;
    match[i].target.bIndex = -1;
    match[i].target.eIndex = -1;
    match[i].target.st[0]  = 0.0;
    match[i].target.st[1]  = 0.0;
    match[i].target.st[2]  = 0.0;
  }

  /* set things up in parameter (or 3) space */
  for (ib = 0; ib < fit->tgt->nBodys; ib++) {
    for (k = i = 0; i < fit->tgt->bodys[ib].nElems; i++) {
      t = fit->tgt->bodys[ib].elems[i].tIndex - 1;
      if (fit->tgt->types[t].nmat == 0) {
        n  = fit->tgt->types[t].nref;
        st = fit->tgt->types[t].gst;
      } else {
        n  = fit->tgt->types[t].nmat;
        st = fit->tgt->types[t].matst;
      }
      for (j = 0; j < n; j++, k++) {
        match[i].target.bIndex = ib+1;
        match[i].target.eIndex = i +1;
        if (dim == 3) {
          match[i].target.st[0] = st[0];
          match[i].target.st[1] = st[1];
          match[i].target.st[2] = st[2];
          stat = aim_InterpolIndex(*fit->aimFPTR, fit->tindx, fit->tgt, "xyz",
                                   ib+1, i+1, &st[3*j], 3, fit->prms_tgt, pos);
        } else {
          match[i].target.st[0] = st[0];
          if (dim == 2) match[i].target.st[1] = st[1];
          stat = aim_InterpolIndex(*fit->aimFPTR, fit->tindx, fit->tgt, "param",
                                   ib+1, i+1, &st[dim*j], dim, fit->prms_tgt,
                                   pos);
        }
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Warning: %d/%d aim_Interpolation %d = %d (match)!\n",
                 k, npts, fit->tindx, stat);
          continue;
        }
        stat = aim_LocateElIndex(*fit->aimFPTR, fit->sindx, fit->src,
                                 fit->prms_src, pos,
                                 &match[i].source.bIndex, &match[i].source.eIndex,
                                 match[i].source.st);
        if (stat != CAPS_SUCCESS)
          printf(" CAPS Warning: %d/%d aim_LocateElement = %d (match)!\n",
                 i, npts, stat);
      }
    }
  }

  fit->mat  = match;
  fit->nmat = npts;

  return CAPS_SUCCESS;
}


static int
caps_Conserve(capsConFit *fit, const char *bname, int dim)
{
  int    i, j, stat, *elems;
  double fopt, *ref, *ftgt, *tmp;
  FILE   *fp = NULL;

#ifdef DEBUG
  fp    = stdout;
#endif
  elems = (int *) EG_alloc(2*fit->npts*sizeof(int));
  if (elems == NULL) {
    printf(" CAPS Error: Malloc on %d element indices (caps_getData)!\n",
           fit->npts);
    return EGADS_MALLOC;
  }
  ref = (double *) EG_alloc(((dim+1)*fit->npts+fit->nrank)*sizeof(double));
  if (ref == NULL) {
    EG_free(elems);
    printf(" CAPS Error: Malloc on %d element %d references (caps_getData)!\n",
           fit->npts, dim);
    return EGADS_MALLOC;
  }
  ftgt = &ref[dim*fit->npts];
  tmp  = &ftgt[fit->npts];
  for (i = 0; i < fit->npts; i++) {
    elems[2*i  ] = 0;
    elems[2*i+1] = 0;
    stat     = aim_LocateElIndex(*fit->aimFPTR, fit->sindx, fit->src,
                                 fit->prms_src, &fit->prms_tgt[dim*i],
                                 &elems[2*i], &elems[2*i+1], &ref[dim*i]);
    if (stat != CAPS_SUCCESS)
      printf(" CAPS Warning: %d/%d aim_LocateElement = %d (caps_getData)!\n",
             i, fit->npts, stat);
  }

  stat = 0;
  for (j = 0; j < fit->nrank; j++) {
    fit->irank = j;
    for (i = 0; i < fit->npts; i++) {
      if (elems[2*i] == 0) continue;
      stat = aim_InterpolIndex(*fit->aimFPTR, fit->sindx, fit->src, fit->name,
                               elems[2*i], elems[2*i+1], &ref[dim*i], fit->nrank,
                               fit->data_src, tmp);
      if (stat != CAPS_SUCCESS)
        printf(" CAPS Warning: %d/%d aim_Interpolation = %d (caps_getData)!\n",
               i, fit->npts, stat);
      ftgt[i] = tmp[j];
    }
    stat = caps_conjGrad(obj_bar, fit, fit->npts, ftgt, 1e-6, fp, &fopt);
    if (stat != CAPS_SUCCESS) break;

    if (j == 0)
      printf(" CAPS Info: Bound '%s' Normalized Integrated '%s'\n",
             bname, fit->name);
    printf("            Rank %d: src = %le, tgt = %le, diff = %le\n",
           j, fit->area_src, fit->area_tgt, fabs(fit->area_src-fit->area_tgt));
  }

  EG_free(ref);
  EG_free(elems);
  return stat;
}


int
caps_getTriangles(capsObject *vobject, int *nGtris,            int **gtris,
                                       int *nGsegs, /*@null@*/ int **gsegs,
                                       int *nDtris,            int **dtris,
                                       int *nDsegs, /*@null@*/ int **dsegs)
{
  int           i, j, ib, n, status, ret, ntris, nsegs, eType, *tris, *segs;
  CAPSLONG      sNum;
  capsObject    *pobject;
  capsProblem   *problem;
  capsVertexSet *vertexset;
  capsDiscr     *discr;
  capsJrnl      args[8];

  *nGtris = 0;
  *nGsegs = 0;
  *nDtris = 0;
  *nDsegs = 0;
  *gtris  = NULL;
  *dtris  = NULL;
  if (gsegs != NULL) *gsegs = NULL;
  if (dsegs != NULL) *dsegs = NULL;
  if (vobject              == NULL)         return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET)    return CAPS_BADTYPE;
  if (vobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(vobject, CAPS_GETTRIANGLES, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem  = (capsProblem *) pobject->blind;

  args[0].type   = jInteger;
  args[1].type   = jPtrFree;
  args[1].length = 0;
  args[2].type   = jInteger;
  args[3].type   = jPtrFree;
  args[3].length = 0;
  args[4].type   = jInteger;
  args[5].type   = jPtrFree;
  args[5].length = 0;
  args[6].type   = jInteger;
  args[7].type   = jPtrFree;
  args[7].length = 0;
  if (problem->dbFlag == 0) {
    status       = caps_jrnlRead(CAPS_GETTRIANGLES, problem, vobject, 8, args,
                                 &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *nGtris = args[0].members.integer;
        *gtris  = args[1].members.pointer;
        *nGsegs = args[2].members.integer;
        if (gsegs != NULL) *gsegs = args[3].members.pointer;
        *nDtris = args[4].members.integer;
        *dtris  = args[5].members.pointer;
        *nDsegs = args[6].members.integer;
        if (dsegs != NULL) *dsegs = args[7].members.pointer;
      }
      return ret;
    }
  }

  status    = CAPS_SUCCESS;
  vertexset = (capsVertexSet *) vobject->blind;
  if (vertexset->discr == NULL) goto vout;
  discr = vertexset->discr;
  if ((discr->nPoints == 0) && (discr->nVerts == 0)) goto vout;

  ntris = nsegs = 0;
  for (ib = 0; ib < discr->nBodys; ib++) {
    for (j = 0; j < discr->bodys[ib].nElems; j++) {
      eType = discr->bodys[ib].elems[j].tIndex;
      if (discr->types[eType-1].tris == NULL) {
        ntris++;
        nsegs += 3;
      } else {
        ntris += discr->types[eType-1].ntri;
        nsegs += discr->types[eType-1].nseg;
      }
    }
  }
  if (ntris != 0) {
    tris = (int *) EG_alloc(3*ntris*sizeof(int));
    if (tris == NULL) {
      status = EGADS_MALLOC;
      goto vout;
    }
    for (ib = 0; ib < discr->nBodys; ib++) {
      for (ntris = j = 0; j < discr->bodys[ib].nElems; j++) {
        eType = discr->bodys[ib].elems[j].tIndex;
        if (discr->types[eType-1].tris == NULL) {
          tris[3*ntris  ] = discr->bodys[ib].elems[j].gIndices[0];
          tris[3*ntris+1] = discr->bodys[ib].elems[j].gIndices[2];
          tris[3*ntris+2] = discr->bodys[ib].elems[j].gIndices[4];
          ntris++;
        } else {
          for (i = 0; i < discr->types[eType-1].ntri; i++, ntris++) {
            n = discr->types[eType-1].tris[3*i  ] - 1;
            tris[3*ntris  ] = discr->bodys[ib].elems[j].gIndices[2*n];
            n = discr->types[eType-1].tris[3*i+1] - 1;
            tris[3*ntris+1] = discr->bodys[ib].elems[j].gIndices[2*n];
            n = discr->types[eType-1].tris[3*i+2] - 1;
            tris[3*ntris+2] = discr->bodys[ib].elems[j].gIndices[2*n];
          }
        }
      }
    }
    segs = NULL;
    if ((gsegs != NULL) && (nsegs != 0)) {
      segs = (int *) EG_alloc(2*nsegs*sizeof(int));
      if (segs == NULL) {
        EG_free(tris);
        status = EGADS_MALLOC;
        goto vout;
      }
      for (ib = 0; ib < discr->nBodys; ib++) {
        for (nsegs = j = 0; j < discr->bodys[ib].nElems; j++) {
          eType = discr->bodys[ib].elems[j].tIndex;
          if (discr->types[eType-1].tris == NULL) {
            segs[2*nsegs  ] = discr->bodys[ib].elems[j].gIndices[0];
            segs[2*nsegs+1] = discr->bodys[ib].elems[j].gIndices[2];
            nsegs++;
            segs[2*nsegs  ] = discr->bodys[ib].elems[j].gIndices[2];
            segs[2*nsegs+1] = discr->bodys[ib].elems[j].gIndices[4];
            nsegs++;
            segs[2*nsegs  ] = discr->bodys[ib].elems[j].gIndices[4];
            segs[2*nsegs+1] = discr->bodys[ib].elems[j].gIndices[0];
            nsegs++;
          } else {
            for (i = 0; i < discr->types[eType-1].nseg; i++, nsegs++) {
              n = discr->types[eType-1].segs[2*i  ] - 1;
              segs[2*nsegs  ] = discr->bodys[ib].elems[j].gIndices[2*n];
              n = discr->types[eType-1].segs[2*i+1] - 1;
              segs[2*nsegs+1] = discr->bodys[ib].elems[j].gIndices[2*n];
            }
          }
        }
      }
      *gsegs  = segs;
      args[3].length = 2*nsegs*sizeof(int);
    }
    *nGsegs = nsegs;
    *nGtris = ntris;
    *gtris  = tris;
    args[1].length = 3*ntris*sizeof(int);
  }
  if ((discr->nDtris == 0) || (discr->dtris == NULL) ||
      (discr->nVerts == 0)) goto vout;

  tris = (int *) EG_alloc(3*discr->nDtris*sizeof(int));
  if (tris == NULL) {
    if (gsegs != NULL) EG_free(*gsegs);
    EG_free(*gtris);
    *nGtris        = 0;
    *gtris         = NULL;
    args[1].length = 0;
    *nGsegs        = 0;
    if (gsegs != NULL) *gsegs = NULL;
    args[3].length = 0;
    status         = EGADS_MALLOC;
    goto vout;
  }
  for (j = 0; j < 3*discr->nDtris; j++) tris[j] = discr->dtris[j];
  *nDtris = discr->nDtris;
  *dtris  = tris;
  args[5].length = 3*discr->nDtris*sizeof(int);
  *nDsegs = discr->nDsegs;

  if (dsegs != NULL) {
    segs = (int *) EG_alloc(2*discr->nDsegs*sizeof(int));
    if (segs == NULL) {
      if (gsegs != NULL) EG_free(*gsegs);
      EG_free(*gtris);
      *nGtris        = 0;
      *gtris         = NULL;
      args[1].length = 0;
      *nGsegs        = 0;
      if (gsegs != NULL) *gsegs = NULL;
      args[3].length = 0;
      EG_free(*dtris);
      *nDtris        = 0;
      *dtris         = NULL;
      args[5].length = 0;
      status         = EGADS_MALLOC;
      goto vout;
    }
    for (j = 0; j < 2*discr->nDsegs; j++) segs[j] = discr->dsegs[j];
    *dsegs  = segs;
    args[7].length = 2*discr->nDsegs*sizeof(int);
  }

vout:
  if (problem->dbFlag == 1) return status;
  args[0].members.integer = *nGtris;
  args[1].members.pointer = *gtris;
  args[2].members.integer = *nGsegs;
  args[3].members.pointer = NULL;
  if (gsegs != NULL) args[3].members.pointer = *gsegs;
  args[4].members.integer = *nDtris;
  args[5].members.pointer = *dtris;
  args[6].members.integer = *nDsegs;
  args[7].members.pointer = NULL;
  if (dsegs != NULL) args[7].members.pointer = *dsegs;
  caps_jrnlWrite(CAPS_GETTRIANGLES, problem, vobject, status, 8, args,
                 problem->sNum, problem->sNum);

  return status;
}


int
caps_outputVertexSet(capsObject *vobject, const char *filename)
{
  int           i, j, k, stat, nGtris, *gtris, nGsegs, nDtris, *dtris, nDsegs;
  FILE          *fp;
  capsObject    *pobject;
  capsProblem   *problem;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;

  if (vobject              == NULL)         return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET)    return CAPS_BADTYPE;
  if (vobject->blind       == NULL)         return CAPS_NULLBLIND;
  if (filename             == NULL)         return CAPS_NULLNAME;
  stat = caps_findProblem(vobject, CAPS_OUTPUTVERTEXSET, &pobject);
  if (stat                 != CAPS_SUCCESS) return stat;
  problem  = (capsProblem *) pobject->blind;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR)   return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    stat = caps_jrnlEnd(problem);
    if (stat != CAPS_CLEAN)                 return CAPS_SUCCESS;
  }

  vertexset = (capsVertexSet *) vobject->blind;
  fp = fopen(filename, "w");
  if (fp                   == NULL)      return CAPS_IOERR;
  stat = caps_getTriangles(vobject, &nGtris, &gtris, &nGsegs, NULL,
                                    &nDtris, &dtris, &nDsegs, NULL);
  if (stat != CAPS_SUCCESS) {
    printf(" caps_outputVertexSet Error: caps_getTriangles = %d!\n", stat);
    fclose(fp);
    return stat;
  }
  printf(" **** writing VertexSet file: %s ****\n", filename);
  fprintf(fp, "%s\n", vobject->parent->name);
  fprintf(fp, "%8d %8d %8d\n", nGtris, nDtris, vertexset->nDataSets);
  for (i = 0; i < nGtris; i++)
    fprintf(fp, "    %8d %8d %8d\n", gtris[3*i  ], gtris[3*i+1], gtris[3*i+2]);
  for (i = 0; i < nDtris; i++)
    fprintf(fp, "    %8d %8d %8d\n", dtris[3*i  ], dtris[3*i+1], dtris[3*i+2]);
  EG_free(gtris);
  EG_free(dtris);
  for (i = 0; i < vertexset->nDataSets; i++) {
    fprintf(fp, "%s\n", vertexset->dataSets[i]->name);
    dataset = (capsDataSet *) vertexset->dataSets[i]->blind;
    fprintf(fp, " %8d %8d\n", dataset->npts, dataset->rank);
    for (j = 0; j < dataset->npts; j++) {
      for (k = 0; k < dataset->rank; k++)
        fprintf(fp, " %lf", dataset->data[j*dataset->rank+k]);
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return CAPS_SUCCESS;
}


static int
caps_fillFieldOut(capsObject *dobject, int *nErr, capsErrs **errors)
{
  int           stat, major, minor, nField, exec, dirty, *ranks, *fInOut;
  char          *intents, *apath, *unitSys, **fnames;
  char          temp[PATH_MAX];
  capsObject    *vobject, *bobject, *pobject, *aobject;
  capsDataSet   *dataset;
  capsVertexSet *vertexset;
  capsAnalysis  *analysis = NULL;
  capsProblem   *problem;

  *nErr     = 0;
  *errors   = NULL;
  dataset   = (capsDataSet *) dobject->blind;
  vobject   = dobject->parent;
  vertexset = (capsVertexSet *) vobject->blind;
  aobject   = vertexset->analysis;
  bobject   = vobject->parent;
  pobject   = bobject->parent;
  problem   = (capsProblem *) pobject->blind;

  if (dataset->ftype != FieldOut) return CAPS_SOURCEERR;

  if (aobject == NULL) {
    snprintf(temp, PATH_MAX, "caps_getData DataSet %s with NULL analysis!",
             dobject->name);
    caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_SOURCEERR;
  }

  /* check to see if analysis is clean or can auto execute */
  stat = caps_analysisInfX(aobject, &apath, &unitSys, &major, &minor, &intents,
                           &nField, &fnames, &ranks, &fInOut, &exec, &dirty);
  if (stat != CAPS_SUCCESS) return stat;
  if (dirty > 0) {
    /* auto execute if available */
    if ((exec == 2) && (dirty < 5)) {
      stat = caps_execX(aobject, nErr, errors);
      if (stat != CAPS_SUCCESS) return stat;
    } else {
      return CAPS_DIRTY;
    }
  }

  analysis = (capsAnalysis *) aobject->blind;

  if (((dobject->last.sNum < analysis->pre.sNum) ||
       (dobject->last.sNum == 0) || (dataset->data == NULL))) {

    /* make sure that VertexSets are up-to-date */
    stat = caps_buildBound(bobject, nErr, errors);
    if ((stat != CAPS_CLEAN) && (stat != CAPS_SUCCESS)) return stat;

    if (vertexset->discr == NULL) return CAPS_SOURCEERR;
    dataset->npts = vertexset->discr->nVerts;
    if (dataset->units != NULL) EG_free(dataset->units);
    dataset->units = NULL;
    if (dataset->npts == 0) dataset->npts = vertexset->discr->nPoints;
    if (dataset->npts == 0) return CAPS_SOURCEERR;
    if (dataset->data != NULL) EG_free(dataset->data);
    dataset->data = (double *) EG_alloc(dataset->npts*
                                        dataset->rank*sizeof(double));
    if (dataset->data == NULL) {
      dataset->npts = 0;
      snprintf(temp, PATH_MAX, "caps_getData %s -- DataSet %s Malloc Error!",
               vertexset->analysis->name, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_SOURCEERR;
    }
    stat = aim_Transfer(problem->aimFPTR, analysis->loadName, vertexset->discr,
                        dobject->name, dataset->npts, dataset->rank,
                        dataset->data, &dataset->units);
    if (stat != CAPS_SUCCESS) {
      EG_free(dataset->data);
      dataset->data = NULL;
      dataset->npts = 0;
      snprintf(temp, PATH_MAX, "caps_getData %s -- aimTransfer %s returns %d!",
               vertexset->analysis->name, dobject->name, stat);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_SOURCEERR;
    } else {
      caps_freeOwner(&dobject->last);
      problem->sNum     += 1;
      dobject->last.sNum = problem->sNum;
      stat = caps_addHistory(dobject, problem);
      if (stat != CAPS_SUCCESS)
        printf(" CAPS Warning: caps_addHistory = %d (caps_getData)\n", stat);
      stat = caps_writeDataSet(dobject);
      if (stat != CAPS_SUCCESS)
        printf(" CAPS Warning: caps_writeDataSet = %d (caps_getData)\n", stat);
    }
    if (caps_unitParse(dataset->units) != CAPS_SUCCESS) {
      snprintf(temp, PATH_MAX, "caps_writeDataSet %s -- DataSet %s Units Error!",
               vertexset->analysis->name, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      EG_free(dataset->units);
      dataset->units = NULL;
      return CAPS_UNITERR;
    }
  }

  return CAPS_SUCCESS;
}


static int
caps_fillFieldIn(capsObject *dobject, int *nErr, capsErrs **errors)
{
  int           i, j, bIndex, eIndex, index, stat, npts;
  int           major, minor, nField, exec, dirty, *ranks, *fInOut;
  double        *values, *params, *param, *data, st[3];
  char          *units, temp[PATH_MAX];
  char          *intents, *apath, *unitSys, **fnames;
  capsObject    *vobject, *bobject, *pobject, *aobject, *foundset, *foundanl;
  capsDataSet   *dataset, *otherset, *ds;
  capsVertexSet *vertexset, *fvs = NULL;
  capsDiscr     *discr, *founddsr = NULL;
  capsAnalysis  *analysis = NULL, *fanal;
  capsBound     *bound;
  capsProblem   *problem;
  capsConFit    fit;

  *nErr     = 0;
  *errors   = NULL;
  dataset   = (capsDataSet *) dobject->blind;
  vobject   = dobject->parent;
  vertexset = (capsVertexSet *) vobject->blind;
  discr     = vertexset->discr;
  aobject   = vertexset->analysis;
  bobject   = vobject->parent;
  bound     = (capsBound *) bobject->blind;
  pobject   = bobject->parent;
  problem   = (capsProblem *) pobject->blind;
  if (aobject != NULL) analysis = (capsAnalysis *) aobject->blind;

  otherset  = NULL;
  foundanl  = NULL;

  /* get link to other DataSet */
  if (dataset->link == NULL) {
    snprintf(temp, PATH_MAX, "FieldIn DataSet '%s' without a link (caps_getData)!",
             dobject->name);
    caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_NULLOBJ;
  }

  foundset = dataset->link;
  if (foundset->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (foundset->type        != DATASET)   return CAPS_BADTYPE;
  if (foundset->blind       == NULL)      return CAPS_NULLBLIND;
  otherset = (capsDataSet *) foundset->blind;
  if (foundset->parent              == NULL)      return CAPS_NULLOBJ;
  if (foundset->parent->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (foundset->parent->type        != VERTEXSET) return CAPS_BADTYPE;
  if (foundset->parent->blind       == NULL)      return CAPS_NULLBLIND;
  fvs = (capsVertexSet *) foundset->parent->blind;

  foundanl = fvs->analysis;
  founddsr = fvs->discr;

  if ((fvs == NULL) || (foundset == NULL) || (otherset == NULL) ||
      (founddsr == NULL)) {
    snprintf(temp, PATH_MAX, "Bound %s -- %s with incomplete linked DataSet!",
             bobject->name, dobject->name);
    caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_SOURCEERR;
  }
  if (foundanl == NULL) {
    snprintf(temp, PATH_MAX, "Bound %s -- Analysis is NULL (caps_getData)!",
             bobject->name);
    caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_BADOBJECT;
  }

  /* check to see if analysis is dirty */
  stat = caps_analysisInfX(foundanl, &apath, &unitSys, &major, &minor, &intents,
                           &nField, &fnames, &ranks, &fInOut, &exec, &dirty);
  if (stat != CAPS_SUCCESS) return stat;

  /* can we auto execute or have we been updated? */
  if (((exec == 2) && (dirty > 0) && (dirty < 5)) ||
      (foundanl->last.sNum > dobject->last.sNum) || (dataset->data == NULL)) {

    fanal = (capsAnalysis *) foundanl->blind;
    if (fanal == NULL) {
      snprintf(temp, PATH_MAX, "Bound %s -- Analysis %s blind NULL (caps_getData)!",
               bobject->name, foundanl->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_NULLVALUE;
    }
    index = aim_Index(problem->aimFPTR, fanal->loadName);
    if (index < 0) return index;

    /* make sure that VertexSets are up-to-date */
    stat = caps_buildBound(bobject, nErr, errors);
    if ((stat != CAPS_CLEAN) && (stat != CAPS_SUCCESS)) return stat;

    npts = discr->nPoints;
    if (discr->nVerts != 0) npts = discr->nVerts;
    params = NULL;
    param  = NULL;
    if (bound->dim == 3) {
      for (j = 0; j < vertexset->nDataSets; j++) {
        if (strcmp("xyz", vertexset->dataSets[j]->name) != 0) continue;
        ds    = (capsDataSet *) vertexset->dataSets[j]->blind;
        param = ds->data;
        break;
      }
      for (j = 0; j < fvs->nDataSets; j++) {
        if (strcmp("xyz", fvs->dataSets[j]->name) != 0) continue;
        ds     = (capsDataSet *) fvs->dataSets[j]->blind;
        params = ds->data;
        break;
      }
    } else {
      if (discr->nVerts != 0) {
        for (j = 0; j < vertexset->nDataSets; j++) {
          if (strcmp("paramd", vertexset->dataSets[j]->name) != 0) continue;
          ds    = (capsDataSet *) vertexset->dataSets[j]->blind;
          param = ds->data;
          break;
        }
      } else {
        for (j = 0; j < vertexset->nDataSets; j++) {
          if (strcmp("param", vertexset->dataSets[j]->name) != 0) continue;
          ds    = (capsDataSet *) vertexset->dataSets[j]->blind;
          param = ds->data;
          break;
        }
      }
      for (j = 0; j < fvs->nDataSets; j++) {
        if (strcmp("param", fvs->dataSets[j]->name) != 0) continue;
        ds     = (capsDataSet *) fvs->dataSets[j]->blind;
        params = ds->data;
        break;
      }
    }

    if (params == NULL) {
      snprintf(temp, PATH_MAX, "Bound %s -- %s without source params (caps_getData)!",
               bobject->name, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_SOURCEERR;
    }
    if (param == NULL) {
      snprintf(temp, PATH_MAX, "Bound %s -- %s without dst params (caps_getData)!",
               bobject->name, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_SOURCEERR;
    }
    if (dataset->data != NULL) EG_free(dataset->data);
    dataset->data = NULL;
    values = (double *) EG_alloc(dataset->rank*npts*sizeof(double));
    if (values == NULL) {
      snprintf(temp, PATH_MAX, "Malloc on %dx%d  Dataset = %s (caps_getData)!",
               npts, dataset->rank, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return EGADS_MALLOC;
    }
    for (i = 0; i < j; i++) values[i] = 0.0;

    stat = caps_getDataX(foundset, &i, &j, &data, &units, nErr, errors);
    if (stat != CAPS_SUCCESS) {
      snprintf(temp, PATH_MAX, "Could not get source %s for FieldIn %s (caps_getData)!",
               foundset->name, dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      EG_free(values);
      return stat;
    }

    /* compute */
    if (otherset->data == NULL) {
      snprintf(temp, PATH_MAX, "Source for %s is NULL (caps_getData)!\n",
               dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      EG_free(values);
      return CAPS_NULLVALUE;
    }
    if (dataset->linkMethod == Interpolate) {
      for (i = 0; i < npts; i++) {
        stat = aim_LocateElIndex(problem->aimFPTR, index, founddsr, params,
                                 &param[bound->dim*i], &bIndex, &eIndex, st);
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Warning: %d/%d aim_LocateElement = %d for %s!\n",
                 i, npts, stat, dobject->name);
          continue;
        }
        stat = aim_InterpolIndex(problem->aimFPTR, index, founddsr,
                                 dobject->name, bIndex, eIndex, st, dataset->rank,
                                 otherset->data, &values[dataset->rank*i]);
        if (stat != CAPS_SUCCESS)
          printf(" CAPS Warning: %d/%d aim_Interpolation = %d for %s!\n",
                 i, npts, stat, dobject->name);
      }
    } else {
      if ((aobject == NULL) || (analysis == NULL)) {
        EG_free(values);
        return CAPS_BADMETHOD;
      }
      fit.name     = dobject->name;
      /*@-immediatetrans@*/
      fit.aimFPTR  = &problem->aimFPTR;
      /*@+immediatetrans@*/
      fit.npts     = npts;
      fit.nrank    = dataset->rank;
      fit.sindx    = index;
      fit.tindx    = aim_Index(problem->aimFPTR, analysis->loadName);
      fit.afact    = 1.e6;
      fit.area_src = 0.0;
      fit.area_tgt = 0.0;
      fit.src      = founddsr;
      fit.prms_src = params;
      fit.data_src = otherset->data;
      fit.tgt      = discr;
      fit.prms_tgt = param;
      fit.data_tgt = values;
      fit.nmat     = 0;
      fit.mat      = NULL;
      stat = caps_Match(&fit, bound->dim);
      if (stat != CAPS_SUCCESS) {
        EG_free(values);
        return stat;
      }
      stat = caps_Conserve(&fit, bobject->name, bound->dim);
      if (stat != CAPS_SUCCESS) {
        if (fit.mat != NULL) EG_free(fit.mat);
        EG_free(values);
        return stat;
      }

      if (fit.mat != NULL) EG_free(fit.mat);
    }

    dataset->data = values;
    dataset->npts = npts;
    if (units != NULL) {
      EG_free(dataset->units);
      dataset->units = EG_strdup(units);
      if (dataset->units == NULL) {
        printf(" CAPS Error: Failed to allocate units!\n");
        return EGADS_MALLOC;
      }
    }
    caps_freeOwner(&dobject->last);
    problem->sNum     += 1;
    dobject->last.sNum = problem->sNum;
    stat = caps_addHistory(dobject, problem);
    if (stat != CAPS_SUCCESS) {
      printf(" CAPS Warning: caps_addHistory = %d (caps_getData)\n", stat);
      stat = caps_writeSerialNum(problem);
      if (stat != CAPS_SUCCESS)
        printf(" CAPS Error: caps_writeSerialNum = %d (caps_getData)\n", stat);
    }
    stat = caps_writeDataSet(dobject);
    if (stat != CAPS_SUCCESS)
      printf(" CAPS Error: caps_writeDataSet = %d (caps_getData)\n", stat);
  }

  return CAPS_SUCCESS;
}


static void
caps_getDataDB(capsDataSet *dataset, int *npts, int *rank, double **data,
               char **units)
{
  *rank  = dataset->rank;
  *npts  = dataset->npts;
  *data  = dataset->data;
  *units = dataset->units;
}


int
caps_getDataX(capsObject *dobject, int *npts, int *rank, double **data,
              char **units, int *nErr, capsErrs **errors)
{
  int           stat;
  char          temp[PATH_MAX];
  capsObject    *vobject, *aobject, *bobject, *foundset, *foundanl;
  capsBound     *bound;
  capsDataSet   *dataset;
  capsVertexSet *vertexset, *fvs;

  *npts   = *rank = 0;
  *data   = NULL;
  *units  = NULL;
  *nErr   = 0;
  *errors = NULL;
  if (dobject              == NULL)      return CAPS_NULLOBJ;
  if (dobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (dobject->type        != DATASET)   return CAPS_BADTYPE;
  if (dobject->blind       == NULL)      return CAPS_NULLBLIND;
  dataset = (capsDataSet *) dobject->blind;
  vobject = dobject->parent;
  if (vobject              == NULL)      return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET) return CAPS_BADTYPE;
  if (vobject->blind       == NULL)      return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vobject->blind;
  aobject   = vertexset->analysis;
  bobject   = vobject->parent;
  if (bobject              == NULL)      return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)     return CAPS_BADTYPE;
  if (bobject->blind       == NULL)      return CAPS_NULLBLIND;
  bound   = (capsBound *) bobject->blind;
  if (bound->state         == Open)      return CAPS_STATEERR;

  /*
   * Sensitivity/BuiltIn - filled in fillVertexSets
   * User                - from explicit calls to setData
   */
  if ((dataset->ftype == GeomSens) || (dataset->ftype == TessSens) ||
      (dataset->ftype == BuiltIn)) {

    /* make sure that VertexSets and Sensitvities are up-to-date */
    stat = caps_buildBound(bobject, nErr, errors);
    if ((stat != CAPS_CLEAN) && (stat != CAPS_SUCCESS)) return stat;

  } else if (dataset->ftype == FieldOut) {

    /* fill in FieldOut DataSet from the AIM */
    stat = caps_fillFieldOut(dobject, nErr, errors);
    if (stat != CAPS_SUCCESS) return stat;

  } else if (dataset->ftype == FieldIn) {

    /* check link to other DataSet */
    if (dataset->link == NULL) {
      snprintf(temp, PATH_MAX, "FieldIn DataSet '%s' without a link!",
               dobject->name);
      caps_makeSimpleErr(dobject, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_NULLOBJ;
    }

    foundset = dataset->link;
    if (foundset->magicnumber         != CAPSMAGIC) return CAPS_BADOBJECT;
    if (foundset->type                != DATASET)   return CAPS_BADTYPE;
    if (foundset->blind               == NULL)      return CAPS_NULLBLIND;
    if (foundset->parent              == NULL)      return CAPS_NULLOBJ;
    if (foundset->parent->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (foundset->parent->type        != VERTEXSET) return CAPS_BADTYPE;
    if (foundset->parent->blind       == NULL)      return CAPS_NULLBLIND;
    fvs = (capsVertexSet *) foundset->parent->blind;

    foundanl = fvs->analysis;

    if ((aobject->last.sNum  == 0) &&
        (foundanl->last.sNum == 0) &&
        (dataset->startup    != NULL)) {
      /* bypass everything because we are in a startup situation */
      *rank  = dataset->rank;
      *npts  = 1;
      *data  = dataset->startup;
      *units = dataset->units;
      return CAPS_SUCCESS;
    }

    /* fill in FieldIn DataSet from a linked FieldOut or User DataSet */
    stat = caps_fillFieldIn(dobject, nErr, errors);
    if (stat != CAPS_SUCCESS) return stat;

  }

  *rank  = dataset->rank;
  *npts  = dataset->npts;
  *data  = dataset->data;
  *units = dataset->units;

  return CAPS_SUCCESS;
}


int
caps_getData(capsObject *dobject, int *npts, int *rank, double **data,
             char **units, int *nErr, capsErrs **errors)
{
  int           ret, stat;
  CAPSLONG      sNum;
  capsObject    *vobject, *bobject, *pobject, *aobject;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsProblem   *problem;
  capsJrnl      args[6];

  if (nErr                 == NULL)      return CAPS_NULLVALUE;
  if (errors               == NULL)      return CAPS_NULLVALUE;
  *npts   = *rank = 0;
  *data   = NULL;
  *units  = NULL;
  *nErr   = 0;
  *errors = NULL;
  if (dobject              == NULL)      return CAPS_NULLOBJ;
  if (dobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (dobject->type        != DATASET)   return CAPS_BADTYPE;
  if (dobject->blind       == NULL)      return CAPS_NULLBLIND;
  dataset = (capsDataSet *) dobject->blind;
  vobject = dobject->parent;
  if (vobject              == NULL)      return CAPS_NULLOBJ;
  if (vobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (vobject->type        != VERTEXSET) return CAPS_BADTYPE;
  if (vobject->blind       == NULL)      return CAPS_NULLBLIND;
  vertexset = (capsVertexSet *) vobject->blind;
  aobject   = vertexset->analysis;
  bobject   = vobject->parent;
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

  if (problem->dbFlag == 1) {
    caps_getDataDB(dataset, npts, rank, data, units);
    return CAPS_SUCCESS;
  }

  if (aobject != NULL) {
    if (aobject->magicnumber  != CAPSMAGIC)          return CAPS_BADOBJECT;
    if (aobject->type         != ANALYSIS)           return CAPS_BADTYPE;
    if (aobject->blind        == NULL)               return CAPS_NULLBLIND;
    analysis = (capsAnalysis *) aobject->blind;
    if ((dataset->ftype == GeomSens) || (dataset->ftype == TessSens) ||
        (dataset->ftype == FieldOut)) {
      if (problem->geometry.sNum > analysis->pre.sNum) return CAPS_MISMATCH;
      if (analysis->pre.sNum     > aobject->last.sNum) return CAPS_MISMATCH;
    }
  }

  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jPointer;
  args[3].type = jString;
  args[4].type = jInteger;
  args[5].type = jErr;
  stat         = caps_jrnlRead(CAPS_GETDATA, problem, dobject, 6, args, &sNum,
                               &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *npts   = args[0].members.integer;
    *rank   = args[1].members.integer;
    *data   = args[2].members.pointer;
    *units  = args[3].members.string;
    *nErr   = args[4].members.integer;
    *errors = args[5].members.errs;
    return ret;
  }

  sNum = problem->sNum;
  if (ret == CAPS_SUCCESS) {
    ret = caps_getDataX(dobject, npts, rank, data, units, nErr, errors);
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
  }
  args[0].members.integer = *npts;
  args[1].members.integer = *rank;
  args[2].length          = *rank*args[0].members.integer*sizeof(double);
  args[2].members.pointer = *data;
  args[3].members.string  = *units;
  args[4].members.integer = *nErr;
  args[5].members.errs    = *errors;
  caps_jrnlWrite(CAPS_GETDATA, problem, dobject, ret, 6, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_getDataSets(capsObject *bobject, const char *dname, int *nobj,
                 capsObject ***dobjs)
{
  int           i, j, n, status, ret;
  CAPSLONG      sNum;
  capsObject    *pobject;
  capsProblem   *problem;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsObject    **objs;
  capsJrnl      args[1];

  *nobj  = 0;
  *dobjs = NULL;
  if (dname                == NULL)         return CAPS_NULLNAME;
  if (bobject              == NULL)         return CAPS_NULLOBJ;
  if (bobject->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (bobject->type        != BOUND)        return CAPS_BADTYPE;
  if (bobject->blind       == NULL)         return CAPS_NULLBLIND;
  status = caps_findProblem(bobject, CAPS_GETDATASETS, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  bound   = (capsBound *)   bobject->blind;

  args[0].type = jObjs;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_GETDATASETS, problem, bobject, 1, args,
                               &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *nobj = args[0].num;
        if (*nobj != 0) {
          objs = (capsObject **) EG_alloc(*nobj*sizeof(capsObject *));
          if (objs == NULL) return EGADS_MALLOC;
          for (i = 0; i < *nobj; i++) objs[i] = args[0].members.objs[i];
          *dobjs = objs;
        }
      }
      return ret;
    }
  }

  status = CAPS_SUCCESS;
  for (n = i = 0; i < bound->nVertexSet; i++) {
    if (bound->vertexSet[i]              == NULL) {
      status = CAPS_NULLOBJ;
      goto gdone;
    }
    if (bound->vertexSet[i]->magicnumber != CAPSMAGIC) {
      status = CAPS_BADOBJECT;
      goto gdone;
    }
    if (bound->vertexSet[i]->type        != VERTEXSET) {
      status = CAPS_BADTYPE;
      goto gdone;
    }
    if (bound->vertexSet[i]->blind       == NULL) {
      status = CAPS_NULLBLIND;
      goto gdone;
    }
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    for (j = 0; j < vertexset->nDataSets; j++)
      if (strcmp(dname, vertexset->dataSets[j]->name) == 0) n++;
  }

  if (n == 0) goto gdone;
  objs = (capsObject **) EG_alloc(n*sizeof(capsObject *));
  if (objs == NULL) {
    status = EGADS_MALLOC;
    goto gdone;
  }

  for (n = i = 0; i < bound->nVertexSet; i++) {
    vertexset = (capsVertexSet *) bound->vertexSet[i]->blind;
    for (j = 0; j < vertexset->nDataSets; j++)
      if (strcmp(dname, vertexset->dataSets[j]->name) == 0) {
        objs[n] = vertexset->dataSets[j];
        n++;
      }
  }

  *nobj  = n;
  *dobjs = objs;

gdone:
  if (problem->dbFlag == 1) return status;
  args[0].num          = *nobj;
  args[0].members.objs = *dobjs;
  caps_jrnlWrite(CAPS_GETDATASETS, problem, bobject, status, 1, args,
                 problem->sNum, problem->sNum);
  return status;
}


int
caps_snDataSets(const capsObject *aobject, int flag, CAPSLONG *sn)
/* flag = 0 - return lowest sn of source
   flag = 1 - also check if the linked analysis is dirty */
{
  int           i, j, k, status, major, minor, nField, exec, dirty;
  int           *ranks, *fInOut;
  char          *intents, *apath, *unitSys, **fnames;
  capsObject    *pobject, *linkanl;
  capsProblem   *problem;
  capsBound     *bound;
  capsVertexSet *vs, *vso;
  capsDataSet   *ds;

  *sn = 0;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (aobject->type        != ANALYSIS)  return CAPS_BADTYPE;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;
  pobject = (capsObject *)  aobject->parent;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  *sn     = -1;

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
        if (ds->link                       == NULL)      continue;

        if (ds->link->magicnumber != CAPSMAGIC)         return CAPS_BADOBJECT;
        if (ds->link->type        != DATASET)           return CAPS_BADTYPE;
        if (ds->link->blind       == NULL)              return CAPS_NULLBLIND;
        if (ds->link->parent              == NULL)      return CAPS_NULLOBJ;
        if (ds->link->parent->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
        if (ds->link->parent->type        != VERTEXSET) return CAPS_BADTYPE;
        if (ds->link->parent->blind       == NULL)      return CAPS_NULLBLIND;
        vso = (capsVertexSet *) ds->link->parent->blind;
        linkanl = vso->analysis;
        if (linkanl == NULL)                           return CAPS_NULLOBJ;

        if (*sn == -1) {
          *sn = linkanl->last.sNum;
        } else {
          if (*sn > linkanl->last.sNum)
            *sn = linkanl->last.sNum;
        }

        if (flag == 1) {
          /* check to see if analysis is dirty */
          status = caps_analysisInfX(linkanl, &apath, &unitSys, &major, &minor,
                                     &intents, &nField, &fnames, &ranks,
                                     &fInOut, &exec, &dirty);
          if (status != CAPS_SUCCESS) return status;

          if (((exec == 2) && (dirty > 0) && (dirty < 5)) ||
              (linkanl->last.sNum > ds->link->last.sNum)  ||
              (ds->link->last.sNum == 0)) {
            *sn = aobject->last.sNum+2;
          }
        }
      }
    }
  }

  if (*sn == -1) {
    *sn = 0;
    return CAPS_NOTFOUND;
  }
  return CAPS_SUCCESS;
}
