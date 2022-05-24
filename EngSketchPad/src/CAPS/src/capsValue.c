/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Value Object Functions
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
#define getcwd   _getcwd
#define PATH_MAX _MAX_PATH
#define snprintf _snprintf
#else
#include <unistd.h>
#include <limits.h>
#endif

#include "udunits.h"

#include "capsBase.h"
#include "capsAIM.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"


/*@-incondefs@*/
extern void ut_free(/*@only@*/ ut_unit* const unit);
/*@+incondefs@*/

extern /*@null@*/ /*@only@*/
       char *EG_strdup(/*@null@*/ const char *str);

extern int  caps_rename(const char *src, const char *dst);
extern int  caps_writeValueObj(capsProblem *problem, capsObject *obj);
extern void caps_jrnlWrite(int funID, capsProblem *problem, capsObject *obj,
                           int status, int nargs, capsJrnl *args, CAPSLONG sNm0,
                           CAPSLONG sNum);
extern int  caps_jrnlEnd(capsProblem *problem);
extern int  caps_jrnlRead(int funID, capsProblem *problem, capsObject *obj,
                          int nargs, capsJrnl *args, CAPSLONG *sNum, int *stat);

extern int  caps_integrateData(const capsObject *vobject, enum capstMethod method,
                               int *rank, double **data, char **units);
extern int  caps_getDataX(const capsObject *vobject, int *npts, int *rank,
                          const double **data, const char **units);
extern void caps_geomOutSensit(capsProblem *problem, int ipmtr, int irow,
                               int icol);
extern int  caps_build(capsObject *pobject, int *nErr, capsErrs **errors);
extern void caps_concatErrs(/*@null@*/ capsErrs *errs, capsErrs **errors);
extern void caps_getAIMerrs(capsAnalysis *analy, int *nErr, capsErrs **errors);
extern int  caps_analysisInfX(const capsObj aobject, char **apath, char **unSys,
                              int *major, int *minor, char **intents,
                              int *nField, char ***fnames, int **ranks,
                              int **fInOut, int *execution, int *status);
extern int  caps_postAnalysiX(capsObject *aobject, int *nErr, capsErrs **errors,
                              int flag);
extern int  caps_unitConvertible(const char *unit1, const char *unit2);
extern int  caps_convert(int count,
                         /*@null@*/ const char  *inUnit, double  *inVal,
                         /*@null@*/ const char *outUnit, double *outVal);
extern int  caps_unitParse(/*@null@*/ const char *unit);
extern int  caps_execX(capsObject *aobject, int *nErr, capsErrs **errors);
extern int  caps_circularAutoExecs(capsObject *asrc,
                                   /*@null@*/ capsObject *aobject);
extern int  caps_transferLinks(capsAnalysis *analysis, int *nErr,
                               capsErrs **errors);


int
caps_checkValueObj(capsObject *object)
{
  const char *name = NULL;
  capsValue  *value;

  if (object              == NULL)        return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC)   return CAPS_BADOBJECT;
  if (object->type        != VALUE)       return CAPS_BADTYPE;
  if (object->blind       == NULL)        return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;

  /* check for valid characters in the name */
  name = object->name;
  if(name == NULL) {
    printf(" CAPS Error: Value Object name is NULL\n");
    return CAPS_NULLNAME;
  }
  /*
    len2 = strlen(name);
    for (j = 0; j < len2; j++) {
      if ( !(((name[j] >= 'a') && (name[j] <= 'z')) ||
             ((name[j] >= 'A') && (name[j] <= 'Z')) ||
              (name[j] == '_') || (name[j] == ':')) )
        return CAPS_BADNAME;
    }
   */

  /* check units */
  if ((value->type == Pointer) || (value->type == PointerMesh)) {
    /* Don't look at AIMptr units -- these are programmer defined */
    if ((value->type == Pointer) && (value->units == NULL)) {
      printf(" CAPS Error: Value '%s' Pointer units is NULL!\n", object->name);
      return CAPS_UNITERR;
    }
  } else {
    if (value->units != NULL) {
      /* Only double and Strings can have units */
      if (value->type == Boolean) {
        printf(" CAPS Error: Value '%s' Integer cannot have units '%s'!\n",
               object->name, value->units);
        return CAPS_UNITERR;
      } else if (value->type == Integer) {
        printf(" CAPS Error: Value '%s' Boolean cannot have units '%s'!\n",
               object->name, value->units);
        return CAPS_UNITERR;
      } else if (value->type == Tuple) {
        printf(" CAPS Error: Value '%s' Tuple cannot have units '%s'!\n",
               object->name, value->units);
        return CAPS_UNITERR;
      }

      if ((value->type == Double) ||
          (value->type == DoubleDeriv)) {
        if (caps_unitParse(value->units) != CAPS_SUCCESS) {
          printf(" CAPS Error: Value '%s' invalid units '%s'!\n",
                 object->name, value->units);
          return CAPS_UNITERR;
        }
      }
    }
  }

  /* set the length */
  value->length = value->ncol*value->nrow;
  if (value->length <= 0 && value->nullVal != IsNull) {
    printf(" CAPS Error: Value '%s' length <= 0 (nrow = %d, ncol = %d) and not IsNull!\n",
           object->name, value->nrow, value->ncol);
    return CAPS_SHAPEERR;
  }

  /* check nullval */
  if (value->type == Integer) {

    if (value->dim > Scalar && value->length > 1) {
      if (value->nullVal == NotAllowed && value->vals.integers == NULL) {
        printf(" CAPS Error: Value '%s' integers is NULL (NotAllowed)!\n",
               object->name);
        return CAPS_NULLVALUE;
      }
    }

  } else if ((value->type == Double) || (value->type == DoubleDeriv)) {

    if (value->dim > Scalar && value->length > 1) {
      if (value->nullVal == NotAllowed && value->vals.reals == NULL) {
        printf(" CAPS Error: Value '%s' reals is NULL (NotAllowed)!\n",
               object->name);
        return CAPS_NULLVALUE;
      }
    }

  } else if (value->type == String) {

    if (value->nullVal == NotAllowed && value->vals.string == NULL) {
      printf(" CAPS Error: Value '%s' string is NULL (NotAllowed)!\n",
             object->name);
      return CAPS_NULLVALUE;
    }
    if (value->nullVal != NotAllowed)
      value->nullVal = (value->vals.string == NULL) ? IsNull : NotNull;

  } else if (value->type == Tuple) {

    if (value->nullVal == NotAllowed && value->vals.tuple == NULL) {
      printf(" CAPS Error: Value '%s' tuple is NULL (NotAllowed)!\n",
             object->name);
      return CAPS_NULLVALUE;
    }
    if (value->nullVal != NotAllowed)
      value->nullVal = (value->vals.tuple == NULL) ? IsNull : NotNull;

    if (value->dim != Vector) {
      printf(" CAPS Error: Value '%s' Tuple dim must be Vector!\n",
             object->name);
      return CAPS_SHAPEERR;
    }

  } else if ((value->type == Pointer) || (value->type == PointerMesh)) {

    if (value->nullVal == NotAllowed && value->vals.AIMptr == NULL) {
      printf(" CAPS Error: Value '%s' AIMptr is NULL (NotAllowed)!\n",
             object->name);
      return CAPS_NULLVALUE;
    }
    if (value->nullVal != NotAllowed)
      value->nullVal = (value->vals.AIMptr == NULL) ? IsNull : NotNull;

  }

  /* look at shapes */
  if (value->dim == Scalar) {
    if (value->length > 1) {
      printf(" CAPS Error: Value '%s' dim = Scalar and length > 1 (length = %d)!\n",
             object->name, value->length);
      return CAPS_SHAPEERR;
    }
  } else if (value->dim == Vector) {
    if ((value->ncol != 1) && (value->nrow != 1)) {
      printf(" CAPS Error: Value '%s' dim = Vector nrow = %d ncol = %d!\n",
             object->name, value->nrow, value->ncol);
      return CAPS_SHAPEERR;
    }
  } else if (value->dim != Array2D) {
    printf(" CAPS Error: Value '%s' Invalid dim = %d!\n",
           object->name, value->dim);
    return CAPS_BADINDEX;
  }

  return CAPS_SUCCESS;
}


static int
caps_updateAnalysisOut(capsObject *vobject, const int funID,
                       int *nErr, capsErrs **errors)
{
  int          i, in, status, major, minor, nField, exec, dirty;
  int          *ranks, *fInOut;
  char         *intents, *apath, *unitSys, **fnames, temp[PATH_MAX];
  capsValue    *value, *valu0, *valIn=NULL;
  capsObject   *pobject, *aobject;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsErrs     *errs = NULL;

  *nErr   = 0;
  *errors = NULL;

  if (vobject->type    != VALUE)         return CAPS_BADTYPE;
  if (vobject->subtype != ANALYSISOUT)   return CAPS_BADTYPE;
  value = (capsValue *) vobject->blind;

  status  = caps_findProblem(vobject, funID, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  /* do we need to update our value? */
  aobject = vobject->parent;
  if (aobject              == NULL)      return CAPS_NULLOBJ;
  if (aobject->blind       == NULL)      return CAPS_NULLBLIND;
  analysis = (capsAnalysis *) aobject->blind;
  if (aobject->parent      == NULL)      return CAPS_NULLOBJ;

  /* check to see if analysis is CLEAN, and auto-execute if possible */
  status = caps_analysisInfX(aobject, &apath, &unitSys, &major, &minor,
                             &intents, &nField, &fnames, &ranks, &fInOut,
                             &exec, &dirty);
  if (status != CAPS_SUCCESS) return status;
  if (dirty > 0) {
    /* auto execute if available */
    if ((exec == 2) && (dirty < 5)) {
      status = caps_execX(aobject, nErr, errors);
      if (status != CAPS_SUCCESS) return status;
    } else {
      return CAPS_DIRTY;
    }
  } else if (analysis->reload == 1) {

    /* fill in any values that have links */
    status = caps_transferLinks(analysis, nErr, errors);
    if (status != CAPS_SUCCESS) return status;

    if (analysis->nAnalysisIn > 0)
      valIn  = (capsValue *) analysis->analysisIn[0]->blind;
/*  analysis->info.inPost = 1;  */
    status   = aim_PostAnalysis(problem->aimFPTR, analysis->loadName,
                                analysis->instStore, &analysis->info, 1,
                                valIn);
    analysis->info.inPost = 0;

    if (*nErr != 0) {
      errs    = *errors;
      *nErr   = 0;
      *errors = NULL;
    }
    caps_getAIMerrs(analysis, nErr, errors);
    caps_concatErrs(errs, errors);
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: PostAnalysis %s (%d) ret = %d (%s -> caps_updateAnalysisOut)!\n",
             analysis->loadName, analysis->info.instance, status, caps_funID[funID]);
      return CAPS_BADINIT;
    }
    analysis->reload = 0;
  }

  valu0 = (capsValue *) analysis->analysisOut[0]->blind;
  in    = value - valu0;
  if (vobject->last.sNum < aobject->last.sNum) {
    if ((value->type == Boolean) || (value->type == Integer)) {
      if (value->length > 1) {
        EG_free(value->vals.integers);
      }
      value->vals.integers = NULL;
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      if (value->length > 1) {
        EG_free(value->vals.reals);
      }
      value->vals.reals = NULL;
      for (i = 0; i < value->nderiv; i++) {
        EG_free(value->derivs[i].name);
        EG_free(value->derivs[i].deriv);
      }
      EG_free(value->derivs);
      value->derivs = NULL;
      value->nderiv = 0;
    } else if (value->type == String) {
      EG_free(value->vals.string);
      value->vals.string = NULL;
    } else if (value->type == Tuple) {
      caps_freeTuple(value->length, value->vals.tuple);
      value->vals.tuple = NULL;
    } else if ((value->type != Pointer) && (value->type != PointerMesh)) {
      return CAPS_BADTYPE;
    }
    caps_freeOwner(&vobject->last);
    vobject->last.sNum = 0;
    status = aim_CalcOutput(problem->aimFPTR, analysis->loadName,
                            analysis->instStore, &analysis->info, in+1, value);
    caps_getAIMerrs(analysis, nErr, errors);
    if (status != CAPS_SUCCESS) return status;
    status = caps_checkValueObj(analysis->analysisOut[in]);
    if (status != CAPS_SUCCESS) {
      snprintf(temp, PATH_MAX, "'%s' CalcOutput Check Value %d (%s)",
               analysis->loadName, status, caps_funID[funID]);
      caps_makeSimpleErr(analysis->analysisOut[in], CERROR, temp,
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return status;
    }
    vobject->last.sNum = aobject->last.sNum;
    status = caps_addHistory(vobject, problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_addHistory = %d (%s)\n",
             status, caps_funID[funID]);
    status = caps_writeValueObj(problem, analysis->analysisOut[in]);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_writeValueObj = %d (%s)\n",
             status, caps_funID[funID]);
  }

  return CAPS_SUCCESS;
}


static int
caps_getValuX(capsObject *object, enum capsvType *type, int *nrow, int *ncol,
              const void **data, const int **partial, const char **units,
              int *nErr, capsErrs **errors)
{
  int         status;
  capsValue   *value;
  capsObject  *pobject, *source, *last;
  capsProblem *problem;

  *nErr   = 0;
  *errors = NULL;
  status  = caps_findProblem(object, 9999, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  source = object;
  do {
    if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
    if (source->type        != VALUE)     return CAPS_BADTYPE;
    if (source->blind       == NULL)      return CAPS_NULLBLIND;
    value = (capsValue *) source->blind;
    if (value->link == object)            return CAPS_CIRCULARLINK;
    last   = source;
    source = value->link;
  } while (value->link != NULL);

  /* do we need to update our value? */
  if (problem->dbFlag == 0)
    if (last->subtype == ANALYSISOUT) {
      status = caps_updateAnalysisOut(last, CAPS_GETVALUE, nErr, errors);
      if (status != CAPS_SUCCESS) return status;
    } else if (last->subtype == GEOMETRYOUT) {
      /* make sure geometry is up-to-date */
      status = caps_build(pobject, nErr, errors);
      if ((status != CAPS_SUCCESS) && (status != CAPS_CLEAN)) return status;
    }

  *type    = value->type;
  *nrow    = value->nrow;
  *ncol    = value->ncol;
  *partial = value->partial;
  *units   = value->units;
  if (value->nullVal == IsNull) *nrow = *ncol = 0;

  /* always return nrow > 1 for Vector */
  if (value->dim == Vector && value->ncol > 1) {
    *nrow = value->ncol;
    *ncol = value->nrow;
  }

  if (data != NULL) {
    *data = NULL;
    if (value->nullVal != IsNull) {
      if ((value->type == Boolean) || (value->type == Integer)) {
        if (value->length == 1) {
          *data = &value->vals.integer;
        } else {
          *data =  value->vals.integers;
        }
      } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
        if (value->length == 1) {
          *data = &value->vals.real;
        } else {
          *data =  value->vals.reals;
        }
      } else if (value->type == String) {
        *data = value->vals.string;
      } else if (value->type == Tuple) {
        *data = value->vals.tuple;
      } else {
        *data = value->vals.AIMptr;
      }
    }
  }

  return CAPS_SUCCESS;
}


int
caps_getValue(capsObject *object, enum capsvType *type, int *nrow, int *ncol,
              const void **data, const int **partial, const char **units,
              int *nErr, capsErrs **errors)
{
  int            i, stat, ret, vlen;
  const char     *string;
  CAPSLONG       sNum;
  capsObject     *pobject;
  capsProblem    *problem;
  capsValue      *value;
  capsJrnl       args[8];
  enum capsvType itype;
  size_t         len;

  if (nErr                == NULL)      return CAPS_NULLVALUE;
  if (errors              == NULL)      return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (nrow                == NULL)      return CAPS_NULLVALUE;
  if (ncol                == NULL)      return CAPS_NULLVALUE;
  if (partial             == NULL)      return CAPS_NULLVALUE;
  if (units               == NULL)      return CAPS_NULLVALUE;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->type        != VALUE)     return CAPS_BADTYPE;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;
  stat  = caps_findProblem(object, CAPS_GETVALUE, &pobject);
  if (stat != CAPS_SUCCESS)             return stat;
  problem = (capsProblem *) pobject->blind;

  /* just return the values if in debug mode */
  if (problem->dbFlag == 1) {
    ret   = caps_getValuX(object, type, nrow, ncol, data, partial, units,
                          nErr, errors);
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
    return ret;
  }

  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jInteger;
  args[3].type = jPointer;
  args[4].type = jPointer;
  args[5].type = jString;
  args[6].type = jInteger;
  args[7].type = jErr;
  if (value->type == Tuple) args[3].type = jTuple;
  stat         = caps_jrnlRead(CAPS_GETVALUE, problem, object, 8, args, &sNum,
                               &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *type    = args[0].members.integer;
    *nrow    = args[1].members.integer;
    *ncol    = args[2].members.integer;
    if (value->type == Tuple) {
      *data  = args[3].members.tuple;
    } else {
      *data  = args[3].members.pointer;
    }
    *partial = args[4].members.pointer;
    *units   = args[5].members.string;
    *nErr    = args[6].members.integer;
    *errors  = args[7].members.errs;
    return ret;
  }

  itype    = Pointer;
  *nrow    = *ncol = 0;
  *data    = NULL;
  *partial = NULL;
  *units   = NULL;
  if (ret == CAPS_SUCCESS) {
    ret   = caps_getValuX(object, &itype, nrow, ncol, data, partial, units,
                          nErr, errors);
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
  }
  args[0].members.integer = itype;
  args[1].members.integer = *nrow;
  args[2].members.integer = *ncol;
  if (itype == Tuple) args[3].members.tuple   = (void *) *data;
/*@-kepttrans@*/
  if (itype != Tuple) args[3].members.pointer = (void *) *data;
/*@+kepttrans@*/
  args[4].members.pointer = (void *) *partial;
  args[5].members.string  = (void *) *units;
  args[6].members.integer = *nErr;
  args[7].members.errs    = *errors;
  args[3].length          = 0;
  args[3].num             = 0;
  args[4].length          = 0;
  vlen                    = *nrow;
  vlen                   *= *ncol;
  if (*data != NULL) {
    if (itype == Tuple) {
      args[3].num = vlen;
    } else {
      len = sizeof(char *);
      if ((itype == Integer) || (itype == Boolean)) {
        len = vlen*sizeof(int);
      } else if ((itype == Double) || (itype == DoubleDeriv)) {
        len = vlen*sizeof(double);
      } else if (itype == String) {
        string = (char *) *data;
        if (vlen == 1) {
          len = (strlen(string)+1)*sizeof(char);
        } else {
          len  = 0;
          for (i = 0; i < vlen; i++)
            len += strlen(string + len) + 1;
          len *= sizeof(char);
        }
      }
      args[3].length = len;
    }
  }
  if (*partial != NULL) {
    len  = *nrow;
    len *= *ncol*sizeof(int);
    args[4].length = len;
  }
  *type = itype;
  caps_jrnlWrite(CAPS_GETVALUE, problem, object, ret, 8, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_makeValueX(capsObject *pobject, const char *vname, enum capssType stype,
                enum capsvType vtype, int nrow, int ncol, const void *data,
                /*@null@*/ int *partial, /*@null@*/ const char *units,
                capsObject **vobj)
{
  int         status, i, vlen;
  char        filename[PATH_MAX], temp[PATH_MAX];
  capsObject  *object, **tmp;
  capsProblem *problem;
  capsValue   *value;
  FILE        *fp;

  vlen    = ncol*nrow;
  problem = (capsProblem *) pobject->blind;

  status  = caps_makeVal(vtype, vlen, data, &value);
  if (status != CAPS_SUCCESS) return status;
  value->nrow = nrow;
  value->ncol = ncol;
  value->dim  = Scalar;
  if ((nrow > 1) || (ncol > 1)) value->dim = Vector;
  if ((nrow > 1) && (ncol > 1)) value->dim = Array2D;

  /* set partial if applicable */
  if (partial != NULL) {
    value->partial = (int*) EG_alloc(vlen*sizeof(int));
    if (value->partial == NULL) {
      if (value->length != 1) {
        if (value->type == Integer) {
          EG_free(value->vals.integers);
        } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
          EG_free(value->vals.reals);
        } else if (value->type == String) {
          EG_free(value->vals.string);
        } else if (value->type == Tuple) {
          caps_freeTuple(value->length, value->vals.tuple);
        }
      }
      EG_free(value);
      return EGADS_MALLOC;
    }
    for (i = 0; i < vlen; i++) value->partial[i] = partial[i];
    value->nullVal = IsPartial;
  }

  /* make the object */
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    if (value->length != 1) {
      if (value->type == Integer) {
        EG_free(value->vals.integers);
      } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
        EG_free(value->vals.reals);
      } else if (value->type == String) {
        EG_free(value->vals.string);
      } else if (value->type == Tuple) {
        caps_freeTuple(value->length, value->vals.tuple);
      }
    }
    EG_free(value->partial);
    EG_free(value);
    return status;
  }
  object->parent = pobject;

  /* save away the object in the problem object */
  if (stype == PARAMETER) {
    value->index = problem->nParam + 1;
    if (problem->params == NULL) {
      problem->params = (capsObject **) EG_alloc(sizeof(capsObject *));
      if (problem->params == NULL) {
        if (value->length != 1) {
          if (value->type == Integer) {
            EG_free(value->vals.integers);
          } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
            EG_free(value->vals.reals);
          } else if (value->type == String) {
            EG_free(value->vals.string);
          } else if (value->type == Tuple) {
            caps_freeTuple(value->length, value->vals.tuple);
          }
        }
        EG_free(value->partial);
        EG_free(value);
        EG_free(object);
        return EGADS_MALLOC;
      }
    } else {
      tmp = (capsObject **) EG_reall( problem->params,
                                     (problem->nParam+1)*sizeof(capsObject *));
      if (tmp == NULL) {
        if (value->length != 1) {
          if (value->type == Integer) {
            EG_free(value->vals.integers);
          } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
            EG_free(value->vals.reals);
          } else if (value->type == String) {
            EG_free(value->vals.string);
          } else if (value->type == Tuple) {
            caps_freeTuple(value->length, value->vals.tuple);
          }
        }
        EG_free(value->partial);
        EG_free(value);
        EG_free(object);
        return EGADS_MALLOC;
      }
      problem->params = tmp;
    }
    problem->params[problem->nParam] = object;
    problem->nParam += 1;
  } else {
    value->index = problem->nUser + 1;
    if (problem->users == NULL) {
      problem->users = (capsObject **) EG_alloc(sizeof(capsObject *));
      if (problem->users == NULL) {
        if (value->length != 1) {
          if (value->type == Integer) {
            EG_free(value->vals.integers);
          } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
            EG_free(value->vals.reals);
          } else if (value->type == String) {
            EG_free(value->vals.string);
          } else if (value->type == Tuple) {
            caps_freeTuple(value->length, value->vals.tuple);
          }
        }
        EG_free(value->partial);
        EG_free(value);
        EG_free(object);
        return EGADS_MALLOC;
      }
    } else {
      tmp = (capsObject **) EG_reall( problem->users,
                                     (problem->nUser+1)*sizeof(capsObject *));
      if (tmp == NULL) {
        if (value->length != 1) {
          if (value->type == Integer) {
            EG_free(value->vals.integers);
          } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
            EG_free(value->vals.reals);
          } else if (value->type == String) {
            EG_free(value->vals.string);
          } else if (value->type == Tuple) {
            caps_freeTuple(value->length, value->vals.tuple);
          }
        }
        EG_free(value->partial);
        EG_free(value);
        EG_free(object);
        return EGADS_MALLOC;
      }
      problem->users = tmp;
    }
    problem->users[problem->nUser] = object;
    problem->nUser += 1;
  }
  problem->sNum    += 1;
  object->last.sNum = problem->sNum;

  /* finish the object off */
  if ((units != NULL) &&
      (strlen(units) > 0)) value->units = EG_strdup(units);
  object->name    = EG_strdup(vname);
  object->type    = VALUE;
  object->subtype = stype;
  object->blind   = value;
/*@-kepttrans@*/
  object->parent  = pobject;
/*@+kepttrans@*/

  /* register for restart */
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\param.txt", problem->root);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx",  problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/param.txt",   problem->root);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",    problem->root);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Warning: Cannot open %s (caps_makeValue)\n", filename);
  } else {
    fprintf(fp, "%d %d\n", problem->nParam, problem->nUser);
    if (problem->params != NULL)
      for (i = 0; i < problem->nParam; i++)
        fprintf(fp, "%s\n", problem->params[i]->name);
    fclose(fp);
    status = caps_rename(temp, filename);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: Cannot rename %s!\n", filename);
  }
  status = caps_writeValueObj(problem, object);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeValueObj = %d (caps_makeValue)\n",
           status);

  *vobj = object;

  return CAPS_SUCCESS;
}


int
caps_makeValue(capsObject *pobject, const char *vname, enum capssType stype,
               enum capsvType vtype, int nrow, int ncol, const void *data,
               /*@null@*/ int *partial, /*@null@*/ const char *units,
               capsObject **vobj)
{
  int         i, stat, ret, vlen;
  CAPSLONG    sNum;
  capsProblem *problem;
  ut_unit     *utunit;
  capsJrnl    args[1];

  if (pobject == NULL)                         return CAPS_NULLOBJ;
  if (vobj    == NULL)                         return CAPS_NULLOBJ;
  *vobj = NULL;
  if (pobject->magicnumber != CAPSMAGIC)       return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)                return CAPS_BADTYPE;
  if (pobject->blind == NULL)                  return CAPS_NULLBLIND;
  if (vname == NULL)                           return CAPS_NULLNAME;
  if ((stype != PARAMETER) && (stype != USER)) return CAPS_BADTYPE;
  if (!((vtype == Double) || (vtype == DoubleDeriv) || (vtype == String))
      && (units != NULL))                      return CAPS_UNITERR;
  vlen = ncol*nrow;
  if (vlen <= 0)                               return CAPS_BADINDEX;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)                    return CAPS_READONLYERR;

  args[0].type = jObject;
  stat         = caps_jrnlRead(CAPS_MAKEVALUE, problem, *vobj, 1, args, &sNum,
                               &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) *vobj = args[0].members.obj;
    return ret;
  }

  /* check for unique name */
  ret = CAPS_SUCCESS;
  if ((stype == PARAMETER) && (problem->params != NULL))
    for (i = 0; i < problem->nParam; i++) {
      if (problem->params[i]       == NULL) continue;
      if (problem->params[i]->name == NULL) continue;
      if (strcmp(problem->params[i]->name, vname) == 0) ret = CAPS_BADNAME;
    }

  /* check the units */
  if ((units != NULL) && (ret == CAPS_SUCCESS)) {
    utunit = ut_parse((ut_system *) problem->utsystem, units, UT_ASCII);
    if (utunit == NULL) ret = CAPS_UNITERR;
    ut_free(utunit);
  }

  sNum = problem->sNum;
  if (ret == CAPS_SUCCESS)
    ret = caps_makeValueX(pobject, vname, stype, vtype, nrow, ncol, data,
                          partial, units, vobj);
  args[0].members.obj = *vobj;
  caps_jrnlWrite(CAPS_MAKEVALUE, problem, *vobj, ret, 1, args, sNum,
                 problem->sNum);

  return ret;
}


static int
caps_setValuX(capsObject *object, enum capsvType vtype, int nrow, int ncol,
              const void *data, /*@null@*/ const int *partial,
              /*@null@*/ const char *units, int *nErr, capsErrs **errors)
{
  int         i, j, k, n, vlen, slen=0, status;
  int         *ints  = NULL;
  double      *reals = NULL;
  char        *str   = NULL;
  capsTuple   *tuple = NULL;
  capsObject  *pobject, *source, *last;
  capsValue   *value;
  capsProblem *problem;

  *nErr   = 0;
  *errors = NULL;
  value   = (capsValue *) object->blind;
  vlen    = nrow*ncol;
  status  = caps_unitConvertible(units, value->units);
  if (status != CAPS_SUCCESS) return status;

  status = caps_findProblem(object, 9999, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  if (data == NULL) {
    if(value->lfixed == Change) vlen = 0;
    else                        vlen = value->length;
  }

  if ((value->type != vtype) && (data != NULL)) {
    if (!(((value->type == Double) || (value->type == DoubleDeriv)) &&
          ((vtype == Integer) || (vtype == Boolean))))
      return CAPS_BADTYPE;
  }
  if (value->nullVal == IsNull ||
      value->nullVal == IsPartial) value->nullVal = NotNull;

  /* are we in range? */
  if (data != NULL) {
    if (value->type == Integer) {
      if (value->limits.ilims[0] != value->limits.ilims[1])
        if (vtype == Integer) {
          ints = (int *) data;
          for (i = 0; i < vlen; i++)
            if ((ints[i] < value->limits.ilims[0]) ||
                (ints[i] > value->limits.ilims[1])) return CAPS_RANGEERR;
        } else {
          reals = (double *) data;
          for (i = 0; i < vlen; i++)
            if ((reals[i] < value->limits.ilims[0]) ||
                (reals[i] > value->limits.ilims[1])) return CAPS_RANGEERR;
        }
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      if (value->limits.dlims[0] != value->limits.dlims[1])
        if (vtype == Integer) {
          ints = (int *) data;
          for (i = 0; i < vlen; i++)
            if ((ints[i] < value->limits.dlims[0]) ||
                (ints[i] > value->limits.dlims[1])) return CAPS_RANGEERR;
        } else {
          reals = (double *) data;
          for (i = 0; i < vlen; i++)
            if ((reals[i] < value->limits.dlims[0]) ||
                (reals[i] > value->limits.dlims[1])) return CAPS_RANGEERR;
        }
    }
  }

  /* check for uniqueness in tuple names */
  if ((vtype == Tuple) && (data != NULL)) {
    tuple = (capsTuple *) data;
    for (i = 0; i < vlen; i++)
      for (j = i+1; j < vlen; j++)
        if (strcmp(tuple[i].name, tuple[j].name) == 0)
          return CAPS_BADVALUE;
  }

  /* remove any existing memory */
  if (vlen != value->length) {
    if ((value->lfixed == Fixed) && (value->type != String) &&
        (value->type != Pointer) && (value->type != Tuple) &&
        (value->type != PointerMesh)) return CAPS_SHAPEERR;
    if ((value->type == Boolean) || (value->type == Integer)) {
      if (vlen > 1) {
        ints = (int *) EG_alloc(vlen*sizeof(int));
        if (ints == NULL) return EGADS_MALLOC;
      }
      if (value->length > 1) EG_free(value->vals.integers);
      if (ints != NULL) value->vals.integers = ints;
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      if (vlen > 1) {
        reals = (double *) EG_alloc(vlen*sizeof(double));
        if (reals == NULL) return EGADS_MALLOC;
      }
      if (value->length > 1) EG_free(value->vals.reals);
      if (reals != NULL) value->vals.reals = reals;
    } else if (value->type == Tuple) {
      if (vlen > 0) {
        status = caps_makeTuple(vlen, &tuple);
        if ((status != CAPS_SUCCESS) || (tuple == NULL)) {
          if (tuple == NULL) status = CAPS_NULLVALUE;
          return status;
        }
      }
      caps_freeTuple(value->length, value->vals.tuple);
      value->vals.tuple = tuple;
    }

    if ((value->nullVal == IsPartial) || (partial != NULL)) {
      ints = (int *) EG_alloc(vlen*sizeof(int));
      if (ints == NULL) return EGADS_MALLOC;
      EG_free(value->partial);
      value->partial = ints;
      for (i = 0; i < vlen; i++) value->partial[i] = NotNull;
      value->nullVal = IsPartial;
    }

    value->length = vlen;
  } else {
    if (value->type == Tuple) {
      status = caps_makeTuple(vlen, &tuple);
      if ((status != CAPS_SUCCESS) || (tuple == NULL)) {
        if (tuple == NULL) status = CAPS_NULLVALUE;
        return status;
      }
      caps_freeTuple(value->length, value->vals.tuple);
      value->vals.tuple = tuple;
    }
  }

  /* set the values */
  if (data == NULL) {
    value->nullVal = IsNull;
  } else {
    if ((value->type == Boolean) || (value->type == Integer)) {
      ints = (int *) data;
      if (vlen == 1) {
        value->vals.integer = ints[0];
      } else {
        for (i = 0; i < vlen; i++) value->vals.integers[i] = ints[i];
      }
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      if ((vtype == Double) || (vtype == DoubleDeriv)) {
        reals = (double *) data;
        if (vlen == 1) {
          value->vals.real = reals[0];
          reals = &value->vals.real;
        } else {
          for (i = 0; i < vlen; i++) value->vals.reals[i] = reals[i];
          reals = value->vals.reals;
       }
      } else {
        ints = (int *) data;
        if (vlen == 1) {
          value->vals.real = ints[0];
          reals = &value->vals.real;
        } else {
          for (i = 0; i < vlen; i++) value->vals.reals[i] = ints[i];
          reals = value->vals.reals;
        }
      }
      status = caps_convert(vlen, units, reals, value->units, reals);
      if (status != CAPS_SUCCESS) return status;
    } else if (value->type == String) {
      /* always re-allocate strings */
      EG_free(value->vals.string);
      value->vals.string = NULL;
      str = (char *) data;
      for (slen = i = 0; i < vlen; i++)
        slen += strlen(str + slen)+1;
      value->vals.string = (char *) EG_alloc(slen*sizeof(char));
      if (value->vals.string == NULL) return EGADS_MALLOC;
      for (i = 0; i < slen; i++) value->vals.string[i] = str[i];
    } else if (value->type == Tuple) {
      tuple = (capsTuple *) data;
      for (i = 0; i < vlen; i++) {
        value->vals.tuple[i].name  = EG_strdup(tuple[i].name);
        value->vals.tuple[i].value = EG_strdup(tuple[i].value);
        if ((tuple[i].name  != NULL) &&
            (value->vals.tuple[i].name  == NULL)) return EGADS_MALLOC;
        if ((tuple[i].value != NULL) &&
            (value->vals.tuple[i].value == NULL)) return EGADS_MALLOC;
      }
    } else {
  /*@-kepttrans@*/
      value->vals.AIMptr = (void *) data;
  /*@+kepttrans@*/
    }
    if (partial != NULL) {
      if (value->partial == NULL) {
        ints = (int *) EG_alloc(vlen*sizeof(int));
        if (ints == NULL) return EGADS_MALLOC;
        value->partial = ints;
      }
      for (i = 0; i < vlen; i++) value->partial[i] = partial[i];
      value->nullVal = IsPartial;
    } else {
      EG_free(value->partial);
      value->partial = NULL;
    }
  }

  if (value->dim == Vector) {
    if (value->ncol == 1) value->nrow = vlen;
    else                  value->ncol = vlen;
  } else {
    value->nrow = nrow;
    value->ncol = ncol;
  }

  if (object->subtype != USER) {
    caps_freeOwner(&object->last);
    problem->sNum    += 1;
    object->last.sNum = problem->sNum;
    status = caps_addHistory(object, problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_addHistory = %d (caps_setValue)\n", status);
  }

  /* apply GeometryIn changes now -- direct and linkages */
  if ((object->subtype == GEOMETRYIN) && (data != NULL)) {
    if (vlen == 1) {
      reals = &value->vals.real;
    } else {
      reals = value->vals.reals;
    }
    for (n = k = 0; k < value->nrow; k++)
      for (j = 0; j < value->ncol; j++, n++) {
        if (value->partial != NULL) {
          status = ocsmSetValuD(problem->modl, value->pIndex, k+1, j+1,
                                value->partial[n] == NotNull ? reals[n] : -HUGEQ);
        } else {
          status = ocsmSetValuD(problem->modl, value->pIndex, k+1, j+1,
                                reals[n]);
        }
        if (status != SUCCESS) {
          printf(" CAPS Error: Cant change %s[%d,%d] = %d (caps_setValue)!\n",
                   object->name, k+1, j+1, status);
          return CAPS_BADVALUE;
        }
      }
  } else if (data != NULL) {
    for (i = 0; i < problem->nGeomIn; i++) {
      source = problem->geomIn[i];
      last   = NULL;
      do {
        if (source->magicnumber != CAPSMAGIC)  break;
        if (source->type        != VALUE)      break;
        if (source->blind       == NULL)       break;
        value = (capsValue *) source->blind;
        if (value->link == problem->geomIn[i]) break;
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      if (last != object) continue;
      /* we hit our object from a GeometryIn link */
      source = problem->geomIn[i];
      value  = (capsValue *) source->blind;
      if ((value->nrow != nrow) || (value->ncol != ncol)) {
        printf(" CAPS Warning: Shape problem with link %s %s (caps_setValue)!\n",
               object->name, source->name);
        continue;
      }
      if ((vtype == Double) || (vtype == DoubleDeriv)) {
        reals = (double *) data;
        if (vlen == 1) {
          value->vals.real = reals[0];
          reals = &value->vals.real;
        } else {
          for (i = 0; i < vlen; i++) value->vals.reals[i] = reals[i];
          reals = value->vals.reals;
        }
      } else {
        ints = (int *) data;
        if (vlen == 1) {
          value->vals.real = ints[0];
          reals = &value->vals.real;
        } else {
          for (i = 0; i < vlen; i++) value->vals.reals[i] = ints[i];
          reals = value->vals.reals;
        }
      }
      status = caps_convert(vlen, units, reals, value->units, reals);
      if (status != CAPS_SUCCESS) return status;
      for (n = k = 0; k < nrow; k++)
        for (j = 0; j < ncol; j++, n++) {
          if (value->nullVal == IsPartial) {
            status = ocsmSetValuD(problem->modl, value->pIndex, k+1, j+1,
                                  value->partial[n] == NotNull ? reals[n] : -HUGEQ);
          } else {
            status = ocsmSetValuD(problem->modl, value->pIndex, k+1, j+1,
                                  reals[n]);
          }
          if (status != SUCCESS) {
            printf(" CAPS Error: Cant change %s[%d,%d] = %d (caps_setValue)!\n",
                   object->name, k+1, j+1, status);
            return CAPS_BADVALUE;
          }
        }
      caps_freeOwner(&source->last);
      source->last       = object->last;
      source->last.pname = EG_strdup(object->last.pname);
      source->last.pID   = EG_strdup(object->last.pID);
      source->last.user  = EG_strdup(object->last.user);
    }
  }

  /* register the change for restarts */
  if (object->subtype != USER) {
    status = caps_writeValueObj(problem, object);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_writeValueObj = %d (caps_setValue)\n",
             status);
  }

  return CAPS_SUCCESS;
}


int
caps_setValue(capsObject *object, enum capsvType vtype, int nrow, int ncol,
              /*@null@*/ const void *data, /*@null@*/ const int *partial,
              /*@null@*/ const char *units, int *nErr, capsErrs **errors)
{
  int         i, stat, ret, vlen;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsValue   *value;
  capsProblem *problem;
  capsJrnl    args[2];

  if (nErr                == NULL)        return CAPS_NULLVALUE;
  if (errors              == NULL)        return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if (object              == NULL)        return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC)   return CAPS_BADOBJECT;
  if (object->type        != VALUE)       return CAPS_BADTYPE;
  if (object->subtype     == GEOMETRYOUT) return CAPS_BADTYPE;
  if (object->subtype     == ANALYSISOUT) return CAPS_BADTYPE;
  if (object->blind       == NULL)        return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;
  if (object->subtype     == GEOMETRYIN)
    if (value->gInType    == 2)           return CAPS_READONLYERR;
  if (value->link         != NULL)        return CAPS_LINKERR;
  if (data == NULL) {
    if (value->nullVal == NotAllowed)     return CAPS_NULLVALUE;
  } else {
    vlen = nrow*ncol;
    if (vlen <= 0)                        return CAPS_RANGEERR;
    if (value->sfixed == Fixed) {
      if (value->dim == Scalar) {
        if (vlen > 1)                     return CAPS_SHAPEERR;
      } else if (value->dim == Vector) {
        if ((ncol != 1) && (nrow != 1))   return CAPS_SHAPEERR;
      }
    }
    if ((value->lfixed == Fixed) && (vlen != value->length))
                                          return CAPS_SHAPEERR;
    if (partial != NULL) {
      for (i = 0; i < vlen; i++)
        if (!(partial[i] == NotNull ||
              partial[i] == IsNull))        return CAPS_NULLVALUE;
    }
  }

  stat = caps_findProblem(object, CAPS_SETVALUE, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  if (( object->subtype == GEOMETRYIN) &&
      (pobject->subtype == STATIC    ))   return CAPS_READONLYERR;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)               return CAPS_READONLYERR;

  args[0].type = jInteger;
  args[1].type = jErr;
  stat         = caps_jrnlRead(CAPS_SETVALUE, problem, object, 2, args, &sNum,
                               &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *nErr   = args[0].members.integer;
    *errors = args[1].members.errs;
    return ret;
  }

  sNum = problem->sNum;
  ret  = caps_setValuX(object, vtype, nrow, ncol, data, partial, units, nErr,
                       errors);
  args[0].members.integer = *nErr;
  args[1].members.errs    = *errors;
  caps_jrnlWrite(CAPS_SETVALUE, problem, object, ret, 2, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_getLimits(const capsObj object, enum capsvType *vtype, const void **limits,
               const char **units)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsValue   *value;
  capsProblem *problem;
  capsJrnl    args[3];

  if (limits              == NULL)      return CAPS_NULLVALUE;
  *limits = NULL;
  if (units               == NULL)      return CAPS_NULLVALUE;
  *units = NULL;
  if (vtype               == NULL)      return CAPS_NULLVALUE;
  *vtype = 0;

  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->type        != VALUE)     return CAPS_BADTYPE;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;
  if ((value->type != Integer) && (value->type != Double) &&
      (value->type != DoubleDeriv))     return CAPS_BADTYPE;

  status = caps_findProblem(object, CAPS_GETLIMITS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jInteger;
  args[1].type = jPointer;
  args[2].type = jString;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_GETLIMITS, problem, object, 3, args, &sNum,
                               &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      *vtype  = args[0].members.integer;
      *limits = args[1].members.pointer;
      *units  = args[2].members.string;
      return ret;
    }
  }

  *vtype         = value->type;
  *units         = value->units;
  args[1].length = 0;
  if (value->type == Integer) {
    if (value->limits.ilims[0] == value->limits.ilims[1]) goto complete;
    *limits = value->limits.ilims;
    args[1].length = 2*sizeof(int);
  } else {
    if (value->limits.dlims[0] == value->limits.dlims[1]) goto complete;
    *limits = value->limits.dlims;
    args[1].length = 2*sizeof(double);
  }

complete:
  if (problem->dbFlag == 1) return CAPS_SUCCESS;
  args[0].members.integer =          *vtype;
  args[1].members.pointer = (void *) *limits;
  args[2].members.string  = (char *) *units;
  caps_jrnlWrite(CAPS_GETLIMITS, problem, object, CAPS_SUCCESS, 3, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}


int
caps_setLimits(capsObject *object, enum capsvType vtype, void *limits,
               const char *units, int *nErr, capsErrs **errors)
{
  int         *ints, i, status;
  double      *realp, reals[2];
  capsValue   *value;
  capsObject  *pobject;
  capsProblem *problem;

  if (nErr                 == NULL)       return CAPS_NULLVALUE;
  if (errors               == NULL)       return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if  (object              == NULL)       return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if  (object->type        != VALUE)      return CAPS_BADTYPE;
  if ((object->subtype     != USER) && (object->subtype != PARAMETER))
                                          return CAPS_BADTYPE;
  if  (object->blind       == NULL)       return CAPS_NULLBLIND;
  value   = (capsValue *) object->blind;
  status  = caps_findProblem(object, CAPS_SETLIMITS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)               return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  if (limits == NULL) {
    value->limits.ilims[0] = value->limits.ilims[1] = 0;
    value->limits.dlims[0] = value->limits.dlims[1] = 0.0;
    return CAPS_SUCCESS;
  }
  if (value->type != vtype) {
    if (!( ((value->type == Double) || (value->type == DoubleDeriv)) &&
            (vtype == Integer) )) return CAPS_BADTYPE;
  }

  status = caps_unitConvertible(units, value->units);
  if (status != CAPS_SUCCESS) return status;

  if (value->type == Integer) {
    if ( units != NULL ) return CAPS_UNITERR;
    ints = (int *) limits;
    if (ints[0] >= ints[1]) return CAPS_RANGEERR;
    if (value->length == 1) {
      if ((value->vals.integer < ints[0]) ||
          (value->vals.integer > ints[1])) return CAPS_RANGEERR;
    } else {
      for (i = 0; i < value->length; i++)
        if ((value->vals.integers[i] < ints[0]) ||
            (value->vals.integers[i] > ints[1])) return CAPS_RANGEERR;
    }
    value->limits.ilims[0] = ints[0];
    value->limits.ilims[1] = ints[1];
  } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
    if (vtype == Integer) {
      ints = (int *) limits;
      reals[0] = ints[0];
      reals[1] = ints[1];
    } else {
      realp = (double *) limits;
      reals[0] = realp[0];
      reals[1] = realp[1];
    }
    if (reals[0] >= reals[1]) return CAPS_RANGEERR;

    status = caps_convert(2, units, reals, value->units, reals);
    if (status != CAPS_SUCCESS) return status;

    if (value->length == 1) {
      if ((value->vals.real < reals[0]) ||
          (value->vals.real > reals[1])) return CAPS_RANGEERR;
    } else {
      for (i = 0; i < value->length; i++)
        if ((value->vals.reals[i] < reals[0]) ||
            (value->vals.reals[i] > reals[1])) return CAPS_RANGEERR;
    }
    value->limits.dlims[0] = reals[0];
    value->limits.dlims[1] = reals[1];
  } else {
    return CAPS_BADTYPE;
  }

  return CAPS_SUCCESS;
}


int
caps_getValueProps(capsObject *object, int *dim, int *pmtr,
                   enum capsFixed *lfixed, enum capsFixed *sfixed,
                   enum capsNull *nval)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsValue   *value;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[5];

  if (dim    == NULL) return CAPS_NULLVALUE;
  if (pmtr   == NULL) return CAPS_NULLVALUE;
  if (lfixed == NULL) return CAPS_NULLVALUE;
  if (sfixed == NULL) return CAPS_NULLVALUE;
  if (nval   == NULL) return CAPS_NULLVALUE;

  *dim    = *pmtr   = 0;
  *lfixed = *sfixed = Fixed;
  if (object              == NULL)      return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (object->type        != VALUE)     return CAPS_BADTYPE;
  if (object->blind       == NULL)      return CAPS_NULLBLIND;
  value  = (capsValue *) object->blind;
  status = caps_findProblem(object, CAPS_GETVALUEPROPS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jInteger;
  args[3].type = jInteger;
  args[4].type = jInteger;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_GETVALUEPROPS, problem, object, 5, args,
                               &sNum, &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      *dim    = args[0].members.integer;
      *pmtr   = args[1].members.integer;
      *lfixed = args[2].members.integer;
      *sfixed = args[3].members.integer;
      *nval   = args[4].members.integer;
      return ret;
    }
  }

  args[0].members.integer = *dim    = value->dim;
  args[1].members.integer = *pmtr   = value->gInType;
  args[2].members.integer = *lfixed = value->lfixed;
  args[3].members.integer = *sfixed = value->sfixed;
  args[4].members.integer = *nval   = value->nullVal;
  if (problem->dbFlag == 1) return CAPS_SUCCESS;

  caps_jrnlWrite(CAPS_GETVALUEPROPS, problem, object, CAPS_SUCCESS, 5, args,
                 problem->sNum, problem->sNum);

  return CAPS_SUCCESS;
}


int
caps_setValueProps(capsObject *object, int dim, enum capsFixed lfixed,
                   enum capsFixed sfixed, enum capsNull nval,
                   int *nErr, capsErrs **errors)
{
  int         status;
  capsValue   *value;
  capsObject  *pobject;
  capsProblem *problem;

  if (nErr                 == NULL)       return CAPS_NULLVALUE;
  if (errors               == NULL)       return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if  (object              == NULL)       return CAPS_NULLOBJ;
  if  (object->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  if  (object->type        != VALUE)      return CAPS_BADTYPE;
  if ((object->subtype != PARAMETER) &&
      (object->subtype != USER))          return CAPS_BADTYPE;
  if  (object->blind == NULL)             return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;
  if ((nval == NotAllowed) &&
      (value->nullVal == IsNull))         return CAPS_NULLVALUE;
  status  = caps_findProblem(object, CAPS_SETVALUEPROPS, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)               return CAPS_READONLYERR;

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)             return CAPS_SUCCESS;
  }

  /* check dimension */
  if (dim == 0) {
    if ((value->ncol != 1) || (value->nrow != 1)) return CAPS_SHAPEERR;
  } else if (dim == 1) {
    if ((value->ncol >  1) && (value->nrow >  1)) return CAPS_SHAPEERR;
  } else if (dim != 2) {
    return CAPS_RANGEERR;
  }

  value->dim     = dim;
  value->lfixed  = lfixed;
  value->sfixed  = sfixed;
  if ((nval == IsNull)  && (value->nullVal == NotNull)) return CAPS_SUCCESS;
  if ((nval == NotNull) && (value->nullVal == IsNull))  return CAPS_SUCCESS;
  value->nullVal = nval;

  return CAPS_SUCCESS;
}


int
caps_convertValue(capsObject *object, double inp, const char *units,
                  double *outp)
{
  int         status, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsValue   *value;
  capsProblem *problem;
  capsJrnl    args[1];

  if (object              == NULL)        return CAPS_NULLOBJ;
  if (object->magicnumber != CAPSMAGIC)   return CAPS_BADOBJECT;
  if (object->type        != VALUE)       return CAPS_BADTYPE;
  if (object->blind       == NULL)        return CAPS_NULLBLIND;
  value = (capsValue *) object->blind;
  if (units               == NULL)        return CAPS_UNITERR;
  if (value->units        == NULL)        return CAPS_UNITERR;
  status = caps_findProblem(object, CAPS_CONVERTVALUE, &pobject);
  if (status != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag == 1)               return CAPS_READONLYERR;

  args[0].type = jDouble;
  status       = caps_jrnlRead(CAPS_CONVERTVALUE, problem, object, 1, args,
                               &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    *outp = args[0].members.real;
    return ret;
  }

  ret = caps_convert(1, units, &inp, value->units, outp);
  args[0].members.real = *outp;
  caps_jrnlWrite(CAPS_CONVERTVALUE, problem, object, ret, 1, args,
                 problem->sNum, problem->sNum);

  return ret;
}


int
caps_dupValues(capsValue *val1, capsValue *val2)
{
  int i, j, len, status;

  /* this function assumes that val2 has no allocated data */
  val2->type    = val1->type;
  val2->length  = val1->length;
  val2->dim     = val1->dim;
  val2->nrow    = val1->nrow;
  val2->ncol    = val1->ncol;
  val2->lfixed  = val1->lfixed;
  val2->sfixed  = val1->sfixed;
  val2->nullVal = val1->nullVal;
  val2->index   = val1->index;
  val2->pIndex  = val1->pIndex;
  val2->gInType = val1->gInType;
  val2->partial = NULL;
  val2->nderiv    = 0;
  val2->derivs    = NULL;
  if (val1->partial != NULL) {
    val2->partial = (int *) EG_alloc(val1->length*sizeof(int));
    if (val2->partial == NULL) return EGADS_MALLOC;
    for (i = 0; i < val1->length; i++) val2->partial[i] = val1->partial[i];
  }
  if ((val1->nderiv != 0) && (val1->derivs != NULL)) {
    val2->derivs = (capsDeriv *) EG_alloc(val1->nderiv*sizeof(capsDeriv));
    if (val2->derivs == NULL) return EGADS_MALLOC;
    for (i = 0; i < val1->nderiv; i++) {
      val2->derivs[i].name  = EG_strdup(val1->derivs[i].name);
      val2->derivs[i].len_wrt  = val1->derivs[i].len_wrt;
      val2->derivs[i].deriv = NULL;
    }
    for (i = 0; i < val1->nderiv; i++) {
      if (val2->derivs[i].name == NULL) return EGADS_MALLOC;
      val2->derivs[i].deriv = (double *) EG_alloc(
                           val1->length*val1->derivs[i].len_wrt*sizeof(double));
      if (val2->derivs[i].deriv == NULL) return EGADS_MALLOC;
      for (j = 0; j < val1->length*val1->derivs[i].len_wrt; j++)
        val2->derivs[i].deriv[j] = val1->derivs[i].deriv[j];
    }
    val2->nderiv = val1->nderiv;
  }

  /* set the actual value(s) */
  switch (val1->type) {
    case Boolean:
      if (val1->length == 1) {
        val2->vals.integer  = val1->vals.integer;
      } else {
    	val2->vals.integers = NULL;
    	if (val1->vals.integers == NULL) break;
        val2->vals.integers = (int *) EG_alloc(val1->length*sizeof(int));
        if (val2->vals.integers == NULL) return EGADS_MALLOC;
        for (i = 0; i < val1->length; i++)
          val2->vals.integers[i] = val1->vals.integers[i];
      }
      break;

    case Integer:
      if (val1->length == 1) {
        val2->vals.integer  = val1->vals.integer;
      } else {
        val2->vals.integers = NULL;
        if (val1->vals.integers != NULL) {
          val2->vals.integers = (int *) EG_alloc(val1->length*sizeof(int));
          if (val2->vals.integers == NULL) return EGADS_MALLOC;
          for (i = 0; i < val1->length; i++)
            val2->vals.integers[i] = val1->vals.integers[i];
        }
      }
      val2->limits.ilims[0] = val1->limits.ilims[0];
      val2->limits.ilims[1] = val1->limits.ilims[1];
      break;

    case Double:
      if (val1->length == 1) {
        val2->vals.real  = val1->vals.real;
      } else {
        val2->vals.reals = NULL;
        if (val1->vals.reals != NULL) {
          val2->vals.reals = (double *) EG_alloc(val1->length*sizeof(double));
          if (val2->vals.reals == NULL) return EGADS_MALLOC;
          for (i = 0; i < val1->length; i++)
            val2->vals.reals[i] = val1->vals.reals[i];
        }
      }
      val2->limits.dlims[0] = val1->limits.dlims[0];
      val2->limits.dlims[1] = val1->limits.dlims[1];
      break;

    case String:
      val2->vals.string = NULL;
      if ((val1->vals.string != NULL)  && (strlen(val1->vals.string) > 0)) {
        j = strlen(val1->vals.string) + 1;
        val2->vals.string = (char *) EG_alloc(j*sizeof(char));
        if (val2->vals.string == NULL) return EGADS_MALLOC;
        for (i = 0; i < j; i++)
          val2->vals.string[i] = val1->vals.string[i];
      }
      break;

    case Tuple:
      val2->vals.tuple = NULL;
      if (val1->vals.tuple == NULL)
        break;
      status = caps_makeTuple(val1->length, &val2->vals.tuple);
      if (status != CAPS_SUCCESS)   return status;
      if (val2->vals.tuple == NULL) return EGADS_MALLOC;        /* for splint */
      for (i = 0; i < val1->length; i++) {
        val2->vals.tuple[i].name  = EG_strdup(val1->vals.tuple[i].name);
        val2->vals.tuple[i].value = EG_strdup(val1->vals.tuple[i].value);
        if ((val1->vals.tuple[i].name  != NULL) &&
            (val2->vals.tuple[i].name  == NULL)) return EGADS_MALLOC;
        if ((val1->vals.tuple[i].value != NULL) &&
            (val2->vals.tuple[i].value == NULL)) return EGADS_MALLOC;
      }
      break;

    case Pointer:
    case PointerMesh:
      val2->vals.AIMptr = val1->vals.AIMptr;

  }

  /* copy units */
  val2->units = NULL;
  if (val1->units != NULL) {
    len = strlen(val1->units) + 1;
    val2->units = (char *) EG_alloc(len*sizeof(char));
    if (val2->units == NULL) return EGADS_MALLOC;
    for (i = 0; i < len; i++) val2->units[i] = val1->units[i];
  }

  /* copy mesh writer */
  val2->meshWriter = NULL;
  if (val1->meshWriter != NULL) {
    len = strlen(val1->meshWriter) + 1;
    val2->meshWriter = (char *) EG_alloc(len*sizeof(char));
    if (val2->meshWriter == NULL) return EGADS_MALLOC;
    for (i = 0; i < len; i++) val2->meshWriter[i] = val1->meshWriter[i];
  }

  /* linkage */
  val2->link       = val1->link;
  val2->linkMethod = val1->linkMethod;

  return CAPS_SUCCESS;
}


static int
caps_compatValues(capsValue *val1, capsValue *val2, capsProblem *problem)
{
  int     status;
  ut_unit *utunit1, *utunit2;

  /* check units */
  if ((val1->units != NULL) && (val2->units == NULL)) return CAPS_UNITERR;
  if ((val1->units == NULL) && (val2->units != NULL)) return CAPS_UNITERR;
  if ((val1->units != NULL) && (val2->units != NULL)) {
    if (strcmp(val1->units, val2->units) != 0) {
      utunit1 = ut_parse((ut_system *) problem->utsystem, val1->units, UT_ASCII);
      utunit2 = ut_parse((ut_system *) problem->utsystem, val2->units, UT_ASCII);
      status  = ut_are_convertible(utunit1, utunit2);
      ut_free(utunit1);
      ut_free(utunit2);
      if (status == 0) return CAPS_UNITERR;
    }
  }

  /* check type */
  if (val1->type != val2->type)                   return CAPS_BADTYPE;

  /* check shape */
  if (val2->lfixed == Fixed)
    if (val1->length != val2->length)             return CAPS_SHAPEERR;
  if (val2->sfixed == Fixed) {
    if (val2->dim    >  val1->dim)                return CAPS_SHAPEERR;
    if (val2->nrow   != val1->nrow)               return CAPS_SHAPEERR;
    if (val2->ncol   != val1->ncol)               return CAPS_SHAPEERR;
  } else {
    if (val2->dim == Scalar) {
      if (val1->length != 1)                      return CAPS_SHAPEERR;
    } else if (val2->dim == Vector) {
      if ((val1->ncol != 1) && (val1->nrow != 1)) return CAPS_SHAPEERR;
    }
  }

  return CAPS_SUCCESS;
}


int
caps_transferValueX(capsObject *source, enum capstMethod method,
                    capsObject *target, int *nErr, capsErrs **errors)
{
  int            status, vlen, rank, ncol, nrow;
  char           *iunits;
  double         *ireals;
  const void     *data;
  const int      *partial;
  const char     *units;
  const double   *reals;
  capsValue      *value, *sval;
  capsProblem    *problem;
  capsObject     *pobject, *last;
  enum capsvType vtype;

  if  (nErr                == NULL)         return CAPS_NULLVALUE;
  if  (errors              == NULL)         return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if  (source              == NULL)         return CAPS_NULLOBJ;
  if  (source->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if ((source->type        != VALUE) &&
      (source->type        != DATASET))     return CAPS_BADTYPE;
  if  (source->blind       == NULL)         return CAPS_NULLBLIND;
  if  (target              == NULL)         return CAPS_NULLOBJ;
  if  (target->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if  (target->type        != VALUE)        return CAPS_BADTYPE;
  if  (target->subtype     == GEOMETRYOUT)  return CAPS_BADTYPE;
  if  (target->subtype     == ANALYSISOUT)  return CAPS_BADTYPE;
  if  (target->blind       == NULL)         return CAPS_NULLBLIND;
  value   = (capsValue *)   target->blind;
  status  = caps_findProblem(target, 9999, &pobject);
  if (status               != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  if (source->type == VALUE) {

    status = caps_getValuX(source, &vtype, &nrow, &ncol, &data, &partial,
                           &units, nErr, errors);
    if (status != CAPS_SUCCESS) return status;
    last = value->link;
    value->link = NULL;
    status = caps_setValuX(target, vtype, nrow, ncol, data, partial, units,
                           nErr, errors);
    value->link = last;
    if (status != CAPS_SUCCESS) return status;

  } else {

    if (method == Copy) {

      sval = (capsValue *) source->blind;
      if (sval->nullVal == IsNull) return CAPS_NULLVALUE;
      status = caps_getDataX(source, &vlen, &rank, &reals, &units);
      if (status != CAPS_SUCCESS) return status;
      /* check for compatibility */
      status = caps_makeVal(Double, vlen*rank, reals, &sval);
      if (status != CAPS_SUCCESS) return status;
      sval->units = (char *) units;
      sval->ncol  = vlen;
      sval->nrow  = rank;
      if (rank != 1) sval->dim = Array2D;
      status = caps_compatValues(sval, value, problem);
      if (status != CAPS_SUCCESS) {
        if (sval->length > 1) EG_free(sval->vals.reals);
        EG_free(sval);
        return status;
      }
      value->ncol = vlen;
      value->nrow = rank;
      if (rank != 1) value->dim = Array2D;
      last = value->link;
      value->link = NULL;
      status = caps_setValuX(target, Double, rank, vlen, reals, NULL, units,
                             nErr, errors);
      value->link = last;
      if (sval->length > 1) EG_free(sval->vals.reals);
      EG_free(sval);
      if (status != CAPS_SUCCESS) return status;

    } else {

      if (source->subtype == UNCONNECTED) return CAPS_BADMETHOD;

      /* Integrate / weighted average */
      status = caps_integrateData(source, method, &rank, &ireals, &iunits);
      if (status != CAPS_SUCCESS) return status;
      status = caps_makeVal(Double, rank, ireals, &sval);
      if (status != CAPS_SUCCESS) {
        EG_free(ireals);
        EG_free(iunits);
        return status;
      }
      sval->units = iunits;
      sval->nrow  = rank;
      sval->ncol  = 1;
      status = caps_compatValues(sval, value, problem);
      EG_free(sval->units);
      if (sval->length > 1) EG_free(sval->vals.reals);
      EG_free(sval);
      if (status != CAPS_SUCCESS) {
        EG_free(ireals);
        return status;
      }
      last = value->link;
      value->link = NULL;
      status = caps_setValuX(target, Double, rank, 1, ireals, NULL, iunits,
                             nErr, errors);
      value->link = last;
      EG_free(ireals);
      if (status != CAPS_SUCCESS) return status;
    }

  }

  /* mark owner */
  caps_freeOwner(&target->last);
  problem->sNum    += 1;
  target->last.sNum = problem->sNum;
  status = caps_addHistory(target, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_transferValues)\n",
           status);

  status = caps_writeValueObj(problem, target);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeValueObj = %d (caps_transferValues)\n",
           status);

  /* invalidate any link if exists */
  value->linkMethod = Copy;
  value->link       = NULL;

  return CAPS_SUCCESS;
}


int
caps_transferValues(capsObject *source, enum capstMethod method,
                    capsObject *target, int *nErr, capsErrs **errors)
{
  int         stat, ret;
  CAPSLONG    sNum;
  capsObject  *pobject;
  capsProblem *problem;
  capsJrnl    args[2];

  if  (nErr                == NULL)         return CAPS_NULLVALUE;
  if  (errors              == NULL)         return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if  (source              == NULL)         return CAPS_NULLOBJ;
  if  (source->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if ((source->type        != VALUE) &&
      (source->type        != DATASET))     return CAPS_BADTYPE;
  if  (source->blind       == NULL)         return CAPS_NULLBLIND;
  if  (target              == NULL)         return CAPS_NULLOBJ;
  if  (target->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if  (target->type        != VALUE)        return CAPS_BADTYPE;
  if  (target->subtype     == GEOMETRYOUT)  return CAPS_BADTYPE;
  if  (target->subtype     == ANALYSISOUT)  return CAPS_BADTYPE;
  if  (target->blind       == NULL)         return CAPS_NULLBLIND;
  stat    = caps_findProblem(target, CAPS_TRANSFERVALUES, &pobject);
  if (stat                 != CAPS_SUCCESS) return stat;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag      == 1)            return CAPS_READONLYERR;

  args[0].type = jInteger;
  args[1].type = jErr;
  stat         = caps_jrnlRead(CAPS_TRANSFERVALUES, problem, target, 2, args,
                               &sNum, &ret);
  if (stat == CAPS_JOURNALERR) return stat;
  if (stat == CAPS_JOURNAL) {
    *nErr   = args[0].members.integer;
    *errors = args[1].members.errs;
    return ret;
  }

  if (ret == CAPS_SUCCESS) {
    ret   = caps_transferValueX(source, method, target, nErr, errors);
    *nErr = 0;
    if (*errors != NULL) *nErr = (*errors)->nError;
  }
  args[0].members.integer = *nErr;
  args[1].members.errs    = *errors;
  caps_jrnlWrite(CAPS_TRANSFERVALUES, problem, target, ret, 2, args, sNum,
                 problem->sNum);

  return ret;
}


int
caps_linkValue(/*@null@*/ capsObject *link, enum capstMethod method,
               capsObject *target, int *nErr, capsErrs **errors)
{
  int          status, rank, vlen;
  const double *reals;
  const char   *units;
  capsValue    *value, *sval;
  capsObject   *source, *pobject;
  capsProblem  *problem;

  if (nErr                == NULL)         return CAPS_NULLVALUE;
  if (errors              == NULL)         return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  if (target              == NULL)         return CAPS_NULLOBJ;
  if (target->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (target->type        != VALUE)        return CAPS_BADTYPE;
  if (target->subtype     == GEOMETRYOUT)  return CAPS_BADTYPE;
  if (target->subtype     == ANALYSISOUT)  return CAPS_BADTYPE;
  if (target->blind       == NULL)         return CAPS_NULLBLIND;
  value  = (capsValue *)   target->blind;
  status = caps_findProblem(target, CAPS_LINKVALUE, &pobject);
  if (status              != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  if (problem->dbFlag     == 1)            return CAPS_READONLYERR;

  if (target->type == VALUE) {
    if (target->subtype == GEOMETRYIN) {
      if (pobject->subtype == STATIC)      return CAPS_READONLYERR;
      if (method != Copy)                  return CAPS_BADMETHOD;
    }
  }

  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR)  return CAPS_JOURNALERR;
  if (problem->stFlag == 4) {
    status = caps_jrnlEnd(problem);
    if (status != CAPS_CLEAN)              return CAPS_SUCCESS;
  }

  if (link == NULL) {
    /* mark owner */
    caps_freeOwner(&target->last);
    problem->sNum    += 1;
    target->last.sNum = problem->sNum;
    status = caps_addHistory(target, problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_addHistory = %d (caps_linkValue)\n", status);
    /* remove existing link */
    value->linkMethod = Copy;
    value->link       = NULL;
    status = caps_writeSerialNum(problem);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_writePSerialNum = %d (caps_linkValue)\n",
             status);
    if (target->type == VALUE) {
      status = caps_writeValueObj(problem, target);
      if (status != CAPS_SUCCESS)
        printf(" CAPS Warning: caps_writeValueObj = %d (caps_linkValue)\n",
               status);
    }
    return CAPS_SUCCESS;
  }

  /* look at link */
  if (link->magicnumber != CAPSMAGIC)       return CAPS_BADOBJECT;
  if (link->type        == VALUE) {

    if (target->subtype == USER)            return CAPS_BADTYPE;
    if (method != Copy)                     return CAPS_BADMETHOD;
    source = link;
    do {
      if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (source->type        != VALUE)     return CAPS_BADTYPE;
      if (source->blind       == NULL)      return CAPS_NULLBLIND;
      sval = (capsValue *) source->blind;
      if (sval->link          == target)    return CAPS_CIRCULARLINK;
      source = sval->link;
    } while (sval->link != NULL);

    if ((sval->type == Pointer) || (sval->type == PointerMesh)) {
      /* must have the same "units" */
      if  ((sval->units  != NULL) || (value->units != NULL)) {
        if (sval->units  == NULL) return CAPS_UNITERR;
        if (value->units == NULL) return CAPS_UNITERR;
        if (strcmp(value->units, sval->units) != 0) return CAPS_UNITERR;
      }
    } else {
      /* make sure we are compatible */
      status = caps_compatValues(sval, value, problem);
      if (status != CAPS_SUCCESS) return status;
    }

  } else if (link->type == DATASET) {

    status = caps_getDataX(link, &vlen, &rank, &reals, &units);
    if (status != CAPS_SUCCESS) return status;

    /* check for compatibility */
    status = caps_makeVal(Double, vlen*rank, reals, &sval);
    if (status != CAPS_SUCCESS) return status;
    sval->units = (char *) units;
    if (rank != 1) {
      sval->dim  = Array2D;
      sval->nrow = rank;
      sval->ncol = vlen;
    }
    status = caps_compatValues(sval, value, problem);
    if (sval->length > 1) EG_free(sval->vals.reals);
    EG_free(sval);
    if (status != CAPS_SUCCESS) return status;

  } else {
    return CAPS_BADTYPE;
  }

  /* set the link */
  value->linkMethod = method;
  value->link       = link;

  /* check for circular auto execution links */
  if (target->subtype == ANALYSISIN) {
    status = caps_circularAutoExecs(target, NULL);
    if (status != CAPS_SUCCESS) {
      value->linkMethod = Copy;
      value->link       = NULL;
      return CAPS_CIRCULARLINK;
    }
  }

  caps_freeOwner(&target->last);
  problem->sNum    += 1;
  target->last.sNum = problem->sNum;
  status = caps_addHistory(target, problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_addHistory = %d (caps_linkValue)\n", status);
  status = caps_writeSerialNum(problem);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeSerialNum = %d (caps_linkValue)!\n",
           status);
  status = caps_writeValueObj(problem, target);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeValueObj = %d (caps_linkValue)!\n",
           status);

  return status;
}


int
caps_hasDeriv(capsObject *vobj, int *nderiv, char ***names,
              int *nErr, capsErrs **errors)
{
  int         i, status, ret;
  char        **namex, **namey;
  CAPSLONG    sNum;
  capsValue   *value;
  capsObject  *pobject, *source, *last;
  capsProblem *problem;
  capsJrnl    args[1];

  *nErr = 0;
  *errors = NULL;
  *nderiv  = 0;
  *names = NULL;
  if (vobj              == NULL)        return CAPS_NULLOBJ;
  if (vobj->magicnumber != CAPSMAGIC)   return CAPS_BADOBJECT;
  if (vobj->type        != VALUE)       return CAPS_BADTYPE;
  if (vobj->blind       == NULL)        return CAPS_NULLBLIND;
  status = caps_findProblem(vobj, CAPS_HASDERIV, &pobject);
  if (status            != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jStrings;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_HASDERIV, problem, vobj, 1, args, &sNum,
                               &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      if (ret == CAPS_SUCCESS) {
        *nderiv = args[0].num;
        namey   = args[0].members.strings;
        namex   = (char **) EG_alloc(args[0].num*sizeof(char *));
        if (namex == NULL) {
          ret = EGADS_MALLOC;
        } else {
          for (i = 0; i < args[0].num; i++) namex[i] = namey[i];
        }
        *names = namex;
      }
      return ret;
    }
  }

  args[0].num             = 0;
  args[0].members.strings = NULL;
  source = vobj;
  do {
    if (source->magicnumber != CAPSMAGIC) {
      status = CAPS_BADOBJECT;
      goto herr;
    }
    if (source->type        != VALUE)     {
      status = CAPS_BADTYPE;
      goto herr;
    }
    if (source->blind       == NULL)      {
      status = CAPS_NULLBLIND;
      goto herr;
    }
    value = (capsValue *) source->blind;
    if (value->link == vobj)              {
      status = CAPS_CIRCULARLINK;
      goto herr;
    }
    last   = source;
    source = value->link;
  } while (value->link != NULL);

  /* do we need to update our value? */
  if (problem->dbFlag == 0)
    if (last->subtype == ANALYSISOUT) {
      status = caps_updateAnalysisOut(last, CAPS_HASDERIV, nErr, errors);
      if (status != CAPS_SUCCESS) goto herr;
    } else if (vobj->subtype == GEOMETRYOUT) {
      /* make sure geometry is up-to-date */
      status = caps_build(pobject, nErr, errors);
      if ((status != CAPS_SUCCESS) && (status != CAPS_CLEAN)) goto herr;
    }

  /* check the value object after it has been updated */
  value = (capsValue *) vobj->blind;
  if (value->type != DoubleDeriv) {
    status = CAPS_BADTYPE;
    goto herr;
  }
  if (value->nderiv == 0) {
    status = CAPS_SUCCESS;
    goto herr;
  }

  status = EGADS_MALLOC;
  namex  = (char **) EG_alloc(value->nderiv*sizeof(char *));
  if (namex != NULL) {
    status = CAPS_SUCCESS;
    for (i = 0; i < value->nderiv; i++) namex[i] = value->derivs[i].name;
    args[0].num             = *nderiv = value->nderiv;
    args[0].members.strings = *names  = namex;
  }

herr:
  if (problem->dbFlag == 0)
    caps_jrnlWrite(CAPS_HASDERIV, problem, vobj, status, 1, args, problem->sNum,
                   problem->sNum);

  return status;
}


static int
caps_registerGIN(capsProblem *problem, const char *fullname)
{
  int        i, j, i1, i2, len, len_wrt, index, irow, icol;
  char       name[MAX_NAME_LEN];
  capsValue  *value, *vout;
  capsRegGIN *regGIN;
  capsDeriv  *derivs;

  len  = strlen(fullname);
  irow = icol = 0;
  for (i = 0; i < len; i++) {
    name[i] = fullname[i];
    if (fullname[i] == '[') break;
  }
  name[i] = '\0';
  i1 = i2 = 0;
  if (i != len) {
    for (j = i+1; j < len; j++)
      if (fullname[j] == ',') break;
    if (j == len) {
      sscanf(&fullname[i+1], "%d", &i1);
    } else {
      sscanf(&fullname[i+1], "%d,%d", &i1, &i2);
    }
  }

  for (index = 1; index <= problem->nGeomIn; index++)
    if (strcmp(problem->geomIn[index-1]->name,name) == 0) break;
  if (index > problem->nGeomIn) {
    printf(" CAPS Error: Object Not Found: %s (caps_getDeriv)!\n", name);
    return CAPS_NOTFOUND;
  }
  value = (capsValue *) problem->geomIn[index-1]->blind;
  if (value == NULL) {
    printf(" CAPS Error: Object: %s has NULL pointer (caps_getDeriv)!\n", name);
    return CAPS_NULLOBJ;
  }

  if (i1 > 0) {
    if ((value->nrow == 1) && (value->ncol == 1)) {
      printf(" CAPS Error: Object: %s has no index (caps_getDeriv)!\n", name);
      return CAPS_BADINDEX;
    }
    if (((i1  > 0) && (i2 <= 0)) ||
        ((i1 <= 0) && (i2  > 0))) {
      printf(" CAPS Error: Object: both indices in %s must be positive (caps_getDeriv)!\n",
             fullname);
      return CAPS_BADINDEX;
    }
    if (i2 != 0) {
      irow = i1;
      icol = i2;
    } else {
      if (value->nrow == 1) {
        irow = 1;
        icol = i1;
      } else {
        irow = i1;
        icol = 1;
      }
    }
    if ((irow < 1) || (irow > value->nrow)) {
      printf(" CAPS Error: Object: %s irow = %d [1-%d] (caps_getDeriv)!\n",
             name, irow, value->nrow);
      return CAPS_BADINDEX;
    }
    if ((icol < 1) || (icol > value->ncol)) {
      printf(" CAPS Error: Object: %s icol = %d [1-%d] (caps_getDeriv)!\n",
             name, icol, value->ncol);
      return CAPS_BADINDEX;
    }
  }

  /* are we already in the list? */
  for (i = 0; i < problem->nRegGIN; i++)
    if ((problem->regGIN[i].index == index) &&
        (problem->regGIN[i].irow  == irow)  &&
        (problem->regGIN[i].icol  == icol)) return CAPS_SUCCESS;

  /* make room */
  len = problem->nRegGIN + 1;
  regGIN = (capsRegGIN *) EG_reall(problem->regGIN, len*sizeof(capsRegGIN));
  if (regGIN == NULL) {
    printf(" CAPS Error: Cant ReAlloc registry for %s (caps_getDeriv)!\n",
           fullname);
    return EGADS_MALLOC;
  }
  problem->regGIN = regGIN;

  len_wrt = 1;
  if (irow == 0 && icol == 0) len_wrt = value->length;

  for (i = 0; i < problem->nGeomOut; i++) {
    if (problem->geomOut[i]              == NULL)      continue;
    if (problem->geomOut[i]->magicnumber != CAPSMAGIC) continue;
    if (problem->geomOut[i]->blind       == NULL)      continue;
    vout = (capsValue *) problem->geomOut[i]->blind;

    derivs = (capsDeriv *) EG_reall(vout->derivs, len*sizeof(capsDeriv));
    if (derivs == NULL) {
      printf(" CAPS Error: Cant Malloc dots for %s (caps_getDeriv)!\n",
             problem->geomOut[i]->name);
      return EGADS_MALLOC;
    }
    vout->derivs = derivs;

    vout->derivs[len-1].name    = NULL;
    vout->derivs[len-1].len_wrt = len_wrt;
    vout->derivs[len-1].deriv   = NULL;
    vout->nderiv                = len;
  }
  for (i = 0; i < problem->nGeomOut; i++) {
    if (problem->geomOut[i]              == NULL)      continue;
    if (problem->geomOut[i]->magicnumber != CAPSMAGIC) continue;
    if (problem->geomOut[i]->blind       == NULL)      continue;
    vout = (capsValue *) problem->geomOut[i]->blind;
    vout->derivs[len-1].name = EG_strdup(fullname);
  }

  problem->regGIN[len-1].name  = EG_strdup(fullname);
  problem->regGIN[len-1].index = index;
  problem->regGIN[len-1].irow  = irow;
  problem->regGIN[len-1].icol  = icol;
  problem->nRegGIN             = len;

  return CAPS_SUCCESS;
}


int
caps_getDeriv(capsObject *vobj, const char *name, int *len, int *len_wrt,
              double **deriv, int *nErr, capsErrs **errors)
{
  int         i, ret, ipmtr, irow, icol, status, nbody, buildTo, builtTo;
  int         irs, ire, ics, ice, index, outLevel;
  CAPSLONG    sNum;
  capsValue   *valueOut, *value_wrt;
  capsProblem *problem;
  capsObject  *pobject;
  capsJrnl    args[5];
  size_t      length;
  modl_T      *MODL;

  if (nErr              == NULL)         return CAPS_NULLVALUE;
  if (errors            == NULL)         return CAPS_NULLVALUE;
  if (len               == NULL)         return CAPS_NULLVALUE;
  if (len_wrt           == NULL)         return CAPS_NULLVALUE;
  if (deriv             == NULL)         return CAPS_NULLVALUE;
  *nErr   = 0;
  *errors = NULL;
  *len    = *len_wrt = 0;
  *deriv    = NULL;
  if (vobj              == NULL)         return CAPS_NULLOBJ;
  if (vobj->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (vobj->type        != VALUE)        return CAPS_BADTYPE;
  if (vobj->blind       == NULL)         return CAPS_NULLBLIND;
  if (name == NULL)                      return CAPS_NULLNAME;
  status = caps_findProblem(vobj, CAPS_GETDERIV, &pobject);
  if (status            != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;

  args[0].type = jInteger;
  args[1].type = jInteger;
  args[2].type = jPointer;
  args[3].type = jInteger;
  args[4].type = jErr;
  if (problem->dbFlag == 0) {
    status     = caps_jrnlRead(CAPS_GETDERIV, problem, vobj, 5, args, &sNum,
                               &ret);
    if (status == CAPS_JOURNALERR) return status;
    if (status == CAPS_JOURNAL) {
      *len     =            args[0].members.integer;
      *len_wrt =            args[1].members.integer;
      *deriv   = (double *) args[2].members.pointer;
      *nErr    =            args[3].members.integer;
      *errors  =            args[4].members.errs;
      return ret;
    }
  }

  sNum   = problem->sNum;

  /* do we need to update our value? */
  if (problem->dbFlag == 0)
    if (vobj->subtype == ANALYSISOUT) {
      status = caps_updateAnalysisOut(vobj, CAPS_GETDERIV, nErr, errors);
      if (status != CAPS_SUCCESS) return status;
    } else if (vobj->subtype == GEOMETRYOUT) {
      /* make sure geometry is up-to-date */
      status = caps_build(pobject, nErr, errors);
      if ((status != CAPS_SUCCESS) && (status != CAPS_CLEAN)) goto finis;

      status = caps_registerGIN(problem, name);
      if (status != CAPS_SUCCESS) goto finis;
    }

  /* check the type after updating */
  valueOut = (capsValue *) vobj->blind;
  if (valueOut->type != DoubleDeriv) {
    status = CAPS_BADTYPE;
    goto finis;
  }

  for (i = 0; i < valueOut->nderiv; i++) {
    if (strcmp(name, valueOut->derivs[i].name) != 0) continue;
    if ((valueOut->derivs[i].deriv == NULL) && (vobj->subtype == GEOMETRYOUT)) {
      index = problem->regGIN[i].index;
      irow  = problem->regGIN[i].irow;
      icol  = problem->regGIN[i].icol;

      if (problem->geomIn[index-1] == NULL) return CAPS_NULLOBJ;
      value_wrt = (capsValue *) problem->geomIn[index-1]->blind;
      ipmtr = value_wrt->pIndex;

      if (irow == 0 && icol == 0) {
        irs = ics = 1;
        ire = value_wrt->nrow; ice = value_wrt->ncol;
      } else {
        irs = irow; ics = icol;
        ire = irow; ice = icol;
      }

      for (irow = irs; irow <= ire; irow++) {
        for (icol = ics; icol <= ice; icol++) {
          /* clear all then set */
          status = ocsmSetDtime(problem->modl, 0);
          if (status != SUCCESS) goto finis;
          status = ocsmSetVelD(problem->modl, 0,     0,    0,    0.0);
          if (status != SUCCESS) goto finis;
          status = ocsmSetVelD(problem->modl, ipmtr, irow, icol, 1.0);
          if (status != SUCCESS) goto finis;
          buildTo = 0;
          nbody   = 0;
          outLevel = ocsmSetOutLevel(0);
          if (problem->outLevel > 0)
            printf(" CAPS Info: Building sensitivity information for: %s[%d,%d]\n",
                   problem->geomIn[index-1]->name, irow, icol);
          status  = ocsmBuild(problem->modl, buildTo, &builtTo, &nbody, NULL);
          ocsmSetOutLevel(outLevel);
    /*    printf(" CAPS Info: parameter %d %d %d sensitivity status = %d\n",
                 open, irow, icol, status);  */
          fflush(stdout);
          if (status != SUCCESS) goto finis;
          caps_geomOutSensit(problem, ipmtr, irow, icol);
          MODL = (modl_T*)problem->modl;
          if (MODL->dtime != 0 && problem->outLevel > 0)
            printf(" CAPS Info: Sensitivity finite difference used for: %s[%d,%d]\n",
                   problem->geomIn[index-1]->name, irow, icol);
        }
      }
    }
    *len     = valueOut->length;
    *len_wrt = valueOut->derivs[i].len_wrt;
    *deriv   = valueOut->derivs[i].deriv;
    status = CAPS_SUCCESS;
    goto finis;
  }
  status = CAPS_NOTFOUND;

finis:
  if (problem->dbFlag == 1) return status;
  length  = *len;
  length *= *len_wrt*sizeof(double);
  args[0].members.integer = *len;
  args[1].members.integer = *len_wrt;
  args[2].members.pointer = *deriv;
  args[2].length          = length;
  args[3].members.integer = *nErr;
  args[4].members.errs    = *errors;
  caps_jrnlWrite(CAPS_GETDERIV, problem, vobj, status, 5, args, sNum,
                 problem->sNum);

  return status;
}
