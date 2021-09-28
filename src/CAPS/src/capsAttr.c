/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Attribute Functions
 *
 *      Copyright 2014-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "capsBase.h"

extern /*@null@*/ /*@only@*/ char *EG_strdup(/*@null@*/ const char *str);

extern int  caps_writeObject(capsObject *cobj);
extern void caps_jrnlWrite(capsProblem *problem, capsObject *obj, int status,
                           int nargs, capsJrnl *args, CAPSLONG sNum0,
                           CAPSLONG sNum);
extern int  caps_jrnlRead(capsProblem *problem, capsObject *obj, int nargs,
                          capsJrnl *args, CAPSLONG *sNum, int *status);

int
caps_attrByName(capsObject *cobj, char *name, capsObject **attr)
{
  int            i, status, ret, atype, len;
  CAPSLONG       sNum;
  const void     *data;
  enum capsvType type;
  capsValue      *value;
  capsObject     *object, *pobject;
  capsProblem    *problem;
  capsJrnl       args[1];

  *attr = NULL;
  if (cobj              == NULL)         return CAPS_NULLOBJ;
  if (name              == NULL)         return CAPS_NULLNAME;
  if (cobj->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (cobj->attrs       == NULL)         return CAPS_NOTFOUND;
  status = caps_findProblem(cobj, CAPS_ATTRBYNAME, &pobject);
  if (status            != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  
  for (i = 0; i < cobj->attrs->nattrs; i++)
    if (strcmp(name, cobj->attrs->attrs[i].name) == 0) break;
  if (i == cobj->attrs->nattrs) return CAPS_NOTFOUND;
  
  status = caps_makeVal(Integer, 1, &i, &value);
  if (status != CAPS_SUCCESS) return status;
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(value);
    return status;
  }
  object->type    = VALUE;
  object->subtype = USER;
  object->parent  = pobject;
  object->blind   = value;
  
  args[0].type        = jValObj;
  args[0].members.obj = object;
  status = caps_jrnlRead(problem, cobj, 1, args, &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) {
      object->name = EG_strdup(name);
      *attr = object;
    } else {
      caps_delete(object);
    }
    return ret;
  }
  EG_free(object->blind);
  object->blind = NULL;

  sNum  = problem->sNum;
  len   = cobj->attrs->attrs[i].length;
  atype = cobj->attrs->attrs[i].type;
  if (atype == ATTRINT) {
    type = Integer;
    if (len == 1) {
      data = &cobj->attrs->attrs[i].vals.integer;
    } else {
      data =  cobj->attrs->attrs[i].vals.integers;
    }
  } else if (atype == ATTRREAL) {
    type = Double;
    if (len == 1) {
      data = &cobj->attrs->attrs[i].vals.real;
    } else {
      data =  cobj->attrs->attrs[i].vals.reals;
    }
  } else {
    type = String;
    data = cobj->attrs->attrs[i].vals.string;
  }
  status = caps_makeVal(type, len, data, &value);
  if (status != CAPS_SUCCESS) {
    caps_delete(object);
    goto adone;
  }
  if ((len == 1) || (type == String)) {
    value->nrow = 1;
    value->ncol = 1;
    value->dim  = Scalar;
  } else {
    value->nrow = len;
    value->ncol = 1;
    value->dim  = Vector;
  }
  object->name  = EG_strdup(name);
  object->blind = value;
  
  *attr  = object;
  status = CAPS_SUCCESS;

adone:
  args[0].members.obj = *attr;
  caps_jrnlWrite(problem, cobj, status, 1, args, sNum, problem->sNum);
  
  return status;
}


int
caps_attrByIndex(capsObject *cobj, int in, capsObject **attr)
{
  int            status, ret, atype, len;
  CAPSLONG       sNum;
  const void     *data;
  enum capsvType type;
  capsValue      *value;
  capsObject     *object, *pobject;
  capsProblem    *problem;
  capsJrnl       args[1];
  
  *attr = NULL;
  if (cobj              == NULL)              return CAPS_NULLOBJ;
  if (cobj->magicnumber != CAPSMAGIC)         return CAPS_BADOBJECT;
  if (cobj->attrs       == NULL)              return CAPS_NOTFOUND;
  if ((in < 1) || (in > cobj->attrs->nattrs)) return CAPS_BADINDEX;
  status = caps_findProblem(cobj, CAPS_ATTRBYINDEX, &pobject);
  if (status            != CAPS_SUCCESS)      return status;
  problem = (capsProblem *) pobject->blind;
  
  ret    = 0;
  status = caps_makeVal(Integer, 1, &ret, &value);
  if (status != CAPS_SUCCESS) return status;
  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(value);
    return status;
  }
  object->type    = VALUE;
  object->subtype = USER;
  object->parent  = pobject;
  object->blind   = value;
  
  args[0].type        = jValObj;
  args[0].members.obj = object;
  status = caps_jrnlRead(problem, cobj, 1, args, &sNum, &ret);
  if (status == CAPS_JOURNALERR) return status;
  if (status == CAPS_JOURNAL) {
    if (ret == CAPS_SUCCESS) {
      object->name = EG_strdup(cobj->attrs->attrs[in-1].name);
      *attr = object;
    } else {
      caps_delete(object);
    }
    return ret;
  }
  EG_free(object->blind);
  object->blind = NULL;

  len   = cobj->attrs->attrs[in-1].length;
  atype = cobj->attrs->attrs[in-1].type;
  if (atype == ATTRINT) {
    type = Integer;
    if (len == 1) {
      data = &cobj->attrs->attrs[in-1].vals.integer;
    } else {
      data =  cobj->attrs->attrs[in-1].vals.integers;
    }
  } else if (atype == ATTRREAL) {
    type = Double;
    if (len == 1) {
      data = &cobj->attrs->attrs[in-1].vals.real;
    } else {
      data =  cobj->attrs->attrs[in-1].vals.reals;
    }
  } else {
    type = String;
    data = cobj->attrs->attrs[in-1].vals.string;
  }
  status = caps_makeVal(type, len, data, &value);
  if (status != CAPS_SUCCESS) {
    caps_delete(object);
    goto idone;
  }
  if ((len == 1) || (type == String)) {
    value->nrow = 1;
    value->ncol = 1;
    value->dim  = Scalar;
  } else {
    value->nrow = len;
    value->ncol = 1;
    value->dim  = Vector;
  }
  object->name  = EG_strdup(cobj->attrs->attrs[in-1].name);
  object->blind = value;
  
  *attr  = object;
  status = CAPS_SUCCESS;

idone:
  args[0].members.obj = *attr;
  caps_jrnlWrite(problem, cobj, status, 1, args, sNum, problem->sNum);
  
  return status;
}


int
caps_setAttr(capsObject *cobj, /*@null@*/ const char *aname, capsObject *aval)
{
  int          atype, i, status, slen, find = -1;
  const int    *ints  = NULL;
  const double *reals = NULL;
  const char   *name;
  egAttr       *attr;
  egAttrs      *attrs;
  capsObject   *pobject;
  capsProblem  *problem;
  capsValue    *value;
  
  if (cobj              == NULL)         return CAPS_NULLOBJ;
  if (cobj->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (aval              == NULL)         return CAPS_NULLVALUE;
  if (aval->type        != VALUE)        return CAPS_BADTYPE;
  if (aval->blind       == NULL)         return CAPS_NULLBLIND;
  name = aname;
  if (name == NULL) name = aval->name;
  if (name              == NULL)         return CAPS_NULLNAME;
  status = caps_findProblem(cobj, CAPS_SETATTR, &pobject);
  if (status            != CAPS_SUCCESS) return status;
  value   = aval->blind;
  problem = (capsProblem *) pobject->blind;
  
  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4)               return CAPS_SUCCESS;
  
  if (value->type == Integer) {
    if (value->dim == 2) return CAPS_BADRANK;
    atype = ATTRINT;
    if (value->length == 1) {
      ints = &value->vals.integer;
    } else {
      ints =  value->vals.integers;
    }
  } else if (value->type == Double) {
    if (value->dim == 2) return CAPS_BADRANK;
    atype = ATTRREAL;
    if (value->length == 1) {
      reals = &value->vals.real;
    } else {
      reals =  value->vals.reals;
    }
  } else if (value->type == String){
    atype = ATTRSTRING;
  } else {
    return CAPS_BADVALUE;
  }

  attrs = (egAttrs *) cobj->attrs;
  
  if (attrs != NULL)
    for (i = 0; i < attrs->nattrs; i++)
      if (strcmp(attrs->attrs[i].name, name) == 0) {
        find = i;
        break;
      }
  
  if ((find != -1) && (attrs != NULL)) {
    
    /* an existing attribute -- reset the values */
    
    if (attrs->attrs[find].type == ATTRINT) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.integers);
    } else if (attrs->attrs[find].type == ATTRREAL) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.reals);
    } else {
      EG_free(attrs->attrs[find].vals.string);
    }
    
  } else {
    
    if (attrs == NULL) {
      attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
      if (attrs == NULL) return EGADS_MALLOC;
      attrs->nattrs = 0;
      attrs->attrs  = NULL;
      attrs->nseqs  = 0;
      attrs->seqs   = NULL;
      cobj->attrs   = attrs;
    }
    if (attrs->attrs == NULL) {
      attr = (egAttr *) EG_alloc((attrs->nattrs+1)*sizeof(egAttr));
    } else {
      attr = (egAttr *) EG_reall(attrs->attrs,
                                 (attrs->nattrs+1)*sizeof(egAttr));
    }
    if (attr == NULL) return EGADS_MALLOC;
    attrs->attrs = attr;
    find = attrs->nattrs;
    attrs->attrs[find].vals.string = NULL;
    attrs->attrs[find].name        = EG_strdup(name);
    if (attrs->attrs[find].name == NULL) return EGADS_MALLOC;
    attrs->nattrs += 1;
  }
  
  attrs->attrs[find].type   = atype;
  attrs->attrs[find].length = value->length;
  if (atype == ATTRINT) {
    if (ints != NULL) {
      if (value->length == 1) {
        attrs->attrs[find].vals.integer = ints[0];
      } else {
        attrs->attrs[find].vals.integers = (int *)
                                           EG_alloc(value->length*sizeof(int));
        if (attrs->attrs[find].vals.integers == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < value->length; i++)
            attrs->attrs[find].vals.integers[i] = ints[i];
        }
      }
    }
  } else if (atype == ATTRREAL) {
    if (reals != NULL) {
      if (value->length == 1) {
        attrs->attrs[find].vals.real = reals[0];
      } else {
        attrs->attrs[find].vals.reals = (double *)
                                        EG_alloc(value->length*sizeof(double));
        if (attrs->attrs[find].vals.reals == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < value->length; i++)
            attrs->attrs[find].vals.reals[i] = reals[i];
        }
      }
    }
  } else {
    for (slen = i = 0; i < value->length; i++)
      slen += strlen(value->vals.string + slen) + 1;
    EG_free(attrs->attrs[find].vals.string);
    attrs->attrs[find].vals.string = (char *) EG_alloc(slen*sizeof(char));
    if (attrs->attrs[find].vals.string == NULL) return EGADS_MALLOC;
    for (i = 0; i < slen; i++)
      attrs->attrs[find].vals.string[i] = value->vals.string[i];
  }
  
  status = caps_writeObject(cobj);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeObject = %d (caps_setAttr)\n", status);
  
  return CAPS_SUCCESS;
}


int
caps_deleteAttr(capsObject *cobj, /*@null@*/ char *name)
{
  int         i, status, find = -1;
  egAttrs     *attrs;
  capsObject  *pobject;
  capsProblem *problem;

  if (cobj              == NULL)         return CAPS_NULLOBJ;
  if (cobj->magicnumber != CAPSMAGIC)    return CAPS_BADOBJECT;
  if (cobj->attrs       == NULL)         return CAPS_NOTFOUND;
  status = caps_findProblem(cobj, CAPS_DELETEATTR, &pobject);
  if (status            != CAPS_SUCCESS) return status;
  problem = (capsProblem *) pobject->blind;
  
  /* ignore if restarting */
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag == 4)               return CAPS_SUCCESS;
  
  /* delete all? */
  if (name == NULL) {
    caps_freeAttrs(&cobj->attrs);
    status = caps_writeObject(cobj);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_writeObject = %d (caps_deleteAttr)\n", status);
    return CAPS_SUCCESS;
  }
  
  /* find the attribute */
  attrs  = cobj->attrs;
  for (i = 0; attrs->nattrs; i++)
    if (strcmp(name, attrs->attrs[i].name) == 0) {
      find = i;
      break;
    }
  
  if (find == -1) return CAPS_NOTFOUND;
 
  /* remove it! */
  EG_free(attrs->attrs[find].name);
  if (attrs->attrs[find].type == ATTRINT) {
    if (attrs->attrs[find].length > 1)
      EG_free(attrs->attrs[find].vals.integers);
  } else if (attrs->attrs[find].type == ATTRREAL) {
    if (attrs->attrs[find].length > 1)
      EG_free(attrs->attrs[find].vals.reals);
  } else {
    EG_free(attrs->attrs[find].vals.string);
  }
  for (i = find+1; i < attrs->nattrs; i++)
    attrs->attrs[i-1] = attrs->attrs[i];
  attrs->nattrs -= 1;
  
  status = caps_writeObject(cobj);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeObject = %d (caps_deleteAttr)\n", status);
  
  return CAPS_SUCCESS;
}
