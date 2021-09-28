/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Problem Object Functions
 *
 *      Copyright 2014-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#ifndef WIN32
#include <strings.h>
#include <unistd.h>
#include <limits.h>
#else
#define strcasecmp  stricmp
#define snprintf   _snprintf
#define getcwd     _getcwd
#define PATH_MAX   _MAX_PATH
#endif

#include "udunits.h"

#include "capsBase.h"
#include "capsAIM.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"
#include "udp.h"

static int    CAPSnLock   = 0;
static char **CAPSlocks   = NULL;
static int    CAPSextSgnl = 1;

extern /*@null@*/
       void *caps_initUnits();
extern void  caps_initFunIDs();
extern int   caps_freeError(/*@only@*/ capsErrs *errs);
extern int   caps_Aprx1DFree(/*@only@*/ capsAprx1D *approx);
extern int   caps_Aprx2DFree(/*@only@*/ capsAprx2D *approx);
extern int   caps_checkDiscr(capsDiscr *discr, int l, char *line);
extern int   caps_filter(capsProblem *problem, capsAnalysis *analysis);
extern int   caps_statFile(const char *path);
extern int   caps_rmFile(const char *path);
extern int   caps_rmDir(const char *path);
extern int   caps_mkDir(const char *path);
extern int   caps_cpDir(const char *src, const char *dst);
extern int   caps_rename(const char *src, const char *dst);


void
caps_rmLock(void)
{
  int i;
  
  for (i = 0; i < CAPSnLock; i++)
    if (CAPSlocks[i] != NULL)
      caps_rmFile(CAPSlocks[i]);

  for (i = 0; i < CAPSnLock; i++)
    if (CAPSlocks[i] != NULL)
     EG_free(CAPSlocks[i]);
  
  EG_free(CAPSlocks);
  CAPSnLock = 0;
  CAPSlocks = NULL;
  
  /* cleanup udp storage */
  ocsmFree(NULL);
}


void
caps_rmLockOnClose(const char* root)
{
  int i, j;

  for (i = 0; i < CAPSnLock; i++) {
    if ((CAPSlocks[i] != NULL) &&
        (strncmp(CAPSlocks[i], root, strlen(root)) == 0)) {
      caps_rmFile(CAPSlocks[i]);
      EG_free(CAPSlocks[i]);
      for (j = i; j < CAPSnLock-1; j++)
        CAPSlocks[j] = CAPSlocks[j+1];
      CAPSnLock--;
      break;
    }
  }
}


void
caps_externSignal()
{
  if (CAPSextSgnl == 1) {
    atexit(caps_rmLock);
    caps_initFunIDs();
  }

  CAPSextSgnl = 0;
}


static void
caps_intHandler(int sig)
{
  caps_rmLock();
  (void) signal(sig, SIG_DFL);
  raise(sig);
}


static void
caps_initSignals()
{
  if (CAPSextSgnl <= 0) return;
  atexit(caps_rmLock);
  (void) signal(SIGSEGV, caps_intHandler);
  (void) signal(SIGINT,  caps_intHandler);
#ifdef WIN32
  (void) signal(SIGABRT, caps_intHandler);
#else
/*@-unrecog@*/
  (void) signal(SIGHUP,  caps_intHandler);
  (void) signal(SIGBUS,  caps_intHandler);
/*@+unrecog@*/
#endif
  caps_initFunIDs();

  CAPSextSgnl = -1;
}


int
caps_isNameOK(const char *name)
{
  int i, n;
  
  if (name == NULL) return CAPS_NULLNAME;
  
  n = strlen(name);
  for (i = 0; i < n; i++)
    if ((name[i] == '/') || (name[i] == '\\') || (name[i] == ':') ||
        (name[i] == ' ') || (name[i] ==  '@') || (name[i] == '!') ||
        (name[i] == '(') || (name[i] ==  ')')) return CAPS_BADNAME;
  
  return CAPS_SUCCESS;
}


static void
caps_prunePath(char *path)
{
  int i, j, k, len, hit;
  
  /* first remove any "./" */
  do {
    hit = 0;
    len = strlen(path);
    for (i = 1; i < len-2; i++)
      if ( (path[i-1] != '.') && (path[i] == '.') &&
          ((path[i+1] == '/') || (path[i+1] == '\\'))) break;
    if (i == len-2) continue;
    hit++;
    for (j = i; j < len-1; j++) path[j] = path[j+2];
  } while (hit != 0);
  
  /* remove the levels */
  do {
    hit = 0;
    len = strlen(path);
    for (i = 0; i < len-3; i++)
      if ((path[i] == '.') && (path[i+1] == '.') &&
          ((path[i+2] == '/') || (path[i+2] == '\\'))) break;
    if (i == len-3) return;
    hit++;
    while (i+3*hit < len-3)
      if ((path[i+3*hit] == '.') && (path[i+1+3*hit] == '.') &&
          ((path[i+2+3*hit] == '/') || (path[i+2+3*hit] == '\\'))) {
        hit++;
      } else {
        break;
      }
    k = i + 3*hit;
    hit++;
    
    /* backup */
    for (j = i; j >= 0; j--)
      if ((path[j] == '/') || (path[j] == '\\')) {
        hit--;
        if (hit == 0) break;
      }
    if (j < 0) {
      printf(" CAPS Warning: Invalid path = %s\n", path);
      return;
    }
    for (i = k-1; i <= len; j++, i++) path[j] = path[i];
    hit++;
  
  } while (hit != 0);
  
}


void
caps_freeValue(capsValue *value)
{
  int i;

  if (value == NULL) return;

  if (value->units      != NULL) EG_free(value->units);
  if (value->meshWriter != NULL) EG_free(value->meshWriter);
  if ((value->type == Boolean) || (value->type == Integer)) {
    if (value->length > 1) EG_free(value->vals.integers);
  } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
    if (value->length > 1) EG_free(value->vals.reals);
  } else if (value->type == String) {
    EG_free(value->vals.string);
  } else if (value->type == Tuple) {
    caps_freeTuple(value->length, value->vals.tuple);
  }
  if (value->partial != NULL) EG_free(value->partial);
  if (value->derivs == NULL) return;
  for (i = 0; i < value->nderiv; i++) {
    if (value->derivs[i].name != NULL) EG_free(value->derivs[i].name);
    if (value->derivs[i].deriv  != NULL) EG_free(value->derivs[i].deriv);
  }
  EG_free(value->derivs);
}


void
caps_freeFList(capsObject *obj)
{
  int       i;
  capsFList *flist, *next;
  
  if (obj->flist == NULL) return;
  
  flist = obj->flist;
  do {
    if (flist->type == jTuple) {
      caps_freeTuple(flist->num, flist->member.tuple);
    } else if (flist->type == jPointer) {
      EG_free(flist->member.pointer);
    } else if (flist->type == jOwn) {
      caps_freeOwner(&flist->member.own);
    } else if (flist->type == jOwns) {
      for (i = 0; i < flist->num; i++) caps_freeOwner(&flist->member.owns[i]);
      EG_free(flist->member.owns);
    } else if (flist->type == jStrings) {
      for (i = 0; i < flist->num; i++) EG_free(flist->member.strings[i]);
      EG_free(flist->member.strings);
    } else if (flist->type == jEgos) {
      EG_deleteObject(flist->member.model);
    } else {
      printf(" CAPS Internal: caps_freeFList type = %d\n", flist->type);
    }
    next  = flist->next;
    EG_free(flist);
    flist = next;
  } while (flist != NULL);
  
  obj->flist = NULL;
}


int
caps_hierarchy(capsObject *obj, char **full)
{
  int        i, n, pos, index, len = 0;
  char       *path, number[5];
  capsObject *object;
  capsBound  *bound;
  capsValue  *value;
  static char type[7] = {'U', 'P', 'V', 'A', 'B', 'S', 'D'};
  static char sub[11] = {'N', 'S', 'P', 'I', 'O', 'P', 'U', 'I', 'O', 'C', 'N'};
  
  *full = NULL;
  if (obj == NULL) return CAPS_SUCCESS;
  if (obj->type == PROBLEM) {
    path = (char *) EG_alloc(2*sizeof(char));
    if (path == NULL) return EGADS_MALLOC;
    path[0] = '.';
    path[1] =  0;
    *full = path;
    return CAPS_SUCCESS;
  }

  object = obj;
  do {
    if ((object->type == VALUE) || (object->type == BOUND)) {
      len += 8;
    } else {
      len += strlen(object->name) + 4;
    }
    object = object->parent;
  } while (object->type != PROBLEM);
  
  path = (char *) EG_alloc(len*sizeof(char));
  if (path == NULL) return EGADS_MALLOC;
  
  pos    = len;
  object = obj;
  do {
    if ((object->type == VALUE) || (object->type == BOUND)) {
      n  = 4;
    } else {
      n  = strlen(object->name);
    }
    pos -= n+4;
    path[pos  ] = type[object->type];
    path[pos+1] = sub[object->subtype];
    path[pos+2] = '-';
    if ((object->type == VALUE) || (object->type == BOUND)) {
      if (object->type == VALUE) {
        value = (capsValue *) object->blind;
        index = value->index;
      } else {
        bound = (capsBound *) object->blind;
        index = bound->index;
      }
      snprintf(number, 5, "%4.4d", index);
      for (i = 0; i < 4; i++) path[pos+i+3] = number[i];
    } else {
      for (i = 0; i < n; i++) path[pos+i+3] = object->name[i];
    }
    path[pos+n+3] = '/';
    object = object->parent;
  } while (object->type != PROBLEM);
  path[len-1] = 0;
  
  *full = path;
  return CAPS_SUCCESS;
}


static int
caps_string2obj(capsProblem *problem, /*@null@*/ const char *full,
                capsObject **object)
{
  int           i, j, in, it, is, pos, len, index;
  char          name[PATH_MAX];
  capsObject    *obj;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vs;
  static char   type[7] = {'U', 'P', 'V', 'A', 'B', 'S', 'D'};
  static char   sub[11] = {'N', 'S', 'P', 'I', 'O', 'P', 'U', 'I', 'O', 'C', 'N'};
  
  *object = NULL;
  if (full == NULL) return CAPS_SUCCESS;
  len = strlen(full);
  if ((len == 1) && (full[0] == '.')) {
    *object = problem->mySelf;
    return CAPS_SUCCESS;
  }
  
  obj = NULL;
  pos = 0;
  do {
    for (it = 0; it < 7; it++)
      if (full[pos] == type[it]) break;
    if (it == 7) {
      printf(" CAPS Error: type %c not found (caps_string2obj)\n", full[pos]);
      return CAPS_BADOBJECT;
    }
    pos++;
    for (is = 0; is < 11; is++)
      if (full[pos] == sub[is]) break;
    if (is == 12) {
      printf(" CAPS Error: subtype %c not found (caps_string2obj)\n", full[pos]);
      return CAPS_BADOBJECT;
    }
    pos++;
    if (full[pos] != '-') {
      printf(" CAPS Error: %c not a seperator (caps_string2obj)\n", full[pos]);
      return CAPS_BADOBJECT;
    }
    pos++;
    in = pos;
    for (; pos < len; pos++)
      if ((full[pos] == '/') || (full[pos] == '\\')) break;
    
    if ((it == 2) || (it == 4)) {
      
      /* look at index */
      index = atoi(&full[in]);
      if (obj == NULL) {
        /* in the problem */
        if (it == 2) {
          if (is == 2) {
            if (index > problem->nParam) {
              printf(" CAPS Error: Bad index %d for pValue (caps_string2obj)\n",
                     index);
              return CAPS_BADINDEX;
            }
            obj = problem->params[index-1];
          } else if (is == 3) {
            if (index > problem->nGeomIn) {
              printf(" CAPS Error: Bad index %d for giValue (caps_string2obj)\n",
                     index);
              return CAPS_BADINDEX;
            }
            obj = problem->geomIn[index-1];
          } else if (is == 4) {
            if (index > problem->nGeomOut) {
              printf(" CAPS Error: Bad index %d for goValue (caps_string2obj)\n",
                     index);
              return CAPS_BADINDEX;
            }
            obj = problem->geomOut[index-1];
          } else {
            printf(" CAPS Error: incorrect sub %d for Value (caps_string2obj)\n",
                   is);
            return CAPS_BADOBJECT;
          }
        } else {
          for (i = 0; i < problem->nBound; i++) {
            if (problem->bounds[i] == NULL) continue;
            bound = (capsBound *) problem->bounds[i]->blind;
            if (bound == NULL) continue;
            if (bound->index == index) {
              obj = problem->bounds[i];
              break;
            }
          }
          if (i == problem->nBound) {
            printf(" CAPS Error: Bad index %d for Bound (caps_string2obj)\n",
                   index);
            return CAPS_BADINDEX;
          }
        }
      } else {
        /* in an object */
        if (it == 4) {
          printf(" CAPS Error: Bad Bound child %s (caps_string2obj)\n", full);
          return CAPS_BADOBJECT;
        }
        if (obj->type != ANALYSIS) {
          printf(" CAPS Error: Bad Value child %s (caps_string2obj)\n", full);
          return CAPS_BADOBJECT;
        }
        analysis = (capsAnalysis *) obj->blind;
        if (analysis == NULL) {
          printf(" CAPS Error: NULL Analysis %s (caps_string2obj)\n", full);
          return CAPS_BADOBJECT;
        }
        if (is == 3) {
          if (index > analysis->nAnalysisIn) {
            printf(" CAPS Error: Bad index %d for aiValue (caps_string2obj)\n",
                   index);
            return CAPS_BADINDEX;
          }
          obj = analysis->analysisIn[index-1];
        } else if (is == 4) {
          if (index > analysis->nAnalysisOut) {
            printf(" CAPS Error: Bad index %d for aoValue (caps_string2obj)\n",
                   index);
            return CAPS_BADINDEX;
          }
          obj = analysis->analysisOut[index-1];
        } else {
          printf(" CAPS Error: Incorrect sub %d for Value (caps_string2obj)\n",
                 is);
          return CAPS_BADOBJECT;
        }
      }
      
    } else {
      
      /* name is from in to pos-1 */
      j = 0;
      for (i = in; i < pos; i++, j++) name[j] = full[i];
      name[j] = 0;
      if (obj == NULL) {
        /* in Problem */
        if (it != 3) {
          printf(" CAPS Error: Incorrect type %d for Problem (caps_string2obj)\n",
                 it);
          return CAPS_BADOBJECT;
        }
        for (i = 0; i < problem->nAnalysis; i++)
          if (strcmp(name, problem->analysis[i]->name) == 0) {
            obj = problem->analysis[i];
            break;
          }
        if (i == problem->nAnalysis) {
          printf(" CAPS Error: Analysis %s Not Found (caps_string2obj)\n", name);
          return CAPS_NOTFOUND;
        }
      } else {
        /* in Bound or VertexSet? */
        if (obj->type == BOUND) {
          if (it != 5) {
            printf(" CAPS Error: Bad Bound child %s (caps_string2obj)\n", full);
            return CAPS_BADOBJECT;
          }
          bound = (capsBound *) obj->blind;
          if (bound == NULL) {
            printf(" CAPS Error: NULL Bound %s (caps_string2obj)\n", full);
            return CAPS_BADOBJECT;
          }
          for (i = 0; i < bound->nVertexSet; i++)
            if (strcmp(name, bound->vertexSet[i]->name) == 0) {
              obj = bound->vertexSet[i];
              break;
            }
          if (i == bound->nVertexSet) {
            printf(" CAPS Error: VertexSet %s Not Found (caps_string2obj)\n",
                   name);
            return CAPS_NOTFOUND;
          }
        } else if (obj->type == VERTEXSET) {
          if (it != 6) {
            printf(" CAPS Error: Bad VertexSet child %s (caps_string2obj)\n",
                   full);
            return CAPS_BADOBJECT;
          }
          vs = (capsVertexSet *) obj->blind;
          if (vs == NULL) {
            printf(" CAPS Error: NULL VertexSet %s (caps_string2obj)\n", full);
            return CAPS_BADOBJECT;
          }
          for (i = 0; i < vs->nDataSets; i++)
            if (strcmp(name, vs->dataSets[i]->name) == 0) {
              obj = vs->dataSets[i];
              break;
            }
          if (i == vs->nDataSets) {
            printf(" CAPS Error: DataSet %s Not Found (caps_string2obj)\n",
                   name);
            return CAPS_NOTFOUND;
          }
        } else {
          printf(" CAPS Error: Incorrect type %d for child (caps_string2obj)\n",
                 it);
          return CAPS_BADOBJECT;
        }
      }
    }
    
    pos++;
  } while (pos < len);
  
  *object = obj;
  
  return CAPS_SUCCESS;
}


static int
caps_writeDoubles(FILE *fp, int len, /*@null@*/ const double *reals)
{
  size_t n;

  if (len < 0) return CAPS_BADINDEX;
  if ((reals == NULL) && (len != 0)) return CAPS_NULLVALUE;
  n = fwrite(&len, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if (reals == NULL) return CAPS_SUCCESS;
  n = fwrite(reals, sizeof(double), len, fp);
  if (n != len) return CAPS_IOERR;

  return CAPS_SUCCESS;
}


static int
caps_writeString(FILE *fp, /*@null@*/ const char *string)
{
  int    len = 0;
  size_t n;

  if (string != NULL) len = strlen(string) + 1;
  n = fwrite(&len, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if (string == NULL) return CAPS_SUCCESS;
  n = fwrite(string, sizeof(char), len, fp);
  if (n != len) return CAPS_IOERR;

  return CAPS_SUCCESS;
}


static int
caps_writeStrings(FILE *fp, int len, /*@null@*/ const char *string)
{
  int    slen = 0, i;
  size_t n;

  if (string != NULL)
    for (i = 0; i < len; i++)
      slen += strlen(string + slen) + 1;
  n = fwrite(&slen, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if (string == NULL) return CAPS_SUCCESS;
  n = fwrite(string, sizeof(char), slen, fp);
  if (n != slen) return CAPS_IOERR;

  return CAPS_SUCCESS;
}


static int
caps_writeTuple(FILE *fp, int len, enum capsNull nullVal,
                /*@null@*/ capsTuple *tuple)
{
  int i, stat;

  if (len < 0) return CAPS_BADINDEX;
  if ((tuple == NULL) && (len != 0) && (nullVal != IsNull))
    return CAPS_NULLVALUE;

  if (tuple == NULL) return CAPS_SUCCESS;
  for (i = 0; i < len; i++) {
    stat = caps_writeString(fp, tuple[i].name);
    if (stat != CAPS_SUCCESS) return stat;
    stat = caps_writeString(fp, tuple[i].value);
    if (stat != CAPS_SUCCESS) return stat;
  }

  return CAPS_SUCCESS;
}


static int
caps_writeOwn(FILE *fp, capsOwn writer, capsOwn own)
{
  int    i, stat, nLines;
  size_t n;

  nLines = 0;
  if (own.lines != NULL) nLines = own.nLines;
  n = fwrite(&nLines, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if (own.lines != NULL)
    for (i = 0; i < nLines; i++) {
      stat = caps_writeString(fp, own.lines[i]);
      if (stat != CAPS_SUCCESS) return stat;
    }

  if (own.pname == NULL) {
    stat = caps_writeString(fp, writer.pname);
  } else {
    stat = caps_writeString(fp, own.pname);
  }
  if (stat != CAPS_SUCCESS) return stat;
  if (own.pID == NULL) {
    stat = caps_writeString(fp, writer.pID);
  } else {
    stat = caps_writeString(fp, own.pID);
  }
  if (stat != CAPS_SUCCESS) return stat;
  if (own.user == NULL) {
    stat = caps_writeString(fp, writer.user);
  } else {
    stat = caps_writeString(fp, own.user);
  }
  if (stat != CAPS_SUCCESS) return stat;

  n = fwrite(own.datetime, sizeof(short), 6, fp);
  if (n != 6) return CAPS_IOERR;

  n = fwrite(&own.sNum, sizeof(CAPSLONG), 1, fp);
  if (n != 1) return CAPS_IOERR;

  return CAPS_SUCCESS;
}


static int
caps_writeHistory(FILE *fp, const capsObject *obj)
{
  int    i, j, nHistory, nLines, stat;
  size_t n;
  
  nHistory = obj->nHistory;
  if (obj->history == NULL) nHistory = 0;
  n = fwrite(&nHistory, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if (nHistory == 0) return CAPS_SUCCESS;
  
  if (obj->history != NULL)
    for (j = 0; j < nHistory; j++) {
      nLines = 0;
      if (obj->history[j].lines != NULL) nLines = obj->history[j].nLines;
      n = fwrite(&nLines, sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
      if (obj->history[j].lines != NULL)
        for (i = 0; i < nLines; i++) {
          stat = caps_writeString(fp, obj->history[j].lines[i]);
          if (stat != CAPS_SUCCESS) return stat;
        }

    stat = caps_writeString(fp, obj->history[j].pname);
    if (stat != CAPS_SUCCESS) return stat;
    stat = caps_writeString(fp, obj->history[j].pID);
    if (stat != CAPS_SUCCESS) return stat;
    stat = caps_writeString(fp, obj->history[j].user);
    if (stat != CAPS_SUCCESS) return stat;

    n = fwrite(obj->history[j].datetime, sizeof(short), 6, fp);
    if (n != 6) return CAPS_IOERR;

    n = fwrite(&obj->history[j].sNum, sizeof(CAPSLONG), 1, fp);
    if (n != 1) return CAPS_IOERR;
  }
  
  return CAPS_SUCCESS;
}


static int
caps_writeAttrs(FILE *fp, egAttrs *attrs)
{
  int    nattr, i, stat;
  size_t n;
  egAttr *attr;

  nattr = 0;
  if (attrs != NULL) nattr = attrs->nattrs;
  n = fwrite(&nattr, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if ((nattr == 0) || (attrs == NULL)) return CAPS_SUCCESS;

  attr = attrs->attrs;
  for (i = 0; i < nattr; i++) {
    n = fwrite(&attr[i].type,   sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    n = fwrite(&attr[i].length, sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    stat = caps_writeString(fp, attr[i].name);
    if (stat != CAPS_SUCCESS) return CAPS_IOERR;
    if (attr[i].type == ATTRINT) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = fwrite(&attr[i].vals.integer, sizeof(int),              1, fp);
      } else if (attr[i].length > 1) {
        n = fwrite(attr[i].vals.integers, sizeof(int), attr[i].length, fp);
      }
      if (n != attr[i].length) return CAPS_IOERR;
    } else if (attr[i].type == ATTRREAL) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = fwrite(&attr[i].vals.real, sizeof(double),              1, fp);
      } else if (attr[i].length > 1) {
        n = fwrite(attr[i].vals.reals, sizeof(double), attr[i].length, fp);
      }
      if (n != attr[i].length) return CAPS_IOERR;
    } else {
      stat = caps_writeStrings(fp, attr[i].length, attr[i].vals.string);
      if (stat != CAPS_SUCCESS) return CAPS_IOERR;
    }
  }

  return CAPS_SUCCESS;
}


static int
caps_writeValue(FILE *fp, capsOwn writer, capsObject *obj)
{
  int       i, j, stat;
  char      *name;
  size_t    n;
  capsValue *value;
  
  stat = caps_writeHistory(fp, obj);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeOwn(fp, writer, obj->last);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeAttrs(fp, obj->attrs);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeString(fp, obj->name);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_hierarchy(obj, &name);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, name);
  EG_free(name);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  
  value = (capsValue *) obj->blind;
  if (value == NULL) return CAPS_NULLVALUE;

  n = fwrite(&value->type,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->length,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->dim,     sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->nrow,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->ncol,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->lfixed,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->sfixed,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->nullVal, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->index,   sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->pIndex,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->gInType, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&value->nderiv,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if (value->type == Integer) {
    n = fwrite(value->limits.ilims, sizeof(int),    2, fp);
    if (n != 2) return CAPS_IOERR;
  } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
    n = fwrite(value->limits.dlims, sizeof(double), 2, fp);
    if (n != 2) return CAPS_IOERR;
  }

  stat = caps_writeString(fp, value->units);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, value->meshWriter);
  if (stat != CAPS_SUCCESS) return stat;
  name = NULL;
  if (value->link != NULL) {
    stat = caps_hierarchy(value->link, &name);
    if (stat != CAPS_SUCCESS) return stat;
  }
  stat = caps_writeString(fp, name);
  EG_free(name);
  if (stat != CAPS_SUCCESS) return stat;
  n = fwrite(&value->linkMethod, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if ((value->length == 1)     && (value->type != String) &&
      (value->type != Pointer) && (value->type != Tuple)  &&
      (value->type != PointerMesh)) {
    if ((value->type == Double) || (value->type == DoubleDeriv)) {
      n = fwrite(&value->vals.real, sizeof(double), 1, fp);
      if (n != 1) return CAPS_IOERR;
    } else {
      n = fwrite(&value->vals.integer, sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
    }
  } else {
    if ((value->type == Pointer) || (value->type == PointerMesh)) {
      /* what do we do? */
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      n = fwrite(value->vals.reals, sizeof(double), value->length, fp);
      if (n != value->length) return CAPS_IOERR;
    } else if (value->type == String) {
      stat = caps_writeStrings(fp, value->length, value->vals.string);
      if (stat != CAPS_SUCCESS) return stat;
    } else if (value->type == Tuple) {
      stat = caps_writeTuple(fp, value->length, value->nullVal,
                             value->vals.tuple);
      if (stat != CAPS_SUCCESS) return stat;
    } else {
      n = fwrite(value->vals.integers, sizeof(int), value->length, fp);
      if (n != value->length) return CAPS_IOERR;
    }
  }

  if (value->nullVal == IsPartial) {
    n = fwrite(value->partial, sizeof(int), value->length, fp);
    if (n != value->length) return CAPS_IOERR;
  }

  if (value->nderiv != 0) {
    for (i = 0; i < value->nderiv; i++) {
      stat = caps_writeString(fp, value->derivs[i].name);
      if (stat != CAPS_SUCCESS) return stat;
      n = fwrite(&value->derivs[i].rank,  sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
      j = value->length*value->derivs[i].rank;
      if (value->derivs[i].deriv == NULL) j = 0;
      n = fwrite(&j,                    sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
      if (j != 0) {
        n = fwrite(value->derivs[i].deriv, sizeof(double), j, fp);
        if (n != j) return CAPS_IOERR;
      }
    }
  }

  return CAPS_SUCCESS;
}


int
caps_writeValueObj(capsProblem *problem, capsObject *valobj)
{
  int  status;
  char filename[PATH_MAX], temp[PATH_MAX], *full;
  FILE *fp;
  
  status = caps_hierarchy(valobj, &full);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Warning: caps_hierarchy = %d\n", status);
    return status;
  }
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\%s", problem->root, full);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/%s",   problem->root, full);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",   problem->root);
#endif
  EG_free(full);

  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  status = caps_writeValue(fp, problem->writer, valobj);
  fclose(fp);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot write %s!\n", filename);
    return status;
  }
  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }
  
  return CAPS_SUCCESS;
}


int
caps_dumpGeomVals(capsProblem *problem, int flag)
{
  int  i, status;
  char filename[PATH_MAX], current[PATH_MAX], temp[PATH_MAX];
  FILE *fp;
  
#ifdef WIN32
  snprintf(current, PATH_MAX, "%s\\capsRestart", problem->root);
  snprintf(temp,    PATH_MAX, "%s\\xxTempxx",    current);
#else
  snprintf(current, PATH_MAX, "%s/capsRestart",  problem->root);
  snprintf(temp,    PATH_MAX, "%s/xxTempxx",     current);
#endif

  if ((flag == 0) || (flag == 1))
    for (i = 0; i < problem->nGeomIn; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VI-%4.4d", current, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/VI-%4.4d",  current, i+1);
#endif
      fp = fopen(temp, "wb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s!\n", filename);
        return CAPS_DIRERR;
      }
      status = caps_writeValue(fp, problem->writer, problem->geomIn[i]);
      fclose(fp);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot write %s!\n", filename);
        return status;
      }
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot rename %s!\n", filename);
        return status;
      }
    }
  
  if ((flag == 0) || (flag == 2))
    for (i = 0; i < problem->nGeomOut; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VO-%4.4d", current, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/VO-%4.4d",  current, i+1);
#endif
      fp = fopen(temp, "wb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s!\n", filename);
        return CAPS_DIRERR;
      }
      status = caps_writeValue(fp, problem->writer, problem->geomOut[i]);
      fclose(fp);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot write %s!\n", filename);
        return status;
      }
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot rename %s!\n", filename);
        return status;
      }
    }
  
  return CAPS_SUCCESS;
}


static int
caps_writeAnalysis(FILE *fp, capsProblem *problem, capsObject *aobject)
{
  int          i, stat;
  capsAnalysis *analysis;
  size_t       n;
  
  analysis = (capsAnalysis *) aobject->blind;
  
  stat = caps_writeHistory(fp, aobject);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeOwn(fp, problem->writer, aobject->last);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeOwn(fp, problem->writer, analysis->pre);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeAttrs(fp, aobject->attrs);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_writeString(fp, aobject->name);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, analysis->loadName);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, analysis->path);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, analysis->unitSys);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_writeString(fp, analysis->intents);
  if (stat != CAPS_SUCCESS) return stat;
  n = fwrite(&analysis->major,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&analysis->minor,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&analysis->autoexec, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fwrite(&analysis->nField,   sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  for (i = 0; i < analysis->nField; i++) {
    stat = caps_writeString(fp, analysis->fields[i]);
    if (stat != CAPS_SUCCESS) return stat;
  }
  if (analysis->nField != 0) {
    n = fwrite(analysis->ranks,  sizeof(int), analysis->nField, fp);
    if (n != analysis->nField) return CAPS_IOERR;
    n = fwrite(analysis->fInOut, sizeof(int), analysis->nField, fp);
    if (n != analysis->nField) return CAPS_IOERR;
  }
  
  return CAPS_SUCCESS;
}


static int
caps_writeAnalysisObj(capsProblem *problem, capsObject *aobject)
{
  int  status;
  char filename[PATH_MAX], current[PATH_MAX], temp[PATH_MAX];
  FILE *fp;

#ifdef WIN32
  snprintf(current,  PATH_MAX, "%s\\capsRestart\\AN-%s",
           problem->root, aobject->name);
  snprintf(temp,     PATH_MAX, "%s\\xxTempxx", current);
  snprintf(filename, PATH_MAX, "%s\\analysis", current);
#else
  snprintf(current,  PATH_MAX, "%s/capsRestart/AN-%s",
           problem->root, aobject->name);
  snprintf(temp,     PATH_MAX, "%s/xxTempxx",  current);
  snprintf(filename, PATH_MAX, "%s/analysis",  current);
#endif
  
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  status = caps_writeAnalysis(fp, problem, aobject);
  fclose(fp);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot write Analysis %s!\n", filename);
    return status;
  }
  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }
  
  return CAPS_SUCCESS;
}
  

int
caps_dumpAnalysis(capsProblem *problem, capsObject *aobject)
{
  int          i, status;
  char         filename[PATH_MAX], current[PATH_MAX], temp[PATH_MAX];
  capsAnalysis *analysis;
  FILE         *fp;
  
  analysis = (capsAnalysis *) aobject->blind;
#ifdef WIN32
  snprintf(current,  PATH_MAX, "%s\\capsRestart\\AN-%s",
           problem->root, aobject->name);
  snprintf(temp,     PATH_MAX, "%s\\xxTempxx", current);
  snprintf(filename, PATH_MAX, "%s\\analysis", current);
#else
  snprintf(current,  PATH_MAX, "%s/capsRestart/AN-%s",
           problem->root, aobject->name);
  snprintf(temp,     PATH_MAX, "%s/xxTempxx",  current);
  snprintf(filename, PATH_MAX, "%s/analysis",  current);
#endif
  
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  status = caps_writeAnalysis(fp, problem, aobject);
  fclose(fp);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot write Analysis %s!\n", filename);
    return status;
  }
  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }
  
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\analysis.txt", current);
  snprintf(temp,     PATH_MAX, "%s\\xxTempxx",     current);
#else
  snprintf(filename, PATH_MAX, "%s/analysis.txt",  current);
  snprintf(temp,     PATH_MAX, "%s/xxTempxx",      current);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  fprintf(fp, "%d %d\n", analysis->nAnalysisIn, analysis->nAnalysisOut);
  if (analysis->analysisIn != NULL)
    for (i = 0; i < analysis->nAnalysisIn; i++)
      fprintf(fp, "%s\n", analysis->analysisIn[i]->name);
  if (analysis->analysisOut != NULL)
    for (i = 0; i < analysis->nAnalysisOut; i++)
      fprintf(fp, "%s\n", analysis->analysisOut[i]->name);
  fclose(fp);
  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }
 
  if (analysis->analysisIn != NULL)
    for (i = 0; i < analysis->nAnalysisIn; i++) {
  #ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VI-%4.4d", current, i+1);
  #else
      snprintf(filename, PATH_MAX, "%s/VI-%4.4d",  current, i+1);
  #endif
      fp = fopen(temp, "wb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s!\n", filename);
        return CAPS_DIRERR;
      }
      status = caps_writeValue(fp, problem->writer, analysis->analysisIn[i]);
      fclose(fp);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot write %s!\n", filename);
        return status;
      }
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot rename %s!\n", filename);
        return status;
      }
    }
  
  if (analysis->analysisOut != NULL)
    for (i = 0; i < analysis->nAnalysisOut; i++) {
  #ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VO-%4.4d", current, i+1);
  #else
      snprintf(filename, PATH_MAX, "%s/VO-%4.4d",  current, i+1);
  #endif
      fp = fopen(temp, "wb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s!\n", filename);
        return CAPS_DIRERR;
      }
      status = caps_writeValue(fp, problem->writer, analysis->analysisOut[i]);
      fclose(fp);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot write %d %s!\n", status, filename);
        return status;
      }
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Cannot rename %s!\n", filename);
        return status;
      }
    }
  
  return CAPS_SUCCESS;
}


int
caps_writeDataSet(capsObject *dobject)
{
  int         stat;
  char        filename[PATH_MAX], temp[PATH_MAX], *full, *name=NULL;
  capsObject  *pobject;
  capsProblem *problem;
  capsDataSet *ds;
  size_t      n;
  FILE        *fp;

  stat = caps_findProblem(dobject, 9999, &pobject);
  if (stat != CAPS_SUCCESS) return stat;
  
  ds      = (capsDataSet *) dobject->blind;
  problem = (capsProblem *) pobject->blind;
  stat    = caps_hierarchy(dobject, &full);
  if (stat != CAPS_SUCCESS) return stat;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\%s", problem->root, full);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/%s",   problem->root, full);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",   problem->root);
#endif
  EG_free(full);
  
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  
  stat = caps_writeHistory(fp, dobject);
  if (stat != CAPS_SUCCESS) goto wrterrs;
  stat = caps_writeOwn(fp, problem->writer, dobject->last);
  if (stat != CAPS_SUCCESS) goto wrterrs;
  stat = caps_writeAttrs(fp, dobject->attrs);
  if (stat != CAPS_SUCCESS) goto wrterrs;
  stat = caps_writeString(fp, dobject->name);
  if (stat != CAPS_SUCCESS) goto wrterrs;

  n = fwrite(&ds->ftype,  sizeof(int), 1, fp);
  if (n != 1) goto wrterrs;
  n = fwrite(&ds->npts,   sizeof(int), 1, fp);
  if (n != 1) goto wrterrs;
  n = fwrite(&ds->rank,   sizeof(int), 1, fp);
  if (n != 1) goto wrterrs;
  stat = caps_writeString(fp, ds->units);
  if (stat != CAPS_SUCCESS) goto wrterrs;

  stat = caps_writeDoubles(fp, ds->npts*ds->rank, ds->data);
  if (stat != CAPS_SUCCESS) goto wrterrs;

  name = NULL;
  if (ds->link != NULL) {
    stat = caps_hierarchy(ds->link, &name);
    if (stat != CAPS_SUCCESS) return stat;
  }
  stat = caps_writeString(fp, name);
  EG_free(name);
  if (stat != CAPS_SUCCESS) return stat;
  n = fwrite(&ds->linkMethod,  sizeof(int), 1, fp);
  if (n != 1) goto wrterrs;
  fclose(fp);

  stat = caps_rename(temp, filename);
  if (stat != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return stat;
  }

  return CAPS_SUCCESS;
  
wrterrs:
  fclose(fp);
  return CAPS_IOERR;
}


int
caps_dumpBound(capsObject *pobject, capsObject *bobject)
{
  int           i, j, status;
  char          filename[PATH_MAX], temp[PATH_MAX], *full;
  capsProblem   *problem;
  capsBound     *bound;
  capsVertexSet *vs;
  FILE          *fp;
  
  problem = (capsProblem *) pobject->blind;
  bound   = (capsBound *)   bobject->blind;
  
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\vSets.txt",
           problem->root, bound->index);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\xxTempxx",
           problem->root, bound->index);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d/vSets.txt",
           problem->root, bound->index);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/BN-%4.4d/xxTempxx",
           problem->root, bound->index);
#endif
  fp = fopen(temp, "w");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  fprintf(fp, "%d\n", bound->nVertexSet);
  if (bound->vertexSet != NULL)
    for (i = 0; i < bound->nVertexSet; i++)
      fprintf(fp, "%d %s\n", bound->vertexSet[i]->subtype,
                             bound->vertexSet[i]->name);
  fclose(fp);
  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }
  
  if (bound->vertexSet != NULL)
    for (i = 0; i < bound->nVertexSet; i++) {
      status = caps_hierarchy(bound->vertexSet[i], &full);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Warning: caps_hierarchy = %d\n", status);
        continue;
      }
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\%s", problem->root, full);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/%s",   problem->root, full);
#endif
      status = caps_mkDir(filename);
      if (status != CAPS_SUCCESS) {
        printf(" CAPS Warning: Cannot mkdir %s (caps_dumpBound)\n", filename);
        EG_free(full);
        continue;
      }
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\%s\\dSets.txt",
               problem->root, full);
      snprintf(temp,     PATH_MAX, "%s\\capsRestart\\%s\\xxTempxx",
               problem->root, full);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/%s/dSets.txt",
               problem->root, full);
      snprintf(temp,     PATH_MAX, "%s/capsRestart/%s/xxTempxx",
               problem->root, full);
#endif
      EG_free(full);
      vs = (capsVertexSet *) bound->vertexSet[i]->blind;
      if (vs == NULL) continue;
      
      fp = fopen(temp, "w");
      if (fp == NULL) {
        printf(" CAPS Warning: Cannot open %s!\n", filename);
        continue;
      }
      fprintf(fp, "%d\n", vs->nDataSets);
      if (vs->dataSets != NULL)
        for (j = 0; j < vs->nDataSets; j++) {
          fprintf(fp, "%s\n", vs->dataSets[j]->name);
          status = caps_writeDataSet(vs->dataSets[j]);
          if (status != CAPS_SUCCESS) {
            printf(" CAPS Warning: Writing %s = %d\n",
                   vs->dataSets[j]->name, status);
            continue;
          }
        }
      fclose(fp);
      status = caps_rename(temp, filename);
      if (status != CAPS_SUCCESS)
        printf(" CAPS Warning: Cannot rename %s!\n", filename);
    }

  return CAPS_SUCCESS;
}


int
caps_writeProblem(const capsObject *pobject)
{
  int         i, ivec[2], stat;
  char        filename[PATH_MAX], temp[PATH_MAX];
  capsProblem *problem;
  size_t      n;
  FILE        *fp;

  problem = (capsProblem *) pobject->blind;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\Problem",  problem->root);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/Problem",    problem->root);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",   problem->root);
#endif
  
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  ivec[0] = CAPSMAJOR;
  ivec[1] = CAPSMINOR;

  n = fwrite(ivec,               sizeof(int),      2, fp);
  if (n != 2) goto writerr;
  n = fwrite(&pobject->subtype,  sizeof(int),      1, fp);
  if (n != 1) goto writerr;
  stat = caps_writeHistory(fp, pobject);
  if (stat != CAPS_SUCCESS) goto writerr;
  stat = caps_writeOwn(fp, problem->writer, pobject->last);
  if (stat != CAPS_SUCCESS) goto writerr;
  stat = caps_writeAttrs(fp, pobject->attrs);
  if (stat != CAPS_SUCCESS) goto writerr;
  stat = caps_writeString(fp, pobject->name);
  if (stat != CAPS_SUCCESS) goto writerr;
  stat = caps_writeString(fp, problem->phName);
  if (stat != CAPS_SUCCESS) goto writerr;
  stat = caps_writeOwn(fp, problem->writer, problem->geometry);
  if (stat != CAPS_SUCCESS) goto writerr;
  n = fwrite(&problem->sNum,      sizeof(CAPSLONG), 1, fp);
  if (n != 1) goto writerr;
#ifdef WIN32
  n = fwrite(&problem->jpos,      sizeof(__int64),  1, fp);
#else
  n = fwrite(&problem->jpos,      sizeof(long),     1, fp);
#endif
  if (n != 1) goto writerr;
  n = fwrite(&problem->outLevel,  sizeof(int),      1, fp);
  if (n != 1) goto writerr;
  n = fwrite(&problem->nEGADSmdl, sizeof(int),      1, fp);
  if (n != 1) goto writerr;
  i = problem->nRegGIN;
  if (problem->regGIN == NULL) i = 0;
  n = fwrite(&i,                  sizeof(int),      1, fp);
  if (n != 1) goto writerr;

  if (problem->regGIN != NULL)
    for (i = 0; i < problem->nRegGIN; i++) {
      stat = caps_writeString(fp, problem->regGIN[i].name);
      if (stat != CAPS_SUCCESS) goto writerr;
      n = fwrite(&problem->regGIN[i].index, sizeof(int), 1, fp);
      if (n != 1) goto writerr;
      n = fwrite(&problem->regGIN[i].irow,  sizeof(int), 1, fp);
      if (n != 1) goto writerr;
      n = fwrite(&problem->regGIN[i].icol,  sizeof(int), 1, fp);
      if (n != 1) goto writerr;
    }
  
  fclose(fp);
  stat = caps_rename(temp, filename);
#ifdef WIN32
  if (stat != CAPS_SUCCESS) {
    Sleep(100);                   /* VMware filesystem sync issue? */
    stat = caps_rename(temp, filename);
  }
#endif
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: Cannot rename %s!\n", filename);
  
  return CAPS_SUCCESS;
  
writerr:
  fclose(fp);
  return CAPS_IOERR;
}


int
caps_writeVertexSet(capsObject *vobject)
{
  int           i, dim, status;
  char          filename[PATH_MAX], temp[PATH_MAX], *full;
  capsObject    *pobject;
  capsProblem   *problem;
  capsVertexSet *vs;
  size_t        n;
  FILE          *fp;

  status = caps_findProblem(vobject, 9999, &pobject);
  if (status != CAPS_SUCCESS) return status;
  
  vs      = (capsVertexSet *) vobject->blind;
  problem = (capsProblem *)   pobject->blind;
  status  = caps_hierarchy(vobject, &full);
  if (status != CAPS_SUCCESS) return status;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\%s\\vs", problem->root, full);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\xxTempxx", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/%s/vs",    problem->root, full);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/xxTempxx",   problem->root);
#endif
  EG_free(full);
  
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }

  status = caps_writeHistory(fp, vobject);
  if (status != CAPS_SUCCESS) goto wrterr;
  status = caps_writeOwn(fp, problem->writer, vobject->last);
  if (status != CAPS_SUCCESS) goto wrterr;
  status = caps_writeAttrs(fp, vobject->attrs);
  if (status != CAPS_SUCCESS) goto wrterr;
  status = caps_writeString(fp, vobject->name);
  if (status != CAPS_SUCCESS) goto wrterr;

  /* write needed discr info */
  dim = 0;
  if (vs->discr != NULL) dim = vs->discr->dim;
  n = fwrite(&dim, sizeof(int), 1, fp);
  if (n != 1) goto wrterr;
  if (vs->analysis != NULL) dim = 0;

  if (dim != 0) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
    n = fwrite(&vs->discr->nVerts, sizeof(int), 1, fp);
    if (n != 1) goto wrterr;
    for (i = 0; i < vs->discr->nVerts; i++) {
      n = fwrite(&vs->discr->verts[3*i], sizeof(double), 3, fp);
      if (n != 3) goto wrterr;
    }
#endif
/*@+nullderef@*/
  }
  fclose(fp);

  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }

  return CAPS_SUCCESS;
  
wrterr:
  fclose(fp);
  return CAPS_IOERR;
}


int
caps_writeBound(const capsObject *bobject)
{
  int         i, status;
  char        filename[PATH_MAX], temp[PATH_MAX];
  capsObject  *pobject;
  capsProblem *problem;
  capsBound   *bound;
  size_t      n;
  FILE        *fp;
  
  bound   = (capsBound *)   bobject->blind;
  pobject = bobject->parent;
  problem = (capsProblem *) pobject->blind;

#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\bound",
           problem->root, bound->index);
  snprintf(temp,     PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\xxTempxx",
           problem->root, bound->index);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d/bound",
           problem->root, bound->index);
  snprintf(temp,     PATH_MAX, "%s/capsRestart/BN-%4.4d/xxTempxx",
           problem->root, bound->index);
#endif
  fp = fopen(temp, "wb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  
  status = caps_writeHistory(fp, bobject);
  if (status != CAPS_SUCCESS) goto werror;
  status = caps_writeOwn(fp, problem->writer, bobject->last);
  if (status != CAPS_SUCCESS) goto werror;
  status = caps_writeAttrs(fp, bobject->attrs);
  if (status != CAPS_SUCCESS) goto werror;
  status = caps_writeString(fp, bobject->name);
  if (status != CAPS_SUCCESS) goto werror;
  n = fwrite(&bound->dim,     sizeof(int),    1, fp);
  if (n != 1) goto werror;
  n = fwrite(&bound->state,   sizeof(int),    1, fp);
  if (n != 1) goto werror;
  n = fwrite(&bound->plimits, sizeof(double), 4, fp);
  if (n != 4) goto werror;
  n = fwrite(&bound->iBody,   sizeof(int),    1, fp);
  if (n != 1) goto werror;
  n = fwrite(&bound->iEnt,    sizeof(int),    1, fp);
  if (n != 1) goto werror;
  n = fwrite(&bound->index,   sizeof(int),    1, fp);
  if (n != 1) goto werror;

  i = 0;
  if (bound->curve != NULL) i = bound->curve->nrank;
  n = fwrite(&i, sizeof(int), 1, fp);
  if (n != 1) goto werror;
  if (bound->curve != NULL) {
    n = fwrite(&bound->curve->periodic, sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = fwrite(&bound->curve->nts,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = 2*bound->curve->nts*bound->curve->nrank;
    status = caps_writeDoubles(fp, n, bound->curve->interp);
    if (status != CAPS_SUCCESS) goto werror;
    n = fwrite(bound->curve->trange,    sizeof(double), 2, fp);
    if (n != 2) goto werror;
    n = fwrite(&bound->curve->ntm,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = 2*bound->curve->ntm;
    status = caps_writeDoubles(fp, n, bound->curve->tmap);
    if (status != CAPS_SUCCESS) goto werror;
  }

  i = 0;
  if (bound->surface != NULL) i = bound->surface->nrank;
  n = fwrite(&i, sizeof(int), 1, fp);
  if (n != 1) goto werror;
  if (bound->surface != NULL) {
    n = fwrite(&bound->surface->periodic, sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = fwrite(&bound->surface->nus,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = fwrite(&bound->surface->nvs,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = 4*bound->surface->nus*bound->surface->nvs*bound->surface->nrank;
    status = caps_writeDoubles(fp, n, bound->surface->interp);
    if (status != CAPS_SUCCESS) goto werror;
    n = fwrite(bound->surface->urange,    sizeof(double), 2, fp);
    if (n != 2) goto werror;
    n = fwrite(bound->surface->vrange,    sizeof(double), 2, fp);
    if (n != 2) goto werror;
    n = fwrite(&bound->surface->num,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = fwrite(&bound->surface->nvm,      sizeof(int),    1, fp);
    if (n != 1) goto werror;
    n = 8*bound->surface->num*bound->surface->nvm;
    status = caps_writeDoubles(fp, n, bound->surface->uvmap);
    if (status != CAPS_SUCCESS) goto werror;
  }
  fclose(fp);

  status = caps_rename(temp, filename);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Cannot rename %s!\n", filename);
    return status;
  }

  return CAPS_SUCCESS;
  
werror:
  fclose(fp);
  return CAPS_IOERR;
}


int
caps_writeObject(capsObject *object)
{
  int         status = CAPS_SUCCESS;
  capsObject  *pobject;
  capsProblem *problem;
  
  if (object->type == PROBLEM) {
    status = caps_writeProblem(object);
  } else if (object->type == VALUE) {
    status = caps_findProblem(object, 9999, &pobject);
    if (status != CAPS_SUCCESS) return status;
    problem = (capsProblem *) pobject->blind;
    status  = caps_writeValueObj(problem, object);
  } else if (object->type == ANALYSIS) {
    status = caps_findProblem(object, 9999, &pobject);
    if (status != CAPS_SUCCESS) return status;
    problem = (capsProblem *) pobject->blind;
    status  = caps_writeAnalysisObj(problem, object);
  } else if (object->type == BOUND) {
    status = caps_writeBound(object);
  } else if (object->type == VERTEXSET) {
    status = caps_writeVertexSet(object);
  } else if (object->type == DATASET) {
    status = caps_writeDataSet(object);
  }
  
  return status;
}


static int
caps_writeErrs(FILE *fp, capsErrs *errs)
{
  int    i, j, nError, stat;
  char   *full;
  size_t n;
  
  nError = 0;
  if (errs != NULL) nError = errs->nError;
  n = fwrite(&nError, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if ((nError == 0) || (errs == NULL)) return CAPS_SUCCESS;
  
  for (i = 0; i < nError; i++) {
    full = NULL;
    if (errs->errors[i].errObj != NULL) {
      stat = caps_hierarchy(errs->errors[i].errObj, &full);
      if (stat != CAPS_SUCCESS)
        printf(" CAPS Warning: caps_hierarchy = %d (caps_writeErrs)\n", stat);
    }
    stat = caps_writeString(fp, full);
    if (stat != CAPS_SUCCESS) {
      printf(" CAPS Warning: caps_writeString = %d (caps_writeErrs)\n", stat);
      return CAPS_IOERR;
    }
    n = fwrite(&errs->errors[i].eType,  sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    n = fwrite(&errs->errors[i].index,  sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    n = fwrite(&errs->errors[i].nLines, sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    for (j = 0; j < errs->errors[i].nLines; j++) {
      stat = caps_writeString(fp, errs->errors[i].lines[j]);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Warning: %d caps_writeString = %d (caps_writeErrs)\n",
               j, stat);
        return CAPS_IOERR;
      }
    }
  }
  
  return CAPS_SUCCESS;
}


void
caps_jrnlWrite(capsProblem *problem, capsObject *obj, int status,
               int nargs, capsJrnl *args, CAPSLONG sNum0, CAPSLONG sNum)
{
  int       i, j, stat;
  char      filename[PATH_MAX], *full;
  capsFList *flist;
  size_t    n;

  if (problem->jrnl   == NULL)            return;
  if (problem->stFlag == CAPS_JOURNALERR) return;
  
  n = fwrite(&problem->funID, sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jwrterr;
  n = fwrite(&sNum0,          sizeof(CAPSLONG), 1, problem->jrnl);
  if (n != 1) goto jwrterr;
  n = fwrite(&status,         sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jwrterr;
/*
  printf(" *** Writing from FunID = %s  status = %d ***\n",
         caps_funID[problem->funID], status);
*/
  if (status == CAPS_SUCCESS) {
    n = fwrite(&obj->last.sNum, sizeof(CAPSLONG), 1, problem->jrnl);
    if (n != 1) goto jwrterr;
    /* cleanup? */
    if (obj->flist != NULL) {
      flist = (capsFList *) obj->flist;
      if (obj->last.sNum > flist->sNum) caps_freeFList(obj);
    }
    
    for (i = 0; i < nargs; i++)
      switch (args[i].type) {
          
        case jInteger:
          n = fwrite(&args[i].members.integer, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          break;
          
        case jDouble:
          n = fwrite(&args[i].members.real, sizeof(double), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          break;
          
        case jString:
          stat = caps_writeString(problem->jrnl, args[i].members.string);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_writeString = %d!\n", stat);
            goto jwrterr;
          }
          break;
          
        case jStrings:
          n = fwrite(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          for (j = 0; j < args[i].num; j++) {
            stat = caps_writeString(problem->jrnl, args[i].members.strings[j]);
            if (stat != CAPS_SUCCESS) {
              printf(" CAPS Warning: Journal %d caps_writeString = %d!\n",
                     j, stat);
              goto jwrterr;
            }
          }
          break;
          
        case jTuple:
          n = fwrite(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          stat = caps_writeTuple(problem->jrnl, args[i].num, NotNull,
                                 args[i].members.tuple);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_writeTuple = %d!\n", stat);
            goto jwrterr;
          }
          break;
          
        case jPointer:
        case jPtrFree:
          n = fwrite(&args[i].length, sizeof(size_t), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          if (args[i].length != 0) {
            n = fwrite(args[i].members.pointer, sizeof(char), args[i].length,
                       problem->jrnl);
            if (n != args[i].length) goto jwrterr;
          }
          break;
          
        case jObject:
          stat = caps_hierarchy(args[i].members.obj, &full);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_hierarchy = %d!\n", stat);
            goto jwrterr;
          }
          stat = caps_writeString(problem->jrnl, full);
          EG_free(full);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Jrnl caps_writeString Obj = %d!\n", stat);
            goto jwrterr;
          }
          break;
          
        case jObjs:
          n = fwrite(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          if (args[i].num != 0) {
            for (j = 0; j < args[i].num; j++) {
              stat = caps_hierarchy(args[i].members.objs[j], &full);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Journal caps_hierarchy = %d!\n", stat);
                goto jwrterr;
              }
              stat = caps_writeString(problem->jrnl, full);
              EG_free(full);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Jrnl caps_writeString Obj = %d!\n", stat);
                goto jwrterr;
              }
            }
          }
          break;
          
        case jValObj:
          stat = caps_writeValue(problem->jrnl, problem->writer,
                                 args[i].members.obj);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Jrnl caps_writeValue = %d!\n", stat);
            goto jwrterr;
          }
          break;
          
        case jErr:
          stat = caps_writeErrs(problem->jrnl, args[i].members.errs);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_writeErrs = %d!\n", stat);
            goto jwrterr;
          }
          break;
          
        case jOwn:
          stat = caps_writeOwn(problem->jrnl, problem->writer,
                               args[i].members.own);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_Own = %d!\n", stat);
            goto jwrterr;
          }
          break;

        case jOwns:
          n = fwrite(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jwrterr;
          if (args[i].num != 0) {
            for (j = 0; j < args[i].num; j++) {
              stat = caps_writeOwn(problem->jrnl, problem->writer,
                                   args[i].members.owns[j]);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Journal caps_Owns %d = %d!\n", j, stat);
                goto jwrterr;
              }
            }
          }
          break;

        case jEgos:
          if (args[i].members.model == NULL) {
            j = -1;
            n = fwrite(&j, sizeof(int), 1, problem->jrnl);
            if (n != 1) goto jwrterr;
          } else {
            n = fwrite(&problem->nEGADSmdl, sizeof(int), 1, problem->jrnl);
            if (n != 1) goto jwrterr;
#ifdef WIN32
            snprintf(filename, PATH_MAX, "%s\\capsRestart\\model%4.4d.egads",
                     problem->root, problem->nEGADSmdl);
#else
            snprintf(filename, PATH_MAX, "%s/capsRestart/model%4.4d.egads",
                     problem->root, problem->nEGADSmdl);
#endif
            stat = EG_saveModel(args[i].members.model, filename);
            if (stat != CAPS_SUCCESS)
              printf(" CAPS Warning: EG_saveModel = %d (caps_jrnlWrite)!\n",
                     stat);
            if (args[i].members.model->oclass == MODEL)
             EG_deleteObject(args[i].members.model);
            problem->nEGADSmdl++;
          }
        }
      }
  
  n = fwrite(&sNum,           sizeof(CAPSLONG), 1, problem->jrnl);
  if (n != 1) goto jwrterr;
  n = fwrite(&problem->funID, sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jwrterr;
  if (status == CAPS_SUCCESS) {
#ifdef WIN32
    problem->jpos = _ftelli64(problem->jrnl);
#else
    problem->jpos = ftell(problem->jrnl);
#endif
  }
  fflush(problem->jrnl);
  
  stat = caps_writeProblem(problem->mySelf);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeProblem = %d (caps_jrnlWrite)!\n", stat);
  
  return;
  
jwrterr:
  printf(" CAPS ERROR: Writing Journal File -- disabled (funID = %d)\n",
         problem->funID);
  fclose(problem->jrnl);
  problem->jrnl = NULL;
}


static void
caps_sizeCB(void *modl, int ipmtr, int nrow, int ncol)
{
  int         index, j, k, n, status;
  double      *reals, dot;
  modl_T      *MODL;
  capsObject  *object;
  capsValue   *value;
  capsProblem *problem;

  MODL    = (modl_T *) modl;
  problem = (capsProblem *) MODL->userdata;
  if ((problem->funID != CAPS_SETVALUE) &&
      (problem->funID != CAPS_READPARAMETERS)) {
    printf(" CAPS Internal: caps_sizeCB called from funID = %d!\n",
           problem->funID);
    return;
  }
  if (problem->funID == CAPS_READPARAMETERS) return;

  /* find our GeometryIn index */
  for (index = 0; index < problem->nGeomIn; index++) {
    object = problem->geomIn[index];
    if (object        == NULL) continue;
    if (object->blind == NULL) continue;
    value  = (capsValue *) object->blind;
    if (value->pIndex == ipmtr) break;
  }
  if (index == problem->nGeomIn) {
    printf(" CAPS Warning: cant find ocsm ipmtr = %d (caps_sizeCB)!\n", ipmtr);
    return;
  }
  object = problem->geomIn[index];
  value  = (capsValue *) object->blind;
/*
  printf(" ###### caps_sizeCB = %s  %d -> %d   %d -> %d ######\n",
         object->name, value->nrow, nrow, value->ncol, ncol);
 */

  reals = NULL;
  if (nrow*ncol != 1) {
    reals = (double *) EG_alloc(nrow*ncol*sizeof(double));
    if (reals == NULL) {
      printf(" CAPS Warning: %s resize %d %d Malloc(caps_sizeCB)\n",
             object->name, nrow, ncol);
      return;
    }
  }
  if (value->length != 1) EG_free(value->vals.reals);
  value->length = nrow*ncol;
  value->nrow   = nrow;
  value->ncol   = ncol;
  if (value->length != 1) value->vals.reals = reals;

  reals = value->vals.reals;
  if (value->length == 1) reals = &value->vals.real;
  for (n = k = 0; k < nrow; k++)
    for (j = 0; j < ncol; j++, n++)
      ocsmGetValu(problem->modl, ipmtr, k+1, j+1, &reals[n], &dot);

  value->dim = Scalar;
  if ((ncol > 1) && (nrow > 1)) {
    value->dim = Array2D;
  } else if ((ncol > 1) || (nrow > 1)) {
    value->dim = Vector;
  }

  caps_freeOwner(&object->last);
  object->last.sNum = problem->sNum;
  caps_fillDateTime(object->last.datetime);
  status = caps_writeValueObj(problem, object);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeValueObj = %d\n", status);
}


static int
caps_readString(FILE *fp, char **string)
{
  int    len;
  size_t n;

  *string = NULL;
  n = fread(&len, sizeof(int), 1, fp);
  if (n   != 1) return CAPS_IOERR;
  if (len <  0) return CAPS_IOERR;
  if (len == 0) return CAPS_SUCCESS;

  *string = (char *) EG_alloc(len*sizeof(char));
  if (*string == NULL) return EGADS_MALLOC;

  n = fread(*string, sizeof(char), len, fp);
  if (n != len) {
    EG_free(*string);
    *string = NULL;
    return CAPS_IOERR;
  }

  return CAPS_SUCCESS;
}


static int
caps_readStrings(FILE *fp, int len, char **string)
{
  int    slen;
  size_t n;

  *string = NULL;
  n = fread(&slen, sizeof(int), 1, fp);
  if (n   != 1) return CAPS_IOERR;
  if (len <  0) return CAPS_IOERR;
  if (len == 0) return CAPS_SUCCESS;

  *string = (char *) EG_alloc(slen*sizeof(char));
  if (*string == NULL) return EGADS_MALLOC;

  n = fread(*string, sizeof(char), slen, fp);
  if (n != slen) {
    EG_free(*string);
    *string = NULL;
    return CAPS_IOERR;
  }

  return CAPS_SUCCESS;
}


static int
caps_readDoubles(FILE *fp, int *len, double **reals)
{
  size_t n;

  *reals = NULL;
  n = fread(len, sizeof(int), 1, fp);
  if (n    != 1) return CAPS_IOERR;
  if (*len <  0) return CAPS_IOERR;
  if (*len == 0) return CAPS_SUCCESS;

  *reals = (double *) EG_alloc(*len*sizeof(double));
  if (*reals == NULL) return EGADS_MALLOC;

  n = fread(*reals, sizeof(double), *len, fp);
  if (n != *len) {
    EG_free(*reals);
    *reals = NULL;
    return CAPS_IOERR;
  }

  return CAPS_SUCCESS;
}


static int
caps_readTuple(FILE *fp, int len, enum capsNull nullVal, capsTuple **tuple)
{
  int       i, stat;
  capsTuple *tmp;

  *tuple = NULL;
  if (nullVal == IsNull) return CAPS_SUCCESS;
  stat = caps_makeTuple(len, &tmp);
  if (stat != CAPS_SUCCESS) return stat;

  for (i = 0; i < len; i++) {
    stat = caps_readString(fp, &tmp[i].name);
    if (stat != CAPS_SUCCESS) {
      caps_freeTuple(len, tmp);
      return stat;
    }
    stat = caps_readString(fp, &tmp[i].value);
    if (stat != CAPS_SUCCESS) {
      caps_freeTuple(len, tmp);
      return stat;
    }
  }
  *tuple = tmp;

  return CAPS_SUCCESS;
}


static int
caps_readOwn(FILE *fp, capsOwn *own)
{
  int    i, stat;
  size_t n;

  own->nLines = 0;
  own->lines  = NULL;
  own->pname  = NULL;
  own->pID    = NULL;
  own->user   = NULL;
  
  n = fread(&own->nLines, sizeof(int), 1, fp);
  if (n != 1) {
    caps_freeOwner(own);
    return CAPS_IOERR;
  }
  if (own->nLines != 0) {
    own->lines = (char **) EG_alloc(own->nLines*sizeof(char *));
    if (own->lines == NULL) {
      own->nLines = 0;
      caps_freeOwner(own);
      return EGADS_MALLOC;
    }
    for (i = 0; i < own->nLines; i++) {
      stat = caps_readString(fp, &own->lines[i]);
      if (stat != CAPS_SUCCESS) return stat;
    }
  }

  stat = caps_readString(fp, &own->pname);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_readString(fp, &own->pID);
  if (stat != CAPS_SUCCESS) {
    caps_freeOwner(own);
    return stat;
  }
  stat = caps_readString(fp, &own->user);
  if (stat != CAPS_SUCCESS) {
    caps_freeOwner(own);
    return stat;
  }

  n = fread(own->datetime, sizeof(short), 6, fp);
  if (n != 6) {
    caps_freeOwner(own);
    return CAPS_IOERR;
  }

  n = fread(&own->sNum, sizeof(CAPSLONG), 1, fp);
  if (n != 1) {
    caps_freeOwner(own);
    return CAPS_IOERR;
  }

  return CAPS_SUCCESS;
}


static int
caps_readHistory(FILE *fp, capsObject *obj)
{
  int    i, j, stat;
  size_t n;
  
  n = fread(&obj->nHistory, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if (obj->nHistory == 0) return CAPS_SUCCESS;
  
  obj->history = (capsOwn *) EG_alloc(obj->nHistory*sizeof(capsOwn));
  if (obj->history == NULL) {
    obj->nHistory = 0;
    return EGADS_MALLOC;
  }
  
  for (j = 0; j < obj->nHistory; j++) {
    obj->history[j].nLines = 0;
    obj->history[j].lines  = NULL;
    obj->history[j].pname  = NULL;
    obj->history[j].pID    = NULL;
    obj->history[j].user   = NULL;
  }
  
  for (j = 0; j < obj->nHistory; j++) {
    n = fread(&obj->history[j].nLines, sizeof(int), 1, fp);
    if (n != 1) return CAPS_IOERR;
    if (obj->history[j].nLines != 0) {
      obj->history[j].lines = (char **)
                                EG_alloc(obj->history[j].nLines*sizeof(char *));
      if (obj->history[j].lines == NULL) {
        obj->history[j].nLines = 0;
        return EGADS_MALLOC;
      }
      for (i = 0; i < obj->history[j].nLines; i++) {
        stat = caps_readString(fp, &obj->history[j].lines[i]);
        if (stat != CAPS_SUCCESS) return stat;
      }
    }

    stat = caps_readString(fp, &obj->history[j].pname);
    if (stat != CAPS_SUCCESS) return stat;
    stat = caps_readString(fp, &obj->history[j].pID);
    if (stat != CAPS_SUCCESS) return stat;
    stat = caps_readString(fp, &obj->history[j].user);
    if (stat != CAPS_SUCCESS) return stat;

    n = fread(obj->history[j].datetime, sizeof(short), 6, fp);
    if (n != 6) return CAPS_IOERR;

    n = fread(&obj->history[j].sNum, sizeof(CAPSLONG), 1, fp);
    if (n != 1) return CAPS_IOERR;
  }
 
  return CAPS_SUCCESS;
}


static int
caps_readAttrs(FILE *fp, egAttrs **attrx)
{
  int     nattr, i, stat;
  size_t  n;
  egAttr  *attr;
  egAttrs *attrs;

  *attrx = NULL;
  n = fread(&nattr, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if (nattr == 0) return CAPS_SUCCESS;

  attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
  if (attrs == NULL) return EGADS_MALLOC;
  attr  = (egAttr *)  EG_alloc(nattr*sizeof(egAttr));
  if (attr == NULL) {
    EG_free(attrs);
    return EGADS_MALLOC;
  }
  attrs->nattrs = nattr;
  attrs->attrs  = attr;
  attrs->nseqs  = 0;
  attrs->seqs   = NULL;
  for (i = 0; i < nattr; i++) {
    attr[i].name   = NULL;
    attr[i].length = 1;
    attr[i].type   = ATTRINT;
  }

/*@-mustfreefresh@*/
  for (i = 0; i < nattr; i++) {
    n = fread(&attr[i].type,   sizeof(int), 1, fp);
    if (n != 1) {
      caps_freeAttrs(&attrs);
      return CAPS_IOERR;
    }
    n = fread(&attr[i].length, sizeof(int), 1, fp);
    if (n != 1) {
      caps_freeAttrs(&attrs);
      return CAPS_IOERR;
    }
    stat = caps_readString(fp, &attr[i].name);
    if (stat != CAPS_SUCCESS) {
      caps_freeAttrs(&attrs);
      return CAPS_IOERR;
    }
    if (attr[i].type == ATTRINT) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = fread(&attr[i].vals.integer, sizeof(int),              1, fp);
      } else if (attr[i].length > 1) {
        attr[i].vals.integers = (int *) EG_alloc(attr[i].length*sizeof(int));
        if (attr[i].vals.integers == NULL) {
          caps_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = fread(attr[i].vals.integers, sizeof(int), attr[i].length, fp);
      }
      if (n != attr[i].length) {
        caps_freeAttrs(&attrs);
        return CAPS_IOERR;
      }
    } else if (attr[i].type == ATTRREAL) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = fread(&attr[i].vals.real, sizeof(double),              1, fp);
      } else if (attr[i].length > 1) {
        attr[i].vals.reals = (double *) EG_alloc(attr[i].length*sizeof(double));
        if (attr[i].vals.reals == NULL) {
          caps_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = fread(attr[i].vals.reals, sizeof(double), attr[i].length, fp);
      }
      if (n != attr[i].length) return CAPS_IOERR;
    } else {
      stat = caps_readStrings(fp, attr[i].length, &attr[i].vals.string);
      if (stat != CAPS_SUCCESS) return CAPS_IOERR;
    }
  }
/*@+mustfreefresh@*/

  *attrx = attrs;
  return CAPS_SUCCESS;
}


static int
caps_readValue(FILE *fp, capsProblem *problem, capsObject *obj)
{
  int       i, stat;
  char      *name;
  size_t    n;
  capsValue *value;
  
  value = (capsValue *) obj->blind;
 
  stat = caps_readHistory(fp, obj);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_readOwn(fp, &obj->last);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_readAttrs(fp, &obj->attrs);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_readString(fp, &obj->name);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;
  stat = caps_readString(fp, &name);
  EG_free(name);
  if (stat != CAPS_SUCCESS) return CAPS_IOERR;

  n = fread(&value->type,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->length,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->dim,     sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->nrow,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->ncol,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->lfixed,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->sfixed,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->nullVal, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->index,   sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->pIndex,  sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->gInType, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  n = fread(&value->nderiv,    sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if (value->type == Integer) {
    n = fread(value->limits.ilims, sizeof(int),    2, fp);
    if (n != 2) return CAPS_IOERR;
  } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
    n = fread(value->limits.dlims, sizeof(double), 2, fp);
    if (n != 2) return CAPS_IOERR;
  }

  stat = caps_readString(fp, &value->units);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_readString(fp, &value->meshWriter);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_readString(fp, &name);
  if (stat != CAPS_SUCCESS) return stat;
  stat = caps_string2obj(problem, name, &value->link);
  EG_free(name);
  if (stat != CAPS_SUCCESS) return stat;
  n = fread(&value->linkMethod, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;

  if ((value->length == 1)     && (value->type != String) &&
      (value->type != Pointer) && (value->type != Tuple)  &&
      (value->type != PointerMesh)) {
    if ((value->type == Double) || (value->type == DoubleDeriv)) {
      n = fread(&value->vals.real, sizeof(double), 1, fp);
      if (n != 1) return CAPS_IOERR;
    } else {
      n = fread(&value->vals.integer, sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
    }
  } else {
    if ((value->type == Pointer) || (value->type == PointerMesh)) {
      /* what do we do? */
    } else if ((value->type == Double) || (value->type == DoubleDeriv)) {
      value->vals.reals = EG_alloc(value->length*sizeof(double));
      if (value->vals.reals == NULL) return EGADS_MALLOC;
      n = fread(value->vals.reals, sizeof(double), value->length, fp);
      if (n != value->length) return CAPS_IOERR;
    } else if (value->type == String) {
      stat = caps_readStrings(fp, value->length, &value->vals.string);
      if (stat != CAPS_SUCCESS) return stat;
    } else if (value->type == Tuple) {
      value->vals.tuple = NULL;
      if (value->length != 0) {
        stat = caps_readTuple(fp, value->length, value->nullVal,
                              &value->vals.tuple);
        if (stat != CAPS_SUCCESS) return stat;
      }
    } else {
      value->vals.integers = EG_alloc(value->length*sizeof(int));
      if (value->vals.integers == NULL) return EGADS_MALLOC;
      n = fread(value->vals.integers, sizeof(int), value->length, fp);
      if (n != value->length) return CAPS_IOERR;
    }
  }

  if (value->nullVal == IsPartial) {
    value->partial = EG_alloc(value->length*sizeof(int));
    if (value->partial == NULL) return EGADS_MALLOC;
    n = fread(value->partial, sizeof(int), value->length, fp);
    if (n != value->length) return CAPS_IOERR;
  }

  if (value->nderiv != 0) {
    value->derivs = (capsDeriv *) EG_alloc(value->nderiv*sizeof(capsDeriv));
    if (value->derivs == NULL) return EGADS_MALLOC;
    for (i = 0; i < value->nderiv; i++) {
      value->derivs[i].name = NULL;
      value->derivs[i].rank = 0;
      value->derivs[i].deriv  = NULL;
    }
    for (i = 0; i < value->nderiv; i++) {
      stat = caps_readString(fp, &value->derivs[i].name);
      if (stat != CAPS_SUCCESS) return stat;
      n = fread(&value->derivs[i].rank, sizeof(int), 1, fp);
      if (n != 1) return CAPS_IOERR;
      value->derivs[i].deriv = (double *) EG_alloc(
                           value->length*value->derivs[i].rank*sizeof(double));
      if (value->derivs[i].deriv == NULL) return EGADS_MALLOC;
      n = fread(value->derivs[i].deriv, sizeof(double),
                value->length*value->derivs[i].rank, fp);
      if (n != value->length*value->derivs[i].rank) return CAPS_IOERR;
    }
  }

  return CAPS_SUCCESS;
}


static int
caps_readErrs(capsProblem *problem, capsErrs **errx)
{
  int      i, j, nError, stat;
  char     *full;
  capsErrs *errs;
  size_t   n;
  FILE     *fp;
  
  *errx = NULL;
  fp    = problem->jrnl;
  n     = fread(&nError, sizeof(int), 1, fp);
  if (n != 1) return CAPS_IOERR;
  if (nError == 0) return CAPS_SUCCESS;
  
  errs = (capsErrs *) EG_alloc(sizeof(capsErrs));
  if (errs == NULL) return EGADS_MALLOC;
  errs->nError = nError;
  errs->errors = (capsError *) EG_alloc(nError*sizeof(capsError));
  if (errs->errors == NULL) {
    EG_free(errs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nError; i++) {
    errs->errors[i].nLines = 0;
    errs->errors[i].lines  = NULL;
  }
  
  for (i = 0; i < nError; i++) {
    stat = caps_readString(fp, &full);
    if (stat != CAPS_SUCCESS) {
      printf(" CAPS Warning: caps_readString = %d (caps_readErrs)\n", stat);
      caps_freeError(errs);
      return CAPS_IOERR;
    }
    stat = caps_string2obj(problem, full, &errs->errors[i].errObj);
    EG_free(full);
    if (stat != CAPS_SUCCESS)
      printf(" CAPS Warning: caps_string2obj = %d (caps_readErrs)\n", stat);

    n = fread(&errs->errors[i].eType,  sizeof(int), 1, fp);
    if (n != 1) goto readErr;
    n = fread(&errs->errors[i].index,  sizeof(int), 1, fp);
    if (n != 1) goto readErr;
    n = fread(&errs->errors[i].nLines, sizeof(int), 1, fp);
    if (n != 1) goto readErr;
    errs->errors[i].lines = (char **) EG_alloc(errs->errors[i].nLines*
                                               sizeof(char *));
    if (errs->errors[i].lines == NULL) {
      caps_freeError(errs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < errs->errors[i].nLines; j++) {
      stat = caps_readString(fp, &errs->errors[i].lines[j]);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Warning: %d caps_readString = %d (caps_readErrs)\n",
               j, stat);
        goto readErr;
      }
    }
  }
  
  *errx = errs;
  return CAPS_SUCCESS;
  
readErr:
/*@-dependenttrans@*/
  caps_freeError(errs);
/*@+dependenttrans@*/
  return CAPS_IOERR;
}


int
caps_jrnlRead(capsProblem *problem, capsObject *obj, int nargs, capsJrnl *args,
              CAPSLONG *serial, int *status)
{
  int       i, j, k, stat, funID;
  char      filename[PATH_MAX], *full;
  CAPSLONG  sNum0, sNum, objSN;
  capsFList *flist;
  size_t    n;
#ifdef WIN32
  __int64   fpos;
#else
  long      fpos;
#endif
  
  *serial = 0;
  *status = CAPS_SUCCESS;
  if (problem->jrnl   == NULL)            return CAPS_SUCCESS;
  if (problem->stFlag == CAPS_JOURNALERR) return CAPS_JOURNALERR;
  if (problem->stFlag != 4)               return CAPS_SUCCESS;
  
#ifdef WIN32
  fpos = _ftelli64(problem->jrnl);
#else
  fpos = ftell(problem->jrnl);
#endif
  
  /* are we at the last success? */
  if (fpos >= problem->jpos) {
    printf(" CAPS Info: Hit last success -- going live!\n");
    problem->stFlag = 0;
    fclose(problem->jrnl);
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal", problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal",   problem->root);
#endif
/*@-dependenttrans@*/
    problem->jrnl = fopen(filename, "ab");
/*@+dependenttrans@*/
    if (problem->jrnl == NULL) {
      printf(" CAPS Error: Cannot open %s (caps_readJrnl)\n", filename);
      return CAPS_DIRERR;
    }
#ifdef WIN32
    _fseeki64(problem->jrnl, problem->jpos, SEEK_SET);
#else
    fseek(problem->jrnl,     problem->jpos, SEEK_SET);
#endif
    return CAPS_SUCCESS;
  }
  
  /* lets get our record */
  n = fread(&funID,   sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jreaderr;
  if (funID != problem->funID) {
    printf(" CAPS Fatal: FunID = %s, should be '%s'!\n", caps_funID[funID], caps_funID[problem->funID]);
    goto jreadfatal;
  }
  n = fread(&sNum0,   sizeof(CAPSLONG), 1, problem->jrnl);
  if (n != 1) goto jreaderr;
  n = fread(status,   sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jreaderr;
/*
  printf(" *** Reading from FunID = %s  status = %d ***\n", caps_funID[funID], *status);
*/
  if (*status == CAPS_SUCCESS) {
    n = fread(&objSN, sizeof(CAPSLONG), 1, problem->jrnl);
    if (n != 1) goto jreaderr;
    /* cleanup? */
    if (obj != NULL)
      if (obj->flist != NULL) {
        flist = (capsFList *) obj->flist;
        if (objSN > flist->sNum) caps_freeFList(obj);
      }
    
    for (i = 0; i < nargs; i++)
      switch (args[i].type) {
          
        case jInteger:
          n = fread(&args[i].members.integer, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          break;
          
        case jDouble:
          n = fread(&args[i].members.real, sizeof(double), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          break;
          
        case jString:
          stat = caps_readString(problem->jrnl, &args[i].members.string);
          if (stat == CAPS_IOERR) goto jreaderr;
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_readString = %d!\n", stat);
            goto jreadfatal;
          }
          if (obj != NULL) {
            flist = (capsFList *) EG_alloc(sizeof(capsFList));
            if (flist == NULL) {
              printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
            } else {
              flist->type           = jPointer;
              flist->num            = 1;
              flist->member.pointer = args[i].members.string;
              flist->sNum           = objSN;
              flist->next           = NULL;
              if (obj->flist != NULL) flist->next = obj->flist;
              obj->flist            = flist;
            }
          }
          break;
          
        case jStrings:
          n = fread(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          args[i].members.strings = (char **) EG_alloc(args[i].num*
                                                       sizeof(char *));
          if (args[i].members.strings == NULL) {
            printf(" CAPS Warning: Journal strings Malloc Error!\n");
            goto jreadfatal;
          }
          for (j = 0; j < args[i].num; j++) {
            stat = caps_readString(problem->jrnl, &args[i].members.strings[j]);
            if (stat != CAPS_SUCCESS) {
              printf(" CAPS Warning: Jrnl %d caps_readString Str = %d!\n",
                     j, stat);
              for (k = 0; k < j; k++) EG_free(args[i].members.strings[k]);
              EG_free(args[i].members.strings);
              args[i].members.strings = NULL;
              goto jreaderr;
            }
          }
          if (obj != NULL) {
            flist = (capsFList *) EG_alloc(sizeof(capsFList));
            if (flist == NULL) {
              printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
            } else {
              flist->type           = jStrings;
              flist->num            = args[i].num;
              flist->member.strings = args[i].members.strings;
              flist->sNum           = objSN;
              flist->next           = NULL;
              if (obj->flist != NULL) flist->next = obj->flist;
              obj->flist            = flist;
            }
          }
          break;
          
        case jTuple:
          n = fread(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          stat = caps_readTuple(problem->jrnl, args[i].num, NotNull,
                                &args[i].members.tuple);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_readTuple = %d!\n", stat);
            goto jreaderr;
          }
          if (obj != NULL) {
            flist = (capsFList *) EG_alloc(sizeof(capsFList));
            if (flist == NULL) {
              printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
            } else {
              flist->type         = jTuple;
              flist->num          = args[i].num;
              flist->member.tuple = args[i].members.tuple;
              flist->sNum         = objSN;
              flist->next         = NULL;
              if (obj->flist != NULL) flist->next = obj->flist;
              obj->flist          = flist;
            }
          }
          break;
          
        case jPointer:
        case jPtrFree:
          n = fread(&args[i].length, sizeof(size_t), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          args[i].members.pointer = NULL;
          if (args[i].length != 0) {
            args[i].members.pointer = (char *) EG_alloc(args[i].length*
                                                        sizeof(char));
            if (args[i].members.pointer == NULL) {
              printf(" CAPS Warning: Journal Pointer Malloc Error!\n");
              goto jreadfatal;
            }
            n = fread(args[i].members.pointer, sizeof(char), args[i].length,
                      problem->jrnl);
            if (n != args[i].length) {
              EG_free(args[i].members.pointer);
              args[i].members.pointer = NULL;
              goto jreaderr;
            }
            if ((obj != NULL) && (args[i].type == jPointer)) {
              flist = (capsFList *) EG_alloc(sizeof(capsFList));
              if (flist == NULL) {
                printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
              } else {
                flist->type           = jPointer;
                flist->num            = 1;
                flist->member.pointer = args[i].members.pointer;
                flist->sNum           = objSN;
                flist->next           = NULL;
                if (obj->flist != NULL) flist->next = obj->flist;
                obj->flist            = flist;
              }
            }
          }
          break;
          
        case jObject:
          stat = caps_readString(problem->jrnl, &full);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Jrnl caps_readString Obj = %d!\n", stat);
            goto jreaderr;
          }
          stat = caps_string2obj(problem, full, &args[i].members.obj);
          EG_free(full);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_string2obj = %d!\n", stat);
            goto jreaderr;
          }
          break;
          
        case jObjs:
          n = fread(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          if (args[i].num != 0) {
            args[i].members.objs = (capsObject **) EG_alloc(args[i].num*
                                                          sizeof(capsObject *));
            if (args[i].members.objs == NULL) {
              printf(" CAPS Warning: Journal Objects Malloc Error!\n");
              goto jreadfatal;
            }
            for (j = 0; j < args[i].num; j++) {
              stat = caps_readString(problem->jrnl, &full);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Jrnl caps_readString Obj = %d!\n", stat);
                EG_free(args[i].members.objs);
                args[i].members.objs = NULL;
                goto jreaderr;
              }
              stat = caps_string2obj(problem, full, &args[i].members.objs[j]);
              EG_free(full);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Journal caps_string2obj = %d!\n", stat);
                EG_free(args[i].members.objs);
                args[i].members.objs = NULL;
                goto jreaderr;
              }
            }
            if (obj != NULL) {
              flist = (capsFList *) EG_alloc(sizeof(capsFList));
              if (flist == NULL) {
                printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
              } else {
                flist->type           = jPointer;
                flist->num            = 1;
                flist->member.pointer = args[i].members.objs;
                flist->sNum           = objSN;
                flist->next           = NULL;
                if (obj->flist != NULL) flist->next = obj->flist;
                obj->flist            = flist;
              }
            }
          }
          break;
          
        case jValObj:
          stat = caps_readValue(problem->jrnl, problem, args[i].members.obj);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Jrnl caps_readValue = %d!\n", stat);
            caps_delete(args[i].members.obj);
            goto jreaderr;
          }
          break;
          
        case jErr:
          stat = caps_readErrs(problem, &args[i].members.errs);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_readErrs = %d!\n", stat);
            goto jreaderr;
          }
          break;
          
        case jOwn:
          stat = caps_readOwn(problem->jrnl, &args[i].members.own);
          if (stat != CAPS_SUCCESS) {
            printf(" CAPS Warning: Journal caps_Own = %d!\n", stat);
            goto jreaderr;
          }
          if (obj != NULL) {
            flist = (capsFList *) EG_alloc(sizeof(capsFList));
            if (flist == NULL) {
              printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
            } else {
              flist->type       = jOwn;
              flist->member.own = args[i].members.own;
              flist->sNum       = objSN;
              flist->next       = NULL;
              if (obj->flist != NULL) flist->next = obj->flist;
              obj->flist          = flist;
            }
          }
          break;
        
        case jOwns:
          n = fread(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          if (args[i].num != 0) {
            args[i].members.owns = (capsOwn *) EG_alloc(args[i].num*
                                                        sizeof(capsOwn));
            if (args[i].members.owns == NULL) {
              printf(" CAPS Warning: Journal Owner Malloc Error!\n");
              goto jreadfatal;
            }
            for (j = 0; j < args[i].num; j++) {
              stat = caps_readOwn(problem->jrnl, &args[i].members.owns[j]);
              if (stat != CAPS_SUCCESS) {
                printf(" CAPS Warning: Journal caps_Owns %d = %d!\n", j, stat);
                EG_free(args[i].members.owns);
                args[i].members.owns = NULL;
                goto jreaderr;
              }
            }
            if (obj != NULL) {
              flist = (capsFList *) EG_alloc(sizeof(capsFList));
              if (flist == NULL) {
                printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
              } else {
                flist->type        = jOwns;
                flist->num         = args[i].num;
                flist->member.owns = args[i].members.owns;
                flist->sNum        = objSN;
                flist->next        = NULL;
                if (obj->flist != NULL) flist->next = obj->flist;
                obj->flist         = flist;
              }
            }
          }
          break;

        case jEgos:
          args[i].members.model = NULL;
          n = fread(&args[i].num, sizeof(int), 1, problem->jrnl);
          if (n != 1) goto jreaderr;
          if (args[i].num != -1) {
#ifdef WIN32
            snprintf(filename, PATH_MAX, "%s\\capsRestart\\model%4.4d.egads",
                     problem->root, args[i].num);
#else
            snprintf(filename, PATH_MAX, "%s/capsRestart/model%4.4d.egads",
                     problem->root, args[i].num);
#endif
            stat = EG_loadModel(problem->context, 1, filename,
                                &args[i].members.model);
            if (stat != CAPS_SUCCESS) {
              printf(" CAPS Warning: EG_loadModel = %d (caps_jrnlRead)!\n",
                     stat);
              goto jreaderr;
            }
            if (obj != NULL) {
              flist = (capsFList *) EG_alloc(sizeof(capsFList));
              if (flist == NULL) {
                printf(" CAPS Warning: Cannot Allocate Journal Free List!\n");
              } else {
                flist->type         = jEgos;
                flist->num          = args[i].num;
                flist->member.model = args[i].members.model;
                flist->sNum         = objSN;
                flist->next         = NULL;
                if (obj->flist != NULL) flist->next = obj->flist;
                obj->flist         = flist;
              }
            }
          }
          
      }
  }
  
  n = fread(&sNum,  sizeof(CAPSLONG), 1, problem->jrnl);
  if (n != 1) goto jreaderr;
  n = fread(&funID, sizeof(int),      1, problem->jrnl);
  if (n != 1) goto jreaderr;
  if (funID != problem->funID) {
    printf(" CAPS Fatal: Ending FunID = %d, should be %d!\n",
           funID, problem->funID);
    goto jreadfatal;
  }
  if (sNum > problem->sNum) {
    printf(" CAPS Info: Hit ending serial number -- going live!\n");
    problem->stFlag = 0;
    fclose(problem->jrnl);
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal", problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal",   problem->root);
#endif
/*@-dependenttrans@*/
    problem->jrnl = fopen(filename, "ab");
/*@+dependenttrans@*/
    if (problem->jrnl == NULL) {
      printf(" CAPS Error: Cannot open %s (caps_readJrnl)\n", filename);
      return CAPS_DIRERR;
    }
  }

  *serial = sNum;
  return CAPS_JOURNAL;
  
jreaderr:
  printf(" CAPS Info: Incomplete Journal Record @ %d -- going live!\n",
         problem->funID);
  problem->stFlag = 0;
  fclose(problem->jrnl);
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal",   problem->root);
#endif
/*@-dependenttrans@*/
  problem->jrnl = fopen(filename, "ab");
/*@+dependenttrans@*/
  if (problem->jrnl == NULL) {
    printf(" CAPS Error: Cannot open %s (caps_readJrnl)\n", filename);
    return CAPS_DIRERR;
  }
#ifdef WIN32
  _fseeki64(problem->jrnl, fpos, SEEK_SET);
#else
  fseek(problem->jrnl,     fpos, SEEK_SET);
#endif
  return CAPS_SUCCESS;
  
jreadfatal:
  fclose(problem->jrnl);
  problem->jrnl   = NULL;
  problem->stFlag = CAPS_JOURNALERR;
  return CAPS_JOURNALERR;
}


static int
caps_readInitObj(capsObject **obj, int type, int subtype, char *name,
                 capsObject *parent)
{
  int        status;
  capsObject *object;
  
  status = caps_makeObject(obj);
  if (status != CAPS_SUCCESS) return status;

  object              = *obj;
  object->magicnumber = CAPSMAGIC;
  object->type        = type;
  object->subtype     = subtype;
  object->name        = EG_strdup(name);
  object->attrs       = NULL;
  object->blind       = NULL;
  object->flist       = NULL;
  object->parent      = parent;
  
  return CAPS_SUCCESS;
}


static int
caps_readDataSet(capsProblem *problem, const char *base, capsObject *dobject)
{
  int         i, stat;
  char        filename[PATH_MAX], *name = NULL;
  capsDataSet *ds;
  size_t      n;
  FILE        *fp;
  
  ds = (capsDataSet *) EG_alloc(sizeof(capsDataSet));
  if (ds == NULL) return EGADS_MALLOC;
  ds->ftype      = BuiltIn;
  ds->npts       = 0;
  ds->rank       = 0;
  ds->data       = NULL;
  ds->units      = NULL;
  ds->startup    = NULL;
  ds->linkMethod = Interpolate;
  ds->link       = NULL;
  dobject->blind = ds;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\DN-%s", base, dobject->name);
#else
  snprintf(filename, PATH_MAX, "%s/DN-%s",  base, dobject->name);
#endif
  
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  
  stat = caps_readHistory(fp, dobject);
  if (stat != CAPS_SUCCESS) goto rderrs;
  stat = caps_readOwn(fp, &dobject->last);
  if (stat != CAPS_SUCCESS) goto rderrs;
  stat = caps_readAttrs(fp, &dobject->attrs);
  if (stat != CAPS_SUCCESS) goto rderrs;
  stat = caps_readString(fp, &name);
  if (stat != CAPS_SUCCESS) goto rderrs;

  n = fread(&ds->ftype,  sizeof(int), 1, fp);
  if (n != 1) goto rderrs;
  n = fread(&ds->npts,   sizeof(int), 1, fp);
  if (n != 1) goto rderrs;
  n = fread(&ds->rank,   sizeof(int), 1, fp);
  if (n != 1) goto rderrs;
  stat = caps_readString(fp, &ds->units);
  if (stat != CAPS_SUCCESS) goto rderrs;
  
  stat = caps_readDoubles(fp, &i, &ds->data);
  if (stat != CAPS_SUCCESS) goto rderrs;
  if (i != ds->npts*ds->rank) {
    printf(" CAPS Error: %s len mismatch %d %d (caps_readDataSet)!\n",
           name, i, ds->npts*ds->rank);
    goto rderrs;
  }
  EG_free(name);
  name = NULL;
  
  stat = caps_readString(fp, &name);
  if (stat != CAPS_SUCCESS) goto rderrs;
  if (name != NULL) {
    stat = caps_string2obj(problem, name, &ds->link);
    EG_free(name);
    name = NULL;
    if (stat != CAPS_SUCCESS) goto rderrs;
  }
  n = fread(&ds->linkMethod, sizeof(int), 1, fp);
  if (n != 1) goto rderrs;
  fclose(fp);

  return CAPS_SUCCESS;
  
rderrs:
  if (name != NULL) EG_free(name);
  fclose(fp);
  return CAPS_IOERR;
}


static int
caps_readVertexSet(capsProblem *problem, capsObject *bobject,
                   capsObject *vobject)
{
  int           i, dim, status;
  char          base[PATH_MAX], filename[PATH_MAX], cstype, *name;
  capsBound     *bound;
  capsAnalysis  *analysis;
  capsVertexSet *vs;
  capsObject    *aobject;
  size_t        n;
  FILE          *fp;

  bound  = (capsBound *)     bobject->blind;
  vs     = (capsVertexSet *) vobject->blind;
  cstype = 'U';
  if (vobject->subtype == CONNECTED) {
    cstype   = 'C';
    aobject  = vs->analysis;
    analysis = (capsAnalysis *) aobject->blind;
    vs->discr->aInfo     = &analysis->info;
    vs->discr->instStore =  analysis->instStore;
  }
  
#ifdef WIN32
  snprintf(base, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\S%c-%s",
           problem->root, bound->index, cstype, vobject->name);
  snprintf(filename, PATH_MAX, "%s\\vs", base);
#else
  snprintf(base, PATH_MAX, "%s/capsRestart/BN-%4.4d/S%c-%s",
           problem->root, bound->index, cstype, vobject->name);
  snprintf(filename, PATH_MAX, "%s/vs",  base);
#endif
  
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  
  status = caps_readHistory(fp, vobject);
  if (status != CAPS_SUCCESS) goto readerr;
  status = caps_readOwn(fp, &vobject->last);
  if (status != CAPS_SUCCESS) goto readerr;
  status = caps_readAttrs(fp, &vobject->attrs);
  if (status != CAPS_SUCCESS) goto readerr;
  status = caps_readString(fp, &name);
  if (status != CAPS_SUCCESS) goto readerr;
  EG_free(name);

  /* read needed discr info */
  n = fread(&dim, sizeof(int), 1, fp);
  if (n != 1) goto readerr;
  vs->discr->dim = dim;
  if (vobject->subtype == CONNECTED) dim = 0;

  if (dim != 0) {
    n = fread(&vs->discr->nVerts, sizeof(int), 1, fp);
    if (n != 1) goto readerr;
    vs->discr->verts = (double *) EG_alloc(3*vs->discr->nVerts*sizeof(double));
    if (vs->discr->verts == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    for (i = 0; i < vs->discr->nVerts; i++) {
      n = fread(&vs->discr->verts[3*i], sizeof(double), 3, fp);
      if (n != 3) goto readerr;
    }
  }
  fclose(fp);
  
  /* load the DataSets */
  for (i = 0; i < vs->nDataSets; i++) {
    status = caps_readDataSet(problem, base, vs->dataSets[i]);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: %s caps_readDataSet = %d\n",
             vs->dataSets[i]->name, status);
      return status;
    }
  }
  
  return CAPS_SUCCESS;
  
readerr:
  fclose(fp);
  return CAPS_IOERR;
}


static int
caps_readInitDSets(capsProblem *problem, capsObject *bobject,
                   capsObject *vobject)
{
  int           i, status;
  char          filename[PATH_MAX], name[PATH_MAX], cstype;
  capsBound     *bound;
  capsVertexSet *vs;
  FILE          *fp;

  cstype = 'U';
  if (vobject->subtype == CONNECTED) cstype = 'C';
  
  bound = (capsBound *)     bobject->blind;
  vs    = (capsVertexSet *) vobject->blind;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\S%c-%s\\dSets.txt",
           problem->root, bound->index, cstype, vobject->name);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d/S%c-%s/dSets.txt",
           problem->root, bound->index, cstype, vobject->name);
#endif

  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s (caps_readInitDSets)!\n", filename);
    return CAPS_DIRERR;
  }
  fscanf(fp, "%d", &vs->nDataSets);
  vs->dataSets = (capsObject **) EG_alloc(vs->nDataSets*sizeof(capsObject *));
  if (vs->dataSets == NULL) {
    fclose(fp);
    return EGADS_MALLOC;
  }
  for (i = 0; i < vs->nDataSets; i++) vs->dataSets[i] = NULL;
  for (i = 0; i < vs->nDataSets; i++) {
    fscanf(fp, "%s", name);
    status = caps_readInitObj(&vs->dataSets[i], DATASET, NONE, name, vobject);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: %s caps_readInitObj = %d (caps_readInitDSets)!\n",
             name, status);
      fclose(fp);
      return status;
    }
  }
  fclose(fp);
  
  return CAPS_SUCCESS;
}


static int
caps_readInitVSets(capsProblem *problem, capsObject *bobject)
{
  int           i, j, type, status;
  char          filename[PATH_MAX], name[PATH_MAX];
  capsObject    *aobject;
  capsBound     *bound;
  capsVertexSet *vs;
  FILE          *fp;
  
  bound = (capsBound *) bobject->blind;
  
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\vSets.txt",
           problem->root, bound->index);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d/vSets.txt",
           problem->root, bound->index);
#endif
  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s (caps_readInitVSets)!\n", filename);
    return CAPS_DIRERR;
  }
  fscanf(fp, "%d", &bound->nVertexSet);
  bound->vertexSet = (capsObject **)
                     EG_alloc(bound->nVertexSet*sizeof(capsObject *));
  if (bound->vertexSet == NULL) {
    fclose(fp);
    return EGADS_MALLOC;
  }
  for (i = 0; i < bound->nVertexSet; i++) bound->vertexSet[i] = NULL;
  for (i = 0; i < bound->nVertexSet; i++) {
    fscanf(fp, "%d %s", &type, name);
    aobject = NULL;
    if (type == CONNECTED) {
      for (j = 0; j < problem->nAnalysis; j++)
        if (strcmp(name, problem->analysis[j]->name) == 0) {
          aobject = problem->analysis[j];
          break;
        }
      if (aobject == NULL) {
        fclose(fp);
        printf(" CAPS Error: Analysis %s Not Found (caps_readInitVSets)\n",
               name);
        return CAPS_NOTFOUND;
      }
    }
    status = caps_readInitObj(&bound->vertexSet[i], VERTEXSET, type, name,
                              bobject);
    if (status != CAPS_SUCCESS) {
      fclose(fp);
      printf(" CAPS Error: caps_readInitObj = %d (caps_readInitVSets)\n",
             status);
      return status;
    }
    vs = (capsVertexSet *) EG_alloc(sizeof(capsVertexSet));
    if (vs == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    vs->analysis  = aobject;
    vs->nDataSets = 0;
    vs->dataSets  = NULL;
    vs->discr     = (capsDiscr *) EG_alloc(sizeof(capsDiscr));
    if (vs->discr == NULL) {
      EG_free(vs);
      fclose(fp);
      return EGADS_MALLOC;
    }
    vs->discr->dim             = 0;
    vs->discr->instStore       = NULL;
    vs->discr->aInfo           = NULL;
    vs->discr->nVerts          = 0;
    vs->discr->verts           = NULL;
    vs->discr->celem           = NULL;
    vs->discr->nDtris          = 0;
    vs->discr->dtris           = NULL;
    vs->discr->nPoints         = 0;
    vs->discr->nTypes          = 0;
    vs->discr->types           = NULL;
    vs->discr->nBodys          = 0;
    vs->discr->bodys           = NULL;
    vs->discr->tessGlobal      = NULL;
    vs->discr->ptrm            = NULL;
    bound->vertexSet[i]->blind = vs;
    
    status = caps_readInitDSets(problem, bobject, bound->vertexSet[i]);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Bound %d caps_readInitDSets = %d (caps_open)!\n",
             i, status);
      fclose(fp);
      return status;
    }
    
  }
  fclose(fp);
  
  return CAPS_SUCCESS;
}


static int
caps_readBound(capsObject *bobject)
{
  int         i, status;
  char        filename[PATH_MAX];
  capsObject  *pobject;
  capsProblem *problem;
  capsBound   *bound;
  size_t      n;
  FILE        *fp;
  
  pobject = bobject->parent;
  problem = (capsProblem *) pobject->blind;
  bound   = (capsBound *)   bobject->blind;

#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\BN-%4.4d\\bound",
           problem->root, bound->index);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/BN-%4.4d/bound",
           problem->root, bound->index);
#endif
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open bound %s!\n", filename);
    return CAPS_DIRERR;
  }
  
  status = caps_readHistory(fp, bobject);
  if (status != CAPS_SUCCESS) goto rerror;
  status = caps_readOwn(fp, &bobject->last);
  if (status != CAPS_SUCCESS) goto rerror;
  status = caps_readAttrs(fp, &bobject->attrs);
  if (status != CAPS_SUCCESS) goto rerror;
  status = caps_readString(fp, &bobject->name);
  if (status != CAPS_SUCCESS) goto rerror;
  n = fread(&bound->dim,     sizeof(int),    1, fp);
  if (n != 1) goto rerror;
  n = fread(&bound->state,   sizeof(int),    1, fp);
  if (n != 1) goto rerror;
  n = fread(&bound->plimits, sizeof(double), 4, fp);
  if (n != 4) goto rerror;
  n = fread(&bound->iBody,   sizeof(int),    1, fp);
  if (n != 1) goto rerror;
  n = fread(&bound->iEnt,    sizeof(int),    1, fp);
  if (n != 1) goto rerror;
  n = fread(&bound->index,   sizeof(int),    1, fp);
  if (n != 1) goto rerror;

  n = fread(&i,              sizeof(int),    1, fp);
  if (n != 1) goto rerror;
  if (i != 0) {
    bound->curve = (capsAprx1D *) EG_alloc(sizeof(capsAprx1D));
    if (bound->curve == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    bound->curve->nrank = i;
  }
  
  if (bound->curve != NULL) {
    n = fread(&bound->curve->periodic, sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    n = fread(&bound->curve->nts,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    status = caps_readDoubles(fp, &i, &bound->curve->interp);
    if (status != CAPS_SUCCESS) goto rerror;
    n = fread(bound->curve->trange,    sizeof(double),   2, fp);
    if (n != 2) goto rerror;
    n = fread(&bound->curve->ntm,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    status = caps_readDoubles(fp, &i, &bound->curve->tmap);
    if (status != CAPS_SUCCESS) goto rerror;
  }

  n = fread(&i,                          sizeof(int),      1, fp);
  if (n != 1) goto rerror;
  if (i != 0) {
    bound->surface = (capsAprx2D *) EG_alloc(sizeof(capsAprx2D));
    if (bound->surface == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    bound->surface->nrank = i;
  }
  if (bound->surface != NULL) {
    n = fread(&bound->surface->periodic, sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    n = fread(&bound->surface->nus,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    n = fread(&bound->surface->nvs,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    status = caps_readDoubles(fp, &i, &bound->surface->interp);
    if (status != CAPS_SUCCESS) goto rerror;
    n = fread(bound->surface->urange,    sizeof(double),   2, fp);
    if (n != 2) goto rerror;
    n = fread(bound->surface->vrange,    sizeof(double),   2, fp);
    if (n != 2) goto rerror;
    n = fread(&bound->surface->num,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    n = fread(&bound->surface->nvm,      sizeof(int),      1, fp);
    if (n != 1) goto rerror;
    status = caps_readDoubles(fp, &i, &bound->surface->uvmap);
    if (status != CAPS_SUCCESS) goto rerror;
  }
  fclose(fp);
  
  /* get VertexSets */
  for (i = 0; i < bound->nVertexSet; i++) {
    status = caps_readVertexSet(problem, bobject, bound->vertexSet[i]);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: caps_readVertexSet = %d (caps_readBound)\n", status);
      return status;
    }
  }

  return CAPS_SUCCESS;
  
rerror:
  fclose(fp);
  return CAPS_IOERR;
}


static int
caps_readAnalysis(capsProblem *problem, capsObject *aobject)
{
  int          i, stat, eFlag, nIn, nOut, nField, *ranks, *fInOut;
  size_t       n;
  char         **fields, base[PATH_MAX], filename[PATH_MAX], *name;
  void         *instStore;
  capsAnalysis *analysis;
  FILE         *fp;

  analysis = (capsAnalysis *) aobject->blind;
  
#ifdef WIN32
  snprintf(base,     PATH_MAX, "%s\\capsRestart\\AN-%s",
           problem->root, aobject->name);
  snprintf(filename, PATH_MAX, "%s\\analysis", base);
#else
  snprintf(base,     PATH_MAX, "%s/capsRestart/AN-%s",
           problem->root, aobject->name);
  snprintf(filename, PATH_MAX, "%s/analysis", base);
#endif
  
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
    return CAPS_DIRERR;
  }
  
  stat = caps_readHistory(fp, aobject);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readOwn(fp, &aobject->last);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readOwn(fp, &analysis->pre);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readAttrs(fp, &aobject->attrs);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readString(fp, &name);
  if (stat != CAPS_SUCCESS) goto raerr;
  EG_free(name);
  
  stat = caps_readString(fp, &analysis->loadName);
  if (stat != CAPS_SUCCESS) goto raerr;
  if (analysis->loadName == NULL) {
    fclose(fp);
    return CAPS_NULLNAME;
  }
  stat = caps_readString(fp, &analysis->path);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readString(fp, &analysis->unitSys);
  if (stat != CAPS_SUCCESS) goto raerr;
  stat = caps_readString(fp, &analysis->intents);
  if (stat != CAPS_SUCCESS) goto raerr;
  n = fread(&analysis->major,    sizeof(int), 1, fp);
  if (n != 1) goto raerr;
  n = fread(&analysis->minor,    sizeof(int), 1, fp);
  if (n != 1) goto raerr;
  n = fread(&analysis->autoexec, sizeof(int), 1, fp);
  if (n != 1) goto raerr;
  n = fread(&analysis->nField,   sizeof(int), 1, fp);
  if (n != 1) goto raerr;
  if (analysis->nField != 0) {
    analysis->fields = (char **) EG_alloc(analysis->nField*sizeof(char *));
    if (analysis->fields == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    for (i = 0; i < analysis->nField; i++) {
      stat = caps_readString(fp, &analysis->fields[i]);
      if (stat != CAPS_SUCCESS) goto raerr;
    }
    analysis->ranks = (int *) EG_alloc(analysis->nField*sizeof(int));
    if (analysis->ranks == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    n = fread(analysis->ranks,  sizeof(int), analysis->nField, fp);
    if (n != analysis->nField) goto raerr;
    analysis->fInOut = (int *) EG_alloc(analysis->nField*sizeof(int));
    if (analysis->fInOut == NULL) {
      fclose(fp);
      return EGADS_MALLOC;
    }
    n = fread(analysis->fInOut, sizeof(int), analysis->nField, fp);
    if (n != analysis->nField) goto raerr;
  }
  fclose(fp);
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\%s", problem->root, analysis->path);
#else
  snprintf(filename, PATH_MAX, "%s/%s",  problem->root, analysis->path);
#endif
  analysis->fullPath = EG_strdup(filename);

  /* try to load the AIM */
  eFlag     = 0;
  nField    = 0;
  fields    = NULL;
  ranks     = NULL;
  fInOut    = NULL;
  instStore = NULL;
  stat      = aim_Initialize(&problem->aimFPTR, analysis->loadName, &eFlag,
                             analysis->unitSys, &analysis->info,
                             &analysis->major, &analysis->minor, &nIn, &nOut,
                             &nField, &fields, &ranks, &fInOut, &instStore);
  if (stat < CAPS_SUCCESS) {
    if (fields != NULL) {
      for (i = 0; i < nField; i++) EG_free(fields[i]);
      EG_free(fields);
    }
    EG_free(ranks);
    EG_free(fInOut);
    return stat;
  }
  if (nIn <= 0) {
    if (fields != NULL) {
      for (i = 0; i < nField; i++) EG_free(fields[i]);
      EG_free(fields);
    }
    EG_free(ranks);
    EG_free(fInOut);
    return CAPS_BADINIT;
  }
  analysis->info.instance = stat;
  analysis->instStore     = instStore;
  analysis->eFlag         = eFlag;
  if ((analysis->autoexec == 1) && (eFlag == 0)) analysis->autoexec = 0;
  if (fields != NULL) {
    for (i = 0; i < nField; i++) EG_free(fields[i]);
    EG_free(fields);
  }
  EG_free(ranks);
  EG_free(fInOut);
  if (nField != analysis->nField) {
    printf(" CAPS Error: %s # Fields = %d -- from file = %d (caps_open)!\n",
           analysis->loadName, nField, analysis->nField);
    return CAPS_MISMATCH;
  }
  if (nIn != analysis->nAnalysisIn) {
    printf(" CAPS Error: %s # Inputs = %d -- from file = %d (caps_open)!\n",
           analysis->loadName, nIn, analysis->nAnalysisIn);
    return CAPS_MISMATCH;
  }
  if (nOut != analysis->nAnalysisOut) {
    printf(" CAPS Error: %s # Outputs = %d -- from file = %d (caps_open)!\n",
           analysis->loadName, nOut, analysis->nAnalysisOut);
    return CAPS_MISMATCH;
  }

  /* read the value objects */
  
  if (analysis->analysisIn != NULL) {
    for (i = 0; i < analysis->nAnalysisIn; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VI-%4.4d", base, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/VI-%4.4d",  base, i+1);
#endif
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
        return CAPS_DIRERR;
      }
      stat = caps_readValue(fp, problem, analysis->analysisIn[i]);
      fclose(fp);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: %s AnalysisIn %d/%d readValue = %d (caps_open)!\n",
               aobject->name, i+1, analysis->nAnalysisIn, stat);
        return stat;
      }
    }
  }
  
  if (analysis->analysisOut != NULL) {
    for (i = 0; i < analysis->nAnalysisOut; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\VO-%4.4d", base, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/VO-%4.4d",  base, i+1);
#endif
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
        return CAPS_DIRERR;
      }
      stat = caps_readValue(fp, problem, analysis->analysisOut[i]);
      fclose(fp);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: %s AnalysisOut %d/%d readValue = %d (caps_open)!\n",
               aobject->name, i+1, analysis->nAnalysisOut, stat);
        return stat;
      }
    }
  }
  
  return CAPS_SUCCESS;
  
raerr:
  fclose(fp);
  return CAPS_IOERR;
}


static int
caps_readState(capsObject *pobject)
{
  int          i, j, ivec[2], nIn, nOut, stat;
  char         filename[PATH_MAX], name[PATH_MAX], *phName;
  capsProblem  *problem;
  capsAnalysis *analysis;
  capsBound    *bound;
  capsValue    *value;
  size_t       n;
  FILE         *fp, *fptxt = NULL;

  problem = (capsProblem *) pobject->blind;
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\Problem", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/Problem",   problem->root);
#endif
  
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    printf(" CAPS Error: Cannot open %s!\n", filename);
    return CAPS_DIRERR;
  }
  ivec[0] = CAPSMAJOR;
  ivec[1] = CAPSMINOR;

  n = fread(ivec,               sizeof(int),      2, fp);
  if (n != 2) goto readerr;
#ifdef DEBUG
  printf(" CAPS Info: Reading files written by Major = %d Minor = %d\n",
         ivec[0], ivec[1]);
#endif
  n = fread(&pobject->subtype,  sizeof(int),      1, fp);
  if (n != 1) goto readerr;
  stat = caps_readHistory(fp, pobject);
  if (stat != CAPS_SUCCESS) goto readerr;
  stat = caps_readOwn(fp, &pobject->last);
  if (stat != CAPS_SUCCESS) goto readerr;
  stat = caps_readAttrs(fp, &pobject->attrs);
  if (stat != CAPS_SUCCESS) goto readerr;
  if (pobject->name != NULL) EG_free(pobject->name);
  stat = caps_readString(fp, &pobject->name);
  if (stat != CAPS_SUCCESS) goto readerr;
  stat = caps_readString(fp, &phName);
  if (stat != CAPS_SUCCESS) goto readerr;
  if (phName != NULL) EG_free(phName);
  stat = caps_readOwn(fp, &problem->geometry);
  if (stat != CAPS_SUCCESS) goto readerr;
  n = fread(&problem->sNum,      sizeof(CAPSLONG), 1, fp);
  if (n != 1) goto readerr;
#ifdef WIN32
  n = fread(&problem->jpos,      sizeof(__int64),  1, fp);
#else
  n = fread(&problem->jpos,      sizeof(long),     1, fp);
#endif
  if (n != 1) goto readerr;
  n = fread(&problem->outLevel,  sizeof(int),      1, fp);
  if (n != 1) goto readerr;
  n = fread(&problem->nEGADSmdl, sizeof(int),      1, fp);
  if (n != 1) goto readerr;
  n = fread(&problem->nRegGIN,   sizeof(int),      1, fp);
  if (n != 1) goto readerr;

  for (i = 0; i < problem->nRegGIN; i++) {
    stat = caps_readString(fp, &problem->regGIN[i].name);
    if (stat != CAPS_SUCCESS) goto readerr;
    n = fread(&problem->regGIN[i].index, sizeof(int), 1, fp);
    if (n != 1) goto readerr;
    n = fread(&problem->regGIN[i].irow,  sizeof(int), 1, fp);
    if (n != 1) goto readerr;
    n = fread(&problem->regGIN[i].icol,  sizeof(int), 1, fp);
    if (n != 1) goto readerr;
  }
  fclose(fp);
  
  /* get the number of objects */

#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\param.txt", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/param.txt",   problem->root);
#endif
  fptxt = fopen(filename, "r");
  if (fptxt != NULL) {
    fscanf(fptxt, "%d", &problem->nParam);
    fclose(fptxt);
  }
  
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\geom.txt", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/geom.txt",   problem->root);
#endif
  fptxt = fopen(filename, "r");
  if (fptxt != NULL) {
    fscanf(fptxt, "%d %d", &problem->nGeomIn, &problem->nGeomOut);
    fclose(fptxt);
  }
  
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\analy.txt", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/analy.txt",   problem->root);
#endif
  fptxt = fopen(filename, "r");
  if (fptxt != NULL) {
    fscanf(fptxt, "%d", &problem->nAnalysis);
    if (problem->nAnalysis > 0) {
      problem->analysis = (capsObject **)
                          EG_alloc(problem->nAnalysis*sizeof(capsObject *));
      if (problem->analysis == NULL) return EGADS_MALLOC;
      for (i = 0; i < problem->nAnalysis; i++) problem->analysis[i] = NULL;
      for (i = 0; i < problem->nAnalysis; i++) {
        stat = caps_readInitObj(&problem->analysis[i], ANALYSIS, NONE, NULL,
                                problem->mySelf);
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Error: Analysis %d caps_readInitObj = %d (caps_open)!\n",
                 i, stat);
          fclose(fptxt);
          return stat;
        }
      }
      for (i = 0; i < problem->nAnalysis; i++) {
        fscanf(fptxt, "%d %d %s", &nIn, &nOut, name);
        analysis = (capsAnalysis *) EG_alloc(sizeof(capsAnalysis));
        if (analysis == NULL) {
          fclose(fptxt);
          return EGADS_MALLOC;
        }
        problem->analysis[i]->name  = EG_strdup(name);
        problem->analysis[i]->blind = analysis;
        analysis->loadName          = NULL;
        analysis->fullPath          = NULL;
        analysis->path              = NULL;
        analysis->unitSys           = NULL;
        analysis->major             = CAPSMAJOR;
        analysis->minor             = CAPSMINOR;
        analysis->instStore         = NULL;
        analysis->intents           = NULL;
        analysis->autoexec          = 0;
        analysis->eFlag             = 0;
        analysis->nField            = 0;
        analysis->fields            = NULL;
        analysis->ranks             = NULL;
        analysis->fInOut            = NULL;
        analysis->nAnalysisIn       = nIn;
        analysis->analysisIn        = NULL;
        analysis->nAnalysisOut      = nOut;
        analysis->analysisOut       = NULL;
        analysis->nBody             = 0;
        analysis->bodies            = NULL;
        analysis->nTess             = 0;
        analysis->tess              = NULL;
        analysis->pre.nLines        = 0;
        analysis->pre.lines         = NULL;
        analysis->pre.pname         = NULL;
        analysis->pre.pID           = NULL;
        analysis->pre.user          = NULL;
        analysis->pre.sNum          = 0;
        analysis->info.magicnumber         = CAPSMAGIC;
        analysis->info.instance            = stat;
        analysis->info.problem             = problem;
        analysis->info.analysis            = analysis;
        analysis->info.pIndex              = 0;
        analysis->info.irow                = 0;
        analysis->info.icol                = 0;
        analysis->info.errs.nError         = 0;
        analysis->info.errs.errors         = NULL;
        analysis->info.wCntxt.aimWriterNum = 0;
        for (j = 0; j < 6; j++) analysis->pre.datetime[j] = 0;
        
        /* make the analysis value objects */
        if (nIn > 0) {
          analysis->analysisIn = (capsObject **)
                                 EG_alloc(nIn*sizeof(capsObject *));
          if (analysis->analysisIn == NULL) {
            fclose(fptxt);
            return EGADS_MALLOC;
          }
          for (j = 0; j < nIn; j++) analysis->analysisIn[j] = NULL;
          value = (capsValue *) EG_alloc(nIn*sizeof(capsValue));
          if (value == NULL) {
            fclose(fptxt);
            EG_free(analysis->analysisIn);
            analysis->analysisIn = NULL;
            return EGADS_MALLOC;
          }
          for (j = 0; j < nIn; j++) {
            value[j].length          = value[j].nrow = value[j].ncol = 1;
            value[j].type            = Integer;
            value[j].dim             = value[j].pIndex = 0;
            value[j].index           = j+1;
            value[j].lfixed          = value[j].sfixed = Fixed;
            value[j].nullVal         = NotAllowed;
            value[j].units           = NULL;
            value[j].meshWriter      = NULL;
            value[j].link            = NULL;
            value[j].vals.reals      = NULL;
            value[j].limits.dlims[0] = value[j].limits.dlims[1] = 0.0;
            value[j].linkMethod      = Copy;
            value[j].gInType         = 0;
            value[j].partial         = NULL;
            value[j].nderiv            = 0;
            value[j].derivs            = NULL;
          }
          for (j = 0; j < nIn; j++) {
            stat = caps_readInitObj(&analysis->analysisIn[j], VALUE, ANALYSISIN,
                                    NULL, problem->analysis[i]);
            if (stat != CAPS_SUCCESS) {
              printf(" CAPS Error: aIn %d caps_readInitObj = %d (caps_open)!\n",
                     j, stat);
              fclose(fptxt);
/*@-kepttrans@*/
              EG_free(value);
/*@+kepttrans@*/
              EG_free(analysis->analysisIn);
              analysis->analysisIn = NULL;
              return stat;
            }
/*@-immediatetrans@*/
            analysis->analysisIn[j]->blind = &value[j];
/*@+immediatetrans@*/
          }
          analysis->analysisIn[0]->blind = value;
        }
        if (nOut > 0) {
          analysis->analysisOut = (capsObject **)
                                  EG_alloc(nOut*sizeof(capsObject *));
          if (analysis->analysisOut == NULL) {
            fclose(fptxt);
            return EGADS_MALLOC;
          }
          for (j = 0; j < nOut; j++) analysis->analysisOut[j] = NULL;
          value = (capsValue *) EG_alloc(nIn*sizeof(capsValue));
          if (value == NULL) {
            fclose(fptxt);
            EG_free(analysis->analysisOut);
            analysis->analysisOut = NULL;
            return EGADS_MALLOC;
          }
          for (j = 0; j < nOut; j++) {
            value[j].length          = value[j].nrow = value[j].ncol = 1;
            value[j].type            = Integer;
            value[j].dim             = value[j].pIndex = 0;
            value[j].index           = j+1;
            value[j].lfixed          = value[j].sfixed = Fixed;
            value[j].nullVal         = NotAllowed;
            value[j].units           = NULL;
            value[j].meshWriter      = NULL;
            value[j].link            = NULL;
            value[j].vals.reals      = NULL;
            value[j].limits.dlims[0] = value[j].limits.dlims[1] = 0.0;
            value[j].linkMethod      = Copy;
            value[j].gInType         = 0;
            value[j].partial         = NULL;
            value[j].nderiv            = 0;
            value[j].derivs            = NULL;
          }
          for (j = 0; j < nOut; j++) {
            stat = caps_readInitObj(&analysis->analysisOut[j], VALUE, ANALYSISOUT,
                                    NULL, problem->analysis[i]);
            if (stat != CAPS_SUCCESS) {
              printf(" CAPS Error: aOut %d caps_readInitObj = %d (caps_open)!\n",
                     j, stat);
              fclose(fptxt);
/*@-kepttrans@*/
              EG_free(value);
/*@+kepttrans@*/
              EG_free(analysis->analysisOut);
              analysis->analysisOut = NULL;
              return stat;
            }
/*@-immediatetrans@*/
            analysis->analysisOut[j]->blind = &value[j];
/*@+immediatetrans@*/
          }
          analysis->analysisOut[0]->blind = value;
        }
      }
    }
    fclose(fptxt);
  }

#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\bound.txt", problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/bound.txt",   problem->root);
#endif
  fptxt = fopen(filename, "r");
  if (fptxt != NULL) {
    fscanf(fptxt, "%d %d", &problem->nBound, &problem->mBound);
    if (problem->nBound > 0) {
      problem->bounds = (capsObject **)
                        EG_alloc(problem->nBound*sizeof(capsObject *));
      if (problem->bounds == NULL) {
        fclose(fptxt);
        return EGADS_MALLOC;
      }
      for (i = 0; i < problem->nBound; i++) problem->bounds[i] = NULL;
      for (i = 0; i < problem->nBound; i++) {
        fscanf(fptxt, "%d %s", &j, name);
        bound = (capsBound *) EG_alloc(sizeof(capsBound));
        if (bound == NULL) {
          fclose(fptxt);
          return EGADS_MALLOC;
        }
        bound->dim        = 0;
        bound->state      = Empty;
        bound->lunits     = NULL;
        bound->plimits[0] = bound->plimits[1] = bound->plimits[2] =
                            bound->plimits[3] = 0.0;
        bound->geom       = NULL;
        bound->iBody      = 0;
        bound->iEnt       = 0;
        bound->curve      = NULL;
        bound->surface    = NULL;
        bound->index      = j;
        bound->nVertexSet = 0;
        bound->vertexSet  = NULL;

        stat = caps_readInitObj(&problem->bounds[i], BOUND, NONE, NULL,
                                problem->mySelf);
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Error: Bound %d caps_readInitObj = %d (caps_open)!\n",
                 i, stat);
          EG_free(bound);
          fclose(fptxt);
          return stat;
        }
        problem->bounds[i]->blind = bound;
        stat = caps_readInitVSets(problem, problem->bounds[i]);
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Error: Bound %d caps_readInitVSets = %d (caps_open)!\n",
                 i, stat);
          fclose(fptxt);
          return stat;
        }
      }
    }
    fclose(fptxt);
  }
  
  /* make all of the rest of the Objects */
  
  if (problem->nParam > 0) {
    problem->params = (capsObject **)
                      EG_alloc(problem->nParam*sizeof(capsObject *));
    if (problem->params == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nParam; i++) problem->params[i] = NULL;
    for (i = 0; i < problem->nParam; i++) {
      stat = caps_makeVal(Integer, 1, &i, &value);
      if (stat != CAPS_SUCCESS) return stat;
      stat = caps_readInitObj(&problem->params[i], VALUE, PARAMETER, NULL,
                              problem->mySelf);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: Param %d caps_readInitObj = %d (caps_open)!\n",
               i, stat);
/*@-kepttrans@*/
        EG_free(value);
/*@+kepttrans@*/
        return stat;
      }
/*@-kepttrans@*/
      problem->params[i]->blind = value;
/*@+kepttrans@*/
    }
  }
  
  if (problem->nGeomIn > 0) {
    problem->geomIn = (capsObject **)
                      EG_alloc(problem->nGeomIn*sizeof(capsObject *));
    if (problem->geomIn == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nGeomIn; i++) problem->geomIn[i] = NULL;
    value = (capsValue *) EG_alloc(problem->nGeomIn*sizeof(capsValue));
    if (value == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nGeomIn; i++) {
      value[i].nrow            = 1;
      value[i].ncol            = 1;
      value[i].type            = Double;
      value[i].dim             = Scalar;
      value[i].index           = i+1;
      value[i].pIndex          = 0;
      value[i].lfixed          = value[i].sfixed = Fixed;
      value[i].nullVal         = NotAllowed;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].length          = 1;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv            = 0;
      value[i].derivs            = NULL;
    }
    for (i = 0; i < problem->nGeomIn; i++) {
      stat = caps_readInitObj(&problem->geomIn[i], VALUE, GEOMETRYIN, NULL,
                              problem->mySelf);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: GeomIn %d caps_readInitObj = %d (caps_open)!\n",
               i, stat);
        EG_free(value);
        EG_free(problem->geomIn);
        problem->geomIn = NULL;
        return stat;
      }
/*@-immediatetrans@*/
      problem->geomIn[i]->blind = &value[i];
/*@+immediatetrans@*/
    }
  }
  
  if (problem->nGeomOut > 0) {
    problem->geomOut = (capsObject **)
                       EG_alloc(problem->nGeomOut*sizeof(capsObject *));
    if (problem->geomOut == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nGeomOut; i++) problem->geomOut[i] = NULL;
    value = (capsValue *) EG_alloc(problem->nGeomOut*sizeof(capsValue));
    if (value == NULL) return EGADS_MALLOC;
    for (i = 0; i < problem->nGeomOut; i++) {
      value[i].length          = 1;
      value[i].type            = DoubleDeriv;
      value[i].nrow            = 1;
      value[i].ncol            = 1;
      value[i].dim             = Scalar;
      value[i].index           = i+1;
      value[i].pIndex          = 0;
      value[i].lfixed          = value[i].sfixed = Change;
      value[i].nullVal         = IsNull;
      value[i].units           = NULL;
      value[i].meshWriter      = NULL;
      value[i].link            = NULL;
      value[i].vals.reals      = NULL;
      value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
      value[i].linkMethod      = Copy;
      value[i].gInType         = 0;
      value[i].partial         = NULL;
      value[i].nderiv            = 0;
      value[i].derivs            = NULL;
    }
    for (i = 0; i < problem->nGeomOut; i++) {
      stat = caps_readInitObj(&problem->geomOut[i], VALUE, GEOMETRYOUT, NULL,
                              problem->mySelf);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: GeomOut %d caps_readInitObj = %d (caps_open)!\n",
               i, stat);
        EG_free(value);
        EG_free(problem->geomOut);
        problem->geomOut = NULL;
        return stat;
      }
/*@-immediatetrans@*/
      problem->geomOut[i]->blind = &value[i];
/*@+immediatetrans@*/
    }
  }

  /* fill top-level objects */
  
  if (problem->params != NULL) {
    for (i = 0; i < problem->nParam; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\VP-%4.4d",
               problem->root, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/VP-%4.4d",
               problem->root, i+1);
#endif
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
        return CAPS_DIRERR;
      }
      stat = caps_readValue(fp, problem, problem->params[i]);
      fclose(fp);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: parameter %d/%d readValue = %d (caps_open)!\n",
               i+1, problem->nParam, stat);
        return stat;
      }
    }
  }
  
  if (problem->geomIn != NULL) {
    for (i = 0; i < problem->nGeomIn; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\VI-%4.4d",
               problem->root, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/VI-%4.4d",
               problem->root, i+1);
#endif
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
        return CAPS_DIRERR;
      }
      stat = caps_readValue(fp, problem, problem->geomIn[i]);
      fclose(fp);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: geomIn %d/%d readValue = %d (caps_open)!\n",
               i+1, problem->nGeomIn, stat);
        return stat;
      }
    }
  }
  
  if (problem->geomOut != NULL) {
    for (i = 0; i < problem->nGeomOut; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\VO-%4.4d",
               problem->root, i+1);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/VO-%4.4d",
               problem->root, i+1);
#endif
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        printf(" CAPS Error: Cannot open %s (caps_open)!\n", filename);
        return CAPS_DIRERR;
      }
      stat = caps_readValue(fp, problem, problem->geomOut[i]);
      fclose(fp);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: geomOut %d/%d readValue = %d (caps_open)!\n",
               i+1, problem->nGeomOut, stat);
        return stat;
      }
    }
  }
  
  if (problem->analysis != NULL) {
    for (i = 0; i < problem->nAnalysis; i++) {
      stat = caps_readAnalysis(problem, problem->analysis[i]);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: Analysis %d/%d readAnalysis = %d (caps_open)!\n",
               i+1, problem->nAnalysis, stat);
        return stat;
      }
    }
  }
  
  if (problem->bounds != NULL) {
    for (i = 0; i < problem->nBound; i++) {
      stat = caps_readBound(problem->bounds[i]);
      if (stat != CAPS_SUCCESS) {
        printf(" CAPS Error: bound index %d -- readBound = %d (caps_open)!\n",
               i+1, stat);
        return stat;
      }
    }
  }
  
  return CAPS_SUCCESS;
  
readerr:
  fclose(fp);
  return CAPS_IOERR;
}


int
caps_close(capsObject *pobject, int complete, /*@null@*/ const char *phName)
{
  int          i, j, stat, state, npts;
  char         path[PATH_MAX];
  ego          model, body;
  capsProblem  *problem;
  capsAnalysis *analysis;
  FILE         *fp;

  if (pobject == NULL)                   return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)          return CAPS_BADTYPE;
  if (pobject->blind == NULL)            return CAPS_NULLBLIND;
  stat    = caps_writeProblem(pobject);
  if (stat != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeProblem = %d (caps_close)!\n", stat);
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_CLOSE;
  if (problem->jrnl != NULL) fclose(problem->jrnl);

  if (complete == 1) {
#ifdef WIN32
    snprintf(path, PATH_MAX, "%s\\capsClosed", problem->root);
#else
    snprintf(path, PATH_MAX, "%s/capsClosed",  problem->root);
#endif
    fp = fopen(path, "w");
    if (fp == NULL) {
      printf(" CAPS Warning: Failed to open capsClosed!\n");
    } else {
      fclose(fp);
    }
  }
  
  if ((phName != NULL) && (complete == 1)) {
    if (problem->phName == NULL) {
#ifdef WIN32
      snprintf(path, PATH_MAX, "%s\\..\\%s", problem->root, phName);
#else
      snprintf(path, PATH_MAX, "%s/../%s",   problem->root, phName);
#endif
      caps_prunePath(path);
      stat = caps_statFile(path);
      if (stat != EGADS_NOTFOUND) {
        printf(" CAPS Warning: %s is not empty -- not renamed!\n", path);
      } else {
        stat = caps_rename(problem->root, path);
        if (stat != CAPS_SUCCESS) {
          printf(" CAPS Warning: Cannot rename %s!\n", path);
          return stat;
        }
      }
    } else {
      printf("CAPS Warning: New Phase Name not available for nonScratch!\n");
    }
  }

  /* remove the lock file */
  caps_rmLockOnClose(problem->root);

  if (problem->phName != NULL) EG_free(problem->phName);
  if (problem->lunits != NULL) {
    for (i = 0; i < problem->nBodies; i++)
      if (problem->lunits[i] != NULL) EG_free(problem->lunits[i]);
    EG_free(problem->lunits);
  }
  caps_freeFList(pobject);
  caps_freeOwner(&problem->writer);
  caps_freeOwner(&problem->geometry);

  /* deal with geometry */
  if (problem->modl != NULL) {
    if (pobject->subtype == PARAMETRIC) {
      if (problem->analysis != NULL) {
        for (i = 0; i < problem->nAnalysis; i++) {
          analysis = (capsAnalysis *) problem->analysis[i]->blind;
          if (analysis == NULL) continue;
          if (analysis->tess != NULL) {
            for (j = 0; j < analysis->nTess; j++)
              if (analysis->tess[j] != NULL) {
                /* delete bodies in the tessellation not on the OpenCSM stack */
                body = NULL;
                if (j >= analysis->nBody) {
                  stat = EG_statusTessBody(analysis->tess[j], &body, &state,
                                           &npts);
                  if (stat <  EGADS_SUCCESS) return stat;
                  if (stat == EGADS_OUTSIDE) return CAPS_SOURCEERR;
                }
                EG_deleteObject(analysis->tess[j]);
                if (body != NULL) EG_deleteObject(body);
                analysis->tess[j] = NULL;
              }
            EG_free(analysis->tess);
            analysis->tess   = NULL;
            analysis->nTess  = 0;
          }
        }
      }
      if (problem->stFlag != 1) {
        /* close up OpenCSM */
        ocsmFree(problem->modl);
      }
      if (problem->bodies != NULL) EG_free(problem->bodies);
    } else {
      if (problem->stFlag != 2) {
        model = (ego) problem->modl;
        EG_deleteObject(model);
      }
    }
  }

  /* deal with the CAPS Problem level Objects */
  caps_freeValueObjects(1, problem->nParam,   problem->params);
  caps_freeValueObjects(0, problem->nGeomIn,  problem->geomIn);
  caps_freeValueObjects(0, problem->nGeomOut, problem->geomOut);

  if (problem->bounds != NULL) {
    while (problem->nBound > 0) caps_delete(problem->bounds[0]);
    EG_free(problem->bounds);
  }

  if (problem->regGIN != NULL) {
    for (i = 0; i < problem->nRegGIN; i++) EG_free(problem->regGIN[i].name);
    EG_free(problem->regGIN);
  }

  /* close up the AIMs */
  if (problem->analysis != NULL) {
    for (i = 0; i < problem->nAnalysis; i++) {
      caps_freeFList(problem->analysis[i]);
      analysis = (capsAnalysis *) problem->analysis[i]->blind;
      caps_freeAnalysis(0, analysis);
      caps_freeOwner(&problem->analysis[i]->last);
      caps_freeAttrs(&problem->analysis[i]->attrs);
      problem->analysis[i]->magicnumber = 0;
      EG_free(problem->analysis[i]->name);
      problem->analysis[i]->name = NULL;
      EG_free(problem->analysis[i]);
    }
    EG_free(problem->analysis);
  }

  /* close up the AIM subsystem */
  aim_cleanupAll(&problem->aimFPTR);

  /* close up EGADS and free the problem */
  if (problem->root     != NULL) EG_free(problem->root);
  if ((problem->context != NULL) && (problem->stFlag != 2))
    EG_close(problem->context);
  EG_free(problem->execs);
  EG_free(problem);

  /* cleanup and invalidate the object */
  caps_freeAttrs(&pobject->attrs);
  caps_freeOwner(&pobject->last);
  pobject->magicnumber = 0;
  EG_free(pobject->name);
  pobject->name = NULL;
  EG_free(pobject);

  return CAPS_SUCCESS;
}


int
caps_phaseState(const char *prPath, /*@null@*/ const char *phName, int *bFlag)
{
  int        i, n, len, status;
  char       root[PATH_MAX], current[PATH_MAX];
  const char *prName;

  *bFlag = 0;
  if (prPath == NULL) return CAPS_NULLNAME;
  if (phName != NULL) {
    if (strcmp(phName, "Scratch") == 0) {
      printf(" CAPS Error: Cannot use the phase Scratch!\n");
      return CAPS_BADNAME;
    }
    len = strlen(phName);
    for (i = 0; i < len; i++)
      if ((phName[i] == '/') || (phName[i] == '\\')) {
        printf(" CAPS Error: Cannot use slashes in phase name: %s\n", phName);
        return CAPS_BADNAME;
      }
  }
  
  /* set up our path and name */
  len = strlen(prPath);
  for (i = len-1; i > 0; i--)
    if ((prPath[i] == '/') || (prPath[i] == '\\')) break;
  if (i != 0) i++;
  prName = &prPath[i];
  status = caps_isNameOK(prName);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: %s is not a valid Problem Name!\n", prName);
    return status;
  }
  n = -1;
#ifdef WIN32
  /* do we have a Windows drive? */
  if (prPath[1] == ':') {
    int drive, oldrive;
    if (prPath[0] >= 97) {
      drive = prPath[0] - 96;
    } else {
      drive = prPath[0] - 64;
    }
    oldrive = _getdrive();
    status  = _chdrive(drive);
    if (status == -1) {
      printf(" CAPS Error: Cannot change drive to %c!\n", prPath[0]);
      return CAPS_DIRERR;
    }
    (void) _chdrive(oldrive);
    n = 0;
  }
#endif
  if (n == -1)
    if ((prPath[0] == '/') || (prPath[0] == '\\')) n = 0;
  if (n ==  0) {
    /* absolute path */
    status = caps_statFile(prPath);
    if (status == EGADS_SUCCESS) {
      printf(" CAPS Info: %s lands on a flat file!\n", prPath);
      return CAPS_DIRERR;
    } else if (status == EGADS_NOTFOUND) {
      printf(" CAPS Info: Path %s does not exist\n", prPath);
      return status;
    }
#ifdef WIN32
    if (prPath[1] == ':') {
      if (phName == NULL) {
        snprintf(root, PATH_MAX, "%s\\Scratch", prPath);
      } else {
        snprintf(root, PATH_MAX, "%s\\%s", prPath, phName);
      }
    } else {
      if (phName == NULL) {
        snprintf(root, PATH_MAX, "%c:%s\\Scratch", _getdrive()+64, prPath);
      } else {
        snprintf(root, PATH_MAX, "%c:%s\\%s", _getdrive()+64, prPath, phName);
      }
    }
#else
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s/Scratch", prPath);
    } else {
      snprintf(root, PATH_MAX, "%s/%s", prPath, phName);
    }
#endif
  } else {
    /* relative path -- make it absolute */
    (void) getcwd(current, PATH_MAX);
#ifdef WIN32
    snprintf(root, PATH_MAX, "%s\\%s", current, prPath);
#else
    snprintf(root, PATH_MAX, "%s/%s",  current, prPath);
#endif
    status = caps_statFile(root);
    if (status == EGADS_SUCCESS) {
      printf(" CAPS Info: Path %s lands on a flat file\n", root);
      return status;
    } else if (status == EGADS_NOTFOUND) {
      printf(" CAPS Info: Path %s does not exist\n", root);
      return status;
    }
#ifdef WIN32
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s\\%s\\Scratch", current, prPath);
    } else {
      snprintf(root, PATH_MAX, "%s\\%s\\%s", current, prPath, phName);
    }
#else
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s/%s/Scratch", current, prPath);
    } else {
      snprintf(root, PATH_MAX, "%s/%s/%s", current, prPath, phName);
    }
#endif
  }
  caps_prunePath(root);
  status = caps_statFile(root);
  if (status == EGADS_SUCCESS) {
    printf(" CAPS Info: Path %s lands on a flat file\n", root);
    return status;
  } else if (status == EGADS_NOTFOUND) {
    printf(" CAPS Info: Path %s does not exist\n", root);
    return status;
  }
  
#ifdef WIN32
  snprintf(current, PATH_MAX, "%s\\capsLock", root);
#else
  snprintf(current, PATH_MAX, "%s/capsLock",  root);
#endif
  status = caps_statFile(current);
  if (status == EGADS_SUCCESS) *bFlag += 1;

#ifdef WIN32
  snprintf(current, PATH_MAX, "%s\\capsClosed", root);
#else
  snprintf(current, PATH_MAX, "%s/capsClosed",  root);
#endif
  status = caps_statFile(current);
  if (status == EGADS_SUCCESS) *bFlag += 2;
  
  return CAPS_SUCCESS;
}


int
caps_journalState(const capsObject *pobject)
{
  capsProblem *problem;
  
  if (pobject == NULL)                   return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)          return CAPS_BADTYPE;
  if (pobject->blind == NULL)            return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  
  return problem->stFlag;
}


int
caps_open(const char *prPath, /*@null@*/ const char *phName, int flag,
          void *ptr, int outLevel, capsObject **pobject,
          int *nErr, capsErrs **errors)
{
  int           i, j, k, n, len, status, oclass, mtype, *senses, idot = 0;
  int           type, nattr, ngIn, ngOut, buildTo, builtTo, ibody;
  int           nrow, ncol, nbrch, npmtr, nbody, ivec[2];
  char          byte, *units, *env, root[PATH_MAX], current[PATH_MAX];
  char          name[MAX_NAME_LEN];
  char          filename[PATH_MAX], temp[PATH_MAX], **tmp, line[129];
  CAPSLONG      fileLen, ret, sNum;
  double        dot, lower, upper, data[4], *reals;
  ego           model, ref, *childs;
  capsObject    *object, *aobject, *objs;
  capsProblem   *problem;
  capsAnalysis  *analysis;
  capsBound     *bound;
  capsVertexSet *vertexset;
  capsDataSet   *dataset;
  capsValue     *value, *valIn = NULL;
  modl_T        *MODL = NULL;
  size_t        nw;
  FILE          *fp;
  const char    *aname, *astring, *prName, *fname = NULL;
  const int     *aints;
  const double  *areals;

  if (nErr   == NULL) return CAPS_NULLVALUE;
  if (errors == NULL) return CAPS_NULLVALUE;
  *pobject = NULL;
  *nErr    = 0;
  *errors  = NULL;
  if (prPath == NULL) return CAPS_NULLNAME;
  len = strlen(prPath);
  for (i = 0; i < len; i++)
    if (prPath[i] == ' ') {
      caps_makeSimpleErr(NULL, CERROR, "Cannot use spaces in path (caps_open):",
                         prPath, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_BADNAME;
    }
  if (phName != NULL) {
    if (strcmp(phName, "Scratch") == 0) {
      caps_makeSimpleErr(NULL, CERROR, "Cannot use the phase Scratch (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_BADNAME;
    }
    len = strlen(phName);
    for (i = 0; i < len; i++)
      if ((phName[i] == '/') || (phName[i] == '\\')) {
        caps_makeSimpleErr(NULL, CERROR,
                           "Cannot use slashes in phase name (caps_open):",
                           phName, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return CAPS_BADNAME;
      }
  }
  for (i = 0; i < PATH_MAX-1; i++) filename[i] = ' ';
  filename[PATH_MAX-1] = 0;

  if (flag == 0) {
    /* filename input */
    fname = (const char *) ptr;
    if (fname == NULL) return CAPS_NULLNAME;
    n = -1;
#ifdef WIN32
    /* do we have a Windows drive? */
    if (fname[1] == ':') {
      int id, oldrive;
      if (fname[0] >= 97) {
        id = fname[0] - 96;
      } else {
        id = fname[0] - 64;
      }
      oldrive = _getdrive();
      status  = _chdrive(id);
      if (status == -1) {
        snprintf(temp, PATH_MAX, "Cannot change drive to %c (caps_open)!",
                 fname[0]);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return CAPS_DIRERR;
      }
      (void) _chdrive(oldrive);
      n = 0;
    }
#endif
    if (n == -1)
      if ((fname[0] == '/') || (fname[0] == '\\')) n = 0;
    if (n ==  0) {
      /* absolute path */
#ifdef WIN32
      if (fname[1] == ':') {
        snprintf(filename, PATH_MAX, "%s", fname);
      } else {
        snprintf(filename, PATH_MAX, "%c:%s", _getdrive()+64, fname);
      }
#else
      snprintf(filename, PATH_MAX, "%s", fname);
#endif
    } else {
      /* relative path -- make it absolute */
      (void) getcwd(current, PATH_MAX);
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\%s", current, fname);
#else
      snprintf(filename, PATH_MAX, "%s/%s",  current, fname);
#endif
    }
    caps_prunePath(filename);
  } else if (flag == 1) {
    MODL   = (modl_T *) ptr;
  } else if (flag == 2) {
    model  = (ego) ptr;
    status = EG_getTopology(model, &ref, &oclass, &mtype, data, &len, &childs,
                            &senses);
    if (status != EGADS_SUCCESS) return status;
  } else if (flag == 3) {
    if (phName == NULL) {
      caps_makeSimpleErr(NULL, CERROR,
                         "Cannot start with a NULL PhaseName (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    }
    fname = (const char *) ptr;
    if (fname == NULL) return CAPS_NULLNAME;
  } else if (flag != 4) {
    /* other values of flag */
    return CAPS_NOTIMPLEMENT;
  }

  /* does file exist? */
  fileLen = 0;
  if (flag == 0) {
    fp = fopen(filename, "rb");
    if (fp == NULL) return CAPS_NOTFOUND;
    do {
      ret = fread(&byte, sizeof(char), 1, fp);
      fileLen++;
    } while (ret != 0);
    fclose(fp);
    if (fileLen >  0) fileLen--;
    if (fileLen == 0) return CAPS_BADVALUE;

    /* find the file extension */
    len = strlen(filename);
    for (idot = len-1; idot > 0; idot--)
      if (filename[idot] == '.') break;
    if (idot == 0) return CAPS_BADNAME;
/*@-unrecog@*/
    if ((strcasecmp(&filename[idot],".csm")   != 0) &&
        (strcasecmp(&filename[idot],".egads") != 0)) return CAPS_BADTYPE;
/*@+unrecog@*/
  }
  
  /* set up our path and name */
  len = strlen(prPath);
  for (i = len-1; i > 0; i--)
    if ((prPath[i] == '/') || (prPath[i] == '\\')) break;
  if (i != 0) i++;
  prName = &prPath[i];
  status = caps_isNameOK(prName);
  if (status != CAPS_SUCCESS) {
    caps_makeSimpleErr(NULL, CERROR, "Not a valid Problem Name (caps_open):",
                       prName, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return status;
  }
  n = -1;
#ifdef WIN32
  /* do we have a Windows drive? */
  if (prPath[1] == ':') {
    int drive, oldrive;
    if (prPath[0] >= 97) {
      drive = prPath[0] - 96;
    } else {
      drive = prPath[0] - 64;
    }
    oldrive = _getdrive();
    status  = _chdrive(drive);
    if (status == -1) {
      snprintf(temp, PATH_MAX, "Cannot change drive to %c (caps_open)!",
               prPath[0]);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    }
    (void) _chdrive(oldrive);
    n = 0;
  }
#endif
  if (n == -1)
    if ((prPath[0] == '/') || (prPath[0] == '\\')) n = 0;
  if (n ==  0) {
    /* absolute path */
    status = caps_statFile(prPath);
    if (status == EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Lands on a flat file (caps_open):",
                         prPath, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    } else if (status == EGADS_NOTFOUND) {
      status = caps_mkDir(prPath);
      if (status != EGADS_SUCCESS) {
        caps_makeSimpleErr(NULL, CERROR, "Cannot mkDir (caps_open):",
                           prPath, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return status;
      }
    }
#ifdef WIN32
    if (prPath[1] == ':') {
      if (phName == NULL) {
        snprintf(root, PATH_MAX, "%s\\Scratch", prPath);
      } else {
        snprintf(root, PATH_MAX, "%s\\%s", prPath, phName);
        if (flag == 3)
          snprintf(filename, PATH_MAX, "%s\\%s", prPath, fname);
      }
    } else {
      if (phName == NULL) {
        snprintf(root, PATH_MAX, "%c:%s\\Scratch", _getdrive()+64, prPath);
      } else {
        snprintf(root, PATH_MAX, "%c:%s\\%s", _getdrive()+64, prPath, phName);
        if (flag == 3)
          snprintf(filename, PATH_MAX, "%c:%s\\%s", _getdrive()+64, prPath,
                   fname);
      }
    }
#else
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s/Scratch", prPath);
    } else {
      snprintf(root, PATH_MAX, "%s/%s", prPath, phName);
      if (flag == 3)
        snprintf(filename, PATH_MAX, "%s/%s", prPath, fname);
    }
#endif
  } else {
    /* relative path -- make it absolute */
    (void) getcwd(current, PATH_MAX);
#ifdef WIN32
    snprintf(root, PATH_MAX, "%s\\%s", current, prPath);
#else
    snprintf(root, PATH_MAX, "%s/%s",  current, prPath);
#endif
    status = caps_statFile(root);
    if (status == EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Path lands on a flat file (caps_open):",
                         root, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    } else if (status == EGADS_NOTFOUND) {
      status = caps_mkDir(root);
      if (status != EGADS_SUCCESS) {
        caps_makeSimpleErr(NULL, CERROR, "Cannot make Path (caps_open):",
                           root, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return status;
      }
    }
#ifdef WIN32
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s\\%s\\Scratch", current, prPath);
    } else {
      snprintf(root, PATH_MAX, "%s\\%s\\%s", current, prPath, phName);
      if (flag == 3)
        snprintf(filename, PATH_MAX, "%s\\%s\\%s", current, prPath, fname);
    }
#else
    if (phName == NULL) {
      snprintf(root, PATH_MAX, "%s/%s/Scratch", current, prPath);
    } else {
      snprintf(root, PATH_MAX, "%s/%s/%s", current, prPath, phName);
      if (flag == 3)
        snprintf(filename, PATH_MAX, "%s/%s/%s", current, prPath, fname);
    }
#endif
  }
  caps_prunePath(root);
  /* not a continuation -- must not have a directory (unless Scratch)! */
  if (flag == 3) {
    caps_prunePath(filename);
    nw = strlen(filename)+1;
    memcpy(temp, filename, nw);
    nw = strlen(root)+1;
    memcpy(filename, root, nw);
    nw = strlen(temp)+1;
    memcpy(root, temp, nw);
    status = caps_statFile(root);
    if (status == EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Lands on a flat file (caps_open):",
                         root, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    } else if (status != EGADS_NOTFOUND) {
      caps_makeSimpleErr(NULL, CERROR, "Path already exists (caps_open):",
                         root, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return EGADS_EXISTS;
    }
#ifdef WIN32
    snprintf(current, PATH_MAX, "%s\\capsClosed", filename);
#else
    snprintf(current, PATH_MAX, "%s/capsClosed",  filename);
#endif
    status = caps_statFile(current);
    if (status != EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Not closed (caps_open):",
                         filename, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    }
    status = caps_cpDir(filename, root);
    if (status != EGADS_SUCCESS) {
      snprintf(temp, PATH_MAX, "Copy directory = %d (caps_open)", status);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return status;
    }
    /* remove closed & lock markers */
#ifdef WIN32
    snprintf(current, PATH_MAX, "%s\\capsClosed", root);
#else
    snprintf(current, PATH_MAX, "%s/capsClosed",  root);
#endif
    status = caps_rmFile(current);
    if ((status != EGADS_SUCCESS) && (status != EGADS_NOTFOUND))
      printf(" CAPS Warning: Cannot remove Closed file!\n");
#ifdef WIN32
    snprintf(current, PATH_MAX, "%s\\capsLock", root);
#else
    snprintf(current, PATH_MAX, "%s/capsLock",  root);
#endif
    status = caps_rmFile(current);
    if ((status != EGADS_SUCCESS) && (status != EGADS_NOTFOUND))
      printf(" CAPS Warning: Cannot remove Lock file (caps_open)\n");
  } else if (flag == 4) {
    /* should not be closed */
#ifdef WIN32
    snprintf(current, PATH_MAX, "%s\\capsClosed", root);
#else
    snprintf(current, PATH_MAX, "%s/capsClosed",  root);
#endif
    status = caps_statFile(current);
    if (status != EGADS_NOTFOUND) {
      caps_makeSimpleErr(NULL, CERROR,
                         "Found Closed file on continuation (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_EXISTS;
    }
  } else {
    status = caps_statFile(root);
    if (status == EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Lands on a flat file (caps_open):",
                         root, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return CAPS_DIRERR;
    } else if (status != EGADS_NOTFOUND) {
      if (phName != NULL) {
        caps_makeSimpleErr(NULL, CERROR, "Already exists (caps_open):",
                           root, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return EGADS_EXISTS;
      } else {
#ifdef WIN32
        snprintf(current, PATH_MAX, "%s\\capsLock", root);
#else
        snprintf(current, PATH_MAX, "%s/capsLock",  root);
#endif
        if (caps_statFile(current) != EGADS_NOTFOUND) {
          caps_makeSimpleErr(NULL, CERROR, "Lock file exists (caps_open)!",
                      "If you believe that no one else is running here remove:",
                             current, errors);
          if (*errors != NULL) *nErr = (*errors)->nError;
          return CAPS_DIRERR;
        }
        caps_rmDir(root);
      }
    }
    status = caps_mkDir(root);
    if (status != EGADS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Cannot mkDir (caps_open)!",
                         root, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      return status;
    }
  }
  
  /* do we have a lock file? */
#ifdef WIN32
  snprintf(current, PATH_MAX, "%s\\capsLock", root);
#else
  snprintf(current, PATH_MAX, "%s/capsLock",  root);
#endif
  if (caps_statFile(current) != EGADS_NOTFOUND) {
    caps_makeSimpleErr(NULL, CERROR, "Lock file exists (caps_open)!",
                      "If you believe that no one else is running here remove:",
                       current, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_DIRERR;
  }
  tmp = (char **) EG_reall(CAPSlocks, (CAPSnLock+1)*sizeof(char *));
  if (tmp == NULL) {
    caps_makeSimpleErr(NULL, CERROR, "ReAllocating Locks storage (caps_open)!",
                       NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return EGADS_MALLOC;
  }
  CAPSlocks = tmp;
  CAPSlocks[CAPSnLock] = EG_strdup(current);
  CAPSnLock++;
  fp = fopen(current, "w");
  if (fp == NULL) {
    caps_makeSimpleErr(NULL, CERROR, "Cannot open Lock file (caps_open)!",
                       NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    return CAPS_DIRERR;
  }
  fclose(fp);
  caps_initSignals();

  /* check the outLevel env (mostly used for testing) */
  env = getenv("CAPS_OUTLEVEL");
  if (env != NULL && (env[0] == '0' ||
                      env[0] == '1' ||
                      env[0] == '2')) {
    outLevel = atoi(env);
  }

  /* we are OK -- make the Problem */

  problem = (capsProblem *) EG_alloc(sizeof(capsProblem));
  if (problem == NULL) return EGADS_MALLOC;

  /* initialize all members */
  problem->signature       = NULL;
  problem->context         = NULL;
  problem->utsystem        = NULL;
  problem->root            = EG_strdup(root);
  problem->phName          = NULL;
  problem->stFlag          = flag;
  problem->jrnl            = NULL;
  problem->outLevel        = outLevel;
  problem->funID           = CAPS_OPEN;
  problem->modl            = NULL;
  problem->nParam          = 0;
  problem->params          = NULL;
  problem->nGeomIn         = 0;
  problem->geomIn          = NULL;
  problem->nGeomOut        = 0;
  problem->geomOut         = NULL;
  problem->nAnalysis       = 0;
  problem->analysis        = NULL;
  problem->mBound          = 0;
  problem->nBound          = 0;
  problem->bounds          = NULL;
  problem->geometry.nLines = 0;
  problem->geometry.lines  = NULL;
  problem->geometry.pname  = NULL;
  problem->geometry.pID    = NULL;
  problem->geometry.user   = NULL;
  problem->geometry.sNum   = 0;
  for (j = 0; j < 6; j++) problem->geometry.datetime[j] = 0;
  problem->nBodies         = 0;
  problem->bodies          = NULL;
  problem->lunits          = NULL;
  problem->nEGADSmdl       = 0;
  problem->nRegGIN         = 0;
  problem->regGIN          = NULL;
  problem->nExec           = 0;
  problem->execs           = NULL;
  problem->sNum            = problem->writer.sNum = 1;
  problem->jpos            = 0;
  problem->writer.nLines   = 0;
  problem->writer.lines    = NULL;
  problem->writer.pname    = EG_strdup(prName);
  caps_getStaticStrings(&problem->signature, &problem->writer.pID,
                        &problem->writer.user);
  for (j = 0; j < 6; j++) problem->writer.datetime[j] = 0;
  problem->aimFPTR.aim_nAnal = 0;
  for (j = 0; j < MAXANAL; j++) {
    problem->aimFPTR.aimName[j]     = NULL;
    problem->aimFPTR.aimDLL[j]      = NULL;
    problem->aimFPTR.aimInit[j]     = NULL;
    problem->aimFPTR.aimDiscr[j]    = NULL;
    problem->aimFPTR.aimFreeD[j]    = NULL;
    problem->aimFPTR.aimLoc[j]      = NULL;
    problem->aimFPTR.aimInput[j]    = NULL;
    problem->aimFPTR.aimPAnal[j]    = NULL;
    problem->aimFPTR.aimPost[j]     = NULL;
    problem->aimFPTR.aimOutput[j]   = NULL;
    problem->aimFPTR.aimCalc[j]     = NULL;
    problem->aimFPTR.aimXfer[j]     = NULL;
    problem->aimFPTR.aimIntrp[j]    = NULL;
    problem->aimFPTR.aimIntrpBar[j] = NULL;
    problem->aimFPTR.aimIntgr[j]    = NULL;
    problem->aimFPTR.aimIntgrBar[j] = NULL;
    problem->aimFPTR.aimBdoor[j]    = NULL;
    problem->aimFPTR.aimClean[j]    = NULL;
  }
  if (phName != NULL) problem->phName = EG_strdup(phName);

  status = caps_makeObject(&object);
  if (status != CAPS_SUCCESS) {
    EG_free(problem);
    return status;
  }
  problem->mySelf = object;
  object->type    = PROBLEM;
  object->blind   = problem;

  problem->utsystem = caps_initUnits();
  if (problem->utsystem == NULL) {
    caps_close(object, 0, NULL);
    return CAPS_UNITERR;
  }

  /* get EGADS context or open up EGADS */
  if (flag == 2) {
    status = EG_getContext(model, &problem->context);
  } else {
    status = EG_open(&problem->context);
  }
  if (status != EGADS_SUCCESS) {
    caps_close(object, 0, NULL);
    return status;
  }
  if (problem->context == NULL) {
    caps_close(object, 0, NULL);
    return EGADS_NOTCNTX;
  }

  /* load the CAPS state */
  if (flag > 2) {
    status = caps_readState(object);
    if (status != CAPS_SUCCESS) {
      caps_close(object, 0, NULL);
      return status;
    }
  
    /* reload the Geometry */
    if (object->subtype == PARAMETRIC) {
    
#ifdef WIN32
      snprintf(current, PATH_MAX, "%s\\capsRestart.cpc", root);
#else
      snprintf(current, PATH_MAX, "%s/capsRestart.cpc",  root);
#endif
      if (problem->outLevel != 1) ocsmSetOutLevel(problem->outLevel);
      status = ocsmLoad(current, &problem->modl);
      if (status < SUCCESS) {
        printf(" CAPS Error: Cannot ocsmLoad %s (caps_open)!\n", current);
        caps_close(object, 0, NULL);
        return status;
      }
      MODL = (modl_T *) problem->modl;
      if (MODL == NULL) {
        caps_makeSimpleErr(NULL, CERROR, "Cannot get OpenCSM MODL (caps_open)!",
                           NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return CAPS_NOTFOUND;
      }
      MODL->context   = problem->context;
/*@-kepttrans@*/
      MODL->userdata  = problem;
/*@+kepttrans@*/
      MODL->tessAtEnd = 0;
      status = ocsmRegSizeCB(problem->modl, caps_sizeCB);
      if (status != SUCCESS)
        printf(" CAPS Warning: ocsmRegSizeCB = %d (caps_open)!", status);
      /* check that Branches are properly ordered */
      status = ocsmCheck(problem->modl);
      if (status < SUCCESS) {
        snprintf(temp, PATH_MAX, "ocsmCheck = %d (caps_open)!", status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return status;
      }
      fflush(stdout);
      /* reset the GeomIns */
      if (problem->geomIn != NULL)
        for (i = 0; i < problem->nGeomIn; i++) {
          if (problem->geomIn[i] == NULL) continue;
          value = (capsValue *) problem->geomIn[i]->blind;
          reals = value->vals.reals;
          if (value->length == 1) reals = &value->vals.real;
          status = ocsmGetPmtr(problem->modl, value->pIndex, &type,
                               &nrow, &ncol, name);
          if (status != SUCCESS) {
            snprintf(temp, PATH_MAX, "ocsmGetPmtr %d fails with %d (caps_open)!",
                     value->pIndex, status);
            caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
            if (*errors != NULL) *nErr = (*errors)->nError;
            caps_close(object, 0, NULL);
            return status;
          }
          if ((ncol != value->ncol) || (nrow != value->nrow)) {
            snprintf(temp, PATH_MAX, "%s ncol = %d %d, nrow = %d %d (caps_open)!",
                     name, ncol, value->ncol, nrow, value->nrow);
            caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
            if (*errors != NULL) *nErr = (*errors)->nError;
            caps_close(object, 0, NULL);
            return CAPS_MISMATCH;
          }
          /* flip storage
          for (n = j = 0; j < ncol; j++)
            for (k = 0; k < nrow; k++, n++) {  */
          for (n = k = 0; k < nrow; k++)
            for (j = 0; j < ncol; j++, n++) {
              status = ocsmSetValuD(problem->modl, value->pIndex,
                                    k+1, j+1, reals[n]);
              if (status != SUCCESS) {
                snprintf(temp, PATH_MAX, "%d ocsmSetValuD[%d,%d] fails with %d!",
                         value->pIndex, k+1, j+1, status);
                caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
                if (*errors != NULL) *nErr = (*errors)->nError;
                caps_close(object, 0, NULL);
                return status;
              }
            }
        }

      buildTo = 0;                          /* all */
      nbody   = 0;
      status  = ocsmBuild(problem->modl, buildTo, &builtTo, &nbody, NULL);
      fflush(stdout);
      if (status != SUCCESS) {
        snprintf(temp, PATH_MAX, "ocsmBuild to %d fails with %d (caps_open)!",
                 builtTo, status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return status;
      }
      nbody = 0;
      for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;
        if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
        nbody++;
      }

      if (nbody > 0) {
        problem->bodies = (ego *)   EG_alloc(nbody*sizeof(ego));
        problem->lunits = (char **) EG_alloc(nbody*sizeof(char *));
        if ((problem->bodies != NULL) && (problem->lunits != NULL)) {
          problem->nBodies = nbody;
          i = 0;
          for (ibody = 1; ibody <= MODL->nbody; ibody++) {
            if (MODL->body[ibody].onstack != 1) continue;
            if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
            problem->bodies[i] = MODL->body[ibody].ebody;
            caps_fillLengthUnits(problem, problem->bodies[i],
                                 &problem->lunits[i]);
            i++;
          }
        } else {
          if (problem->lunits != NULL) EG_free(problem->lunits);
          problem->lunits = NULL;
          snprintf(temp, PATH_MAX, "Malloc on %d Body (caps_open)!\n", nbody);
          caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
          if (*errors != NULL) *nErr = (*errors)->nError;
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
      }
      status = ocsmInfo(problem->modl, &nbrch, &npmtr, &nbody);
      if (status != SUCCESS) {
        caps_close(object, 0, NULL);
        snprintf(temp, PATH_MAX, "ocsmInfo returns %d (caps_open)!", status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        return status;
      }
      for (ngIn = ngOut = i = 0; i < npmtr; i++) {
        status = ocsmGetPmtr(problem->modl, i+1, &type, &nrow, &ncol, name);
        if (status != SUCCESS) {
          caps_close(object, 0, NULL);
          return status;
        }
        if (type == OCSM_OUTPMTR) ngOut++;
        if (type == OCSM_DESPMTR) ngIn++;
        if (type == OCSM_CFGPMTR) ngIn++;
        if (type == OCSM_CONPMTR) ngIn++;
      }
      if (ngIn != problem->nGeomIn) {
        snprintf(temp, PATH_MAX, "# Design Vars = %d -- from %s = %d (caps_open)!",
                 ngIn, filename, problem->nGeomIn);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return CAPS_MISMATCH;
      }
      if (ngOut != problem->nGeomOut) {
        snprintf(temp, PATH_MAX, "# Geometry Outs = %d -- from %s = %d (caps_open)!",
                 ngOut, filename, problem->nGeomOut);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return CAPS_MISMATCH;
      }
      /* check geomOuts */
      if (problem->geomOut != NULL)
        for (i = j = 0; j < npmtr; j++) {
          ocsmGetPmtr(problem->modl, j+1, &type, &nrow, &ncol, name);
          if (type != OCSM_OUTPMTR) continue;
          if (strcmp(name, problem->geomOut[i]->name) != 0) {
            snprintf(temp, PATH_MAX, "%d Geometry Outs %s != %s (caps_open)!",
                   i+1, name, problem->geomOut[i]->name);
            caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
            if (*errors != NULL) *nErr = (*errors)->nError;
            caps_close(object, 0, NULL);
            return CAPS_MISMATCH;
          }
          i++;
        }

    } else {
      
      /* Problem is static */
      if (problem->outLevel != 1)
        EG_setOutLevel(problem->context, problem->outLevel);
#ifdef WIN32
      snprintf(current, PATH_MAX, "%s\\capsRestart.egads", root);
#else
      snprintf(current, PATH_MAX, "%s/capsRestart.egads",  root);
#endif
      status = EG_loadModel(problem->context, 1, current, &model);
      if (status != EGADS_SUCCESS) {
        snprintf(temp, PATH_MAX, "%s EG_loadModel = %d (caps_open)!",
                 current, status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return status;
      }
      problem->modl = model;
      status = EG_getTopology(model, &ref, &oclass, &mtype, data,
                              &problem->nBodies, &problem->bodies, &senses);
      if (status != EGADS_SUCCESS) {
        snprintf(temp, PATH_MAX, "%s EG_getTopology = %d (caps_open)!",
                 current, status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return status;
      }
      nbody = problem->nBodies;
      /* get length units if exist */
      if (problem->nBodies > 0) {
        problem->lunits = (char **) EG_alloc(problem->nBodies*sizeof(char *));
        if ((problem->lunits != NULL) && (problem->bodies != NULL))
          for (i = 0; i < problem->nBodies; i++)
            caps_fillLengthUnits(problem, problem->bodies[i],
                                 &problem->lunits[i]);
      }
    }
    
    /* set the bodies for the AIMs */
    if ((nbody > 0) && (problem->bodies != NULL) &&
        (problem->analysis != NULL))
      for (i = 0; i < problem->nAnalysis; i++) {
        analysis = (capsAnalysis *) problem->analysis[i]->blind;
        if (analysis == NULL) continue;
        status   = caps_filter(problem, analysis);
        if (status != CAPS_SUCCESS)
          printf(" CAPS Warning: %s caps_filter = %d (caps_open)!\n",
                 problem->analysis[i]->name, status);
      }

    /* get the capsDiscr structures */
    if (problem->bounds != NULL)
      for (i = 0; i < problem->nBound; i++) {
        if (problem->bounds[i]              == NULL)      continue;
        if (problem->bounds[i]->magicnumber != CAPSMAGIC) continue;
        if (problem->bounds[i]->type        != BOUND)     continue;
        if (problem->bounds[i]->blind       == NULL)      continue;
        bound = (capsBound *) problem->bounds[i]->blind;
        for (j = 0; j < bound->nVertexSet; j++) {
          if (bound->vertexSet[j]              == NULL)      continue;
          if (bound->vertexSet[j]->magicnumber != CAPSMAGIC) continue;
          if (bound->vertexSet[j]->type        != VERTEXSET) continue;
          if (bound->vertexSet[j]->blind       == NULL)      continue;
          vertexset = (capsVertexSet *) bound->vertexSet[j]->blind;
          if (vertexset->analysis != NULL)
            if (vertexset->analysis->blind != NULL) {
              analysis = (capsAnalysis *) vertexset->analysis->blind;
              vertexset->discr->dim       = bound->dim;
              vertexset->discr->instStore = analysis->instStore;
              status = aim_Discr(problem->aimFPTR, analysis->loadName,
                                 problem->bounds[i]->name, vertexset->discr);
              if (status != CAPS_SUCCESS) {
                aim_FreeDiscr(problem->aimFPTR, analysis->loadName,
                              vertexset->discr);
                snprintf(temp, PATH_MAX, "Bound = %s, Analysis = %s aimDiscr = %d",
                         problem->bounds[i]->name, analysis->loadName, status);
                caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
                if (*errors != NULL) *nErr = (*errors)->nError;
                caps_close(object, 0, NULL);
                return status;
              } else {
                /* check the validity of the discretization just returned */
                status = caps_checkDiscr(vertexset->discr, 129, line);
                if (status != CAPS_SUCCESS) {
                  snprintf(temp, PATH_MAX, "Bound = %s, Analysis = %s chkDiscr=%d",
                           problem->bounds[i]->name, analysis->loadName, status);
                  caps_makeSimpleErr(NULL, CERROR, temp, line, NULL, errors);
                  if (*errors != NULL) *nErr = (*errors)->nError;
                  aim_FreeDiscr(problem->aimFPTR, analysis->loadName,
                                vertexset->discr);
                  caps_close(object, 0, NULL);
                  return status;
                }
              }
              /* do the counts match? */
              if (vertexset->nDataSets > 0) {
                dataset = (capsDataSet *) vertexset->dataSets[0]->blind;
                if (dataset->npts != vertexset->discr->nPoints) {
                  snprintf(temp, PATH_MAX, "DataSet = %s, npts = %d %d!",
                           vertexset->dataSets[0]->name, dataset->npts,
                           vertexset->discr->nPoints);
                  caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
                  if (*errors != NULL) *nErr = (*errors)->nError;
                  aim_FreeDiscr(problem->aimFPTR, analysis->loadName,
                                vertexset->discr);
                  caps_close(object, 0, NULL);
                  return CAPS_MISMATCH;
                }
              }
              if (vertexset->nDataSets > 1) {
                dataset = (capsDataSet *) vertexset->dataSets[1]->blind;
                if (dataset->npts != vertexset->discr->nVerts) {
                  snprintf(temp, PATH_MAX, "DataSet = %s, npts = %d %d!",
                           vertexset->dataSets[1]->name, dataset->npts,
                           vertexset->discr->nVerts);
                  caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
                  if (*errors != NULL) *nErr = (*errors)->nError;
                  aim_FreeDiscr(problem->aimFPTR, analysis->loadName,
                                vertexset->discr);
                  caps_close(object, 0, NULL);
                  return CAPS_MISMATCH;
                }
              }
            }
        }
      }

  } else if ((flag == 1) || (strcasecmp(&filename[idot],".csm") == 0)) {
    object->subtype     = PARAMETRIC;
    object->name        = EG_strdup(prName);
    object->last.nLines = 0;
    object->last.lines  = NULL;
    object->last.pname  = EG_strdup(prName);
    object->last.sNum   = problem->sNum;
    caps_getStaticStrings(&problem->signature, &object->last.pID,
                          &object->last.user);

    if (problem->outLevel != 1) ocsmSetOutLevel(problem->outLevel);

    if (flag == 0) {
      /* do an OpenCSM load */
      status = ocsmLoad((char *) filename, &problem->modl);
      if (status < SUCCESS) {
        snprintf(temp, PATH_MAX, "Cannot Load %s (caps_open)!", filename);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return status;
      }
      MODL = (modl_T *) problem->modl;
    } else {
      problem->modl = MODL;
    }
    if (MODL == NULL) {
      caps_makeSimpleErr(NULL, CERROR, "Cannot get OpenCSM MODL (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_NOTFOUND;
    }
    MODL->context   = problem->context;
/*@-kepttrans@*/
    MODL->userdata  = problem;
/*@+kepttrans@*/
    MODL->tessAtEnd = 0;
    status = ocsmRegSizeCB(problem->modl, caps_sizeCB);
    if (status != SUCCESS)
      printf(" CAPS Warning: ocsmRegSizeCB = %d (caps_open)!\n", status);

    env = getenv("DUMPEGADS");
    if (env != NULL) {
      MODL->dumpEgads = 1;
      MODL->loadEgads = 1;
    }

    /* check that Branches are properly ordered */
    status = ocsmCheck(problem->modl);
    if (status < SUCCESS) {
      snprintf(temp, PATH_MAX, "ocsmCheck = %d (caps_open)!", status);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return status;
    }
    fflush(stdout);

    status = ocsmInfo(problem->modl, &nbrch, &npmtr, &nbody);
    if (status != SUCCESS) {
      snprintf(temp, PATH_MAX, "ocsmInfo returns %d (caps_open)!", status);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return status;
    }

    /* count the GeomIns and GeomOuts */
    for (ngIn = ngOut = i = 0; i < npmtr; i++) {
      status = ocsmGetPmtr(problem->modl, i+1, &type, &nrow, &ncol, name);
      if (status != SUCCESS) {
        caps_close(object, 0, NULL);
        return status;
      }
      if (type == OCSM_OUTPMTR) ngOut++;
      if (type == OCSM_DESPMTR) ngIn++;
      if (type == OCSM_CFGPMTR) ngIn++;
      if (type == OCSM_CONPMTR) ngIn++;
    }

    /* allocate the objects for the geometry inputs */
    if (ngIn != 0) {
      problem->geomIn = (capsObject **) EG_alloc(ngIn*sizeof(capsObject *));
      if (problem->geomIn == NULL) {
        caps_close(object, 0, NULL);
        return EGADS_MALLOC;
      }
      for (i = 0; i < ngIn; i++) problem->geomIn[i] = NULL;
      value = (capsValue *) EG_alloc(ngIn*sizeof(capsValue));
      if (value == NULL) {
        caps_close(object, 0, NULL);
        return EGADS_MALLOC;
      }
      for (i = j = 0; j < npmtr; j++) {
        ocsmGetPmtr(problem->modl, j+1, &type, &nrow, &ncol, name);
        if ((type != OCSM_DESPMTR) && (type != OCSM_CFGPMTR) &&
            (type != OCSM_CONPMTR)) continue;
        if ((nrow == 0) || (ncol == 0)) continue;    /* ignore string pmtrs */
        value[i].nrow            = nrow;
        value[i].ncol            = ncol;
        value[i].type            = Double;
        value[i].dim             = Scalar;
        value[i].index           = i+1;
        value[i].pIndex          = j+1;
        value[i].lfixed          = value[i].sfixed = Fixed;
        value[i].nullVal         = NotAllowed;
        value[i].units           = NULL;
        value[i].meshWriter      = NULL;
        value[i].link            = NULL;
        value[i].vals.reals      = NULL;
        value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
        value[i].linkMethod      = Copy;
        value[i].length          = value[i].nrow*value[i].ncol;
        if ((ncol > 1) && (nrow > 1)) {
          value[i].dim           = Array2D;
        } else if ((ncol > 1) || (nrow > 1)) {
          value[i].dim           = Vector;
        }
        value[i].gInType         = 0;
        if (type == OCSM_CFGPMTR) value[i].gInType = 1;
        if (type == OCSM_CONPMTR) value[i].gInType = 2;
        value[i].partial         = NULL;
        value[i].nderiv            = 0;
        value[i].derivs            = NULL;

        status = caps_makeObject(&objs);
        if (status != CAPS_SUCCESS) {
          EG_free(value);
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
        if (i == 0) objs->blind = value;
/*@-kepttrans@*/
        objs->parent       = object;
/*@+kepttrans@*/
        objs->name         = NULL;
        objs->type         = VALUE;
        objs->subtype      = GEOMETRYIN;
        objs->last.sNum    = 1;
/*@-immediatetrans@*/
        objs->blind        = &value[i];
/*@+immediatetrans@*/
        problem->geomIn[i] = objs;
        i++;
      }
      problem->nGeomIn = ngIn;
      for (i = 0; i < ngIn; i++) {
        ocsmGetPmtr(problem->modl, value[i].pIndex, &type, &nrow, &ncol, name);
        if (nrow*ncol > 1) {
          reals = (double *) EG_alloc(nrow*ncol*sizeof(double));
          if (reals == NULL) {
            caps_close(object, 0, NULL);
            return EGADS_MALLOC;
          }
          value[i].vals.reals = reals;
        } else {
          reals = &value[i].vals.real;
        }
        problem->geomIn[i]->name = EG_strdup(name);
/*      flip storage
        for (n = j = 0; j < ncol; j++)
          for (k = 0; k < nrow; k++, n++) {  */
        for (n = k = 0; k < nrow; k++)
          for (j = 0; j < ncol; j++, n++) {
            status = ocsmGetValu(problem->modl, value[i].pIndex, k+1, j+1,
                                 &reals[n], &dot);
            if (status != SUCCESS) {
              caps_close(object, 0, NULL);
              return status;
            }
          }
        if (type == OCSM_CFGPMTR) continue;
        status = ocsmGetBnds(problem->modl, value[i].pIndex, 1, 1,
                             &lower, &upper);
        if (status != SUCCESS) continue;
        if ((lower != -HUGEQ) || (upper != HUGEQ)) {
          value[i].limits.dlims[0] = lower;
          value[i].limits.dlims[1] = upper;
        }
      }
    }

    /* allocate the objects for the geometry outputs */
    if (ngOut != 0) {
      units = NULL;
      if (problem->lunits != NULL) units = problem->lunits[problem->nBodies-1];
      problem->geomOut = (capsObject **) EG_alloc(ngOut*sizeof(capsObject *));
      if (problem->geomOut == NULL) {
        caps_close(object, 0, NULL);
        return EGADS_MALLOC;
      }
      for (i = 0; i < ngOut; i++) problem->geomOut[i] = NULL;
      value = (capsValue *) EG_alloc(ngOut*sizeof(capsValue));
      if (value == NULL) {
        caps_close(object, 0, NULL);
        return EGADS_MALLOC;
      }
      for (i = j = 0; j < npmtr; j++) {
        ocsmGetPmtr(problem->modl, j+1, &type, &nrow, &ncol, name);
        if (type != OCSM_OUTPMTR) continue;
        value[i].length          = 1;
        value[i].type            = DoubleDeriv;
        value[i].nrow            = 1;
        value[i].ncol            = 1;
        value[i].dim             = Scalar;
        value[i].index           = i+1;
        value[i].pIndex          = j+1;
        value[i].lfixed          = value[i].sfixed = Change;
        value[i].nullVal         = IsNull;
        value[i].units           = NULL;
        value[i].meshWriter      = NULL;
        value[i].link            = NULL;
        value[i].vals.reals      = NULL;
        value[i].limits.dlims[0] = value[i].limits.dlims[1] = 0.0;
        value[i].linkMethod      = Copy;
        value[i].gInType         = 0;
        value[i].partial         = NULL;
        value[i].nderiv            = 0;
        value[i].derivs            = NULL;
        caps_geomOutUnits(name, units, &value[i].units);

        status = caps_makeObject(&objs);
        if (status != CAPS_SUCCESS) {
          for (k = 0; k < i; k++)
            if (value[k].length > 1) EG_free(value[k].vals.reals);
          EG_free(value);
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
/*@-kepttrans@*/
        objs->parent                   = object;
/*@+kepttrans@*/
        objs->name                     = EG_strdup(name);
        objs->type                     = VALUE;
        objs->subtype                  = GEOMETRYOUT;
        objs->last.sNum                = 0;
/*@-immediatetrans@*/
        objs->blind                    = &value[i];
/*@+immediatetrans@*/
        problem->geomOut[i]            = objs;
        problem->geomOut[i]->last.sNum = problem->sNum;
        i++;
      }
      problem->nGeomOut = ngOut;
    }

    /* write an OpenCSM checkpoint file */
#ifdef WIN32
    snprintf(current, PATH_MAX, "%s\\capsRestart.cpc", root);
#else
    snprintf(current, PATH_MAX, "%s/capsRestart.cpc",  root);
#endif
    caps_rmFile(current);
    status = ocsmSave(problem->modl, current);
    if (status != CAPS_SUCCESS) {
      caps_close(object, 0, NULL);
      return status;
    }
    fflush(stdout);

  } else if ((flag == 2) || (strcasecmp(&filename[idot],".egads") == 0)) {
    object->subtype     = STATIC;
    object->name        = EG_strdup(prName);
    object->last.nLines = 0;
    object->last.lines  = NULL;
    object->last.pname  = EG_strdup(prName);
    object->last.sNum   = problem->sNum;
    caps_getStaticStrings(&problem->signature, &object->last.pID,
                          &object->last.user);
    if (flag == 0) {
      status = EG_loadModel(problem->context, 1, filename, &model);
      if (status != EGADS_SUCCESS) {
        caps_close(object, 0, NULL);
        return status;
      }
    }
    problem->modl = model;
    status = EG_getTopology(model, &ref, &oclass, &mtype, data,
                            &problem->nBodies, &problem->bodies, &senses);
    if (status != EGADS_SUCCESS) {
      caps_close(object, 0, NULL);
      return status;
    }
    /* get length units if exist */
    if (problem->nBodies > 0) {
      problem->lunits = (char **) EG_alloc(problem->nBodies*sizeof(char *));
      if ((problem->lunits != NULL) && (problem->bodies != NULL))
        for (i = 0; i < problem->nBodies; i++)
          caps_fillLengthUnits(problem, problem->bodies[i],
                               &problem->lunits[i]);
    }

    /* get parameter values if saved by OpenCSM */
    status = EG_attributeNum(model, &nattr);
    if ((status == EGADS_SUCCESS) && (nattr != 0)) {
      for (ngIn = ngOut = i = 0; i < nattr; i++) {
        status = EG_attributeGet(model, i+1, &aname, &type, &len, &aints,
                                 &areals, &astring);
        if (status != EGADS_SUCCESS) continue;
        if (type   != ATTRREAL)      continue;
        if (strncmp(aname, "_outpmtr_", 9) == 0) ngOut++;
        if (strncmp(aname, "_despmtr_", 9) == 0) ngIn++;
        if (strncmp(aname, "_cfgpmtr_", 9) == 0) ngIn++;
      }

      /* allocate the objects for the geometry inputs */
      if (ngIn != 0) {
        problem->geomIn = (capsObject **) EG_alloc(ngIn*sizeof(capsObject *));
        if (problem->geomIn == NULL) {
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
        for (i = 0; i < ngIn; i++) problem->geomIn[i] = NULL;
        value = (capsValue *) EG_alloc(ngIn*sizeof(capsValue));
        if (value == NULL) {
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
        for (i = j = 0; j < nattr; j++) {
          status = EG_attributeGet(model, j+1, &aname, &type, &len, &aints,
                                   &areals, &astring);
          if (status != EGADS_SUCCESS) continue;
          if (type   != ATTRREAL)      continue;
          if ((strncmp(aname, "_despmtr_", 9) != 0) &&
              (strncmp(aname, "_cfgpmtr_", 9) != 0)) continue;
          value[i].nrow             = len;
          value[i].ncol             = 1;
          value[i].type             = Double;
          value[i].dim              = Scalar;
          value[i].index            = value[i].pIndex = j+1;
          value[i].lfixed           = value[i].sfixed = Fixed;
          value[i].nullVal          = NotAllowed;
          value[i].units            = NULL;
          value[i].meshWriter       = NULL;
          value[i].link             = NULL;
          value[i].vals.reals       = NULL;
          value[i].limits.dlims[0]  = value[i].limits.dlims[1] = 0.0;
          value[i].linkMethod       = Copy;
          value[i].gInType          = (strncmp(aname, "_cfgpmtr_", 9) == 0) ? 1 : 0;
          value[i].partial          = NULL;
          value[i].nderiv             = 0;
          value[i].derivs             = NULL;
          value[i].length           = len;
          if (len > 1) value[i].dim = Vector;

          status = caps_makeObject(&objs);
          if (status != CAPS_SUCCESS) {
            EG_free(value);
            caps_close(object, 0, NULL);
            return EGADS_MALLOC;
          }
          if (i == 0) objs->blind = value;
          /*@-kepttrans@*/
          objs->parent       = object;
          /*@+kepttrans@*/
          objs->name         = NULL;
          objs->type         = VALUE;
          objs->subtype      = GEOMETRYIN;
          objs->last.sNum    = 1;
          /*@-immediatetrans@*/
          objs->blind        = &value[i];
          /*@+immediatetrans@*/
          problem->geomIn[i] = objs;
          i++;
        }
        problem->nGeomIn = ngIn;
        for (i = 0; i < ngIn; i++) {
          EG_attributeGet(model, value[i].pIndex, &aname, &type, &len, &aints,
                          &areals, &astring);
          if (len > 1) {
            reals = (double *) EG_alloc(len*sizeof(double));
            if (reals == NULL) {
              caps_close(object, 0, NULL);
              return EGADS_MALLOC;
            }
            value[i].vals.reals = reals;
          } else {
            reals = &value[i].vals.real;
          }
          problem->geomIn[i]->name = EG_strdup(&aname[9]);
          for (j = 0; j < len; j++) reals[j] = areals[j];
        }
      }

      /* allocate the objects for the geometry outputs */
      if (ngOut != 0) {
        problem->geomOut = (capsObject **) EG_alloc(ngOut*sizeof(capsObject *));
        if (problem->geomOut == NULL) {
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
        for (i = 0; i < ngOut; i++) problem->geomOut[i] = NULL;
        value = (capsValue *) EG_alloc(ngOut*sizeof(capsValue));
        if (value == NULL) {
          caps_close(object, 0, NULL);
          return EGADS_MALLOC;
        }
        for (i = j = 0; j < nattr; j++) {
          status = EG_attributeGet(model, j+1, &aname, &type, &len, &aints,
                                   &areals, &astring);
          if (status != EGADS_SUCCESS) continue;
          if (type   != ATTRREAL)      continue;
          if (strncmp(aname, "_outpmtr_", 9) != 0) continue;
          value[i].nrow             = len;
          value[i].ncol             = 1;
          value[i].type             = Double;
          value[i].dim              = Scalar;
          value[i].index            = value[i].pIndex = j+1;
          value[i].lfixed           = value[i].sfixed = Fixed;
          value[i].nullVal          = NotAllowed;
          value[i].units            = NULL;
          value[i].meshWriter       = NULL;
          value[i].link             = NULL;
          value[i].vals.reals       = NULL;
          value[i].limits.dlims[0]  = value[i].limits.dlims[1] = 0.0;
          value[i].linkMethod       = Copy;
          value[i].gInType          = 0;
          value[i].partial          = NULL;
          value[i].nderiv             = 0;
          value[i].derivs             = NULL;
          value[i].length           = len;
          if (len > 1) value[i].dim = Vector;

          status = caps_makeObject(&objs);
          if (status != CAPS_SUCCESS) {
            EG_free(value);
            caps_close(object, 0, NULL);
            return EGADS_MALLOC;
          }
          if (i == 0) objs->blind = value;
          /*@-kepttrans@*/
          objs->parent        = object;
          /*@+kepttrans@*/
          objs->name          = NULL;
          objs->type          = VALUE;
          objs->subtype       = GEOMETRYOUT;
          objs->last.sNum     = 1;
          /*@-immediatetrans@*/
          objs->blind         = &value[i];
          /*@+immediatetrans@*/
          problem->geomOut[i] = objs;
          i++;
        }
        problem->nGeomOut = ngOut;
        for (i = 0; i < ngOut; i++) {
          EG_attributeGet(model, value[i].pIndex, &aname, &type, &len, &aints,
                          &areals, &astring);
          if (len > 1) {
            reals = (double *) EG_alloc(len*sizeof(double));
            if (reals == NULL) {
              caps_close(object, 0, NULL);
              return EGADS_MALLOC;
            }
            value[i].vals.reals = reals;
          } else {
            reals = &value[i].vals.real;
          }
          problem->geomOut[i]->name = EG_strdup(&aname[9]);
          for (j = 0; j < len; j++) reals[j] = areals[j];
        }
      }
      
#ifdef WIN32
      snprintf(current, PATH_MAX, "%s\\capsRestart.egads", root);
#else
      snprintf(current, PATH_MAX, "%s/capsRestart.egads",  root);
#endif
      caps_rmFile(current);
      status = EG_saveModel(model, current);
      if (status != EGADS_SUCCESS)
        printf(" CAPS Warning: Cannot save EGADS file = %d\n", status);
      
      /* set geometry serial number */
      problem->geometry.sNum = 1;
    }

  } else {
    if (flag == 0) {
      snprintf(temp, PATH_MAX,
               "Start Flag = %d  filename = %s  NOT initialized (caps_open)!",
               flag, filename);
    } else {
      snprintf(temp, PATH_MAX, "Start Flag = %d NOT initialized (caps_open)!",
               flag);
    }
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_BADINIT;
  }
  
  /* phase start -- use capsRestart & reset journal */
  if (flag == 3) {
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal.txt",
             problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal.txt",
             problem->root);
#endif
    fp = fopen(filename, "w");
    if (fp == NULL) {
      snprintf(temp, PATH_MAX, "Cannot open %s on Phase (caps_open)",
               filename);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_DIRERR;
    }
    fprintf(fp, "%d %d\n", CAPSMAJOR, CAPSMINOR);
    env = getenv("ESP_ARCH");
    if (env == NULL) {
      snprintf(temp, PATH_MAX, "ESP_ARCH env variable is not set! (caps_open)");
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_JOURNALERR;
    }
    fprintf(fp, "%s\n", env);
    env = getenv("CASREV");
    if (env == NULL) {
      snprintf(temp, PATH_MAX, "CASREV env variable is not set! (caps_open)");
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_JOURNALERR;
    }
    fprintf(fp, "%s\n", env);
    fclose(fp);
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal", problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal",   problem->root);
#endif
    status = caps_rmFile(filename);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Warning: Cannot delete %s (caps_open)!\n", filename);
/*@-dependenttrans@*/
    problem->jrnl = fopen(filename, "wb");
/*@+dependenttrans@*/
    if (problem->jrnl == NULL) {
      snprintf(temp, PATH_MAX, "Cannot open %s on Phase (caps_open)", filename);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_DIRERR;
    }
    i  = CAPS_OPEN;
    nw = fwrite(&i,   sizeof(int),       1, problem->jrnl);
    if (nw != 1) goto owrterr;
    ret = problem->sNum;
    nw  = fwrite(&ret, sizeof(CAPSLONG), 1, problem->jrnl);
    if (nw != 1) goto owrterr;
    i   = CAPS_SUCCESS;
    nw  = fwrite(&i,   sizeof(int),      1, problem->jrnl);
    if (nw != 1) goto owrterr;
  
    ivec[0] = CAPSMAJOR;
    ivec[1] = CAPSMINOR;
    nw = fwrite(ivec, sizeof(int),       2, problem->jrnl);
    if (nw != 2) goto owrterr;
  
    ret = problem->sNum;
    nw  = fwrite(&ret, sizeof(CAPSLONG), 1, problem->jrnl);
    if (nw != 1) goto owrterr;
    i   = CAPS_OPEN;
    nw  = fwrite(&i,   sizeof(int),      1, problem->jrnl);
    if (nw != 1) goto owrterr;
    fflush(problem->jrnl);
    
    /* cleanup Model files */
    for (i = 0; i < problem->nEGADSmdl; i++) {
#ifdef WIN32
      snprintf(filename, PATH_MAX, "%s\\capsRestart\\model%4.4d.egads",
               problem->root, i);
#else
      snprintf(filename, PATH_MAX, "%s/capsRestart/model%4.4d.egads",
               problem->root, i);
#endif
      status = caps_rmFile(filename);
      if (status != EGADS_SUCCESS)
        printf(" CAPS Warning: Cannot remove file: %s\n", filename);
    }
    problem->nEGADSmdl = 0;
    
    /* call post on all analysis objects & get AIMpointers */
    ret = 0;
    do {
      n    = -1;
      sNum = problem->sNum+1;
      for (i = 0; i < problem->nAnalysis; i++) {
/*@-nullderef@*/
        aobject  = problem->analysis[i];
/*@+nullderef@*/
        if (aobject  == NULL) continue;
        analysis = (capsAnalysis *) aobject->blind;
        if (analysis == NULL) continue;
        if (analysis->pre.sNum == 0)   continue;
        if (aobject->last.sNum > sNum) continue;
        if (aobject->last.sNum <= ret) continue;
        n    = i;
        sNum = aobject->last.sNum;
      }
      if (n == -1) continue;
/*@-nullderef@*/
      aobject  = problem->analysis[n];
/*@+nullderef@*/
      analysis = (capsAnalysis *) aobject->blind;
      if (analysis->nAnalysisIn > 0)
        valIn  = (capsValue *) analysis->analysisIn[0]->blind;
      status   = aim_PostAnalysis(problem->aimFPTR, analysis->loadName,
                                  analysis->instStore, &analysis->info, 1,
                                  valIn);
      if (status != CAPS_SUCCESS) {
        snprintf(temp, PATH_MAX, "Post Analysis on %s (%d) returns %d (caps_open)!",
                 analysis->loadName, analysis->info.instance, status);
        caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return CAPS_BADINIT;
      }
      for (i = 0; i < analysis->nAnalysisOut; i++) {
        value = analysis->analysisOut[i]->blind;
        if  (value == NULL) continue;
        if ((value->type != Pointer) && (value->type != PointerMesh)) continue;
        status = aim_CalcOutput(problem->aimFPTR, analysis->loadName,
                              analysis->instStore, &analysis->info, i+1, value);
        if (status != CAPS_SUCCESS)
          printf(" CAPS Warning: aim_CalcOutput = %d (caps_open)\n", status);
      }
      ret = sNum;
    } while (n != -1);

    problem->jpos = 0;
    status = caps_writeProblem(object);
    if (status != CAPS_SUCCESS) {
      caps_close(object, 0, NULL);
      return status;
    }
    
    *pobject = object;
    return CAPS_SUCCESS;
  }
  
  /* continuation -- use capsRestart & read from journal */
  if (flag == 4) {
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal.txt",
             problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal.txt",
             problem->root);
#endif
    fp = fopen(filename, "r");
    if (fp == NULL) {
      snprintf(temp, PATH_MAX, "Cannot open %s on Continuation (caps_open)",
               filename);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_BADINIT;
    } else {
      fscanf(fp, "%d %d", &ivec[0], &ivec[1]);
      env   = getenv("ESP_ARCH");
      fscanf(fp, "%s", temp);
      units = getenv("CASREV");
      fscanf(fp, "%s", filename);
      fclose(fp);
      i = 0;
      if ((env   == NULL) && (strlen(temp)     != 0)) i = 1;
      if ((units == NULL) && (strlen(filename) != 0)) i = 1;
      if (i == 0) {
        if (strcmp(env,   temp)     != 0) i = 1;
        if (strcmp(units, filename) != 0) i = 1;
      }
      if ((ivec[0] != CAPSMAJOR) || (ivec[1] != CAPSMINOR) || (i == 1)) {
        snprintf(temp, PATH_MAX, "Journal from CAPS %d.%d and running %d.%d!",
                 ivec[0], ivec[1], CAPSMAJOR, CAPSMINOR);
        snprintf(current, PATH_MAX, "Architecture %s vs %s\n", temp,     env);
        snprintf(root,    PATH_MAX, "OpenCASCADE  %s vs %s\n", filename, units);
        caps_makeSimpleErr(NULL, CERROR, temp, current, root, errors);
        if (*errors != NULL) *nErr = (*errors)->nError;
        caps_close(object, 0, NULL);
        return CAPS_BADINIT;
      }
    }
    
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal", problem->root);
#else
    snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal",   problem->root);
#endif
/*@-dependenttrans@*/
    problem->jrnl = fopen(filename, "rb");
/*@+dependenttrans@*/
    if (problem->jrnl == NULL) {
      printf(" CAPS Error: Cannot open %s for read (caps_open)!\n", filename);
      caps_close(object, 0, NULL);
      return CAPS_DIRERR;
    }
    nw = fread(&i,   sizeof(int),       1, problem->jrnl);
    if (nw != 1) goto owrterr;
    if (i != CAPS_OPEN) {
      caps_makeSimpleErr(NULL, CERROR, "Journal Sequence Fail 0 (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_IOERR;
    }
    nw = fread(&ret, sizeof(CAPSLONG),  1, problem->jrnl);
    if (nw != 1) goto owrterr;
    nw = fread(&i,   sizeof(int),       1, problem->jrnl);
    if (nw != 1) goto owrterr;
    if (i != CAPS_SUCCESS) {
      caps_makeSimpleErr(NULL, CERROR, "Journal Sequence Fail 1 (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_IOERR;
    }
    nw = fread(ivec, sizeof(int),      2, problem->jrnl);
    if (nw != 2) goto owrterr;
    if ((ivec[0] != CAPSMAJOR) || (ivec[1] != CAPSMINOR)) {
      snprintf(temp, PATH_MAX, "Journal Sequence Fail  %d %d (caps_open)!",
               ivec[0], ivec[1]);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_IOERR;
    }
    nw = fread(&ret, sizeof(CAPSLONG), 1, problem->jrnl);
    if (nw != 1) goto owrterr;
    nw = fread(&i,   sizeof(int),      1, problem->jrnl);
    if (nw != 1) goto owrterr;
    if (i != CAPS_OPEN) {
      caps_makeSimpleErr(NULL, CERROR, "Journal Sequence Fail 2 (caps_open)!",
                         NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_IOERR;
    }
    
    *pobject = object;
    return CAPS_SUCCESS;
  }
  
  /* start up */
  problem->writer.sNum = problem->sNum;
  caps_fillDateTime(problem->writer.datetime);
  
  /* setup preliminary restart data */
#ifdef WIN32
  snprintf(current, PATH_MAX, "%s\\capsRestart", root);
#else
  snprintf(current, PATH_MAX, "%s/capsRestart",  root);
#endif
  status = caps_statFile(current);
  if (status == EGADS_SUCCESS) {
    snprintf(temp, PATH_MAX, "%s is a flat file (caps_open)", current);
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_DIRERR;
  } else if (status == EGADS_NOTFOUND) {
    status = caps_mkDir(current);
    if (status != EGADS_SUCCESS) {
      snprintf(temp, PATH_MAX, "Cannot mkDir %s (caps_open)", current);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return status;
    }
    /* populate the geometry info */
#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\geom.txt", current);
    snprintf(temp,     PATH_MAX, "%s\\xxTempxx", current);
#else
    snprintf(filename, PATH_MAX, "%s/geom.txt",  current);
    snprintf(temp,     PATH_MAX, "%s/xxTempxx",  current);
#endif
    fp = fopen(temp, "w");
    if (fp == NULL) {
      snprintf(temp, PATH_MAX, "Cannot open %s (caps_open)\n", filename);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return CAPS_DIRERR;
    }
    fprintf(fp, "%d %d\n", problem->nGeomIn, problem->nGeomOut);
    if (problem->geomIn != NULL)
      for (i = 0; i < problem->nGeomIn; i++)
        fprintf(fp, "%s\n", problem->geomIn[i]->name);
    if (problem->geomOut != NULL)
      for (i = 0; i < problem->nGeomOut; i++)
        fprintf(fp, "%s\n", problem->geomOut[i]->name);
    fclose(fp);
    status = caps_rename(temp, filename);
    if (status != CAPS_SUCCESS) {
      snprintf(temp, PATH_MAX, "Cannot rename %s (caps_open)!\n", filename);
      caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
      if (*errors != NULL) *nErr = (*errors)->nError;
      caps_close(object, 0, NULL);
      return status;
    }
    status = caps_dumpGeomVals(problem, 0);
    if (status != CAPS_SUCCESS) {
      caps_close(object, 0, NULL);
      return CAPS_DIRERR;
    }
  }
  status = caps_writeProblem(object);
  if (status != CAPS_SUCCESS) {
    caps_close(object, 0, NULL);
    return status;
  }
  
  /* open journal file */
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsRestart\\capsJournal.txt",
           problem->root);
#else
  snprintf(filename, PATH_MAX, "%s/capsRestart/capsJournal.txt",
           problem->root);
#endif
  fp = fopen(filename, "w");
  if (fp == NULL) {
    snprintf(temp, PATH_MAX, "Cannot open %s on Phase (caps_open)", filename);
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_DIRERR;
  }
  fprintf(fp, "%d %d\n", CAPSMAJOR, CAPSMINOR);
  env = getenv("ESP_ARCH");
  if (env == NULL) {
    snprintf(temp, PATH_MAX, "ESP_ARCH env variable is not set! (caps_open)");
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_JOURNALERR;
  }
  fprintf(fp, "%s\n", env);
  env = getenv("CASREV");
  if (env == NULL) {
    snprintf(temp, PATH_MAX, "CASREV env variable is not set! (caps_open)");
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_JOURNALERR;
  }
  fprintf(fp, "%s\n", env);
  fclose(fp);
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\capsJournal", current);
#else
  snprintf(filename, PATH_MAX, "%s/capsJournal",  current);
#endif
/*@-dependenttrans@*/
  problem->jrnl = fopen(filename, "wb");
/*@+dependenttrans@*/
  if (problem->jrnl == NULL) {
    snprintf(temp, PATH_MAX, "Cannot open %s (caps_open)", filename);
    caps_makeSimpleErr(NULL, CERROR, temp, NULL, NULL, errors);
    if (*errors != NULL) *nErr = (*errors)->nError;
    caps_close(object, 0, NULL);
    return CAPS_DIRERR;
  }
  i   = CAPS_OPEN;
  nw  = fwrite(&i,   sizeof(int),      1, problem->jrnl);
  if (nw != 1) goto owrterr;
  ret = 0;
  nw  = fwrite(&ret, sizeof(CAPSLONG), 1, problem->jrnl);
  if (nw != 1) goto owrterr;
  i   = CAPS_SUCCESS;
  nw  = fwrite(&i,   sizeof(int),      1, problem->jrnl);
  if (nw != 1) goto owrterr;
  
  ivec[0] = CAPSMAJOR;
  ivec[1] = CAPSMINOR;
  nw  = fwrite(ivec, sizeof(int),      2, problem->jrnl);
  if (nw != 2) goto owrterr;
  
  ret = problem->sNum;
  nw  = fwrite(&ret, sizeof(CAPSLONG), 1, problem->jrnl);
  if (nw != 1) goto owrterr;
  i   = CAPS_OPEN;
  nw  = fwrite(&i,   sizeof(int),      1, problem->jrnl);
  if (nw != 1) goto owrterr;
  fflush(problem->jrnl);
  
  *pobject = object;
  return CAPS_SUCCESS;
  
owrterr:
  fclose(problem->jrnl);
  printf(" CAPS Error: IO error on journal file (caps_open)!\n");
  caps_makeSimpleErr(NULL, CERROR, "IO error on journal file (caps_open)!",
                     NULL, NULL, errors);
  if (*errors != NULL) *nErr = (*errors)->nError;
  caps_close(object, 0, NULL);
  return CAPS_IOERR;
}


int
caps_build(capsObject *pobject, int *nErr, capsErrs **errors)
{
  int          i, j, k, m, n, stat, status, gstatus, state, npts;
  int          type, nrow, ncol, buildTo, builtTo, nbody, ibody;
  char         name[MAX_NAME_LEN], error[129], vstr[MAX_STRVAL_LEN], *units;
  double       dot, *values;
  ego          body;
  modl_T       *MODL;
  capsValue    *value;
  capsProblem  *problem;
  capsAnalysis *analy;
  capsObject   *source, *object, *last;
  capsErrs     *errs;

  *nErr    = 0;
  *errors  = NULL;
  if (pobject              == NULL)      return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type        != PROBLEM)   return CAPS_BADOBJECT;
  if (pobject->blind       == NULL)      return CAPS_NULLBLIND;
  problem  = (capsProblem *)  pobject->blind;

  /* nothing to do for static geometry */
  if (pobject->subtype == STATIC)        return CAPS_CLEAN;

  /* do we need new geometry? */
  /* check for dirty geometry inputs */
  gstatus = 0;
  for (i = 0; i < problem->nGeomIn; i++) {
    source = object = problem->geomIn[i];
    do {
      if (source->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
      if (source->type        != VALUE)     return CAPS_BADTYPE;
      if (source->blind       == NULL)      return CAPS_NULLBLIND;
      value = (capsValue *) source->blind;
      if (value->link == object)            return CAPS_CIRCULARLINK;
      last   = source;
      source = value->link;
    } while (value->link != NULL);
    if (last->last.sNum > problem->geometry.sNum) {
      gstatus = 1;
      break;
    }
  }
  if ((gstatus == 0) && (problem->geometry.sNum > 0)) return CAPS_CLEAN;

  /* generate new geometry */
  MODL           = (modl_T *) problem->modl;
  MODL->context  = problem->context;
  MODL->userdata = problem;

  /* Note: Geometry In already updated */

  /* have OpenCSM do the rebuild */
  if (problem->bodies != NULL) {
    for (i = 0; i < problem->nBodies; i++)
      if (problem->lunits[i] != NULL) EG_free(problem->lunits[i]);
    /* remove old bodies & tessellations for all Analyses */
    for (i = 0; i < problem->nAnalysis; i++) {
      analy = (capsAnalysis *) problem->analysis[i]->blind;
      if (analy == NULL) continue;
      if (analy->tess != NULL) {
        for (j = 0; j < analy->nTess; j++)
          if (analy->tess[j] != NULL) {
            /* delete the body in the tessellation if not on the OpenCSM stack */
            body = NULL;
            if (j >= analy->nBody) {
              stat = EG_statusTessBody(analy->tess[j], &body, &state, &npts);
              if (stat <  EGADS_SUCCESS) return stat;
              if (stat == EGADS_OUTSIDE) return CAPS_SOURCEERR;
            }
            EG_deleteObject(analy->tess[j]);
            if (body != NULL) EG_deleteObject(body);
            analy->tess[j] = NULL;
          }
        EG_free(analy->tess);
        analy->tess   = NULL;
        analy->nTess  = 0;
      }
      if (analy->bodies != NULL) {
        EG_free(analy->bodies);
        analy->bodies = NULL;
        analy->nBody  = 0;
      }
      /* remove tracked sensitivity calculations */
      analy->info.pIndex = 0;
      analy->info.irow   = 0;
      analy->info.icol   = 0;
    }
    EG_free(problem->bodies);
    EG_free(problem->lunits);
    problem->nBodies       = 0;
    problem->bodies        = NULL;
    problem->lunits        = NULL;
    problem->geometry.sNum = 0;
  }
  buildTo = 0;                          /* all */
  nbody   = 0;
  status  = ocsmBuild(problem->modl, buildTo, &builtTo, &nbody, NULL);
  fflush(stdout);
  if (status != SUCCESS) {
    caps_makeSimpleErr(pobject, CERROR,
                       "caps_build Error: ocsmBuild fails!",
                       NULL, NULL, errors);
    if (*errors != NULL) {
      errs  = *errors;
      *nErr = errs->nError;
    }
    return status;
  }
  nbody = 0;
  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    nbody++;
  }

  units = NULL;
  if (nbody > 0) {
    problem->lunits = (char **) EG_alloc(nbody*sizeof(char *));
    problem->bodies = (ego *)   EG_alloc(nbody*sizeof(ego));
    if ((problem->bodies == NULL) || (problem->lunits == NULL)) {
      if (problem->bodies != NULL) EG_free(problem->bodies);
      if (problem->lunits != NULL) EG_free(problem->lunits);
      for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;
        if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
        EG_deleteObject(MODL->body[ibody].ebody);
      }
      caps_makeSimpleErr(pobject, CERROR,
                         "caps_build: Error on Body memory allocation!",
                         NULL, NULL, errors);
      if (*errors != NULL) {
        errs  = *errors;
        *nErr = errs->nError;
      }
      return EGADS_MALLOC;
    }
    problem->nBodies = nbody;
    i = 0;
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
      if (MODL->body[ibody].onstack != 1) continue;
      if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
      problem->bodies[i] = MODL->body[ibody].ebody;
      caps_fillLengthUnits(problem, problem->bodies[i], &problem->lunits[i]);
      i++;
    }
    units = problem->lunits[nbody-1];
  }
  caps_freeOwner(&problem->geometry);
  problem->sNum         += 1;
  problem->geometry.sNum = problem->sNum;
  caps_fillDateTime(problem->geometry.datetime);
  
  /* clear sensitivity registry */
  if (problem->regGIN != NULL) {
    for (i = 0; i < problem->nRegGIN; i++) EG_free(problem->regGIN[i].name);
    EG_free(problem->regGIN);
    problem->regGIN  = NULL;
    problem->nRegGIN = 0;
  }

  /* get geometry outputs */
  for (i = 0; i < problem->nGeomOut; i++) {
    if (problem->geomOut[i]->magicnumber != CAPSMAGIC) continue;
    if (problem->geomOut[i]->type        != VALUE)     continue;
    if (problem->geomOut[i]->blind       == NULL)      continue;
    value = (capsValue *) problem->geomOut[i]->blind;
    if (value->derivs != NULL) {
      /* clear all dots */
      for (j = 0; j < value->nderiv; j++) {
        if (value->derivs[j].name != NULL) EG_free(value->derivs[j].name);
        if (value->derivs[j].deriv  != NULL) EG_free(value->derivs[j].deriv);
      }
      EG_free(value->derivs);
      value->derivs = NULL;
      value->nderiv = 0;
    }
    if (value->type == String) {
      EG_free(value->vals.string);
      value->vals.string = NULL;
    } else {
      if (value->length != 1)
        EG_free(value->vals.reals);
      value->vals.reals = NULL;
    }
    if (value->partial != NULL) {
      EG_free(value->partial);
      value->partial = NULL;
    }
    status = ocsmGetPmtr(problem->modl, value->pIndex, &type, &nrow, &ncol,
                         name);
    if (status != SUCCESS) {
      snprintf(error, 129, "Cannot get info on Output %s",
               problem->geomOut[i]->name);
      caps_makeSimpleErr(problem->geomOut[i], CERROR,
                         "caps_build Error: ocsmGetPmtr fails!",
                         error, NULL, errors);
      if (*errors != NULL) {
        errs  = *errors;
        *nErr = errs->nError;
      }
      return status;
    }
    if (strcmp(name,problem->geomOut[i]->name) != 0) {
      snprintf(error, 129, "Cannot Geom Output[%d] %s != %s",
               i, problem->geomOut[i]->name, name);
      caps_makeSimpleErr(problem->geomOut[i], CERROR,
                         "caps_build Error: ocsmGetPmtr MisMatch!",
                         error, NULL, errors);
      if (*errors != NULL) {
        errs  = *errors;
        *nErr = errs->nError;
      }
      return CAPS_MISMATCH;
    }
    if ((nrow == 0) || (ncol == 0)) {
      status = ocsmGetValuS(problem->modl, value->pIndex, vstr);
      if (status != SUCCESS) {
        snprintf(error, 129, "Cannot get string on Output %s",
                 problem->geomOut[i]->name);
        caps_makeSimpleErr(problem->geomOut[i], CERROR,
                           "caps_build Error: ocsmGetValuSfails!",
                           error, NULL, errors);
        if (*errors != NULL) {
          errs  = *errors;
          *nErr = errs->nError;
        }
        return status;
      }
      value->nullVal     = NotNull;
      value->type        = String;
      value->length      = 1;
      value->nrow        = 1;
      value->ncol        = 1;
      value->dim         = Scalar;
      value->vals.string = EG_strdup(vstr);
      if (value->vals.string == NULL) value->nullVal = IsNull;
    } else {
      value->nullVal = NotNull;
      value->type    = DoubleDeriv;
      value->nrow    = nrow;
      value->ncol    = ncol;
      value->length  = nrow*ncol;
      value->dim     = Scalar;
      if ((nrow > 1) || (ncol > 1)) value->dim = Vector;
      if ((nrow > 1) && (ncol > 1)) value->dim = Array2D;
      if (value->length == 1) {
        values = &value->vals.real;
      } else {
        values = (double *) EG_alloc(value->length*sizeof(double));
        if (values == NULL) {
          value->nullVal = IsNull;
          snprintf(error, 129, "length = %d doubles for %s",
                   value->length, problem->geomOut[i]->name);
          caps_makeSimpleErr(problem->geomOut[i], CERROR,
                             "caps_build Error: Memory Allocation fails!",
                             error, NULL, errors);
          if (*errors != NULL) {
            errs  = *errors;
            *nErr = errs->nError;
          }
          return EGADS_MALLOC;
        }
        value->vals.reals = values;
      }
      /* flip storage
      for (n = m = j = 0; j < ncol; j++)
        for (k = 0; k < nrow; k++, n++) {  */
      for (n = m = k = 0; k < nrow; k++)
        for (j = 0; j < ncol; j++, n++) {
          status = ocsmGetValu(problem->modl, value->pIndex, k+1, j+1,
                               &values[n], &dot);
          if (status != SUCCESS) {
            snprintf(error, 129, "irow = %d icol = %d on %s",
                     k+1, j+1, problem->geomOut[i]->name);
            caps_makeSimpleErr(problem->geomOut[i], CERROR,
                               "caps_build Error: Output ocsmGetValu fails!",
                               error, NULL, errors);
            if (*errors != NULL) {
              errs  = *errors;
              *nErr = errs->nError;
            }
            return status;
          }
          if (values[n] == -HUGEQ) m++;
        }
      if (m != 0) {
        value->nullVal = IsNull;
        if (m != nrow*ncol) {
          value->partial = (int *) EG_alloc(nrow*ncol*sizeof(int));
          if (value->partial == NULL) {
            snprintf(error, 129, "nrow = %d ncol = %d on %s",
                     nrow, ncol, problem->geomOut[i]->name);
            caps_makeSimpleErr(problem->geomOut[i], CERROR,
                               "caps_build Error: Alloc of partial fails!",
                               error, NULL, errors);
            if (*errors != NULL) {
              errs  = *errors;
              *nErr = errs->nError;
            }
            return EGADS_MALLOC;
          }
          for (n = k = 0; k < nrow; k++)
            for (j = 0; j < ncol; j++, n++) {
              value->partial[n] = NotNull;
              if (values[n] == -HUGEQ) value->partial[n] = IsNull;
            }
          value->nullVal = IsPartial;
        }
      }
    }

    if (value->units != NULL) EG_free(value->units);
    value->units = NULL;
    caps_geomOutUnits(name, units, &value->units);

    caps_freeOwner(&problem->geomOut[i]->last);
    problem->geomOut[i]->last.sNum = problem->sNum;
    caps_fillDateTime(problem->geomOut[i]->last.datetime);
  }
  
  /* update the value objects for restart */
  status = caps_writeProblem(pobject);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_writeProblem = %d (caps_build)\n", status);
  status = caps_dumpGeomVals(problem, 2);
  if (status != CAPS_SUCCESS)
    printf(" CAPS Warning: caps_dumpGeomVals = %d (caps_build)\n", status);

  if (problem->nBodies == 0)
    printf(" CAPS Warning: No bodies generated (caps_build)!\n");
  
  return CAPS_SUCCESS;
}


int
caps_outLevel(capsObject *pobject, int outLevel)
{
  int          old;
  capsProblem  *problem;

  if (pobject == NULL)                   return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)          return CAPS_BADTYPE;
  if (pobject->blind == NULL)            return CAPS_NULLBLIND;
  if ((outLevel < 0) || (outLevel > 2))  return CAPS_RANGEERR;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_OUTLEVEL;

  if (pobject->subtype == PARAMETRIC) {
    ocsmSetOutLevel(outLevel);
    old = problem->outLevel;
  } else {
    old = EG_setOutLevel(problem->context, outLevel);
  }
  if (old >= 0) problem->outLevel = outLevel;

  return old;
}


int
caps_getRootPath(capsObject *pobject, const char **root)
{
  capsProblem  *problem;

  *root = NULL;
  if (pobject == NULL)                   return CAPS_NULLOBJ;
  if (pobject->magicnumber != CAPSMAGIC) return CAPS_BADOBJECT;
  if (pobject->type != PROBLEM)          return CAPS_BADTYPE;
  if (pobject->blind == NULL)            return CAPS_NULLBLIND;
  problem = (capsProblem *) pobject->blind;
  problem->funID = CAPS_GETROOTPATH;

  *root = problem->root;

  return CAPS_SUCCESS;
}
