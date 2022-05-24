/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Object Output Utility
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "caps.h"
#include <string.h>

#ifdef WIN32
#define PATH_MAX _MAX_PATH
#else
#include <limits.h>
#endif

#define DEBUG

typedef struct {
    char           *name;
    enum capsvType vtype;
    int            nrow;
    int            ncol;
    bool           nulldata;
    int            *iValue;
    double         *dValue;
    char           *sValue;
    capsTuple      *tValue;
    int            ndot;
    char           **derivNames;
    double         **derivs;
    int            *lens;
    int            *len_wrts;
} ValueData;

typedef struct {
    char      *id;
    ValueData *analysisIn;
    ValueData *analysisOut;
    ValueData *analysisDynO;
    int       inCount;
    int       outCount;
    int       dynOCount;
    bool      dirty;
} AIMData;

typedef struct {
    char *sourceAim;
    char *targetAim;
    char *sourceVar;
    char *targetVar;
} ValLinkData;

typedef struct {
    char *sourceAim;
    char *targetAim;
    char *sourceVar;
    char *targetVar;
    char *bound;
} GeomLinkData;



static int
printValues(capsObj pobject, capsObj object, int indent)
{
  int            i, j, status, len, nerr, dim, nrow, ncol, range, rank, pmtr;
  int            nLines, *ints;
  char           *str, *phase, *pname, *pID, *uID, **dotnames, **lines;
  const char     *units;
  const int      *partial;
  double         *reals;
  const void     *data;
  short          datetime[6];
  CAPSLONG       sNum;
  enum capsoType otype;
  enum capssType stype;
  enum capsvType vtype;
  enum capsFixed lfixed, sfixed;
  enum capsNull  nval;
  capsTuple      *tuple;
  capsObj        link, parent;
  capsErrs       *errors;
  capsOwn        own;

  status = caps_info(object, &pname, &otype, &stype, &link, &parent, &own);
  if (status != CAPS_SUCCESS) return status;
  status = caps_ownerInfo(pobject, own, &phase, &pname, &pID, &uID, &nLines,
                          &lines, datetime, &sNum);
  if (status != CAPS_SUCCESS) return status;
  if ((sNum == 0) && (stype != USER)) {
    for (i = 0; i < indent; i++) printf(" ");
    printf(" value = UNITITIALIZED\n");
    return CAPS_SUCCESS;
  }
  status = caps_getValue(object, &vtype, &nrow, &ncol, &data, &partial, &units,
                         &nerr, &errors);
  if (errors != NULL) caps_freeError(errors);
  if (status != CAPS_SUCCESS) return status;
  status = caps_getValueProps(object, &dim, &pmtr, &lfixed, &sfixed, &nval);
  if (status != CAPS_SUCCESS) return status;
  len   = nrow*ncol;
  range = 0;
  for (i = 0; i < indent; i++) printf(" ");
  if (len == 1) {
    printf(" value =");
  } else {
    printf(" value (%dx%d) =", nrow, ncol);
  }
  if (data == NULL) {
    printf(" NULL with len = %d", len);
    if (vtype == Tuple) printf("\n");
  } else {
    if ((vtype == Boolean) || (vtype == Integer)) {
      ints = (int *) data;
      for (i = 0; i < len; i++) printf(" %d", ints[i]);
      if (vtype == Integer) {
        data   = NULL;
        status = caps_getLimits(object, &vtype, &data, &units);
        if ((data != NULL) && (status == CAPS_SUCCESS)) range = 1;
      }
    } else if ((vtype == Double) || (vtype == DoubleDeriv)) {
      reals = (double *) data;
      for (i = 0; i < len; i++) printf(" %lf", reals[i]);
      data   = NULL;
      status = caps_getLimits(object, &vtype, &data, &units);
      if ((data != NULL) && (status == CAPS_SUCCESS)) range = 2;
    } else if (vtype == String) {
      str = (char *) data;
      printf(" %s", str);
    } else if (vtype == Tuple) {
      printf("\n");
      tuple = (capsTuple *) data;
      for (j = 0; j < len; j++) {
        for (i = 0; i < indent+2; i++) printf(" ");
        printf("%d: %s -> %s\n", j+1, tuple[j].name, tuple[j].value);
      }
    } else if ((vtype == Pointer) || (vtype == PointerMesh)) {
#if defined(WIN32) && defined(_OCC64)
      printf(" %llx", (long long) data);
#else
      printf(" %lx", (long) data);
#endif
    } else {
      return CAPS_BADTYPE;
    }
  }
  if (vtype == Tuple) return CAPS_SUCCESS;
  
  printf(" %s", units);
  if (range == 1) {
    ints = (int *) data;
    if (ints  != NULL) printf(" lims=[%d-%d]", ints[0], ints[1]);
  } else if (range == 2) {
    reals = (double *) data;
    if (reals != NULL) printf(" lims=[%lf-%lf]", reals[0], reals[1]);
  }
  printf("\n");
  if (vtype != DoubleDeriv) return CAPS_SUCCESS;

  status = caps_hasDeriv(object, &dim, &dotnames,
                         &nerr, &errors);
  if (errors != NULL) caps_freeError(errors);
  if (status != CAPS_SUCCESS) return status;
  
  for (i = 0; i < dim; i++) {
    status = caps_getDeriv(object, dotnames[i], &len, &rank, &reals,
                         &nerr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) continue;
    for (j = 0; j < indent; j++) printf(" ");
    printf(" dot %2d: %s  rank = %d\n  ", i+1, dotnames[i], rank);
    if (reals != NULL) {
      for (j = 0; j < indent+2; j++) printf(" ");
      for (j = 0; j < len*rank; j++) printf(" %lf", reals[j]);
    }
    printf("\n");
  }
  EG_free(dotnames);
  
  return CAPS_SUCCESS;
}


void
caps_printObjects(capsObj pobject, capsObj object, int indent)
{
  int            i, j, k, status, nBody, nTess, nParam, nGeomIn, nGeomOut;
  int            nAnalysis, nBound, nAnalIn, nAnalOut, nConnect, nUnConnect;
  int            nAttr, nHist, nDataSet, npts, rank, nErr, nLines, nAnalDynO;
  char           *name, *units, *phase, *pname, *pID, *userID, *oname, **lines;
  double         *data;
  short          datetime[6];
  ego            *eobjs;
  CAPSLONG       sn, lsn;
  capsObj        link, parent, obj, attr;
  capsOwn        own, *hist;
  capsErrs       *errors;
  enum capsoType type, otype;
  enum capssType subtype, osubtype;
  static char    *oType[9]  = { "BODIES", "ATTRIBUTES", "UNUSED", "PROBLEM",
                                "VALUE", "ANALYSIS", "BOUND", "VERTEXSET",
                                "DATASET" };
  static char    *sType[12] = { "NONE", "STATIC", "PARAMETRIC", "GEOMETRYIN",
                                "GEOMETRYOUT", "PARAMETER", "USER",
                                "ANALYSISIN", "ANALYSISOUT", "CONNECTED",
                                "UNCONNECTED", "ANALYSISDYNO" };

#ifdef TESTHIERACHY
  {
    extern int caps_findProblem(const capsObject *object, int funID,
                                      capsObject **pobject);
    extern int caps_hierarchy(capsObject *obj, char **full);
    extern int caps_string2obj(capsProblem *problm, /*@null@*/ const char *full,
                               capsObject **object);
    
    char        *full;
    capsProblem *problem;
    capsObject  *pobject, *tmp;
    
    status = caps_findProblem(object, 9999, &pobject);
    if (status != CAPS_SUCCESS) {
      printf(" caps_findProblem = %d\n", status);
      exit(1);
    }
    problem = pobject->blind;
    status  = caps_hierarchy(object, &full);
    if (status != CAPS_SUCCESS) {
      printf(" caps_hierarchy = %d\n", status);
      exit(1);
    }
    if (full != NULL) {
      status = caps_string2obj(problem, full, &tmp);
      EG_free(full);
      if (status != CAPS_SUCCESS) {
        printf(" caps_string2obj = %d\n", status);
        exit(1);
      }
      if (tmp != object) {
        printf(" tmp = %lx, object = %lx\n", (long) tmp, (long) object);
        exit(1);
      }
    }
  }
#endif
  
  status = caps_info(object, &name, &type, &subtype, &link, &parent, &own);
  if (status != CAPS_SUCCESS) {
#if defined(WIN32) && defined(_OCC64)
    printf(" CAPS Error: Object %llx returns %d from caps_info!\n",
           (long long) object, status);
#else
    printf(" CAPS Error: Object %lx returns %d from caps_info!\n",
           (long) object, status);
#endif
    return;
  }
  status = caps_size(object, ATTRIBUTES, NONE, &nAttr, &nErr, &errors);
  if (errors != NULL) caps_freeError(errors);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Object %s returns %d from caps_size(Attribute)!\n",
           name, status);
    return;
  }
  for (i = 0; i < indent; i++) printf(" ");
  printf(" %s has type %s, subtype %s with %d attributes\n",
         name, oType[type+2], sType[subtype], nAttr);
  /* output owner */
  status = caps_ownerInfo(pobject, own, &phase, &pname, &pID, &userID, &nLines,
                          &lines, datetime, &sn);
  if (status != CAPS_SUCCESS) {
    printf(" CAPS Error: Object %s returns %d from caps_ownerInfo!\n",
           name, status);
  } else {
    if (datetime[0] != 0) {
      for (i = 0; i < indent; i++) printf(" ");
#ifdef WIN32
      printf("   last: %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %lld\n", pname,
             pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
             datetime[4], datetime[5], sn);
#else
      printf("   last: %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %ld\n", pname,
             pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
             datetime[4], datetime[5], sn);
#endif
    }
  }
  
  /* output attributes */
  for (j = 1; j <= nAttr; j++) {
    status = caps_attrByIndex(object, j, &attr);
    if (status != CAPS_SUCCESS) {
      for (i = 0; i < indent+2; i++) printf(" ");
      printf(" CAPS Error: Object %s Attr %d ret=%d from caps_attrByIndex!\n",
             name, j, status);
      continue;
    }
    caps_printObjects(pobject, attr, indent+2);
  }
  
  /* output history */
  status = caps_getHistory(object, &nHist, &hist);
  if (status == CAPS_SUCCESS)
    for (j = 0; j < nHist; j++) {
      status = caps_ownerInfo(pobject, hist[j], &phase, &pname, &pID, &userID,
                              &nLines, &lines, datetime, &sn);
      if (status == CAPS_SUCCESS) {
        for (i = 0; i < indent; i++) printf(" ");
        printf("   hist: %s", phase);
        if (nLines > 0) printf(" -> %s", lines[0]);
        printf("\n");
        for (i = 0; i < indent; i++) printf(" ");
#ifdef WIN32
        printf("         %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %lld\n", pname,
               pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
               datetime[4], datetime[5], sn);
#else
        printf("         %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %ld\n", pname,
               pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
               datetime[4], datetime[5], sn);
#endif
        for (k = 1; k < nLines; k++) {
          for (i = 0; i < indent; i++) printf(" ");
          printf("         %s\n", lines[k]);
          for (i = 0; i < indent; i++) printf(" ");
#ifdef WIN32
          printf("         %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %lld\n", pname,
                 pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
                 datetime[4], datetime[5], sn);
#else
          printf("         %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %ld\n", pname,
                 pID, userID, datetime[0], datetime[1], datetime[2], datetime[3],
                 datetime[4], datetime[5], sn);
#endif
        }
      }
    }
  
  if (type == PROBLEM) {
    status = caps_size(object, BODIES, NONE, &nBody, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(Body)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VALUE, PARAMETER, &nParam, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(Parameter)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VALUE, GEOMETRYIN, &nGeomIn, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(GeomIn)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VALUE, GEOMETRYOUT, &nGeomOut, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(GeomOut)!\n",
             name, status);
      return;
    }
    status = caps_size(object, ANALYSIS, NONE, &nAnalysis, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(Analysis)!\n",
             name, status);
      return;
    }
    status = caps_size(object, BOUND, NONE, &nBound, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(Bound)!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d Bodies, %d Parameters, %d GeomIns, %d GeomOuts,",
           nBody, nParam, nGeomIn, nGeomOut);
    printf(" %d Analyses, %d Bounds\n", nAnalysis, nBound);

    if (nParam > 0) {
      printf("\n");
      for (i = 0; i < nParam; i++) {
        status = caps_childByIndex(object, VALUE, PARAMETER, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(Param,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
  
    if (nGeomIn > 0) {
      printf("\n");
      for (i = 0; i < nGeomIn; i++) {
        status = caps_childByIndex(object, VALUE, GEOMETRYIN, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(GeomIn,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    if (nGeomOut > 0) {
      printf("\n");
      for (i = 0; i < nGeomOut; i++) {
        status = caps_childByIndex(object, VALUE, GEOMETRYOUT, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(GeomOut,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
    if (nAnalysis > 0) {
      printf("\n");
      for (i = 0; i < nAnalysis; i++) {
        status = caps_childByIndex(object, ANALYSIS, NONE, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(Analysis,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
    if (nBound > 0) {
      printf("\n");
      for (i = 0; i < nBound; i++) {
        status = caps_childByIndex(object, BOUND, NONE, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(Bound,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }

  } else if (type == VALUE) {
    if (link != NULL) {
      for (i = 0; i < indent+2; i++) printf(" ");
      printf(" linked to %s\n", link->name);
    } else {
      status = printValues(pobject, object, indent+2);
      if (status != CAPS_SUCCESS) {
        for (i = 0; i < indent; i++) printf(" ");
        printf(" CAPS Error: printVal returns %d!\n", status);
      }
    }
    
  } else if (type == ANALYSIS) {
    status = caps_size(object, VALUE, ANALYSISIN, &nAnalIn, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(AnalysisIn)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VALUE, ANALYSISOUT, &nAnalOut, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(AnalysisIn)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VALUE, ANALYSISDYNO, &nAnalDynO, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(AnalysisIn)!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d AnalysisIns, %d AnalysisOuts, %d AnalysisDynOs\n",
           nAnalIn, nAnalOut, nAnalDynO);

    status = caps_getBodies(object, &nBody, &eobjs, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_getBodies!\n",
             name, status);
      return;
    }
    status = caps_getTessels(object, &nTess, &eobjs, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_getTessels!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d Bodies, %d Tessellations\n", nBody, nTess);

    if (nAnalIn > 0) {
      printf("\n");
      for (i = 0; i < nAnalIn; i++) {
        status = caps_childByIndex(object, VALUE, ANALYSISIN, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(AnalIn,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
    if (nAnalOut > 0) {
      printf("\n");
      for (i = 0; i < nAnalOut; i++) {
        status = caps_childByIndex(object, VALUE, ANALYSISOUT, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(AnalOut,%d)!\n",
                 name, status, i+1);
          return;
        }
        sn     = 0;
        status = caps_info(obj, &oname, &otype, &osubtype, &link, &parent, &own);
        if (status == CAPS_SUCCESS)
          caps_ownerInfo(pobject, own, &phase, &pname, &pID, &userID, &nLines,
                         &lines, datetime, &sn);
        if (sn == 0) continue;
        
        caps_printObjects(pobject, obj, indent+2);

        status = caps_info(obj, &oname, &otype, &osubtype, &link, &parent, &own);
        if (status == CAPS_SUCCESS) {
          status = caps_ownerInfo(pobject, own, &phase, &pname, &pID, &userID,
                                  &nLines, &lines, datetime, &lsn);
          if (status == CAPS_SUCCESS)
            if (lsn != sn) {
              for (j = 0; j < indent+2; j++) printf(" ");
#ifdef WIN32
              printf("   lazy: %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %lld\n",
                     pname, pID, userID, datetime[0], datetime[1], datetime[2],
                     datetime[3], datetime[4], datetime[5], lsn);
#else
              printf("   lazy: %s %s %s  %d/%02d/%02d  %02d:%02d:%02d  %ld\n",
                     pname, pID, userID, datetime[0], datetime[1], datetime[2],
                     datetime[3], datetime[4], datetime[5], lsn);
#endif
              for (k = 0; k < nLines; k++) {
                for (j = 0; j < indent; j++) printf(" ");
                printf("         %s\n", lines[k]);
              }
            }
        }
      }
    }
    
    if (nAnalDynO > 0) {
      printf("\n");
      for (i = 0; i < nAnalDynO; i++) {
        status = caps_childByIndex(object, VALUE, ANALYSISDYNO, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(AnalDynO,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
  } else if (type == BOUND) {
    status = caps_size(object, VERTEXSET, CONNECTED, &nConnect, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(VSconnected)!\n",
             name, status);
      return;
    }
    status = caps_size(object, VERTEXSET, UNCONNECTED, &nUnConnect,
                       &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(VSconnected)!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d Connecteds, %d UnConnecteds\n", nConnect, nUnConnect);

    if (nConnect > 0) {
      printf("\n");
      for (i = 0; i < nConnect; i++) {
        status = caps_childByIndex(object, VERTEXSET, CONNECTED, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(Connect,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
    if (nUnConnect > 0) {
      printf("\n");
      for (i = 0; i < nUnConnect; i++) {
        status = caps_childByIndex(object, VERTEXSET, UNCONNECTED, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(UnConnect,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
  
  } else if (type == VERTEXSET) {
    status = caps_size(object, DATASET, NONE, &nDataSet, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s returns %d from caps_size(DataSet)!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d DataSets\n", nDataSet);
    
    if (nDataSet > 0) {
      printf("\n");
      for (i = 0; i < nDataSet; i++) {
        status = caps_childByIndex(object, DATASET, NONE, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(DataSet,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(pobject, obj, indent+2);
      }
    }
    
  } else if (type == DATASET) {
    status = caps_getData(object, &npts, &rank, &data, &units, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
      printf(" CAPS Error: Object %s return=%d from caps_getData!\n",
             name, status);
      return;
    }
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d points, rank=%d, units=%s\n", npts, rank, units);
    
  }

}


static /*@null@*/ char *valueString(ValueData *valObj)
{
    char *varStr, *valueStr, *derivStr, *tmpStr;
    int i, j, k, strIndex;
    capsTuple *currentTuple;
    
    varStr = (char *) malloc((70 + strlen(valObj->name))* sizeof(char));
    if (varStr == NULL) return NULL;

    sprintf(varStr, "{ \"name\": \"%s\",\n", valObj->name);
    strcat(varStr, "\t\"value\": ");
    strIndex = 0;
    if (!valObj->nulldata) {
        if (valObj->vtype == String || valObj->vtype == Pointer ||
            valObj->vtype == PointerMesh) {
            valueStr = (char *) malloc((4 + strlen(valObj->sValue))*sizeof(char));
            if (valueStr == NULL) {
                free(varStr);
                return NULL;
            }
            sprintf(valueStr, "\"%s\"", valObj->sValue);
            free(valObj->sValue);
        } else {
            valueStr = (char *) malloc((10+valObj->nrow*valObj->ncol*10)*
                                       sizeof(char));
            if (valueStr == NULL) {
                free(varStr);
                return NULL;
            }
            if (valObj->nrow > 1) {
                strIndex += sprintf(&valueStr[strIndex], "[");
            }
            for (i = 0; i < valObj->nrow; i++) {
                if (i > 0) {
                    strIndex += sprintf(&valueStr[strIndex], ",\n");
                }
                if (valObj->ncol > 1) {
                    strIndex += sprintf(&valueStr[strIndex], "[");
                }
                for (j = 0; j < valObj->ncol; j++) {
                    if (j > 0) {
                        strIndex += sprintf(&valueStr[strIndex], ", ");
                    }
                    if ((valObj->vtype == Boolean) || (valObj->vtype == Integer)) {
                        strIndex += sprintf(&valueStr[strIndex], "%d",
                                            valObj->iValue[i * valObj->ncol + j]);
                    } else if ((valObj->vtype == Double) ||
                               (valObj->vtype == DoubleDeriv)) {
                        strIndex += sprintf(&valueStr[strIndex], "%f",
                                            valObj->dValue[i * valObj->ncol + j]);
                    } else {
#ifndef __clang_analyzer__
                        currentTuple = &valObj->tValue[i*valObj->ncol+j];
                        /* check if tuple value is string without quotes */
                        tmpStr = realloc(valueStr, (strlen(currentTuple->name)  +
                                                    strlen(currentTuple->value) +
                                                    strIndex+15));
                        if (tmpStr == NULL) {
                            free(valueStr);
                            free(varStr);
                            return NULL;
                        }
                        valueStr = tmpStr;
                        if ((currentTuple->value[0] != '{' &&
                            currentTuple->value[0] != '[') &&
                            currentTuple->value[0] != '"') {
                            strIndex += sprintf(&valueStr[strIndex], "{\"%s\": \"%s\"}",
                                                currentTuple->name,
                                                currentTuple->value);
                        } else {
                            strIndex += sprintf(&valueStr[strIndex], "{\"%s\": %s}",
                                                currentTuple->name,
                                                currentTuple->value);
                        }
#endif
                    }                    
                }
                if (valObj->ncol > 1) {
                    strIndex += sprintf(&valueStr[strIndex], "]");
                }
            }
            if (valObj->nrow > 1) {
                sprintf(&valueStr[strIndex], "]");
            }

            if ((valObj->vtype == Boolean) || (valObj->vtype == Integer)) {
                 free(valObj->iValue);
            } else if ((valObj->vtype == Double) || (valObj->vtype == DoubleDeriv)) {
                free(valObj->dValue);
            } else {
                free(valObj->tValue);
            }      
        }

        tmpStr = (char *) realloc(varStr,
                                  (strlen(valueStr)+strlen(varStr)+30)*sizeof(char));
        if (tmpStr == NULL) {
            free(valueStr);
            free(varStr);
            return NULL;
        }
        varStr = tmpStr;
        strcat(varStr, valueStr);
        strcat(varStr, ",");
        free(valueStr);
    } else {
        strcat(varStr, "null,\n");
    }

    if ((!valObj->nulldata) && (valObj->vtype == DoubleDeriv)) {
        derivStr = (char *) malloc(50 * sizeof(char));
        if (derivStr == NULL) {
            free(varStr);
            return NULL;
        }
        strIndex = sprintf(derivStr, " \"deriv\": {");
        for (i = 0; i < valObj->ndot; i++) {
            tmpStr = realloc(derivStr, 
                             (valObj->lens[i] * valObj->len_wrts[i] * 30 +
                              strlen(derivStr)+30)*sizeof(char));
            if (tmpStr == NULL) {
                free(derivStr);
                free(varStr);
                return NULL;
            }
            derivStr = tmpStr;
            if (i > 0) {
                strIndex += sprintf(&derivStr[strIndex], ",");
            }
            strIndex += sprintf(&derivStr[strIndex], "\n  \"%s\": ",
                                valObj->derivNames[i]);
            if (valObj->lens[i] > 1) {
                strIndex += sprintf(&derivStr[strIndex], "[");
            }
            
            for (j=0; j<valObj->lens[i]; j++) {
                if (j > 0) {
                    strIndex += sprintf(&derivStr[strIndex], ",\n");
                }
                if (valObj->len_wrts[i] > 1) {
                    strIndex += sprintf(&derivStr[strIndex], "[");
                }
                for (k=0; k<valObj->len_wrts[i]; k++) {
                    if (k > 0) {
                        strIndex += sprintf(&derivStr[strIndex], ", ");
                    }
                    strIndex += sprintf(&derivStr[strIndex], "%f",
                                        valObj->derivs[i][j*valObj->lens[i]+k]);
                }
                if (valObj->len_wrts[i] > 1) {
                    strIndex += sprintf(&derivStr[strIndex], "]");
                }
            }
            if (valObj->lens[i] > 1) {
                strIndex += sprintf(&derivStr[strIndex], "]");
            }
        }
        sprintf(&derivStr[strIndex], "}");
        tmpStr = (char *) realloc(varStr,
                                  (strlen(varStr)+strlen(derivStr)+10)*sizeof(char));
        if (tmpStr == NULL) {
            free(derivStr);
            free(varStr);
            return NULL;
        }
        varStr = tmpStr;
        strcat(varStr, derivStr);
        free(derivStr);
        if (valObj->derivNames != NULL) free(valObj->derivNames);
        if (valObj->derivs     != NULL) free(valObj->derivs);
        if (valObj->lens       != NULL) free(valObj->lens);
        if (valObj->len_wrts   != NULL) free(valObj->len_wrts);
    } else {
        strcat(varStr, "  \"deriv\": null");
    }

    strcat(varStr, "}");
    return varStr;
}


static int getValueData(capsObj *valueObj, ValueData *valObj)
{
    int            status, ncol, nrow, dim, pmtr, valLen, *ints;
    int            i, nErr, ndot, len, len_wrt;
    double         *reals, *deriv;
    const char     *units;
    char           *str, **names;
    const int      *partial;
    const void     *data;
    capsErrs       *errors;
    capsTuple      *tuple;
    enum capsvType vtype;
    enum capsFixed lfixed, sfixed;
    enum capsNull  nval;

    status = caps_getValue(*valueObj, &vtype, &nrow, &ncol, &data, &partial,
                           &units, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Object %s ret=%d from caps_getValue(Value)!\n",
                valObj->name, status);
        valObj->nulldata = true;
        return status;
    }
    status = caps_getValueProps(*valueObj, &dim, &pmtr, &lfixed, &sfixed, &nval);
    if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Object %s ret=%d from caps_getValueProps(Value)!\n",
               valObj->name, status);
        return status;
    }

    valObj->nrow     = nrow;
    valObj->ncol     = ncol;
    valObj->vtype    = vtype;
    valObj->nulldata = (data == NULL);
    
    valLen = nrow * ncol;
    if (data != NULL) {
        if ((vtype == Boolean) || (vtype == Integer)) {
            ints = (int *) data;
            valObj->iValue = (int *) malloc(valLen * sizeof(int));
            if (valObj->iValue == NULL) return EGADS_MALLOC;
            for (i = 0; i < valLen; i++) {
                valObj->iValue[i] = ints[i];
            }
        } else if ((vtype == Double) || (vtype == DoubleDeriv)) {
            reals = (double *) data;
            valObj->dValue = (double *) malloc(valLen * sizeof(double));
            if (valObj->dValue == NULL) return EGADS_MALLOC;
            for (i = 0; i < valLen; i++) {
                valObj->dValue[i] = reals[i];
            }
            // if DoubleDeriv, get derivative values
            if (vtype == DoubleDeriv) {
                valObj->ndot = 0;
                status = caps_hasDeriv(*valueObj, &ndot, &names, &nErr, &errors);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Object %s ret=%d from caps_hasDeriv(Value)!\n",
                            valObj->name, status);
                }
                valObj->ndot       = ndot;
                valObj->derivNames = NULL;
                valObj->derivs     = NULL;
                valObj->lens       = NULL;
                valObj->len_wrts   = NULL;
                if (ndot > 0) {
                    valObj->derivNames = (char **) malloc(ndot * sizeof(char *));
                    valObj->derivs     = (double **) malloc(ndot * sizeof(double *));
                    valObj->lens       = (int *) malloc(ndot * sizeof(int));
                    valObj->len_wrts   = (int *) malloc(ndot * sizeof(int));

                    if ((valObj->derivNames == NULL) || (valObj->derivs   == NULL) || 
                        (valObj->lens       == NULL) || (valObj->len_wrts == NULL)) {
                        if (valObj->derivNames != NULL) free(valObj->derivNames);
                        if (valObj->derivs     != NULL) free(valObj->derivs);
                        if (valObj->lens       != NULL) free(valObj->lens);
                        if (valObj->len_wrts   != NULL) free(valObj->len_wrts);
                        return EGADS_MALLOC;
                    }
                    for (i = 0; i < ndot; i++) {
                        status = caps_getDeriv(*valueObj, names[i], &len,
                                               &len_wrt, &deriv,
                                              &nErr, &errors);
                        if (status != CAPS_SUCCESS) {
                            printf(" CAPS Error: Object %s ret=%d from caps_getDeriv(Value_%d)!\n",
                                    valObj->name, status, i);
                            continue;
                        }
                        printf("  %s: %dx%d\n", names[i], len, len_wrt);
                        valObj->derivNames[i] = names[i];
                        valObj->derivs[i]     = deriv;
                        valObj->lens[i]       = len;
                        valObj->len_wrts[i]   = len_wrt;
                    }
                }
                
            }
        } else if (vtype == String) {
            str = (char *) data;
            valObj->sValue = (char *) malloc((strlen(str)+1) * sizeof(char));
            if (valObj->sValue == NULL) return EGADS_MALLOC;
            strcpy(valObj->sValue, str);
        } else if (vtype == Tuple) {
            tuple = (capsTuple *) data;
            valObj->tValue = (capsTuple *) malloc(valLen * sizeof(capsTuple));
            if (valObj->tValue == NULL) return EGADS_MALLOC;
            for (i = 0; i < valLen; i++) {
                valObj->tValue[i] = tuple[i];
            }
        } else if (vtype == Pointer || vtype == PointerMesh) {
            valObj->sValue = (char *) malloc(10 * sizeof(char));
            if (valObj->sValue == NULL) return EGADS_MALLOC;
            strcpy(valObj->sValue, "pointer");
        } else {
            /* if not of above types, set to null */
            valObj->nulldata = true;
        }
    }
  
    return CAPS_SUCCESS;
}


int caps_outputObjects(capsObj problemObj, /*@null@*/ char **stream)
{
    int              status, nErr, i, j, k, nGpts, nDpts, nConnect, nDataSet;
    int              nAnalysis, nBound, nAnalIn, nAnalOut, nDynAnalOut;
    int              major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char             *aName, *vName, *bName, *dName, *name, *env = NULL;
    char             *analysisPath = NULL, *unitSystem, *intents, **fnames;
    char             *varStr, *tmpStr, filename[PATH_MAX], *jsonText = NULL;
    capsObj          link, parent, dataLink;
    capsObj          analysisObj, valueObj, vertexObj, dataSetObj, boundObject;
    capsObj          bObj, aObj;
    capsOwn          own;
    capsErrs         *errors;
    enum capsoType   type;
    enum capssType   subtype;
    enum capsfType   ftype;
    enum capsdMethod dmethod;
    AIMData          *aims          = NULL;
    ValLinkData      *valLinks      = NULL, *tmpValLink;
    GeomLinkData     *geomLinks     = NULL, *tmpGeomLink;
    int              valLinksIndex  = 0;
    int              geomLinksIndex = 0;
    FILE             *fptr;

#ifdef DEBUG
    printf("\n In caps_outputObjects:\n");
#endif

    if (stream == NULL) {
        /* for now key off CAPS_FLOW (no spaces in path!) */
        env = getenv("CAPS_FLOW");
        if (env == NULL) {
            printf(" CAPS_Error: CAPS_FLOW not in the environment\n");
            return CAPS_NOTIMPLEMENT;
        }
        /* build up our filename from CAPS_FLOW */
        k = strlen(env);
        for (j = k-1; j >= 0; j--)
            if (env[j] == ' ') break;
        if (j < 0) {
            printf(" CAPS_Error: Bad Environment: %s\n", env);
            return CAPS_BADNAME;
        }
        snprintf(filename, PATH_MAX, "%s", &env[j+1]);
        k = strlen(filename);
        for (i = k-1; i >= 0; i--)
            if (filename[i] == '.') break;
        if (i < 0) {
            printf(" CAPS_Error: Bad environment: %s\n", filename);
            return CAPS_BADNAME;
        }
        filename[i  ] = '_';
        filename[i+1] = 'd';
        filename[i+2] = 'a';
        filename[i+3] = 't';
        filename[i+4] = 'a';
        filename[i+5] = '.';
        filename[i+6] = 'j';
        filename[i+7] = 's';
        filename[i+8] =  0;
#ifdef DEBUG
        printf("   JavaScript filename: %s\n", filename);
#endif
    } else {
        *stream = NULL;
    }
  
    /* make CAPS static */
    status = caps_debug(problemObj);
    if (status != 1) {
        printf(" CAPS_Error: CAPS not in debug mode\n");
        return CAPS_STATEERR;
    }

    /* Get number of AIMS */
    status = caps_size(problemObj, ANALYSIS, NONE, &nAnalysis, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Problem Obj returns %d from caps_size(Analysis)!\n",
               status);
        goto cleanup;
    }
#ifdef DEBUG
    printf("   %d Analyses found\n", nAnalysis);
#endif
    if (nAnalysis == 0) {
        status = CAPS_STATEERR;
        goto cleanup;
    }

    status = caps_size(problemObj, BOUND, NONE, &nBound, &nErr, &errors);
    if (errors != NULL) caps_freeError(errors);
    if (status != CAPS_SUCCESS) {
        printf(" CAPS Error: Problem Obj returns %d from caps_size(Bound)!\n",
                status);
        goto cleanup;
    }
#ifdef DEBUG
    printf("   %d Bounds found\n", nBound);
#endif

    aims = (AIMData *) malloc(nAnalysis * sizeof(AIMData));
    if (aims == NULL) {
        printf(" CAPS Error: Cannot allocate %d AIM storage!\n", nAnalysis);
        status = EGADS_MALLOC;
        goto cleanup;
    }
    for (i = 0; i < nAnalysis; i++) {
        aims[i].id           = NULL;
        aims[i].inCount      = 0;
        aims[i].outCount     = 0;
        aims[i].dynOCount    = 0;
        aims[i].analysisIn   = NULL;
        aims[i].analysisOut  = NULL;
        aims[i].analysisDynO = NULL;
    }

    /* loop through each analysis (AIM) */
    for (i = 0; i < nAnalysis; i++) {
        status = caps_childByIndex(problemObj, ANALYSIS, NONE, i+1,
                                   &analysisObj);
        if (status != CAPS_SUCCESS) {
            printf(" CAPS Error: Problem Obj ret=%d from caps_cBI(Analy,%d)!\n",
                   status, i+1);
            goto cleanup;
        }
        status = caps_info(analysisObj, &aName, &type, &subtype, &link,
                           &parent, &own);
        if (status != CAPS_SUCCESS) {
            printf(" CAPS Error: Analy Obj %d ret=%d from caps_info!\n",
                   i+1, status);
            goto cleanup;
        }
        /* analysisInfo -> dirty: color Box differently, don't populate output */
        status = caps_analysisInfo(analysisObj, &analysisPath, &unitSystem,
                                   &major, &minor, &intents, &nFields, &fnames,
                                   &ranks, &fInOut, &exec, &dirty);
        if (status != CAPS_SUCCESS) {
            printf(" CAPS Error: Analy Obj %d ret=%d from caps_analysisInfo!\n",
                   i+1, status);
            goto cleanup;
        }
        if (dirty != 0) {
            aims[i].dirty = true;
        } else {
            aims[i].dirty = false;
        }
        /* get # inputs */
        status = caps_size(analysisObj, VALUE, ANALYSISIN, &nAnalIn,
                           &nErr, &errors);
        if (errors != NULL) caps_freeError(errors);
        if (status != CAPS_SUCCESS) {
            printf(" CAPS Error: Obj %s ret=%d from caps_size(AnalysisIn)!\n",
                    aName, status);
            goto cleanup;
        }
        /* get # outputs */
        status = caps_size(analysisObj, VALUE, ANALYSISOUT, &nAnalOut,
                           &nErr, &errors);
        if (errors != NULL) caps_freeError(errors);
        if (status != CAPS_SUCCESS) {
            printf(" CAPS Error: Obj %s ret=%d from caps_size(AnalysisOut)!\n",
                   aName, status);
            goto cleanup;
        }
        /* get # dynamic outputs */
        if (dirty == 0) {
            status = caps_size(analysisObj, VALUE, ANALYSISDYNO, &nDynAnalOut,
                               &nErr, &errors);
            if (errors != NULL) caps_freeError(errors);
            if (status != CAPS_SUCCESS) {
                printf(" CAPS Error: Obj %s ret=%d from caps_size(AnalysisDynO)!\n",
                       aName, status);
                goto cleanup;
            }
        } else {
            nDynAnalOut = 0;
        }
        aims[i].id        = aName;
        aims[i].inCount   = nAnalIn;
        aims[i].outCount  = nAnalOut;
        aims[i].dynOCount = nDynAnalOut;

        if (nAnalIn > 0) {
            aims[i].analysisIn = (ValueData *) malloc(nAnalIn*sizeof(ValueData));
            if (aims[i].analysisIn == NULL) {
                printf(" CAPS Error: Name Malloc for %s -- %d AnalysisIn!\n",
                       aName, nAnalIn);
                status = EGADS_MALLOC;
                goto cleanup;
            }
            for (j = 0; j < nAnalIn; j++) {
                /* get value object */
                status = caps_childByIndex(analysisObj, VALUE, ANALYSISIN, j+1,
                                           &valueObj);
               if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Object %s ret=%d from caps_child(AnalIn,%d)!\n",
                           aName, status, j+1);
                    goto cleanup;
                }
                status = caps_info(valueObj, &vName, &type, &subtype, &dataLink,
                                   &parent, &own);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Object %s ret=%d from caps_info(AnalIn,%d)!\n",
                           vName, status, j+1);
                    goto cleanup;
                }
                aims[i].analysisIn[j].name = vName;
                status = getValueData(&valueObj, &aims[i].analysisIn[j]);
                if (status != CAPS_SUCCESS) goto cleanup;
                /* if linked to another aim */
                if (dataLink != NULL) {
                    status = caps_info(dataLink, &name, &type, &subtype, &link,
                                       &parent, &own);
                    if (status != CAPS_SUCCESS) {
                        printf(" CAPS Error: Object %s ret=%d from caps_info(link)!\n",
                               vName, status);
                        goto cleanup;
                    }
                    tmpValLink = (ValLinkData *) realloc(valLinks,
                                       (valLinksIndex+1) * sizeof(ValLinkData));
                    if (tmpValLink == NULL) {
                        printf(" CAPS Error: Link Malloc for %s -- %d Links!\n",
                               aName, valLinksIndex+1);
                        status = EGADS_MALLOC;
                        goto cleanup;
                    }
                    valLinks = tmpValLink;
                    valLinks[valLinksIndex].sourceAim = parent->name;
/*@-kepttrans@*/
                    valLinks[valLinksIndex].targetAim = aName;
                    valLinks[valLinksIndex].sourceVar = name;
                    valLinks[valLinksIndex].targetVar = vName;
/*@+kepttrans@*/
                    valLinksIndex++;
                }
            }
        }
        if (nAnalOut > 0) {
            aims[i].analysisOut = (ValueData *) malloc(nAnalOut*sizeof(ValueData));
            if (aims[i].analysisOut == NULL) {
                printf(" CAPS Error: Name Malloc for %s -- %d AnalysisOut!\n",
                       aName, nAnalOut);
                status = EGADS_MALLOC;
                goto cleanup;
            }
            for (j = 0; j < nAnalOut; j++) {
                status = caps_childByIndex(analysisObj, VALUE, ANALYSISOUT, j+1,
                                           &valueObj);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %s ret=%d from caps_child(AnalOut,%d)!\n",
                           aName, status, j+1);
                    goto cleanup;
                }
                status = caps_info(valueObj, &vName, &type, &subtype, &dataLink,
                                   &parent, &own);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %s ret=%d from caps_info(AnalOut,%d)!\n",
                           vName, status, j+1);
                    goto cleanup;
                }
/*@-kepttrans@*/
                aims[i].analysisOut[j].name = vName;
/*@+kepttrans@*/
                if (aims[i].dirty) {
                    aims[i].analysisOut[j].nulldata = true;
                } else {
                    status = getValueData(&valueObj, &aims[i].analysisOut[j]);
                    if (status != CAPS_SUCCESS) goto cleanup;
                }
            }
        }
        if (nDynAnalOut > 0) {
            aims[i].analysisDynO = (ValueData *) malloc(nDynAnalOut*sizeof(ValueData));
            if (aims[i].analysisDynO == NULL) {
                printf(" CAPS Error: Name Malloc for %s -- %d AnalysisDynOut!\n",
                       aName, nDynAnalOut);
                status = EGADS_MALLOC;
                goto cleanup;
            }
            for (j = 0; j < nDynAnalOut; j++) {
                status = caps_childByIndex(analysisObj, VALUE, ANALYSISDYNO, j+1,
                                           &valueObj);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %s ret=%d from caps_child(AnalDynOut,%d)!\n",
                           aName, status, j+1);
                    goto cleanup;
                }
                status = caps_info(valueObj, &vName, &type, &subtype, &dataLink,
                                   &parent, &own);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %s ret=%d from caps_info(AnalDynOut,%d)!\n",
                           vName, status, j+1);
                    goto cleanup;
                }
/*@-kepttrans@*/
                aims[i].analysisDynO[j].name = vName;
/*@+kepttrans@*/
                status = getValueData(&valueObj, &aims[i].analysisDynO[j]);
                if (status != CAPS_SUCCESS) goto cleanup;
            }
        }
    }

    /* loop through bounds */
    if (nBound > 0) {
        for (i = 0; i < nBound; i++) {
            status = caps_childByIndex(problemObj, BOUND, NONE, i+1, &boundObject);
            if (status != CAPS_SUCCESS) {
                printf(" CAPS Error: Problem Object ret=%d from caps_child(Bound,%d)!\n",
                       status, i+1);
                goto cleanup;
            }
            status = caps_info(boundObject, &bName, &type, &subtype, &link,
                               &parent, &own);
            if (status != CAPS_SUCCESS) {
                printf(" CAPS Error: Bound Object ret=%d from caps_info,%d)!\n",
                       status, i+1);
                goto cleanup;
            }
            status = caps_size(boundObject, VERTEXSET, CONNECTED, &nConnect,
                               &nErr, &errors);
            if (errors != NULL) caps_freeError(errors);
            if (status != CAPS_SUCCESS) {
                printf(" CAPS Error: Bound %s ret=%d from caps_size(VSconnected)!\n",
                       bName, status);
                goto cleanup;
            }

            if (nConnect == 0) continue;
            for (j = 0; j < nConnect; j++) {
                status = caps_childByIndex(boundObject, VERTEXSET, CONNECTED,
                                           j+1, &vertexObj);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %s ret=%d from caps_child(VSConnected,%d)!\n",
                           bName, status, j+1);
                    goto cleanup;
                }
                
                status = caps_vertexSetInfo(vertexObj, &nGpts, &nDpts, &bObj,
                                            &aObj);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: Obj %d ret=%d from caps_vertesSetInfo!\n",
                           j+1, status);
                    goto cleanup;
                }
                status = caps_size(vertexObj, DATASET, NONE, &nDataSet, &nErr,
                                   &errors);
                if (status != CAPS_SUCCESS) {
                    printf(" CAPS Error: VSObj %d ret=%d from caps_size(DataSet)!\n",
                           j+1, status);
                    goto cleanup;
                }
                if (nDataSet == 0) continue;
                for (k = 0; k < nDataSet; k++) {
                    status = caps_childByIndex(vertexObj, DATASET, NONE, k+1,
                                               &dataSetObj);
                    if (status != CAPS_SUCCESS) {
                        printf(" CAPS Error: DSObj %d ret=%d from CBI(DataSet)!\n",
                               k+1, status);
                        goto cleanup;
                    }
                    status = caps_dataSetInfo(dataSetObj, &ftype, &dataLink,
                                              &dmethod);
                    if (status != CAPS_SUCCESS) {
                        printf(" CAPS Error: DSObj %d ret=%d from dataSetInfo!\n",
                               k+1, status);
                        goto cleanup;
                    }
                    if (ftype == FieldIn) {
                        status = caps_info(dataSetObj, &dName, &type, &subtype,
                                           &link, &parent, &own);
                        if (status != CAPS_SUCCESS) {
                            printf(" CAPS Error: DSObj %d ret=%d from caps_info!\n",
                                   k+1, status);
                            goto cleanup;
                        }
                        status = caps_info(dataLink, &name, &type, &subtype,
                                           &link, &parent, &own);
                        if (status != CAPS_SUCCESS) {
                            printf(" CAPS Error: DSLnk %d ret=%d from caps_info!\n",
                                   k+1, status);
                            goto cleanup;
                        }
                        tmpGeomLink = (GeomLinkData *) realloc(geomLinks,
                                     (geomLinksIndex+1) * sizeof(GeomLinkData));
                        if (tmpGeomLink == NULL) {
                            printf(" CAPS Error: Link Malloc for %s -- %d Xfers!\n",
                                   bName, geomLinksIndex+1);
                            status = EGADS_MALLOC;
                            goto cleanup;
                        }
                        geomLinks = tmpGeomLink;
                        geomLinks[geomLinksIndex].sourceAim = parent->name;
                        geomLinks[geomLinksIndex].targetAim = aObj->name;
/*@-kepttrans@*/
                        geomLinks[geomLinksIndex].sourceVar = name;
                        geomLinks[geomLinksIndex].targetVar = dName;
                        geomLinks[geomLinksIndex].bound     = bName;
/*@+kepttrans@*/
                        geomLinksIndex++;
                    }
                } 
            }
        }
    }
    
    /* make the stream */
    status   = EGADS_MALLOC;
    jsonText = (char *) EG_alloc(30 * sizeof(char));
    if (jsonText == NULL) goto cleanup;
    sprintf(jsonText, "dataJSON = {\n\"aims\": [\n");

    for (i = 0; i < nAnalysis; i++) {
        char *idStr, *inStr, *outStr, *dynOStr, *aimStr;
        char dirtyStr[20];

        idStr   = (char *) malloc((20 + strlen(aims[i].id)   * sizeof(char)));
        inStr   = (char *) malloc( 20 * sizeof(char));
        outStr  = (char *) malloc( 20 * sizeof(char));
        dynOStr = (char *) malloc( 20 * sizeof(char));
        if ((idStr   == NULL) || (inStr == NULL) || (outStr == NULL) ||
            (dynOStr == NULL)) {
            if (idStr   != NULL) free(idStr);
            if (inStr   != NULL) free(inStr);
            if (outStr  != NULL) free(outStr);
            if (dynOStr != NULL) free(dynOStr);
            goto cleanup;
        }

        sprintf(idStr,    "\"id\": \"%s\"",  aims[i].id);
        sprintf(dirtyStr, "\"dirty\": %s",   aims[i].dirty ? "true" : "false");
        sprintf(inStr,    "\"inVars\" : [");
        
        for (j = 0; j < aims[i].inCount; j++) {
            if (j > 0) strcat(inStr, ",\n");
            varStr = valueString(&aims[i].analysisIn[j]);
            if (varStr == NULL) {
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            tmpStr = (char *) realloc(inStr, (strlen(inStr)+strlen(varStr)+10)*
                                      sizeof(char));
            if (tmpStr == NULL) {
                free(varStr);
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            inStr = tmpStr;
            strcat(inStr, varStr);
            free(varStr);
        }
        strcat(inStr, "]");

        sprintf(outStr, "\"outVars\" : [");
        for (j = 0; j < aims[i].outCount; j++) {
            if (j > 0) strcat(outStr, ",\n");
            varStr = valueString(&aims[i].analysisOut[j]);
            if (varStr == NULL) {
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            
            tmpStr = (char *) realloc(outStr, (strlen(outStr)+strlen(varStr)+10)*
                                      sizeof(char));
            if (tmpStr == NULL) {
                free(varStr);
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            outStr = tmpStr;
            strcat(outStr, varStr);
            free(varStr);
        }
        strcat(outStr, "]");

        sprintf(dynOStr, "\"dynOutVars\" : [");
        for (j = 0; j < aims[i].dynOCount; j++) {
            if (j > 0) strcat(dynOStr, ",\n");
            varStr = valueString(&aims[i].analysisDynO[j]);
            if (varStr == NULL) {
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            tmpStr = (char *) realloc(dynOStr, (strlen(dynOStr)+strlen(varStr)+10)*
                                      sizeof(char));
            if (tmpStr == NULL) {
                free(varStr);
                free(idStr);
                free(inStr);
                free(outStr);
                free(dynOStr);
                goto cleanup;
            }
            dynOStr = tmpStr;
            strcat(dynOStr, varStr);
            free(varStr);
        }
        strcat(dynOStr, "]");

        aimStr = (char *) malloc((strlen(idStr) +strlen(dirtyStr)+strlen(inStr) +
                                  strlen(outStr)+strlen(dynOStr)+20)*sizeof(char));
        if (aimStr == NULL) {
            free(idStr);
            free(inStr);
            free(outStr);
            free(dynOStr);
            goto cleanup;
        }
        sprintf(aimStr, "{ %s,\n %s,\n %s,\n %s,\n %s}",
                idStr, dirtyStr, inStr, outStr, dynOStr);
        
        tmpStr = (char *) EG_reall(jsonText, (strlen(jsonText)+strlen(aimStr) +
                                              30) * sizeof(char));
        if (tmpStr == NULL) {
            free(aimStr);
            free(idStr);
            free(inStr);
            free(outStr);
            free(dynOStr);
            goto cleanup;
        }
        jsonText = tmpStr;
        if (i > 0) strcat(jsonText, ",\n");
        strcat(jsonText, aimStr);
        free(aimStr);
        free(idStr);
        free(inStr);
        free(outStr);
        free(dynOStr);
    }

/*@-nullderef@*/
    strcat(jsonText, "],\n \"valLinks\": [\n");
    for (i = 0; i < valLinksIndex; i++) {
        char *sourceStr, *targetStr, *dataStr, *linkStr;

        sourceStr = (char *) malloc((strlen(valLinks[i].sourceAim) + 20) *
                                    sizeof(char));
        targetStr = (char *) malloc((strlen(valLinks[i].targetAim) + 20) *
                                    sizeof(char));
        dataStr   = (char *) malloc((strlen(valLinks[i].sourceVar) +
                                     strlen(valLinks[i].targetVar) + 60) *
                                    sizeof(char));
        if ((sourceStr == NULL) || (targetStr == NULL) || (dataStr == NULL)) {
            if (sourceStr != NULL) free(sourceStr);
            if (targetStr != NULL) free(targetStr);
            if (dataStr   != NULL) free(dataStr);
            goto cleanup;
        }
        sprintf(sourceStr, "\"source\": \"%s\"", valLinks[i].sourceAim);
        sprintf(targetStr, "\"target\": \"%s\"", valLinks[i].targetAim);
        sprintf(dataStr,   "\"data\": [{\"sourceVar\": \"%s\", \"targetVar\": \"%s\"}]",
                valLinks[i].sourceVar, valLinks[i].targetVar);

        linkStr = (char *) malloc((strlen(sourceStr) + strlen(targetStr) +
                                   strlen(dataStr) + 20) * sizeof(char));
        if (linkStr == NULL) {
            free(sourceStr);
            free(targetStr);
            free(dataStr);
            goto cleanup;
        }
        sprintf(linkStr, "{ %s,\n %s,\n %s }", sourceStr, targetStr, dataStr);

        tmpStr = (char *) EG_reall(jsonText, (strlen(jsonText)+strlen(linkStr)+
                                              30) * sizeof(char));
        if (tmpStr == NULL) {
            free(linkStr);
            free(sourceStr);
            free(targetStr);
            free(dataStr);
            goto cleanup;
        }
        jsonText = tmpStr;
        if (i > 0) strcat(jsonText, ",\n");
        strcat(jsonText, linkStr);
        free(linkStr);
        free(sourceStr);
        free(targetStr);
        free(dataStr);
    }
    
    strcat(jsonText, "],\n \"geomLinks\": [\n");
    for (i = 0; i < geomLinksIndex; i++) {
        char *sourceStr, *targetStr, *dataStr, *linkStr;
      
        sourceStr = (char *) malloc((strlen(geomLinks[i].sourceAim) + 20) *
                                    sizeof(char));
        targetStr = (char *) malloc((strlen(geomLinks[i].targetAim) + 20) *
                                    sizeof(char));
        dataStr   = (char *) malloc((strlen(geomLinks[i].sourceVar) +
                                     strlen(geomLinks[i].targetVar) + 80) *
                                    sizeof(char));
        if ((sourceStr == NULL) || (targetStr == NULL) || (dataStr == NULL)) {
            if (sourceStr != NULL) free(sourceStr);
            if (targetStr != NULL) free(targetStr);
            if (dataStr   != NULL) free(dataStr);
            goto cleanup;
        }
        sprintf(sourceStr, "\"source\": \"%s\"", geomLinks[i].sourceAim);
        sprintf(targetStr, "\"target\": \"%s\"", geomLinks[i].targetAim);
        sprintf(dataStr,   "\"data\": [{\"bound\": \"%s\", \"sourceVar\": \"%s\", \"targetVar\": \"%s\"}]",
                geomLinks[i].bound, geomLinks[i].sourceVar, geomLinks[i].targetVar);
        
        linkStr = (char *) malloc((strlen(sourceStr) + strlen(targetStr) +
                                   strlen(dataStr)   + 20) * sizeof(char));
        if (linkStr == NULL) {
            free(sourceStr);
            free(targetStr);
            free(dataStr);
            goto cleanup;
        }
        sprintf(linkStr, "{ %s,\n %s,\n %s }", sourceStr, targetStr, dataStr);
        
        tmpStr = (char *) EG_reall(jsonText, (strlen(jsonText)+strlen(linkStr)+
                                              30) * sizeof(char));
        if (tmpStr == NULL) {
            free(linkStr);
            free(sourceStr);
            free(targetStr);
            free(dataStr);
            goto cleanup;
        }
        jsonText = tmpStr;
        if (i > 0) strcat(jsonText, ",\n");
        strcat(jsonText, linkStr);
        free(sourceStr);
        free(targetStr);
        free(dataStr);
        free(linkStr);
    }
/*@+nullderef@*/
    strcat(jsonText, "] \n }");
    
    if (stream == NULL) {
        fptr = fopen(filename, "w");
        if (fptr == NULL) {
            printf(" CAPS Error: error opening file %s\n", filename);
            status = CAPS_IOERR;
            goto cleanup;
        }
        fprintf(fptr, "%s", jsonText);
        fclose(fptr);
  
        // start browser and wait until done
        system(env);
    }
#ifdef DEBUG
    printf("\n");
#endif
    status = CAPS_SUCCESS;
  
cleanup:
    if ((stream != NULL) && (status == CAPS_SUCCESS)) {
        *stream = jsonText;
    } else {
        EG_free(jsonText);
    }
    if (geomLinks != NULL) free(geomLinks);
    if (valLinks  != NULL) free(valLinks);
    if (aims != NULL) {
        for (i = 0; i < nAnalysis; i++) {
            if (aims[i].analysisDynO != NULL) free(aims[i].analysisDynO);
            if (aims[i].analysisOut  != NULL) free(aims[i].analysisOut);
            if (aims[i].analysisIn   != NULL) free(aims[i].analysisIn);
        }
        free(aims);
    }
  
    /* return of caps_debug(problemObj) should be 0 */
    i = caps_debug(problemObj);
    if (i != 0) {
        printf(" CAPS_Error: CAPS still in debug mode!\n");
        if ((stream != NULL) && (status == CAPS_SUCCESS)) {
            EG_free(*stream);
            *stream = NULL;
        }
        return CAPS_STATEERR;
    }

    return status;
}
