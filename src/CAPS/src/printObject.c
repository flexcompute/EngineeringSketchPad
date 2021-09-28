/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Object Output Utility
 *
 *      Copyright 2014-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "caps.h"


extern void caps_printObjects(capsObj object, int indent);


static int
printValues(capsObj object, int indent)
{
  int            i, j, status, len, nerr, dim, nrow, ncol, range, rank, pmtr;
  int            nLines, *ints;
  char           *str, *pname, *pID, *uID, **dotnames, **lines;
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
  status = caps_ownerInfo(own, &pname, &pID, &uID, &nLines, &lines, datetime,
                          &sNum);
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
caps_printObjects(capsObj object, int indent)
{
  int            i, j, k, status, nBody, nParam, nGeomIn, nGeomOut;
  int            nAnalysis, nBound, nAnalIn, nAnalOut, nConnect, nUnConnect;
  int            nAttr, nDataSet, npts, rank, nErr, nLines;
  char           *name, *units, *pname, *pID, *userID, *oname, **lines;
  double         *data;
  short          datetime[6];
  CAPSLONG       sn, lsn;
  capsObj        link, parent, obj, attr;
  capsOwn        own;
  capsErrs       *errors;
  enum capsoType type, otype;
  enum capssType subtype, osubtype;
  static char    *oType[9]  = { "BODIES", "ATTRIBUTES", "UNUSED", "PROBLEM",
                                "VALUE", "ANALYSIS", "BOUND", "VERTEXSET",
                                "DATASET" };
  static char    *sType[11] = { "NONE", "STATIC", "PARAMETRIC", "GEOMETRYIN",
                                "GEOMETRYOUT", "PARAMETER", "USER",
                                "ANALYSISIN", "ANALYSISOUT", "CONNECTED",
                                "UNCONNECTED" };

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
  status = caps_ownerInfo(own, &pname, &pID, &userID, &nLines, &lines, datetime,
                          &sn);
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
      for (k = 0; k < nLines; k++) {
        for (i = 0; i < indent; i++) printf(" ");
        printf("         %s\n", lines[k]);
      }
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
    caps_printObjects(attr, indent+2);
    status = caps_delete(attr);
    if (status != CAPS_SUCCESS)
      printf(" CAPS Error: Object %s returns %d from caps_delete(Attr)!\n",
             attr->name, status);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
      }
    }

  } else if (type == VALUE) {
    if (link != NULL) {
      for (i = 0; i < indent+2; i++) printf(" ");
      printf(" linked to %s\n", link->name);
    } else {
      status = printValues(object, indent+2);
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
    for (i = 0; i < indent; i++) printf(" ");
    printf("   %d AnalysisIns, %d AnalysisOuts\n", nAnalIn, nAnalOut);

    if (nAnalIn > 0) {
      printf("\n");
      for (i = 0; i < nAnalIn; i++) {
        status = caps_childByIndex(object, VALUE, ANALYSISIN, i+1, &obj);
        if (status != CAPS_SUCCESS) {
          printf(" CAPS Error: Object %s ret=%d from caps_child(AnalIn,%d)!\n",
                 name, status, i+1);
          return;
        }
        caps_printObjects(obj, indent+2);
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
          caps_ownerInfo(own, &pname, &pID, &userID, &nLines, &lines, datetime,
                         &sn);
        if (sn == 0) continue;
        
        caps_printObjects(obj, indent+2);

        status = caps_info(obj, &oname, &otype, &osubtype, &link, &parent, &own);
        if (status == CAPS_SUCCESS) {
          status = caps_ownerInfo(own, &pname, &pID, &userID, &nLines, &lines,
                                  datetime, &lsn);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
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
        caps_printObjects(obj, indent+2);
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
