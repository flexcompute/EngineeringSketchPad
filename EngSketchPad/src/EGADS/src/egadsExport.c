/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Export a Model (via a string) for use in egadsLite
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egads.h"
/*@-redef@*/
typedef int    INT_;
typedef INT_   INT_3D[3];
typedef double DOUBLE_2D[2];
/*@+redef@*/
#include "uvmap_struct.h"

  extern int EG_outLevel( const egObject *object );
  extern int EG_flattenBSpline( egObject *object, egObject **result );


typedef struct {
  int       nobjs;                /* number in the map */
  egObject **objs;                /* vector of egos for map */
} egMap;


typedef struct {
  egMap pcurves;
  egMap curves;
  egMap surfaces;
} egGeoMap;


typedef struct {
  void    *data;
  size_t  ptr;
  size_t  size;
} stream_T;


#define CHUNK 10000
  
static int
Fwrite(void *data, size_t size, int nitems, stream_T *stream)
{
  void *temp_data;

  while (stream->ptr + size*nitems > stream->size) {
    stream->size += CHUNK;
    temp_data = EG_reall(stream->data, stream->size);
    if (temp_data == NULL) return -1;
    stream->data = temp_data;
  }

  memcpy(&(((char *)stream->data)[stream->ptr]), data, size*nitems);
  stream->ptr += size*nitems;

  return nitems;
}


static void
Fclose(stream_T *stream)
{
  if (stream->data != NULL) EG_free(stream->data);
}


static int
EG_writeString(stream_T *fp, /*@null@*/ const char *string)
{
  int    len = 0;
  size_t n;
  
  if (string != NULL) len = strlen(string) + 1;
  n = Fwrite(&len, sizeof(int), 1, fp);
  if (n != 1) return EGADS_WRITERR;
  
  if (string == NULL) return EGADS_SUCCESS;
  n = Fwrite((void *) string, sizeof(char), len, fp);
  if (n != len) return EGADS_WRITERR;
  
  return EGADS_SUCCESS;
}


static int
EG_writeAttrs(stream_T *fp, egAttrs *attrs)
{
  int     nattr, i, stat;
  size_t  n;
  egAttr  *attr;
  
  nattr = 0;
  if (attrs != NULL) {
    for (i = 0; i < attrs->nattrs; i++)
      if (attrs->attrs[i].type != ATTRPTR) nattr++;
  }
  n = Fwrite(&nattr, sizeof(int), 1, fp);
  if  (n     != 1) return EGADS_WRITERR;
  if ((nattr == 0) || (attrs == NULL)) return EGADS_SUCCESS;
  
  attr = attrs->attrs;
  for (i = 0; i < nattr; i++) {
    n = Fwrite(&attr[i].type,   sizeof(int), 1, fp);
    if (n != 1) return EGADS_WRITERR;
    n = Fwrite(&attr[i].length, sizeof(int), 1, fp);
    if (n != 1) return EGADS_WRITERR;
    stat = EG_writeString(fp, attr[i].name);
    if (stat != EGADS_SUCCESS) return EGADS_WRITERR;
    if (attr[i].type == ATTRINT) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = Fwrite(&attr[i].vals.integer, sizeof(int),              1, fp);
      } else if (attr[i].length > 1) {
        n = Fwrite(attr[i].vals.integers, sizeof(int), attr[i].length, fp);
      }
      if (n != attr[i].length) return EGADS_WRITERR;
    } else if ((attr[i].type == ATTRREAL) || (attr[i].type == ATTRCSYS)) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = Fwrite(&attr[i].vals.real, sizeof(double),              1, fp);
      } else if (attr[i].length > 1) {
        n = Fwrite(attr[i].vals.reals, sizeof(double), attr[i].length, fp);
      }
      if (n != attr[i].length) return EGADS_WRITERR;
    } else if (attr[i].type == ATTRSTRING) {
      stat = EG_writeString(fp, attr[i].vals.string);
      if (stat != EGADS_SUCCESS) return EGADS_WRITERR;
    }
  }

  return EGADS_SUCCESS;
}


static int
EG_lookAtMap(egObject *object, int oclass, const egGeoMap *maps, int flag)
{
  int   i;
  egMap map;
  
  if (oclass == PCURVE) {
    map = maps->pcurves;
  } else if (oclass == CURVE) {
    map = maps->curves;
  } else if (oclass == SURFACE) {
    map = maps->surfaces;
  } else {
    printf(" Bad Geometry type: %d\n", oclass);
    return EGADS_GEOMERR;
  }
  for (i = 0; i < map.nobjs; i++)
    if (object == map.objs[i]) return i+1;
  
  for (i = 0; i < map.nobjs; i++)
    if (EG_isSame(object, map.objs[i]) == EGADS_SUCCESS) return i+1;

  if (flag == 0)
    printf(" Geometry type: %d -- Not found in %d objs!\n", oclass, map.nobjs);
  return EGADS_NOTFOUND;
}


static int
EG_populateGeom(egObject *bobject, egGeoMap *maps)
{
  int      i, j, outLevel, stat, n, oclass, mtype, nchild, npcrv, ncrv, nsurf;
  int      *senses, *ivec;
  double   data[4], *rvec;
  egObject *obj, **objs, *ref, **children;
  
  npcrv    = ncrv = nsurf = 0;
  outLevel = EG_outLevel(bobject);
  
  /* only see pcurves in loops */
  stat = EG_getBodyTopos(bobject, NULL, LOOP, &n, &objs);
  if (stat != EGADS_SUCCESS) return stat;
  for (i = 0; i < n; i++) {
    obj  = objs[i];
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      return stat;
    }
    if (ref == NULL) continue;
    npcrv += nchild;
    nsurf++;
    /* do surfaces reference surfaces? */
    do {
      obj = ref;
      stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) {
        EG_free(objs);
        return stat;
      }
      if (ivec != NULL) EG_free(ivec);
      EG_free(rvec);
      if (ref != NULL)
        if (ref->oclass == CURVE) ref = NULL;
      if (ref != NULL) nsurf++;
    } while (ref != NULL);
    /* do any of the pcurves have reference geometry? */
    for (j = 0; j < nchild; j++) {
      obj = children[j+nchild];
      do {
        stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) break;
        if (ivec != NULL) EG_free(ivec);
        EG_free(rvec);
        if (ref  != NULL) npcrv++;
        obj = ref;
      } while (obj != NULL);
    }
  }
  if (npcrv != 0) {
    maps->pcurves.objs = (egObject **) EG_alloc(npcrv*sizeof(egObject *));
    if (maps->pcurves.objs == NULL) {
      EG_free(objs);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++) {
      obj = objs[i];
      EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                     &senses);
      if (ref == NULL) continue;
      for (j = 0; j < nchild; j++) {
        obj  = children[j+nchild];
        stat = EG_lookAtMap(obj, PCURVE, maps, 1);
        if (stat == EGADS_NOTFOUND) {
          maps->pcurves.objs[maps->pcurves.nobjs] = obj;
          maps->pcurves.nobjs++;
        }
        do {
          stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
          if (stat != EGADS_SUCCESS) break;
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          if (ref  != NULL) {
            stat = EG_lookAtMap(ref, PCURVE, maps, 1);
            if (stat == EGADS_NOTFOUND) {
              maps->pcurves.objs[maps->pcurves.nobjs] = ref;
              maps->pcurves.nobjs++;
            }
          }
          obj = ref;
        } while (obj != NULL);
      }
    }
    if (outLevel > 0)
      printf("   PCurve Map: %d (alloc %d)!\n", maps->pcurves.nobjs, npcrv);
  }
  EG_free(objs);

  /* look for Surfaces in Loops (already counted) and Faces */
  stat = EG_getBodyTopos(bobject, NULL, FACE, &n, &objs);
  if (stat != EGADS_SUCCESS) return stat;
  for (i = 0; i < n; i++) {
    obj  = objs[i];
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      return stat;
    }
    if (ref == NULL) continue;
    nsurf++;
    /* does the surface have reference geometry? */
    while (ref != NULL) {
      obj  = ref;
      stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) break;
      if (ivec != NULL) EG_free(ivec);
      EG_free(rvec);
      if (ref != NULL) {
        if (ref->oclass == CURVE) {
          ncrv++;
        } else if ((ref->mtype == EXTRUSION) || (ref->mtype == REVOLUTION)) {
          ncrv++;
        } else {
          nsurf++;
        }
      }
    }
  }
  if (nsurf != 0) {
    maps->surfaces.objs = (egObject **) EG_alloc(nsurf*sizeof(egObject *));
    if (maps->surfaces.objs == NULL) {
      EG_free(objs);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++) {
      obj = objs[i];
      EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                     &senses);
      if (ref == NULL) continue;
      stat = EG_lookAtMap(ref, SURFACE, maps, 1);
      if (stat == EGADS_NOTFOUND) {
        maps->surfaces.objs[maps->surfaces.nobjs] = ref;
        maps->surfaces.nobjs++;
        do {
          obj = ref;
          EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          if (ref != NULL)
              if (ref->oclass == CURVE) ref = NULL;
          if (ref != NULL) {
            stat = EG_lookAtMap(ref, SURFACE, maps, 1);
            if (stat == EGADS_NOTFOUND) {
              maps->surfaces.objs[maps->surfaces.nobjs] = ref;
              maps->surfaces.nobjs++;
            }
          }
        } while (ref != NULL);
      }
    }
    EG_free(objs);
    EG_getBodyTopos(bobject, NULL, LOOP, &n, &objs);
    for (i = 0; i < n; i++) {
      obj = objs[i];
      EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                     &senses);
      if (ref == NULL) continue;
      do {
        obj = ref;
        EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
        if (ivec != NULL) EG_free(ivec);
        EG_free(rvec);
        if (ref != NULL)
          if (ref->oclass == CURVE) ref = NULL;
        if (ref != NULL) {
          stat = EG_lookAtMap(ref, SURFACE, maps, 1);
          if (stat == EGADS_NOTFOUND) {
            maps->surfaces.objs[maps->surfaces.nobjs] = ref;
            maps->surfaces.nobjs++;
          }
        }
      } while (ref != NULL);
    }
    if (outLevel > 0)
      printf("   Surface Map: %d (alloc %d)!\n", maps->surfaces.nobjs, nsurf);
  }
  EG_free(objs);
  
  /* look for Curves in Edges & as reference for some Surfaces */
  
  stat = EG_getBodyTopos(bobject, NULL, EDGE, &n, &objs);
  if (stat != EGADS_SUCCESS) return stat;
  for (i = 0; i < n; i++) {
    obj  = objs[i];
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      return stat;
    }
    if ((ref == NULL) || (mtype == DEGENERATE)) continue;
    ncrv++;
    /* does the Curve have reference geometry? */
    do {
      obj  = ref;
      stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) break;
      if (ivec != NULL) EG_free(ivec);
      EG_free(rvec);
      if (ref  != NULL) ncrv++;
    } while (ref != NULL);
  }
  if (ncrv != 0) {
    maps->curves.objs = (egObject **) EG_alloc(ncrv*sizeof(egObject *));
    if (maps->curves.objs == NULL) {
      EG_free(objs);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++) {
      obj = objs[i];
      EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                     &senses);
      if ((ref == NULL) || (mtype == DEGENERATE)) continue;
      stat = EG_lookAtMap(ref, CURVE, maps, 1);
      if (stat == EGADS_NOTFOUND) {
        maps->curves.objs[maps->curves.nobjs] = ref;
        maps->curves.nobjs++;
      }
      /* does the Curve have reference geometry? */
      do {
        obj  = ref;
        stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) break;
        if (ivec != NULL) EG_free(ivec);
        EG_free(rvec);
        if (ref  != NULL) {
          stat = EG_lookAtMap(ref, CURVE, maps, 1);
          if (stat == EGADS_NOTFOUND) {
            maps->curves.objs[maps->curves.nobjs] = ref;
            maps->curves.nobjs++;
          }
        }
      } while (ref != NULL);
    }
    EG_free(objs);
    EG_getBodyTopos(bobject, NULL, FACE, &n, &objs);
    for (i = 0; i < n; i++) {
      obj  = objs[i];
      EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                     &senses);
      /* does the surface have curve reference geometry? */
      if (ref == NULL) continue;
      do {
        obj  = ref;
        stat = EG_getGeometry(obj, &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) break;
        if (ivec != NULL) EG_free(ivec);
        EG_free(rvec);
        if (ref != NULL) {
          if (ref->oclass == CURVE) {
            stat = EG_lookAtMap(ref, CURVE, maps, 1);
            if (stat == EGADS_NOTFOUND) {
              maps->curves.objs[maps->curves.nobjs] = ref;
              maps->curves.nobjs++;
            }
          } else if ((ref->mtype == EXTRUSION) || (ref->mtype == REVOLUTION)) {
            stat = EG_lookAtMap(ref, CURVE, maps, 1);
            if (stat == EGADS_NOTFOUND) {
              maps->curves.objs[maps->curves.nobjs] = ref;
              maps->curves.nobjs++;
            }
          }
        }
      } while (ref != NULL);
    }

    if (outLevel > 0)
      printf("   Curve Map: %d (alloc %d)!\n", maps->curves.nobjs, ncrv);
  }
  EG_free(objs);

  
  return EGADS_SUCCESS;
}


static int
EG_writeGeometry(egObject *gobject, const egGeoMap *maps, stream_T *fp)
{
  int      stat, nreal, nint, iref, oclass, mtype, *ivec;
  double   *rvec;
  egObject *robject, *bspline;
  size_t   n;

  stat = EG_getGeometry(gobject, &oclass, &mtype, &robject, &ivec, &rvec);
  if (stat != EGADS_SUCCESS) return stat;

  iref = nreal = nint = 0;
  if (oclass == PCURVE) {
    switch (mtype) {
        
      case LINE:
        nreal = 4;
        break;
        
      case CIRCLE:
        nreal = 7;
        break;
        
      case ELLIPSE:
        nreal = 8;
        break;
        
      case PARABOLA:
        nreal = 7;
        break;
        
      case HYPERBOLA:
        nreal = 8;
        break;
        
      case TRIMMED:
        iref = EG_lookAtMap(robject, PCURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 2;
        break;
        
      case BEZIER:
        nint  = 3;
        nreal = 2*ivec[2];
        if ((ivec[0]&2) != 0) nreal += ivec[2];
        break;
        
      case BSPLINE:
        nint  = 4;
        nreal = ivec[3] + 2*ivec[2];
        if ((ivec[0]&2) != 0) nreal += ivec[2];
        if ((ivec[0]&4) != 0) {
          printf(" EGADS Warning: Periodic PCurve!\n");
        }
        break;
        
      case OFFSET:
        iref = EG_lookAtMap(robject, PCURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 1;
        break;
    }
    
  } else if (oclass == CURVE) {
    switch (mtype) {
        
      case LINE:
        nreal = 6;
        break;
        
      case CIRCLE:
        nreal = 10;
        break;
        
      case ELLIPSE:
        nreal = 11;
        break;
        
      case PARABOLA:
        nreal = 10;
        break;
        
      case HYPERBOLA:
        nreal = 11;
        break;
        
      case TRIMMED:
        iref = EG_lookAtMap(robject, CURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 2;
        break;
        
      case BEZIER:
        nint  = 3;
        nreal = 3*ivec[2];
        if ((ivec[0]&2) != 0) nreal += ivec[2];
        break;
        
      case BSPLINE:
        nint  = 4;
        nreal = ivec[3] + 3*ivec[2];
        if ((ivec[0]&2) != 0) nreal += ivec[2];
        if ((ivec[0]&4) != 0) {
/*        printf(" note: Converting Periodic Curve!\n");  */
          stat = EG_flattenBSpline(gobject, &bspline);
          EG_free(ivec);
          EG_free(rvec);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_flattenBSpline = %d\n", stat);
            return stat;
          }
          stat = EG_getGeometry(bspline, &oclass, &mtype, &robject,
                                &ivec, &rvec);
          if (stat != EGADS_SUCCESS) return stat;
          nint  = 4;
          nreal = ivec[3] + 3*ivec[2];
          if ((ivec[0]&2) != 0) nreal += ivec[2];
          ivec[0] |= 1;           /* mark that this has been flattened */
/*        printf("     periodic flag = %d\n", ivec[0]);  */
          EG_deleteObject(bspline);
        }
        break;
        
      case OFFSET:
        iref = EG_lookAtMap(robject, CURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 4;
        break;
    }

  } else {
    /* surface */
    switch (mtype) {
        
      case PLANE:
        nreal = 9;
        break;
        
      case SPHERICAL:
        nreal = 10;
        break;
        
      case CONICAL:
        nreal = 14;
        break;
        
      case CYLINDRICAL:
        nreal = 13;
        break;
        
      case TOROIDAL:
        nreal = 14;
        break;
        
      case REVOLUTION:
        iref = EG_lookAtMap(robject, CURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        iref  = -iref;
        nreal = 6;
        break;
        
      case EXTRUSION:
        iref = EG_lookAtMap(robject, CURVE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        iref  = -iref;
        nreal = 3;
        break;
        
      case TRIMMED:
        iref = EG_lookAtMap(robject, SURFACE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 4;
        break;
        
      case BEZIER:
        nint  = 5;
        nreal = 3*ivec[2]*ivec[4];
        if ((ivec[0]&2) != 0) nreal += ivec[2]*ivec[4];
        break;
        
      case BSPLINE:
        nint  = 7;
        nreal = ivec[3] + ivec[6] + 3*ivec[2]*ivec[5];
        if ((ivec[0]&2) != 0) nreal += ivec[2]*ivec[5];
        if ((ivec[0]&12) != 0) {
/*        printf(" note: Converting Periodic Surface %d!\n", ivec[0]);  */
          stat = EG_flattenBSpline(gobject, &bspline);
          EG_free(ivec);
          EG_free(rvec);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_flattenBSpline = %d\n", stat);
            return stat;
          }
          stat = EG_getGeometry(bspline, &oclass, &mtype, &robject,
                                &ivec, &rvec);
          if (stat != EGADS_SUCCESS) return stat;
          nreal = ivec[3] + ivec[6] + 3*ivec[2]*ivec[5];
          if ((ivec[0]&2) != 0) nreal += ivec[2]*ivec[5];
          ivec[0] |= 1;           /* mark that this has been flattened */
 /*       printf("     periodic flag = %d\n", ivec[0]);  */
          EG_deleteObject(bspline);
        }
        break;
        
      case OFFSET:
        iref = EG_lookAtMap(robject, SURFACE, maps, 0);
        if (iref < EGADS_SUCCESS) {
          if (ivec != NULL) EG_free(ivec);
          EG_free(rvec);
          return iref;
        }
        nreal = 1;
        break;
    }
  }
  if (nreal == 0) {
    printf(" OCLASS = %d, MTYPE = %d not found!\n", oclass, mtype);
    if (ivec != NULL) EG_free(ivec);
    EG_free(rvec);
    return EGADS_GEOMERR;
  }
  
  n = Fwrite(&iref,  sizeof(int), 1, fp);
  if (n != 1) {
    if (ivec != NULL) EG_free(ivec);
    EG_free(rvec);
    return EGADS_WRITERR;
  }
  n = Fwrite(&nint,  sizeof(int), 1, fp);
  if (n != 1) {
    if (ivec != NULL) EG_free(ivec);
    EG_free(rvec);
    return EGADS_WRITERR;
  }
  n = Fwrite(&nreal, sizeof(int), 1, fp);
  if (n != 1) {
    if (ivec != NULL) EG_free(ivec);
    EG_free(rvec);
    return EGADS_WRITERR;
  }
  
  if (nint != 0) {
    n = Fwrite(ivec, sizeof(int), nint, fp);
    if (n != nint) {
      if (ivec != NULL) EG_free(ivec);
      EG_free(rvec);
      return EGADS_WRITERR;
    }
  }
  n = Fwrite(rvec, sizeof(double), nreal, fp);
  if (n != nreal) {
    if (ivec != NULL) EG_free(ivec);
    EG_free(rvec);
    return EGADS_WRITERR;
  }

  if (ivec != NULL) EG_free(ivec);
  EG_free(rvec);
  return EGADS_SUCCESS;
}


static int
EG_exportBody(egObject *bobject, stream_T *fp)
{
  int      i, n, m, no, nchild, stat, oclass, mtype, iref, ntypes[8], *senses;
  int      outLevel;
  double   tol, data[4], bbox[6], lims[4];
  egGeoMap maps;
  egObject **objs, *obj, *ref, **children;

  if (bobject == NULL)               return EGADS_NULLOBJ;
  if (bobject->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bobject->oclass != BODY)       return EGADS_NOTMODEL;

  outLevel = EG_outLevel(bobject);
  m = bobject->mtype;
  n = Fwrite(&m, sizeof(int), 1, fp);
  if (n != 1) return EGADS_WRITERR;

  stat = EG_getBodyTopos(bobject, NULL, NODE, &ntypes[3], &objs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_free(objs);
  stat = EG_getBodyTopos(bobject, NULL, EDGE, &ntypes[4], &objs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_free(objs);
  stat = EG_getBodyTopos(bobject, NULL, LOOP, &ntypes[5], &objs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_free(objs);
  ntypes[6] = ntypes[7] = 0;
  if (bobject->mtype != WIREBODY) {
    stat = EG_getBodyTopos(bobject, NULL, FACE, &ntypes[6], &objs);
    if (stat != EGADS_SUCCESS) return stat;
    EG_free(objs);
    if (bobject->mtype != FACEBODY) {
      stat = EG_getBodyTopos(bobject, NULL, SHELL, &ntypes[7], &objs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_free(objs);
    }
  }

  /* fill up our geometric entries */
  maps.pcurves.nobjs = maps.curves.nobjs = maps.surfaces.nobjs = 0;
  maps.pcurves.objs  = maps.curves.objs  = maps.surfaces.objs  = NULL;
  stat = EG_populateGeom(bobject, &maps);
  if (stat != EGADS_SUCCESS) {
    EG_free(maps.pcurves.objs);
    EG_free(maps.curves.objs);
    EG_free(maps.surfaces.objs);
    return stat;
  }
  ntypes[0] = maps.pcurves.nobjs;
  ntypes[1] = maps.curves.nobjs;
  ntypes[2] = maps.surfaces.nobjs;

  n = Fwrite(ntypes, sizeof(int), 8, fp);
  if (n != 8) {
    EG_free(maps.pcurves.objs);
    EG_free(maps.curves.objs);
    EG_free(maps.surfaces.objs);
    return EGADS_WRITERR;
  }
  
  /* pcurves */
  if (outLevel > 0)
    printf(" Writing PCurves...\n");
  for (i = 0; i < maps.pcurves.nobjs; i++) {
/*@-nullderef@*/
    obj = maps.pcurves.objs[i];
/*@+nullderef@*/
    n = Fwrite(&obj->mtype, sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    stat = EG_writeGeometry(obj, &maps, fp);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  
  /* curves */
  if (outLevel > 0)
    printf(" Writing Curves...\n");
  for (i = 0; i < maps.curves.nobjs; i++) {
/*@-nullderef@*/
    obj = maps.curves.objs[i];
/*@+nullderef@*/
    n = Fwrite(&obj->mtype, sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    stat = EG_writeGeometry(obj, &maps, fp);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  
  /* surfaces */
  if (outLevel > 0)
    printf(" Writing Surfaces...\n");
  for (i = 0; i < maps.surfaces.nobjs; i++) {
/*@-nullderef@*/
    obj = maps.surfaces.objs[i];
/*@+nullderef@*/
    n = Fwrite(&obj->mtype, sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    stat = EG_writeGeometry(obj, &maps, fp);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  
  /* nodes */
  EG_getBodyTopos(bobject, NULL, NODE, &no, &objs);
  if (outLevel > 0)
    printf(" Writing %d Nodes...\n", no);
  for (i = 0; i < no; i++) {
/*@-nullderef@*/
    obj  = objs[i];
/*@+nullderef@*/
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat == EGADS_SUCCESS) stat = EG_getTolerance(obj, &tol);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    n = Fwrite(data, sizeof(double), 3, fp);
    if (n != 3) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(&tol, sizeof(double), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  EG_free(objs);

  /* edges */
  EG_getBodyTopos(bobject, NULL, EDGE, &no, &objs);
  if (outLevel > 0)
    printf(" Writing %d Edges...\n", no);
  for (i = 0; i < no; i++) {
/*@-nullderef@*/
    obj  = objs[i];
/*@+nullderef@*/
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat == EGADS_SUCCESS) stat = EG_getTolerance  (obj, &tol);
    if (stat == EGADS_SUCCESS) stat = EG_getBoundingBox(obj, bbox);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    iref = 0;
    if (mtype != DEGENERATE) {
      if (ref->mtype == BSPLINE)
        if (EG_getRange(ref, lims, &n) == EGADS_SUCCESS)
          if (n != 0)
            if ((data[0] < lims[0]) || (data[1] > lims[1])) {
              printf(" EGADS Problem: Edge #%d BSPLINE limits out of Period!\n",
                     i+1);
              printf("                %lf %lf   %lf %lf\n",
                     data[0], data[1], lims[0], lims[1]);
            }
      iref = EG_lookAtMap(ref, CURVE, &maps, 0);
    }
    if (iref < EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return iref;
    }
    n = Fwrite(&mtype, sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(&iref,  sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    iref = EG_indexBodyTopo(bobject, children[0]);
    if (iref <= EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return iref;
    }
    n = Fwrite(&iref,  sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    if (nchild == 2) {
      iref = EG_indexBodyTopo(bobject, children[1]);
      if (iref <= EGADS_SUCCESS) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return iref;
      }
    }
    n = Fwrite(&iref,  sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    
    n = Fwrite(data, sizeof(double), 2, fp);
    if (n != 2) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(bbox, sizeof(double), 6, fp);
    if (n != 6) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(&tol,  sizeof(double), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  EG_free(objs);
  
  /* loops */
  EG_getBodyTopos(bobject, NULL, LOOP, &no, &objs);
  if (outLevel > 0)
    printf(" Writing %d Loops...\n", no);
  for (i = 0; i < no; i++) {
/*@-nullderef@*/
    obj  = objs[i];
/*@+nullderef@*/
    stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                          &senses);
    if (stat == EGADS_SUCCESS) stat = EG_getBoundingBox(obj, bbox);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
    n = Fwrite(&mtype,  sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(&nchild, sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    iref = 0;
    if (ref != NULL) iref = EG_lookAtMap(ref, SURFACE, &maps, 0);
    if (iref < EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return iref;
    }
    n = Fwrite(&iref,   sizeof(int), 1, fp);
    if (n != 1) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(bbox, sizeof(double), 6, fp);
    if (n != 6) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    n = Fwrite(senses, sizeof(int), nchild, fp);
    if (n != nchild) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return EGADS_WRITERR;
    }
    /* edges */
    for (m = 0; m < nchild; m++) {
      iref = EG_indexBodyTopo(bobject, children[m]);
      if (iref <= EGADS_SUCCESS) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return iref;
      }
      n = Fwrite(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
    }
    /* pcurves */
    if (ref != NULL)
      for (m = 0; m < nchild; m++) {
        iref = EG_lookAtMap(children[m+nchild], PCURVE, &maps, 0);
        if (iref < EGADS_SUCCESS) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return iref;
        }
        n = Fwrite(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
      }
    
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      EG_free(maps.pcurves.objs);
      EG_free(maps.curves.objs);
      EG_free(maps.surfaces.objs);
      return stat;
    }
  }
  EG_free(objs);

  if (bobject->mtype != WIREBODY) {
    /* faces */
    EG_getBodyTopos(bobject, NULL, FACE, &no, &objs);
    if (outLevel > 0)
      printf(" Writing %d Faces...\n", no);
    for (i = 0; i < no; i++) {
/*@-nullderef@*/
      obj  = objs[i];
/*@+nullderef@*/
      stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild, &children,
                            &senses);
      if (stat == EGADS_SUCCESS) stat = EG_getTolerance  (obj, &tol);
      if (stat == EGADS_SUCCESS) stat = EG_getBoundingBox(obj, bbox);
      if (stat != EGADS_SUCCESS) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return stat;
      }
      if (ref != NULL)
        if (ref->mtype == BSPLINE)
          if (EG_getRange(ref, lims, &n) == EGADS_SUCCESS)
            if (n != 0)
              if ((data[0] < lims[0]) || (data[1] > lims[1]) ||
                  (data[2] < lims[2]) || (data[3] > lims[3])) {
                printf(" EGADS Problem: Face #%d BSPLINE limits out of Period!\n",
                       i+1);
                printf("                %lf %lf %lf %lf   %lf %lf %lf %lf\n",
                       data[0], data[1], data[2], data[3],
                       lims[0], lims[1], lims[2], lims[3]);
              }
      n = Fwrite(&mtype,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      n = Fwrite(&nchild, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      iref = 0;
      if (ref != NULL) iref = EG_lookAtMap(ref, SURFACE, &maps, 0);
      if (iref < EGADS_SUCCESS) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return iref;
      }
      n = Fwrite(&iref,   sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      n = Fwrite(data, sizeof(double), 4, fp);
      if (n != 4) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      n = Fwrite(bbox, sizeof(double), 6, fp);
      if (n != 6) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      n = Fwrite(&tol,  sizeof(double), 1, fp);
      if (n != 1) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      n = Fwrite(senses, sizeof(int), nchild, fp);
      if (n != nchild) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
      /* loops */
      for (m = 0; m < nchild; m++) {
        iref = EG_indexBodyTopo(bobject, children[m]);
        if (iref <= EGADS_SUCCESS) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return iref;
        }
        n = Fwrite(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
      }

      stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
      if (stat != EGADS_SUCCESS) {
        EG_free(objs);
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return stat;
      }
    }
    EG_free(objs);

    if (bobject->mtype != FACEBODY) {
      /* shells */
      EG_getBodyTopos(bobject, NULL, SHELL, &no, &objs);
      if (outLevel > 0)
        printf(" Writing %d Shells...\n", no);
      for (i = 0; i < no; i++) {
/*@-nullderef@*/
        obj  = objs[i];
/*@+nullderef@*/
        stat = EG_getTopology(obj, &ref, &oclass, &mtype, data, &nchild,
                              &children, &senses);
        if (stat == EGADS_SUCCESS) stat = EG_getBoundingBox(obj, bbox);
        if (stat != EGADS_SUCCESS) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return stat;
        }
        n = Fwrite(&mtype,  sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
        n = Fwrite(&nchild, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
        n = Fwrite(bbox, sizeof(double), 6, fp);
        if (n != 6) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
        /* faces */
        for (m = 0; m < nchild; m++) {
          iref = EG_indexBodyTopo(bobject, children[m]);
          if (iref <= EGADS_SUCCESS) {
            EG_free(objs);
            EG_free(maps.pcurves.objs);
            EG_free(maps.curves.objs);
            EG_free(maps.surfaces.objs);
            return iref;
          }
          n = Fwrite(&iref, sizeof(int), 1, fp);
          if (n != 1) {
            EG_free(objs);
            EG_free(maps.pcurves.objs);
            EG_free(maps.curves.objs);
            EG_free(maps.surfaces.objs);
            return EGADS_WRITERR;
          }
        }
        
        stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
        if (stat != EGADS_SUCCESS) {
          EG_free(objs);
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return stat;
        }
      }
      EG_free(objs);
    }
    
  }

  /* finish off the body */
  
  stat = EG_getTopology(bobject, &ref, &oclass, &mtype, data, &nchild,
                        &children, &senses);
  if (stat == EGADS_SUCCESS) stat = EG_getBoundingBox(bobject, bbox);
  if (stat != EGADS_SUCCESS) {
    EG_free(maps.pcurves.objs);
    EG_free(maps.curves.objs);
    EG_free(maps.surfaces.objs);
    return stat;
  }
  if (ntypes[7] != 0) {
    if (senses == NULL) {
      n = 1;
      for (i = 0; i < nchild; i++) {
        if (Fwrite(&n, sizeof(int), 1, fp) != 1) {
          EG_free(maps.pcurves.objs);
          EG_free(maps.curves.objs);
          EG_free(maps.surfaces.objs);
          return EGADS_WRITERR;
        }
      }
    } else {
      n = Fwrite(senses, sizeof(int), nchild, fp);
      if (n != nchild) {
        EG_free(maps.pcurves.objs);
        EG_free(maps.curves.objs);
        EG_free(maps.surfaces.objs);
        return EGADS_WRITERR;
      }
    }
  }
  n = Fwrite(bbox, sizeof(double), 6, fp);
  if (n != 6) {
    EG_free(maps.pcurves.objs);
    EG_free(maps.curves.objs);
    EG_free(maps.surfaces.objs);
    return EGADS_WRITERR;
  }
  stat = EG_writeAttrs(fp, (egAttrs *) bobject->attrs);

  EG_free(maps.pcurves.objs);
  EG_free(maps.curves.objs);
  EG_free(maps.surfaces.objs);

  return stat;
}


static int
EG_exportTess(egObject *tess, stream_T *fp)
{
  int          status, len;
  int          ntri, iedge, iface, nedge, nface;
  const double *pxyz  = NULL, *puv    = NULL, *pt    = NULL;
  const int    *ptype = NULL, *pindex = NULL, *ptris = NULL, *ptric = NULL;
  egTessel     *btess;
  egObject     *body;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;

  /* get the body from tessellation */
  btess = (egTessel *) tess->blind;
  body  = btess->src;

  /* get the sizes */
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE,  &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE,  &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }

  /* write the number of edges and faces */
  if (Fwrite(&nedge, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  if (Fwrite(&nface, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;

  /* write out the edge tessellation */
  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess, iedge+1, &len, &pxyz, &pt);
    if (status != EGADS_SUCCESS) return status;
    if (Fwrite(&len, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    if (len == 0) continue;
/*@-nullpass@*/
    if (Fwrite((void *) pxyz, sizeof(double), 3*len, fp) != 3*len)
      return EGADS_WRITERR;
    if (Fwrite((void *) pt,   sizeof(double),   len, fp) !=   len)
      return EGADS_WRITERR;
/*@+nullpass@*/
  }

  /* write out face tessellations */
  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess, iface+1, &len, &pxyz, &puv, &ptype, &pindex,
                            &ntri, &ptris, &ptric);
    if (status != EGADS_SUCCESS) return status;
    if (Fwrite(&len,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    if (Fwrite(&ntri, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    if ((len == 0) || (ntri == 0)) continue;
    ntri *= 3;
/*@-nullpass@*/
    if (Fwrite((void *) pxyz,  sizeof(double), 3*len, fp) != 3*len)
      return EGADS_WRITERR;
    if (Fwrite((void *) puv,   sizeof(double), 2*len, fp) != 2*len)
      return EGADS_WRITERR;
    if (Fwrite((void *) ptris, sizeof(int),     ntri, fp) !=  ntri)
      return EGADS_WRITERR;
/*@+nullpass@*/
  }

  /* write out the tessellation attributes */
  return EG_writeAttrs(fp, (egAttrs *) tess->attrs);
}


static int
EG_uvmapExport(void *uvmap, int *trmap, stream_T *fp)
{
  int          i, ival, trmp = 0, msrch = 0;
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  if (uvstruct->mdef  != 1) return EGADS_UVMAP;
  if (uvstruct->ndef  != 1) return EGADS_UVMAP;
  if (uvstruct->msrch != NULL) msrch = 1;
  if (trmap           != NULL) trmp  = 1;
  ival = (int) uvstruct->isrch;
  if (Fwrite(&ival,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  ival = (int) uvstruct->ibface;
  if (Fwrite(&ival,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  ival = (int) uvstruct->nbface;
  if (Fwrite(&ival,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  ival = (int) uvstruct->nnode;
  if (Fwrite(&ival,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  if (Fwrite(&msrch, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  if (Fwrite(&trmp,  sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
  
  if (Fwrite(&uvstruct->idibf[1], sizeof(int), uvstruct->nbface, fp) !=
      uvstruct->nbface) return EGADS_WRITERR;
  
  for (i = 1; i <= uvstruct->nbface; i++) {
    if (Fwrite(uvstruct->inibf[i],  sizeof(int), 3, fp) != 3)
      return EGADS_WRITERR;
    if (Fwrite(uvstruct->ibfibf[i], sizeof(int), 3, fp) != 3)
      return EGADS_WRITERR;
  }
        
  for (i = 1; i <= uvstruct->nnode; i++)
    if (Fwrite(uvstruct->u[i], sizeof(double), 2, fp) != 2)
      return EGADS_WRITERR;
  
  if (uvstruct->msrch != NULL)
    if (Fwrite(&uvstruct->msrch[1], sizeof(int), uvstruct->nbface, fp) !=
        uvstruct->nbface) return EGADS_WRITERR;

  if (trmap != NULL)
    if (Fwrite(trmap, sizeof(int), uvstruct->nbface, fp) != uvstruct->nbface)
      return EGADS_WRITERR;
    
  return EGADS_SUCCESS;
}


static int
EG_exportEBody(egObject *EBody, stream_T *fp)
{
  int      i, j, k, stat, nds[2];
  double   real = 0.0;
  egObject *body, *obj;
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  
  if (EBody == NULL)               return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)      return EGADS_NOTTESS;
  
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_exportEBody)!\n");
    return EGADS_NOTFOUND;
  }
  if (ebody->done == 0) {
    printf(" EGADS Error: EBody not finialized (EG_exportEBody)!\n");
    return EGADS_EFFCTOBJ;
  }
  body = ebody->ref;
  if (Fwrite(&EBody->mtype,         sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  stat = EG_writeAttrs(fp, (egAttrs *) EBody->attrs);
  if (stat != EGADS_SUCCESS) return stat;
  
  if (Fwrite(&ebody->eedges.nobjs,  sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  if (Fwrite(&ebody->eloops.nobjs,  sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  if (Fwrite(&ebody->efaces.nobjs,  sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  if (Fwrite(&ebody->eshells.nobjs, sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  if (Fwrite(&ebody->nedge,         sizeof(int),    1, fp) != 1)
    return EGADS_WRITERR;
  if (Fwrite(&ebody->angle,         sizeof(double), 1, fp) != 1)
    return EGADS_WRITERR;

  if (body->mtype == SOLIDBODY)
    if (Fwrite(ebody->senses, sizeof(int), ebody->eshells.nobjs, fp) !=
        ebody->eshells.nobjs) return EGADS_WRITERR;
  
  /* source Edges */
  for (j = 0; j < ebody->nedge; j++) {
    nds[0] = EG_indexBodyTopo(body, ebody->edges[j].edge);
    if (nds[0] <= EGADS_SUCCESS) {
      printf(" EGADS Error: Source Edge = %d (EG_exportEBody)!\n", nds[0]);
      return EGADS_TOPOERR;
    }
    if (Fwrite(&nds[0],                sizeof(int),    1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&ebody->edges[j].curve, sizeof(int),    1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&ebody->edges[j].npts,  sizeof(int),    1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(ebody->edges[j].dstart, sizeof(double), 3, fp) != 3)
      return EGADS_WRITERR;
    if (Fwrite(ebody->edges[j].dend,   sizeof(double), 3, fp) != 3)
      return EGADS_WRITERR;
    if (Fwrite(ebody->edges[j].ts,     sizeof(double), ebody->edges[j].npts,
               fp) != ebody->edges[j].npts) return EGADS_WRITERR;
  }
  
  /* EEdges */
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    obj    = ebody->eedges.objs[i];
    eedge  = (egEEdge *) obj->blind;
    nds[0] = EG_indexBodyTopo(body, eedge->nodes[0]);
    nds[1] = EG_indexBodyTopo(body, eedge->nodes[1]);
    if ((nds[0] <= EGADS_SUCCESS) || (nds[1] <= EGADS_SUCCESS)) {
      printf(" EGADS Error: Node indices = %d %d (EG_exportEBody)!\n",
             nds[0], nds[1]);
      return EGADS_TOPOERR;
    }
    if (Fwrite(&obj->mtype,   sizeof(short),  1, fp) != 1) return EGADS_WRITERR;
    if (Fwrite(&eedge->nsegs, sizeof(int),    1, fp) != 1) return EGADS_WRITERR;
    if (Fwrite(nds,           sizeof(int),    2, fp) != 2) return EGADS_WRITERR;
    if (Fwrite(eedge->trange, sizeof(double), 2, fp) != 2) return EGADS_WRITERR;
    for (j = 0; j < eedge->nsegs; j++) {
      nds[0] = 0;
      if (eedge->segs[j].nstart != NULL)
        nds[0] = EG_indexBodyTopo(body, eedge->segs[j].nstart);
      if (nds[0] < EGADS_SUCCESS) {
        printf(" EGADS Error: Node start = %d (EG_exportEBody)!\n", nds[0]);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&eedge->segs[j].iedge,  sizeof(int),    1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eedge->segs[j].sense,  sizeof(int),    1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&nds[0],                sizeof(int),    1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eedge->segs[j].tstart, sizeof(double), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eedge->segs[j].tend,   sizeof(double), 1, fp) != 1)
        return EGADS_WRITERR;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) return stat;
  }
  
  /* ELoops */
  for (i = 0; i < ebody->eloops.nobjs; i++) {
    obj   = ebody->eloops.objs[i];
    eloop = (egELoop *) obj->blind;
    if (Fwrite(&obj->mtype, sizeof(short), 1, fp) != 1) return EGADS_WRITERR;
    if (eloop == NULL) {
      nds[0] = nds[1] = 0;
      if (Fwrite(nds,   sizeof(int),    2, fp) != 2) return EGADS_WRITERR;
      if (Fwrite(&real, sizeof(double), 1, fp) != 1) return EGADS_WRITERR;
      continue;
    }
    if (Fwrite(&eloop->eedges.nobjs, sizeof(int),    1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eloop->nedge,        sizeof(int),    1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eloop->area,         sizeof(double), 1, fp) != 1)
      return EGADS_WRITERR;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: EEdge %d index = %d in ELoop %d (EG_exportEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&k, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    }
    if (Fwrite(eloop->senses, sizeof(int), eloop->eedges.nobjs, fp) !=
        eloop->eedges.nobjs) return EGADS_WRITERR;
    for (j = 0; j < eloop->nedge; j++) {
      k = EG_indexBodyTopo(body, eloop->edgeUVs[j].edge);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Edge %d index = %d in ELoop %d (EG_exportEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&k,                       sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eloop->edgeUVs[j].sense, sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eloop->edgeUVs[j].npts,  sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(eloop->edgeUVs[j].iuv,    sizeof(int), eloop->edgeUVs[j].npts,
                 fp) != eloop->edgeUVs[j].npts) return EGADS_WRITERR;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) return stat;
  }
  
  /* EFaces */
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    obj   = ebody->efaces.objs[i];
    eface = (egEFace *) obj->blind;
    if (Fwrite(&obj->mtype,          sizeof(short), 1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eface->npatch,       sizeof(int),   1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eface->eloops.nobjs, sizeof(int),   1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eface->last,         sizeof(int),   1, fp) != 1)
      return EGADS_WRITERR;
    if (eface->npatch != 1) {
      stat = EG_uvmapExport(eface->uvmap, eface->trmap, fp);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d  uvmapExport = %d (EG_exportEBody)!\n",
               i+1, stat);
        return stat;
      }
    } else {
      if (Fwrite(eface->range, sizeof(double), 4, fp) != 4)
        return EGADS_WRITERR;
    }
    for (j = 0; j < eface->eloops.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eface->eloops.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Loop %d index = %d in EFace %d (EG_exportEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&k, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    }
    if (Fwrite(eface->senses, sizeof(int), eface->eloops.nobjs, fp) !=
        eface->eloops.nobjs) return EGADS_WRITERR;
    for (j = 0; j < eface->npatch; j++) {
      k = EG_indexBodyTopo(body, eface->patches[j].face);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d index = %d in EFace %d (EG_exportEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&k,                          sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eface->patches[j].start,    sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eface->patches[j].nuvs,     sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eface->patches[j].ndeflect, sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(&eface->patches[j].ntris,    sizeof(int), 1, fp) != 1)
        return EGADS_WRITERR;
      if (Fwrite(eface->patches[j].uvtris,    sizeof(int),
                 3*eface->patches[j].ntris, fp) != 3*eface->patches[j].ntris)
        return EGADS_WRITERR;
      if (Fwrite(eface->patches[j].uvs,       sizeof(double),
                 2*eface->patches[j].nuvs, fp) != 2*eface->patches[j].nuvs)
        return EGADS_WRITERR;
      if (Fwrite(eface->patches[j].deflect,   sizeof(double),
                 3*eface->patches[j].ndeflect, fp) !=
          3*eface->patches[j].ndeflect) return EGADS_WRITERR;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) return stat;
  }
  
  /* EShells */
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    obj    = ebody->eshells.objs[i];
    eshell = (egEShell *) obj->blind;
    if (Fwrite(&obj->mtype,           sizeof(short), 1, fp) != 1)
      return EGADS_WRITERR;
    if (Fwrite(&eshell->efaces.nobjs, sizeof(int),   1, fp) != 1)
      return EGADS_WRITERR;
    for (j = 0; j < eshell->efaces.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eshell->efaces.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d index = %d in EFace %d (EG_exportEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      if (Fwrite(&k, sizeof(int), 1, fp) != 1) return EGADS_WRITERR;
    }
    stat = EG_writeAttrs(fp, (egAttrs *) obj->attrs);
    if (stat != EGADS_SUCCESS) return stat;
  }
  
  return EGADS_SUCCESS;
}


int
EG_exportModel(ego mobject, size_t *nbytes, char **stream)
{
  int      i, j, n, oclass, mtype, nbody, *senses, rev[2] = {1, 1};
  double   bbox[6];
  egObject *ref, **bodies;
  egTessel *btess;
  egEBody  *ebody;
  stream_T myStream;
  stream_T *fp = &(myStream);

  /* default returns */
  *nbytes = 0;
  *stream = NULL;

  if (mobject == NULL)               return EGADS_NULLOBJ;
  if (mobject->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject->oclass != MODEL)      return EGADS_NOTMODEL;

  i = EG_getTopology(mobject, &ref, &oclass, &mtype, NULL, &nbody, &bodies,
                     &senses);
  if (i != EGADS_SUCCESS) return i;
  i = EG_getBoundingBox(mobject, bbox);
  if (i != EGADS_SUCCESS) return i;
  
  fp->size = CHUNK;
  fp->ptr  = 0;
  fp->data = EG_alloc(fp->size);

  /* put header */
  i = MAGIC;
  n = Fwrite(&i,        sizeof(int),    1, fp);
  if (n != 1) {
    Fclose(fp);
    return EGADS_WRITERR;
  }
  n = Fwrite(rev,       sizeof(int),    2, fp);
  if (n != 2) {
    Fclose(fp);
    return EGADS_WRITERR;
  }

  n = Fwrite(bbox,      sizeof(double), 6, fp);
  if (n != 6) {
    Fclose(fp);
    return EGADS_WRITERR;
  }
  n = Fwrite(&nbody,    sizeof(int),    1, fp);
  if (n != 1) {
    Fclose(fp);
    return EGADS_WRITERR;
  }
  i = EG_writeAttrs(fp, (egAttrs *) mobject->attrs);
  if (i != EGADS_SUCCESS) {
    Fclose(fp);
    return i;
  }

  /* write all of the bodies */
  for (n = 0; n < nbody; n++) {
    i = EG_exportBody(bodies[n], fp);
    if (i == EGADS_SUCCESS) continue;
    /* errorred out -- cleanup */
    Fclose(fp);
    return i;
  }
  
  /* write possible tessellation and EBody Objects */
  n = Fwrite(&mtype,    sizeof(int),    1, fp);
  if (n != 1) {
    Fclose(fp);
    return EGADS_WRITERR;
  }
  for (i = nbody; i < mtype; i++) {
    oclass = bodies[i]->oclass;
    if (oclass == TESSELLATION) {
      btess = (egTessel *) bodies[i]->blind;
      ref   = btess->src;
    } else if (oclass == EBODY) {
      ebody = (egEBody *)  bodies[i]->blind;
      ref   = ebody->ref;
    } else {
      printf(" Export Error: %d Entry in Model has class = %d!\n", i+1, oclass);
      Fclose(fp);
      return EGADS_NOTBODY;
    }
    j = Fwrite(&oclass, sizeof(int),    1, fp);
    if (j != 1) {
      Fclose(fp);
      return EGADS_WRITERR;
    }
    for (n = 0; n < i; n++)
      if (ref == bodies[n]) break;
    if (n == nbody) {
      printf(" Export Error: %d Entry in Model cannot find Body!\n", i+1);
      Fclose(fp);
      return EGADS_NOTBODY;
    }
    n++;
    j = Fwrite(&n,      sizeof(int),    1, fp);
    if (j != 1) {
      Fclose(fp);
      return EGADS_WRITERR;
    }
      
    if (bodies[i]->oclass == TESSELLATION) {
      n = EG_exportTess(bodies[i], fp);
      if (n != EGADS_SUCCESS) {
        /* errorred out -- cleanup */
        Fclose(fp);
        return n;
      }
    } else if (bodies[i]->oclass == EBODY) {
      n = EG_exportEBody(bodies[i], fp);
      if (n != EGADS_SUCCESS) {
        /* errorred out -- cleanup */
        Fclose(fp);
        return n;
      }
    }
  }

  /* return results */
  *nbytes = fp->ptr;
  *stream = fp->data;

  return EGADS_SUCCESS;
}


#ifdef STANDALONE
int main(int argc, char *argv[])
{
  size_t nbytes;
  ego    context, model;
  char   *stream;
  FILE   *fp;
  
  if (argc != 3) {
    printf(" Usage: writeLite modelFile liteFile\n\n");
    exit(EXIT_FAILURE);
  }
  /* initialize */
  printf(" EG_open          = %d\n", EG_open(&context));
  printf(" EG_loadModel     = %d  %s\n", EG_loadModel(context, 0, argv[1],
                                                      &model), argv[1]);

  printf(" EG_exportModel   = %d\n", EG_exportModel(model, &nbytes, &stream));
  
  fp = fopen(argv[2], "wb");
  if (fp == NULL) exit(EXIT_FAILURE);
  fwrite(stream, sizeof(char), nbytes, fp);
  fclose(fp);
  EG_free(stream);

  printf(" EG_deleteObject  = %d\n", EG_deleteObject(model));
  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}
#endif
