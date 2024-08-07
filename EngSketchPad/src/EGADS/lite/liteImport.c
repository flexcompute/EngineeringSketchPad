/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Import a Model from EGADS (via a string)
 *
 *      Copyright 2011-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
/*@-redef@*/
typedef int    INT_;
typedef INT_   INT_3D[3];
typedef double DOUBLE_2D[2];
/*@+redef@*/
#include "uvmap_struct.h"
#include "liteClasses.h"
#include "liteDevice.h"


/* #define DEBUG */
#define FULLATTR


  extern int EG_close(egObject *context);
  extern int EG_initTessBody(egObject *object, egObject **tess);
  extern int EG_setTessEdge(egObject *tess, int eIndex, int len,
                            const double *xyz, const double *t);
  extern int EG_setTessFace(egObject *tess, int fIndex, int len,
                            const double *xyz, const double *uv, int ntri,
                            const int *tris);
  extern int EG_statusTessBody(egObject *tess, egObject **body, int *state,
                               int *npts);
  extern int EG_objectBodyTopo(const egObject *body, int oclass, int index,
                               egObject **obj);
  extern int EG_effectNeighbor(egEFace *eface);
#ifdef __NVCC__
  extern int EG_evaluateDev(const egObject *geom_d, const double *param,
                            double *ev);
  extern int EG_attributeRetDev(const egObject *obj_d, const char *name,
                                int *atype, int *len, int **ints,
                                double **reals, char **str);
  extern int EG_addStrAttrDev(egObject *obj_d, const char *name,
                              const char *str);
#else
  extern int EG_evaluate(const egObject *geom, const double *param, double *ev);
  extern int EG_attributeRet(const egObject *obj, const char *name, int *atype,
                             int *len, /*@null@*/ const int **ints,
                                       /*@null@*/ const double **reals,
                                       /*@null@*/ const char **str);
#endif


typedef struct {
  void   *data;
  size_t ptr;
  size_t size;
  int    swap;
} stream_T;



static void
swap(void *buffer, size_t size)
{
  char *buf, save;

  buf = (char *) buffer;
  if (size == 2) {
    save   = buf[1];
    buf[1] = buf[0];
    buf[0] = save;
  } else if (size == 4) {
    save   = buf[3];
    buf[3] = buf[0];
    buf[0] = save;
    save   = buf[2];
    buf[2] = buf[1];
    buf[1] = save;
  } else {
    save   = buf[7];
    buf[7] = buf[0];
    buf[0] = save;
    save   = buf[6];
    buf[6] = buf[1];
    buf[1] = save;
    save   = buf[5];
    buf[5] = buf[2];
    buf[2] = save;
    save   = buf[4];
    buf[4] = buf[3];
    buf[3] = save;
  }
}


static int
Fread(void *data, size_t size, int nitems, stream_T *stream)
{
  int  i;
  char *buf;

  buf = (char *) data;
  memcpy(data, &(((char *) stream->data)[stream->ptr]), size*nitems);
  if ((size != sizeof(char)) && (stream->swap == 1))
    for (i = 0; i < nitems; i++) swap(&buf[i*size], size);
  stream->ptr += size*nitems;

  return nitems;
}


#ifdef FULLATTR
static void
EG_attrBuildSeq(egAttrs *attrs)
{
  int       i, j, l, n, snum, *hit, nospace = 0, nseqs = 0;
  char      *root, *newname;
  egAttr    *attr;
  egAttrSeq *seqs = NULL, *tmp;
  egAttrs   attrs_, *attrs_h = &attrs_;

  EG_GET_ATTRS(attrs_h, attrs);

  hit = (int *) EG_alloc(attrs_h->nattrs*sizeof(int));
  if (hit == NULL) {
    printf(" EGADS Internal: Malloc on %d attributes!\n", attrs_h->nattrs);
    return;
  }
  for (i = 0; i < attrs_h->nattrs; i++) hit[i] = 0;

  /* build the sequence structure */
  for (i = 0; i < attrs_h->nattrs; i++) {
    if (hit[i] != 0) continue;
    if (attrs->attrs[i].name == NULL) continue;
    l = strlen(attrs->attrs[i].name);
    for (n = j = 0; j < l; j++)
      if (attrs->attrs[i].name[j] == 32) n++;
    if (n == 0) {
      /* only the first can have no space */
      nospace = 1;
    } else if (n >  1) {
      printf(" EGADS Internal: More than a single space (%d) in an Attr name!\n",
             n);
      continue;
    }
    /* make the root name */
    n    = l;
    root = attrs->attrs[i].name;
    if (nospace == 0) {
      root = EG_strdup(attrs->attrs[i].name);
      if (root == NULL) {
        printf(" EGADS Internal: Null root on %s!\n", attrs->attrs[i].name);
        EG_free(hit);
        return;
      }
      for (n = 0; n < l; n++)
        if (root[n] == 32) {
          root[n] = 0;
          break;
        }
    }

    /* count the members */
    snum = hit[i] = 1;
    for (j = i+1; j < attrs->nattrs; j++) {
      if (hit[j] != 0) continue;
      if (attrs->attrs[j].name == NULL) continue;
      if (strlen(attrs->attrs[j].name) < n) continue;
      if (attrs->attrs[j].name[n] == 32)
        if (strncmp(attrs->attrs[i].name, attrs->attrs[j].name, n-1) == 0)
          snum++;
    }

    /* only a single member */
    if (snum == 1) {
      if (nospace == 1) continue;
      /* remove existing seq number */
      EG_free(attrs->attrs[i].name);
      attrs->attrs[i].name = root;
      continue;
    }

    /* build the sequence */
    if (nseqs == 0) {
      seqs = (egAttrSeq *) EG_alloc(sizeof(egAttrSeq));
      if (seqs == NULL) {
        EG_free(root);
        printf(" EGADS Internal: Malloc on Base Sequence!\n");
        continue;
      }
    } else {
      tmp = (egAttrSeq *) EG_reall(seqs, (nseqs+1)*sizeof(egAttrSeq));
      if (tmp == NULL) {
        EG_free(root);
        printf(" EGADS Internal: Malloc on %d Sequence!\n", nseqs+1);
        continue;
      }
      seqs = tmp;
    }
    seqs[nseqs].attrSeq = (int *) EG_alloc(snum*sizeof(int));
    if (seqs[nseqs].attrSeq == NULL) {
      EG_free(root);
      printf(" EGADS Internal: Malloc on %d Attr Seq Pointers!\n", snum);
      continue;
    }
    if (nospace == 1) root = EG_strdup(attrs->attrs[i].name);
    seqs[nseqs].nSeq = snum;
    seqs[nseqs].root = root;

    /* load the sequence */
    seqs[nseqs].attrSeq[0] = i;
    snum = 1;
    for (j = i+1; j < attrs->nattrs; j++) {
      if (hit[j] != 0) continue;
      if (attrs->attrs[j].name == NULL) continue;
      if (strlen(attrs->attrs[j].name) < n) continue;
      if (attrs->attrs[j].name[n] == 32)
        if (strncmp(attrs->attrs[i].name, attrs->attrs[j].name, n-1) == 0) {
          seqs[nseqs].attrSeq[snum] = j;
          snum++;
          hit[j] = 1;
        }
    }

    /* check/correct the sequence numbers */
    for (j = 0; j < snum; j++) {
      attr = &attrs->attrs[seqs[nseqs].attrSeq[j]];
      if ((j != 0) || (nospace == 0)) {
        l    = 0;
        sscanf(&attr->name[n], "%d", &l);
#ifdef DEBUG
        if (l == j+1) printf(" seq = %d, oldname = %s\n", j+1, attr->name);
#endif
        if (l == j+1) continue;
      }
      newname = (char *) EG_alloc((n+8)*sizeof(char));
      if (newname == NULL) {
        printf(" EGADS Internal: Malloc on name %s!\n", root);
        continue;
      }
      snprintf(newname, n+8, "%s %d", root, j+1);
      EG_free(attr->name);
      attr->name = newname;
#ifdef DEBUG
      printf(" seq = %d, newname = %s\n", j+1, newname);
#endif
    }
#ifdef DEBUG
    printf(" seq %d: root = %s,  snum = %d\n",
           nseqs, seqs[nseqs].root, seqs[nseqs].nSeq);
    for (j = 0; j < seqs[nseqs].nSeq; j++) {
      attr = &attrs->attrs[seqs[nseqs].attrSeq[j]];
      printf(" %d: %d  %s\n", j+1, seqs[nseqs].attrSeq[j], attr->name);
    }
#endif
    nseqs++;
  }
  EG_free(hit);

  attrs_h->nseqs = nseqs;
  attrs_h->seqs  = seqs;
  EG_SET_ATTRS(attrs, attrs_h);
}
#endif


static int
EG_addStrAttr(egObject *obj, const char *name, const char *str)
{
  int      i, length, find = -1;
  char     *temp;
  egAttr   *attr;
  egAttrs  *attrs;
  egObject obj_,   *obj_h   = &obj_;
  egAttrs  attrs_, *attrs_h = &attrs_;
  egAttr   attr_,  *attr_h  = &attr_;

  if (obj == NULL)               return EGADS_NULLOBJ;
  EG_GET_OBJECT(obj_h, obj);
  if (obj_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj_h->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj_h->oclass == NIL)        return EGADS_EMPTY;
  if (obj_h->oclass == REFERENCE)  return EGADS_REFERCE;

  if ((name == NULL) || (str == NULL)) {
      printf(" EGADS Internal: NULL Name/Value (EG_addStrAttr)!\n");
    return EGADS_NONAME;
  }
  length = strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] <= ' ') {
      length = 0;
      break;
    }
  if (length == 0) {
    printf(" EGADS Internal: BAD Name (EG_addStrAttr)!\n");
    return EGADS_INDEXERR;
  }
  attrs = (egAttrs *) obj_h->attrs;

  if (attrs != NULL) {
    EG_NEW(&temp, char, length+2);
    if (temp != NULL) {
      EG_GET_ATTRS(attrs_h, attrs);
      for (i = 0; i < attrs_h->nattrs; i++) {
        EG_GET_STRN(temp, &(attrs_h->attrs[i].name), length+1);
        if (strcmp(temp,name) == 0) {
          find = i;
          break;
        }
      }
      EG_FREE(temp);
    }
  }

  if ((find != -1) && (attrs != NULL)) {

    /* an existing attribute -- reset the values */

    EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));
    if (attr_h->type == ATTRINT) {
      if (attr_h->length != 1)
        EG_FREE(attr_h->vals.integers);
    } else if (attr_h->type == ATTRREAL) {
      if (attr_h->length != 1)
        EG_FREE(attr_h->vals.reals);
    } else if (attr_h->type == ATTRCSYS) {
      EG_FREE(attr_h->vals.reals);
    } else if (attr_h->type == ATTRSTRING) {
      EG_FREE(attr_h->vals.string);
    }

  } else {

    if (attrs == NULL) {
      EG_NEW(&attrs, egAttrs, 1);
      if (attrs == NULL) {
        printf(" EGADS Internal: Attrs MALLOC for %s (EG_addStrAttr)!\n",
               name);
        return EGADS_MALLOC;
      }
      attrs_h->nattrs = 0;
      attrs_h->attrs  = NULL;
      attrs_h->nseqs  = 0;
      attrs_h->seqs   = NULL;
/*@-nullret@*/
      EG_SET_ATTRS(attrs, attrs_h);
/*@+nullret@*/
      EG_SET_OBJECT_PTR(&(obj->attrs), &attrs);
    }
    EG_GET_ATTRS(attrs_h, attrs);
    if (attrs_h->attrs == NULL) {
      EG_NEW(&attr, egAttr, (attrs_h->nattrs+1));
    } else {
      EG_REALLOC(&attr, attrs_h->attrs, egAttr, attrs_h->nattrs,
                 (attrs_h->nattrs+1));
    }
    if (attr == NULL) {
      printf(" EGADS Internal: Attr MALLOC for %s (EG_addStrAttr)!\n",
             name);
      return EGADS_MALLOC;
    }
    attrs_h->attrs = attr;
    find = attrs_h->nattrs;
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));
    attr_h->vals.string = NULL;
    EG_SET_STR(&(attr_h->name), name);
    if (attr_h->name == NULL) return EGADS_MALLOC;
    EG_SET_ATTR(&(attrs_h->attrs[find]), attr_h);
    attrs_h->nattrs += 1;
  }

  EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));

  attr_h->type        = ATTRSTRING;
  attr_h->length      = 0;
  EG_SET_STR(&(attr_h->vals.string), str);
  if (attr_h->vals.string != NULL) {
    attr_h->length = strlen(str);
  }
  EG_SET_ATTR(&(attrs_h->attrs[find]), attr_h);

  EG_SET_ATTRS(attrs, attrs_h);

  return EGADS_SUCCESS;
}


static void
EG_freeAttrs(egAttrs **attrx)
{
  int     i;
  egAttrs *attrs;
  egAttrs attrs_, *attrs_h = &attrs_;
  egAttr  attr_,  *attr_h  = &attr_;
  egAttrs *nil = NULL;

  attrs = *attrx;
  if (attrs == NULL) return;
  EG_GET_ATTRS(attrs_h, *attrx);
  EG_SET_ATTRS_PTR(attrx, nil);

  /* remove any attributes */
  for (i = 0; i < attrs_h->nseqs; i++) {
    EG_FREE(attrs_h->seqs[i].root);
    EG_FREE(attrs_h->seqs[i].attrSeq);
  }
  if (attrs_h->seqs != NULL) EG_FREE(attrs_h->seqs);
  for (i = 0; i < attrs_h->nattrs; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
    if (attr_.name != NULL) EG_FREE(attr_.name);
    if (attr_.type == ATTRINT) {
      if (attr_.length > 1) EG_FREE(attr_.vals.integers);
    } else if ((attr_.type == ATTRREAL) ||
               (attr_.type == ATTRCSYS)) {
      if (attr_.length > 1) EG_FREE(attr_.vals.reals);
    } else if (attr_.type == ATTRSTRING) {
      EG_FREE(attr_.vals.string);
    }
  }
  EG_FREE(attrs_h->attrs);
  EG_FREE(attrs);
}


static int
EG_readString(stream_T *fp, char **string)
{
  int    len;
  size_t n;
  char   *string_h;
  int    status = EGADS_SUCCESS;

  *string = string_h = NULL;
  n = Fread(&len, sizeof(int), 1, fp);
  if (n   != 1) return EGADS_READERR;
  if (len <  0) return EGADS_READERR;
  if (len == 0) return EGADS_SUCCESS;

  string_h = (char *) EG_alloc(len*sizeof(char));
  if (string_h == NULL) return EGADS_MALLOC;

  n = Fread(string_h, sizeof(char), len, fp);
  if (n != len) {
    EG_free(string_h);
    *string = NULL;
    return EGADS_READERR;
  }
  EG_SET_STR(&(string[0]), string_h);

  EG_free(string_h);

  return status;
}


static int
EG_readAttrs(stream_T *fp, egAttrs **attrx)
{
  int     nattr, i, status;
#ifdef FULLATTR
  int     len, j, nspace = 0;
#endif
  size_t  n;
  egAttr  *attr;
  egAttrs *attrs;
  egAttrs attrs_, *attrs_h = &attrs_;
  egAttr  attr_, *attr_h = &attr_;
  void    *temp;

  *attrx = NULL;
  n = Fread(&nattr, sizeof(int), 1, fp);
  if (n     != 1) return EGADS_READERR;
  if (nattr == 0) return EGADS_SUCCESS;

  EG_NEW(&attrs, egAttrs, 1);
  if (attrs == NULL) return EGADS_MALLOC;
  EG_NEW(&attr, egAttr, nattr);
  if (attr == NULL) {
    EG_FREE(attrs);
    return EGADS_MALLOC;
  }
  attrs_h->nattrs = nattr;
  attrs_h->attrs  = attr;
  attrs_h->nseqs  = 0;
  attrs_h->seqs   = NULL;
/*@-nullret@*/
  EG_SET_ATTRS(attrs, attrs_h);
/*@+nullret@*/
  for (i = 0; i < nattr; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
    attr_.name   = NULL;
    attr_.length = 1;
    attr_.type   = ATTRINT;
    EG_SET_ATTR(&(attrs_h->attrs[i]), attr_h);
  }

/*@-mustfreefresh@*/
  for (i = 0; i < nattr; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));

    n = Fread(&attr_.type,   sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    n = Fread(&attr_.length, sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    status = EG_readString(fp, &attr_.name);
    if (status != EGADS_SUCCESS) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
#ifdef FULLATTR
    if (attr_.name != NULL) {
      len = strlen(attr_.name);
      for (j = 0; j < len; j++)
        if (attr_.name[j] == 32) nspace++;
    }
#endif
    if (attr_.type == ATTRINT) {
      n = attr_.length;
      if (attr_.length == 1) {
        n = Fread(&attr_.vals.integer, sizeof(int), 1, fp);
      } else if (attr_.length > 1) {
        EG_NEW(&attr_.vals.integers, int, attr_.length);
        if (attr_.vals.integers == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        temp = EG_alloc(attr_.length*sizeof(int));
        if (temp == NULL) {
          EG_FREE(attr_.vals.integers);
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread((int *) temp, sizeof(int), attr_.length, fp);
        EG_COPY(attr_.vals.integers, temp, int, attr_.length);
        EG_free(temp);
      }
      if (n != attr_.length) {
        EG_FREE(attr_.vals.integers);
        EG_freeAttrs(&attrs);
        return EGADS_READERR;
      }
    } else if ((attr_.type == ATTRREAL) || (attr_.type == ATTRCSYS)) {
      n = attr_.length;
      if (attr_.length == 1) {
        n = Fread(&attr_.vals.real, sizeof(double), 1, fp);
      } else if (attr_.length > 1) {
        EG_NEW(&attr_.vals.reals, double, attr_.length);
        if (attr_.vals.reals == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        temp = EG_alloc(attr_.length*sizeof(double));
        if (temp == NULL) {
          EG_FREE(attr_.vals.reals);
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread((double *) temp, sizeof(double), attr_.length, fp);
        EG_COPY(attr_.vals.reals, temp, double, attr_.length);
        EG_free(temp);
      }
      if (n != attr_.length) {
        EG_FREE(attr_.vals.reals);
        EG_freeAttrs(&attrs);
        return EGADS_READERR;
      }
    } else {
      status = EG_readString(fp, &attr_.vals.string);
      if (status != EGADS_SUCCESS) return EGADS_READERR;
    }
    EG_SET_ATTR(&(attrs_h->attrs[i]), attr_h);
  }
/*@+mustfreefresh@*/

#ifdef FULLATTR
  /* sequences exist! */
  if (nspace != 0) EG_attrBuildSeq(attrs);
#endif

  *attrx = attrs;
  return EGADS_SUCCESS;
}


static int
EG_readGeometry(liteGeometry *lgeom, int *iref, stream_T *fp)
{
  int          n, nhead, ndata;
  liteGeometry lgeom_, *lgeom_h = &lgeom_;
  void         *temp;

  EG_GET_GEOM(lgeom_h, lgeom);

  *iref           = 0;
  lgeom_h->ref    = NULL;
  lgeom_h->header = NULL;
  lgeom_h->data   = NULL;
/*@-nullret@*/
  EG_SET_GEOM(lgeom, lgeom_h);
/*@+nullret@*/
  n = Fread(iref,   sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&nhead, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ndata, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;

  if (nhead != 0) {
    EG_NEW(&(lgeom_h->header), int, nhead);
    EG_COPY(&(lgeom->header), &(lgeom_h->header), int *, 1);
    if (lgeom_h->header == NULL) return EGADS_MALLOC;
    temp = EG_alloc(nhead*sizeof(int));
    if (temp == NULL) return EGADS_MALLOC;
    n = Fread((int *) temp, sizeof(int), nhead, fp);
    EG_COPY(lgeom_h->header, temp, int, nhead);
    EG_free(temp);
    if (n != nhead) return EGADS_READERR;
  }
  EG_NEW(&(lgeom_h->data), double, ndata);
  EG_COPY(&(lgeom->data), &(lgeom_h->data), double *, 1);
  if (lgeom_h->data == NULL) return EGADS_MALLOC;
  temp = EG_alloc(ndata*sizeof(double));
  if (temp == NULL) return EGADS_MALLOC;
  n = Fread((double *) temp, sizeof(double), ndata, fp);
  EG_COPY(lgeom_h->data, temp, double, ndata);
  EG_free(temp);
  if (n != ndata) return EGADS_READERR;

  return EGADS_SUCCESS;
}


static int
EG_readBody(egObject *context, egObject *mobject, int bindex, stream_T *fp)
{
  int          i, j, m, n, stat, mtype, iref, ntypes[8];
  double       t, d, x0[2], x1[2], data[6];
  int          atype, alen;
  const int    *ints;
  const double *reals;
  const char   *str;
  egObject     *obj, *bobj, *pcobj;
  liteGeometry *lgeom;
  liteNode     *lnode;
  liteEdge     *ledge;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;
  egObject     context_, *context_h = &context_;
  egObject     mobject_, *mobject_h = &mobject_;
  egObject     bobj_, *bobj_h = &bobj_;
  liteModel    lmodel_, *lmodel_h = &lmodel_;
  liteBody     lbody_, *lbody_h = &lbody_;
  void         *nil = NULL;

  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (mobject == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(mobject_h, mobject);
  if (mobject_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject_h->oclass != MODEL)      return EGADS_NOTMODEL;
  lmodel = (liteModel *) mobject_h->blind;
  EG_GET_MODEL(lmodel_h, lmodel);

  n = Fread(&mtype, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(ntypes, sizeof(int), 8, fp);
  if (n != 8) return EGADS_READERR;
  EG_NEW(&lbody, liteBody, 1);
  if (lbody == NULL) return EGADS_MALLOC;
  EG_GET_BODY(lbody_h, lbody);
  lbody_h->pcurves.objs   = NULL;
  lbody_h->curves.objs    = NULL;
  lbody_h->surfaces.objs  = NULL;
  lbody_h->nodes.objs     = NULL;
  lbody_h->edges.objs     = NULL;
  lbody_h->loops.objs     = NULL;
  lbody_h->faces.objs     = NULL;
  lbody_h->shells.objs    = NULL;
  lbody_h->senses         = NULL;
  lbody_h->pcurves.nobjs  = 0;
  lbody_h->curves.nobjs   = 0;
  lbody_h->surfaces.nobjs = 0;
  lbody_h->nodes.nobjs    = 0;
  lbody_h->edges.nobjs    = 0;
  lbody_h->loops.nobjs    = 0;
  lbody_h->faces.nobjs    = 0;
  lbody_h->shells.nobjs   = 0;

#ifdef DEBUG
  printf(" Reading Body #%d: %d %d %d %d %d %d %d %d\n", bindex+1, ntypes[0],
         ntypes[1], ntypes[2], ntypes[3], ntypes[4], ntypes[5], ntypes[6],
         ntypes[7]);
#endif

  stat = EG_makeObject(context, &bobj);
  if (stat != EGADS_SUCCESS) {
    EG_FREE(lbody);
    return stat;
  }
  EG_GET_OBJECT(bobj_h, bobj);
  bobj_h->oclass = BODY;
  bobj_h->mtype  = mtype;
  bobj_h->blind  = lbody;
  bobj_h->topObj = mobject;
  EG_SET_OBJECT(&bobj, bobj_h);
  EG_SET_OBJECT_PTR(&(lmodel_h->bodies[bindex]), &bobj);

  /* make all of the objects */
  if (ntypes[0] > 0) {
    EG_NEW(&(lbody_h->pcurves.objs), egObject *, ntypes[0]);
    if (lbody_h->pcurves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[0]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->pcurves.objs[i]), &nil);
    lbody_h->pcurves.nobjs = ntypes[0];
    for (i = 0; i < ntypes[0]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->pcurves.objs[i]), &obj);
    }
  }
  if (ntypes[1] > 0) {
    EG_NEW(&(lbody_h->curves.objs), egObject *, ntypes[1]);
    if (lbody_h->curves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[1]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->curves.objs[i]), &nil);
    lbody_h->curves.nobjs = ntypes[1];
    for (i = 0; i < ntypes[1]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->curves.objs[i]), &obj);
    }
  }
  if (ntypes[2] > 0) {
    EG_NEW(&(lbody_h->surfaces.objs), egObject *, ntypes[2]);
    if (lbody_h->surfaces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[2]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->surfaces.objs[i]), &nil);
    lbody_h->surfaces.nobjs = ntypes[2];
    for (i = 0; i < ntypes[2]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->surfaces.objs[i]), &obj);
    }
  }
  if (ntypes[3] > 0) {
    EG_NEW(&(lbody_h->nodes.objs), egObject *, ntypes[3]);
    if (lbody_h->nodes.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[3]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->nodes.objs[i]), &nil);
    lbody_h->nodes.nobjs = ntypes[3];
    for (i = 0; i < ntypes[3]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->nodes.objs[i]), &obj);
    }
  }
  if (ntypes[4] > 0) {
    EG_NEW(&(lbody_h->edges.objs), egObject *, ntypes[4]);
    if (lbody_h->edges.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[4]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->edges.objs[i]), &nil);
    lbody_h->edges.nobjs = ntypes[4];
    for (i = 0; i < ntypes[4]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->edges.objs[i]), &obj);
    }
  }
  if (ntypes[5] > 0) {
    EG_NEW(&(lbody_h->loops.objs), egObject *, ntypes[5]);
    if (lbody_h->loops.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[5]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->loops.objs[i]), &nil);
    lbody_h->loops.nobjs = ntypes[5];
    for (i = 0; i < ntypes[5]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->loops.objs[i]), &obj);
    }
  }
  if (ntypes[6] > 0) {
    EG_NEW(&(lbody_h->faces.objs), egObject *, ntypes[6]);
    if (lbody_h->faces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[6]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->faces.objs[i]), &nil);
    lbody_h->faces.nobjs = ntypes[6];
    for (i = 0; i < ntypes[6]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->faces.objs[i]), &obj);
    }
  }
  if (ntypes[7] > 0) {
    EG_NEW(&(lbody_h->shells.objs), egObject *, ntypes[7]);
    if (lbody_h->shells.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[7]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->shells.objs[i]), &nil);
    lbody_h->shells.nobjs = ntypes[7];
    for (i = 0; i < ntypes[7]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->shells.objs[i]), &obj);
    }
  }
/*@-nullret@*/
  EG_SET_BODY(lbody, lbody_h);
/*@+nullret@*/

#ifdef DEBUG
  printf(" Reading %d PCurves...\n", ntypes[0]);
#endif
  /* pcurves */
  if (lbody_h->pcurves.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    int header_h[4];
    double data_h[6];
    for (i = 0; i < lbody_h->pcurves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
      EG_GET_GEOM(lgeom_h, lgeom);
      if (iref != 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->pcurves.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->pcurves.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = PCURVE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
      if (mtype != BSPLINE) continue;
      EG_GET_HEADERS_AT(header_h, lgeom_h->header, 4);
      for (n = 2; n < header_h[2]; n++) {
        m     = header_h[3] + 2*n - 4;
        EG_GET_DATA_AT(data_h, &(lgeom_h->data[m]), 6);
        x0[0] = data_h[  2] - data_h[  0];
        x0[1] = data_h[  3] - data_h[  1];
        d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
        if (d != 0.0) {
          x0[0] /= d;
          x0[1] /= d;
        }
        x1[0] = data_h[  4] - data_h[  2];
        x1[1] = data_h[  5] - data_h[  3];
        d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
        if (d != 0.0) {
          x1[0] /= d;
          x1[1] /= d;
        }
        d = x0[0]*x1[0] + x0[1]*x1[1];
        if (d < -0.95) {
#ifdef DEBUG
          double last_data_h;
          EG_GET_DATA_AT(&last_data_h, &(lgeom_h->data[header_h[3]-1]), 1);
          printf(" EGADS Info: PCurve %d dot flip at %d/%d (%lf) -- %lf %lf!\n",
                 i, n-2, header_h[2]-2, d,
                 data_h[0], last_data_h);
#endif
          stat = EG_addStrAttr(obj, ".Bad", "CPrev");
          if (stat != EGADS_SUCCESS)
            printf("             EG_addStrAttr CPrev= %d\n", stat);
        }
      }
    }
  }

  /* curves */
#ifdef DEBUG
  printf(" Reading %d Curves...\n", ntypes[1]);
#endif
  if (lbody_h->curves.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->curves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
      if (iref != 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->curves.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->curves.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = CURVE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* surfaces */
#ifdef DEBUG
  printf(" Reading %d Surfaces...\n", ntypes[2]);
#endif
  if (lbody_h->surfaces.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->surfaces.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
/*@-nullderef@*/
      if (iref < 0) {
        if(lgeom->ref == NULL) {
          EG_GET_GEOM(lgeom_h, lgeom);
          if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
          if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
          EG_FREE(lgeom);
          return EGADS_READERR;
        }
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->curves.objs[-iref-1]));
      }
/*@+nullderef@*/
      if (iref > 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->surfaces.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->surfaces.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = SURFACE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* nodes */
#ifdef DEBUG
  printf(" Reading %d Nodes...\n", ntypes[3]);
#endif
  if (lbody_h->nodes.objs != NULL) {
    liteNode lnode_, *lnode_h = &lnode_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->nodes.nobjs; i++) {
      EG_NEW(&lnode, liteNode, 1);
      if (lnode == NULL) return EGADS_MALLOC;
      n = Fread(lnode_h->xyz,  sizeof(double), 3, fp);
      if (n != 3) {
        EG_FREE(lnode);
        return EGADS_READERR;
      }
      n = Fread(&lnode_h->tol, sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(lnode);
        return EGADS_READERR;
      }
      EG_SET_NODE(lnode, lnode_h);
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->nodes.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = NODE;
      obj_h->mtype  = 0;
      obj_h->blind  = lnode;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* edges */
#ifdef DEBUG
  printf(" Reading %d Edges...\n", ntypes[4]);
#endif
  if (lbody_h->edges.objs != NULL) {
    liteEdge ledge_, *ledge_h = &ledge_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->edges.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&ledge, liteEdge, 1);
      if (ledge == NULL) return EGADS_MALLOC;
      EG_GET_EDGE(ledge_h, ledge);
      ledge_h->curve    = NULL;
      ledge_h->nodes[0] = NULL;
      ledge_h->nodes[1] = NULL;

      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
#ifndef __clang_analyzer__
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->curve),
                                       &(lbody_h->curves.objs[iref-1]));
#endif
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
#ifndef __clang_analyzer__
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->nodes[0]),
                                       &(lbody_h->nodes.objs[iref-1]));
#endif
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
#ifndef __clang_analyzer__
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->nodes[1]),
                                       &(lbody_h->nodes.objs[iref-1]));
#endif
/*@+nullderef@*/

      n = Fread(ledge_h->trange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
      n = Fread(ledge_h->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
      n = Fread(&ledge_h->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullret@*/
      EG_SET_EDGE(ledge, ledge_h);
/*@+nullret@*/

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->edges.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = EDGE;
      obj_h->mtype  = mtype;
      obj_h->blind  = ledge;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* loops */
#ifdef DEBUG
  printf(" Reading %d Loops...\n", ntypes[5]);
#endif
  if (lbody_h->loops.objs != NULL) {
    liteLoop lloop_, *lloop_h = &lloop_;
    int      *itemp;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->loops.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lloop, liteLoop, 1);
      if (lloop == NULL) return EGADS_MALLOC;
      lloop_h->nedges  = m;
      lloop_h->surface = NULL;
      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(lloop);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(lloop_h->surface),
                          &(lbody_h->surfaces.objs[iref-1]));
#endif
/*@+nullderef@*/
        m *= 2;
      }
      n = Fread(lloop_h->bbox,  sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lloop);
        return EGADS_READERR;
      }
      if (lloop_h->nedges != 0) {
        egObject **otemp;
        EG_NEW(&(lloop_h->senses), int, lloop_h->nedges);
        if (lloop_h->senses == NULL) {
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }

        itemp = (int *) EG_alloc(lloop_h->nedges*sizeof(int));
        if (itemp == NULL) {
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        n = Fread((int *) itemp,  sizeof(int), lloop_h->nedges, fp);
        if (n != lloop_h->nedges) {
          EG_free(itemp);
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_READERR;
        }
        EG_COPY(lloop_h->senses, itemp, int, lloop_h->nedges);
        EG_free(itemp);

        EG_NEW(&(lloop_h->edges), egObject *, m);
        if (lloop_h->edges == NULL) {
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        otemp = (egObject **) EG_alloc(m*sizeof(egObject *));
        if (otemp == NULL) {
          EG_FREE(lloop_h->edges);
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        for (j = 0; j < m; j++) {
          n = Fread(&iref, sizeof(int), 1, fp);
          if (n != 1) {
            EG_free(otemp);
            EG_FREE(lloop_h->edges);
            EG_FREE(lloop_h->senses);
            EG_FREE(lloop);
            return EGADS_READERR;
          }
          if (j < lloop_h->nedges) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
            EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->edges.objs[iref-1]));
#endif
/*@+nullderef@*/
          } else {
/*@-nullderef@*/
#ifndef __clang_analyzer__
            EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->pcurves.objs[iref-1]));
#endif
/*@+nullderef@*/
          }
        }
        EG_COPY(lloop_h->edges, otemp, egObject *, m);
        if (lloop_h->surface != NULL) {
          liteEdge ledge_, *ledge_h = &ledge_;
          for (n = 0; n < lloop_h->nedges; n++) {
            EG_GET_OBJECT_PTR(&ledge, &(otemp[n]->blind));
            EG_GET_EDGE(ledge_h, ledge);
            pcobj = otemp[lloop_h->nedges+n];
#ifdef __NVCC__
            stat  = EG_attributeRetDev(pcobj, ".Bad", &atype, &alen,
                                       NULL, NULL, NULL);
            if (stat != EGADS_SUCCESS) continue;
            EG_evaluateDev(pcobj, &ledge_h->trange[0], data);
#else
            stat  = EG_attributeRet(pcobj, ".Bad", &atype, &alen,
                                    &ints, &reals, &str);
            if (stat != EGADS_SUCCESS) continue;
            EG_evaluate(pcobj, &ledge_h->trange[0], data);
#endif
            d     = sqrt(data[2]*data[2] + data[3]*data[3]);
            x0[0] = x0[1] = 0.0;
            if (d != 0.0) {
              x0[0] = data[2]/d;
              x0[1] = data[3]/d;
            }
            for (j = 1; j < 1000; j++) {
              t = ledge_h->trange[0]+j*(ledge_h->trange[1]-ledge_h->trange[0])/999.;
#ifdef __NVCC__
              EG_evaluateDev(pcobj, &t, data);
#else
              EG_evaluate(pcobj, &t, data);
#endif
              d = sqrt(data[2]*data[2] + data[3]*data[3]);
              x1[0] = x1[1] = 0.0;
              if (d != 0.0) {
                x1[0] = data[2]/d;
                x1[1] = data[3]/d;
              }
              if (x0[0]*x1[0] + x0[1]*x1[1] < -0.95) {
#ifdef __NVCC__
                stat = EG_addStrAttrDev(pcobj, ".Bad", "fold");
#else
                stat = EG_addStrAttr(pcobj, ".Bad", "fold");
#endif
                if (stat != EGADS_SUCCESS)
                  printf(" EGADS Info: EG_addStrAttr fold= %d\n", stat);
              }
              x0[0] = x1[0];
              x0[1] = x1[1];
            }
          }
        }
        EG_free(otemp);
      }
/*@-nullret@*/
      EG_SET_LOOP(lloop, lloop_h);
/*@+nullret@*/

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->loops.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = LOOP;
      obj_h->mtype  = mtype;
      obj_h->blind  = lloop;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* faces */
#ifdef DEBUG
  printf(" Reading %d Faces...\n", ntypes[6]);
#endif
  if (lbody_h->faces.objs != NULL) {
    liteFace lface_, *lface_h = &lface_;
#ifdef DEBUG
    liteLoop lloop_, *lloop_h = &lloop_;
    liteLoop *lloop_d;
#endif
    int      *itemp;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->faces.nobjs; i++) {
      egObject **otemp;
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lface, liteFace, 1);
      if (lface == NULL) return EGADS_MALLOC;
      lface_h->nloops  = m;
      lface_h->surface = NULL;
      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(lface_h->surface),
                          &(lbody_h->surfaces.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
      n = Fread(lface_h->urange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(lface_h->vrange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(lface_h->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(&lface_h->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      EG_NEW(&(lface_h->senses), int, lface_h->nloops);
      if (lface_h->senses == NULL) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      itemp = (int *) EG_alloc(lface_h->nloops*sizeof(int));
      if (itemp == NULL) {
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_MALLOC;
      }
      n = Fread((int *) itemp,  sizeof(int), lface_h->nloops, fp);
      if (n != lface_h->nloops) {
        EG_free(itemp);
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_READERR;
      }
      EG_COPY(lface_h->senses, itemp, int, lface_h->nloops);
      EG_free(itemp);
      EG_NEW(&(lface_h->loops), egObject *, lface_h->nloops);
      if (lface_h->loops == NULL) {
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_MALLOC;
      }
      otemp = (egObject **) EG_alloc(lface_h->nloops*sizeof(egObject *));
      for (j = 0; j < lface_h->nloops; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(otemp);
          EG_FREE(lface_h->loops);
          EG_FREE(lface_h->senses);
          EG_FREE(lface);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->loops.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
/*@-nullpass@*/
      EG_COPY(lface_h->loops, otemp, egObject *, lface_h->nloops);
/*@+nullpass@*/
      EG_free(otemp);

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->faces.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = FACE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lface;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
#ifdef DEBUG
      for (j = 0; j < lface_h->nloops; j++) {
        EG_GET_OBJECT_PTR(&(obj), &(lface_h->loops[j]));
        EG_GET_OBJECT_PTR(&(obj), &(obj->blind));
        lloop_d = (liteLoop *) obj;
        EG_GET_LOOP(lloop_h, lloop_d);
        if (lloop_h->surface == NULL) continue;
        for (n = 0; n < lloop_h->nedges; n++) {
          EG_GET_OBJECT_PTR(&(obj), &(lloop_h->edges[lloop_h->nedges+n]));
#ifdef __NVCC__
          stat  = EG_attributeRetDev(obj, ".Bad", &atype, &alen,
                                     NULL, NULL, &str);
#else
          stat = EG_attributeRet(obj, ".Bad", &atype,
                                 &alen, &ints, &reals, &str);
#endif
          if ((stat == EGADS_SUCCESS) && (atype == ATTRSTRING)) {
            if (strcmp(str, "fold") == 0)
              printf(" EGADS Info: Body %d Face %d Loop#%d/Edge#%d Bad PCurve -- %s!\n",
                     bindex+1, i+1, j+1, n+1, str);
#ifdef __NVCC__
            EG_free(str);
#endif
          }
        }
      }
#endif
/*@-nullret@*/
      EG_SET_FACE(lface, lface_h);
/*@+nullret@*/
    }
  }

  /* shells */
#ifdef DEBUG
  printf(" Reading %d Shells...\n", ntypes[7]);
#endif
  if (lbody_h->shells.objs != NULL) {
    liteShell lshell_, *lshell_h = &lshell_;
    egObject  obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->shells.nobjs; i++) {
      egObject  **otemp;
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lshell, liteShell, 1);
      if (lshell == NULL) return EGADS_MALLOC;
      lshell_h->nfaces = m;
      n = Fread(lshell_h->bbox, sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lshell);
        return EGADS_READERR;
      }
      EG_NEW(&(lshell_h->faces), egObject *, m);
      if (lshell_h->faces == NULL) {
        EG_FREE(lshell);
        return EGADS_MALLOC;
      }
      otemp = (egObject **) EG_alloc(m*sizeof(egObject *));
      if (otemp == NULL) {
        EG_FREE(lshell_h->faces);
        EG_FREE(lshell);
        return EGADS_MALLOC;
      }
      for (j = 0; j < m; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(otemp);
          EG_FREE(lshell_h->faces);
          EG_FREE(lshell);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->faces.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
      EG_COPY(lshell_h->faces, otemp, egObject *, m);
      EG_free(otemp);
      EG_SET_SHELL(lshell, lshell_h);

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->shells.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = SHELL;
      obj_h->mtype  = mtype;
      obj_h->blind  = lshell;
      obj_h->topObj = bobj;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* finish off the body */

  if (lbody_h->shells.nobjs != 0) {
    int *itemp;
    EG_NEW(&(lbody_h->senses), int, lbody_h->shells.nobjs);
    if (lbody_h->senses == NULL) return EGADS_MALLOC;
    itemp = (int *) EG_alloc(lbody_h->shells.nobjs*sizeof(int));
    if (itemp == NULL) {
      EG_FREE(lbody_h->senses);
      return EGADS_MALLOC;
    }
    n = Fread(itemp, sizeof(int), lbody_h->shells.nobjs, fp);
    if (n != lbody_h->shells.nobjs) {
      EG_free(itemp);
      EG_FREE(lbody_h->senses);
      return EGADS_READERR;
    }
    EG_COPY(lbody_h->senses, itemp, int, lbody_h->shells.nobjs);
    EG_free(itemp);
  }
  n = Fread(lbody_h->bbox, sizeof(double), 6, fp);
  if (n != 6) return EGADS_READERR;
  EG_GET_OBJECT(bobj_h, bobj);
  stat = EG_readAttrs(fp, (egAttrs **) &bobj_h->attrs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_SET_OBJECT(&bobj, bobj_h);
/*@-nullret@*/
  EG_SET_BODY(lbody, lbody_h);
/*@+nullret@*/

  return EGADS_SUCCESS;
}


/* this function needs CUDA attention! */
static int
EG_uvmapImport(void **uvmap, int **trmap, double *range, stream_T *fp)
{
  int          i, stat, msrch, trmp, *map;
  uvmap_struct *uvstruct;

  *uvmap = NULL;
  *trmap = NULL;
  uvstruct = (uvmap_struct *) EG_alloc(sizeof(uvmap_struct));
  if (uvstruct == NULL) {
    printf(" EGADS Error: Failed to allocate UVmap (EG_uvmapImport)!\n");
    return EGADS_MALLOC;
  }
  uvstruct->ndef   = 1;
  uvstruct->mdef   = 1;
  uvstruct->idef   = 1;
  uvstruct->isrch  = 0;
  uvstruct->ibface = 0;
  uvstruct->nbface = 0;
  uvstruct->nnode  = 0;

  uvstruct->idibf  = NULL;
  uvstruct->msrch  = NULL;
  uvstruct->inibf  = NULL;
  uvstruct->ibfibf = NULL;
  uvstruct->u      = NULL;

  stat = EGADS_READERR;
  if (Fread(&uvstruct->isrch,  sizeof(int), 1, fp) != 1) goto errOut;
  if (Fread(&uvstruct->ibface, sizeof(int), 1, fp) != 1) goto errOut;
  if (Fread(&uvstruct->nbface, sizeof(int), 1, fp) != 1) goto errOut;
  if (Fread(&uvstruct->nnode,  sizeof(int), 1, fp) != 1) goto errOut;
  if (Fread(&msrch,            sizeof(int), 1, fp) != 1) goto errOut;
  if (Fread(&trmp,             sizeof(int), 1, fp) != 1) goto errOut;

  uvstruct->idibf = (int *) EG_alloc((uvstruct->nbface+1)*sizeof(int));
  if (uvstruct->idibf == NULL) {
    printf(" EGADS Error: malloc %d id (EG_uvmapImport)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto errOut;
  }
  uvstruct->idibf[0] = 0;
  if (Fread(&uvstruct->idibf[1], sizeof(int), uvstruct->nbface, fp) !=
      uvstruct->nbface) goto errOut;

  uvstruct->inibf = (INT_3D *) EG_alloc((uvstruct->nbface+1)*sizeof(INT_3D));
  if (uvstruct->inibf == NULL) {
    printf(" EGADS Error: malloc %d ini (EG_uvmapImport)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto errOut;
  }
  uvstruct->ibfibf = (INT_3D *) EG_alloc((uvstruct->nbface+1)*sizeof(INT_3D));
  if (uvstruct->ibfibf == NULL) {
    printf(" EGADS Error: malloc %d ibf (EG_uvmapImport)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto errOut;
  }
  uvstruct->inibf[0][0]  = uvstruct->inibf[0][1]  = uvstruct->inibf[0][2]  = 0;
  uvstruct->ibfibf[0][0] = uvstruct->ibfibf[0][1] = uvstruct->ibfibf[0][2] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    if (Fread(uvstruct->inibf[i],  sizeof(int), 3, fp) != 3) goto errOut;
    if (Fread(uvstruct->ibfibf[i], sizeof(int), 3, fp) != 3) goto errOut;
  }
  uvstruct->u = (DOUBLE_2D *) EG_alloc((uvstruct->nnode+1)*sizeof(DOUBLE_2D));
  if (uvstruct->u == NULL) {
    printf(" EGADS Error: malloc %d u (EG_uvmapImport)!\n", uvstruct->nnode);
    stat = EGADS_MALLOC;
    goto errOut;
  }
  uvstruct->u[0][0] = uvstruct->u[0][1] = 0.0;
  for (i = 1; i <= uvstruct->nnode; i++)
    if (Fread(uvstruct->u[i], sizeof(double), 2, fp) != 2) goto errOut;

  range[0] = range[1] = uvstruct->u[1][0];
  range[2] = range[3] = uvstruct->u[1][1];
  for (i = 2; i <= uvstruct->nnode; i++) {
    if (uvstruct->u[i][0] < range[0]) range[0] = uvstruct->u[i][0];
    if (uvstruct->u[i][0] > range[1]) range[1] = uvstruct->u[i][0];
    if (uvstruct->u[i][1] < range[2]) range[2] = uvstruct->u[i][1];
    if (uvstruct->u[i][1] > range[3]) range[3] = uvstruct->u[i][1];
  }

  if (msrch == 1) {
    uvstruct->msrch = (int *) EG_alloc((uvstruct->nbface+1)*sizeof(int));
    if (uvstruct->msrch == NULL) {
      printf(" EGADS Error: malloc %d msrch (EG_uvmapRead)!\n",
             uvstruct->nbface);
      stat = EGADS_MALLOC;
      goto errOut;
    }
    uvstruct->msrch[0] = 0;
    if (Fread(&uvstruct->msrch[1], sizeof(int), uvstruct->nbface, fp) !=
        uvstruct->nbface) goto errOut;
  }

  if (trmp == 1) {
    map = (int *) EG_alloc(uvstruct->nbface*sizeof(int));
    if (map == NULL) {
      printf(" EGADS Error: malloc %d trmap (EG_uvmapRead)!\n",
             uvstruct->nbface);
      stat = EGADS_MALLOC;
      goto errOut;
    }
    if (Fread(map, sizeof(int), uvstruct->nbface, fp) != uvstruct->nbface)
      goto errOut;
    *trmap = map;
  }

  *uvmap = uvstruct;
  return EGADS_SUCCESS;

errOut:
  EG_free(uvstruct->idibf);
  EG_free(uvstruct->msrch);
  EG_free(uvstruct->inibf);
  EG_free(uvstruct->ibfibf);
  EG_free(uvstruct->u);
/*@-dependenttrans@*/
  EG_free(uvstruct);
/*@+dependenttrans@*/
  return stat;
}


static int
EG_readEBody(egObject *context, egObject *mobject, int bindex, int iref,
             stream_T *fp)
{
  int       index, i, j, n, stat, mtype, nsegs, nds[2];
  double    area;
  egEBody   *ebody;
  egEShell  *eshell;
  egEFace   *eface;
  egELoop   *eloop;
  egEEdge   *eedge;
  egObject  bobj_,    *bobj_h    = &bobj_,    *tobj, *bobj, *body;
  egObject  context_, *context_h = &context_;
  egObject  mobject_, *mobject_h = &mobject_;
  liteModel lmodel_,  *lmodel_h  = &lmodel_,  *lmodel;

  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (mobject == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(mobject_h, mobject);
  if (mobject_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject_h->oclass != MODEL)      return EGADS_NOTMODEL;
  lmodel = (liteModel *) mobject_h->blind;
  EG_GET_MODEL(lmodel_h, lmodel);
  body   = lmodel_h->bodies[iref-1];

  n = Fread(&mtype, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;

  EG_NEW(&ebody, egEBody, 1);
  if (ebody == NULL) return EGADS_MALLOC;
  ebody->ref           = body;
  ebody->eedges.objs   = NULL;
  ebody->eloops.objs   = NULL;
  ebody->efaces.objs   = NULL;
  ebody->eshells.objs  = NULL;
  ebody->senses        = NULL;
  ebody->eedges.nobjs  = 0;
  ebody->eloops.nobjs  = 0;
  ebody->efaces.nobjs  = 0;
  ebody->eshells.nobjs = 0;
  ebody->angle         = 0.0;
  ebody->done          = 1;
  ebody->nedge         = 0;
  ebody->edges         = NULL;

  stat = EG_makeObject(context, &bobj);
  if (stat != EGADS_SUCCESS) {
    EG_FREE(ebody);
    return stat;
  }
  EG_GET_OBJECT(bobj_h, bobj);
  stat = EG_readAttrs(fp, (egAttrs **) &bobj_h->attrs);
  if (stat != EGADS_SUCCESS) {
    EG_FREE(ebody);
    return stat;
  }
  bobj_h->oclass = EBODY;
  bobj_h->mtype  = mtype;
  bobj_h->blind  = ebody;
  bobj_h->topObj = mobject;
  EG_SET_OBJECT(&bobj, bobj_h);
  EG_SET_OBJECT_PTR(&(lmodel_h->bodies[bindex]), &bobj);

  /***** needs attention for CUDA *****/
  n = Fread(&ebody->eedges.nobjs,  sizeof(int),    1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ebody->eloops.nobjs,  sizeof(int),    1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ebody->efaces.nobjs,  sizeof(int),    1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ebody->eshells.nobjs, sizeof(int),    1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ebody->nedge,         sizeof(int),    1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ebody->angle,         sizeof(double), 1, fp);
  if (n != 1) return EGADS_READERR;

  if (mtype == SOLIDBODY) {
    ebody->senses = (int *) EG_alloc(ebody->eshells.nobjs*sizeof(int));
    if (ebody->senses == NULL) {
      printf(" EGADS Error: Cannot allocate %d Senses (EG_importEBody)!\n",
             ebody->eshells.nobjs);
      return EGADS_MALLOC;
    }
    n = Fread(ebody->senses,       sizeof(int), ebody->eshells.nobjs, fp);
    if (n != ebody->eshells.nobjs) return EGADS_READERR;
  }

  /* source Edges */
  ebody->edges = (egEdVert *) EG_alloc(ebody->nedge*sizeof(egEdVert));
  if (ebody->edges == NULL) {
    printf(" EGADS Error: Cannot allocate %d EdVerts (EG_importEBody)!\n",
           ebody->nedge);
    return EGADS_MALLOC;
  }
  for (j = 0; j < ebody->nedge; j++) {
    ebody->edges[j].edge = NULL;
    ebody->edges[j].npts = 0;
    ebody->edges[j].ts   = NULL;
  }
  for (j = 0; j < ebody->nedge; j++) {
    n = Fread(&index,                 sizeof(int),    1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&ebody->edges[j].curve, sizeof(int),    1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&ebody->edges[j].npts,  sizeof(int),    1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(ebody->edges[j].dstart, sizeof(double), 3, fp);
    if (n != 3) return EGADS_READERR;
    n = Fread(ebody->edges[j].dend,   sizeof(double), 3, fp);
    if (n != 3) return EGADS_READERR;
    ebody->edges[j].ts = EG_alloc(ebody->edges[j].npts*sizeof(double));
    if (ebody->edges[j].ts == NULL) {
      printf(" EGADS Error: Edge %d Cannot allocate %d ts (EG_importEBody)!\n",
             j+1, ebody->edges[j].npts);
      return EGADS_MALLOC;
    }
    n = Fread(ebody->edges[j].ts,     sizeof(double), ebody->edges[j].npts, fp);
    if (n != ebody->edges[j].npts) return EGADS_READERR;
    stat = EG_objectBodyTopo(body, EDGE, index, &ebody->edges[j].edge);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Source Edge = %d (EG_exportEBody)!\n", index);
      return EGADS_TOPOERR;
    }
  }

  /* EEdges */
  ebody->eedges.objs = (egObject **)
                       EG_alloc(ebody->eedges.nobjs*sizeof(egObject *));
  if (ebody->eedges.objs == NULL) {
    printf(" EGADS Error: Cannot allocate %d EEdges (EG_importEBody)!\n",
           ebody->eedges.nobjs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eedges.nobjs; i++) ebody->eedges.objs[i] = NULL;
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    n = Fread(&mtype, sizeof(short),  1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&nsegs, sizeof(int),    1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(nds,    sizeof(int),    2, fp);
    if (n != 2) return EGADS_READERR;
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EEdge %d/%d Object (EG_importEBody)!\n",
             i+1, ebody->eedges.nobjs);
      return stat;
    }
    tobj->oclass          = EEDGE;
    tobj->mtype           = mtype;
    tobj->topObj          = bobj;
    ebody->eedges.objs[i] = tobj;
    eedge = (egEEdge *) EG_alloc(sizeof(egEEdge));
    if (eedge == NULL) {
      printf(" EGADS Error: Malloc on %d EEdge blind (EG_importEBody)!\n", i+1);
      return EGADS_MALLOC;
    }
    tobj->blind   = eedge;
    n = Fread(eedge->trange, sizeof(double), 2, fp);
    if (n != 2) return EGADS_READERR;
    eedge->sedges = ebody->edges;
    stat = EG_objectBodyTopo(body, NODE, nds[0], &eedge->nodes[0]);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EEdge = %d First Node (EG_importEBody)!\n", i+1);
      return EGADS_TOPOERR;
    }
    stat = EG_objectBodyTopo(body, NODE, nds[1], &eedge->nodes[1]);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EEdge = %d Last Node (EG_importEBody)!\n", i+1);
      return EGADS_TOPOERR;
    }
    eedge->nsegs = nsegs;
    eedge->segs  = (egEEseg *) EG_alloc(eedge->nsegs*sizeof(egEEseg));
    if (eedge->segs == NULL) {
      printf(" EGADS Error: Malloc EEdge %d nsegs = %d (EG_importEBody)!\n",
             i+1, eedge->nsegs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eedge->nsegs; j++) eedge->segs[j].nstart = NULL;
    for (j = 0; j < eedge->nsegs; j++) {
      n = Fread(&eedge->segs[j].iedge,  sizeof(int),    1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eedge->segs[j].sense,  sizeof(int),    1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&index,                 sizeof(int),    1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eedge->segs[j].tstart, sizeof(double), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eedge->segs[j].tend,   sizeof(double), 1, fp);
      if (n != 1) return EGADS_READERR;
      if (index != 0) {
        stat = EG_objectBodyTopo(body, NODE, index, &eedge->segs[j].nstart);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: EEdge = %d Start Node %d (EG_importEBody)!\n",
                 i+1, j+1);
          return EGADS_TOPOERR;
        }
      }
    }
    stat = EG_readAttrs(fp, (egAttrs **) &tobj->attrs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: readAttrs = %d  EEdge = %d (EG_importEBody)!\n",
             stat, i+1);
      return stat;
    }
  }

  /* ELoops */
  ebody->eloops.objs = (egObject **)
                       EG_alloc(ebody->eloops.nobjs*sizeof(egObject *));
  if (ebody->eloops.objs == NULL) {
    printf(" EGADS Error: Cannot allocate %d ELoops (EG_importEBody)!\n",
           ebody->eloops.nobjs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eloops.nobjs; i++) ebody->eloops.objs[i] = NULL;
  for (i = 0; i < ebody->eloops.nobjs; i++) {
    n = Fread(&mtype, sizeof(short),  1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(nds,    sizeof(int),    2, fp);
    if (n != 2) return EGADS_READERR;
    n = Fread(&area,  sizeof(double), 1, fp);
    if (n != 1) return EGADS_READERR;
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make ELoop %d/%d Object (EG_importEBody)!\n",
             i+1, ebody->eloops.nobjs);
      return stat;
    }
    tobj->oclass          = ELOOPX;
    tobj->mtype           = mtype;
/*@-kepttrans@*/
    tobj->topObj          = bobj;
/*@+kepttrans@*/
    ebody->eloops.objs[i] = tobj;
    tobj->blind           = NULL;
    if ((nds[0] == 0) && (nds[1] == 0)) continue;
    eloop = (egELoop *) EG_alloc(sizeof(egELoop));
    if (eloop == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop blind (EG_importEBody)!\n", i+1);
      return EGADS_MALLOC;
    }
    tobj->blind         = eloop;
    eloop->eedges.nobjs = nds[0];
    eloop->nedge        = nds[1];
    eloop->edgeUVs      = NULL;
    eloop->senses       = NULL;
    eloop->area         = area;
    eloop->eedges.objs  = (egObject **) EG_alloc(eloop->eedges.nobjs*
                                                 sizeof(egObject *));
    if (eloop->eedges.objs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d Objects (EG_importEBody)!\n",
             i+1, eloop->eedges.nobjs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->eedges.nobjs; j++) eloop->eedges.objs[j] = NULL;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      n = Fread(&index, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      stat = EG_objectBodyTopo(bobj, EEDGE, index, &eloop->eedges.objs[j]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: ELoop = %d Edge %d  %d  (EG_importEBody)!\n",
               i+1, j+1, index);
        return EGADS_TOPOERR;
      }
    }
    eloop->senses = (int *) EG_alloc(eloop->eedges.nobjs*sizeof(int));
    if (eloop->senses == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d senses (EG_importEBody)!\n",
             i+1, eloop->eedges.nobjs);
      return EGADS_MALLOC;
    }
    n = Fread(eloop->senses, sizeof(int), eloop->eedges.nobjs, fp);
    if (n != eloop->eedges.nobjs) return EGADS_READERR;

    eloop->edgeUVs = (egEdgeUV *) EG_alloc(eloop->nedge*sizeof(egEdgeUV));
    if (eloop->edgeUVs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d edgeUVs (EG_importEBody)!\n",
             i+1, eloop->nedge);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->nedge; j++) {
      n = Fread(&index,                   sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eloop->edgeUVs[j].sense, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eloop->edgeUVs[j].npts,  sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      eloop->edgeUVs[j].iuv = (int *) EG_alloc(eloop->edgeUVs[j].npts*
                                               sizeof(int));
      if (eloop->edgeUVs == NULL) {
        printf(" EGADS Error: Malloc on %d ELoop %d iUVs %d (EG_importEBody)!\n",
               i+1, j+1, eloop->edgeUVs[j].npts);
        return EGADS_MALLOC;
      }
      n = Fread(eloop->edgeUVs[j].iuv, sizeof(int), eloop->edgeUVs[j].npts, fp);
      if (n != eloop->edgeUVs[j].npts) return EGADS_READERR;
      stat = EG_objectBodyTopo(body, EDGE, index, &eloop->edgeUVs[j].edge);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: ELoop = %d EdgeUV %d  %d  (EG_importEBody)!\n",
               i+1, j+1, index);
        return EGADS_TOPOERR;
      }
    }
    stat = EG_readAttrs(fp, (egAttrs **) &tobj->attrs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: readAttrs = %d  ELoop = %d (EG_importEBody)!\n",
             stat, i+1);
      return stat;
    }
  }

  /* EFaces */
  ebody->efaces.objs = (egObject **)
                       EG_alloc(ebody->efaces.nobjs*sizeof(egObject *));
  if (ebody->efaces.objs == NULL) {
    printf(" EGADS Error: Cannot allocate %d EFaces (EG_importEBody)!\n",
           ebody->efaces.nobjs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->efaces.nobjs; i++) ebody->efaces.objs[i] = NULL;
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    n = Fread(&mtype, sizeof(short), 1, fp);
    if (n != 1) return EGADS_READERR;
    eface = (egEFace *) EG_alloc(sizeof(egEFace));
    if (eface == NULL) {
      printf(" EGADS Error: Malloc on %d EFace blind (EG_importEBody)!\n", i+1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EFace %d/%d Object (EG_importEBody)!\n",
             i+1, ebody->efaces.nobjs);
      EG_free(eface);
      return stat;
    }
    ebody->efaces.objs[i] = tobj;
    tobj->oclass          = EFACE;
    tobj->mtype           = mtype;
/*@-kepttrans@*/
    tobj->topObj          = bobj;
/*@+kepttrans@*/
    tobj->blind           = eface;
    eface->trmap          = NULL;
    eface->uvmap          = NULL;
    eface->patches        = NULL;
    eface->eloops.objs    = NULL;
    eface->senses         = NULL;
    eface->patches        = NULL;
    n = Fread(&eface->npatch,       sizeof(int),   1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&eface->eloops.nobjs, sizeof(int),   1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&eface->last,         sizeof(int),   1, fp);
    if (n != 1) return EGADS_READERR;
    if (eface->npatch != 1) {
      stat = EG_uvmapImport(&eface->uvmap, &eface->trmap, eface->range, fp);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d  uvmapImport = %d (EG_importEBody)!\n",
               i+1, stat);
        return stat;
      }
    } else {
      n = Fread(eface->range, sizeof(double), 4, fp);
      if (n != 4) return EGADS_READERR;
    }
    eface->eloops.objs = (egObject **) EG_alloc(eface->eloops.nobjs*
                                                sizeof(egObject *));
    if (eface->eloops.objs == NULL) {
      printf(" EGADS Error: Malloc on %d EFace %d Objects (EG_importEBody)!\n",
             i+1, eface->eloops.nobjs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->eloops.nobjs; j++) eface->eloops.objs[j] = NULL;
    for (j = 0; j < eface->eloops.nobjs; j++) {
      n = Fread(&index, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      stat = EG_objectBodyTopo(bobj, ELOOPX, index, &eface->eloops.objs[j]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace = %d ELoop %d  %d  (EG_importEBody)!\n",
               i+1, j+1, index);
        return EGADS_TOPOERR;
      }
    }
    eface->senses = (int *) EG_alloc(eface->eloops.nobjs*sizeof(int));
    if (eface->senses == NULL) {
      printf(" EGADS Error: Malloc on %d EFace senses %d (EG_importEBody)!\n",
             i+1, eface->eloops.nobjs);
      return EGADS_MALLOC;
    }
    n = Fread(eface->senses, sizeof(int), eface->eloops.nobjs, fp);
    if (n != eface->eloops.nobjs) return EGADS_READERR;

    eface->patches = (egEPatch *) EG_alloc(abs(eface->npatch)*sizeof(egEPatch));
    if (eface->patches == NULL) {
      printf(" EGADS Error: Malloc on %d Patch Object %d (EG_importEBody)!\n",
             i+1, abs(eface->npatch));
      return EGADS_MALLOC;
    }
    for (j = 0; j < abs(eface->npatch); j++) {
      eface->patches[j].uvtris  = NULL;
      eface->patches[j].uvtric  = NULL;
      eface->patches[j].uvs     = NULL;
      eface->patches[j].deflect = NULL;
      eface->patches[j].tol     = -1.0;
    }
    for (j = 0; j < abs(eface->npatch); j++) {
      n = Fread(&index,                      sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eface->patches[j].start,    sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eface->patches[j].nuvs,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eface->patches[j].ndeflect, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&eface->patches[j].ntris,    sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      eface->patches[j].uvtris = (int *) EG_alloc(3*eface->patches[j].ntris*
                                                  sizeof(int));
      if (eface->patches[j].uvtris == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvtris %d (EG_importEBody)!\n",
               i+1, eface->patches[j].ntris);
        return EGADS_MALLOC;
      }
      n = Fread(eface->patches[j].uvtris,    sizeof(int),
                3*eface->patches[j].ntris, fp);
      if (n != 3*eface->patches[j].ntris) return EGADS_READERR;
      eface->patches[j].uvs = (double *) EG_alloc(2*eface->patches[j].nuvs*
                                                  sizeof(double));
      if (eface->patches[j].uvs == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvs %d (EG_importEBody)!\n",
               i+1, eface->patches[j].nuvs);
        return EGADS_MALLOC;
      }
      n = Fread(eface->patches[j].uvs,       sizeof(double),
                2*eface->patches[j].nuvs, fp);
      if (n != 2*eface->patches[j].nuvs) return EGADS_READERR;
      eface->patches[j].deflect = (double *)
                          EG_alloc(3*eface->patches[j].ndeflect*sizeof(double));
      if (eface->patches[j].deflect == NULL) {
        printf(" EGADS Error: Malloc on %d Patch deflect %d (EG_importEBody)!\n",
               i+1, eface->patches[j].ndeflect);
        return EGADS_MALLOC;
      }
      n = Fread(eface->patches[j].deflect,   sizeof(double),
                3*eface->patches[j].ndeflect, fp);
      if (n != 3*eface->patches[j].ndeflect) return EGADS_READERR;

      stat = EG_objectBodyTopo(body, FACE, index, &eface->patches[j].face);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace = %d Patch %d  %d  (EG_importEBody)!\n",
               i+1, j+1, index);
        return EGADS_TOPOERR;
      }
    }
    if (eface->npatch == 1) {
      stat = EG_effectNeighbor(eface);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: %d Patch effectNeighbor = %d (EG_importEBody)!\n",
               i+1, stat);
        return EGADS_MALLOC;
      }
    }
    stat = EG_readAttrs(fp, (egAttrs **) &tobj->attrs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: readAttrs = %d  EFace = %d (EG_importEBody)!\n",
             stat, i+1);
      return stat;
    }
  }

  /* EShells */
  ebody->eshells.objs = (egObject **) EG_alloc(ebody->eshells.nobjs*
                                               sizeof(egObject *));
  if (ebody->eshells.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EShell Objects (EG_importEBody)!\n",
           ebody->eshells.nobjs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eshells.nobjs; i++) ebody->eshells.objs[i] = NULL;
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    n = Fread(&mtype, sizeof(short), 1, fp);
    if (n != 1) return EGADS_READERR;

    eshell = (egEShell *) EG_alloc(sizeof(egEShell));
    if (eshell == NULL) {
      printf(" EGADS Error: Malloc on %d EShell blind (EG_importEBody)!\n", i+1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EShell %d/%d Object (EG_importEBody)!\n",
             i+1, ebody->eshells.nobjs);
      EG_free(eshell);
      return stat;
    }
    tobj->oclass           = ESHELL;
    tobj->mtype            = mtype;
/*@-kepttrans@*/
    tobj->topObj           = bobj;
/*@+kepttrans@*/
    ebody->eshells.objs[i] = tobj;
    tobj->blind            = eshell;
    n = Fread(&eshell->efaces.nobjs, sizeof(int),   1, fp);
    if (n != 1) return EGADS_READERR;
    eshell->efaces.objs    = (egObject **) EG_alloc(eshell->efaces.nobjs*
                                                    sizeof(egObject *));
    if (eshell->efaces.objs == NULL) {
      printf(" EGADS Error: Malloc on %d EShell %d Objects (EG_importEBody)!\n",
             i+1, eshell->efaces.nobjs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eshell->efaces.nobjs; j++) eshell->efaces.objs[j] = NULL;
    for (j = 0; j < eshell->efaces.nobjs; j++) {
      n = Fread(&index, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      stat = EG_objectBodyTopo(bobj, EFACE, index, &eshell->efaces.objs[j]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EShell = %d EFace %d  %d  (EG_importEBody)!\n",
               i+1, j+1, index);
        return EGADS_TOPOERR;
      }
    }
    stat = EG_readAttrs(fp, (egAttrs **) &tobj->attrs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: readAttrs = %d  EShell = %d (EG_importEBody)!\n",
             stat, i+1);
      return stat;
    }
  }

  /***** ************************ *****/

  return EGADS_SUCCESS;
}


static int
EG_readTess(egObject *mobject, int bindex, int iref, stream_T *fp)
{
  int       i, n, stat, nedge, nface, len, ntri, *tris;
  double    *xyzs, *uvs, *ts;
  egObject  bobj_,    *bobj_h    = &bobj_,    *bobj, *ref;
  egObject  mobject_, *mobject_h = &mobject_;
  liteModel lmodel_,  *lmodel_h  = &lmodel_,  *lmodel;

  if (mobject == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(mobject_h, mobject);
  if (mobject_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject_h->oclass != MODEL)      return EGADS_NOTMODEL;
  lmodel = (liteModel *) mobject_h->blind;
  EG_GET_MODEL(lmodel_h, lmodel);

  n = Fread(&nedge, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&nface, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;

  stat = EG_initTessBody(lmodel_h->bodies[iref-1], &bobj);
  if (stat != EGADS_SUCCESS) return stat;
  EG_GET_OBJECT(bobj_h, bobj);
  bobj_h->topObj = mobject;
  EG_SET_OBJECT(&bobj, bobj_h);
  EG_SET_OBJECT_PTR(&(lmodel_h->bodies[bindex]), &bobj);

  /* get Edges */
  for (i = 0; i < nedge; i++) {
    n = Fread(&len, sizeof(int), 1, fp);
    if (n != 1) return EGADS_READERR;
#ifdef DEBUG
    printf(" Reading Edge Tess %d len = %d\n", i+1, len);
#endif
    if (len == 0) continue;
    xyzs = (double *) EG_alloc(3*len*sizeof(double));
    ts   = (double *) EG_alloc(  len*sizeof(double));
    if ((xyzs == NULL) || (ts == NULL)) {
      printf(" EGADS Error: Alloc %d Edge %d (EG_importModel)!\n", i+1, len);
      return EGADS_MALLOC;
    }
    n = Fread(xyzs, sizeof(double), 3*len, fp);
    if (n != 3*len) {
      EG_free(ts);
      EG_free(xyzs);
      return EGADS_READERR;
    }
    n = Fread(ts,   sizeof(double),   len, fp);
    if (n !=   len) {
      EG_free(ts);
      EG_free(xyzs);
      return EGADS_READERR;
    }
    stat = EG_setTessEdge(bobj, i+1, len, xyzs, ts);
    EG_free(ts);
    EG_free(xyzs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Edge %d EG_setTessEdge = %d (EG_importModel)!\n",
             i+1, stat);
      return stat;
    }
  }

  /* get Faces */
  for (i = 0; i < nface; i++) {
    n = Fread(&len,  sizeof(int), 1, fp);
    if (n != 1) return EGADS_READERR;
    n = Fread(&ntri, sizeof(int), 1, fp);
    if (n != 1) return EGADS_READERR;
#ifdef DEBUG
    printf(" Reading Face Tess %d len = %d  ntri = %d\n", i+1, len, ntri);
#endif
    if ((len == 0) || (ntri == 0)) continue;
    xyzs = (double *) EG_alloc( 3*len*sizeof(double));
    uvs  = (double *) EG_alloc( 2*len*sizeof(double));
    tris = (int *)    EG_alloc(3*ntri*sizeof(int));
    if ((xyzs == NULL) || (uvs == NULL) || (tris == NULL)) {
      printf(" EGADS Error: Alloc %d Face %d %d (EG_importModel)!\n",
             i+1, len, ntri);
      return EGADS_MALLOC;
    }
    n = Fread(xyzs, sizeof(double), 3*len, fp);
    if (n !=  3*len) {
      EG_free(tris);
      EG_free(uvs);
      EG_free(xyzs);
      return EGADS_READERR;
    }
    n = Fread(uvs,  sizeof(double), 2*len, fp);
    if (n !=  2*len) {
      EG_free(tris);
      EG_free(uvs);
      EG_free(xyzs);
      return EGADS_READERR;
    }
    n = Fread(tris, sizeof(int),   3*ntri, fp);
    if (n != 3*ntri) {
      EG_free(tris);
      EG_free(uvs);
      EG_free(xyzs);
      return EGADS_READERR;
    }
    stat = EG_setTessFace(bobj, i+1, len, xyzs, uvs, ntri, tris);
    EG_free(tris);
    EG_free(uvs);
    EG_free(xyzs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d EG_setTessFace = %d (EG_importModel)!\n",
             i+1, stat);
      return stat;
    }
  }

  EG_GET_OBJECT(bobj_h, bobj);
  stat = EG_readAttrs(fp, (egAttrs **) &bobj_h->attrs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_SET_OBJECT(&bobj, bobj_h);

  stat = EG_statusTessBody(bobj, &ref, &i, &len);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: %d EG_statusTessBody = %d (EG_importModel)!\n",
           bindex, stat);
  } else {
    if (i != 1)
      printf(" EGADS Warning: Body %d Tess NOT closed (EG_importModel)!\n",
             bindex);
  }

  return EGADS_SUCCESS;
}


static int
EG_freeLiteModel(/*@null@*/ liteModel *model)
{
  egObject *bodies_d = NULL;

  if (NULL != model) {
    EG_COPY(&(bodies_d), &(model->bodies), egObject *, 1);
    if (NULL != bodies_d) EG_FREE(bodies_d);
    EG_FREE(model);
  }
  return EGADS_SUCCESS;
}


static int
EG_allocLiteModel(int nbody, double bbox[6], liteModel **model)
{
  int       i;
  liteModel obj_, *obj_h = &obj_;
  liteModel *obj = NULL;
  void      *nil = NULL;

  obj_h->nbody   = nbody;
  obj_h->bbox[0] = bbox[0]; obj_h->bbox[1] = bbox[1]; obj_h->bbox[2] = bbox[2];
  obj_h->bbox[3] = bbox[3]; obj_h->bbox[4] = bbox[4]; obj_h->bbox[5] = bbox[5];
  obj_h->bodies  = NULL;
  if (nbody > 0) {
    EG_NEW((void **)&(obj_h->bodies), egObject *, obj_h->nbody);
    if (obj_h->bodies == NULL) goto modelCleanup;
  }

  EG_NEW((void **)&(obj), liteModel, 1);
  if (obj == NULL) goto modelCleanup;
  EG_COPY(obj, obj_h, liteModel, 1);
  if (obj == NULL) goto modelCleanup;

  for (i = 0; i < obj_h->nbody; ++i) {
    egObject *ptr;
/*@-nullderef@*/
    EG_COPY(&(obj_h->bodies[i]), &(nil), void *, 1);
    EG_GET_OBJECT_PTR(&ptr, &(obj_h->bodies[i]));
/*@+nullderef@*/
    if (ptr != NULL) goto modelCleanup;
  }

  *model = obj;
  return EGADS_SUCCESS;

modelCleanup:
/*@-dependenttrans@*/
  EG_FREE(obj);
/*@+dependenttrans@*/
  EG_FREE(obj_h->bodies);
  return EGADS_MALLOC;
}


/*** also needs CUDA attention ***/
static int
EG_reallocLiteModel(egObject *mobject)
{
  int       i, mtype, nbody;
  egObject  **bodies;
  liteModel *lmodel;

  mtype  = mobject->mtype;
  lmodel = (liteModel *) mobject->blind;
  nbody  = lmodel->nbody;
  bodies = (egObject **) EG_reall(lmodel->bodies, mtype*sizeof(egObject *));
  if (bodies == NULL) return EGADS_MALLOC;
  lmodel->bodies = bodies;

  for (i = nbody; i < mtype; i++) lmodel->bodies[i] = NULL;

  return EGADS_SUCCESS;
}


int
EG_importModel(egObject *context, const size_t nbytes, const char *stream,
               egObject **model)
{
  int       i, j, n, oclass, mtype, iref, rev[2];
  liteModel *lmodel = NULL;
  egObject  obj_, *obj_h = &obj_;
  egObject  *obj;
  egObject  context_;
  egObject  *context_h = &context_;
  stream_T  myStream;
  stream_T  *fp = &myStream;
  int       nbody;
  double    bbox[6];

  *model = NULL;

  EG_GET_OBJECT(context_h, context);

  if (context_h == NULL)               return EGADS_NULLOBJ;
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (context_h->topObj != NULL)       return EGADS_EXISTS;

  fp->size = nbytes;
  fp->ptr  = 0;
  fp->data = (void *) stream;
  fp->swap = 0;

  /* get header */
  n = Fread(&i,     sizeof(int),    1, fp);
  if (n != 1) {
    return EGADS_READERR;
  }
  if (i != MAGIC) {
    swap(&i, sizeof(int));
    if (i != MAGIC) {
      printf(" EGADS Error: Not a EGADS Lite file!\n");
      return EGADS_READERR;
    }
    fp->swap = 1;
  }
  n = Fread(rev,    sizeof(int),    2, fp);
  if (n != 2) {
    return EGADS_READERR;
  }
  if (rev[0] != 1) {
    printf(" EGADS Error: EGADS Lite file revision = %d %d!\n", rev[0], rev[1]);
    return EGADS_READERR;
  }

  n = Fread(bbox,   sizeof(double), 6, fp);
  if (n != 6) {
    return EGADS_READERR;
  }
  n = Fread(&nbody, sizeof(int),    1, fp);
  if (n != 1) {
    return EGADS_READERR;
  }
  i = EG_allocLiteModel(nbody, bbox, &lmodel);
  if (i != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_allocLiteModel = %d!\n", i);
    return i;
  }

  i = EG_makeObject(context, &obj);
  if (i != EGADS_SUCCESS) {
    EG_freeLiteModel(lmodel);
    printf(" EGADS Error: makeObject on Model = %d!\n", i);
    return i;
  }
  EG_GET_OBJECT(obj_h, obj);
  obj_h->oclass = MODEL;
  obj_h->mtype  = 0;
  obj_h->blind  = lmodel;
  i = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
  if (i != EGADS_SUCCESS) {
    EG_close(context);
    return i;
  }
/*@-nullret@*/
  EG_SET_OBJECT(&obj, obj_h);
/*@+nullret@*/

  /* get all of the bodies */
  for (n = 0; n < nbody; n++) {
    i = EG_readBody(context, obj, n, fp);
    if (i == EGADS_SUCCESS) continue;
    /* errorred out -- cleanup */
    EG_close(context);
    return i;
  }

  if (rev[1] != 0) {
    n = Fread(&mtype, sizeof(int), 1, fp);
    if (n != 1) {
      mtype = 0;
      printf(" EGADS Info: Cannot read extended Model data!\n");
    }
    if (mtype > nbody) {
      obj->mtype = obj_h->mtype = mtype;
      i = EG_reallocLiteModel(obj);
      if (i != EGADS_SUCCESS) {
        /* errorred out -- cleanup */
        EG_close(context);
        return i;
      }
      for (j = nbody; j < mtype; j++) {
        n = Fread(&oclass, sizeof(int), 1, fp);
        if (n != 1) {
          EG_close(context);
          return EGADS_READERR;
        }
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_close(context);
          return EGADS_READERR;
        }
        if (oclass == TESSELLATION) {
          i = EG_readTess(obj, j, iref, fp);
          if (i != EGADS_SUCCESS) {
            /* errorred out -- cleanup */
            printf(" Import Error: %d Tess Entry in Model has status = %d!\n",
                   j+1, i);
            EG_close(context);
            return i;
          }
        } else if (oclass == EBODY) {
          i = EG_readEBody(context, obj, j, iref, fp);
          if (i != EGADS_SUCCESS) {
            /* errorred out -- cleanup */
            printf(" Import Error: %d EBody Entry in Model has status = %d!\n",
                   j+1, i);
            EG_close(context);
            return i;
          }
        } else {
          printf(" Import Error: %d Entry in Model has class = %d!\n",
                 j+1, oclass);
          EG_close(context);
          return EGADS_NOTTOPO;
        }
      }
    }
  }

  EG_SET_OBJECT_PTR(&(context->topObj), &obj);
  *model = obj;

  return EGADS_SUCCESS;
}
