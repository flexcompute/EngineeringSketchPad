/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *      Tesselate EGADS bodies and write:
 *         CART3D unstructured surface triangulation tri file
 *                with FaceID in Component field
 *         TESS   information about vertex/geometry ownership
 *         PLOT3D structured surface grid file from quadding algorithm
 *         PLOT3D structured surface grid file from direct uv parameters space 
 *                evaluation (untrimmed)
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include <math.h>
#include <string.h>

//#define INFACEOCC

#ifdef INFACEOCC
extern int EG_inFaceOCC( const egObject *face, double tol, const double *uv );
#endif
extern int aflr4egads( ego model, ego *tesses );

#ifdef WIN32
#define ULONG    unsigned long long
#define snprintf _snprintf
#else
#define ULONG    unsigned long
#endif

//#define DEBUG

#define FUZZ   1.e-14
#define PI     3.1415926535897931159979635
#define MAXUVS 1024                        /* max Us or Vs for untrimmed data */
#define MAXSEG 1023
#define CROSS(a,b,c)  a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                      a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                      a[2] = (b[0]*c[1]) - (b[1]*c[0])


typedef struct {
  int    nu;
  int    nv;
  double us[MAXUVS];
  double vs[MAXUVS];
} UVmap;


typedef struct parmSeg {
  struct parmSeg *prev;          /* for front -- seg connection or NULL */
  double         parms[2];       /* for front -- dir is second value (-1,1) */
  double         alen;           /* front only */
  double         size;
  struct parmSeg *next;
} parmSeg;


/* Global vars */

static char   version[5] = "1.17";
static int    aflr4  = 0;
static int    wrtqud = 0;
static int    wrtuv  = 0;
static double mxang  = 0.0;
static double mxedg  = 0.0;
static double mxchd  = 0.0;
static double ggf    = 1.2;


/* front handling routines for geometric growth smoothing */

static void addSeg(int *nSeg, parmSeg *segs, parmSeg seg)
{
  parmSeg *link;
  
  if (*nSeg >= MAXSEG) {
    printf(" ERROR: No more room for Segments!\n");
    exit(EXIT_FAILURE);
  }
  segs[*nSeg] = seg;
  link = seg.prev;
  if (link != NULL) {
    if (link->next != NULL) {
      printf(" ERROR: Double hit for prev/next!\n");
      exit(EXIT_FAILURE);
    }
    link->next = &segs[*nSeg];
  }
  link = seg.next;
  if (link != NULL) {
    if (link->prev != NULL) {
      printf(" ERROR: Double hit for next/prev!\n");
      exit(EXIT_FAILURE);
    }
    link->prev = &segs[*nSeg];
  }
  *nSeg += 1;
}


static void addFront(parmSeg **frontsx, parmSeg **poolx, double params[2],
                     double alen, double size, /*@null@*/ parmSeg *connect)
{
  parmSeg *seg, *fronts, *pool, *link;
  
  fronts = *frontsx;
  pool   = *poolx;
  
  /* get a free segment */
  if (pool != NULL) {
    seg  = pool;
    pool = seg->next;
  } else {
    seg = (parmSeg *) EG_alloc(sizeof(parmSeg));
    if (seg == NULL) {
      printf(" ERROR: allocating Front Segment!\n");
      exit(EXIT_FAILURE);
    }
  }
  seg->prev     = connect;
  seg->parms[0] = params[0];
  seg->parms[1] = params[1];
  seg->alen     = alen;
  seg->size     = size;
  seg->next     = NULL;
  if (fronts == NULL) {
    fronts = seg;
  } else {
    link = fronts;
    while (link->next != NULL) link = link->next;
    link->next = seg;
  }
  
  *frontsx = fronts;
  *poolx   = pool;
}


static void delFront(parmSeg **frontsx, parmSeg **poolx, parmSeg *front)
{
  int     hit = 0;
  parmSeg *fronts, *pool, *link, *last;
  
  fronts = *frontsx;
  pool   = *poolx;
  
  if (front == NULL) {
    printf(" ERROR: NULL Front to delete!\n");
    exit(EXIT_FAILURE);
  }
  
  /* delete a segment */
  if (front == fronts) {
    fronts = front->next;
    hit++;
  } else {
    link = last = fronts;
    while (link->next != NULL) {
      last = link;
      link = link->next;
      if (link == front) {
        last->next = front->next;
        hit++;
        break;
      }
    }
  }
  if (hit == 0) {
    printf(" ERROR: Front not found to delete!\n");
    exit(EXIT_FAILURE);
  }
  
  front->next = NULL;
  if (pool == NULL) {
    pool = front;
  } else {
    link = pool;
    while (link->next != NULL) link = link->next;
    link->next = front;
  }
  
  *frontsx = fronts;
  *poolx   = pool;
}


static void freeFront(/*@null@*/ parmSeg *fronts, /*@null@*/ parmSeg *pool)
{
  parmSeg *link, *last;
  
  if (fronts != NULL) {
    link = fronts;
    while (link != NULL) {
      last = link;
      link = link->next;
      EG_free(last);
    }
  }
  
  if (pool != NULL) {
    link = pool;
    while (link != NULL) {
      last = link;
      link = link->next;
      EG_free(last);
    }
  }
}


/*
 * 	calculates and returns a complete Body tessellation
 *                 note that this is consistent with the EGADS global numbers
 *
 * 	where:	body	- ego of a body tessellation
 *              nface   - number of faces in the body
 *              face    - the EGADS Face objects
 * 		nvert	- Number of vertices (returned)
 * 		verts	- coordinates (returned) 3*nverts in len -- freeable
 * 		ntriang	- number of triangles (returned)
 * 		triang	- triangle indices (returned) 3*ntriang in len
 *			  freeable
 *              comp    - Cart3D component ID per triangle -- freeable
 */

static int
bodyTessellation(ego tess, int nface, int *nvert, double **verts,
                 int *ntriang, int **triang, int **comp)

{
  int          status, i, j, k, base, npts, ntri, *tri, *table, *compID;
  int          plen, tlen;
  const int    *tris, *tric, *ptype, *pindex;
  double       *xyzs;
  const double *points, *uv;

  *nvert  = *ntriang = 0;
  *verts  = NULL;
  *triang = NULL;
  *comp   = NULL;
  npts = ntri = 0;

  for (i = 1; i <= nface; i++) {
    status = EG_getTessFace(tess, i, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) {
      printf(" Face %d: EG_getTessFace status = %d (bodyTessellation)!\n",
             i, status);
    } else {
      npts += plen;
      ntri += tlen;
    }
  }

  /* get the memory associated with the points */

  table = (int *) EG_alloc(2*npts*sizeof(int));
  if (table == NULL) {
    printf(" Error: Can not allocate node table (bodyTessellation)!\n");
    return EGADS_MALLOC;
  }

  xyzs = (double *) EG_alloc(3*npts*sizeof(double));
  if (xyzs == NULL) {
    printf(" Error: Can not allocate XYZs (bodyTessellation)!\n");
    EG_free(table);
    return EGADS_MALLOC;
  }

  /* zipper up the edges -- a Face at a time */

  npts = 0;
  for (j = 1; j <= nface; j++) {
    status = EG_getTessFace(tess, j, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) continue;

    for (i = 0; i < plen; i++) {
      table[2*npts  ] = ptype[i];
      table[2*npts+1] = pindex[i];
      xyzs[3*npts  ]  = points[3*i  ];
      xyzs[3*npts+1]  = points[3*i+1];
      xyzs[3*npts+2]  = points[3*i+2];

      /* for non-interior pts -- try to match with others */

      if (ptype[i] != -1)
        for (k = 0; k < npts; k++)
          if ((table[2*k]==table[2*npts]) && (table[2*k+1]==table[2*npts+1])) {
            table[2*npts  ] = k;
            table[2*npts+1] = 0;
            break;
          }
      
      npts++;
    }

  }

  /* fill up the whole triangle list -- a Face at a time */

  tri = (int *) EG_alloc(3*ntri*sizeof(int));
  if (tri == NULL) {
    printf(" Error: Can not allocate triangles (bodyTessellation)!\n");
    EG_free(xyzs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  compID = (int *) EG_alloc(ntri*sizeof(int));
  if (compID == NULL) {
    printf(" Error: Can not allocate components (bodyTessellation)!\n");
    EG_free(tri);
    EG_free(xyzs);
    EG_free(table);
    return EGADS_MALLOC;
  }

  ntri = base = 0;
  for (j = 1; j <= nface; j++) {

    /* get the face tessellation and store it away */
    status = EG_getTessFace(tess, j, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) continue;

    for (i = 0; i < tlen; i++, ntri++) {

      k = tris[3*i  ] + base;

      if (table[2*k-1] == 0) {
        tri[3*ntri  ] = table[2*k-2] + 1;
      } else {
        tri[3*ntri  ] = k;
      }

      k = tris[3*i+1] + base;

      if (table[2*k-1] == 0) {
        tri[3*ntri+1] = table[2*k-2] + 1;
      } else {
        tri[3*ntri+1] = k;
      }

      k = tris[3*i+2] + base;

      if (table[2*k-1] == 0) {
        tri[3*ntri+2] = table[2*k-2] + 1;
      } else {
        tri[3*ntri+2] = k;
      }

      compID[ntri] = j;
    }
    base += plen;
  }


  /*  remove the unused points -- crunch the point list
   *  NOTE: the returned pointer verts has the full length (not realloc'ed)
   */
  for (i = 0; i <   npts; i++) table[i] = 0;
  for (i = 0; i < 3*ntri; i++) table[tri[i]-1]++;
  for (plen = i = 0; i < npts; i++) {
    if (table[i] == 0) continue;
    xyzs[3*plen  ] = xyzs[3*i  ];
    xyzs[3*plen+1] = xyzs[3*i+1];
    xyzs[3*plen+2] = xyzs[3*i+2];
    plen++;
    table[i] = plen;
  }

  /* reset the triangle indices */
  for (i = 0; i < 3*ntri; i++) {
    k      = tri[i]-1;
    tri[i] = table[k];
  }

  EG_free(table);

  *nvert   = plen;
  *verts   = xyzs;
  *ntriang = ntri;
  *triang  = tri;
  *comp    = compID;

  return EGADS_SUCCESS;

}


static int strroot(char *s)
{

  /* Return string index location in string s by locating last . */

  int i, iroot, len;

  len   = strlen(s);
  iroot = len-1;

  i = len-1;
  while ( s[i] != '.' && i > 0 ) {
    i--;
  }
  if ( i >= 0 ) iroot = i;

  return iroot;
}


static void writeAttr(FILE *fp, ego obj, /*@null@*/ char *filter)
{
  int          i, j, n, stat, nattr, atype, alen, flen;
  const int    *ints;
  const double *reals;
  const char   *name, *str;
  
  stat = EG_attributeNum(obj, &nattr);
  if (stat != EGADS_SUCCESS) return;
  
  for (n = i = 0; i < nattr; i++) {
    stat = EG_attributeGet(obj, i+1, &name, &atype, &alen, &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    if ((atype != ATTRINT) && (atype != ATTRREAL) && (atype != ATTRSTRING))
      continue;
    if (filter != NULL) {
      flen = strlen(filter);
      if (flen > strlen(name)) continue;
      for (j = 0; j < flen; j++)
        if (filter[j] != name[j]) break;
      if (j != flen) continue;
    }
    n++;
  }
  fprintf(fp, " %6d\n", n);
  if (n == 0) return;

  for (i = 0; i < nattr; i++) {
    stat = EG_attributeGet(obj, i+1, &name, &atype, &alen, &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    if ((atype != ATTRINT) && (atype != ATTRREAL) && (atype != ATTRSTRING))
      continue;
    if (filter != NULL) {
      flen = strlen(filter);
      if (flen > strlen(name)) continue;
      for (j = 0; j < flen; j++)
        if (filter[j] != name[j]) break;
      if (j != flen) continue;
    }
    if (atype == ATTRSTRING) alen = strlen(str);
    fprintf(fp, " %6d %6d %s\n", atype, alen, name);
    if (atype == ATTRSTRING) {
      fprintf(fp, " %s\n", str);
    } else if (atype == ATTRREAL) {
      for (j = 0; j < alen; j++) {
        fprintf(fp, " %20.13le", reals[j]);
        if ((j+1)%4 == 0) fprintf(fp,"\n");
      }
      if (j%4 != 0) fprintf(fp, "\n");
    } else {
      for (j = 0; j < alen; j++) {
        fprintf(fp, " %10d", ints[j]);
        if ((j+1)%8 == 0) fprintf(fp,"\n");
      }
      if (j%8 != 0) fprintf(fp, "\n");
    }
  }

}


static void getNodeSpacing(ego tess, int nedge, ego *edges, double *spacing)
{
  int          i, j, m, n, nvert, oclass, mtype, status, *senses;
  double       dist, range[2];
  const double *xyzs, *ts;
  ego          body, geom, *objs;
  
  /* compute the smallest segment spacing touching a Node from the Edge
     tessellations */
   
  status = EG_statusTessBody(tess, &body, &n, &nvert);
  if (status != EGADS_SUCCESS) {
    printf(" EG_statusTessBody = %d (getNodeSpacing)!\n", status);
    return;
  }
  
  for (i = 0; i < nedge; i++) {
    status = EG_getTopology(edges[i], &geom, &oclass, &mtype, range, &n, &objs,
                            &senses);
    if (status != EGADS_SUCCESS) {
      printf(" %d: EG_getTopology = %d (getNodeSpacing)!\n", i+1, status);
      continue;
    }
    if (mtype == DEGENERATE) continue;
    
    status = EG_getTessEdge(tess, i+1, &m, &xyzs, &ts);
    if (status != EGADS_SUCCESS) {
      printf(" %d: EG_getTessEdge = %d (getNodeSpacing)!\n", i+1, status);
      continue;
    }
    if (m == 0) continue;
    j = EG_indexBodyTopo(body, objs[0]);
    if (j <= EGADS_SUCCESS) {
      printf(" %d: EG_indexBodyTopo 0 = %d (getNodeSpacing)!\n", i+1, status);
      continue;
    }
    dist = sqrt((xyzs[0]-xyzs[3])*(xyzs[0]-xyzs[3]) +
                (xyzs[1]-xyzs[4])*(xyzs[1]-xyzs[4]) +
                (xyzs[2]-xyzs[5])*(xyzs[2]-xyzs[5]));
    if (spacing[j-1] == 0.0) {
      spacing[j-1] = dist;
    } else {
      if (dist < spacing[j-1]) spacing[j-1] = dist;
    }
    
    if (mtype == TWONODE) j = EG_indexBodyTopo(body, objs[1]);
    if (j <= EGADS_SUCCESS) {
      printf(" %d: EG_indexBodyTopo 1 = %d (getNodeSpacing)!\n", i+1, status);
      continue;
    }
    dist = sqrt((xyzs[3*m-3]-xyzs[3*m-6])*(xyzs[3*m-3]-xyzs[3*m-6]) +
                (xyzs[3*m-2]-xyzs[3*m-5])*(xyzs[3*m-2]-xyzs[3*m-5]) +
                (xyzs[3*m-1]-xyzs[3*m-4])*(xyzs[3*m-1]-xyzs[3*m-4]));
    if (spacing[j-1] == 0.0) {
      spacing[j-1] = dist;
    } else {
      if (dist < spacing[j-1]) spacing[j-1] = dist;
    }
  }
}


static void fillIn(double parm, double delta, int *n, double *parms)
{
  int i, ii, hit;
  
  /* adjust the parameter sequence by inserting Node positions and spacing
     before, after or both */
  
  for (i = 0; i < *n; i++)
    if (parms[i] > parm) break;
  
  ii = i-1;
  if (ii == -1) ii = 0;
  if (i  == *n) i  = *n-1;
#ifdef DEBUG
  printf("         %lf %lf %lf", parms[ii], parm, parms[i]);
#endif
  /* do we insert the value? */
  if ((fabs(parm-parms[ii]) <= delta) || (fabs(parm-parms[i]) <= delta)) {
    hit = i;
    if (fabs(parm-parms[ii]) <= delta) hit = ii;
#ifdef DEBUG
    printf(" Node @ %d", hit);
#endif
  } else {
    if (*n == MAXUVS) return;
    hit = i;
    for (i = *n; i >= hit; i--) parms[i] = parms[i-1];
    parms[hit] = parm;
    *n += 1;
#ifdef DEBUG
    printf(" Add @ %d", hit);
#endif
  }
  
  /* look forward */
  if (hit != *n-1)
    if (parms[hit+1]-parms[hit] > 2.0*delta) {
      if (*n == MAXUVS) return;
      for (i = *n; i >= hit+1; i--) parms[i] = parms[i-1];
      parms[hit+1] = parms[hit] + delta;
      *n += 1;
#ifdef DEBUG
      printf(" Add+");
#endif
    }
  
  /* look back */
  if (hit != 0)
    if (parms[hit]-parms[hit-1] > 2.0*delta) {
      if (*n == MAXUVS) return;
      for (i = *n; i >= hit+1; i--) parms[i] = parms[i-1];
      parms[hit] -= delta;
      *n += 1;
#ifdef DEBUG
      printf(" Add-");
#endif
    }
#ifdef DEBUG
  printf("\n");
#endif
}


static void insertNodeSpacing(ego tess, ego face, int iface, double *spacing,
                              UVmap *map)
{
  int          i, m, n, status;
  double       udist, vdist, result[18], *u1, *v1;
  const int    *pindex, *ptype, *tris, *tric;
  const double *xyzs, *uvs;
  
  u1 = &result[3];
  v1 = &result[6];
  
  status = EG_getTessFace(tess, iface, &m, &xyzs, &uvs, &ptype, &pindex,
                          &n, &tris, &tric);
  if (status != EGADS_SUCCESS) {
    printf(" %d: EG_getTessFace = %d (insertNodeSpacing)!\n", iface, status);
    return;
  }
  for (i = 0; i < m; i++) {
    if (ptype[i] != 0) continue;
#ifdef DEBUG
    printf(" Face %d: Node = %d, spacing = %lf   uv = %lf %lf\n",
           iface, pindex[i], spacing[pindex[i]-1], uvs[2*i  ], uvs[2*i+1]);
#endif
    status = EG_evaluate(face, &uvs[2*i], result);
    if (status != EGADS_SUCCESS) {
      printf(" %d: EG_evaluate = %d (insertNodeSpacing)!\n", iface, status);
      return;
    }
    udist = sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]);
    if (udist != 0.0) {
      udist = spacing[pindex[i]-1]/udist;
      fillIn(uvs[2*i  ], udist, &map->nu, map->us);
      if (map->nu == MAXUVS) break;
    }
    vdist = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    if (vdist != 0.0) {
      vdist = spacing[pindex[i]-1]/vdist;
      fillIn(uvs[2*i+1], vdist, &map->nv, map->vs);
    }
#ifdef DEBUG
    printf("         udelta = %lf   vdelta = %lf\n", udist, vdist);
#endif
  }
}


static void dumpFront(int nSeg, parmSeg *segs, parmSeg *front)
{
  int     i;
  parmSeg *link;
  
  link = front;
  while (link != NULL) {
    printf(" F %lx  parms = %lf %lf,  len = %lf\n    size = %lf  connect = %lx\n",
           (ULONG) link, link->parms[0], link->parms[1], link->alen, link->size,
           (ULONG) link->prev);
    link = link->next;
  }

  for (i = 0; i < nSeg; i++) {
    link = &segs[i];
    printf(" S %lx  parms = %lf %lf,  size = %lf\n    prev = %lx  next = %lx\n",
           (ULONG) link, link->parms[0], link->parms[1], link->size,
           (ULONG) link->prev, (ULONG) link->next);
  }
}


static void patchSeg(int nSeg, parmSeg *segs, parmSeg **frontx,
                     parmSeg **poolx, parmSeg *fill)
{
  int     i;
  parmSeg *front, *pool, *last, *link;
  
  /* look in our current segments */
  for (i = 0; i < nSeg; i++) {
    if (fabs(fill->parms[0] - segs[i].parms[1]) < FUZZ) fill->prev = &segs[i];
    if (fabs(fill->parms[1] - segs[i].parms[0]) < FUZZ) fill->next = &segs[i];
  }
  
  front = *frontx;
  pool  = *poolx;

  link = front;
  while (link != NULL) {
    last = link;
    link = link->next;
    if (last->parms[1] < 0) {
      if (fabs(fill->parms[1] - last->parms[0]) < FUZZ)
        delFront(&front, &pool, last);
    } else {
      if (fabs(fill->parms[0] - last->parms[0]) < FUZZ)
        delFront(&front, &pool, last);
    }
  }
  
  *frontx = front;
  *poolx  = pool;
}


static void smoothParm(int *np, double *parms, double *r, double *q, double fact)
{
  int     i, j, k, m, n, mark[MAXUVS], nSeg = 0;
  double  smals, ave, alen, dist, qi, pi, size, fra, params[2], al[MAXUVS];
  parmSeg segs[MAXSEG], seg, *link, *smallest, *other;
  parmSeg *first = NULL, *front = NULL, *pool = NULL;
  
  n = *np;
  if (n <= 2) return;
  al[0] = 0.0;
  for (i = 0; i < n-1; i++) {
    mark[i] = 1;               /* available */
    al[i+1] = al[i] + r[i];    /* compute segment arclength */
  }
  
  /* get the smallest & average segment */
  smals = ave = r[0];
  j     = 0;
  for (i = 1; i < n-1; i++) {
    ave += r[i];
    if (r[i] < smals) {
      j     = i;
      smals = r[i];
    }
  }
  ave /= (n-1);
  if (fabs(smals-ave)/ave < 0.1) return;
  
  /* add the first segment */
  seg.prev     = NULL;
  seg.parms[0] = parms[j];
  seg.parms[1] = parms[j+1];
  seg.alen     = 0.0;
  seg.size     = r[j];
  seg.next     = NULL;
  addSeg(&nSeg, segs, seg);
  mark[j]      = -1;   /* used */

  /* initialize the front */
  if (j == 0) {
    params[0] = parms[j+1];
    params[1] = +1.0;
    addFront(&front, &pool, params, al[1],   r[0],   &segs[0]);
    first     = &segs[0];
    params[0] = parms[n-1];
    params[1] = -1.0;
    addFront(&front, &pool, params, al[n-1], r[n-2], NULL);
  } else if (j == n-2) {
    params[0] = parms[0];
    params[1] = +1.0;
    addFront(&front, &pool, params, al[0],   r[0],   NULL);
    params[0] = parms[n-2];
    params[1] = -1.0;
    addFront(&front, &pool, params, al[j],   r[j],   &segs[0]);
  } else {
    params[0] = parms[0];
    params[1] = +1.0;
    addFront(&front, &pool, params, al[0],   r[0],   NULL);
    params[0] = parms[j];
    params[1] = -1.0;
    addFront(&front, &pool, params, al[j],   r[j],   &segs[0]);
    params[0] = parms[j+1];
    params[1] = +1.0;
    addFront(&front, &pool, params, al[j+1], r[j],   &segs[0]);
    params[0] = parms[n-1];
    params[1] = -1.0;
    addFront(&front, &pool, params, al[n-1], r[n-2], NULL);
  }
  
  do {
    smals *= fact;
    /* first look at available original segments */
    j      = -1;
    for (i = 0; i < n-1; i++)
      if (mark[i] == 1)
        if (r[i] < smals)
          if (j == -1) {
            ave = r[i];
            j   = i;
          } else {
            if (ave > r[i]) {
              ave = r[i];
              j   = i;
            }
          }
#ifdef DEBUG
    printf(" next = %d  ave = %lf  small = %lf  nSeg = %d\n",
           j, ave, smals, nSeg);
#endif
    
    if (j != -1) {
      /* add the next segment */
      seg.prev     = NULL;
      seg.parms[0] = parms[j];
      seg.parms[1] = parms[j+1];
      seg.alen     = 0.0;
      seg.size     = r[j];
      seg.next     = NULL;
      patchSeg(nSeg, segs, &front, &pool, &seg);
      addSeg(&nSeg, segs, seg);
      if (seg.parms[0] == parms[0]) first = &segs[nSeg-1];
      mark[j]      = -1;   /* used */
      smals        = r[j];
      if ((seg.prev == NULL) && (j != 0)) {
        params[0] = parms[j];
        params[1] = -1.0;
        addFront(&front, &pool, params, al[j],   r[j], &segs[nSeg-1]);
      }
      if ((seg.next == NULL) && (j != n-2)) {
        params[0] = parms[j+1];
        params[1] = +1.0;
        addFront(&front, &pool, params, al[j+1], r[j], &segs[nSeg-1]);
      }
      continue;
    }
    if (front == NULL) continue;
      
    /* now find the smallest front segment */
    smallest = link = front;
    while (link != NULL) {
      if (smallest->size > link->size) smallest = link;
      link = link->next;
    }
  
    /* look ahead for fronts colliding */
    link  = front;
    ave   = al[n-1];
    other = NULL;
    while (link != NULL) {
      if (smallest->parms[1]*link->parms[1] < 0.0) {
        if (smallest->parms[1] < 0) {
          if (link->parms[0] < smallest->parms[0]) {
            alen = fabs(smallest->alen - link->alen);
            if (alen < ave) {
              other = link;
              ave   = alen;
            }
          }
        } else {
          if (link->parms[0] > smallest->parms[0]) {
            alen = fabs(smallest->alen - link->alen);
            if (alen < ave) {
              other = link;
              ave   = alen;
            }
          }
        }
      }
      link = link->next;
    }
    if (other != NULL) {
      dist = 0.5*(other->size + smallest->size);
      m    = ave/dist + .49;
      alen = smallest->size/other->size;
      if (alen < 1.0) alen = 1.0/alen;
#ifdef DEBUG
      printf("   alen = %lf  delta = %lf   divs = %d   ratio = %lf\n",
             ave, dist, m, alen);
#endif
      if (m == 0) m = 1;
      if (alen > 10.0) m = 0;
      if ((m > 0) && (m < 10)) {          /* 10 is arbitrary */
        if (smallest->parms[0] > other->parms[0]) {
          link     = smallest;
          smallest = other;
          other    = link;
        }
        
        /* mark original segments */
        alen = smallest->alen;
        for (j = 0; j < n-1; j++)
          if ((alen >= al[j]) && (alen <= al[j+1])) break;
        if (j == n-1) {
          printf(" ERROR: Cannot Interpolate 0  %lf  %lf -- collide!\n",
                 alen, al[n-1]);
          exit(EXIT_FAILURE);
        }
        dist = (alen - al[j])/(al[j+1] - al[j]);
        if ((fabs(dist-1.0) < FUZZ) && (mark[j] == -1) && (j != n-2)) {
          j++;
          dist = 0.0;
        }
        alen = other->alen;
        for (k = 0; k < n-1; k++)
          if ((alen >= al[k]) && (alen <= al[k+1])) break;
        if (k == n-1) {
          if (fabs(alen-al[n-1]) > FUZZ) {
            printf(" ERROR: Cannot Interpolate 1  %lf  %lf -- collide!\n",
                   alen, al[n-1]);
            exit(EXIT_FAILURE);
          } else {
            k--;
          }
        }
        dist = (alen - al[k])/(al[k+1] - al[k]);
        if ((fabs(dist) < FUZZ) && (mark[k] == -1) && (k != 0)) {
          k--;
          dist = 1.0;
        }
        for (i = j; i <= k; i++) {
          if (mark[i] == -1) {
            printf(" ERROR: Used Segment %d -- %d %d -- collide!\n", i, j, k);
            exit(EXIT_FAILURE);
          }
          mark[i] = 0;
        }
        
        /* make new segments -- scale parameter by arcLength */
        alen = 0.0;
        qi   = smallest->size;
        for (i = 0; i < m; i++) {
          ave   = i+1;
          ave  /= m;
          size  = (1.0-ave)*smallest->size + ave*other->size;
          alen += 0.5*(qi + size);
          qi    = size;
        }
        qi       = smallest->size;
        pi       = smallest->parms[0];
        seg.alen = dist = 0.0;
        for (i = 0; i < m; i++) {
          ave          = i+1;
          ave         /= m;
          size         = (1.0-ave)*smallest->size + ave*other->size;
          dist        += 0.5*(qi + size);
          fra          = dist/alen;
          seg.size     = 0.5*(qi + size);
          seg.parms[0] = pi;
          seg.parms[1] = (1.0-fra)*smallest->parms[0] + fra*other->parms[0];
          seg.prev     = NULL;
          seg.next     = NULL;
          qi           = size;
          pi           = seg.parms[1];
          patchSeg(nSeg, segs, &front, &pool, &seg);
          addSeg(&nSeg, segs, seg);
          if (seg.parms[0] == parms[0]) first = &segs[nSeg-1];
        }
        smals /= fact;
        continue;
      }
    }
    
    /* add a single segment from the front's smallest */
    smals = fact*smallest->size;
    alen  = smallest->parms[1]*smals + smallest->alen;
    for (j = 0; j < n-1; j++)
      if ((alen >= al[j]) && (alen <= al[j+1])) break;
    if (j == n-1) {
      printf(" ERROR: Cannot Interpolate  alen = %lf [%lf %lf]  %lx!\n",
             alen, al[0], al[n-1], (ULONG) smallest);
      dumpFront(nSeg, segs, front);
      exit(EXIT_FAILURE);
    }
    dist = (alen - al[j])/(al[j+1] - al[j]);
    if (mark[j] == -1) {
      printf(" ERROR: Used Segment -- dist = %lf!\n", dist);
      exit(EXIT_FAILURE);
    }
    mark[j]   = 0;
    qi        = q[j] + dist*(q[j+1] - q[j]);
    pi        = smals/qi;
    seg.prev  = NULL;
    seg.alen  = 0.0;
    seg.size  = smals;
    seg.next  = NULL;
    params[1] = smallest->parms[1];
    if (params[1] > 0) {
      seg.parms[0] =             smallest->parms[0];
      seg.parms[1] = params[0] = smallest->parms[0] + pi;
    } else {
      seg.parms[0] = params[0] = smallest->parms[0] - pi;
      seg.parms[1] =             smallest->parms[0];
    }
    patchSeg(nSeg, segs, &front, &pool, &seg);
    addSeg(&nSeg, segs, seg);
    if (seg.parms[0] == parms[0]) first = &segs[nSeg-1];
    addFront(&front, &pool, params, alen, smals, &segs[nSeg-1]);
    
  } while (front != NULL);
  
  /* reset our parameter sequence from the linked list */
  if (first == NULL) {
    printf(" ERROR: No First Segment!\n");
    dumpFront(nSeg, segs, front);
    exit(EXIT_FAILURE);
  }
  link = first;
  i    = 0;
  pi   = parms[n-1];
  while (link != NULL) {
    if (i == 0) parms[0] = link->parms[0];
    parms[i+1] = link->parms[1];
    i++;
    link = link->next;
  }
  if ((i != nSeg) || (parms[i] != pi)) {
    printf(" ERROR: Finialization -- %d %d  %lf %lf!\n", i, nSeg, parms[i], pi);
    dumpFront(nSeg, segs, front);
    exit(EXIT_FAILURE);
  }
  *np = nSeg+1;
  
  freeFront(front, pool);
}


static void smoothMap(ego face, double factor, UVmap *map)
{
  int     i, j, stat;
  double  uv[2], *u1, *v1, r[MAXUVS], q[MAXUVS], result[18];
  double  dp[3*MAXUVS], xyzs[3*MAXUVS];
  
  if (factor < 1.0) return;
  u1 = &result[3];
  v1 = &result[6];
  
  /* look in the U direction */
  for (i = 0; i < map->nu; i++) r[i] = q[i] = 0.0;
  for (j = 0; j < map->nv; j++) {
    uv[1] = map->vs[j];
    for (i = 0; i < map->nu; i++) {
      uv[0] = map->us[i];
      stat  = EG_evaluate(face, uv, result);
      if (stat != EGADS_SUCCESS) {
        printf(" smoothMap: Fill U EG_evaluate = %d!\n", stat);
        return;
      }
      xyzs[3*i  ] = result[0];
      xyzs[3*i+1] = result[1];
      xyzs[3*i+2] = result[2];
      dp[3*i  ]   = u1[0];
      dp[3*i+1]   = u1[1];
      dp[3*i+2]   = u1[2];
    }
    q[0] += sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    for (i = 1; i < map->nu; i++) {
      q[i  ] += sqrt(dp[3*i  ]*dp[3*i  ] + dp[3*i+1]*dp[3*i+1] +
                     dp[3*i+2]*dp[3*i+2]);
      r[i-1] += sqrt((xyzs[3*i-3]-xyzs[3*i  ])*(xyzs[3*i-3]-xyzs[3*i  ]) +
                     (xyzs[3*i-2]-xyzs[3*i+1])*(xyzs[3*i-2]-xyzs[3*i+1]) +
                     (xyzs[3*i-1]-xyzs[3*i+2])*(xyzs[3*i-1]-xyzs[3*i+2]));
    }
  }
  q[0] /= map->nv;
  for (i = 0; i < map->nu-1; i++) {
    q[i+1] /= map->nv;
    r[i]   /= map->nv;
  }
#ifdef DEBUG
  printf(" nUs = %d\n %lf             %lf\n", map->nu, map->us[0], q[0]);
  for (i = 1; i < map->nu; i++)
    printf(" %lf  %lf  %lf\n", map->us[i], r[i-1], q[i]);
#endif
  smoothParm(&map->nu, map->us, r, q, factor);
#ifdef DEBUG
  for (i = 0; i < map->nu; i++) r[i] = q[i] = 0.0;
  for (j = 0; j < map->nv; j++) {
    uv[1] = map->vs[j];
    for (i = 0; i < map->nu; i++) {
      uv[0] = map->us[i];
      stat  = EG_evaluate(face, uv, result);
      if (stat != EGADS_SUCCESS) {
        printf(" smoothMap: Fill U EG_evaluate = %d!\n", stat);
        return;
      }
      xyzs[3*i  ] = result[0];
      xyzs[3*i+1] = result[1];
      xyzs[3*i+2] = result[2];
      dp[3*i  ]   = u1[0];
      dp[3*i+1]   = u1[1];
      dp[3*i+2]   = u1[2];
    }
    q[0] += sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    for (i = 1; i < map->nu; i++) {
      q[i  ] += sqrt(dp[3*i  ]*dp[3*i  ] + dp[3*i+1]*dp[3*i+1] +
                     dp[3*i+2]*dp[3*i+2]);
      r[i-1] += sqrt((xyzs[3*i-3]-xyzs[3*i  ])*(xyzs[3*i-3]-xyzs[3*i  ]) +
                     (xyzs[3*i-2]-xyzs[3*i+1])*(xyzs[3*i-2]-xyzs[3*i+1]) +
                     (xyzs[3*i-1]-xyzs[3*i+2])*(xyzs[3*i-1]-xyzs[3*i+2]));
    }
  }
  q[0] /= map->nv;
  for (i = 0; i < map->nu-1; i++) {
    q[i+1] /= map->nv;
    r[i]   /= map->nv;
  }
  printf(" nUs = %d\n %lf             %lf\n", map->nu, map->us[0], q[0]);
  r[map->nu-1] = r[map->nu-2];
  for (i = 1; i < map->nu; i++)
    printf(" %lf  %lf  %lf   %lf\n", map->us[i], r[i-1], q[i], r[i]/r[i-1]);
#endif
  
  /* look in the V direction */
  for (j = 0; j < map->nv; j++) r[j] = q[j] = 0.0;
  for (i = 0; i < map->nu; i++) {
    uv[0] = map->us[i];
    for (j = 0; j < map->nv; j++) {
      uv[1] = map->vs[j];
      stat  = EG_evaluate(face, uv, result);
      if (stat != EGADS_SUCCESS) {
        printf(" smoothMap: Fill V EG_evaluate = %d!\n", stat);
        return;
      }
      xyzs[3*j  ] = result[0];
      xyzs[3*j+1] = result[1];
      xyzs[3*j+2] = result[2];
      dp[3*j  ]   = v1[0];
      dp[3*j+1]   = v1[1];
      dp[3*j+2]   = v1[2];
    }
    q[0] += sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    for (j = 1; j < map->nv; j++) {
      q[j  ] += sqrt(dp[3*j  ]*dp[3*j  ] + dp[3*j+1]*dp[3*j+1] +
                     dp[3*j+2]*dp[3*j+2]);
      r[j-1] += sqrt((xyzs[3*j-3]-xyzs[3*j  ])*(xyzs[3*j-3]-xyzs[3*j  ]) +
                     (xyzs[3*j-2]-xyzs[3*j+1])*(xyzs[3*j-2]-xyzs[3*j+1]) +
                     (xyzs[3*j-1]-xyzs[3*j+2])*(xyzs[3*j-1]-xyzs[3*j+2]));
    }
  }
  q[0] /= map->nu;
  for (j = 0; j < map->nv-1; j++) {
    q[j+1] /= map->nu;
    r[j]   /= map->nu;
  }
#ifdef DEBUG
  printf(" nVs = %d\n %lf             %lf\n", map->nv, map->vs[0], q[0]);
  for (j = 1; j < map->nv; j++)
    printf(" %lf  %lf  %lf\n", map->vs[j], r[j-1], q[j]);
#endif
  smoothParm(&map->nv, map->vs, r, q, factor);
#ifdef DEBUG
  for (j = 0; j < map->nv; j++) r[j] = q[j] = 0.0;
  for (i = 0; i < map->nu; i++) {
    uv[0] = map->us[i];
    for (j = 0; j < map->nv; j++) {
      uv[1] = map->vs[j];
      stat  = EG_evaluate(face, uv, result);
      if (stat != EGADS_SUCCESS) {
        printf(" smoothMap: Fill V EG_evaluate = %d!\n", stat);
        return;
      }
      xyzs[3*j  ] = result[0];
      xyzs[3*j+1] = result[1];
      xyzs[3*j+2] = result[2];
      dp[3*j  ]   = v1[0];
      dp[3*j+1]   = v1[1];
      dp[3*j+2]   = v1[2];
    }
    q[0] += sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    for (j = 1; j < map->nv; j++) {
      q[j  ] += sqrt(dp[3*j  ]*dp[3*j  ] + dp[3*j+1]*dp[3*j+1] +
                     dp[3*j+2]*dp[3*j+2]);
      r[j-1] += sqrt((xyzs[3*j-3]-xyzs[3*j  ])*(xyzs[3*j-3]-xyzs[3*j  ]) +
                     (xyzs[3*j-2]-xyzs[3*j+1])*(xyzs[3*j-2]-xyzs[3*j+1]) +
                     (xyzs[3*j-1]-xyzs[3*j+2])*(xyzs[3*j-1]-xyzs[3*j+2]));
    }
  }
  q[0] /= map->nu;
  for (j = 0; j < map->nv-1; j++) {
    q[j+1] /= map->nu;
    r[j]   /= map->nu;
  }
  printf(" nVs = %d\n %lf             %lf\n", map->nv, map->vs[0], q[0]);
  r[map->nv-1] = r[map->nv-2];
  for (j = 1; j < map->nv; j++)
    printf(" %lf  %lf  %lf   %lf\n", map->vs[j], r[j-1], q[j], r[j]/r[j-1]);
#endif
}


static void updateMap(ego face, double mxedg, double sag, double angle, UVmap *map)
{
  int    i, j, max, stat, cnt;
  double d, dist, dot, last, uv[2], xyz[3], *u1, *v1, result[18], xyzs[3*MAXUVS];
  
  dot = cos(PI*angle/180.0);
  u1  = &result[3];
  v1  = &result[6];
  
  /* sag -- look in U direction */
  do {
    max  = -1;
    dist = 0.0;
    for (j = 0; j < map->nv; j++) {
      uv[1] = map->vs[j];
      for (i = 0; i < map->nu; i++) {
        uv[0] = map->us[i];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill U EG_evaluate = %d!\n", stat);
          return;
        }
        xyzs[3*i  ] = result[0];
        xyzs[3*i+1] = result[1];
        xyzs[3*i+2] = result[2];
      }
      for (i = 1; i < map->nu; i++) {
        xyz[0] = 0.5*(xyzs[3*i-3]  + xyzs[3*i  ]);
        xyz[1] = 0.5*(xyzs[3*i-2]  + xyzs[3*i+1]);
        xyz[2] = 0.5*(xyzs[3*i-1]  + xyzs[3*i+2]);
        uv[0]  = 0.5*(map->us[i-1] + map->us[i]);
        stat   = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Half U EG_evaluate = %d!\n", stat);
          return;
        }
        d = sqrt((xyz[0]-result[0])*(xyz[0]-result[0]) +
                 (xyz[1]-result[1])*(xyz[1]-result[1]) +
                 (xyz[2]-result[2])*(xyz[2]-result[2]));
        if (d <= sag) continue;
        if (d > dist) {
          max  = i;
          dist = d;
        }
      }
    }
    if (map->nu == MAXUVS) break;
    if (max != -1) {
      for (i = map->nu; i >= max; i--) map->us[i] = map->us[i-1];
      map->us[max] = 0.5*(map->us[max+1] + map->us[max-1]);
      map->nu++;
    }
  } while (max != -1);
  
  /* sag -- look in V direction */
  do {
    max  = -1;
    dist = 0.0;
    for (i = 0; i < map->nu; i++) {
      uv[0] = map->us[i];
      for (j = 0; j < map->nv; j++) {
        uv[1] = map->vs[j];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill V EG_evaluate = %d!\n", stat);
          return;
        }
        xyzs[3*j  ] = result[0];
        xyzs[3*j+1] = result[1];
        xyzs[3*j+2] = result[2];
      }
      for (j = 1; j < map->nv; j++) {
        xyz[0] = 0.5*(xyzs[3*j-3]  + xyzs[3*j  ]);
        xyz[1] = 0.5*(xyzs[3*j-2]  + xyzs[3*j+1]);
        xyz[2] = 0.5*(xyzs[3*j-1]  + xyzs[3*j+2]);
        uv[1]  = 0.5*(map->vs[j-1] + map->vs[j]);
        stat   = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Half V EG_evaluate = %d!\n", stat);
          return;
        }
        d = sqrt((xyz[0]-result[0])*(xyz[0]-result[0]) +
                 (xyz[1]-result[1])*(xyz[1]-result[1]) +
                 (xyz[2]-result[2])*(xyz[2]-result[2]));
        if (d <= sag) continue;
        if (d > dist) {
          max  = j;
          dist = d;
        }
      }
    }
    if (map->nv == MAXUVS) break;
    if (max != -1) {
      for (j = map->nv; j >= max; j--) map->vs[j] = map->vs[j-1];
      map->vs[max] = 0.5*(map->vs[max+1] + map->vs[max-1]);
      map->nv++;
    }
  } while (max != -1);
  
  /* angle -- look in U direction */
  last = -1.0;
  cnt  = 0;
  do {
    max  = -1;
    dist = 1.0;
    for (j = 0; j < map->nv; j++) {
      uv[1] = map->vs[j];
      for (i = 0; i < map->nu; i++) {
        uv[0] = map->us[i];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill Ua EG_evaluate = %d!\n", stat);
          return;
        }
        CROSS(xyz, u1, v1);
        d = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
        if (d == 0) d = 1.0;
        xyzs[3*i  ] = xyz[0]/d;
        xyzs[3*i+1] = xyz[1]/d;
        xyzs[3*i+2] = xyz[2]/d;
      }
      for (i = 1; i < map->nu; i++) {
        if (sqrt(xyzs[3*i-3]*xyzs[3*i-3] + xyzs[3*i-2]*xyzs[3*i-2] +
                 xyzs[3*i-1]*xyzs[3*i-1]) == 0.0) continue;
        if (sqrt(xyzs[3*i  ]*xyzs[3*i  ] + xyzs[3*i+1]*xyzs[3*i+1] +
                 xyzs[3*i+2]*xyzs[3*i+2]) == 0.0) continue;
        d = xyzs[3*i  ]*xyzs[3*i-3] + xyzs[3*i+1]*xyzs[3*i-2] +
            xyzs[3*i+2]*xyzs[3*i-1];
        if (d >= dot) continue;
        if (d < -0.2) continue;     /* there is some problem - C0? */
        if (d < dist) {
          max  = i;
          dist = d;
        }
      }
    }
    if (map->nu == MAXUVS) break;
    if (max != -1) {
/*    printf(" Umax = %d (%d)  %lf (%lf)\n", max, map->nu, dist, dot);  */
      for (i = map->nu; i >= max; i--) map->us[i] = map->us[i-1];
      map->us[max] = 0.5*(map->us[max+1] + map->us[max-1]);
      map->nu++;
    }
    if (dist >= last) cnt++;
    last = dist;
  } while ((max != -1) && (cnt < 10)) ;
  
  /* angle -- look in V direction */
  last = -1.0;
  cnt  = 0;
  do {
    max  = -1;
    dist = 1.0;
    for (i = 0; i < map->nu; i++) {
      uv[0] = map->us[i];
      for (j = 0; j < map->nv; j++) {
        uv[1] = map->vs[j];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill Va EG_evaluate = %d!\n", stat);
          return;
        }
        CROSS(xyz, u1, v1);
        d = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
        if (d == 0) d = 1.0;
        xyzs[3*j  ] = xyz[0]/d;
        xyzs[3*j+1] = xyz[1]/d;
        xyzs[3*j+2] = xyz[2]/d;
      }
      for (j = 1; j < map->nv; j++) {
        if (sqrt(xyzs[3*j-3]*xyzs[3*j-3] + xyzs[3*j-2]*xyzs[3*j-2] +
                 xyzs[3*j-1]*xyzs[3*j-1]) == 0.0) continue;
        if (sqrt(xyzs[3*j  ]*xyzs[3*j  ] + xyzs[3*j+1]*xyzs[3*j+1] +
                 xyzs[3*j+2]*xyzs[3*j+2]) == 0.0) continue;
        d = xyzs[3*j  ]*xyzs[3*j-3] + xyzs[3*j+1]*xyzs[3*j-2] +
            xyzs[3*j+2]*xyzs[3*j-1];
        if (d >= dot) continue;
        if (d < -0.2) continue;     /* there is some problem - C0? */
        if (d < dist) {
          max  = j;
          dist = d;
        }
      }
    }
    if (map->nv == MAXUVS) break;
    if (max != -1) {
/*    printf(" Vmax = %d (%d)  %lf (%lf)\n", max, map->nv, dist, dot);  */
      for (j = map->nv; j >= max; j--) map->vs[j] = map->vs[j-1];
      map->vs[max] = 0.5*(map->vs[max+1] + map->vs[max-1]);
      map->nv++;
    }
    if (dist >= last) cnt++;
    last = dist;
  } while ((max != -1) && (cnt < 10)) ;
  
  /* maxedge -- look in U direction */
  do {
    max  = -1;
    dist = 0.0;
    for (j = 0; j < map->nv; j++) {
      uv[1] = map->vs[j];
      for (i = 0; i < map->nu; i++) {
        uv[0] = map->us[i];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill U EG_evaluate = %d!\n", stat);
          return;
        }
        xyzs[3*i  ] = result[0];
        xyzs[3*i+1] = result[1];
        xyzs[3*i+2] = result[2];
      }
      for (i = 1; i < map->nu; i++) {
        xyz[0] = xyzs[3*i-3] - xyzs[3*i  ];
        xyz[1] = xyzs[3*i-2] - xyzs[3*i+1];
        xyz[2] = xyzs[3*i-1] - xyzs[3*i+2];
        d = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
        if (d <= mxedg) continue;
        if (d > dist) {
          max  = i;
          dist = d;
        }
      }
    }
    if (map->nu == MAXUVS) break;
    if (max != -1) {
      for (i = map->nu; i >= max; i--) map->us[i] = map->us[i-1];
      map->us[max] = 0.5*(map->us[max+1] + map->us[max-1]);
      map->nu++;
    }
  } while (max != -1);
  
  /* sag -- look in V direction */
  do {
    max  = -1;
    dist = 0.0;
    for (i = 0; i < map->nu; i++) {
      uv[0] = map->us[i];
      for (j = 0; j < map->nv; j++) {
        uv[1] = map->vs[j];
        stat  = EG_evaluate(face, uv, result);
        if (stat != EGADS_SUCCESS) {
          printf(" updateMap: Fill V EG_evaluate = %d!\n", stat);
          return;
        }
        xyzs[3*j  ] = result[0];
        xyzs[3*j+1] = result[1];
        xyzs[3*j+2] = result[2];
      }
      for (j = 1; j < map->nv; j++) {
        xyz[0] = xyzs[3*j-3] - xyzs[3*j  ];
        xyz[1] = xyzs[3*j-2] - xyzs[3*j+1];
        xyz[2] = xyzs[3*j-1] - xyzs[3*j+2];
        d = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
        if (d <= mxedg) continue;
        if (d > dist) {
          max  = j;
          dist = d;
        }
      }
    }
    if (map->nv == MAXUVS) break;
    if (max != -1) {
      for (j = map->nv; j >= max; j--) map->vs[j] = map->vs[j-1];
      map->vs[max] = 0.5*(map->vs[max+1] + map->vs[max-1]);
      map->nv++;
    }
  } while (max != -1);

}


static int parse_args(const int argc, char *argv[], char *ifile)
{

    int i;

    printf("********** egads2cgt version %s **********\n", version);

    for (i = 1; i <= argc-1; i++) {

      if (strcmp("-i", argv[i]) == 0) {
	 i++;
         strcpy(ifile, argv[i]);
      } else if (strcmp("-aflr4", argv[i]) == 0) {
        aflr4 = 1;
      } else if (strcmp("-q", argv[i]) == 0) {
 	 wrtqud = 1;
      } else if (strcmp("-uv", argv[i]) == 0) {
 	 wrtuv = 1;
      } else if (strcmp("-maxa", argv[i]) == 0) {
	 i++;
         sscanf(argv[i], "%lf", &mxang);
      } else if (strcmp("-maxe", argv[i]) == 0) {
	 i++;
         sscanf(argv[i], "%lf", &mxedg);
      } else if (strcmp("-maxc", argv[i]) == 0) {
	 i++;
         sscanf(argv[i], "%lf", &mxchd);
      } else if (strcmp("-ggf", argv[i]) == 0) {
        i++;
        sscanf(argv[i], "%lf", &ggf);
      } else {
	printf("Usage: egads2cgt [argument list] (defaults are in parenthesis)\n");
        printf("   -i input geometry filename <*.egads, *.stp, *.igs> ()\n");
        printf("   -aflr4 use AFLR4 for surface triangulation\n");
        printf("   -q write structured patches from quadding scheme to plot3d surface grid file\n");
        printf("   -uv write structured patches from uv evaluation to plot3d surface grid file\n");
        printf("   -ggf  factr <geometric growth factor with -uv for isocline smoothing> (1.2)\n");
        printf("   -maxa mxang <Max allow dihedral angle (deg)> (15.0)\n");
        printf("   -maxe mxedg <Max allow edge length> (0.025 * size)\n");
        printf("   -maxc mxchd <Max allow chord-height tolerance> (0.001 * size)\n");
        return 1;
      }
    }

    return 0;
}


int main(int argc, char *argv[])

{
  int          i, j, k, m, mm, n, status, oclass, mtype, nbody, nvert, ntriang;
  int          nnode, nedge, nloop, nface, dim, ptype, pindex;
  int          qnverts, npatch, npatchtot, nn, ip, ii, jj, iroot, iper, count;
  const int    *qptype, *qpindex, *pvindex, *pvbounds;
  int          *triang, *comp, *nipatch, *njpatch, *senses;
  char         ifile[120], filename[120], rootname[120], trifilename[120];
  char         p3dfilename[120], uvfilename[120];
  const char   *OCCrev;
  double       params[3], box[6], size, *verts, range[4], uv[2], result[18];
  double       qparams[3], alen, tol, *ndist;
  const double *qxyz, *quv;
  FILE         *fp;
  ego          context, model, geom, *bodies, tess, body;
  ego          *dum, *nodes, *edges, *loops, *faces, *tesses;
  UVmap        *map;

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d %s\n\n", i, j, OCCrev);
  
  /* get arguments */
  status = parse_args(argc, argv, ifile);
  if (status != 0) exit(status);
  
  printf("mxang:    %lf\n", mxang);
  printf("mxedg:    %lf\n", mxedg);
  printf("mxchd:    %lf\n", mxchd);
  if (wrtuv == 1) printf("ggf:      %lf\n", ggf);

  /* Build filenames from rootname */
  iroot = strroot(ifile);
  printf("iroot:    %d\n", iroot);
  strncpy(rootname, ifile, iroot);
  rootname[iroot] = 0;
  printf("rootname: %s\n", rootname);

  /* Initialize */
  status = EG_open(&context);
  if (status != EGADS_SUCCESS) {
    printf(" EG_open = %d!\n\n", status);
    return 1;
  }
  status = EG_loadModel(context, 0, ifile, &model);
  if (status != EGADS_SUCCESS) {
    printf(" EG_loadModel = %d\n\n", status);
    EG_close(context);
    return 1;
  }
  status = EG_getBoundingBox(model, box);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getBoundingBox = %d\n\n", status);
    EG_deleteObject(model);
    EG_close(context);
    return 1;
  }

  size = sqrt((box[0]-box[3])*(box[0]-box[3]) + (box[1]-box[4])*(box[1]-box[4]) +
              (box[2]-box[5])*(box[2]-box[5]));

  /* Get all Bodies */
  status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &triang);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n\n", status);
    EG_deleteObject(model);
    EG_close(context);
    return 1;
  }
  printf(" Number of Bodies = %d\n\n", nbody);

  /* Set tesselation parameters */
  if        (mxedg == 0.0) {
     params[0] = 0.025*size;
  } else if (mxedg >  0.0) {
     params[0] = mxedg;
  } else if (mxedg <  0.0) {
     params[0] = fabs(mxedg*size);
  }

  if        (mxchd == 0.0) {
     params[1] =  0.001*size;
  } else if (mxchd >  0.0) {
     params[1] = mxchd;
  } else if (mxchd <  0.0) {
     params[1] = fabs(mxchd*size);
  }
  if        (mxang == 0.0) {
     params[2] = 15.0;
  } else {
     params[2] = fabs(mxang);
  }
  printf(" Tess params: %lf %lf %lf\n\n", params[0], params[1], params[2]);
  tesses = (ego *) EG_alloc(nbody*sizeof(ego));
  if (tesses == NULL) {
      printf(" Alloc error on %d EGOs\n\n", nbody);
      EG_deleteObject(model);
      EG_close(context);
      return 1;
  }
  for (i = 0; i < nbody; i++) tesses[i] = NULL;
  
  /* use AFLR4? */
  if (aflr4 == 1) {
    status = aflr4egads(model, tesses);
    if (status != EGADS_SUCCESS) {
      printf(" aflr4egads = %d\n\n", status);
      for (i = 0; i < nbody; i++)
        if (tesses[i] != NULL) EG_deleteObject(tesses[i]);
      EG_free(tesses);
      EG_deleteObject(model);
      EG_close(context);
      return 1;
    }
    printf(" Tessellations completed!\n\n");
  }

  /* ---------------- Tesselate each body ---------------- */

  for (i = 0; i < nbody; i++) {
    snprintf(filename, 120, "%s.%3.3d.p3d", rootname, i+1);

    /* Add attribute */
    qparams[0] = 0.15;
    qparams[1] = 10.0;
    qparams[2] = 0.0;
    status = EG_attributeAdd(bodies[i], ".qParams", ATTRREAL, 3, NULL, qparams,
                             NULL);
    if (status != EGADS_SUCCESS) {
      printf(" Body %d: EG_attributeAdd = %d\n", i+1, status);
    }

    mtype  = 0;
    status = EG_getTopology(bodies[i], &geom, &oclass, &mtype, NULL, &j,
                            &dum, &triang);
    if (status != EGADS_SUCCESS) {
      printf(" Body %d: EG_getTopology = %d\n", i+1, status);
      continue;
    }

    status = EG_tolerance(bodies[i], &tol);
    if (status != EGADS_SUCCESS) {
      printf(" Body %d: EG_tolerance = %d\n", i+1, status);
      continue;
    }
    if        (mtype == WIREBODY) {
      printf(" Body %2d: Type = WireBody  -- tolerance = %le\n", i+1, tol);
    } else if (mtype == FACEBODY) {
      printf(" Body %2d: Type = FaceBody  -- tolerance = %le\n", i+1, tol);
    } else if (mtype == SHEETBODY) {
      printf(" Body %2d: Type = SheetBody -- tolerance = %le\n", i+1, tol);
    } else {
      printf(" Body %2d: Type = SolidBody -- tolerance = %le\n", i+1, tol);
    }

    if (aflr4 == 0) {
      status = EG_makeTessBody(bodies[i], params, &tess);
      if (status != EGADS_SUCCESS) {
        printf(" EG_makeTessBody %d = %d\n", i, status);
        continue;
      }
      printf("          Tessellation completed!\n");
      tesses[i] = tess;
    } else {
      tess = tesses[i];
    }
    
    status = EG_getBodyTopos(bodies[i], NULL, NODE, &nnode, &nodes);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos N %d = %d\n", i, status);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, EDGE, &nedge, &edges);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos E %d = %d\n", i, status);
      EG_free(nodes);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, LOOP, &nloop, &loops);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos L %d = %d\n", i, status);
      EG_free(edges);
      EG_free(nodes);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos F %d = %d\n", i, status);
      EG_free(loops);
      EG_free(edges);
      EG_free(nodes);
      continue;
    }
    
    /* ------------------- Write out tess owner file ----------------------- */

    status = EG_statusTessBody(tess, &body, &n, &nvert);
    if (status != EGADS_SUCCESS) {
      printf(" EG_statusTessBody = %d!\n", status);
      EG_free(faces);
      EG_free(loops);
      EG_free(edges);
      EG_free(nodes);
      continue;
    }
    comp = (int *) EG_alloc(nnode*sizeof(int));
    if (comp == NULL) {
      printf(" Cannot Allocate %d ints!\n", nnode);
      EG_free(faces);
      EG_free(loops);
      EG_free(edges);
      EG_free(nodes);
      continue;
    }
    for (n = 1; n <= nvert; n++) {
      status = EG_getGlobal(tess, n, &ptype, &pindex, NULL);
      if (status != EGADS_SUCCESS) {
        printf(" Vert %d: EG_getGlobal = %d\n", n, status);
        continue;
      }
      if (ptype != 0) continue;
      comp[pindex-1] = n;
    }

    /* write it out */
    snprintf(trifilename, 120, "%s.%3.3d.tess", rootname, i+1);
    fp = fopen(trifilename, "w");
    if (fp == NULL) {
      printf(" Cannot Open file %s -- NO FILE WRITTEN\n", trifilename);
      EG_free(comp);
      EG_free(faces);
      EG_free(loops);
      EG_free(edges);
      EG_free(nodes);
      continue;
    }
    printf(" Writing EGADS tess file: %s\n", trifilename);
    
    /* header */
    fprintf(fp, " %6d %6d %6d %6d %6d\n", mtype, nnode, nedge, nloop, nface);
    writeAttr(fp, bodies[i], NULL);

    /* Nodes */
    for (n = 1; n <= nnode; n++) {
      fprintf(fp, " %6d %6d\n", n, comp[n-1]);
      writeAttr(fp, nodes[n-1], NULL);
    }
    EG_free(comp);
    EG_free(nodes);
    
    /* Edges */
    for (n = 1; n <= nedge; n++) {
      status = EG_getTopology(edges[n-1], &geom, &oclass, &mtype, range, &j,
                              &dum, &triang);
      if (status != EGADS_SUCCESS) {
        printf(" %d: EG_getTopology = %d!\n", n, status);
        fprintf(fp, " %6d %6d\n", n, 0);
        fprintf(fp, " %6d\n", 0);
        continue;
      }
      if (mtype == DEGENERATE) {
        fprintf(fp, " %6d %6d\n", n, 0);
        writeAttr(fp, edges[n-1], NULL);
        continue;
      }

      status = EG_getTessEdge(tess, n, &m, &qxyz, &quv);
      if (status != EGADS_SUCCESS) {
        printf(" %d: EG_getTessEdge = %d!\n", n, status);
        fprintf(fp, " %6d %6d\n", n, 0);
        writeAttr(fp, edges[n-1], NULL);
        continue;
      }

      status = EG_arcLength(edges[n-1], range[0], range[1], &alen);
      if (status != EGADS_SUCCESS) {
        printf(" %d: EG_arcLength = %d!\n", n, status);
      } else {
	if (alen < 0.01)
	  printf(" Edge %2d: arc length = %le, number of pts = %d\n",
                 n, alen, m);
      }

      fprintf(fp, " %6d %6d\n", n, m);
      for (j = 1; j <= m; j++) {
        dim = 0;
        if (j == 1) {
          status = EG_indexBodyTopo(bodies[i], dum[0]);
          if (status > EGADS_SUCCESS) dim = status;
        }
        if (j == m) {
          if (mtype == TWONODE) {
            status = EG_indexBodyTopo(bodies[i], dum[1]);
          } else {
            status = EG_indexBodyTopo(bodies[i], dum[0]);
          }
          if (status > EGADS_SUCCESS) dim = status;
        }
        status = EG_localToGlobal(tess, -n, j, &k);
        if (status != EGADS_SUCCESS)
          printf(" %d/%d: EG_localToGlobal Edge = %d!\n", n, j, status);
        fprintf(fp, " %6d %20.13le %6d  ", k, quv[j-1], dim);
        if (j%2 == 0) fprintf(fp,"\n");
      }
      if ((j-1)%2 != 0) fprintf(fp,"\n");
      writeAttr(fp, edges[n-1], NULL);
    }
    
    /* Loops */
    for (j = 1; j <= nloop; j++) {
      status = EG_getTopology(loops[j-1], &geom, &oclass, &mtype, NULL, &m,
                              &dum, &senses);
      if (status != EGADS_SUCCESS) {
        printf(" Body %d: EG_getTopology L %d = %d\n", i+1, j, status);
        continue;
      }
      fprintf(fp, " %6d %6d\n", j, m);
      for (k = 0; k < m; k++)
        fprintf(fp, " %6d %6d\n",
                EG_indexBodyTopo(bodies[i], dum[k]), senses[k]);
      writeAttr(fp, loops[j-1], NULL);
    }
    EG_free(loops);
    
    /* Faces */
    for (j = 1; j <= nface; j++) {
      status = EG_getTopology(faces[j-1], &geom, &oclass, &mtype, NULL, &mm,
                              &dum, &senses);
      if (status != EGADS_SUCCESS) {
        printf(" Body %d: EG_getTopology F %d = %d\n", i+1, j, status);
        continue;
      }
      printf(" Face %2d: surface type = %d\n", j, geom->mtype);
      status = EG_getTessFace(tess, j, &m, &qxyz, &quv, &qptype, &qpindex,
                              &n, &pvindex, &pvbounds);
      if (status != EGADS_SUCCESS) {
        printf(" %d: EG_getTessFace = %d!\n", j, status);
        continue;
      }
      fprintf(fp, " %6d %6d %6d %6d\n", j, mm, mtype, m);
      for (k = 0; k < mm; k++)
        fprintf(fp, " %6d", EG_indexBodyTopo(bodies[i], dum[k]));
      fprintf(fp, "\n");
      for (k = 1; k <= m; k++) {
        status = EG_localToGlobal(tess, j, k, &n);
        if (status != EGADS_SUCCESS)
          printf(" %d/%d: EG_localToGlobal Face = %d!\n", j, k, status);
        fprintf(fp, " %6d %20.13le %20.13le %6d %6d\n",
                n, quv[2*k-2], quv[2*k-1], qptype[k-1], qpindex[k-1]);
      }
      writeAttr(fp, faces[j-1], NULL);
    }
    fclose(fp);
    
    /* ---------- Write quads file (plot3d format) for each body ----------- */

    if (wrtqud == 1) {
      
      qparams[0] = qparams[1] = qparams[2] = 0.0;
      
      /* First count total number of quad patches to write and write to file */
      npatchtot = 0;
      EG_setOutLevel(context, 0);
      for (n = 0; n < nface; n++) {
        status = EG_makeQuads(tess, qparams, n+1);
        if (EGADS_SUCCESS == status) {
          status    = EG_getQuads(tess, n+1, &qnverts, &qxyz, &quv, &qptype,
                                  &qpindex, &npatch);
          npatchtot = npatchtot + npatch;
        }
      }
      EG_setOutLevel(context, 1);
      printf(" Total number of quad patches %d \n", npatchtot);
      
      if (npatchtot > 0) {
        snprintf(p3dfilename, 120, "%s.%3.3d.p3d", rootname, i+1);
        printf(" Writing PLOT3D Quadded file: %s\n", p3dfilename);
        
        nipatch = njpatch = NULL;
        fp      = fopen(p3dfilename, "w");
        if (fp == NULL) {
          printf(" Error Opening %s!\n", p3dfilename);
          goto bailQuad;
        }
        fprintf(fp, "%d\n", npatchtot);
        
        /* Now get dimensions of all patches and write to file */
        nipatch = (int *) EG_alloc(npatchtot*sizeof(int));
        njpatch = (int *) EG_alloc(npatchtot*sizeof(int));
        if ((nipatch == NULL) || (njpatch == NULL)) {
          printf(" Malloc ERROR on %d Patches!\n", npatchtot);
          goto bailQuad;
        }
        
        EG_setOutLevel(context, 0);
        for (nn = n = 0; n < nface; n++) {
          status = EG_makeQuads(tess, qparams, n+1) ;
          if (EGADS_SUCCESS == status) {
            status = EG_getQuads(tess, n+1, &qnverts, &qxyz, &quv, &qptype,
                                 &qpindex, &npatch);
            for (ip = 0; ip < npatch; ip++, nn++) {
              status = EG_getPatch(tess, n+1, ip+1, &nipatch[nn], &njpatch[nn],
                                   &pvindex, &pvbounds);
            }
          }
        }
        
        for (ip = 0; ip < npatchtot; ip++)
          fprintf(fp, "%d %d %d\n", nipatch[ip], njpatch[ip],1);
        
        /* Write x,y,z of each quad patch */
        for (nn = n = 0; n < nface; n++) {
          
          status = EG_makeQuads(tess, qparams, n+1);
          if (status == EGADS_SUCCESS) {
            status = EG_getQuads(tess, n+1, &qnverts, &qxyz, &quv, &qptype,
                                 &qpindex, &npatch);
            for (ip = 0; ip < npatch; ip++, nn++) {
              status = EG_getPatch(tess, n+1, ip+1, &nipatch[nn], &njpatch[nn],
                                   &pvindex, &pvbounds);
              for (count = jj = 0; jj < njpatch[nn]; jj++) {
                for (ii = 0; ii < nipatch[nn]; ii++) {
                  m = pvindex[jj*nipatch[nn] + ii] - 1;
                  fprintf(fp, "%20.13le ", qxyz[3*m]);
                  count++;
                  if (count%5 == 0) fprintf(fp, "\n");
                }
              }
              fprintf(fp, "\n");
              
              for (count = jj = 0; jj < njpatch[nn]; jj++) {
                for (ii = 0; ii < nipatch[nn]; ii++) {
                  m = pvindex[jj*nipatch[nn] + ii] - 1;
                  fprintf(fp, "%20.13le ", qxyz[3*m+1]);
                  count++;
                  if (count%5 == 0) fprintf(fp, "\n");
                }
              }
              fprintf(fp, "\n");
              
              for (count = jj = 0; jj < njpatch[nn]; jj++) {
                for (ii = 0; ii < nipatch[nn]; ii++) {
                  m = pvindex[jj*nipatch[nn] + ii] - 1;
                  fprintf(fp, "%20.13le ", qxyz[3*m+2]);
                  count++;
                  if (count%5 == 0) fprintf(fp, "\n");
                }
              }
              fprintf(fp, "\n");
              
            }
          }
        }
        EG_setOutLevel(context, 1);
      bailQuad:
        if (nipatch != NULL) EG_free(nipatch);
        if (njpatch != NULL) EG_free(njpatch);
        if (fp      != NULL) fclose(fp);
      }
    }

   /* ----------- Write uv file (plot3d format) for each Body ------------ */

    if (wrtuv == 1) {

      ndist = (double *) EG_alloc(nnode*sizeof(double));
      map   = (UVmap *)  EG_alloc(nface*sizeof(UVmap));
      if ((map == NULL) || (ndist == NULL)) {
        printf(" ERROR Allocating %d Maps!\n", nface);
        if (map   != NULL) EG_free(map);
        if (ndist != NULL) EG_free(ndist);
      } else {
        for (n = 0; n < nnode; n++) ndist[n] = 0.0;
        getNodeSpacing(tess, nedge, edges, ndist);
        for (n = 0; n < nnode; n++)
          printf(" Node %2d: spacing = %lf\n", n+1, ndist[n]);
        snprintf(uvfilename, 120, "%s.%3.3d.uv", rootname, i+1);
        printf(" Writing PLOT3D Untrimmed file: %s\n", uvfilename);
        
        fp = fopen(uvfilename, "w");
        if (fp == NULL) {
          printf(" Error Opening %s!\n", uvfilename);
        } else {
          
          /* write out the uv file a Face at a time */
          for (n = 0; n < nface; n++) {
            status = EG_getRange(faces[n], range, &iper);
            if (status != EGADS_SUCCESS)
              printf(" %d: EG_getRange = %d\n", n+1, status);
            /* start out 9x9 regardless & squeeze in a little */
            map[n].nu = 9;
            map[n].nv = 9;
            for (jj = 0; jj < map[n].nv; jj++)
              map[n].vs[jj] =  range[2] +
                              (range[3] - range[2])*(double)jj/(map[n].nv-1);
            map[n].vs[0] += 1.e-5*(range[3] - range[2]);
            map[n].vs[8] -= 1.e-5*(range[3] - range[2]);
            for (ii = 0; ii < map[n].nu; ii++)
              map[n].us[ii] =  range[0] +
                              (range[1] - range[0])*(double)ii/(map[n].nu-1);
            map[n].us[0] += 1.e-5*(range[1] - range[0]);
            map[n].us[8] -= 1.e-5*(range[1] - range[0]);
            /* enhance based on usual tessellation parameters */
            updateMap(faces[n], params[0], params[1], params[2], &map[n]);
            /* insert Node spacings into both U and V */
            insertNodeSpacing(tess, faces[n], n+1, ndist, &map[n]);
            /* smooth the Us and Vs based on the geometric growth factor */
            smoothMap(faces[n], ggf, &map[n]);
          }
          fprintf(fp, "%d\n", nface);
          for (n = 0; n < nface; n++)
            fprintf(fp,"%d %d %d \n", map[n].nu, map[n].nv, 1);
          
          for (n = 0; n < nface; n++) {

            status = EG_getTopology(faces[n], &geom, &oclass, &mtype, NULL, &mm,
                                    &dum, &senses);
            if (status != EGADS_SUCCESS) {
               printf(" Body %d: EG_getTopology F %d = %d\n", i+1, n+1, status);
               continue;
            }
          
            for (dim = 0; dim < 3; dim++) {
              for (count = jj = 0; jj < map[n].nv; jj++) {
                uv[1] = map[n].vs[jj];
                for (ii = 0; ii < map[n].nu; ii++) {
                  if ( mtype == SFORWARD ) {
                     uv[0] = map[n].us[ii];
                  } else {
                     uv[0] = map[n].us[map[n].nu-ii-1];
                  }
                  status = EG_evaluate(faces[n], uv, result);
                  if (status != EGADS_SUCCESS)
                    printf(" %d: EG_evaluate = %d\n", n+1, status);
                  fprintf(fp, "%20.13le ", result[dim]);
                  count++;
                  if (count%5 == 0) fprintf(fp, "\n");
                }

                if (count%5 != 0) fprintf(fp, "\n");
              }
            }

              for (count = jj = 0; jj < map[n].nv; jj++) {
                uv[1] = map[n].vs[jj];
                for (ii = 0; ii < map[n].nu; ii++) {
                  if (mtype == SFORWARD) {
                     uv[0] = map[n].us[ii];
                  } else {
                     uv[0] = map[n].us[map[n].nu-ii-1];
                  }
#ifdef INFACEOCC
                  status = EG_inFaceOCC(faces[n], 1.5*tol, uv);
#else
                  status = EG_inFace(faces[n], uv);
#endif
                  if (status < 0)
		    printf(" face,ii,jj %d %d %d: inFace = %d\n",
                           n+1, ii, jj, status);
                  fprintf(fp, "%d ", 1-status);
                  count++;
                  if (count%15 == 0) fprintf(fp, "\n");
                }
                if (count%15 != 0) fprintf(fp, "\n");
              }
            
          }
          fclose(fp);
        }
        EG_free(map);
        EG_free(ndist);
      }
    }
    EG_free(edges);
    EG_free(faces);

   /* --------- Zip up the tessellation and write in CART3D format --------- */

    status = bodyTessellation(tess, nface, &nvert, &verts,
                              &ntriang, &triang, &comp);
    if (status != EGADS_SUCCESS) continue;

    /* write it out */
    snprintf(trifilename, 120, "%s.%3.3d.tri", rootname, i+1);

    fp = fopen(trifilename, "w");
    if (fp == NULL) {
      printf(" Can not Open file %s -- NO FILE WRITTEN\n", trifilename);
      EG_free(verts);
      EG_free(triang);
      continue;
    }
    printf(" Writing CART3D tri file: %s\n", trifilename);

    /* header */
    fprintf(fp, "%d  %d\n", nvert, ntriang);

    /* ...vertList     */
    for (j = 0; j < nvert; j++)
      fprintf(fp, " %20.13le %20.13le %20.13le\n",
              verts[3*j  ], verts[3*j+1], verts[3*j+2]);

    /* ...Connectivity */
    for (j = 0; j < ntriang; j++)
      fprintf(fp, "%6d %6d %6d\n", triang[3*j], triang[3*j+1], triang[3*j+2]);

    /* ...Component list*/
    for (j = 0; j < ntriang; j++)
      fprintf(fp, "%6d\n", comp[j]);

    fclose(fp);

    printf("      # verts = %d,  # tris = %d\n\n", nvert, ntriang);

    EG_free(verts);
    EG_free(triang);
    EG_free(comp);
  }

  for (i = 0; i < nbody; i++)
    if (tesses[i] != NULL) EG_deleteObject(tesses[i]);
  EG_free(tesses);
  status = EG_deleteObject(model);
  if (status != EGADS_SUCCESS) printf(" EG_deleteObject = %d\n", status);

  EG_close(context);

  return 0;
}
