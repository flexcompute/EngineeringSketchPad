/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Quad Patch functions
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egadsPatch.h"

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


int
EG_quadEdges(ego tess, quadPatch *patch)
{
  int    i, j, index, stat, nnode, nedge, oclass, mtype, eIndex;
  int    side, off, len, start;
  int    nin[4], ein[2], *senses;
  double d, dist, tol, xyz[3], tlims[2], *t, *ll, *lr, *ur, *ul, *xyzs = NULL;
  ego    body, ref, *objs, *dum, *nodes = NULL, *edges = NULL;
  const double *txyzs, *tts;
  
  for (i = 0; i < 4; i++)
    patch->en[i][0] = patch->en[i][1] = patch->en[i][2] = -1;
  
  /* get corners (Nodes) */
  ll    = &patch->xyzs[0];
  index = (patch->ni-1)*patch->nj;
  lr    = &patch->xyzs[3*index];
  index = (patch->ni-1)*patch->nj + patch->nj-1;
  ur    = &patch->xyzs[3*index];
  index = patch->nj-1;
  ul    = &patch->xyzs[3*index];
  
  /* get EGOs */
  stat = EG_statusTessBody(tess, &body, &i, &index);
  if (stat < EGADS_SUCCESS) {
    printf(" EGADS Error: EG_statusTessBody = %d (EG_quadEdges)!\n", stat);
    return stat;
  }
  stat = EG_getBodyTopos(body, patch->face, NODE, &nnode, &nodes);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Node = %d (EG_quadEdges)!\n", stat);
    return stat;
  }
  if ((nnode > 4) || (nodes == NULL)) {
    printf(" EGADS Error: #Node = %d (EG_quadEdges)!\n", nnode);
    stat = EGADS_TOPOERR;
    goto eCleanup;
  }
  stat = EG_getBodyTopos(body, patch->face, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Node = %d (EG_quadEdges)!\n", stat);
    goto eCleanup;
  }
  if (edges == NULL) {
    printf(" EGADS Error: Edge Objects is NULL (EG_quadEdges)!\n");
    stat = EGADS_TOPOERR;
    goto eCleanup;
  }
  
  /* get the storage */
  len  = patch->ni;
  if (patch->nj > len) len = patch->nj;
  xyzs = (double *) EG_alloc(4*len*sizeof(double));
  if (xyzs == NULL) {
    printf(" EGADS Error: Malloc on %d Edge verts (EG_quadEdges)!\n", len);
    stat = EGADS_MALLOC;
    goto eCleanup;
  }
  t = &xyzs[3*len];
  
  /* match the Nodes */
  for (i = 0; i < nnode; i++) {
    stat = EG_getTopology(nodes[i], &ref, &oclass, &mtype, xyz, &index, &objs,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_getTopology Node %d = %d (EG_quadEdges)!\n",
             i+1, stat);
      goto eCleanup;
    }
    stat = EG_tolerance(nodes[i], &tol);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_tolerance Node %d = %d (EG_quadEdges)!\n",
             i+1, stat);
      goto eCleanup;
    }
    dist   =  1.e300;
    nin[i] = -1;
    d = sqrt((ll[0]-xyz[0])*(ll[0]-xyz[0]) + (ll[1]-xyz[1])*(ll[1]-xyz[1]) +
             (ll[2]-xyz[2])*(ll[2]-xyz[2]));
    if (d < dist) {
      dist   = d;
      nin[i] = 0;
    }
    d = sqrt((lr[0]-xyz[0])*(lr[0]-xyz[0]) + (lr[1]-xyz[1])*(lr[1]-xyz[1]) +
             (lr[2]-xyz[2])*(lr[2]-xyz[2]));
    if (d < dist) {
      dist   = d;
      nin[i] = 1;
    }
    d = sqrt((ur[0]-xyz[0])*(ur[0]-xyz[0]) + (ur[1]-xyz[1])*(ur[1]-xyz[1]) +
             (ur[2]-xyz[2])*(ur[2]-xyz[2]));
    if (d < dist) {
      dist   = d;
      nin[i] = 2;
    }
    d = sqrt((ul[0]-xyz[0])*(ul[0]-xyz[0]) + (ul[1]-xyz[1])*(ul[1]-xyz[1]) +
             (ul[2]-xyz[2])*(ul[2]-xyz[2]));
    if (d < dist) {
      dist   = d;
      nin[i] = 3;
    }
    if (dist > tol) {
      printf(" EGADS Error: Tolerance on Node %d: %le (%le) (EG_quadEdges)!\n",
             i+1, dist, tol);
      stat = EGADS_TOPOERR;
      goto eCleanup;
    }
  }
#ifdef DEBUG
  printf("\n Face %d\n", EG_indexBodyTopo(body, patch->face));
#endif
  
  /* find the Edges */
  for (i = 0; i < nedge; i++) {
    eIndex = EG_indexBodyTopo(body, edges[i]);
    if (eIndex <= EGADS_SUCCESS) {
      printf(" EGADS Error: Edge %d EG_indexBodyTopo = %d (EG_quadEdges)!\n",
             i+1, eIndex);
      stat = EGADS_TOPOERR;
      goto eCleanup;
    }
    stat = EG_tolerance(edges[i], &tol);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_tolerance Edge %d = %d (EG_quadEdges)!\n",
             i+1, stat);
      goto eCleanup;
    }
    stat = EG_getTopology(edges[i], &ref, &oclass, &mtype, tlims, &index, &objs,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_getTopology Edge %3d = %d (EG_quadEdges)!\n",
             eIndex, stat);
      goto eCleanup;
    }
    if (mtype == DEGENERATE) continue;
    if (mtype == ONENODE) {
      printf(" EGADS Error: Edge %3d: ONENODE (EG_quadEdges)!\n", eIndex);
      stat = EGADS_TOPOERR;
      goto eCleanup;
    }
    ein[0] = ein[1] = -1;
    for (j = 0; j < nnode; j++) {
      if (objs[0] == nodes[j]) ein[0] = nin[j];
      if (objs[1] == nodes[j]) ein[1] = nin[j];
    }
    if ((ein[0] == -1) || (ein[1] == -1)) {
      printf(" EGADS Error: Edge %3d Cannot find Node %d %d (EG_quadEdges)!\n",
             eIndex, ein[0], ein[1]);
      stat = EGADS_TOPOERR;
      goto eCleanup;
    }
    side  = -1;
    start = off = 0;
    if ((ein[0] == 0) && (ein[1] == 1)) {
      side  = 0;
      start =  0;
      len   =  patch->ni;
      off   =  patch->nj;
    } else if ((ein[0] == 1) && (ein[1] == 0)) {
      side  = 0;
      start = (patch->ni-1)*patch->nj;
      len   =  patch->ni;
      off   = -patch->nj;
    } else if ((ein[0] == 1) && (ein[1] == 2)) {
      side  = 1;
      start = (patch->ni-1)*patch->nj;
      len   =  patch->nj;
      off   =  1;
    } else if ((ein[0] == 2) && (ein[1] == 1)) {
      side  = 1;
      start = (patch->ni-1)*patch->nj + patch->nj-1;
      len   =  patch->nj;
      off   = -1;
    } else if ((ein[0] == 2) && (ein[1] == 3)) {
      side  = 2;
      start = (patch->ni-1)*patch->nj + patch->nj-1;
      len   =  patch->ni;
      off   = -patch->nj;
    } else if ((ein[0] == 3) && (ein[1] == 2)) {
      side  = 2;
      start =  patch->nj-1;
      len   =  patch->ni;
      off   =  patch->nj;
    } else if ((ein[0] == 3) && (ein[1] == 0)) {
      side  = 3;
      start =  patch->nj-1;
      len   =  patch->nj;
      off   =  -1;
    } else if ((ein[0] == 0) && (ein[1] == 3)) {
      side  = 3;
      start =  0;
      len   =  patch->nj;
      off   =  1;
    }
    if (side != -1) {
      patch->en[side][0] = eIndex;
      patch->en[side][1] = start;
      patch->en[side][2] = off;
    }
#ifdef DEBUG
    printf(" Edge %3d: Side = %d %3d %4d %4d, node indices = %d %d\n",
           eIndex, side, len, start, off, ein[0], ein[1]);
#endif
    
    stat = EG_getTopology(objs[0], &ref, &oclass, &mtype, xyz, &index, &dum,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Edge %d EG_getTopology Node0 = %d (EG_quadEdges)!\n",
             eIndex, stat);
      goto eCleanup;
    }
    xyzs[0] = xyz[0];
    xyzs[1] = xyz[1];
    xyzs[2] = xyz[2];
    t[0]    = tlims[0];
    stat = EG_getTopology(objs[1], &ref, &oclass, &mtype, xyz, &index, &dum,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Edge %d EG_getTopology Node0 = %d (EG_quadEdges)!\n",
             eIndex, stat);
      goto eCleanup;
    }
    xyzs[3*len-3] = xyz[0];
    xyzs[3*len-2] = xyz[1];
    xyzs[3*len-1] = xyz[2];
    t[len-1]      = tlims[1];
    for (j = 1; j < len-1; j++) {
      start += off;
      stat = EG_invEvaluate(edges[i], &patch->xyzs[3*start], &t[j], &xyzs[3*j]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Edge %d invEvaluate %d = %d (EG_quadEdges)!\n",
               eIndex, j, stat);
        goto eCleanup;
      }
      d = sqrt((patch->xyzs[3*start  ]-xyzs[3*j  ])*
               (patch->xyzs[3*start  ]-xyzs[3*j  ]) +
               (patch->xyzs[3*start+1]-xyzs[3*j+1])*
               (patch->xyzs[3*start+1]-xyzs[3*j+1]) +
               (patch->xyzs[3*start+2]-xyzs[3*j+2])*
               (patch->xyzs[3*start+2]-xyzs[3*j+2]));
      if (d > tol) {
        printf(" EGADS Error: Edge %d tolerance %d %le (%le) (EG_quadEdges)!\n",
               eIndex, j, d, tol);
        stat = EGADS_TOPOERR;
        goto eCleanup;
      }
    }
    
    /* save it away */
    stat = EG_getTessEdge(tess, eIndex, &index, &txyzs, &tts);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Edge %d EG_getTessEdge = %d (EG_quadEdges)!\n",
             eIndex, stat);
      goto eCleanup;
    }
    if (index == 0) {
      stat = EG_setTessEdge(tess, eIndex, len, xyzs, t);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Edge %d EG_setTessEdge = %d (EG_quadEdges)!\n",
               eIndex, stat);
        goto eCleanup;
      }
    } else {
      if (index != len) {
        printf(" EGADS Error: Edge %d length = %d %d (EG_quadEdges)!\n",
               eIndex, index, len);
        stat = EGADS_TESSTATE;
        goto eCleanup;
      }
      for (j = 0; j < len; j++) {
        d = sqrt((txyzs[3*j  ]-xyzs[3*j  ])*(txyzs[3*j  ]-xyzs[3*j  ]) +
                 (txyzs[3*j+1]-xyzs[3*j+1])*(txyzs[3*j+1]-xyzs[3*j+1]) +
                 (txyzs[3*j+2]-xyzs[3*j+2])*(txyzs[3*j+2]-xyzs[3*j+2]));
        if (d > tol) {
          printf(" EGADS Error: Edge %d cmp toler %d %le (%le) (EG_quadEdges)!\n",
                 eIndex, j, d, tol);
          stat = EGADS_TOPOERR;
          goto eCleanup;
        }
      }
    }
  }
  
  stat = EGADS_SUCCESS;
  
eCleanup:
  if (xyzs  != NULL) EG_free(xyzs);
  if (edges != NULL) EG_free(edges);
  if (nodes != NULL) EG_free(nodes);
  return stat;
}


int
EG_quadFace(ego tess, quadPatch *patch)
{
  int    i, j, n, index, stat, fIndex, oclass, mtype, nloop, nedge, len, ntri;
  int    eIndex, start, off, rev, stype, aType, aLen, vrt[4];
  int    *senses, *tris = NULL;
  double area, uv[2], x1[3], x2[3], norm[3], fnorm[3], uvBox[4], result[18];
  double *xyzs, *uvs = NULL;
  ego    body, ref, surf, *loops, *edges;
  const int    *ints;
  const double *reals, *txyzs, *tts;
  const char   *str;
  
  stat = EG_statusTessBody(tess, &body, &i, &n);
  if (stat < EGADS_SUCCESS) return stat;
  
  fIndex = EG_indexBodyTopo(body, patch->face);
  if (fIndex <= EGADS_SUCCESS) {
    printf(" EGADS Error: Face EG_indexBodyTopo = %d (EG_quadFace)!\n",
           fIndex);
    stat = fIndex;
    goto fCleanup;
  }
  stat = EG_getTopology(patch->face, &surf, &oclass, &stype, uvBox, &nloop,
                        &loops, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Face %d EG_getTopology Loop = %d (EG_quadFace)!\n",
           fIndex, stat);
    goto fCleanup;
  }
  if (nloop != 1) {
    printf(" EGADS Error: Face %d  nLoops = %d (EG_quadFace)!\n",
           fIndex, nloop);
    stat = EGADS_TOPOERR;
    goto fCleanup;
  }
  stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, NULL, &nedge,
                        &edges, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Face %d EG_getTopology Edges = %d (EG_quadFace)!\n",
           fIndex, stat);
    goto fCleanup;
  }
  for (i = 0; i < nedge-1; i++)
    for (j = i+1; j < nedge; j++)
      if (edges[i] == edges[j]) {
        printf(" EGADS Error: Face %d Edge in Loop twice (EG_quadFace)!\n",
               fIndex);
        stat = EGADS_TOPOERR;
        goto fCleanup;
      }
  
  ntri = 2*(patch->ni-1)*(patch->nj-1);
  uvs  = (double *) EG_alloc(5*patch->ni*patch->nj*sizeof(double));
  if (uvs == NULL) {
    printf(" EGADS Error: Malloc on %dx%d Face UVs (EG_quadFace)!\n",
           patch->ni, patch->nj);
    stat = EGADS_MALLOC;
    goto fCleanup;
  }
  xyzs = &uvs[2*patch->ni*patch->nj];
  tris = (int *) EG_alloc(3*ntri*sizeof(int));
  if (tris == NULL) {
    printf(" EGADS Error: Malloc on %d tris (EG_quadFace)!\n", ntri);
    stat = EGADS_MALLOC;
    goto fCleanup;
  }
  
  /* get UVs at Edges */
  for (i = 0; i < 4; i++) {
    eIndex = patch->en[i][0];
    if (eIndex == -1) {
      printf(" EGADS Error: Face %d/Edge %d No internals (EG_quadFace)!\n",
             fIndex, eIndex);
      stat = EGADS_TOPOERR;
      goto fCleanup;
    }
    start = patch->en[i][1];
    off   = patch->en[i][2];
    stat  = EG_objectBodyTopo(body, EDGE, eIndex, &ref);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d/Edge %d EG_objectBodyTopo = %d (EG_quadFace)!\n",
             fIndex, eIndex, stat);
      goto fCleanup;
    }
    stat = EG_getTessEdge(tess, eIndex, &len, &txyzs, &tts);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d/Edge %d EG_getTessEdge = %d (EG_quadFace)!\n",
             fIndex, eIndex, stat);
      goto fCleanup;
    }
    
    for (j = 0; j < nedge; j++)
      if (edges[j] == ref) break;
    if (j == nedge) {
      printf(" EGADS Error: Face %d -- Edge %d Not Found (EG_quadFace)!\n",
             fIndex, eIndex);
      stat = EGADS_NOTFOUND;
      goto fCleanup;
    }
    for (n = 0; n < len; n++) {
      xyzs[3*start  ] = txyzs[3*n  ];
      xyzs[3*start+1] = txyzs[3*n+1];
      xyzs[3*start+2] = txyzs[3*n+2];
      stat = EG_getEdgeUV(patch->face, ref, senses[j], tts[n], &uvs[2*start]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d EG_getEdgeUV %d %d = %d (EG_quadFace)!\n",
               fIndex, n, j, stat);
        goto fCleanup;
      }
      start += off;
    }
  }
  
  /* get interior UVs */
  for (i = 1; i < patch->ni-1; i++)
    for (j = 1; j < patch->nj-1; j++) {
      index = i*patch->nj + j;
      stat  = EG_invEvaluate(surf, &patch->xyzs[3*index], &uvs[2*index],
                             &xyzs[3*index]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d invEvaluate %d %d = %d (EG_quadFace)!\n",
               fIndex, i, j, stat);
        goto fCleanup;
      }
    }
  
  /* get normal orientation */
  rev    = 1;
  vrt[0] = 0;
  vrt[1] = patch->nj;
  vrt[2] = patch->nj + 1;
  uv[0]  = (uvs[2*vrt[0]  ] + uvs[2*vrt[1]  ] + uvs[2*vrt[2]  ])/3.0;
  uv[1]  = (uvs[2*vrt[0]+1] + uvs[2*vrt[1]+1] + uvs[2*vrt[2]+1])/3.0;
  stat   = EG_evaluate(surf, uv, result);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Face %d evaluate = %d (EG_quadFace)!\n",
           fIndex, stat);
    goto fCleanup;
  }
  x1[0] = result[3];
  x1[1] = result[4];
  x1[2] = result[5];
  x2[0] = result[6];
  x2[1] = result[7];
  x2[2] = result[8];
  CROSS(norm, x1, x2);
  area = sqrt(DOT(norm, norm))*stype;
  if (area == 0.0) {
    printf(" EGADS Error: Face %d Zero Cross eval (EG_quadFace)!\n", fIndex);
    stat = EGADS_DEGEN;
    goto fCleanup;
  }
  norm[0] /= area;
  norm[1] /= area;
  norm[2] /= area;
  x1[0] = xyzs[3*vrt[1]  ] - xyzs[3*vrt[0]  ];
  x2[0] = xyzs[3*vrt[2]  ] - xyzs[3*vrt[0]  ];
  x1[1] = xyzs[3*vrt[1]+1] - xyzs[3*vrt[0]+1];
  x2[1] = xyzs[3*vrt[2]+1] - xyzs[3*vrt[0]+1];
  x1[2] = xyzs[3*vrt[1]+2] - xyzs[3*vrt[0]+2];
  x2[2] = xyzs[3*vrt[2]+2] - xyzs[3*vrt[0]+2];
  CROSS(fnorm, x1, x2);
  area = sqrt(DOT(fnorm, fnorm));
  if (area == 0.0) {
    printf(" EGADS Error: Face %d Zero Cross tri (EG_quadFace)!\n", fIndex);
    stat = EGADS_DEGEN;
    goto fCleanup;
  }
  fnorm[0] /= area;
  fnorm[1] /= area;
  fnorm[2] /= area;
  if (DOT(fnorm,norm) < 0.0) rev = -1;
  
  /* make the tris */
  for (n = i = 0; i < patch->ni-1; i++)
    for (j = 0; j < patch->nj-1; j++, n+=6) {
      if (rev == -1) {
        vrt[3]    =    i *patch->nj + j + 1;
        vrt[2]    = (i+1)*patch->nj + j + 1;
        vrt[1]    = (i+1)*patch->nj + j + 2;
        vrt[0]    =    i *patch->nj + j + 2;
      } else {
        vrt[0]    =    i *patch->nj + j + 1;
        vrt[1]    = (i+1)*patch->nj + j + 1;
        vrt[2]    = (i+1)*patch->nj + j + 2;
        vrt[3]    =    i *patch->nj + j + 2;
      }
      tris[n  ] = vrt[0];
      tris[n+1] = vrt[1];
      tris[n+2] = vrt[2];
      tris[n+3] = vrt[0];
      tris[n+4] = vrt[2];
      tris[n+5] = vrt[3];
    }
  
  /* set the Face tessellation */
  n    = patch->ni*patch->nj;
  stat = EG_setTessFace(tess, fIndex, n, xyzs, uvs, ntri, tris);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Face %d EG_setTessFace = %d (EG_quadFace)!\n",
           fIndex, stat);
    goto fCleanup;
  }
  
  /* make or adjust the .mixed attribute */
  stat = EG_attributeRet(tess, ".mixed", &aType, &aLen, &ints, &reals, &str);
  if (stat == EGADS_NOTFOUND) {
    i = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 6, NULL, NULL, "Mixed");
    if (i != EGADS_SUCCESS)
      printf(" EGADS Warning: Face %d attributeAdd Mixed = %d (EG_quadFace)!\n",
             fIndex, i);
    i = EG_getBodyTopos(body, NULL, FACE, &n, NULL);
    if (i != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d EG_getBodyTopos F = %d (EG_quadFace)!\n",
             fIndex, i);
      stat = i;
      goto fCleanup;
    }
    senses = (int *) EG_alloc(n*sizeof(int));
    if (senses == NULL) {
      printf(" EGADS Error: Face %d Attr Malloc %d (EG_quadFace)!\n",
             fIndex, n);
      stat = EGADS_MALLOC;
      goto fCleanup;
    }
    for (i = 0; i < n; i++) senses[i] = 0;
    senses[fIndex-1] = ntri/2;
    stat = EG_attributeAdd(tess, ".mixed", ATTRINT, n, senses, NULL, NULL);
    EG_free(senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d EG_attributeAdd = %d (EG_quadFace)!\n",
             fIndex, stat);
      goto fCleanup;
    }
  } else if (stat == EGADS_SUCCESS) {
    if (aType != ATTRINT) {
      printf(" EGADS Error: Face %d .mixed aType = %d (EG_quadFace)!\n",
             fIndex, aType);
      stat = EGADS_ATTRERR;
      goto fCleanup;
    }
    senses = (int *) EG_alloc(aLen*sizeof(int));
    if (senses == NULL) {
      printf(" EGADS Error: Face %d Attribute Malloc %d (EG_quadFace)!\n",
             fIndex, aLen);
      stat = EGADS_MALLOC;
      goto fCleanup;
    }
    for (i = 0; i < aLen; i++) senses[i] = ints[i];
    senses[fIndex-1] = ntri/2;
    for (j = i = 0; i < aLen; i++)
      if (senses[i] != 0) j++;
    if (j == aLen) {
      i = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 6, NULL, NULL, "Quad");
      if (i != EGADS_SUCCESS)
        printf(" EGADS Warning: Face %d attributeAdd Mix = %d (EG_quadFace)!\n",
               fIndex, i);
    }
    stat = EG_attributeAdd(tess, ".mixed", ATTRINT, aLen, senses, NULL, NULL);
    EG_free(senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d attributeAdd = %d (EG_quadFace)!\n",
             fIndex, stat);
      goto fCleanup;
    }
  } else {
    printf(" EGADS Error: Face %d EG_attributeRet = %d (EG_quadFace)!\n",
           fIndex, stat);
    goto fCleanup;
  }
  
  stat = EGADS_SUCCESS;
  
fCleanup:
  if (tris != NULL) EG_free(tris);
  if (uvs  != NULL) EG_free(uvs);
  return stat;
}


#ifdef STANDALONE
static void
fillPatch(quadPatch *patch, int i, int j, const double *xyz)
{
  int index;

  index = i*patch->nj + j;
  patch->xyzs[3*index  ] = xyz[0];
  patch->xyzs[3*index+1] = xyz[1];
  patch->xyzs[3*index+2] = xyz[2];
}


static void
fillInterior(quadPatch *patch)
{
  int    i, j, index;
  double xyz[3], xi, et, *xi0, *xin, *xj0, *xjn, *ll, *lr, *ul, *ur;
  
  /* get corners */
  ll    = &patch->xyzs[0];
  index = (patch->ni-1)*patch->nj;
  lr    = &patch->xyzs[3*index];
  index = (patch->ni-1)*patch->nj + patch->nj-1;
  ur    = &patch->xyzs[3*index];
  index = patch->nj-1;
  ul    = &patch->xyzs[3*index];
  
  /* TFI the interior */
  for (i = 1; i < patch->ni-1; i++) {
    xi  = i;
    xi /= patch->ni;
    for (j = 1; j < patch->nj-1; j++) {
      et     = j;
      et    /= patch->nj;
      index  = j;
      xi0    = &patch->xyzs[3*index];
      index  = (patch->ni-1)*patch->nj + j;
      xin    = &patch->xyzs[3*index];
      index  = i*patch->nj;
      xj0    = &patch->xyzs[3*index];
      index  = i*patch->nj + patch->nj-1;
      xjn    = &patch->xyzs[3*index];
      xyz[0] = (1.0-xi)            * xi0[0] +
               (    xi)            * xin[0] +
                          (1.0-et) * xj0[0] +
                          (    et) * xjn[0] -
               (1.0-xi) * (1.0-et) * ll[0]  -
               (1.0-xi) * (    et) * ul[0]  -
               (    xi) * (1.0-et) * lr[0]  -
               (    xi) * (    et) * ur[0];
      xyz[1] = (1.0-xi)            * xi0[1] +
               (    xi)            * xin[1] +
                          (1.0-et) * xj0[1] +
                          (    et) * xjn[1] -
               (1.0-xi) * (1.0-et) * ll[1]  -
               (1.0-xi) * (    et) * ul[1]  -
               (    xi) * (1.0-et) * lr[1]  -
               (    xi) * (    et) * ur[1];
      xyz[2] = (1.0-xi)            * xi0[2] +
               (    xi)            * xin[2] +
                          (1.0-et) * xj0[2] +
                          (    et) * xjn[2] -
               (1.0-xi) * (1.0-et) * ll[2]  -
               (1.0-xi) * (    et) * ul[2]  -
               (    xi) * (1.0-et) * lr[2]  -
               (    xi) * (    et) * ur[2];
      index  = i*patch->nj + j;
      patch->xyzs[3*index  ] = xyz[0];
      patch->xyzs[3*index+1] = xyz[1];
      patch->xyzs[3*index+2] = xyz[2];
    }
  }
}


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int          i, j, k, m, stat, aType, aLen, len, ntri;
  double       data[6], params[3];
  ego          context, bodies[2], tess, model, *objs;
  const int    *ints, *ptype, *pindex, *tris, *tric;
  const double *reals, *xyzs, *uvs;
  const char   *str;
  quadPatch    patches[6];
  
  /* create an EGADS context */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", stat);
    return 1;
  }
  
  /* make a box */
  data[0] = data[1] = data[2] = -1.0;
  data[3] = data[4] = data[5] =  2.0;
  stat = EG_makeSolidBody(context, BOX, data, &bodies[0]);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody box return = %d\n", stat);
    EG_close(context);
    return 1;
  }
  stat = EG_getBodyTopos(bodies[0], NULL, FACE, &j, &objs);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos = %d\n", stat);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  if (j != 6) {
    printf(" Number of Faces = %d!\n", j);
    EG_free(objs);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  for (i = 0; i < 6; i++) {
    patches[i].face = objs[i];
    patches[i].ni   = 0;
    patches[i].nj   = 0;
    patches[i].xyzs = NULL;
    stat = EG_getBodyTopos(bodies[0], patches[i].face, EDGE, &j, NULL);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos %d = %d\n", i+1, stat);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
    if (j != 4) {
      printf(" Face %d: number of Edges = %d\n", i+1, j);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
  }
  EG_free(objs);
  
  params[0] =  0.1;
  params[1] =  0.002;
  params[2] = 15.0;
  stat      = EG_makeTessBody(bodies[0], params, &tess);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTessBody = %d\n", stat);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  
  printf("\n");
  stat = EG_attributeRet(tess, ".tessType", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet tessType = %d\n", stat);
  } else if (aType == ATTRSTRING) {
    printf(" .tessType = %s\n", str);
  }
  
  stat = EG_attributeRet(tess, ".mixed", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet = %d\n", stat);
    EG_deleteObject(tess);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  if ((aType != ATTRINT) || (aLen != 6)) {
    printf(" aType = %d  aLen = %d\n", aType, aLen);
    EG_deleteObject(tess);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  printf(" .mixed    =");
  for (i = 0; i < 6; i++) printf(" %d", ints[i]);
  printf("\n\n");
  
  for (i = 1; i <= 6; i++) {
    stat = EG_getTessFace(tess, i, &len, &xyzs, &uvs, &ptype, &pindex,
                          &ntri, &tris, &tric);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getTessFace %d = %d\n", i, stat);
      for (j = 0; j < i-1; j++) EG_free(patches[j].xyzs);
      EG_deleteObject(tess);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
    printf(" Face %d: npts = %d  ntris = %d\n", i, len, ntri);
    if (ntri/2 != ints[i-1]) {
      printf(" Mismatch on Face %d!\n", i);
      for (j = 0; j < i-1; j++) EG_free(patches[j].xyzs);
      EG_deleteObject(tess);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
    for (k = j = 0; j < len; j++) {
      if (ptype[j] == 0) {
        if (k != 0) {
/*        printf(" Face %d: seg = %d\n", i, k+1);  */
          if (patches[i-1].ni == 0) {
            patches[i-1].ni = k+1;
          } if (patches[i-1].nj == 0) {
            patches[i-1].nj = k+1;
          } else {
            if (k+1 != patches[i-1].ni) {
              printf(" Mismatch Face %d -- first & third %d %d!\n",
                     i, k+1, patches[i-1].ni);
              for (m = 0; m < i-1; m++) EG_free(patches[m].xyzs);
              EG_deleteObject(tess);
              EG_deleteObject(bodies[0]);
              EG_close(context);
              return 1;
            }
          }
        }
        k = 0;
      } else if (ptype[j] < 0) {
/*      printf(" Face %d: seg = %d\n", i, k+1);  */
        if (k+1 != patches[i-1].nj) {
          printf(" Mismatch Face %d -- second & fourth %d %d!\n",
                 i, k+1, patches[i-1].nj);
          for (m = 0; m < i-1; m++) EG_free(patches[m].xyzs);
          EG_deleteObject(tess);
          EG_deleteObject(bodies[0]);
          EG_close(context);
          return 1;
        }
        break;
      }
      k++;
    }
    m = 3*patches[i-1].ni*patches[i-1].nj;
    patches[i-1].xyzs = (double *) EG_alloc(m*sizeof(double));
    if (patches[i-1].xyzs == NULL) {
      printf(" Malloc Error on Face %d -- %d %d!\n",
             i, patches[i-1].ni, patches[i-1].nj);
      for (m = 0; m < i-1; m++) EG_free(patches[m].xyzs);
      EG_deleteObject(tess);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
    for (j = 0; j < m; j++) patches[i-1].xyzs[j] = 0.0;
    
    /* fill in the bounds */
    for (m = k = j = 0; j < len; j++) {
      if (m == 0) {
        fillPatch(&patches[i-1], k, 0, &xyzs[3*j]);
      } else if (m == 1) {
        fillPatch(&patches[i-1], patches[i-1].ni-1, k, &xyzs[3*j]);
      } else if (m == 2) {
        fillPatch(&patches[i-1], patches[i-1].ni-k-1,
                                 patches[i-1].nj-1, &xyzs[3*j]);
      } else {
        if (ptype[j] < 0) break;
        fillPatch(&patches[i-1], 0, patches[i-1].nj-k-1, &xyzs[3*j]);
      }
      if (ptype[j] == 0) {
        if (k != 0) {
          m++;
        }
        k = 0;
      }
      k++;
    }
    
    /* simply fill in the interiors */
    fillInterior(&patches[i-1]);
  }
  EG_deleteObject(tess);
  printf("\n");
  
  /* make the new tessellation */
  stat = EG_initTessBody(bodies[0], &tess);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_initTessBody = %d\n", stat);
    for (j = 0; j < 6; j++) EG_free(patches[j].xyzs);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  
  /* create the Edge discretizations */
  for (i = 1; i <= 6; i++) {
    stat = EG_quadEdges(tess, &patches[i-1]);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadEdges %d = %d\n", i, stat);
      for (j = 0; j < 6; j++) EG_free(patches[j].xyzs);
      EG_deleteObject(tess);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
  }
  
  /* create the Face Tessellations */
  for (i = 1; i <= 6; i++) {
    stat = EG_quadFace(tess, &patches[i-1]);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadFace %d = %d\n", i, stat);
      for (j = 0; j < 6; j++) EG_free(patches[j].xyzs);
      EG_deleteObject(tess);
      EG_deleteObject(bodies[0]);
      EG_close(context);
      return 1;
    }
  }
  for (j = 0; j < 6; j++) EG_free(patches[j].xyzs);
  
  /* finish it off */
  stat = EG_finishTess(tess, params);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_finishTess = %d\n", stat);
    EG_deleteObject(tess);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  
  printf("\n");
  stat = EG_attributeRet(tess, ".tessType", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet tessType = %d\n", stat);
  } else if (aType == ATTRSTRING) {
    printf(" .tessType = %s\n", str);
  }
  
  stat = EG_attributeRet(tess, ".mixed", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet = %d\n", stat);
    EG_deleteObject(tess);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  if ((aType == ATTRINT) && (aLen == 6)) {
    printf(" .mixed    =");
    for (i = 0; i < 6; i++) printf(" %d", ints[i]);
    printf("\n");
  } else {
    printf(" aType = %d  aLen = %d\n", aType, aLen);
  }
  printf("\n");
  
  /* make the Model and write it out */
  bodies[1] = tess;
  stat      = EG_makeTopology(context, NULL, MODEL, 2, NULL, 1, bodies, NULL,
                              &model);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTopology = %d\n", stat);
    EG_deleteObject(tess);
    EG_deleteObject(bodies[0]);
    EG_close(context);
    return 1;
  }
  stat = EG_saveModel(model, "patch.egads");
  if (stat != EGADS_SUCCESS)
    printf(" EG_saveModel = %d\n", stat);
  EG_deleteObject(model);
  
  EG_close(context);
  return 0;
}
#endif
