#ifndef EGADSPATCH_H
#define EGADSPATCH_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Quad Patch function Header
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


typedef struct {
  ego    face;        /* associated Face object */
  int    ni;
  int    nj;
  double *xyzs;       /* 3d coordinates for the mesh (3*ni*nj in length)
                         index = i*nj + j */
  int    en[4][3];    /* eIndex, start, off -- internally set */
} quadPatch;


#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif
  
__ProtoExt__ int EG_quadEdges(ego Tess, quadPatch *patch);
/* where: Tess  -- the open Tessellation Object
          patch -- the quad patch

  note: This needs to be called for each "patch" before EG_quadFace is invoked
 */

__ProtoExt__ int EG_quadFace(ego Tess, quadPatch *patch);
/* where: Tess  -- the open Tessellation Object
          Patch -- the quad patch

  note: This needs to be called for each "patch" before EG_finishTess is called
        This also modifies the ".mixed" attribute on the tessellation object 
        to indicate which faces are quad faces
 */

#ifdef __cplusplus
}
#endif

#endif
