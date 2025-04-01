
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: geomStructures.c,v 1.16 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#include <string.h>
#include <time.h>

#include "memory_util.h"

#include "geomStructures.h"

/**
 * -- c3d_newTriangulation(): Create or resize, one or more
 *                            triangulation structs.  When creating a
 *                            single tri struct, set startComp=0 and
 *                            nComps=1, when creating an array set
 *                            startComp=0 and nComps to number of
 *                            elements, and when resizing the array,
 *                            set an startComp to the last original
 *                            element and nComps to the new total
 *                            size.
 */
int c3d_newTriangulation(p_tsTriangulation *pp_surf, const int startComp,
                         const int nComps)
{
  int i;
  p_tsTriangulation p_comp;

  if (0 >= nComps || 0 > startComp) return -1;
  if (nComps <= startComp) return -1;

  if (0 == startComp) {
    NEW_ARRAY(*pp_surf,tsTriangulation,nComps);
  }
  else {
    RESIZE_ARRAY(pp_surf,tsTriangulation,nComps);
  }

  for (i=startComp; i<nComps; ++i) { /* initialize */
    p_comp = (*pp_surf) + i;

    p_comp->geomName[0] = '\0';
    p_comp->nVerts = p_comp->nTris = p_comp->nVertScalars = p_comp->nTriScalars = 0;

    p_comp->infoCode = (unsigned int)0;
    p_comp->configBBox[0] = p_comp->configBBox[1] = p_comp->configBBox[2] = UNSET;
    p_comp->configBBox[3] = p_comp->configBBox[4] = p_comp->configBBox[5] = UNSET;
    p_comp->a_Tris = NULL; p_comp->a_Verts = NULL;

    p_comp->a_scalar0 = p_comp->a_scalar1 = p_comp->a_scalar2 = NULL;
    p_comp->a_scalar0_t = p_comp->a_scalar1_t = NULL;
    p_comp->a_U = p_comp->a_U_t = NULL;
    p_comp->a_area = NULL;
    p_comp->outwardNormals = TRUE;

    p_comp->nVertData = p_comp->nTriData = 0;
    p_comp->p_vertData = p_comp->p_triData = NULL;

    p_comp->a_dpVerts = NULL;
  }

  return 0;
}

/**
 * -- c3d_allocTriData(): allocate triData info array
 */
int c3d_allocTriData(p_tsTriangulation *pp_surf, const int nTriData)
{
  int i;
  const int idata = (*pp_surf)->nTriData;

  (*pp_surf)->nTriData += nTriData;

  if (NULL == (*pp_surf)->p_triData) {
    NEW_ARRAY((*pp_surf)->p_triData,tsSurfTriX,(*pp_surf)->nTriData);
  }
  else {
    RESIZE_ARRAY(&((*pp_surf)->p_triData),tsSurfTriX,(*pp_surf)->nTriData);
  }

  for (i=idata; i<(*pp_surf)->nTriData; ++i) {
    (*pp_surf)->p_triData[i].name[0] = '\0';
    (*pp_surf)->p_triData[i].dim    = UNSET;
    (*pp_surf)->p_triData[i].offset = UNSET;
    (*pp_surf)->p_triData[i].info   = TRIX_unset;
  }

  return 0;
}

/**
 * -- c3d_allocVertData(): allocate vertData info array
 */
int c3d_allocVertData(p_tsTriangulation *pp_surf, const int nVertData)
{
  int i;
  const int idata = (*pp_surf)->nVertData;

  (*pp_surf)->nVertData += nVertData;

  if (NULL == (*pp_surf)->p_vertData) {
    NEW_ARRAY((*pp_surf)->p_vertData,tsSurfTriX,(*pp_surf)->nVertData);
  }
  else {
    RESIZE_ARRAY(&((*pp_surf)->p_vertData),tsSurfTriX,(*pp_surf)->nVertData);
  }

  for (i=idata; i<(*pp_surf)->nVertData; ++i) {
    (*pp_surf)->p_vertData[i].name[0] = '\0';
    (*pp_surf)->p_vertData[i].dim    = UNSET;
    (*pp_surf)->p_vertData[i].offset = UNSET;
    (*pp_surf)->p_vertData[i].offset = TRIX_unset;
  }

  return 0;
}

/**
 * -- c3d_allocTriangulation(): allocate or re-size space for verts,
 *                              tris and any data associated with the
 *                              triangulation.  The function allocs
 *                              space based on the internal nVerts,
 *                              nTris, nVertData and nTriData
 *                              values.
 */
int c3d_allocTriangulation(p_tsTriangulation *pp_surf)
{
  int i;

  if (NULL == (*pp_surf)->a_Tris) {
    NEW_ARRAY((*pp_surf)->a_Tris, tsTri, (*pp_surf)->nTris);
    memset((*pp_surf)->a_Tris, 0, sizeof(tsTri)*((*pp_surf)->nTris));
  }
  else {
    RESIZE_ARRAY( &((*pp_surf)->a_Tris), tsTri, (*pp_surf)->nTris);
  }

  if (NULL == (*pp_surf)->a_Verts) {
    NEW_ARRAY((*pp_surf)->a_Verts, tsVertex, (*pp_surf)->nVerts);
    memset((*pp_surf)->a_Verts, 0, sizeof(tsVertex)*((*pp_surf)->nVerts));
  }
  else {
    RESIZE_ARRAY(&((*pp_surf)->a_Verts), tsVertex, (*pp_surf)->nVerts);
  }

  if ((*pp_surf)->infoCode & DP_VERTS_CODE ) {
    if (NULL == (*pp_surf)->a_dpVerts) {
      NEW_ARRAY((*pp_surf)->a_dpVerts, tsDPVertex, (*pp_surf)->nVerts);
      memset((*pp_surf)->a_dpVerts, 0, sizeof(tsDPVertex)*((*pp_surf)->nVerts));
    }
    else {
      RESIZE_ARRAY(&((*pp_surf)->a_dpVerts), tsDPVertex, (*pp_surf)->nVerts);
    }
  }

  if (0 != (*pp_surf)->nVertData) {
    if (0 == (*pp_surf)->nVertScalars) {
      for (i=0; i<(*pp_surf)->nVertData; ++i) {
        (*pp_surf)->nVertScalars += (*pp_surf)->p_vertData[i].dim;
      }
    }
    if (NULL == (*pp_surf)->a_scalar0) {
      NEW_ARRAY((*pp_surf)->a_scalar0, double,
                (*pp_surf)->nVerts*((*pp_surf)->nVertScalars));
      memset((*pp_surf)->a_scalar0, 0,
             sizeof(double)*((*pp_surf)->nVerts)*((*pp_surf)->nVertScalars));
    }
    else {
      RESIZE_ARRAY(&((*pp_surf)->a_scalar0), double,
                   (*pp_surf)->nVerts*((*pp_surf)->nVertScalars));
    }
  }

  if (0 != (*pp_surf)->nTriData) {
    if (0 == (*pp_surf)->nTriScalars) {
      for (i=0; i<(*pp_surf)->nTriData; ++i) {
        (*pp_surf)->nTriScalars += (*pp_surf)->p_triData[i].dim;
      }
    }
    if (NULL == (*pp_surf)->a_scalar0_t) {
      NEW_ARRAY((*pp_surf)->a_scalar0_t, double,
                (*pp_surf)->nTris*((*pp_surf)->nTriScalars));
      memset((*pp_surf)->a_scalar0_t, 0,
             sizeof(double)*((*pp_surf)->nTris)*((*pp_surf)->nTriScalars));
    }
    else {
      RESIZE_ARRAY(&((*pp_surf)->a_scalar0_t), double,
                   (*pp_surf)->nTris*((*pp_surf)->nTriScalars));
    }
  }

  return 0;
}

int c3d_resizeScalars(p_tsTriangulation *pp_surf)
{
  int i, hold;

  ASSERT(*pp_surf);

  hold = 0;
  for (i=0; i<(*pp_surf)->nVertData; ++i) {
    hold += (*pp_surf)->p_vertData[i].dim;
  }

  if ((*pp_surf)->nVertScalars < hold) {
    CONT("Cart3D LIB: Increasing TRI nVertScalars from %d to %d\n",
         (*pp_surf)->nVertScalars,hold);
  }

  if ((*pp_surf)->nVertScalars > 0) {
    ASSERT((*pp_surf)->a_scalar0);

    if ((*pp_surf)->nVertScalars < hold) {
      (*pp_surf)->nVertScalars = hold;
      RESIZE_ARRAY(&((*pp_surf)->a_scalar0), double,
                   (*pp_surf)->nVerts*((*pp_surf)->nVertScalars));
    }
  }
  else {
    ASSERT((*pp_surf)->nVertScalars == 0);
    if ((*pp_surf)->nVertScalars < hold) {
      (*pp_surf)->nVertScalars = hold;
      NEW_ARRAY((*pp_surf)->a_scalar0, double,
                (*pp_surf)->nVerts*((*pp_surf)->nVertScalars));
    }
  }

  if ((*pp_surf)->nTriScalars > 0) {

    ASSERT((*pp_surf)->a_scalar0_t);

    hold = 0;
    for (i=0; i<(*pp_surf)->nTriData; ++i) {
      hold += (*pp_surf)->p_triData[i].dim;
    }

    if ((*pp_surf)->nTriScalars < hold) {
      CONT("Cart3D LIB: Increasing TRI nTriScalars from %d to %d\n",
           (*pp_surf)->nTriScalars,hold);
      (*pp_surf)->nTriScalars = hold;
      RESIZE_ARRAY(&((*pp_surf)->a_scalar0_t), double,
                   (*pp_surf)->nTris*((*pp_surf)->nTriScalars));
    }
  }

  return 0;
}

int resizeTriangulation(p_tsTriangulation p_surf, const int *const p_nVerts,
                        const int *const p_nTris)
{
  /*
   * we do not allow changes in nVertScalars or nTriScalars because previous
   * components would need to be updated
   */

  if ( *p_nVerts > 0 &&  *p_nTris > 0 ) {
    p_surf->nVerts += *p_nVerts;
    ASSERT(NULL != p_surf->a_Verts);
    RESIZE_ARRAY(&p_surf->a_Verts, tsVertex, p_surf->nVerts);

    if (p_surf->infoCode & DP_VERTS_CODE) {
      RESIZE_ARRAY(&p_surf->a_dpVerts, tsDPVertex, p_surf->nVerts);
    }

    p_surf->nTris  += *p_nTris;
    ASSERT(NULL != p_surf->a_Tris);
    RESIZE_ARRAY(&p_surf->a_Tris, tsTri, p_surf->nTris);
  }

  if (p_surf->nVertScalars > 0) {
    if (NULL == p_surf->a_scalar0) {
      NEW_ARRAY(p_surf->a_scalar0, double, p_surf->nVerts*p_surf->nVertScalars);
      memset(p_surf->a_scalar0, 0, sizeof(double)*p_surf->nVerts*p_surf->nVertScalars);
    } else {
      RESIZE_ARRAY(p_surf->a_scalar0, double, p_surf->nVerts*p_surf->nVertScalars);
    }
  }

  if (p_surf->nTriScalars > 0) {
    if (NULL == p_surf->a_scalar0_t) {
      NEW_ARRAY(p_surf->a_scalar0_t, double, p_surf->nTris*p_surf->nTriScalars);
      memset(p_surf->a_scalar0_t, 0,
             sizeof(double)*p_surf->nTris*p_surf->nTriScalars);
    } else {
      RESIZE_ARRAY(p_surf->a_scalar0_t, double, p_surf->nTris*p_surf->nTriScalars);
    }
  }

  return 0;
}

p_tsTriangulation deepCopyTriangulation(const p_tsTriangulation p_surf) {
  p_tsTriangulation p_copy;

  NEW(p_copy, tsTriangulation);

  strncpy(p_copy->geomName, p_surf->geomName, FILENAME_LEN);

  p_copy->nVerts       = p_surf->nVerts;
  p_copy->nTris        = p_surf->nTris;
  p_copy->nVertScalars = p_surf->nVertScalars;
  p_copy->nTriScalars  = p_surf->nTriScalars;

  p_copy->configBBox[0] = p_surf->configBBox[0];
  p_copy->configBBox[1] = p_surf->configBBox[1];
  p_copy->configBBox[2] = p_surf->configBBox[2];
  p_copy->configBBox[3] = p_surf->configBBox[3];
  p_copy->configBBox[4] = p_surf->configBBox[4];
  p_copy->configBBox[5] = p_surf->configBBox[5];

  NEW_ARRAY(p_copy->a_Tris, tsTri, p_copy->nTris);
  memcpy(p_copy->a_Tris, p_surf->a_Tris, sizeof(tsTri)*p_copy->nTris);

  NEW_ARRAY(p_copy->a_Verts, tsVertex, p_copy->nVerts);
  memcpy(p_copy->a_Verts, p_surf->a_Verts, sizeof(tsVertex)*p_copy->nVerts);

  if (p_copy->nVertScalars) {
    NEW_ARRAY(p_copy->a_scalar0, double, p_copy->nVerts*p_copy->nVertScalars);
    memcpy(p_copy->a_scalar0, p_surf->a_scalar0,
           sizeof(double)*p_copy->nVerts*p_copy->nVertScalars);
  } else {
    p_copy->a_scalar0 = NULL;
  }

  if (p_copy->nTriScalars) {
    NEW_ARRAY(p_copy->a_scalar0_t, double, p_copy->nTris*p_copy->nTriScalars);
    memcpy(p_copy->a_scalar0_t, p_surf->a_scalar0_t,
           sizeof(double)*p_copy->nTris*p_copy->nTriScalars);
  } else {
    p_copy->a_scalar0_t = NULL;
  }

  return(p_copy);
}

/**
 * -- c3d_freeTriangulation(): destructor
 */
void c3d_freeTriangulation(p_tsTriangulation p_surf, const int verbose) {
  if (p_surf->a_Tris) {
    if (0==p_surf->nTris) {
      printf("freeTriangulation: nTris is zero??\n");
    } else {
      FREE(p_surf->a_Tris, sizeof(tsTri)*p_surf->nTris);
      if (verbose) printf("Freed p_surf->a_Tris\n");
    }
  }

  if (p_surf->a_Verts) {
    if (0==p_surf->nVerts) {
      printf("freeTriangulation: nVerts is zero??\n");
    } else {
      FREE(p_surf->a_Verts, sizeof(tsVertex)*p_surf->nVerts);
      if (verbose) printf("Freed p_surf->a_Verts\n");
    }
  }

  if (p_surf->a_dpVerts) {
    if (0==p_surf->nVerts) {
      printf("freeTriangulation: nVerts is zero??\n");
    } else {
      FREE(p_surf->a_dpVerts, sizeof(tsDPVertex)*p_surf->nVerts);
      if (verbose) printf("Freed p_surf->a_dpVerts\n");
    }
  }

  if (p_surf->a_scalar0) {
    if (p_surf->nVertScalars > 0) {
      FREE(p_surf->a_scalar0, sizeof(double)*p_surf->nVerts*p_surf->nVertScalars);
      if (verbose) printf("Freed p_surf->a_scalar0\n");
    } else {
      printf("freeTriangulation: check nVertScalars value??\n");
    }
  }

  if (p_surf->a_scalar1)
    printf("freeTriangulation: a_scalar1 is NOT NULL ??\n");

  if (p_surf->a_scalar2)
    printf("freeTriangulation: a_scalar2 is NOT NULL ??\n");

  if (p_surf->a_U)
    printf("freeTriangulation: a_U is NOT NULL ??\n");

  if (p_surf->a_scalar0_t) {
    if (p_surf->nTriScalars > 0) {
      FREE(p_surf->a_scalar0_t, sizeof(double)*p_surf->nTris*p_surf->nTriScalars);
      if (verbose) printf("Freed p_surf->a_scalar0_t\n");
    } else {
      printf("freeTriangulation: check nTriScalars value??\n");
    }
  }

  if (p_surf->a_scalar1_t)
    printf("freeTriangulation: a_scalar1_t is NOT NULL ??\n");

  if (p_surf->a_U_t)
    printf("freeTriangulation: a_U_t is NOT NULL ??\n");

  if (p_surf->a_area)
    printf("freeTriangulation: a_area is NOT NULL ??\n");

  p_surf->nVerts = p_surf->nTris = p_surf->nVertScalars = p_surf->nTriScalars = 0;

  p_surf->geomName[0] = '\0';

  if (p_surf->p_vertData) {
    if (0==p_surf->nVertData) {
      printf("freeTriangulation: nVertData is zero??\n");
    } else {
      FREE(p_surf->p_vertData, sizeof(tsSurfTriX)*p_surf->nVertData);
      if (verbose) printf("Freed p_surf->p_vertData\n");
    }
  }
  p_surf->nVertData = 0;

  if (p_surf->p_triData) {
    if (0==p_surf->nTriData) {
      printf("freeTriangulation: nTriData is zero??\n");
    } else {
      FREE(p_surf->p_triData, sizeof(tsSurfTriX)*p_surf->nTriData);
      if (verbose) printf("Freed p_surf->p_triData\n");
    }
  }
  p_surf->nTriData = 0;

  return;
}
