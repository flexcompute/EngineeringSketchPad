/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Internal UVmap Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "UVMAP_LIB.h"

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __PROTO_H_AND_D__
#undef __PROTO_H_AND_D__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#define __PROTO_H_AND_D__   extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#define __PROTO_H_AND_D__ extern
#endif

__PROTO_H_AND_D__ /*@null@*/ /*@out@*/ /*@only@*/ void *EG_alloc(size_t nbytes);
__PROTO_H_AND_D__ /*@null@*/ /*@only@*/ void *EG_reall(/*@null@*/ /*@only@*/
                                       /*@returned@*/ void *ptr, size_t nbytes);
__PROTO_H_AND_D__ void EG_free( /*@null@*/ /*@only@*/ void *pointer );

__PROTO_H_AND_D__ int  EG_inTriExact( double *t1, double *t2, double *t3,
                                      double *p, double *w );



__HOST_AND_DEVICE__ void uvmap_free(void *ptr)
{
  EG_free(ptr);
}


/*@-incondefs@*/
__HOST_AND_DEVICE__
/*@null@*/ void *uvmap_malloc(INT_ *err_flag, size_t size)
{
  void *ptr;

  *err_flag = EGADS_SUCCESS;
  ptr = EG_alloc(size);
  if (ptr == NULL) *err_flag = EGADS_MALLOC;
  return ptr;
}


__HOST_AND_DEVICE__
/*@null@*/ void *uvmap_realloc(INT_ *err_flag, void *ptr, size_t size)
{
  void *tmp;

  *err_flag = EGADS_SUCCESS;
  tmp = EG_reall(ptr, size);
  if (tmp == NULL) *err_flag = EGADS_MALLOC;
  return tmp;
}
/*@+incondefs@*/


__HOST_AND_DEVICE__
void uvmap_struct_free_index (INT_ index, uvmap_struct *uvmap_struct_ptr)
{
  // Free mapping data structure for surface at location index.

  uvmap_struct *struct_ptr;

  struct_ptr = &(uvmap_struct_ptr[0]);

  if (index < 0 || index >= struct_ptr->ndef)
    return;

  struct_ptr = &(uvmap_struct_ptr[index]);

  if (struct_ptr->mdef == 0)
    return;

  uvmap_free (struct_ptr->idibf);
  uvmap_free (struct_ptr->msrch);
  uvmap_free (struct_ptr->inibf);
  uvmap_free (struct_ptr->ibfibf);
  uvmap_free (struct_ptr->u);

  struct_ptr->idef = 0;
  struct_ptr->mdef = 0;
  struct_ptr->isrch = 0;
  struct_ptr->ibface = 0;
  struct_ptr->nbface = 0;

  struct_ptr->idibf = NULL;
  struct_ptr->msrch = NULL;
  struct_ptr->inibf = NULL;
  struct_ptr->ibfibf = NULL;
  struct_ptr->u = NULL;

  return;
}


__HOST_AND_DEVICE__
void uvmap_struct_free (void *ptr)
{
  // Free UV mapping data structure for all surfaces.

  uvmap_struct *struct_ptr;
  uvmap_struct *uvmap_struct_ptr;

  INT_ index, ndef;

  if (ptr == NULL)
    return;

  uvmap_struct_ptr = (uvmap_struct *) ptr;

  struct_ptr = &(uvmap_struct_ptr[0]);

  ndef = struct_ptr->ndef;

  for (index = 0; index < ndef; index++) {
    uvmap_struct_free_index (index, uvmap_struct_ptr);
  }

  uvmap_free (uvmap_struct_ptr);

  return;
}


__HOST_AND_DEVICE__ void
uvmap_struct_set_srch_data(INT_ index, INT_ isrch, INT_ ibface,
                           uvmap_struct *uvmap_struct_ptr)
{
  // Set searching data.

  uvmap_struct *struct_ptr;

  INT_ i;
  INT_ nsrch = 1000000;

  if (uvmap_struct_ptr == NULL)
    return;

  struct_ptr = &(uvmap_struct_ptr[index]);

  if (isrch <= 0 || isrch > nsrch) {

    isrch = 1;

    for (i = 1; i <= struct_ptr->nbface; i++) {
      struct_ptr->msrch[i] = 0;
    }
  }

  struct_ptr->isrch = isrch;
  struct_ptr->ibface = ibface;

  return;
}


__HOST_AND_DEVICE__ void uvmap_error_message(char *Text)
{
  fprintf (stderr, "%s\n", Text);
  fflush (stderr);
  return;
}


__HOST_AND_DEVICE__
INT_ uvmap_struct_find_entry(INT_ idef, INT_ *index,
                             uvmap_struct *uvmap_struct_ptr)
{
  // Find mapping structure index for surface idef.

  uvmap_struct *struct_ptr;

  INT_ ndef;

  if (uvmap_struct_ptr == NULL) {
    uvmap_error_message("*** ERROR 3502 : mapping surface structure not set ***");
    return 3502;
  }

  struct_ptr = &(uvmap_struct_ptr[0]);

  ndef = struct_ptr->ndef;

  *index = 0;

  do {

    struct_ptr = &(uvmap_struct_ptr[*index]);

    (*index)++;
  }
  while (*index < ndef && (struct_ptr->idef != idef || struct_ptr->mdef == 0));

  (*index)--;

  if (struct_ptr->idef != idef || struct_ptr->mdef == 0)
    *index = -1;

  return 0;
}


__HOST_AND_DEVICE__
INT_ uvmap_struct_get_entry(INT_ idef, INT_ *index, INT_ *isrch, INT_ *ibface,
                            INT_ *nbface, INT_ **idibf, INT_ **msrch,
                            INT_3D **inibf, INT_3D **ibfibf, DOUBLE_2D **u,
                            uvmap_struct *uvmap_struct_ptr)
{
  // Get mapping data structure data for surface idef.

  uvmap_struct *struct_ptr;

  INT_ status = 0;

  status = uvmap_struct_find_entry(idef, index, uvmap_struct_ptr);

  if (status == 0 && *index == -1) {
    uvmap_error_message ("*** ERROR 3501 : unable to find mapping surface ***");
    status = 3501;
  }

  if (status)
    return status;

  struct_ptr = &(uvmap_struct_ptr[*index]);

  *isrch = struct_ptr->isrch;
  *ibface = struct_ptr->ibface;
  *nbface = struct_ptr->nbface;

  *idibf = struct_ptr->idibf;
  *msrch = struct_ptr->msrch;
  *inibf = struct_ptr->inibf;
  *ibfibf = struct_ptr->ibfibf;
  *u = struct_ptr->u;

  return 0;
}


__HOST_AND_DEVICE__
INT_ uvmap_find_uv(INT_ idef, double u_[2], void *ptr, INT_ *local_idef,
                   INT_ *ibface, INT_ inode_[3], double s[3])
{
  // Find location of given UV coordinates.

  uvmap_struct *uvmap_struct_ptr;

  INT_      *idibf,  *msrch;
  INT_3D    *ibfibf, *inibf;
  DOUBLE_2D *u;

  INT_ found, ibface_save, index, inode, isrch, j, jbface, k, nbface;
  INT_ status = 0;

  double area[3], area_min, area_sum, du[3][2];
  double smin  = 1.0e-12;
  double smin2 = 0.1;

  uvmap_struct_ptr = (uvmap_struct *) ptr;

  if (uvmap_struct_ptr == NULL) {
    uvmap_error_message ("*** ERROR 3503 mapping surface structure not set ***");
    return 3503;
  }

  // get data from UV mapping data structure

  status = uvmap_struct_get_entry(idef, &index, &isrch, ibface, &nbface,
                                  &idibf, &msrch, &inibf, &ibfibf, &u,
                                  uvmap_struct_ptr);

  if (status)
    return status;

  // save starting tria-face index

  ibface_save = *ibface;

  // area coordinate search loop

  jbface = *ibface;

  do
  {
    *ibface = jbface;

    // set search flag

    msrch[*ibface] = isrch;

    // set UV coordinate deltas

    for (j = 0; j < 3; j++) {
      inode = inibf[*ibface][j];
      for (k = 0; k < 2; k++) {
        du[j][k] = u[inode][k] - u_[k];
      }
    }

    // set areas for area coordinates

    area[0] = du[1][0] * du[2][1] - du[1][1] * du[2][0];
    area[1] = du[2][0] * du[0][1] - du[2][1] * du[0][0];
    area[2] = du[0][0] * du[1][1] - du[0][1] * du[1][0];

    area_sum = area[0] + area[1] + area[2];

    // set minimum area

    area_min = MIN (area[0], area[1]);
    area_min = MIN (area[2], area_min);

    // check if tria-face contains the given UV coordinates

    found = (area_min + smin * area_sum >= 0.0) ? 1 : -2;

    // if not found then set next tria-face to search

    j = 0;

    while (j < 3 && found < -1) {

      if (area[j] < 0.0) {

        jbface = ibfibf[*ibface][j];

        found = (jbface > 0) ? ((msrch[jbface] != isrch) ? -1 : -2) : -3;
      }

      j++;
    }
  }
  while (found == -1);

  // if search is stuck at a boundary tria-face then check with a larger
  // tolerance

  if (found == -3 && smin2 > smin && area_min + smin2 * area_sum >= 0.0)
    found = 1;

  // if not found then do a brute force global search for containing tria-face

  if (found < 0)
  {
    found = -1;

    // loop over tria-faces

    *ibface = 0;

    do {

      (*ibface)++;

      if (msrch[*ibface] != isrch) {

        // set search flag

        msrch[*ibface] = isrch;

        // set UV coordinate deltas

        for (j = 0; j < 3; j++) {
          inode = inibf[*ibface][j];
          for (k = 0; k < 2; k++) {
            du[j][k] = u[inode][k] - u_[k];
          }
        }

        // set areas for area coordinates

        area[0] = du[1][0] * du[2][1] - du[1][1] * du[2][0];
        area[1] = du[2][0] * du[0][1] - du[2][1] * du[0][0];
        area[2] = du[0][0] * du[1][1] - du[0][1] * du[1][0];

        area_sum = area[0] + area[1] + area[2];

        // set minimum area

        area_min = MIN (area[0], area[1]);
        area_min = MIN (area[2], area_min);

        // check if tria-face contains the given UV coordinates
        // if search is at a boundary tria-face then also check with a larger
        // tolerance

        if (area_min + smin * area_sum >= 0.0 ||
            (smin2 > smin && area_min + smin2 * area_sum >= 0.0 &&
             ((area[0] == area_min && ibfibf[*ibface][0] <= 0) ||
              (area[1] == area_min && ibfibf[*ibface][1] <= 0) ||
              (area[2] == area_min && ibfibf[*ibface][2] <= 0))))
          found = 1;
      }
    }
    while (*ibface < nbface && found == -1);
  }

  // set search data in UV mapping data structure

  if (found == -1)
    *ibface = ibface_save;

  isrch++;

  uvmap_struct_set_srch_data(index, isrch, *ibface, uvmap_struct_ptr);

  // if found then for the containing tria-face set nodes/vertices, local
  // surface ID label, and shape-functions

  if (found == 1) {

    inode_[0] = inibf[*ibface][0];
    inode_[1] = inibf[*ibface][1];
    inode_[2] = inibf[*ibface][2];

    if (idibf)
      *local_idef = idibf[*ibface];
    else
      *local_idef = idef;

    s[0] = area[0] / (area[0] + area[1] + area[2]);
    s[1] = area[1] / (area[0] + area[1] + area[2]);
    s[2] = area[2] / (area[0] + area[1] + area[2]);
  }

  // if not found then set return value, default output argument values, and
  // reset tria-face starting index

  else {

    status = -1;

    *ibface = -1;

    inode_[0] = -1;
    inode_[1] = -1;
    inode_[2] = -1;

    *local_idef = -1;

    s[0] = -1.0;
    s[1] = -1.0;
    s[2] = -1.0;
  }

  return status;
}


__HOST_AND_DEVICE__ int
EG_uvmapFindUV(int idef, double uv[2], void *ptr, int *local_idef, int *itria,
               int ivertex[3], double s[3])
{
  // Find location of given UV coordinates.

  INT_ inode[3] = {0, 0, 0};
  INT_ ibface   = 0, local_idef_ = 0;
  int  status   = 0;

  // find location of given UV coordinates

  status = (int) uvmap_find_uv((INT_) idef, uv, ptr,
                               &local_idef_, &ibface, inode, s);

  *itria      = (int) ibface;
  *local_idef = (int) local_idef_;

  ivertex[0]  = (int) inode[0];
  ivertex[1]  = (int) inode[1];
  ivertex[2]  = (int) inode[2];

  // set return value

  if (status > 100000)
    status = EGADS_MALLOC;
  else if (status > 0)
    status = EGADS_UVMAP;
  else if (status == -1)
    status = EGADS_NOTFOUND;
  else
    status = EGADS_SUCCESS;

  return status;
}


__HOST_AND_DEVICE__ static void
EG_triRemap(int *trmap, int itri, int flag, int *verts, double *ws)
{
  int    i1, i2, i3, tris[3];
  double w[3];
  
  if (trmap         == NULL) return;
  if (trmap[itri-1] == 0)    return;
  
  tris[0]  = verts[0];
  tris[1]  = verts[1];
  tris[2]  = verts[2];
  w[0]     = ws[0];
  w[1]     = ws[1];
  w[2]     = ws[2];
  
  i1       =  trmap[itri-1]       & 3;
  i2       = (trmap[itri-1] >> 2) & 3;
  i3       = (trmap[itri-1] >> 4) & 3;
  
  if (flag == 1) {
    /* EGADS tri is source -- uvmap is destination */
    verts[0] = tris[i1-1];
    verts[1] = tris[i2-1];
    verts[2] = tris[i3-1];
    ws[0]    = w[i1-1];
    ws[1]    = w[i2-1];
    ws[2]    = w[i3-1];
  } else {
    /* uvmap is source -- EGADS tri is destination */
    verts[i1-1] = tris[0];
    verts[i2-1] = tris[1];
    verts[i3-1] = tris[2];
    ws[i1-1]    = w[0];
    ws[i2-1]    = w[1];
    ws[i3-1]    = w[2];
  }
}


/*               ********** Exposed Entry Points **********               */


/* return triangle containing the input UV
 *
 * uvmap = pointer to the internal uvmap structure
 * trmap = pointer to triangle map -- can be null
 * uv    = the input target UV (2 in length)
 * fID   = the returned Face ID
 * itri  = the returned index (1-bias) into tris for found location
 * verts = the 3 vertex indices for the triangle
 * ws    = the weights in the triangle for the vertices (3 in len)
 */

__HOST_AND_DEVICE__ int
EG_uvmapLocate(void *uvmap, int *trmap, double *uv, int *fID, int *itri,
               int *verts, double *ws)
{
  int          i, i1, i2, i3, stat, cls = 0;
  double       w[3], neg = 0.0;
  uvmap_struct *uvstruct;
  
  verts[0] = verts[1] = verts[2] = 0;
  ws[0]    = ws[1]    = ws[2]    = 0.0;
  stat = EG_uvmapFindUV(1, uv, uvmap, fID, itri, verts, ws);
  EG_triRemap(trmap, *itri, 0, verts, ws);
  if (stat != EGADS_NOTFOUND) return stat;
  
  /* need to extrapolate */
  uvstruct = (uvmap_struct *) uvmap;
  for (i = 1; i <= uvstruct->nbface; i++) {
    i1   = uvstruct->inibf[i][0];
    i2   = uvstruct->inibf[i][1];
    i3   = uvstruct->inibf[i][2];
    stat = EG_inTriExact(uvstruct->u[i1], uvstruct->u[i2], uvstruct->u[i3],
                         uv, w);
    if (stat == EGADS_SUCCESS) {
      *fID     = uvstruct->idibf[i];
      *itri    = i;
      ws[0]    = w[0];
      ws[1]    = w[1];
      ws[2]    = w[2];
      verts[0] = i1;
      verts[1] = i2;
      verts[2] = i3;
      EG_triRemap(trmap, i, 0, verts, ws);
/*    printf(" EGADS Info: Exact -> Ws = %le %le %le (EG_uvmapLocate)!\n",
             ws[0], ws[1], ws[2]);  */
      return EGADS_SUCCESS;
    }
    if (w[1] < w[0]) w[0] = w[1];
    if (w[2] < w[0]) w[0] = w[2];
    if (cls == 0) {
      cls = i;
      neg = w[0];
    } else {
      if (w[0] > neg) {
        cls = i;
        neg = w[0];
      }
    }
  }
  if (cls == 0) return EGADS_NOTFOUND;
  
  /* extrapolate */
  *fID  = uvstruct->idibf[cls];
  *itri = cls;
  i1    = uvstruct->inibf[cls][0];
  i2    = uvstruct->inibf[cls][1];
  i3    = uvstruct->inibf[cls][2];
  EG_inTriExact(uvstruct->u[i1], uvstruct->u[i2], uvstruct->u[i3], uv, ws);
  verts[0] = i1;
  verts[1] = i2;
  verts[2] = i3;
  EG_triRemap(trmap, cls, 0, verts, ws);
/*
  if (neg < -1.e-4) {
    printf(" EGADS Info: Extrapolate -> Ws = %le %le %le (EG_uvmapLocate)!\n",
           ws[0], ws[1], ws[2]);
    return EGADS_EXTRAPOL;
  }
 */
  return EGADS_SUCCESS;
}


/* return uvmap UV given index
 *
 * uvmap = pointer to the internal uvmap structure
 * index = index into uvmap
 * uv    = uv in uvmap
 */

__HOST_AND_DEVICE__ void
EG_getUVmap(void *uvmap, int index, double *uv)
{
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  uv[0]    = uvstruct->u[index][0];
  uv[1]    = uvstruct->u[index][1];
}


/* return uvmap UV given uv and fID
 *
 * uvmap = pointer to the internal uvmap structure
 * trmap = pointer to triangle map -- can be null
 * fuv   = the input Face UV (2 in length)
 * fuvs  = uvs for the Face
 * tris  = tris indices for fuvs
 * tbeg  = triangle start index
 * tend  = triangle end index
 * uv    = uv in uvmap
 */

__HOST_AND_DEVICE__ int
EG_uv2UVmap(void *uvmap, int *trmap, double *fuv, double *fuvs, int *tris,
            int tbeg, int tend, double *uv)
{
  int          i, j, stat, i1, i2, i3, verts[3], cls = 0;
  double       w[3], neg;
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  for (i = tbeg; i <= tend; i++) {
    j    = i - tbeg;
    i1   = tris[3*j  ] - 1;
    i2   = tris[3*j+1] - 1;
    i3   = tris[3*j+2] - 1;
    stat = EG_inTriExact(&fuvs[2*i1], &fuvs[2*i2], &fuvs[2*i3], fuv, w);
    if (stat == EGADS_SUCCESS) {
      verts[0] = uvstruct->inibf[i][0];
      verts[1] = uvstruct->inibf[i][1];
      verts[2] = uvstruct->inibf[i][2];
      EG_triRemap(trmap, i, 1, verts, w);
      i1    = verts[0];
      i2    = verts[1];
      i3    = verts[2];
      uv[0] = w[0]*uvstruct->u[i1][0] + w[1]*uvstruct->u[i2][0] +
              w[2]*uvstruct->u[i3][0];
      uv[1] = w[0]*uvstruct->u[i1][1] + w[1]*uvstruct->u[i2][1] +
              w[2]*uvstruct->u[i3][1];
      return EGADS_SUCCESS;
    }
    if (w[1] < w[0]) w[0] = w[1];
    if (w[2] < w[0]) w[0] = w[2];
    if (cls == 0) {
      cls = j;
      neg = w[0];
    } else {
      if (w[0] > neg) {
        cls = j;
        neg = w[0];
      }
    }
  }
  if (cls == 0) return EGADS_NOTFOUND;

  /* extrapolate */
  i1 = tris[3*cls  ] - 1;
  i2 = tris[3*cls+1] - 1;
  i3 = tris[3*cls+2] - 1;
  EG_inTriExact(&fuvs[2*i1], &fuvs[2*i2], &fuvs[2*i3], fuv, w);
  verts[0] = uvstruct->inibf[cls+tbeg][0];
  verts[1] = uvstruct->inibf[cls+tbeg][1];
  verts[2] = uvstruct->inibf[cls+tbeg][2];
  EG_triRemap(trmap, cls+tbeg, 1, verts, w);
  i1    = verts[0];
  i2    = verts[1];
  i3    = verts[2];
  uv[0] = w[0]*uvstruct->u[i1][0] + w[1]*uvstruct->u[i2][0] +
          w[2]*uvstruct->u[i3][0];
  uv[1] = w[0]*uvstruct->u[i1][1] + w[1]*uvstruct->u[i2][1] +
          w[2]*uvstruct->u[i3][1];
  
  return EGADS_SUCCESS;
}
