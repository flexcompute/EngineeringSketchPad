#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* libCart3d defaults to WORD_BIT 32 unless told otherwise
 * this makes that assumption explicit and suppresses the warning
 */
#define WORD_BIT 32

#include "capsErrors.h"
#include "c3dio_lib.h"
#include "xddm.h"

#define VALID_TRIX_FILE 1


/* NOTE: needs rewriting for multiple bodies */


int writeTrix(const char *fname, /*@null@*/ p_tsXddm p_xddm, int ibody,
              int nvert, double *xyzs, int nv, /*@null@*/ double **dvar,
              int ntri, int *tris)
{
  int               i, j, k, rc, opts = 0;
  p_tsTriangulation p_surf;
  p_tsXmParent      p_p;

  c3d_newTriangulation(&p_surf, 0, 1);
  p_surf->nVerts = nvert;
  p_surf->nTris  = ntri;
  if (p_xddm != NULL)
    for (i = 0; i < p_xddm->p_parent->nAttr; i++) {
      p_p = p_xddm->p_parent;
      if (strcmp(p_p->p_attr[i].p_name, "ID") == 0) {
        strcpy(p_surf->geomName, p_p->p_attr[i].p_value);
        break;
      }
    }

  if (nv != 0) {                  /* linearization lives at the verts */
    c3d_allocVertData(&p_surf, nv);
    for (i = 0; i < nv; ++i) {
      p_surf->p_vertData[i].dim    = 3;
      p_surf->p_vertData[i].offset = i*3;
      p_surf->p_vertData[i].type   = VTK_Float64;
      p_surf->p_vertData[i].info   = TRIX_shapeLinearization;
    }
  }

  c3d_allocTriData(&p_surf, 2);             /* 2 tags; body and FaceIDs */
  strcpy(p_surf->p_triData[0].name, "IntersectComponents");
  p_surf->p_triData[0].dim    = 1;
  p_surf->p_triData[0].offset = 0;
  p_surf->p_triData[0].type   = VTK_Int16;
  p_surf->p_triData[0].info   = TRIX_componentTag;
  strcpy(p_surf->p_triData[1].name, "GMPtags");
  p_surf->p_triData[1].dim    = 1;
  p_surf->p_triData[1].offset = 1;
  p_surf->p_triData[1].type   = VTK_Int16;
  p_surf->p_triData[1].info   = TRIX_componentTag;

  rc = c3d_allocTriangulation(&p_surf);
  if (rc != 0) {
    printf(" c3d_allocTriangulation failed\n");
    return CAPS_SOURCEERR;
  }

  for (i = 0; i < nvert; i++) {              /* save vert locations */
    p_surf->a_Verts[i].x[0] = (float) xyzs[3*i  ];
    p_surf->a_Verts[i].x[1] = (float) xyzs[3*i+1];
    p_surf->a_Verts[i].x[2] = (float) xyzs[3*i+2];
  }

  for (i = 0; i < ntri; i++) {
    p_surf->a_Tris[i].vtx[0] = tris[4*i  ] - 1;
    p_surf->a_Tris[i].vtx[1] = tris[4*i+1] - 1;
    p_surf->a_Tris[i].vtx[2] = tris[4*i+2] - 1;
  }

  /* -- component numbers -- */
  for (i = 0; i < ntri; i++) {
    p_surf->a_scalar0_t[i]      = ibody;
    p_surf->a_scalar0_t[i+ntri] = tris[4*i+3];
  }

  if ((dvar != NULL) && (p_xddm != NULL))
    for (i = 0; i < nv; i++) { /* actual sensitivities */
      strcpy(p_surf->p_vertData[i].name, p_xddm->a_v[i].p_id);
      for (j = k = 0; k < nvert; k++) {
        p_surf->a_scalar0[j  + i*nvert*3] = dvar[i][3*k  ];
        p_surf->a_scalar0[j+1+ i*nvert*3] = dvar[i][3*k+1];
        p_surf->a_scalar0[j+2+ i*nvert*3] = dvar[i][3*k+2];
        j += 3;
      }
    }

  rc = io_writeSurfTrix(p_surf, 1, fname, opts);
  if (rc != 0) {
    printf(" io_writeSurfTrix failed: %d", rc);
    return CAPS_IOERR;
  }
  c3d_freeTriangulation(p_surf, 0);
  free(p_surf);

  return CAPS_SUCCESS;
}


int readTrix(const char *fname, const char *tag, int nvert, int dim,
             double *data)
{
  int               i, j, k, l, rc, trixOpts, nComps, offset;
  p_tsTriangulation p_surf;

  trixOpts = 0;      /* or TRIX_VERBOSE */

  nComps = 0;        /* ...need to initialize to 0 before calling */
  rc     = io_readSurfTrix(fname, &p_surf, &nComps, "ALL", "ALL", "ALL",
                           trixOpts);
  if (rc == VALID_TRIX_FILE) {
    printf("trix_readSurf: %s is a valid TRIX (VTK) file\n", fname);
    return CAPS_IOERR;
  } else if (rc < 0) {
    printf("trix_readSurf: %s is not a valid TRIX (VTK) file\n", fname);
    return CAPS_IOERR;
  }

/*  printf("\n"); printf("    o Parsed %d components\n", nComps);  */
  for (i = 0; i < nComps; ++i) {
/*
    printf("       Comp %d Name: \"%s\"\n", i, p_surf[i].geomName);
    printf("              nVerts %d nTris %d nVertData %d nTriData %d\n",
           p_surf[i].nVerts, p_surf[i].nTris, p_surf[i].nVertData,
           p_surf[i].nTriData);  */
    if (nvert != p_surf[i].nVerts) {
      printf(" Wrong Number of Vertices!\n");
      for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
      free(p_surf);
      return CAPS_MISMATCH;
    }

/*  printf("       Vertices X,Y,Z:\n");
    for (j=0; j<p_surf[i].nVerts; ++j) {
      printf("          %g %g %g\n",p_surf[i].a_Verts[j].x[X],
             p_surf[i].a_Verts[j].x[Y],p_surf[i].a_Verts[j].x[Z]);
    }
    printf("\n");
    printf("       Connectivity:\n");
    for (j=0; j<p_surf[i].nTris; ++j) {
      printf("          %d %d %d\n",p_surf[i].a_Tris[j].vtx[0],
             p_surf[i].a_Tris[j].vtx[1],p_surf[i].a_Tris[j].vtx[2]);
    }
    printf("\n");
*/
/*
    if (p_surf[i].nTriData > 0) {
      printf("TriData tags:\n");
      for (j=0; j<p_surf[i].nTriData; ++j) {
        printf("         Name = %s DIM = %d Offset = %d Info = %d\n",
               p_surf[i].p_triData[j].name,
               p_surf[i].p_triData[j].dim,
               p_surf[i].p_triData[j].offset,
               (int) p_surf[i].p_triData[j].info);
      }
    }
*/
/*
    if (p_surf[i].nVertData > 0) {
      printf("VertData tags:\n");
      for (j=0; j<p_surf[i].nVertData; ++j) {
        printf("         Name = %s DIM = %d Offset = %d Info = %d\n",
               p_surf[i].p_vertData[j].name,
               p_surf[i].p_vertData[j].dim,
               p_surf[i].p_vertData[j].offset,
               (int) p_surf[i].p_vertData[j].info);
      }
    }
 */
    for (j = 0; j < p_surf[i].nVertData; ++j)
      if (strcmp(tag, p_surf[i].p_vertData[j].name) == 0)
        if (dim == p_surf[i].p_vertData[j].dim) {
          offset = p_surf[i].p_vertData[j].offset;
          for (k = 0; k < nvert; k++)
            for (l = 0; l < dim; l++)
              data[dim*k+l] = p_surf->a_scalar0[dim*k+l  + offset*nvert];
          for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
          free(p_surf);
          return CAPS_SUCCESS;
        } else {
          printf(" DIM mismatch -- Name = %s DIM = %d (%d)\n", tag,
                 p_surf[i].p_vertData[j].dim, dim);
          for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
          free(p_surf);
          return CAPS_MISMATCH;
        }
  }

  for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
  free(p_surf);
  return CAPS_NOTFOUND;
}
