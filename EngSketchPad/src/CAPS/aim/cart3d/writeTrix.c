#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* libCart3d defaults to WORD_BIT 32 unless told otherwise
 * this makes that assumption explicit and suppresses the warning
 */

#define WORD_BIT 32

#include <geomStructures.h>
#define VALID_TRIX_FILE 1

#include "xddm.h"

#include "bodyTess.h"


int writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                  const char *const p_fileName, const int options);

int readSurfTrix(const char *const p_fileName, p_tsTriangulation *pp_config,
                 int *const p_nComps, const char *const p_compName,
                 const char *const p_vertDataNames,
                 const char *const p_triDataNames,
                 const int options);

int writeTrix(const char *fname, int nbody, ego *tess,
              /*@null@*/ p_tsXddm p_xddm, int nv, /*@null@*/ double ***dvar)
{
  int               status = EGADS_SUCCESS;
  int               i, j, k, rc, nvert=0, state, ibody, opts = 0;
  int               nedge, ntri, plen, iov, iot;
  double            *xyzs=NULL;
  const int         *trix, *tric, *ptype, *pindex;
  const double      *points, *uv;
  ego               ref;
  p_tsTriangulation p_surf=NULL;
  p_tsXmParent      p_p;
  int               nface, *tris=NULL;
  verTags           *vtags=NULL;

  rc = c3d_newTriangulation(&p_surf, 0, 1);
  if (rc != 0) {
    printf(" c3d_newTriangulation failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }

  if (p_xddm != NULL)
    for (i = 0; i < p_xddm->p_parent->nAttr; i++) {
      p_p = p_xddm->p_parent;
      if (strcmp(p_p->p_attr[i].p_name, "ID") == 0) {
        strcpy(p_surf->geomName, p_p->p_attr[i].p_value);
        break;
      }
    }

  if (nv != 0) {                  /* linearization lives at the verts */
    rc = c3d_allocVertData(&p_surf, nv);
    if (rc != 0) {
      printf(" c3d_allocVertData failed\n");
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nv; ++i) {
      p_surf->p_vertData[i].dim    = 3;
      p_surf->p_vertData[i].offset = i*3;
      p_surf->p_vertData[i].type   = VTK_Float64;
      p_surf->p_vertData[i].info   = TRIX_shapeLinearization;
    }
  }

  rc = c3d_allocTriData(&p_surf, 2);             /* 2 tags; body and FaceIDs */
  if (rc != 0) {
    printf(" c3d_allocTriData failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }

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

  //p_surf->infoCode += DP_VERTS_CODE;

  p_surf->nVerts = 0;
  p_surf->nTris  = 0;

  for (ibody = 0; ibody < nbody; ibody++) {
    status  = EG_statusTessBody(tess[ibody], &ref, &state, &nvert);
    if (status != EGADS_SUCCESS) goto cleanup;
    p_surf->nVerts += nvert;

    status = EG_getBodyTopos(ref, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (i = 1; i <= nface; i++) {
      status = EG_getTessFace(tess[ibody], i, &plen, &points, &uv, &ptype, &pindex,
                              &ntri, &trix, &tric);
      if (status != EGADS_SUCCESS) goto cleanup;
      p_surf->nTris  += ntri;
    }
  }

  rc = c3d_allocTriangulation(&p_surf);
  if (rc != 0) {
    printf(" c3d_allocTriangulation failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }

  iov = 0;
  iot = 0;
  for (ibody = 0; ibody < nbody; ibody++) {

    status = bodyTess(tess[ibody], &nface, &nedge, &nvert, &xyzs, &vtags, &ntri, &tris);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (i = 0; i < nvert; i++) {              /* save vert locations */
      p_surf->a_Verts[iov+i].x[0] = (float) xyzs[3*i  ];
      p_surf->a_Verts[iov+i].x[1] = (float) xyzs[3*i+1];
      p_surf->a_Verts[iov+i].x[2] = (float) xyzs[3*i+2];
      //surf->a_dpVerts[i].x[0] = xyzs[3*i  ];
      //surf->a_dpVerts[i].x[1] = xyzs[3*i+1];
      //surf->a_dpVerts[i].x[2] = xyzs[3*i+2];
    }

    for (i = 0; i < ntri; i++) {
      p_surf->a_Tris[iot+i].vtx[0]    = iov + tris[4*i  ] - 1;
      p_surf->a_Tris[iot+i].vtx[1]    = iov + tris[4*i+1] - 1;
      p_surf->a_Tris[iot+i].vtx[2]    = iov + tris[4*i+2] - 1;

      /* -- component numbers -- */
      p_surf->a_scalar0_t[iot+i]      = ibody;
      p_surf->a_scalar0_t[iot+ntri+i] = tris[4*i+3];
    }

    if ((dvar != NULL) && (p_xddm != NULL))
      for (i = 0; i < nv; i++) { /* actual sensitivities */
        strcpy(p_surf->p_vertData[i].name, p_xddm->a_v[i].p_id);
        for (j = k = 0; k < nvert; k++) {
          p_surf->a_scalar0[3*iov + j  + i*nvert*3] = dvar[ibody][i][3*k  ];
          p_surf->a_scalar0[3*iov + j+1+ i*nvert*3] = dvar[ibody][i][3*k+1];
          p_surf->a_scalar0[3*iov + j+2+ i*nvert*3] = dvar[ibody][i][3*k+2];
          j += 3;
        }
      }

    iov += nvert;
    iot += ntri;

    EG_free(xyzs); xyzs = NULL;
    EG_free(tris); tris = NULL;
    EG_free(vtags); vtags = NULL;
  }

  rc = writeSurfTrix(p_surf, 1, fname, opts);
  if (rc != 0) {
    printf(" io_writeSurfTrix failed: %d", rc);
    status = EGADS_WRITERR;
    goto cleanup;
  }

  status = EGADS_SUCCESS;

cleanup:
  if (p_surf != NULL) {
    c3d_freeTriangulation(p_surf, 0);
    free(p_surf); /* must use free */
  }

  EG_free(xyzs); xyzs = NULL;
  EG_free(tris); tris = NULL;
  EG_free(vtags); vtags = NULL;

  return status;
}


int readTrix(const char *fname, const char *tag, int *dim,
             double ***data_out)
{
  int               status, i, j, k, l, rc, rank, nvert, trixOpts, nComps, offset;
  double            **data=NULL;
  p_tsTriangulation p_surf;

  *data_out = NULL;

  trixOpts = 0;      /* or TRIX_VERBOSE */

  nComps = 0;        /* ...need to initialize to 0 before calling */
  rc     = readSurfTrix(fname, &p_surf, &nComps, tag, "ALL", "ALL",
                           trixOpts);
  if (rc == VALID_TRIX_FILE) {
    printf("trix_readSurf: %s is a valid TRIX (VTK) file\n", fname);
    return EGADS_READERR;
  } else if (rc < 0) {
    printf("trix_readSurf: %s is not a valid TRIX (VTK) file\n", fname);
    return EGADS_READERR;
  }

  data = (double**)EG_alloc(nComps*sizeof(double*));
  if (data == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < nComps; ++i) {
    data[i] = NULL;
  }

  for (i = 0; i < nComps; ++i) {
    data[i] = (double*)EG_alloc(p_surf[i].nVerts*sizeof(double));
    if (data[i] == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
  }


/*  printf("\n"); printf("    o Parsed %d components\n", nComps);  */
  for (i = 0; i < nComps; ++i) {
/*
    printf("       Comp %d Name: \"%s\"\n", i, p_surf[i].geomName);
    printf("              nVerts %d nTris %d nVertData %d nTriData %d\n",
           p_surf[i].nVerts, p_surf[i].nTris, p_surf[i].nVertData,
           p_surf[i].nTriData);  */
//    if (nvert != p_surf[i].nVerts) {
//      printf(" Wrong Number of Vertices!\n");
//      for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
//      free(p_surf);
//      return CAPS_MISMATCH;
//    }

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

    for (j = 0; j < p_surf[i].nVertData; ++j) {
      if (strcmp(tag, p_surf[i].p_vertData[j].name) == 0) {
        nvert  = p_surf[i].nVerts;
        rank   = p_surf[i].p_vertData[j].dim;
        offset = p_surf[i].p_vertData[j].offset;
        *dim = rank;
        for (k = 0; k < nvert; k++)
          for (l = 0; l < rank; l++)
            data[i][rank*k+l] = p_surf[i].a_scalar0[rank*k+l  + offset*nvert];
      }
    }
  }

  *data_out = data;
  status = EGADS_SUCCESS;

cleanup:
  for (i = 0; i < nComps; ++i) c3d_freeTriangulation(p_surf+i, 1);
  free(p_surf);
  if (status != EGADS_SUCCESS) {
    if (data != NULL) {
      for (i = 0; i < nComps; ++i) {
        EG_free(data[i]);
      }

    }
  }

  return status;
}
