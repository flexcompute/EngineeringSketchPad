/* Common data transfer functions to be used within an AIM */

#include "math.h"

#include "aimUtil.h"

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


/* initialize capsBodyDiscr pointer */
void aim_initBodyDiscr(capsBodyDiscr *discBody)
{
  if (discBody == NULL) return;

  discBody->tess         = NULL;
  discBody->nElems       = 0;
  discBody->elems        = NULL;
  discBody->gIndices     = NULL;
  discBody->dIndices     = NULL;
  discBody->poly         = NULL;
  discBody->globalOffset = 0;
}


/* create a nodal triangle element type */
int aim_nodalTriangleType(capsEleType *eletype)
{
  int status;

  if (eletype == NULL) return CAPS_NULLOBJ;

  status = EGADS_MALLOC;

  eletype->nref  = 3;
  eletype->ndata = 0;            /* data at geom reference positions */
  eletype->ntri  = 1;
  eletype->nmat  = 0;            /* match points at geom ref positions */
  eletype->gst   = NULL;
  eletype->dst   = NULL;
  eletype->matst = NULL;
  eletype->tris  = NULL;

  /*  t
      ^
      |
      3
      | \
      |   \
      |     \
      1-------2---> s
  */

  eletype->tris    = (int *) EG_alloc(3*eletype->ntri*sizeof(int));
  if (eletype->tris == NULL) goto cleanup;
  eletype->tris[0] = 1;
  eletype->tris[1] = 2;
  eletype->tris[2] = 3;

  eletype->gst    = (double *) EG_alloc(2*eletype->nref*sizeof(double));
  if (eletype->gst == NULL) goto cleanup;
  eletype->gst[0] = 0.0;    /* s = 0, t = 0 */
  eletype->gst[1] = 0.0;
  eletype->gst[2] = 1.0;    /* s = 1, t = 0 */
  eletype->gst[3] = 0.0;
  eletype->gst[4] = 0.0;    /* s = 0, t = 1 */
  eletype->gst[5] = 1.0;

  return CAPS_SUCCESS;

cleanup:
  printf(" aimTransferUtil/aim_NodalTriangleType: malloc error!\n");

  EG_free(eletype->gst);
  eletype->gst = NULL;
  EG_free(eletype->dst);
  eletype->dst = NULL;
  EG_free(eletype->matst);
  eletype->matst = NULL;
  EG_free(eletype->tris);
  eletype->tris = NULL;

  return status;
}


/* create a nodal quadrilateral element type */
int aim_nodalQuadType(capsEleType *eletype)
{
  int status;

  if (eletype == NULL) return CAPS_NULLOBJ;

  status = EGADS_MALLOC;

  eletype->nref  = 4;
  eletype->ndata = 0;            /* data at geom reference positions */
  eletype->ntri  = 2;
  eletype->nmat  = 0;            /* match points at geom ref positions */
  eletype->gst   = NULL;
  eletype->dst   = NULL;
  eletype->matst = NULL;
  eletype->tris  = NULL;

  /*  t
      ^
      |
      4-------3
      |     / |
      |   /   |
      | /     |
      1-------2---> s
  */

  eletype->tris    = (int *) EG_alloc(3*eletype->ntri*sizeof(int));
  if (eletype->tris == NULL) goto cleanup;
  eletype->tris[0] = 1;
  eletype->tris[1] = 2;
  eletype->tris[2] = 3;

  eletype->tris[3] = 3;
  eletype->tris[4] = 4;
  eletype->tris[5] = 1;

  eletype->gst    = (double *) EG_alloc(2*eletype->nref*sizeof(double));
  if (eletype->gst == NULL) goto cleanup;
  eletype->gst[0] = 0.0;    /* s = 0, t = 0 */
  eletype->gst[1] = 0.0;
  eletype->gst[2] = 1.0;    /* s = 1, t = 0 */
  eletype->gst[3] = 0.0;
  eletype->gst[4] = 1.0;    /* s = 1, t = 1 */
  eletype->gst[5] = 1.0;
  eletype->gst[6] = 0.0;    /* s = 0, t = 1 */
  eletype->gst[7] = 1.0;

  return CAPS_SUCCESS;

cleanup:
  printf(" aimTransferUtil/aim_nodalQuadType: malloc error!\n");

  EG_free(eletype->gst);
  eletype->gst = NULL;
  EG_free(eletype->dst);
  eletype->dst = NULL;
  EG_free(eletype->matst);
  eletype->matst = NULL;
  EG_free(eletype->tris);
  eletype->tris = NULL;

  return status;
}


/* Create element type for cell centered triangle */
int aim_cellTriangleType(capsEleType *eletype)
{
  int status;

  if (eletype == NULL) return CAPS_NULLOBJ;

  status = EGADS_MALLOC;

  eletype->nref  = 3;
  eletype->ndata = 1;            /* data at 1 reference positions */
  eletype->ntri  = 1;
  eletype->nmat  = 1;            /* match points at 1 ref positions */
  eletype->gst   = NULL;
  eletype->dst   = NULL;
  eletype->matst = NULL;
  eletype->tris  = NULL;

  /*  t
      ^
      |
      3
      | \
      |   \
      |     \
      1-------2---> s
  */

  eletype->tris    = (int *) EG_alloc(3*eletype->ntri*sizeof(int));
  if (eletype->tris == NULL) goto cleanup;
  eletype->tris[0] = 1;
  eletype->tris[1] = 2;
  eletype->tris[2] = 3;

  eletype->gst    = (double *) EG_alloc(2*eletype->nref*sizeof(double));
  if (eletype->gst == NULL) goto cleanup;
  eletype->gst[0] = 0.0;    /* s = 0, t = 0 */
  eletype->gst[1] = 0.0;
  eletype->gst[2] = 1.0;    /* s = 1, t = 0 */
  eletype->gst[3] = 0.0;
  eletype->gst[4] = 0.0;    /* s = 0, t = 1 */
  eletype->gst[5] = 1.0;

  /* reference location of the data */
  eletype->dst    = (double *) EG_alloc(2*eletype->ndata*sizeof(double));
  if (eletype->dst == NULL) goto cleanup;
  eletype->dst[0] = 1./3.;
  eletype->dst[1] = 1./3.;

  eletype->matst    = (double *) EG_alloc(2*eletype->ndata*sizeof(double));
  if (eletype->matst == NULL) goto cleanup;
  eletype->matst[0] = 1./3.;
  eletype->matst[1] = 1./3.;

  return CAPS_SUCCESS;

cleanup:
  printf(" aimTransferUtil/aim_CellTriangleType: malloc error!\n");

  EG_free(eletype->gst);
  eletype->gst = NULL;
  EG_free(eletype->dst);
  eletype->dst = NULL;
  EG_free(eletype->matst);
  eletype->matst = NULL;
  EG_free(eletype->tris);
  eletype->tris = NULL;

  return status;
}


/* Create element type for cell centered quadrilateral */
int aim_cellQuadType(capsEleType *eletype)
{
  int status;

  if (eletype == NULL) return CAPS_NULLOBJ;

  status = EGADS_MALLOC;

  eletype->nref  = 4;
  eletype->ndata = 1;            /* data at 1 reference positions */
  eletype->ntri  = 2;
  eletype->nmat  = 1;            /* match points at 1 ref positions */
  eletype->gst   = NULL;
  eletype->dst   = NULL;
  eletype->matst = NULL;
  eletype->tris  = NULL;

  /*  t
      ^
      |
      4-------3
      |     / |
      |   /   |
      | /     |
      1-------2---> s
  */

  eletype->tris    = (int *) EG_alloc(3*eletype->ntri*sizeof(int));
  if (eletype->tris == NULL) goto cleanup;
  eletype->tris[0] = 1;
  eletype->tris[1] = 2;
  eletype->tris[2] = 3;

  eletype->tris[3] = 3;
  eletype->tris[4] = 4;
  eletype->tris[5] = 1;

  eletype->gst    = (double *) EG_alloc(2*eletype->nref*sizeof(double));
  if (eletype->gst == NULL) goto cleanup;
  eletype->gst[0] = 0.0;    /* s = 0, t = 0 */
  eletype->gst[1] = 0.0;
  eletype->gst[2] = 1.0;    /* s = 1, t = 0 */
  eletype->gst[3] = 0.0;
  eletype->gst[4] = 1.0;    /* s = 1, t = 1 */
  eletype->gst[5] = 1.0;
  eletype->gst[6] = 0.0;    /* s = 0, t = 1 */
  eletype->gst[7] = 1.0;

  /* reference location of the data */
  eletype->dst    = (double *) EG_alloc(2*eletype->ndata*sizeof(double));
  if (eletype->dst == NULL) goto cleanup;
  eletype->dst[0] = 1./2.;
  eletype->dst[1] = 1./2.;

  eletype->matst    = (double *) EG_alloc(2*eletype->ndata*sizeof(double));
  if (eletype->matst == NULL) goto cleanup;
  eletype->matst[0] = 1./2.;
  eletype->matst[1] = 1./2.;

  return CAPS_SUCCESS;

cleanup:
  printf(" aimTransferUtil/aim_CellTriangleType: malloc error!\n");

  EG_free(eletype->gst);
  eletype->gst = NULL;
  EG_free(eletype->dst);
  eletype->dst = NULL;
  EG_free(eletype->matst);
  eletype->matst = NULL;
  EG_free(eletype->tris);
  eletype->tris = NULL;

  return status;
}


/* Use newton's method to solve for the quadrilaterly reference coordiantes */
static int
invEvaluationQuad(const double uvs[], /* (in) uvs that support the geom
                                         (2*npts in length) */
                  const double uv[],  /* (in) the uv position to get st */
                  const int    in[],  /* (in) grid indices */
                        double st[])  /* (out) weighting */
{

  int    i;
  double idet, delta, d, du[2], dv[2], duv[2], dst[2], uvx[2];

  delta  = 100.0;
  for (i = 0; i < 20; i++) {
    uvx[0] = (1.0-st[0])*((1.0-st[1])*uvs[2*in[0]  ]  +
                               st[1] *uvs[2*in[3]  ]) +
                  st[0] *((1.0-st[1])*uvs[2*in[1]  ]  +
                               st[1] *uvs[2*in[2]  ]);
    uvx[1] = (1.0-st[0])*((1.0-st[1])*uvs[2*in[0]+1]  +
                               st[1] *uvs[2*in[3]+1]) +
                  st[0] *((1.0-st[1])*uvs[2*in[1]+1]  +
                               st[1] *uvs[2*in[2]+1]);
    du[0]  = (1.0-st[1])*(uvs[2*in[1]  ] - uvs[2*in[0]  ]) +
                  st[1] *(uvs[2*in[2]  ] - uvs[2*in[3]  ]);
    du[1]  = (1.0-st[0])*(uvs[2*in[3]  ] - uvs[2*in[0]  ]) +
                  st[0] *(uvs[2*in[2]  ] - uvs[2*in[1]  ]);
    dv[0]  = (1.0-st[1])*(uvs[2*in[1]+1] - uvs[2*in[0]+1]) +
                  st[1] *(uvs[2*in[2]+1] - uvs[2*in[3]+1]);
    dv[1]  = (1.0-st[0])*(uvs[2*in[3]+1] - uvs[2*in[0]+1]) +
                  st[0] *(uvs[2*in[2]+1] - uvs[2*in[1]+1]);
    duv[0] = uv[0] - uvx[0];
    duv[1] = uv[1] - uvx[1];
    idet   = du[0]*dv[1] - du[1]*dv[0];
    if (idet == 0.0) break;
    dst[0] = (dv[1]*duv[0] - du[1]*duv[1])/idet;
    dst[1] = (du[0]*duv[1] - dv[0]*duv[0])/idet;
    d      = sqrt(dst[0]*dst[0] + dst[1]*dst[1]);
    if (d >= delta) break;
    delta  = d;
    st[0] += dst[0];
    st[1] += dst[1];
    if (delta < 1.e-8) break;
  }

  if (delta < 1.e-8) return CAPS_SUCCESS;

  return CAPS_NOTFOUND;
}


/* locate an element within the trianglution of an element */
int aim_locateElement(capsDiscr *discr, double *params, double *param,
                      int *bIndex, int *eIndex, double *bary)
{
  int         i, j, k, in[4], itri[3], status;
  int         ib, ibsmall = 0, itsmall = 0, ismall = 0;
  double      we[3], w, smallw = -1.e300;
  capsBodyDiscr *discBody=NULL;
  capsEleType *eletype=NULL;

  if (discr == NULL) return CAPS_NULLOBJ;

  for (ib = 0; ib < discr->nBodys; ib++) {
    discBody = &discr->bodys[ib];
    for (i = 0; i < discBody->nElems; i++) {
      eletype = discr->types + discBody->elems[i].tIndex-1;
      for (j = 0; j < eletype->ntri; j++) {

        itri[0] = eletype->tris[3*j+0]-1;
        itri[1] = eletype->tris[3*j+1]-1;
        itri[2] = eletype->tris[3*j+2]-1;
        in[0]   = discBody->elems[i].gIndices[2*itri[0]] - 1;
        in[1]   = discBody->elems[i].gIndices[2*itri[1]] - 1;
        in[2]   = discBody->elems[i].gIndices[2*itri[2]] - 1;
        status  = EG_inTriExact(&params[2*in[0]], &params[2*in[1]],
                                &params[2*in[2]], param, we);

        if (status == EGADS_SUCCESS) {
          *bIndex = ib+1;
          *eIndex = i+1;
          /* interpolate reference coordinates to bary */
          for (k = 0; k < 2; k++)
            bary[k] = eletype->gst[2*itri[0]+k]*we[0] +
                      eletype->gst[2*itri[1]+k]*we[1] +
                      eletype->gst[2*itri[2]+k]*we[2];

          /* Linear quad */
          if (discr->types[discBody->elems[i].tIndex-1].nref == 4) {
            in[0] = discBody->elems[i].gIndices[0] - 1;
            in[1] = discBody->elems[i].gIndices[2] - 1;
            in[2] = discBody->elems[i].gIndices[4] - 1;
            in[3] = discBody->elems[i].gIndices[6] - 1;
            invEvaluationQuad(params, param, in, bary);
          }

          return CAPS_SUCCESS;
        }

                       w = we[0];
        if (we[1] < w) w = we[1];
        if (we[2] < w) w = we[2];
        if (w > smallw) {
          ismall = i+1;
          itsmall = j;
          ibsmall = ib+1;
          smallw = w;
        }
      }
    }
  }

  /* must extrapolate! */
  if (ismall == 0) return CAPS_NOTFOUND;

  discBody = &discr->bodys[ibsmall-1];
  eletype = discr->types + discBody->elems[ismall-1].tIndex-1;

  itri[0] = eletype->tris[3*itsmall+0]-1;
  itri[1] = eletype->tris[3*itsmall+1]-1;
  itri[2] = eletype->tris[3*itsmall+2]-1;
  in[0]   = discBody->elems[ismall-1].gIndices[2*itri[0]] - 1;
  in[1]   = discBody->elems[ismall-1].gIndices[2*itri[1]] - 1;
  in[2]   = discBody->elems[ismall-1].gIndices[2*itri[2]] - 1;
  EG_inTriExact(&params[2*in[0]], &params[2*in[1]], &params[2*in[2]], param, we);

  *bIndex = ibsmall;
  *eIndex = ismall;
  /* interpolate reference coordinates to bary */
  for (k = 0; k < 2; k++)
    bary[k] = eletype->gst[2*itri[0]+k]*we[0] +
              eletype->gst[2*itri[1]+k]*we[1] +
              eletype->gst[2*itri[2]+k]*we[2];

  return CAPS_SUCCESS;
}


/* Interpolation for a linear triangular element */
static int interpolation_LinearTriangle(capsDiscr *discr, int bIndex, int eIndex,
                                        double *bary, int rank, double *data,
                                        double *result)
{
  int    in[3], i;
  double we[3];
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if ((discr->types[discBody->elems[eIndex-1].tIndex-1].nref  != 3) ||
      (discr->types[discBody->elems[eIndex-1].tIndex-1].ndata != 0))
    return CAPS_BADTYPE;

  /* Linear triangle */
  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  for (i = 0; i < rank; i++)
    result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2];

  return CAPS_SUCCESS;
}


/* Interpolation for a linear quadrilateral element */
static int
interpolation_LinearQuad(capsDiscr *discr, int bIndex, int eIndex, double *bary,
                         int rank, double *data, double *result)
{
  int    in[4], i;
  double we[2];
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if ((discr->types[discBody->elems[eIndex-1].tIndex-1].nref  != 4) ||
      (discr->types[discBody->elems[eIndex-1].tIndex-1].ndata != 0))
    return CAPS_BADTYPE;

  we[0] = bary[0];
  we[1] = bary[1];

  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  in[3] = discBody->elems[eIndex-1].gIndices[6] - 1;
  for (i = 0; i < rank; i++)
    result[i] = (1.0-we[0])*( (1.0-we[1])*data[rank*in[0]+i] +
                                   we[1] *data[rank*in[3]+i] ) +
                     we[0] *( (1.0-we[1])*data[rank*in[1]+i] +
                                   we[1] *data[rank*in[2]+i] );

  return CAPS_SUCCESS;
}


/* Interpolation for cell centered data */
static int
interpolation_CellCenter(capsDiscr *discr, int bIndex, int eIndex, int rank,
                         double *data, double *result)
{
  int ind, i;
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].ndata != 1)
    return CAPS_BADTYPE;

  ind = discBody->elems[eIndex-1].dIndices[0] - 1;
  for (i = 0; i < rank; i++)
    result[i] = data[rank*ind+i];

  return CAPS_SUCCESS;
}


/* Interpolation selector */
int aim_interpolation(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                      double *bary, int rank, double *data, double *result)
{
  int status;
  capsBodyDiscr *discBody=NULL;
  capsEleType *eletype=NULL;

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" aimTransferUtil/aim_interpolation: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" aimTransferUtil/aim_interpolation: name = %s, eIndex = %d [1-%d]!\n",
           name, eIndex, discBody->nElems);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  eletype = &discr->types[discBody->elems[eIndex-1].tIndex-1];
  if (eletype->ndata == 0 ) {
    /* data at the reference locations */
    if (eletype->nref == 3) {
      /* Linear triangle */
      status = interpolation_LinearTriangle(discr, bIndex, eIndex, bary, rank,
                                            data, result);
    } else if (eletype->nref == 4) {
      /* Linear quadrilateral */
      status = interpolation_LinearQuad(discr, bIndex, eIndex, bary, rank, data,
                                        result);
    } else {
      printf(" aimTransferUtil/aim_interpolation: name = %s, eIndex = %d [1-%d], nref not recognized!\n",
             name, eIndex, discBody->nElems);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  } else if (eletype->ndata == 1 ) {
    /* data at cell center */
    status = interpolation_CellCenter(discr, bIndex, eIndex, rank, data, result);
  } else {
    printf(" aimTransferUtil/aim_interpolation: name = %s, tIndex = %d, ndata = %d. Only supports ndata = 0 or ndata = 1!\n",
           name, discBody->elems[eIndex-1].tIndex, eletype->ndata);
    status = CAPS_BADTYPE;
    goto cleanup;
  }

cleanup:
  if (status != CAPS_SUCCESS)
    printf("Premature exit in transferUtils interpolation name = %s, status = %d\n",
           name, status);

  return status;
}


/* Interpolation Bar for a linear triangular element */
static int
interpolateBar_LinearTriangle(capsDiscr *discr, int bIndex, int eIndex,
                              double *bary, int rank,
                              double *r_bar, double *d_bar)
{
  int    in[4], i;
  double we[3];
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 3)
    return CAPS_BADTYPE;

  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  for (i = 0; i < rank; i++) {
/*  result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2]; */
    d_bar[rank*in[0]+i] += we[0]*r_bar[i];
    d_bar[rank*in[1]+i] += we[1]*r_bar[i];
    d_bar[rank*in[2]+i] += we[2]*r_bar[i];
  }

  return CAPS_SUCCESS;
}


/* Interpolation Bar for a linear quadrilateral element */
static int
interpolateBar_LinearQuad(capsDiscr *discr, int bIndex, int eIndex, double *bary,
                                     int rank, double *r_bar, double *d_bar)
{
  int    in[4], i;
  double we[3];
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 4)
    return CAPS_BADTYPE;

  we[0] = bary[0];
  we[1] = bary[1];

  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  in[3] = discBody->elems[eIndex-1].gIndices[6] - 1;

  for (i = 0; i < rank; i++) {
/*  result[i] = (1.0-we[0])*((1.0-we[1])*data[rank*in[0]+i]  +
                                  we[1] *data[rank*in[3]+i]) +
                     we[0] *((1.0-we[1])*data[rank*in[1]+i]  +
                                  we[1] *data[rank*in[2]+i]);  */
    d_bar[rank*in[0]+i] += (1.0-we[0])*(1.0-we[1])*r_bar[i];
    d_bar[rank*in[1]+i] +=      we[0] *(1.0-we[1])*r_bar[i];
    d_bar[rank*in[2]+i] +=      we[0] *     we[1] *r_bar[i];
    d_bar[rank*in[3]+i] += (1.0-we[0])*     we[1] *r_bar[i];
  }

  return CAPS_SUCCESS;
}


/* Interpolation Bar for cell centered data */
static int
interpolateBar_CellCenter(capsDiscr *discr, int bIndex, int eIndex, int rank,
                          double *r_bar, double *d_bar)
{
  int    ind, i;
  capsBodyDiscr *discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].ndata != 1)
    return CAPS_BADTYPE;

  ind = discBody->elems[eIndex-1].dIndices[0] - 1;
  for (i = 0; i < rank; i++) {
/*  result[i] = data[rank*ind+i];  */
    d_bar[rank*ind+i] += r_bar[i];
  }

  return CAPS_SUCCESS;
}


/* Interpolate bar selector */
int aim_interpolateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                       double *bary, int rank, double *r_bar, double *d_bar)
{
  int status;
  capsBodyDiscr *discBody=NULL;
  capsEleType *eletype=NULL;

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" aimTransferUtil/aim_interpolateBar: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" aimTransferUtil/aim_interpolateBar: name = %s, eIndex = %d [1-%d]!\n",
           name, eIndex, discBody->nElems);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  eletype = &discr->types[discBody->elems[eIndex-1].tIndex-1];
  if (eletype->ndata == 0 ) {
    /* data at the reference locations */
    if (eletype->nref == 3) {
      /* Linear triangle */
      status = interpolateBar_LinearTriangle(discr, bIndex, eIndex, bary, rank,
                                             r_bar, d_bar);
    } else if (eletype->nref == 4) {
      /* Linear quad */
      status = interpolateBar_LinearQuad(discr, bIndex, eIndex, bary, rank, r_bar,
                                         d_bar);
    } else {
      printf(" aimTransferUtil/aim_interpolateBar: name = %s, eIndex = %d, nref not recognized!\n",
             name, eIndex);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  } else if (eletype->ndata == 1 ) {
    /* data at cell center */
    status = interpolateBar_CellCenter(discr, bIndex, eIndex, rank, r_bar, d_bar);
  } else {
    printf(" aimTransferUtil/aim_interpolateBar: name = %s, tIndex = %d, ndata = %d. Only supports ndata = 0 or ndata = 1!\n",
           name, discBody->elems[eIndex-1].tIndex, eletype->ndata);
    status = CAPS_BADTYPE;
    goto cleanup;
  }

cleanup:
  if (status != CAPS_SUCCESS)
    printf("Premature exit in transferUtils interpolation name = %s, status = %d\n",
           name, status);
  return status;
}


/* Integration for a linear triangular element */
static int
integration_LinearTriangle(capsDiscr *discr, int bIndex, int eIndex, int rank,
                           double *data, double *result)
{
  int        i, in[3], global[3], ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 3)
    return CAPS_BADTYPE;

  /* Element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  if (data == NULL) {
    *result = 3.0*area;
    return CAPS_SUCCESS;
  }

  for (i = 0; i < rank; i++) {
    result[i] = (data[rank*in[0]+i] +
                 data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area;
  }

cleanup:
  return status;
}


/* Integration for a linear quadrilateral element */
static int
integration_LinearQuad(capsDiscr *discr, int bIndex, int eIndex, int rank,
                       double *data, double *result)
{
  int        i, in[4], global[4], ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area, area2;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  if ((eIndex <= 0) || (eIndex > discBody->nElems)) return CAPS_BADINDEX;

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 4)
    return CAPS_BADTYPE;

  /* Element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  in[3] = discBody->elems[eIndex-1].gIndices[6] - 1;

  /* convert to global indices */
  for (i = 0; i < 4; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[3], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x1, x2, x3);

  area2  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  if (data == NULL) {
    *result = 3.0*area + 3.0*area2;
    return CAPS_SUCCESS;
  }

  for (i = 0; i < rank; i++) {
    result[i] = (data[rank*in[0]+i] +
                 data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area +
                (data[rank*in[0]+i] +
                 data[rank*in[2]+i] +
                 data[rank*in[3]+i])*area2;
  }

cleanup:
  return status;
}


/* Integration for an element with cell centered data */
static int
integration_CellCenter(capsDiscr *discr, int bIndex, int eIndex, int rank,
                       double *data, double *result)
{
  int        i, in[4], global[4], ind, ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area, area2;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  /* Triangle or first half of the quad */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/2.0;      /* 1/2 for area */

  ind = discBody->elems[eIndex-1].dIndices[0] - 1;

  if (data == NULL) {
    *result = area;
  } else {
    for (i = 0; i < rank; i++) {
      result[i] = data[rank*ind+i]*area;
    }
  }

  /* second half of the quad */
  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref == 4) {

    in[3] = discBody->elems[eIndex-1].gIndices[6]-1;

    /* convert to global indices */
    global[3] = discr->tessGlobal[2*in[3]+1];

    /* get coordinates */
    status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz2);
    AIM_STATUS(discr->aInfo, status);
    status = EG_getGlobal(discBody->tess, global[3], &ptype, &pindex, xyz3);
    AIM_STATUS(discr->aInfo, status);

    x1[0] = xyz2[0] - xyz1[0];
    x2[0] = xyz3[0] - xyz1[0];
    x1[1] = xyz2[1] - xyz1[1];
    x2[1] = xyz3[1] - xyz1[1];
    x1[2] = xyz2[2] - xyz1[2];
    x2[2] = xyz3[2] - xyz1[2];

    CROSS(x1, x2, x3);

    area2  = sqrt(DOT(x3, x3))/2.0;      /* 1/2 for area */

    if (data == NULL) {
      *result += area2;
    } else {
      for (i = 0; i < rank; i++) {
        result[i] += data[rank*ind+i]*area2;
      }
    }
  }

cleanup:
  return status;
}


/* Integration selector */
int aim_integration(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                    int rank, double *data, double *result)
{
  int status;
  capsBodyDiscr *discBody=NULL;
  capsEleType *eletype=NULL;

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" aimTransferUtil/aim_integration: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" aimTransferUtil/aim_integration: name = %s, eIndex = %d [1-%d]!\n",
           name, eIndex, discBody->nElems);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  eletype = &discr->types[discBody->elems[eIndex-1].tIndex-1];
  if (eletype->ndata == 0 ) {
    /* data at the reference locations */
    if (eletype->nref == 3) {
      /* Linear triangle */
      status = integration_LinearTriangle(discr, bIndex, eIndex, rank, data,
                                          result);
    } else if (eletype->nref == 4) {
      /* Linear quad */
      status = integration_LinearQuad(discr, bIndex, eIndex, rank, data, result);
    } else {
      printf(" aimTransferUtil/aim_integration: name = %s, eIndex = %d [1-%d], nref not recognized!\n",
             name, eIndex, discBody->nElems);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  } else if (eletype->ndata == 1 ) {
    /* data at cell center */
    status = integration_CellCenter(discr, bIndex, eIndex, rank, data, result);
  } else {
    printf(" aimTransferUtil/aim_integration: name = %s, tIndex = %d, ndata = %d. Only supports ndata = 0 or ndata = 1!\n",
           name, discBody->elems[eIndex-1].tIndex, eletype->ndata);
    status = CAPS_BADTYPE;
    goto cleanup;
  }

cleanup:
  if (status != CAPS_SUCCESS)
    printf("Premature exit in transferUtils interpolation name = %s, status = %d\n",
           name, status);
  return status;
}


/* Integration bar for a linear triangular element */
static int
integrateBar_LinearTriangle(capsDiscr *discr, int bIndex, int eIndex, int rank,
                            double *r_bar, double *d_bar)
{
  int        i, in[3], global[3], ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 3)
    return CAPS_BADTYPE;

  /* Element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  for (i = 0; i < rank; i++) {
/*  result[i] = (data[rank*in[0]+i] +
                 data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area;  */
    d_bar[rank*in[0]+i] += area*r_bar[i];
    d_bar[rank*in[1]+i] += area*r_bar[i];
    d_bar[rank*in[2]+i] += area*r_bar[i];
  }

cleanup:
  return CAPS_SUCCESS;
}


/* Integration bar for a linear quadrilateral element */
static int
integrateBar_LinearQuad(capsDiscr *discr, int bIndex, int eIndex, int rank,
                        double *r_bar, double *d_bar)
{
  int        i, in[4], global[4], ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area, area2;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref != 4)
    return CAPS_BADTYPE;

  /* Element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  in[3] = discBody->elems[eIndex-1].gIndices[6] - 1;

  /* convert to global indices */
  for (i = 0; i < 4; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[3], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area2  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  for (i = 0; i < rank; i++) {

/*  result[i] = (data[rank*in[0]+i] + data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area +
                (data[rank*in[0]+i] + data[rank*in[2]+i] +
                 data[rank*in[3]+i])*area2;  */
    d_bar[rank*in[0]+i] += (area + area2)*r_bar[i];
    d_bar[rank*in[1]+i] +=  area         *r_bar[i];
    d_bar[rank*in[2]+i] += (area + area2)*r_bar[i];
    d_bar[rank*in[3]+i] +=         area2 *r_bar[i];
  }

cleanup:
  return CAPS_SUCCESS;
}


/* Integration bar for an element with cell centered data */
static int
integrateBar_CellCenter(capsDiscr *discr, int bIndex, int eIndex, int rank,
                        double *r_bar, double *d_bar)
{
  int        i, in[4], global[4], ind, ptype, pindex, status;
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area, area2;
  capsBodyDiscr *discBody=NULL;

  discBody = &discr->bodys[bIndex-1];

  /* Element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++) {
    global[i] = discr->tessGlobal[2*in[i]+1];
  }

  /* get coordinates */
  status = EG_getGlobal(discBody->tess, global[0], &ptype, &pindex, xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[1], &ptype, &pindex, xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];

  CROSS(x3, x1, x2);

  area  = sqrt(DOT(x3, x3))/2.0;      /* 1/2 for area */

  ind = discBody->elems[eIndex-1].dIndices[0] - 1;

  for (i = 0; i < rank; i++) {
/*  result[i] = data[rank*ind+i]*area;  */
    d_bar[rank*ind+i] += area*r_bar[i];
  }

  if (discr->types[discBody->elems[eIndex-1].tIndex-1].nref == 4) {

    in[3] = discBody->elems[eIndex-1].gIndices[6]-1;

    /* convert to global indices */
    global[3] = discr->tessGlobal[2*in[3]+1];

    /* get coordinates */
    status = EG_getGlobal(discBody->tess, global[2], &ptype, &pindex, xyz2);
    AIM_STATUS(discr->aInfo, status);
    status = EG_getGlobal(discBody->tess, global[3], &ptype, &pindex, xyz3);
    AIM_STATUS(discr->aInfo, status);

    x1[0] = xyz2[0] - xyz1[0];
    x2[0] = xyz3[0] - xyz1[0];
    x1[1] = xyz2[1] - xyz1[1];
    x2[1] = xyz3[1] - xyz1[1];
    x1[2] = xyz2[2] - xyz1[2];
    x2[2] = xyz3[2] - xyz1[2];

    CROSS(x3, x1, x2);

    area2  = sqrt(DOT(x3, x3))/2.0;      /* 1/2 for area */

    for (i = 0; i < rank; i++) {
/*    result[i] += data[rank*ind+i]*area2;  */
      d_bar[rank*ind+i] += area2*r_bar[i];
    }
  }

cleanup:
  return status;
}


/* Integrate bar selector */
int aim_integrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                     int rank, double *r_bar, double *d_bar)
{
  int status;
  capsBodyDiscr *discBody=NULL;
  capsEleType *eletype=NULL;

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" aimTransferUtil/aim_integrateBar: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" aimTransferUtil/aim_integrateBar: name = %s, eIndex = %d [1-%d]!\n",
           name, eIndex, discBody->nElems);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  eletype = &discr->types[discBody->elems[eIndex-1].tIndex-1];
  if (eletype->ndata == 0 ) {
    /* data at the reference locations */
    if (eletype->nref == 3) {
      /* Linear triangle */
      status = integrateBar_LinearTriangle(discr, bIndex, eIndex, rank, r_bar,
                                           d_bar);
    } else if (eletype->nref == 4) {
      /* Linear quad */
      status = integrateBar_LinearQuad(discr, bIndex, eIndex, rank, r_bar,
                                       d_bar);
    } else {
      printf(" aimTransferUtil/aim_integrateBar: name = %s, eIndex = %d [1-%d], nref not recognized!\n",
             name, eIndex, discBody->nElems);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  } else if (eletype->ndata == 1 ) {
    /* data at cell center */
    status = integrateBar_CellCenter(discr, bIndex, eIndex, rank, r_bar, d_bar);
  } else {
    printf(" aimTransferUtil/aim_integrateBar: name = %s, tIndex = %d, ndata = %d. Only supports ndata = 0 or ndata = 1!\n",
           name, discBody->elems[eIndex-1].tIndex, eletype->ndata);
    status = CAPS_BADTYPE;
    goto cleanup;
  }

cleanup:
  if (status != CAPS_SUCCESS)
    printf("Premature exit in transferUtils integrateBar name = %s, status = %d\n",
           name, status);
  return status;
}
