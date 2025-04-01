/*
 *************************************************************************
 *                                                                       *
 * udpFFDdeform -- udp file to deform b-spline geometry using an FFD box *
 *                                                                       *
 *            Written by Marlena Gomez @ MIT                             *
 *                   and Marshall Galbraith @ MIT                        *
 *            Patterned after .c code written by John Dannenhoffer       *
 *                                                @ Syracuse University  *
 *                                                                       *
 *************************************************************************
 */

/*
 * Copyright (C) 2011/2025 Marlena Gomez and Marshall Galbraith @ MIT
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */
#include <vector>
#include "egads_dot.h"
#include "egadsSplineFit.h" /* EG_spline.1dfit */
#include "spline3d_eval.hpp"

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#define PATH_MAX   _MAX_PATH
#else
#include <unistd.h>
#include <limits.h>
#endif

extern "C" {
 #define NUMUDPINPUTBODYS 2
 #define NUMUDPARGS 3
 #include "udpUtilities.h"
} /* extern "C" */


// n row needs to be number of cp in the box
#define DESIGN_PARAM_ARG(UDP      )                    (UDP)->arg[0]
#define DESIGN_PARAM_VAL(UDP, I, J) ((const double *) ((UDP)->arg[0].val))[(I)*(UDP)->arg[0].ncol+(J)]
#define DESIGN_PARAM_DOT(UDP, I, J) ((const double *) ((UDP)->arg[0].dot))[(I)*(UDP)->arg[0].ncol+(J)]
#define WRITEXYZ_VAL(    UDP      ) ((const int *)    ((UDP)->arg[1].val))[0]
#define PENALTY_ARG(     UDP      )                    (UDP)->arg[2]
#define PENALTY_VAL(     UDP      ) ((double *)       ((UDP)->arg[2].val))[0]
#define PENALTY_DOT(     UDP      ) ((double *)       ((UDP)->arg[2].dot))[0]

// maybe change to 3 so loop over fewer things
#define DESIGN_PARAM_SURREAL(UDP, I, J) SurrealS<1>( DESIGN_PARAM_VAL(UDP, I, J), DESIGN_PARAM_DOT(UDP, I, J))
#define PENALTY_SURREAL(     UDP      ) SurrealS<1>( PENALTY_VAL(     UDP      ), PENALTY_DOT(     UDP      ))

/* data about possible arguments */
static const char  *argNames[NUMUDPARGS] = {"design_params", "writexyz", "penalty"   };
static int          argTypes[NUMUDPARGS] = {ATTRREALSEN,     ATTRINT,    -ATTRREALSEN};
static int          argIdefs[NUMUDPARGS] = {0 ,              0,          0,          };
static double       argDdefs[NUMUDPARGS] = {0.,              0.,         0.          };

extern "C" void EG_getGeometryLen(const egObject *geom, int *nivec, int *nrvec);

extern "C"
int buildBody(ego eOrigBody,                       /* (in) original body to rebuild */
              const int nNewNode, ego *eNewNodes,  /* (in) mixture of orig and new Nodes    (Node BodyIndex ordered) */
              const int nNewCurv, ego *eNewCurvs,  /* (in) mixture of orig and new Curves   (Edge BodyIndex ordered) */
              const int nNewSurf, ego *eNewSurfs,  /* (in) mixture of orig and new Surfaces (Face BodyIndex ordered) */
              ego *eNewBody,                       /* (out) the new body with all geometry replaced */
              int lenmsg,                          /* (in) error message length */
              char *message);                      /* (out) error message */

template<class T> T    DESIGN_PARAM(udp_T *udp, int i, int j);
template<> double      DESIGN_PARAM< double      >(udp_T *udp, int i, int j ) { return DESIGN_PARAM_VAL(udp,i,j); }
template<> SurrealS<1> DESIGN_PARAM< SurrealS<1> >(udp_T *udp, int i, int j ) { return DESIGN_PARAM_SURREAL(udp,i,j); }

template<class T> T    PENALTY(udp_T *udp);
template<> double      PENALTY< double      >( udp_T *udp ) { return PENALTY_VAL(udp); }
template<> SurrealS<1> PENALTY< SurrealS<1> >( udp_T *udp ) { return PENALTY_SURREAL(udp); }

template class SurrealS<1>;
typedef SurrealS<1> SurrealS1;

#define LENMSG 512
#define DEFORMBODY "ffdDeformBody"

template<class T>
struct BSpline {
  int eIndex;
  int ivec[9];
  T *rvec;

  BSpline() : eIndex(0), rvec(NULL) {}
  ~BSpline() { FREE(rvec); }
};

struct ffdBoxData
{
  ego ebody;
  int i_degree;
  int j_degree;
  int k_degree;
  int i_cp;
  int j_cp;
  int k_cp;
  double i_range[2];
  double j_range[2];
  double k_range[2];
  std::vector<double> u, v, w;
  int *ivec;
  double *rvec; // original point distribution, without deformation from design params

  ffdBoxData() :
    ebody(NULL),
    i_degree(0),
    j_degree(0),
    k_degree(0),
    i_cp(0),
    j_cp(0),
    k_cp(0),
    ivec(NULL), rvec(NULL)
  {
    i_range[0] = i_range[1] = 0;
    j_range[0] = j_range[1] = 0;
    k_range[0] = k_range[1] = 0;
  }

  ~ffdBoxData() { delete [] ivec; delete [] rvec; }
};

struct uvArrayData
{
  // for a NODE ivec is {0} and rvec is {x,y,z}
  ego    topo;
  int    eIndex; // 0-based index into enodes, eedges, or efaces of original deform body
  double *uvw;   // u,v,w value in the FFD box, -1 for points that are fixed
  int    *ivec;  // B-spline header
  double *rvec;  // original control point distribution, without perturbation from design params

  int    nCP;    // number of control points
  double *CP;    // reference to control points in rvec
  double *uk;    // reference to u-knots in rvec
  double *vk;    // reference to v-knots in rvec

  uvArrayData() :
    topo(NULL), eIndex(0),
    uvw(NULL), ivec(NULL), rvec(NULL),
    nCP(0), CP(NULL), uk(NULL), vk(NULL) {}

  uvArrayData(const uvArrayData&) = delete;

  uvArrayData(uvArrayData&& copy) :
    topo(copy.topo), eIndex(copy.eIndex),
    uvw(copy.uvw), ivec(copy.ivec), rvec(copy.rvec),
    nCP(copy.nCP), CP(copy.CP),
    uk(copy.uk), vk(copy.vk)
  {
    copy.topo   = NULL;
    copy.eIndex = 0;
    copy.uvw    = NULL;
    copy.ivec   = NULL;
    copy.rvec   = NULL;
    copy.nCP    = 0;
    copy.CP     = NULL;
    copy.uk     = NULL;
    copy.vk     = NULL;
  }

  ~uvArrayData()
  {
    FREE(uvw);
    FREE(ivec);
    FREE(rvec);
    CP = NULL;
    uk = NULL;
    vk = NULL;
  }
};

struct udpDotCache_T
{
  ffdBoxData ffdBox;
  ego source_body;
  ego deformed_body;
  std::vector<uvArrayData> uvNodes;
  std::vector<uvArrayData> uvEdges;
  std::vector<uvArrayData> uvFaces;

  udpDotCache_T() : source_body(NULL), deformed_body(NULL) {}
  ~udpDotCache_T()
  {
    EG_deleteObject(source_body);
    EG_deleteObject(deformed_body);
  }
};

static int freePrivateData(void *data)
{
  udpDotCache_T *cache = (udpDotCache_T*)data;

  delete cache;

  return EGADS_SUCCESS;
}

#define FREEUDPDATA freePrivateData
extern "C" {
/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
} /* extern "C" */


/* minimum p (smoothmin tetrahedron neighbors) / (sum tetrahedron neighbors) * num tet neighbors */
double min_p = 1.0e8;

/*
Smooth approximation of minimum function
Larger values of alpha approach a sharper minimum function
*/
template < class T >
static void smoothmin(const T x, const T y, const double alpha, T& smoothmin)
{
  // The if statements below is a more precise version of
  // T m = min(x,y);
  // smoothmin = m - log( exp(-alpha*(x-m)) + exp(-alpha*(y-m)) )/alpha;

  if ( x < y )
    smoothmin = x - log1p( exp(-alpha*(y-x)) )/alpha;
  else
    smoothmin = y - log1p( exp(-alpha*(x-y)) )/alpha;
}

/*
Return Tetrahedron volume determinant given four points that make up the tetrahedron
Order of points may change sign of the determinant
*/
template < class T >
T tetDeter(std::vector<T> &pa, std::vector<T> &pb, std::vector<T> &pc, std::vector<T> &pd)
{
  T adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  T bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  T det;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;

  det = adz * (bdxcdy - cdxbdy)
      + bdz * (cdxady - adxcdy)
      + cdz * (adxbdy - bdxady);

  return det;

}

/*
Find 6 neighboring control points (above, below, north, south, east, west of current control point in FFD grid)
Neighbors are the next CP over by one in x,y,z directions (total of 6 neighbors for interior)
*/
template < class T >
void findNeighbors(
  const int ivec[],
  const T rvec[],
  const int i,
  const int j,
  const int k,
  std::vector<T>& center,
  std::vector<T>& left,
  std::vector<T>& right,
  std::vector<T>& above,
  std::vector<T>& below,
  std::vector<T>& forward,
  std::vector<T>& backward)
{
  int icp, iknot, jcp, jknot, kcp, kknot, nknot;
  int center_index, above_index, below_index, right_index, left_index, forward_index, backward_index;

  // ideg = ivec[0];
  icp = ivec[1];
  iknot = ivec[2];
  // jdeg = ivec[3];
  jcp = ivec[4];
  jknot = ivec[5];
  // kdeg = ivec[6];
  kcp = ivec[7];
  kknot = ivec[8];

  nknot = iknot+jknot+kknot;

  center_index = nknot + 3 * (i + j * icp + k * icp * jcp);
  right_index  = nknot + 3 * ((i+1) + j * icp + k * icp * jcp);
  left_index   = nknot + 3 * ((i-1) + j * icp + k * icp * jcp);
  above_index  = nknot + 3 * (i + (j+1) * icp + k * icp * jcp);
  below_index  = nknot + 3 * (i + (j-1) * icp + k * icp * jcp);
  forward_index  = nknot + 3 * (i + j * icp + (k+1) * icp * jcp);
  backward_index  = nknot + 3 * (i + j * icp + (k-1) * icp * jcp);

  center[0]    = rvec[center_index  ]; // center
  center[1]    = rvec[center_index+1];
  center[2]    = rvec[center_index+2];

  if (i > 0) {
    left[0]   = rvec[left_index  ]; // negative u
    left[1]   = rvec[left_index+1];
    left[2]   = rvec[left_index+2];
  }

  if (i < icp-1) {
    right[0]  = rvec[right_index  ]; // positive u
    right[1]  = rvec[right_index+1];
    right[2]  = rvec[right_index+2];
  }

  if (j > 0) {
    below[0]   = rvec[below_index  ]; // negative v
    below[1]   = rvec[below_index+1];
    below[2]   = rvec[below_index+2];
  }

  if (j < jcp-1) {
    above[0]  = rvec[above_index  ]; // positive v
    above[1]  = rvec[above_index+1];
    above[2]  = rvec[above_index+2];
  }

  if (k > 0) {
    backward[0]  = rvec[backward_index]; // negative w
    backward[1]  = rvec[backward_index+1];
    backward[2]  = rvec[backward_index+2];
  }

  if (k < kcp-1) {
    forward[0]  = rvec[forward_index  ]; // positive w
    forward[1]  = rvec[forward_index+1];
    forward[2]  = rvec[forward_index+2];
  }
}

/*
Compute the 3D penalty value for the tetrahedron around point (i,j,k)
*/
template < class T >
void computePenalty(
  const int ivec[],
  const T rvec[],
  int i,
  int j,
  int k,
  T& penalty_out)
{
  int            icp, jcp, kcp, index;
  std::vector<T> center(3), above(3), below(3), left(3), right(3), forward(3), backward(3);
  std::vector<T> vol(8), min_vol(4), min_vol2(2);
  T              sum_vol = 0., final_min = 0., penalty = 0., p = 0.;
  bool           isCorner = false, isEdge = false, isFace = false;

  // the constants for the penalty cost function
  double         alpha = 1000., c1 =1.0, c2 = 10.0, c3 = 0.5, c4 = 0.0;

  icp   = ivec[1];
  jcp   = ivec[4];
  kcp   = ivec[7];

  /*
  Definitions of orientation:
  - x-axis going to right
  - y-axis going upward
  - z-axis pointing out of the page

  - left and right neighbors represent negative and positive movement (respectively) in x-direction
  - above and below neighbors represent negative and positive movement (respectively) in y-direction
  - forward and backward neighbors represent negative and positive movement (respectively) in z-direction

  Corners and edges:
  - bottom_back_left means - bottom (y=0), back (z=0), and left (x=0)
  - top_front_right means - top (y=jcp-1), front (z=kcp-1), and right (x=icp-1)
  */


  for (index = 0; index < 3; index++) {
    center[index] = 0.;
    above[index] = 0.;
    below[index] = 0.;
    left[index] = 0.;
    right[index] = 0.;
    forward[index] = 0.;
    backward[index] = 0.;
  }

  findNeighbors(ivec, rvec, i, j, k, center, left, right, above, below,  forward, backward);

  /* i == 0 face of the box */
  if (i == 0) {
    if (j == 0) {
      if (k == 0) {
        isCorner = true;
        vol[0] = tetDeter(center, forward, above, right);
      } else if (k == (kcp -1)) {
        isCorner = true;
        vol[0] = tetDeter(backward, center, above, right);
      } else {
        isEdge = true;
        vol[0] = tetDeter(backward, center, above, right);
        vol[1] = tetDeter(center, forward, above, right);
      }
    } else if (j == jcp-1) {
      if (k == 0) {
        isCorner = true;
        vol[0] = tetDeter(center, forward, right, below);
      } else if (k == (kcp-1)) {
        isCorner = true;
        vol[0] = tetDeter(backward, center, right, below);
      } else {
        isEdge = true;
        vol[0] = tetDeter(backward, center, right, below);
        vol[1] = tetDeter(center, forward, right, below);
      }
    } else if (k == 0) {
      isEdge = true;
      vol[0] = tetDeter(center, forward, above, right);
      vol[1] = tetDeter(center, forward, right, below);
    } else if (k == kcp-1) {
      isEdge = true;
      vol[0] = tetDeter(backward, center, above, right);
      vol[1] = tetDeter(backward, center, right, below);
    } else {
      isFace = true;
      vol[0] = tetDeter(center, forward, above, right);
      vol[1] = tetDeter(center, forward, right, below);
      vol[2] = tetDeter(backward, center, above, right);
      vol[3] = tetDeter(backward, center, right, below);
    }
  }

  /* i == icp-1 face of the box */
  else if (i == (icp-1)) {
    if (j == 0) {
      if (k == 0) {
        isCorner = true;
        vol[0] = tetDeter(center, forward, left, above);
      } else if (k == (kcp -1)) {
        isCorner = true;
        vol[0] = tetDeter(backward, center, left, above);
      } else {
        isEdge = true;
        vol[0] = tetDeter(backward, center, left, above);
        vol[1] = tetDeter(center, forward, left, above);
      }
    } else if (j == jcp-1) {
      if (k == 0) {
        isCorner = true;
        vol[0] = tetDeter(center, forward, below, left);
      } else if (k == (kcp-1)) {
        isCorner = true;
        vol[0] = tetDeter(backward, center, below, left);
      } else {
        isEdge = true;
        vol[0] = tetDeter(backward, center, below, left);
        vol[1] = tetDeter(center, forward, below, left);
      }
    } else if (k == 0) {
      isEdge = true;
      vol[0] = tetDeter(center, forward, left, above);
      vol[1] = tetDeter(center, forward, below, left);
    } else if (k == kcp-1) {
      isEdge = true;
      vol[0] = tetDeter(backward, center, left, above);
      vol[1] = tetDeter(backward, center, below, left);
    } else {
      isFace = true;
      vol[0] = tetDeter(center, forward, left, above);
      vol[1] = tetDeter(center, forward, below, left);
      vol[2] = tetDeter(backward, center, left, above);
      vol[3] = tetDeter(backward, center, below, left);
    }
  }

  else if (j == 0) {
    if (k == 0) {
      isEdge = true;
      vol[0] = tetDeter(center, forward, left, above);
      vol[1] = tetDeter(center, forward, above, right);
    } else if (k == (kcp-1)) {
      isEdge = true;
      vol[0] = tetDeter(backward, center, left, above);
      vol[1] = tetDeter(backward, center, above, right);
    } else {
      isFace = true;
      vol[0] = tetDeter(center, forward, left, above);
      vol[1] = tetDeter(center, forward, above, right);
      vol[2] = tetDeter(backward, center, left, above);
      vol[3] = tetDeter(backward, center, above, right);
    }
  }

  else if (j == (jcp-1)) {
    if (k == 0) {
      isEdge = true;
      vol[0] = tetDeter(center, forward, below, left);
      vol[1] = tetDeter(center, forward, right, below);
    } else if (k == (kcp-1)) {
      isEdge = true;
      vol[0] = tetDeter(backward, center, below, left);
      vol[1] = tetDeter(backward, center, right, below);
    } else {
      isFace = true;
      vol[0] = tetDeter(center, forward, below, left);
      vol[1] = tetDeter(center, forward, right, below);
      vol[2] = tetDeter(backward, center, below, left);
      vol[3] = tetDeter(backward, center, right, below);
    }
  } else if (k == 0) {
    isFace = true;
    vol[0] = tetDeter(center, forward, below, left);
    vol[1] = tetDeter(center, forward, above, right);
    vol[2] = tetDeter(center, forward, right, below);
    vol[3] = tetDeter(center, forward, left, above);
  } else if (k == kcp-1) {
    isFace = true;
    vol[0] = tetDeter(backward, center, below, left);
    vol[1] = tetDeter(backward, center, above, right);
    vol[2] = tetDeter(backward, center, right, below);
    vol[3] = tetDeter(backward, center, left, above);
  }

  /* Interior control points */
  else {
    vol[0] = tetDeter(center, forward, above, right);
    vol[1] = tetDeter(center, forward, right, below);
    vol[2] = tetDeter(center, forward, below, left);
    vol[3] = tetDeter(center, forward, left, above);


    // lower four tets
    vol[4] = tetDeter(backward, center, above, right);
    vol[5] = tetDeter(backward, center, right, below);
    vol[6] = tetDeter(backward, center, below, left);
    vol[7] = tetDeter(backward, center, left, above);

    for (index = 0; index < 8; index++) sum_vol += vol[index];
    for (index = 0; index < 8; index++) vol[index] /= sum_vol;


    smoothmin(vol[0], vol[1], alpha, min_vol[0]);
    smoothmin(vol[2], vol[3], alpha, min_vol[1]);
    smoothmin(vol[4], vol[5], alpha, min_vol[2]);
    smoothmin(vol[6], vol[7], alpha, min_vol[3]);

    smoothmin(min_vol[0], min_vol[1], alpha, min_vol2[0]);
    smoothmin(min_vol[2], min_vol[3], alpha, min_vol2[1]);

    smoothmin(min_vol2[0], min_vol2[1], alpha, final_min);

    p = 8. * final_min; // scale p to 1

  }

  if (isCorner) {
    for (index = 0; index < 2; index++) sum_vol += vol[0];
    vol[0] /= sum_vol;
    p = vol[0];

  } else if (isEdge) {
    for (index = 0; index < 2; index++) sum_vol += vol[index];
    for (index = 0; index < 2; index++) vol[index] /= sum_vol;

    smoothmin(vol[0], vol[1], alpha, final_min);
    p = 2. * final_min; // scale p to 1
  } else if (isFace) {
    for (index = 0; index < 4; index++) sum_vol += vol[index];
    for (index = 0; index < 4; index++) vol[index] /= sum_vol;

    smoothmin(vol[0], vol[1], alpha, min_vol2[0]);
    smoothmin(vol[2], vol[3], alpha, min_vol2[1]);
    smoothmin(min_vol2[0], min_vol2[1], alpha, final_min);
    p = 4. * final_min; // scale p to 1
  }

  penalty = c1 * erfc( (p+c4) * c2) * -c3 * (p-c3);
  penalty_out = penalty;

}

/*
  Raise error if source body is not correctly attributed
  i.e. The body has not gone through step 1 (FFDSplines)
*/
int checkDeformBody(ego context, ego ebody, char *message)
{
  int    status = EGADS_SUCCESS;
  ego    *efaces = NULL, *eedges = NULL;
  int    atype, nface, nedge;
  int    alen, i;
  const int     *pint;
  const double  *pdouble;

  ROUTINE(checkDeformBody);

  /* --------------------------------------------------------------- */

  // check at least one face/edge has an attribute indicating it is perturbed by the ffd box

  status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  for (i = 0; i < nface; i++) {
    status = EG_attributeRet(efaces[i], "ffdID", &atype, &alen, &pint, &pdouble, NULL);
    if (status == EGADS_SUCCESS)
      goto cleanup;
  }

  status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  for (i = 0; i < nedge; i++) {
    status = EG_attributeRet(eedges[i], "ffdID", &atype, &alen, &pint, &pdouble, NULL);
    if (status == EGADS_SUCCESS)
      goto cleanup;
  }

  snprintf(message, LENMSG, "The source body is missing required attributes.\n\
    Please make sure the first step (udpFFDsplines) has been run first.\n");
  status  = EGADS_RANGERR;

cleanup:

  FREE(efaces);
  FREE(eedges);

  return status;
}

/*
Returns EGADS_NOTFOUND (status = -1) if the ebody is not an FFD body.
Returns EGADS_RANGERR (status = -16) if the ebody has attributes for an FFD body, but
the attributes are incorrect or invalid in some way.
Returns EGADS_SUCCESS if the ebody is an FFD body with correct attributes.
Fills in the ffdBoxData structure if successful.
*/
int
checkFFDbox(ego context, ego ebody, ffdBoxData &ffd_box, char *message)
{
  int    status = EGADS_SUCCESS;
  ego    eref, *efaces = NULL, *eloops = NULL, *eedges = NULL;
  double data[18], bbox[6];
  int    oclass, mtype, atype, nface, nloop, nedge, *senses;
  int    alen, i, j, k, iface, nfaceNeighbor;
  bool   uMin_attr = false, vMin_attr = false, wMin_attr = false;
  ego    *efaceNeighbors = NULL, *faceCandidates = NULL;
  const double  *pdouble;

  ROUTINE(checkFFDbox);

  /* --------------------------------------------------------------- */

  status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  /* FFD box must have 6 faces */
  if (nface != 6) {
    status = EGADS_NOTFOUND;
    goto cleanup;
  }

  for (iface = 0; iface < nface; iface++) {

    status = EG_getTopology(efaces[iface], &eref, &oclass, &mtype,
                            data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                            data, &nedge, &eedges, &senses);
    CHECK_STATUS(EG_getTopology);

    faceCandidates = new ego[nface];

    /* fill in the faceCandidates array used in finding opposite pairs of faces */
    for (j = 0; j < nface; j++) {
      faceCandidates[j] = NULL;
      if (iface != j) faceCandidates[j] = efaces[j];
    }

    /* loop over the edges and get the face neighbors of the edge */
    for (j = 0; j < nedge; j++) {
      status = EG_getBodyTopos(ebody, eedges[j], FACE, &nfaceNeighbor, &efaceNeighbors);
      CHECK_STATUS(EG_getBodyTopos);

      for (int k = 0; k < nface; k++) {
        for (int m = 0; m < nfaceNeighbor; m++) {
          if (efaceNeighbors[m] == faceCandidates[k]) faceCandidates[k] = NULL;
        }
      }
      FREE(efaceNeighbors);
    }

    status = EG_attributeRet(efaces[iface], "ffd_uMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      uMin_attr = true;

      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_uMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      ffd_box.i_degree = (int)(size_t) pdouble[0];

      if (ffd_box.i_degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (ffd_box.i_degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Adding specified number of control points in i direction. */
      if (alen == 2) {
        ffd_box.i_cp = (int)(size_t) pdouble[1];
        ffd_box.u.resize(ffd_box.i_cp);
        for (i = 0; i < ffd_box.i_cp; i++) {
          ffd_box.u[i] = double(i) / (ffd_box.i_cp-1);
        }
      } else {
        ffd_box.i_cp = alen-1;
        ffd_box.u.resize(ffd_box.i_cp);

        if (pdouble[1] != 0) {
          snprintf(message, LENMSG, "u-seuqence first value must be zero.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        if (pdouble[alen-1] != 1) {
          snprintf(message, LENMSG, "u-seuqence last value must be one.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        for (int i = 2; i < alen; i++) {
          if (pdouble[i-1] >= pdouble[i]) {
            snprintf(message, LENMSG, "u-seuqence must be in ascending order.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
        }

        for (i = 0; i < ffd_box.i_cp; i++) {
          ffd_box.u[i] = pdouble[i+1];
        }
      }

      /* check degree and number of control points relationship */
      if (ffd_box.i_degree+1 > ffd_box.i_cp ) {
        snprintf(message, LENMSG, "Number of control points is too few for degree of spline in u-direction.\n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Check that the opposite faces of the uMin attribute are not also attributed. */
      for (int k = 0; k < nface; k++) {
        if (faceCandidates[k] != NULL) {
          status = EG_attributeRet(faceCandidates[k], "ffd_vMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
          status = EG_attributeRet(faceCandidates[k], "ffd_wMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }

        }
      }
    }

    status = EG_attributeRet(efaces[iface], "ffd_vMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      vMin_attr = true;

      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_vMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      ffd_box.j_degree = (int)(size_t) pdouble[0];

      if (ffd_box.j_degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (ffd_box.j_degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Adding specified number of control points in j direction. */
      if (alen == 2) {
        ffd_box.j_cp = (int)(size_t) pdouble[1];
        ffd_box.v.resize(ffd_box.j_cp);
        for (j = 0; j < ffd_box.j_cp; j++) {
          ffd_box.v[j] = double(j) / (ffd_box.j_cp-1);
        }
      } else {
        ffd_box.j_cp = alen-1;
        ffd_box.v.resize(ffd_box.j_cp);

        if (pdouble[1] != 0) {
          snprintf(message, LENMSG, "v-seuqence first value must be zero.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        if (pdouble[alen-1] != 1) {
          snprintf(message, LENMSG, "v-seuqence last value must be one.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        for (int i = 2; i < alen; i++) {
          if (pdouble[i-1] >= pdouble[i]) {
            snprintf(message, LENMSG, "v-seuqence must be in ascending order.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
        }

        for (j = 0; j < ffd_box.j_cp; j++) {
          ffd_box.v[j] = pdouble[j+1];
        }
      }

      /* check degree and number of control points relationship */
      if (ffd_box.j_degree+1 > ffd_box.j_cp ) {
        snprintf(message, LENMSG, "Number of control points is too few for degree of spline in vMin direction.\n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }


      /* Check that the opposite faces of the uMin attribute are not also attributed. */
      for (int k = 0; k < nface; k++) {
        if (faceCandidates[k] != NULL) {
          status = EG_attributeRet(faceCandidates[k], "ffd_uMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
          status = EG_attributeRet(faceCandidates[k], "ffd_wMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
        }
      }
    }

    status = EG_attributeRet(efaces[iface], "ffd_wMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      wMin_attr = true;
      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_wMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      ffd_box.k_degree = (int)(size_t) pdouble[0];

      if (ffd_box.k_degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (ffd_box.k_degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Adding specified number of control points in k direction. */
      if (alen == 2) {
        ffd_box.k_cp = (int)(size_t) pdouble[1];
        ffd_box.w.resize(ffd_box.k_cp);
        for (k = 0; k < ffd_box.k_cp; k++) {
          ffd_box.w[k] = double(k) / (ffd_box.k_cp-1);
        }
      } else {
        ffd_box.k_cp = alen-1;
        ffd_box.w.resize(ffd_box.k_cp);

        if (pdouble[1] != 0) {
          snprintf(message, LENMSG, "w-seuqence first value must be zero.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        if (pdouble[alen-1] != 1) {
          snprintf(message, LENMSG, "w-seuqence last value must be one.\n");
          status  = EGADS_RANGERR;
          goto cleanup;
        }
        for (int i = 2; i < alen; i++) {
          if (pdouble[i-1] >= pdouble[i]) {
            snprintf(message, LENMSG, "w-seuqence must be in ascending order.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
        }

        for (k = 0; k < ffd_box.k_cp; k++) {
          ffd_box.w[k] = pdouble[k+1];
        }
      }

      /* check degree and number of control points relationship */
      if (ffd_box.k_degree+1 > ffd_box.k_cp) {
        snprintf(message, LENMSG, "Number of control points is too few for degree of spline in wMin direction.\n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Check that the opposite faces of the uMin attribute are not also attributed. */
      for (int k = 0; k < nface; k++) {
        if (faceCandidates[k] != NULL) {
          status = EG_attributeRet(faceCandidates[k], "ffd_uMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
          status = EG_attributeRet(faceCandidates[k], "ffd_vMin", &atype, &alen, NULL, &pdouble, NULL);
          if (status == 0) {
            snprintf(message, LENMSG, "ffd_uMin, ffd_vMin, ffd_wMin attributes cannot be on opposite sides of the FFD box.\n");
            status  = EGADS_RANGERR;
            goto cleanup;
          }
        }
      }
    }

    delete [] faceCandidates;
  }

  if (uMin_attr && vMin_attr && wMin_attr) {
    status = EG_getBoundingBox(ebody, bbox);
    CHECK_STATUS(EG_getBoundingBox);

    ffd_box.i_range[0] = bbox[0];
    ffd_box.j_range[0] = bbox[1];
    ffd_box.k_range[0] = bbox[2];
    ffd_box.i_range[1] = bbox[3];
    ffd_box.j_range[1] = bbox[4];
    ffd_box.k_range[1] = bbox[5];

    // printf("x range %f %f \n", bbox[0], bbox[3]);
    // printf("y range %f %f \n", bbox[1], bbox[4]);
    // printf("z range %f %f \n", bbox[2], bbox[5]);
    ffd_box.ebody = ebody;

    status = EGADS_SUCCESS;

  } else {
    status  = EGADS_NOTFOUND;
    goto cleanup;
  }


cleanup:

  FREE(efaces);

  return status;
}

/*
Fills in ivec and rvec for ffd_box struct.
Assumes all other values of the ffd_box were already filled in by checkFDDbox.
checkFFDBox has 0 for num_added_i_knots if number of cps was specified instead of not sequence
*/
int
makeFFDBox(ffdBoxData &ffd_box)
{
  int status = EGADS_SUCCESS;

  ROUTINE(makeFFDBox);

  /* --------------------------------------------------------------- */

  ffd_box.ivec = new int[9];

  // TODO: this logic requires uMin to be x direction, vMin to be y, and wMin to be z
  ffd_box.ivec[0] = ffd_box.i_degree;
  ffd_box.ivec[1] = ffd_box.i_cp;
  ffd_box.ivec[2] = ffd_box.ivec[0]+1+ffd_box.ivec[1]; // number of knots in i direction

  ffd_box.ivec[3] = ffd_box.j_degree;
  ffd_box.ivec[4] = ffd_box.j_cp;
  ffd_box.ivec[5] = ffd_box.ivec[3]+1+ffd_box.ivec[4]; // number of knots in j direction

  ffd_box.ivec[6] = ffd_box.k_degree;
  ffd_box.ivec[7] = ffd_box.k_cp;
  ffd_box.ivec[8] = ffd_box.ivec[6]+1+ffd_box.ivec[7]; // number of knots in k direction

  int    ncp   = ffd_box.i_cp*ffd_box.j_cp*ffd_box.k_cp;
  int    uknot = ffd_box.ivec[2];
  int    vknot = ffd_box.ivec[5];
  int    wknot = ffd_box.ivec[8];
  int    nknot = uknot+vknot+wknot;

  ffd_box.rvec = new double[nknot+3*ncp];

  int    i, j, k, rvec_index;

//  double num_u_inner_knots = uknot-2*(ffd_box.i_degree-1);
//  double num_v_inner_knots = vknot-2*(ffd_box.j_degree-1);
//  double num_w_inner_knots = wknot-2*(ffd_box.k_degree-1);

  double x_spacing = (ffd_box.i_range[1] - ffd_box.i_range[0]);
  double y_spacing = (ffd_box.j_range[1] - ffd_box.j_range[0]);
  double z_spacing = (ffd_box.k_range[1] - ffd_box.k_range[0]);

  /* fill in rvec for i direction knot sequence rvec */
  for (i = 0; i < uknot; i++) {
    if (i < ffd_box.i_degree+1)
      ffd_box.rvec[i] = 0.;
    else if (i < uknot-ffd_box.i_degree-1)
      ffd_box.rvec[i] = ffd_box.u[i-ffd_box.i_degree];
    else
      ffd_box.rvec[i] = 1.;
  }

  /* fill in rvec for j direction knot sequence rvec */
  for (j = 0; j < vknot; j++) {
    if (j < ffd_box.j_degree+1)
      ffd_box.rvec[uknot+j] = 0.;
    else if (j < vknot-ffd_box.j_degree-1)
      ffd_box.rvec[uknot+j] = ffd_box.v[j-ffd_box.j_degree];
    else
      ffd_box.rvec[uknot+j] = 1.;
  }

  /* fill in rvec for k direction knot sequence */
  for (k = 0; k < wknot; k++) {
    if (k < ffd_box.k_degree+1)
      ffd_box.rvec[uknot+vknot+k] = 0.;
    else if (k < wknot-ffd_box.k_degree-1)
      ffd_box.rvec[uknot+vknot+k] = ffd_box.w[k-ffd_box.k_degree];
    else
      ffd_box.rvec[uknot+vknot+k] = 1.;
  }

  /* fill in rvec for control point locations */
  for  (k = 0; k < ffd_box.k_cp; k++) {
    for (j = 0; j < ffd_box.j_cp; j++) {
      for (i = 0; i < ffd_box.i_cp; i++) {
        rvec_index = nknot + 3 * (i + j * ffd_box.i_cp + k * ffd_box.i_cp * ffd_box.j_cp);
        ffd_box.rvec[rvec_index  ] = ffd_box.i_range[0] + ffd_box.u[i] * x_spacing;
        ffd_box.rvec[rvec_index+1] = ffd_box.j_range[0] + ffd_box.v[j] * y_spacing;
        ffd_box.rvec[rvec_index+2] = ffd_box.k_range[0] + ffd_box.w[k] * z_spacing;
      }
    }
  }

  return status;
}


int makeUvArray(
  ffdBoxData& ffd_box,
  uvArrayData& uvData,
  char *message)
{
  int    status = EGADS_SUCCESS;
  int    i, j, atype, alen, iUV, iedge, iloop, ideg, icp, jdeg, jcp;
  int    oclass, mtype, nnode, nloop, nloop_edge, *senses;
  double xyz[3], uvw[3], data[4], uv, uvrange[4];
  const char *pstr;
  bool G1 = false;
  ego ecurve, epcurve, esurf, *eloops, *eloop_edges, *enodes;

  ROUTINE(makeUvArray);

  /* --------------------------------------------------------------- */

  MALLOC(uvData.uvw, double, 3*uvData.nCP); // 3 times number of CPs in spline

  /* set an initial value as 'not specified' */
  for (i=0; i < 3*uvData.nCP; i+=3) {
    uvData.uvw[i  ] = -2;
    uvData.uvw[i+1] = -2;
    uvData.uvw[i+2] = -2;
  }


  if (uvData.topo->oclass == NODE) {

    status = EG_attributeRet(uvData.topo, "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
    if (status == EGADS_SUCCESS) {
      /* Fix the Node */
      uvData.uvw[0] = -1;
      uvData.uvw[1] = -1;
      uvData.uvw[2] = -1;
    }

  } else if (uvData.topo->oclass == EDGE) {

    status = EG_attributeRet(uvData.topo, "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
    if (status == EGADS_SUCCESS) {

      for (i=0; i < 3*uvData.nCP; i+=3) {
        /* Fix the Edge */
        uvData.uvw[i+0] = -1;
        uvData.uvw[i+1] = -1;
        uvData.uvw[i+2] = -1;
      }

    } else {

      ideg = uvData.ivec[1];
      icp  = uvData.ivec[2];

      // Check the Nodes of the Edge
      status = EG_getTopology(uvData.topo, &ecurve, &oclass, &mtype, data, &nnode, &enodes, &senses);
      CHECK_STATUS(EG_getTopology);
      if (status == EGADS_SUCCESS) {
        /* Check for G1 continuity at the Edge beg */
        status = EG_attributeRet(enodes[0], "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
        if (status == EGADS_SUCCESS) {
          G1 = strcasecmp(pstr,"G1") == 0;
          for (i=0; i < icp; i++) {
            if (uvData.uk[ideg+i] > data[0]) {
              if (G1) {
                uvData.uvw[3*i+0] = -1;
                uvData.uvw[3*i+1] = -1;
                uvData.uvw[3*i+2] = -1;
              }
              break;
            }
            uvData.uvw[3*i+0] = -1;
            uvData.uvw[3*i+1] = -1;
            uvData.uvw[3*i+2] = -1;
          }
        }

        /* Check for G1 continuity at the Edge end */
        status = EG_attributeRet(enodes[1], "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
        if (status == EGADS_SUCCESS) {
          G1 = strcasecmp(pstr,"G1") == 0;
          for (i=icp-1; i >= 0; i--) {

            if (uvData.uk[i+1] < data[1]) {
              if (G1) {
                uvData.uvw[3*i+0] = -1;
                uvData.uvw[3*i+1] = -1;
                uvData.uvw[3*i+2] = -1;
              }
              break;
            }

            uvData.uvw[3*i+0] = -1;
            uvData.uvw[3*i+1] = -1;
            uvData.uvw[3*i+2] = -1;
          }
        }
      }
    }

  } else if (uvData.topo->oclass == FACE) {

    status = EG_attributeRet(uvData.topo, "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
    if (status == EGADS_SUCCESS) {
      for (i=0; i < 3*uvData.nCP; i+=3) {
        /* Fix the Face */
        uvData.uvw[i+0] = -1;
        uvData.uvw[i+1] = -1;
        uvData.uvw[i+2] = -1;
      }
    } else {

      status = EG_getTopology(uvData.topo, &esurf, &oclass, &mtype, uvrange, &nloop, &eloops, &senses);
      CHECK_STATUS(EG_getTopology);

      for (iloop = 0; iloop < nloop; iloop++) {

        status = EG_getTopology(eloops[iloop], &esurf, &oclass, &mtype, data, &nloop_edge, &eloop_edges, &senses);
        CHECK_STATUS(EG_getTopology);

        for (iedge = 0; iedge < nloop_edge; iedge++) {

          if (eloop_edges[iedge]->mtype == DEGENERATE) {
            status = EG_getTopology(eloop_edges[iedge], &ecurve, &oclass, &mtype, data, &nnode, &enodes, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_attributeRet(enodes[0], "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
          } else {
            status = EG_attributeRet(eloop_edges[iedge], "ffdFixed", &atype, &alen, NULL, NULL, &pstr);
          }

          if (status == EGADS_NOTFOUND) continue;

          if ((atype != ATTRSTRING) ||
              (pstr == NULL) ||
              (strcasecmp(pstr,"G0") != 0 &&
               strcasecmp(pstr,"G1") != 0)) {
            snprintf(message, LENMSG, "'ffdFixed' Edge attribute must be a string 'G0' or 'G1'\n");
            status = EGADS_ATTRERR;
            goto cleanup;
          }

          G1 = strcasecmp(pstr,"G1") == 0;

          epcurve = eloop_edges[nloop_edge+iedge];

          status = EG_isIsoPCurve(epcurve, &iUV, &uv, NULL);
          CHECK_STATUS(EG_isIsoPCurve);

          ideg = uvData.ivec[1];
          icp  = uvData.ivec[2];
          jdeg = uvData.ivec[4];
          jcp  = uvData.ivec[5];

          if (iUV == UISO) {

            if ( fabs(uv - uvrange[0]) < fabs(uv - uvrange[1]) ) {

              // U-min ISO curve
              for (j=0; j < jcp; j++) {
                for (i=0; i < icp; i++) {
                  if (uvData.uk[ideg+i] > uv) {
                    if (G1) {
                      uvData.uvw[3*(j*icp+i)+0] = -1;
                      uvData.uvw[3*(j*icp+i)+1] = -1;
                      uvData.uvw[3*(j*icp+i)+2] = -1;
                    }
                    break;
                  }
                  uvData.uvw[3*(j*icp+i)+0] = -1;
                  uvData.uvw[3*(j*icp+i)+1] = -1;
                  uvData.uvw[3*(j*icp+i)+2] = -1;
                }
              }

            } else {

              // U-max ISO curve
              for (j=0; j < jcp; j++) {
                for (i=icp-1; i >= 0; i--) {
                  if (uvData.uk[i+1] < uv) {
                    if (G1) {
                      uvData.uvw[3*(j*icp+i)+0] = -1;
                      uvData.uvw[3*(j*icp+i)+1] = -1;
                      uvData.uvw[3*(j*icp+i)+2] = -1;
                    }
                    break;
                  }
                  uvData.uvw[3*(j*icp+i)+0] = -1;
                  uvData.uvw[3*(j*icp+i)+1] = -1;
                  uvData.uvw[3*(j*icp+i)+2] = -1;
                }
              }

            }

          } else { // VISO

            if ( fabs(uv - uvrange[2]) < fabs(uv - uvrange[3]) ) {

              // V-min ISO curve
              for (j=0; j < jcp; j++) {
                if (uvData.vk[jdeg+j] > uv) {
                  if (G1) {
                    for (i=0; i < icp; i++) {
                      uvData.uvw[3*(j*icp+i)+0] = -1;
                      uvData.uvw[3*(j*icp+i)+1] = -1;
                      uvData.uvw[3*(j*icp+i)+2] = -1;
                    }
                  }
                  break;
                }
                for (i=0; i < icp; i++) {
                  uvData.uvw[3*(j*icp+i)+0] = -1;
                  uvData.uvw[3*(j*icp+i)+1] = -1;
                  uvData.uvw[3*(j*icp+i)+2] = -1;
                }
              }

            } else {

              // V-max ISO curve
              for (j=jcp-1; j >= 0; j--) {
                if (uvData.vk[j+1] < uv) {
                  if (G1) {
                    for (i=0; i < icp; i++) {
                      uvData.uvw[3*(j*icp+i)+0] = -1;
                      uvData.uvw[3*(j*icp+i)+1] = -1;
                      uvData.uvw[3*(j*icp+i)+2] = -1;
                    }
                  }
                  break;
                }
                for (i=0; i < icp; i++) {
                  uvData.uvw[3*(j*icp+i)+0] = -1;
                  uvData.uvw[3*(j*icp+i)+1] = -1;
                  uvData.uvw[3*(j*icp+i)+2] = -1;
                }
              }

            }
          }
        }
      }
    }
  }


  for (i=0; i < 3*uvData.nCP; i+=3) {

    // skip anything already marked as frozen
    if ( uvData.uvw[i] == -1 ) continue;

    // bspline xyz
    xyz[0] = uvData.CP[i  ];
    xyz[1] = uvData.CP[i+1];
    xyz[2] = uvData.CP[i+2];

    uvw[0] = uvw[1] = uvw[2] = 0.5;
    // ivec and rvec are for the ffd bspline, uvw is returned
    status = spline3dInvEval(ffd_box.ivec, ffd_box.rvec, xyz, uvw);
    if (status == -3) {
      /* outside the box (ffdbox uvw by definition in [0,1] */
      uvData.uvw[i  ] = -1;
      uvData.uvw[i+1] = -1;
      uvData.uvw[i+2] = -1;
      continue;
    }
    CHECK_STATUS(spline3dInvEval);

    /* inside the box */
    uvData.uvw[i  ] = uvw[0];
    uvData.uvw[i+1] = uvw[1];
    uvData.uvw[i+2] = uvw[2];
  }


  status = EGADS_SUCCESS;

cleanup:

  return status;
}


template<class T>
int copyBSplines(
  std::vector< uvArrayData >& uvData,
  std::vector< BSpline<T> >& bspl_data)
{
  int    status = EGADS_SUCCESS;
  std::size_t i, j, rvec_size = 0;

  ROUTINE(copyBSplines);

  /* --------------------------------------------------------------- */

  bspl_data.resize(uvData.size());

  for (i = 0; i < uvData.size(); i++) {

    bspl_data[i].eIndex = uvData[i].eIndex;

    if (uvData[i].topo->oclass == NODE) {

      rvec_size = 3;
    }
    else if (uvData[i].topo->oclass == EDGE) {

      for (j = 0; j < 4; j++)
        bspl_data[i].ivec[j] = uvData[i].ivec[j];

      rvec_size = 3*bspl_data[i].ivec[2] + bspl_data[i].ivec[3];

      if (bspl_data[i].ivec[0] & 2) { // we need to add the bspline weights
        rvec_size += bspl_data[i].ivec[2]; // add length of cps
      }
    }
    else if (uvData[i].topo->oclass == FACE) {
      for (j = 0; j < 7; j++)
        bspl_data[i].ivec[j] = uvData[i].ivec[j];

      rvec_size = bspl_data[i].ivec[3] + bspl_data[i].ivec[6] +   // knots
                  3*(bspl_data[i].ivec[2]*bspl_data[i].ivec[5]);  // control points

      if (bspl_data[i].ivec[0] & 2) { // we need to add the bspline weights
        rvec_size += (bspl_data[i].ivec[2]*bspl_data[i].ivec[5]); // add length of cps
      }
    }

    MALLOC(bspl_data[i].rvec, T, rvec_size);

    for (j = 0; j < rvec_size; j++) {
      bspl_data[i].rvec[j] = uvData[i].rvec[j];
    }
  }

  status = EGADS_SUCCESS;

cleanup:

  return status;
}


/*
Move CPs stored in rvec of uvData according to input FFD box defined by ffdIvec and ffdRvec
Outputs: penalty, gradients of each FFD CP x,y,z
Modifies: uvData
*/
template <class T>
int moveGeomCps(BSpline<T>& ffdBox, uvArrayData &uvData, BSpline<T> & uv_vec)
{

  int status = EGADS_SUCCESS;
  int i, j, k;
  T   uvw[3], eval[30];

  ROUTINE(moveCpsHelper);

  /* --------------------------------------------------------------- */

  if (uvData.topo->oclass == NODE) {
    uvw[0] = uvData.uvw[0];
    uvw[1] = uvData.uvw[1];
    uvw[2] = uvData.uvw[2];

    if (uvw[0] != -1) {
      status = spline3dEval(ffdBox.ivec, ffdBox.rvec, uvw, eval);
      CHECK_STATUS(spline3dEval);

      uv_vec.rvec[0] = eval[0];
      uv_vec.rvec[1] = eval[1];
      uv_vec.rvec[2] = eval[2];
    }

  }

  else if (uvData.topo->oclass == EDGE) {

    const int icp = uv_vec.ivec[2];

    for (k=0; k < 3*icp; k+=3) {
      uvw[0] = uvData.uvw[k  ];
      uvw[1] = uvData.uvw[k+1];
      uvw[2] = uvData.uvw[k+2];

      if (uvw[0] != -1) {
        status = spline3dEval(ffdBox.ivec, ffdBox.rvec, uvw, eval);
        CHECK_STATUS(spline3dEval);

        uv_vec.rvec[uv_vec.ivec[3]+k  ] = eval[0];
        uv_vec.rvec[uv_vec.ivec[3]+k+1] = eval[1];
        uv_vec.rvec[uv_vec.ivec[3]+k+2] = eval[2];
      }
    }
  }

  else if (uvData.topo->oclass == FACE) {
    const int nknot = uv_vec.ivec[3]+uv_vec.ivec[6];
    const int icp = uv_vec.ivec[2];
    const int jcp = uv_vec.ivec[5];
    for (i = 0; i < icp; i ++) {
      for (j = 0; j < jcp; j++) {
        k = 3 * (i + icp * j);
        uvw[0] = uvData.uvw[k  ];
        uvw[1] = uvData.uvw[k+1];
        uvw[2] = uvData.uvw[k+2];

        if (uvw[0] != -1) {
          status = spline3dEval(ffdBox.ivec, ffdBox.rvec, uvw, eval);
          CHECK_STATUS(spline3dEval);

          uv_vec.rvec[nknot+k  ] = eval[0];
          uv_vec.rvec[nknot+k+1] = eval[1];
          uv_vec.rvec[nknot+k+2] = eval[2];
        }
      }
    }
  }

  status = EGADS_SUCCESS;

cleanup:
  return status;
}

/*
Move the control points of the FFD box according to inputed perturbations
And then update the rvec data for the moved surfaces and edges of the source body
*/
template <class T>
int moveCps(
   BSpline<T>& ffdBox,
   std::vector< uvArrayData > &uvNodes,
   std::vector< uvArrayData > &uvEdges,
   std::vector< uvArrayData > &uvFaces,
   std::vector< BSpline<T> > &bspl_nodes,
   std::vector< BSpline<T> > &bspl_edges,
   std::vector< BSpline<T> > &bspl_faces,
   T &penaltyOut)
{
  int    status = EGADS_SUCCESS;
  int i = 0, j = 0, k = 0;
  T total_penalty = 0;
  T penalty;

  ROUTINE(moveCps);

  /* --------------------------------------------------------------- */

  // const int nknot = ffdIvec[2]+ffdIvec[5]+ffdIvec[8]; // num ffd knots
  const int icp = ffdBox.ivec[1];
  const int jcp = ffdBox.ivec[4];
  const int kcp = ffdBox.ivec[7];
  // const int n = icp*jcp*kcp*3;

  min_p = 1e8; // reset so output is min_p of that iter

  for (i = 0; i < (int)uvNodes.size(); i++) {
    status = moveGeomCps(ffdBox, uvNodes[i], bspl_nodes[i]);
    CHECK_STATUS(moveGeomCps);
  }

  for (i = 0; i < (int)uvEdges.size(); i++) {
    status = moveGeomCps(ffdBox, uvEdges[i], bspl_edges[i]);
    CHECK_STATUS(moveGeomCps);
  }

  for (i = 0; i < (int)uvFaces.size(); i++) {
    status = moveGeomCps(ffdBox, uvFaces[i], bspl_faces[i]);
    CHECK_STATUS(moveGeomCps);
  }

  for (i = 0; i < icp; i++) {
    for (j = 0; j < jcp; j++) {
      for (k = 0; k < kcp; k++) {
        computePenalty(ffdBox.ivec, ffdBox.rvec, i, j, k, penalty);
        total_penalty += penalty;
      }
    }
  }

  penaltyOut = total_penalty;

cleanup:
  return status;
}

template <class T>
int zero_dot(int nnode, ego *enodes, int nedge, ego *eedges, int nface, ego *efaces)
{
  return EGADS_SUCCESS;
}

template <>
int zero_dot<SurrealS1>(int nnode, ego *enodes, int nedge, ego *eedges, int nface, ego *efaces)
{
  int             status = EGADS_SUCCESS;
  ego             eref, geom;
  int             oclass, mtype, *senses, nchild, *ivec = NULL, i, ni, nr, j, periodic;
  ego             *echild = NULL;
  double          data[18], tdata[2], tdata_dot[2];
  SurrealS1       node[3];
  SurrealS1 *rvec = NULL;
  double *dvec = NULL;

  ROUTINE(zero_dot);

  /* --------------------------------------------------------------- */

  for (i = 0; i < nnode; i++) {
    status = EG_evaluate(enodes[i], NULL, data);
    CHECK_STATUS(EG_evaluate);
    node[0] = data[0];
    node[1] = data[1];
    node[2] = data[2];

    status = EG_setGeometry_dot(enodes[i], NODE, 0, NULL, node);
    CHECK_STATUS(EG_setGeometry_dot);
  }


  for (i = 0; i < nedge; i++) {
    status = EG_getTopology(eedges[i], &geom, &oclass, &mtype, tdata, &nchild, &echild, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getRange(eedges[i], tdata, &periodic);
    CHECK_STATUS(EG_getRange);

    tdata_dot[0] = 0.;
    tdata_dot[1] = 0.;

    status = EG_setRange_dot(eedges[i], EDGE, tdata, tdata_dot);
    CHECK_STATUS(EG_setRange_dot);

    if (mtype == DEGENERATE) continue;
    do {

      status = EG_getGeometry(geom, &oclass, &mtype, &eref, &ivec, &dvec);
      CHECK_STATUS(EG_getGeometry);

      EG_getGeometryLen(geom, &ni, &nr);

      rvec = new SurrealS1[nr];
      for (j = 0; j < nr; j++) rvec[j] = dvec[j];

      status = EG_setGeometry_dot(geom, oclass, mtype, ivec, rvec);
      CHECK_STATUS(EG_setGeometry_dot);
      geom = eref;

      EG_free(dvec);
      delete [] rvec;

    } while (geom != NULL);
  }

  for (i = 0; i < nface; i++) {
    status = EG_getTopology(efaces[i], &geom, &oclass, &mtype, data, &nchild, &echild, &senses);
    CHECK_STATUS(EG_getTopology);

    do {
      status = EG_getGeometry(geom, &oclass, &mtype, &eref, &ivec, &dvec);
      CHECK_STATUS(EG_getGeometry);

      EG_getGeometryLen(geom, &ni, &nr);

      rvec = new SurrealS1[nr];
      for (j = 0; j < nr; j++) rvec[j] = dvec[j];

      status = EG_setGeometry_dot(geom, oclass, mtype, ivec, rvec);
      CHECK_STATUS(EG_setGeometry_dot);
      geom = eref;

      EG_free(dvec);
      delete [] rvec;

    } while (geom != NULL);
  }

  status = EGADS_SUCCESS;

cleanup:
  return status;

}

/*
originalBody is body to perturb
uvArrayData vectors should be updated with appropriate ivec and rvec so new geometry can be made
*/
template <class T>
int replaceFacesAndCurves(
   ego source_body,
   std::vector< BSpline<T> > &bspl_nodes,
   std::vector< BSpline<T> > &bspl_edges,
   std::vector< BSpline<T> > &bspl_faces,
   ego &deformed_body,
   char *message)
{
  int    i = 0;
  int    status = EGADS_SUCCESS;
  ego    context;
  int    oclass, mtype, *senses, nnode, nloop;
  ego    *efaces = NULL, *eloops = NULL, *eedges = NULL, *enodes = NULL, *eedge_nodes = NULL;
  double data[18];
  double tdata[2];
  T      node[18];
  int    nNewNode, nNewCurv, nNewSurf;
  ego    *eNewNodes=NULL, *eNewCurvs=NULL, *eNewSurfs=NULL;

  ROUTINE(replaceFacesAndCurves);

  /* --------------------------------------------------------------- */

  status = EG_getContext(source_body, &context);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(source_body, NULL, NODE, &nNewNode, &enodes);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(source_body, NULL, EDGE, &nNewCurv, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(source_body, NULL, FACE, &nNewSurf, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  status = zero_dot<T>(nNewNode, enodes, nNewCurv, eedges, nNewSurf, efaces);
  CHECK_STATUS(zero_dot);

  /* get the current Nodes */
  MALLOC(eNewNodes, ego, nNewNode);
  for (i = 0; i < nNewNode; i++) {
    eNewNodes[i] = enodes[i];
  }

  /* get the current Curves */
  MALLOC(eNewCurvs, ego, nNewCurv);
  for (i = 0; i < nNewCurv; i++) {
    status = EG_getTopology(eedges[i], &eNewCurvs[i], &oclass, &mtype, tdata, &nnode, &eedge_nodes, &senses);
    CHECK_STATUS(EG_getTopology);
  }

  /* Create new Nodes */
  for (i = 0; i < (int)bspl_nodes.size(); i++) {
    node[0] = bspl_nodes[i].rvec[0];
    node[1] = bspl_nodes[i].rvec[1];
    node[2] = bspl_nodes[i].rvec[2];

    status = EG_makeTopology(context, NULL, NODE, 0, node, 0, NULL, NULL, &eNewNodes[bspl_nodes[i].eIndex]);
    CHECK_STATUS(EG_makeTopology);
  }

  /* Create new Curves */
  for (i = 0; i < (int)bspl_edges.size(); i++) {
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, bspl_edges[i].ivec, bspl_edges[i].rvec, &eNewCurvs[bspl_edges[i].eIndex]);
    CHECK_STATUS(EG_getTopology);
  }


  // WireBody
  if (nNewSurf == 0) {
    /* build the body */
    status = buildBody(source_body,
                       nNewNode, eNewNodes,
                       nNewCurv, eNewCurvs,
                       0, NULL,
                       &deformed_body,
                       LENMSG, message);
    CHECK_STATUS(buildBody);
    goto cleanup;
  }

  MALLOC(eNewSurfs, ego, nNewSurf);
  for (i = 0; i < nNewSurf; i++) {
    status = EG_getTopology(efaces[i], &eNewSurfs[i], &oclass, &mtype, data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);
  }

  /* make bsurface from the ivec and rvec in builder */
  for (i = 0; i < (int)bspl_faces.size(); i++) {
    status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bspl_faces[i].ivec, bspl_faces[i].rvec, &eNewSurfs[bspl_faces[i].eIndex]);
    CHECK_STATUS(EG_makeGeometry);
  }

  /* build the body */
  status = buildBody(source_body,
                     nNewNode, eNewNodes,
                     nNewCurv, eNewCurvs,
                     nNewSurf, eNewSurfs,
                     &deformed_body,
                     LENMSG, message);
  CHECK_STATUS(buildBody);

  // printf("new body tolerance = %le \n", tolerance);

  status = EGADS_SUCCESS;

cleanup:

  FREE(enodes);
  FREE(eedges);
  FREE(efaces);

  FREE(eNewNodes);
  FREE(eNewCurvs);
  FREE(eNewSurfs);

  return status;
}



/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */
extern "C"
int
udpExecute(
          ego emodel,                  /* (in)  input model */
          ego  *ebody,                 /* (out) Body pointer */
          int  *nMesh,                 /* (out) number of associated meshes */
          char *string[])              /* (out) error message */
{
  int    status = EGADS_SUCCESS;
  ego    context=NULL;
  ego    eref, geom, source_body = NULL, *ebodies = NULL;
  double data[18], tdata[2], penalty;
  int    oclass, mtype, *senses, num_ffd_knots;
  int    identified_boxes = 0;
  ego    *efaces = NULL, *eedges = NULL, *enodes = NULL, *echild;
  int    nchild, nface, nedge, nnode, nbody, atype, alen;
  int    i, k;
  udpDotCache_T *cache;
  const int     *pint;
  char          *message=NULL;
  void    *realloc_temp=NULL;
  udp_T   *udps = *Udps;
  udp_T   *udp = &udps[0];
  FILE    *fp = NULL;

  BSpline<double>                ffd_vec;
  std::vector< BSpline<double> > bspl_nodes;
  std::vector< BSpline<double> > bspl_edges;
  std::vector< BSpline<double> > bspl_faces;

  ROUTINE(udpExecute);

  /* --------------------------------------------------------------- */

#ifdef DEBUG
  printf("udpExecute(emodel=%llx)\n", (long long)emodel);
#endif

  /* default return values */
  *ebody  = NULL;
  *nMesh  = 0;
  *string = NULL;

  /* setup the cache for tracking dot changes */
  cache = new udpDotCache_T; // need c++ allocation
  ffdBoxData& ffdBoxStruct = cache->ffdBox;
  std::vector< uvArrayData >& uvNodes = cache->uvNodes;
  std::vector< uvArrayData >& uvEdges = cache->uvEdges;
  std::vector< uvArrayData >& uvFaces = cache->uvFaces;

  MALLOC(message, char, LENMSG);
  message[0] = '\0';

  //printf("udpExecute(emodel=%p)\n", emodel);
  status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                          data, &nbody, &ebodies, &senses);
  CHECK_STATUS(EG_getTopology);

  if (oclass != MODEL) {
    snprintf(message, LENMSG, "Expecting a Model\n");
    status = EGADS_NOTMODEL;
    goto cleanup;
  } else if (nbody < 2) {
    snprintf(message, LENMSG, "Expecting Model to contain at least two Bodys (not %d)\n", nbody);
    status = EGADS_NOTBODY;
    goto cleanup;
  } else if (DESIGN_PARAM_ARG(udp).ncol != 3) {
    snprintf(message, LENMSG, "Expecting three columns (x,y,z) for design parameter matrix. There are %d columns \n", DESIGN_PARAM_ARG(udp).ncol);
    status = EGADS_RANGERR;
    goto cleanup;
  }

  /* Cache copy of arguments for future use.
   * This also increments numUdp and copies udps[numUdp] to udps[numUdp].
   * Caching should only be performed after checking for valid inputs.
   */
  status = cacheUdp(emodel);
  CHECK_STATUS(cacheUdp);
  udp = &udps[numUdp];
  udp->data = (void*)cache;

  /* make enough room for the penalty */
  PENALTY_ARG(udp).size = 1;
  PENALTY_ARG(udp).nrow = 1;
  PENALTY_ARG(udp).ncol = 1;

  RALLOC(PENALTY_ARG(udp).val, double, PENALTY_ARG(udp).size);
  RALLOC(PENALTY_ARG(udp).dot, double, PENALTY_ARG(udp).size);

  PENALTY_VAL(udp) = 0;
  PENALTY_DOT(udp) = 0;

  /* Look at input bodies */
  status = EG_getContext(emodel, &context);
  CHECK_STATUS(EG_getContext);

#ifdef DEBUG
  printf("udpExecute(context=%llx)\n", (long long)context);
  for (irow = 0; irow < DESIGN_PARAM_ARG(udp).nrow; irow++) {
    for (icol = 0; icol < DESIGN_PARAM_ARG(udp).ncol; icol++) {
      printf("design_param_val(     0,%d,%d) = %f\n", irow, icol, DESIGN_PARAM_VAL(    0,irow,icol));
      printf("design_param_dot( 0,%d,%d) = %f\n", irow, icol, DESIGN_PARAM_DOT(0,irow,icol));
    }
  }
#endif


  /* There should be one body to deform (deform body), and one ffd box body. */
  for (i = 0; i < nbody; i++) {
    status = checkFFDbox(context, ebodies[i], ffdBoxStruct, message);
    if (status == EGADS_SUCCESS)
    {
      if (identified_boxes > 1) {
        snprintf(message, LENMSG, "Too many FFD boxes identified. \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      } else {
        status = makeFFDBox(ffdBoxStruct);
        CHECK_STATUS(makeFFDBox);

        if ((ffdBoxStruct.i_cp) * (ffdBoxStruct.j_cp) * (ffdBoxStruct.k_cp) * 3 != DESIGN_PARAM_ARG(udp).size) {
          snprintf(message, LENMSG, "design parameters array length does not match the number of control points x,y,z values of the ffd box. \
             design parameters = %d, controls points * 3 = %d \n.", ffdBoxStruct.i_cp * ffdBoxStruct.j_cp * ffdBoxStruct.k_cp * 3, DESIGN_PARAM_ARG(udp).size);
          status = EGADS_RANGERR;
          goto cleanup;
        }

        //num_ffd_knots = ffdBoxStruct.ivec[2] + ffdBoxStruct.ivec[5] + ffdBoxStruct.ivec[8];

        identified_boxes++;
      }
    }
    else if (status == EGADS_RANGERR)
    {
      goto cleanup;
    }
    else
    {
      /* The body is not an FFD body */
      if (source_body != NULL) {
        snprintf(message, LENMSG, "More than one deform body identified. \n"
            "FFD box must have one each of ffd_uMin, ffd_vMin, and ffd_wMin attributes on faces. \n");
        goto cleanup;
      } else {
        status = EG_copyObject(ebodies[i], NULL, &source_body);
        CHECK_STATUS(EG_copyObject);

        status = checkDeformBody(context, source_body, message);
        CHECK_STATUS(checkDeformBody);
      }
    }
  }

  status = EG_getBodyTopos(source_body, NULL, FACE, &nface, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(source_body, NULL, EDGE, &nedge, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(source_body, NULL, NODE, &nnode, &enodes);
  CHECK_STATUS(EG_getBodyTopos);

  for (i = 0; i < nnode; i++) {
    status = EG_attributeRet(enodes[i], "ffdID", &atype, &alen, &pint, NULL, NULL);
    if (status == EGADS_SUCCESS) {

      uvNodes.push_back(uvArrayData());
      k = uvNodes.size()-1;
      uvNodes[k].topo = enodes[i];
      uvNodes[k].eIndex = i;

      status = EG_evaluate(enodes[i], NULL, data);
      CHECK_STATUS(EG_getTopology);

      MALLOC(uvNodes[k].rvec, double, 3);

      uvNodes[k].rvec[0] = data[0];
      uvNodes[k].rvec[1] = data[1];
      uvNodes[k].rvec[2] = data[2];

      uvNodes[k].nCP = 1;
      uvNodes[k].CP  = uvNodes[k].rvec;

      status = makeUvArray(ffdBoxStruct, uvNodes[k], message);
      CHECK_STATUS(makeUvArray);
    }
  }

  for (i = 0; i < nedge; i++) {

    status = EG_getTopology(eedges[i], &geom, &oclass, &mtype, tdata, &nchild, &echild, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_attributeRet(eedges[i], "ffdID", &atype, &alen, &pint, NULL, NULL);
    if (status == EGADS_SUCCESS) {
#if 0
      status = EG_attributeDel(eedges[i], NULL);
      CHECK_STATUS(EG_attributeDel);
#endif

      status = EG_attributeAdd(eedges[i], DEFORMBODY, ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);

      status = EG_attributeAdd(eedges[i], "ffdID", ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);

      uvEdges.push_back(uvArrayData());
      k = uvEdges.size()-1;
      uvEdges[k].topo = eedges[i];
      uvEdges[k].eIndex = i;

      status = EG_getGeometry(geom, &oclass, &mtype, &eref, &uvEdges[k].ivec, &uvEdges[k].rvec);
      CHECK_STATUS(EG_getGeometry);

      uvEdges[k].nCP = uvEdges[k].ivec[2];
      uvEdges[k].CP  = uvEdges[k].rvec + uvEdges[k].ivec[3];
      uvEdges[k].uk  = uvEdges[k].rvec;

      if (mtype != BSPLINE) {
        snprintf(message, LENMSG, "Body Edge %d inside FFD box is not a B-Spline!\n", i+1);
        status = EGADS_RANGERR;
        goto cleanup;
      }

      status = makeUvArray(ffdBoxStruct, uvEdges[k], message);
      CHECK_STATUS(makeUvArray);

    } else {
#if 0
      status = EG_attributeDel(eedges[i], NULL);
      CHECK_STATUS(EG_attributeDel);

      status = EG_attributeAdd(eedges[i], DEFORMBODY, ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);
#endif
    }
  }

  for (i = 0; i < nface; i++) {
    status = EG_attributeRet(efaces[i], "ffdID", &atype, &alen, &pint, NULL, NULL);

    if (status == EGADS_SUCCESS) {

      status = EG_getTopology(efaces[i], &geom, &oclass, &mtype, data, &nchild, &echild, &senses);
      CHECK_STATUS(EG_getTopology);

      uvFaces.push_back(uvArrayData());
      k = uvFaces.size()-1;
      uvFaces[k].topo = efaces[i];
      uvFaces[k].eIndex = i;

      status = EG_getGeometry(geom, &oclass, &mtype, &eref, &uvFaces[k].ivec, &uvFaces[k].rvec);
      CHECK_STATUS(EG_getGeometry);

      if (mtype != BSPLINE) { // check the geometry surface and not the face
        snprintf(message, LENMSG, "body has faces in FFD box which are not B-Splines. \n");
        status = EGADS_RANGERR;
        goto cleanup;
      }

      uvFaces[k].nCP = uvFaces[k].ivec[2]*uvFaces[k].ivec[5];
      uvFaces[k].CP  = uvFaces[k].rvec + uvFaces[k].ivec[3]+uvFaces[k].ivec[6];
      uvFaces[k].uk  = uvFaces[k].rvec;
      uvFaces[k].vk  = uvFaces[k].rvec + uvEdges[k].ivec[3];

      status = makeUvArray(ffdBoxStruct, uvFaces[k], message);
      CHECK_STATUS(makeUvArray);
    }
  }


  status = copyBSplines(uvNodes, bspl_nodes);
  CHECK_STATUS(copyBSplines);

  status = copyBSplines(uvEdges, bspl_edges);
  CHECK_STATUS(copyBSplines);

  status = copyBSplines(uvFaces, bspl_faces);
  CHECK_STATUS(copyBSplines);


  for ( i = 0; i < 9; i++) ffd_vec.ivec[i] = ffdBoxStruct.ivec[i];

  num_ffd_knots = ffd_vec.ivec[2] + ffd_vec.ivec[5] + ffd_vec.ivec[8];

  MALLOC(ffd_vec.rvec, double, num_ffd_knots + ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp*3);
  for ( i = 0; i < num_ffd_knots + ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp*3; i++) ffd_vec.rvec[i] = ffdBoxStruct.rvec[i];

  for (i = 0; i < DESIGN_PARAM_ARG(udp).nrow; i++) {
    ffd_vec.rvec[num_ffd_knots + 3 * i    ] += DESIGN_PARAM_VAL(udp, i, 0);
    ffd_vec.rvec[num_ffd_knots + 3 * i + 1] += DESIGN_PARAM_VAL(udp, i, 1);
    ffd_vec.rvec[num_ffd_knots + 3 * i + 2] += DESIGN_PARAM_VAL(udp, i, 2);
  }

  if (WRITEXYZ_VAL(udp) != 0) {
    /* print out the updated FFD control point positions for plotting */
    fp = fopen("ffdbox.xyz", "w");
    if (fp == NULL) {
      snprintf(message, LENMSG, "Failed to open ffdbox.xyz!\n");
      status = EGADS_WRITERR;
      goto cleanup;
    }
    fprintf(fp,"%d 0 xyz\n", ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp);
    for (i=0; i < ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp*3; i+=3) {
      fprintf(fp, "%.18e %.18e %.18e \n", ffd_vec.rvec[num_ffd_knots+i], ffd_vec.rvec[num_ffd_knots+i+1], ffd_vec.rvec[num_ffd_knots+i+2]);
    }
    fclose(fp); fp = NULL;
  }

  status = moveCps(ffd_vec, uvNodes, uvEdges, uvFaces, bspl_nodes, bspl_edges, bspl_faces, penalty);
  CHECK_STATUS(moveCps);

  PENALTY_VAL(udp) = penalty;

  status = replaceFacesAndCurves(source_body,
                                 bspl_nodes, bspl_edges, bspl_faces,
                                 *ebody, message);
  CHECK_STATUS(replaceFacesAndCurves);

#if 0
  status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(*ebody, NULL, NODE, &nnode, &enodes);
  CHECK_STATUS(EG_getBodyTopos);

  /* attribute sourceBody nodes */
  for (i = 0; i < nnode; i++) {
    // delete attributes --> issue with replace faces and _body attributes for edges and nodes
    status = EG_attributeDel(enodes[i], NULL);
    CHECK_STATUS(EG_attributeDel);
  }
#endif

  /* Put unchanging things in the cache */
  cache->source_body = source_body;
  cache->deformed_body = NULL;
  cache = NULL;

  /* remember this model (body) */
  udps[numUdp].ebody = *ebody;

#ifdef DEBUG
  printf("udpExecute -> *ebody=%lx\n", (long)(*ebody));
#endif

cleanup:
  if (strlen(message) > 0) {
    *string = message;
    printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
    printf("ERROR: status = %d \n", status);
    // freePrivateData(&ffdBox);
    FREE(message);
    *string = udpErrorStr(status);
  } else {
    FREE(message);
  }

  delete cache;

  if (fp != NULL) fclose(fp);

  /* remove all the construction parts first */
  int outLevel = EG_setOutLevel(context, 0);
  EG_deleteObject(context);
  EG_setOutLevel(context, outLevel);

  FREE(enodes);
  FREE(eedges);
  FREE(efaces);

  return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                      *
 ************************************************************************
 */
extern "C"
int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
  int    status = EGADS_SUCCESS;
  int    iudp = 0, judp, ipnt, nchild, stride;
  double point[18], point_dot[18];
  ego    eent, *echildren=NULL;
  SurrealS<1> *pts=NULL, *rdata=NULL;
  udpDotCache_T *cache;
  udp_T   *udp = NULL;

  BSpline<SurrealS1>                ffd_vec;
  std::vector< BSpline<SurrealS1> > bspl_nodes;
  std::vector< BSpline<SurrealS1> > bspl_edges;
  std::vector< BSpline<SurrealS1> > bspl_faces;
  SurrealS1 penalty;

  int    i, j;
  int    num_ffd_knots;
  char   *message=NULL;

  ROUTINE(udpSensitivity);

  /* --------------------------------------------------------------- */

  /* check that ebody matches one of the ebodys */
  iudp = 0;
  for (judp = 1; judp <= numUdp; judp++) {
    if (ebody == udps[judp].ebody) {
      iudp = judp;
      break;
    }
  }
  if (iudp <= 0) {
    return EGADS_NOTMODEL;
  }

#ifdef DEBUG
  int n, ndesignparam, npenalty;

  /* get arguments */
  ndesignparam  = udps[iudp].arg[0].size;
  npenalty = udps[iudp].arg[1].size;

  printf("udpFFDPerturb.udpSensitivity\n");
  // TODO change to i,j
  for (n = 0; n < ndesignparam; n++) {
    printf("DESIGN_PARAM_VAL(%d,%d) = %f\n", udp, n, DESIGN_PARAM_VAL(udp,n) );
  }
  for (n = 0; n < ndesignparam; n++) {
    printf("POLYDESIGN_PARAM_DOT(%d,%d) = %f\n", udp, n, DESIGN_PARAM_DOT(udp,n) );
  }
  for (n = 0; n < npenalty; n++) {
    printf("PENALTY_VAL(%d,%d) = %f\n", udp, n, PENALTY_VAL(udp,n) );
  }
  for (n = 0; n < nparam; n++) {
    printf("PENALTY_DOT(%d,%d) = %f\n", udp, n, PENALTY_DOT(udp,n) );
  }
  printf("\n");
#endif

  /* get the cache */
  udp = &udps[iudp];
  cache = (udpDotCache_T*)udp->data;
  ffdBoxData& ffdBoxStruct = cache->ffdBox;

  MALLOC(message, char, LENMSG);
  message[0] = '\0';

  /* build the sensitivity if needed */
  if (udp->ndotchg > 0) {

    for ( i = 0; i < 9; i++) ffd_vec.ivec[i] = ffdBoxStruct.ivec[i];

    num_ffd_knots = ffd_vec.ivec[2] + ffd_vec.ivec[5] + ffd_vec.ivec[8];

    MALLOC(ffd_vec.rvec, SurrealS1, num_ffd_knots + ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp*3);
    for ( i = 0; i < num_ffd_knots + ffdBoxStruct.i_cp*ffdBoxStruct.j_cp*ffdBoxStruct.k_cp*3; i++) ffd_vec.rvec[i] = ffdBoxStruct.rvec[i];

    for (i = 0; i < DESIGN_PARAM_ARG(udp).nrow; i++) {
      ffd_vec.rvec[num_ffd_knots + 3 * i    ] += DESIGN_PARAM_SURREAL(udp, i, 0);
      ffd_vec.rvec[num_ffd_knots + 3 * i + 1] += DESIGN_PARAM_SURREAL(udp, i, 1);
      ffd_vec.rvec[num_ffd_knots + 3 * i + 2] += DESIGN_PARAM_SURREAL(udp, i, 2);
    }

    status = copyBSplines(cache->uvNodes, bspl_nodes);
    CHECK_STATUS(copyBSplines);

    status = copyBSplines(cache->uvEdges, bspl_edges);
    CHECK_STATUS(copyBSplines);

    status = copyBSplines(cache->uvFaces, bspl_faces);
    CHECK_STATUS(copyBSplines);

    status = moveCps(ffd_vec,
                     cache->uvNodes, cache->uvEdges, cache->uvFaces,
                     bspl_nodes, bspl_edges, bspl_faces, penalty);
    CHECK_STATUS(moveCps);

    if (npnt == 0) {
      for (i = 0; i < udps[numUdp].arg[1].nrow; i++) {
        for (j = 0; j < udps[numUdp].arg[1].ncol; j++) {
          // printf("penalty = %.12f \n", penalty.value());
          // printf("penalty deriv = %.12f \n", penalty.deriv());
          PENALTY_VAL(udp) = penalty.value();
          PENALTY_DOT(udp) = penalty.deriv();
        }
      }
      status = EGADS_SUCCESS;
      goto cleanup;
    } else {
      EG_deleteObject(cache->deformed_body);

      status = replaceFacesAndCurves(cache->source_body,
                                     bspl_nodes, bspl_edges, bspl_faces,
                                     cache->deformed_body, message);
      CHECK_STATUS(replaceFacesAndCurves);

      status = EG_copyGeometry_dot(cache->deformed_body, NULL, NULL, ebody);
      CHECK_STATUS(EG_copyGeometry_dot);
    }
  }


  /* find the ego entity */
  if (entType == OCSM_NODE) {
    status = EG_getBodyTopos(ebody, NULL, NODE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 0;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else if (entType == OCSM_EDGE) {
    status = EG_getBodyTopos(ebody, NULL, EDGE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 1;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else if (entType == OCSM_FACE) {
    status = EG_getBodyTopos(ebody, NULL, FACE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 2;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else {
    printf("udpSensitivity: bad entType=%d\n", entType);
    status = EGADS_GEOMERR;
    goto cleanup;
  }


  /* get the velocities from the entity */
  for (ipnt = 0; ipnt < npnt; ipnt++) {
    status = EG_evaluate_dot(eent, &(uvs[stride*ipnt]), NULL, point, point_dot);
    CHECK_STATUS(EG_evaluate_dot);

    /* return the point velocity */
    vels[3*ipnt  ] = point_dot[0];
    vels[3*ipnt+1] = point_dot[1];
    vels[3*ipnt+2] = point_dot[2];
  }

cleanup:
  if (strlen(message) > 0) {
   printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
   printf("ERROR: status = %d \n", status);
  }

  delete [] pts;
  FREE( rdata );
  FREE(echildren);
  FREE(message);

  return status;

}
