/*
 *************************************************************************
 *                                                                       *
 * udpFFDsplines -- udp file to generate b-spline geometry for FFDdeform *
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

/*
CSM file

Defining the FFD box:
- must be a a cube shape (define with box)
- Must have ffd_uMin, ffd_vMin, ffd_wMin attributes for each direction (x,y,z) on faces of the box.
- Each of these attributes must have a value in the format "degree_of_spline; number_of_control_points"
- By default there are 2 control points.
- The outer frame of control points of the FFD box is fixed, and cannot be perturbed.
- The number of control points in the attribute is added to the default 2 control points,
  so a value of 1 for number of control points will result in 3 total control points, 1 of which
  is interior and can be perturbed.
- The degree of the spline must follow the relationship of: degree <= 1 + number_of_control_points.
*/

#include <utility> // std::swap
#include <vector>

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
 #define NUMUDPARGS 2
 #include "udpUtilities.h"
} /* extern "C" */

#define MINCP( UDP) ((int *)((UDP)->arg[0].val))[0]
#define MINDEG(UDP) ((int *)((UDP)->arg[1].val))[0]

/* data about possible arguments */
static const char *argNames[NUMUDPARGS] = {"mincp", "mindeg"};
static int         argTypes[NUMUDPARGS] = {ATTRINT, ATTRINT };
static int         argIdefs[NUMUDPARGS] = {11     , 3       };
static double      argDdefs[NUMUDPARGS] = {0      , 0       };

#include "spline3d_eval.hpp"

#define LENMSG 100
#define DEFORMBODY "ffdDeformBody"
#define FFDID      "ffdID"

struct convert_T {
  int index;
  double range[4];
};

extern "C" {
/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
} /* extern "C" */


extern "C"
int buildBody(ego eOrigBody,                       /* (in) original body to rebuild */
              const int nNewNode, ego *eNewNodes,  /* (in) mixture of orig and new Nodes    (Node BodyIndex ordered) */
              const int nNewCurv, ego *eNewCurvs,  /* (in) mixture of orig and new Curves   (Edge BodyIndex ordered) */
              const int nNewSurf, ego *eNewSurfs,  /* (in) mixture of orig and new Surfaces (Face BodyIndex ordered) */
              ego *eNewBody,                       /* (out) the new body with all geometry replaced */
              int lenmsg,                          /* (in) error message length */
              char *message);                      /* (out) error message */


/*
Returns EGADS_NOTFOUND (status = -1) if the ebody is not an FFD body.
Returns EGADS_RANGERR (status = -16) if the ebody has attributes for an FFD body, but
the attributes are incorrect or invalid in some way.
Returns EGADS_SUCCESS if the ebody is an FFD body with correct attributes.
Fills in the ffdBoxData structure if successful.
*/
static int
checkFFDbox(ego ebody,           /* (in) potential ffd box */
            char *message)       /* (out) error message */
{
  int    status = EGADS_SUCCESS;
  ego    eref, *efaces = NULL, *eloops = NULL, *eedges = NULL;
  double data[18];
  int    degree, ncp;
  int    oclass, mtype, atype, nface, nloop, nedge, *senses;
  int    alen, i, j, nfaceNeighbor;
  bool   uMin_attr = false, vMin_attr = false, wMin_attr = false;
  ego    *efaceNeighbors = NULL, *faceCandidates = NULL;
  const double  *pdouble;

  ROUTINE(checkFFDbox);

  /* --------------------------------------------------------------- */

  status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* FFD box must have 6 faces */
  if (nface != 6) {
    status = EGADS_NOTFOUND;
    goto cleanup;
  }

  for (i = 0; i < nface; i++) {

    status = EG_getTopology(efaces[i], &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
    CHECK_STATUS(EG_getTopology);

    faceCandidates = new ego[nface];

    /* fill in the faceCandidates array used in finding opposite pairs of faces */
    for (j = 0; j < nface; j++) {
      faceCandidates[j] = NULL;
      if (i != j) faceCandidates[j] = efaces[j];
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

    status = EG_attributeRet(efaces[i], "ffd_uMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      uMin_attr = true;

      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_uMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      degree = (int)pdouble[0];

      if (degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Get the number of u-control points */
      if (alen == 2)
        ncp = (int)pdouble[1];
      else
      {
        ncp = alen-1;
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
      }

      /* check degree and number of control points relationship */
      if (degree  + 1 > ncp ) {
        snprintf(message, LENMSG, "Number of control points is too few for degree of spline in uMin direction.\n");
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

    status = EG_attributeRet(efaces[i], "ffd_vMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      vMin_attr = true;

      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_vMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      degree = (int)pdouble[0];

      if (degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Get the number of v-control points */
      if (alen == 2)
        ncp = (int)pdouble[1];
      else
      {
        ncp = alen-1;
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
      }

      /* check degree and number of control points relationship */
      if (degree  + 1 > ncp ) {
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

    status = EG_attributeRet(efaces[i], "ffd_wMin", &atype, &alen, NULL, &pdouble, NULL);
    if (status == 0) {
      wMin_attr = true;
      if (alen < 2) {
        snprintf(message, LENMSG, "ffd_wMin attribute must be at least length two, but has length %d. \n", alen);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      degree = (int)pdouble[0];

      if (degree < 1) {
        snprintf(message, LENMSG, "FFD Degree < 1 \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      if (degree > SPLINE3D_MAXDEG) {
        snprintf(message, LENMSG, "FFD Degree > %d \n", SPLINE3D_MAXDEG);
        status  = EGADS_RANGERR;
        goto cleanup;
      }

      /* Get the number of w-control points */
      if (alen == 2)
        ncp = (int)pdouble[1];
      else
      {
        ncp = alen-1;
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
      }

      /* check degree and number of control points relationship */
      if (degree  + 1 > ncp ) {
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

    status = EGADS_SUCCESS;

  } else {
    //snprintf(message, LENMSG, "FFD box must have one each of ffd_uMin, ffd_vMin, and ffd_wMin attributes on faces. \n");
    status  = EGADS_NOTFOUND;
    //goto cleanup;
  }


cleanup:

  FREE(efaces);

  return status;
}

struct KnotSpan_T
{
  KnotSpan_T() : index(-1), u(0), du(0), nk(0) {}

  double span() const { return du/(nk+1); }

  bool operator<(const KnotSpan_T& du) { return span() < du.span(); }
  bool operator>(const KnotSpan_T& du) { return span() > du.span(); }

  int index;
  double u;  // lower knot value
  double du; // span in the knots
  int nk;    // number of knots in the span
};

static void
bubbleSortSpan(std::vector<KnotSpan_T>& spans)
{
  std::size_t i, j;
  for (i = 0; i < spans.size()-1; i++)
    // Last i elements are already in place
    for (j = 0; j < spans.size()-i-1; j++)
      if (spans[j] > spans[j+1])
        std::swap(spans[j], spans[j+1]);
}

static void
bubbleSortIndex(std::vector<KnotSpan_T>& spans)
{
  std::size_t i, j;
  for (i = 0; i < spans.size()-1; i++)
    // Last i elements are already in place
    for (j = 0; j < spans.size()-i-1; j++)
      if (spans[j].index > spans[j+1].index)
        std::swap(spans[j], spans[j+1]);
}

static int
newKnots(udp_T *udp, const int *ivec, const double *rvec, const double *range,
         int& nU, double **Us)
{
  int status = EGADS_SUCCESS;

  ROUTINE(newKnots);

  /* --------------------------------------------------------------- */

  nU = 0;
  *Us = NULL;

//int ideg  = ivec[1];
  int iknot = ivec[3];

  int imin = 0;
  int imax = iknot-1;

  // Count the number of control points within the range
  // and find the bounds
  for (int i = 0; i < iknot; i++) {
    if (rvec[i] <= range[0]) {
      imin = i;
    } else if (rvec[i] >= range[1]) {
      imax = i;
      break;
    }
  }

  const int nk = imax-imin;
  // Nothing to do if count sufficient
  if (nk >= MINCP(udp)) return EGADS_SUCCESS;

  std::vector<KnotSpan_T> spans(nk);

  for (int i = imin; i < imax; i++) {
    double ul = i == imin   ? range[0] : rvec[i];
    double ur = i == imax-1 ? range[1] : rvec[i+1];
    spans[i-imin].u  = ul;
    spans[i-imin].du = ur - ul;
    spans[i-imin].index = i-imin;
  }

  bubbleSortSpan(spans);

  int nadd = MINCP(udp) - nk;
  int nadded = 0;

  if (nk == 1) {
    spans.back().nk = nadd;
  } else {
    while (nadded < nadd) {
      // add a split to the largest span
      spans.back().nk++;
      nadded++;

      // bubble down the now shortest span
      int j = spans.size()-2;
      while (spans[j+1] < spans[j]) {
        std::swap(spans[j+1], spans[j]);
        j--;
        if (j < 0) break;
      }
    }

    bubbleSortIndex(spans);
  }

  const double KNACC = 2e-12;

  int kk = 0;
  nU = nadd;
  if (rvec[0]       < range[0]-KNACC) { nU+=1; kk+=1;}
  if (rvec[iknot-1] > range[1]+KNACC) { nU+=1; }

  MALLOC(*Us, double, nU);

  for (std::size_t i = 0; i < spans.size(); i++) {
    for (int k = 0; k < spans[i].nk; k++) {
      (*Us)[kk++] = spans[i].u + spans[i].du*(k+1)/(spans[i].nk+1);
    }
  }

  if (rvec[0]       < range[0]-KNACC) { for (int k = 0; k < 1; k++) (*Us)[k]      = range[0]; }
  if (rvec[iknot-1] > range[1]+KNACC) { for (int k = 0; k < 1; k++) (*Us)[nU-1-k] = range[1]; }

//  std::cout << "--> " << range[0] << " " << range[1] << std::endl;
//  for (int i = 0; i < nU; i++)
//    std::cout << i << " " << (*Us)[i] << std::endl;

cleanup:
  return status;
}

/*
deform_body is body to deform
uvArrayData vectors should be updated with appropriate ivec and rvec so new geometry can be made
*/
static int
convertDeformToBSpline(
  udp_T *udp,
  ego context,
  ego deform_body,
  int n_inter_edge,
  convert_T convert_edges[],
  int n_inter_face,
  convert_T convert_faces[],
  int rface_len,          /* (in)  */
  ego *bspline_body,      /* (out) body with parts converted to b-splines */
  char *message)          /* (out) error message */
{
  int    status = EGADS_SUCCESS;

  int    oclass, mtype, *senses, nnode, nloop, nloop_edge, edge_mtype;
  int    i = 0, nU, nV, maxDeg=0, iloop, iedge, inode, atype, alen;
  ego    ecurve, esurf, *efaces = NULL, *eloops = NULL, *eloop_edges = NULL, *eedges = NULL, *enodes = NULL;
  ego    ref;
  ego    new_surface, surfaceOut, new_curve, curveOut, epcurve=NULL;
  int    nNewNode, nNewCurv, nNewSurf;
  ego    *eNewNodes=NULL, *eNewCurvs=NULL, *eNewSurfs=NULL;
  double data[18], tdata[2];
  int    *ivec = NULL;
  double *rvec = NULL;
  double *Uknots=NULL, *Vknots=NULL;
  const char *pstr_edge=NULL, *pstr_face=NULL;

  ROUTINE(convertDeformToBSpline);

  /* --------------------------------------------------------------- */

  maxDeg = MINDEG(udp);

#if 0
  status = EG_getBodyTopos(deform_body, NULL, NODE, &nnode, NULL);
  CHECK_STATUS(EG_getBodyTopos);

  /* special processing for NodeBody */
  if (nnode == 1) {
      MALLOC(eloops, ego, 1);

      tdata[0] = 0;
      tdata[1] = 1;
      sense     = SFORWARD;

      status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, tdata, 1, enodes, &sense, &eedges[0]);
      CHECK_STATUS(EG_makeTopology);

      status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 1, eedges, &sense, eloops);
      CHECK_STATUS(EG_makeTopology);

      status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, eloops, NULL, bspline_body);
      CHECK_STATUS(EG_makeTopology);

      FREE(eloops);

      goto cleanup;
  }
#endif


  status = EG_getBodyTopos(deform_body, NULL, NODE, &nNewNode, &eNewNodes);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(deform_body, NULL, EDGE, &nNewCurv, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(deform_body, NULL, FACE, &nNewSurf, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  /* get the current Curves */
  MALLOC(eNewCurvs, ego, nNewCurv);
  for (i = 0; i < nNewCurv; i++) {
    status = EG_getTopology(eedges[i], &eNewCurvs[i], &oclass, &mtype, tdata, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);
  }

  /* Convert curves to B-Splines if necessary, and replace them in the eNewCurvs array. */
  for (i = 0; i < n_inter_edge; i++) {

    /* Only look at edges we need to replace. */
    if (convert_edges[i].index != -1) {
      /* skip degenerate edges */
      if (eedges[convert_edges[i].index]->mtype == DEGENERATE) {
       continue;
      }

      status = EG_getTopology(eedges[convert_edges[i].index], &ecurve, &oclass, &edge_mtype, tdata, &nnode, &enodes, &senses);
      CHECK_STATUS(EG_getTopology);

      if (ecurve->mtype != BSPLINE) {
        status = EG_convertToBSpline(eedges[convert_edges[i].index], &new_curve);
        CHECK_STATUS(EG_convertToBSpline);
      } else {
        new_curve = ecurve;
      }

      /* Check the degree */
      status = EG_getGeometry(new_curve, &oclass, &mtype, &ref, &ivec, &rvec);
      CHECK_STATUS(EG_getGeometry);

      maxDeg = MAX(maxDeg, ivec[1]);
      FREE(ivec);
      FREE(rvec);

      /* No new curve, so continue to next */
      if (new_curve == ecurve) continue;

      eNewCurvs[convert_edges[i].index] = new_curve;
    }
  }

  if (nNewSurf > 0) {

    /* get the current Surfaces */
    MALLOC(eNewSurfs, ego, nNewSurf);
    for (i = 0; i < nNewSurf; i++) {
      status = EG_getTopology(efaces[i], &eNewSurfs[i], &oclass, &mtype, data, &nloop, &eloops, &senses);
      CHECK_STATUS(EG_getTopology);
    }
    for (i = 0; i < n_inter_face; i++) {

      if (convert_faces[i].index != -1) {

        status = EG_getTopology(efaces[convert_faces[i].index], &esurf, &oclass, &mtype, data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        if (esurf->mtype != BSPLINE) {
          status = EG_convertToBSpline(efaces[convert_faces[i].index], &new_surface);
          CHECK_STATUS(EG_convertToBSpline);
        } else {
          new_surface = esurf;
        }

        /* Check the degree */
        status = EG_getGeometry(new_surface, &oclass, &mtype, &ref, &ivec, &rvec);
        CHECK_STATUS(EG_getGeometry);

        maxDeg = MAX(maxDeg, ivec[1]);
        maxDeg = MAX(maxDeg, ivec[4]);
        FREE(ivec);
        FREE(rvec);

        /* No new surface, so continue to next */
        if (new_surface == esurf) continue;

        /* Store the new surface */
        eNewSurfs[convert_faces[i].index] = new_surface;
      }
    }
  }


  /* Refine b-splines to the maximum degree */
  for (i = 0; i < n_inter_edge; i++) {

    /* Only look at edges we need to replace. */
    if (convert_edges[i].index != -1) {
      /* degenerate edges */
      if (eedges[convert_edges[i].index]->mtype == DEGENERATE) {
       continue;
      }

      status = EG_getTopology(eedges[convert_edges[i].index], &ecurve, &oclass, &edge_mtype, tdata, &nnode, &enodes, &senses);
      CHECK_STATUS(EG_getTopology);

      new_curve = eNewCurvs[convert_edges[i].index];

      /* Adding more knots to edges with too few */
      status = EG_getGeometry(new_curve, &oclass, &mtype, &ref, &ivec, &rvec);
      CHECK_STATUS(EG_getGeometry);

      status = newKnots(udp, ivec, rvec, convert_edges[i].range, nU, &Uknots);
      CHECK_STATUS(newKnots);

      if (nU > 0 || maxDeg > ivec[1]) {
        status = EG_addKnots(new_curve, maxDeg, nU, Uknots, 0, 0, NULL, &curveOut);
        CHECK_STATUS(EG_addKnots);
        FREE(Uknots);

        if (new_curve != ecurve) {
          status = EG_deleteObject(new_curve);
          CHECK_STATUS(EG_deleteObject);
        }

        new_curve = curveOut;
      }
      FREE(ivec);
      FREE(rvec);

      eNewCurvs[convert_edges[i].index] = new_curve;
    }
  }

  if (nNewSurf > 0) {

    // Add knots and increase the degree of the B-splines
    for (i = 0; i < n_inter_face; i++) {

      if (convert_faces[i].index != -1) {

        status = EG_getTopology(efaces[convert_faces[i].index], &esurf, &oclass, &mtype, data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        new_surface = eNewSurfs[convert_faces[i].index];

        // Check the knot count and add more if necessary
        status = EG_getGeometry(new_surface, &oclass, &mtype, &ref, &ivec, &rvec);
        CHECK_STATUS(EG_getGeometry);

        // u-knots
        status = newKnots(udp, ivec, rvec, convert_faces[i].range, nU, &Uknots);
        CHECK_STATUS(newKnots);

        // v-knots
        status = newKnots(udp, ivec+3, rvec+ivec[3], convert_faces[i].range+2, nV, &Vknots);
        CHECK_STATUS(newKnots);

        if (nU > 0 || nV > 0 || maxDeg > ivec[1] || maxDeg > ivec[4]) {
          status = EG_addKnots(new_surface, maxDeg, nU, Uknots, maxDeg, nV, Vknots, &surfaceOut);
          CHECK_STATUS(EG_addKnots);

          if (new_surface != esurf) {
            status = EG_deleteObject(new_surface);
            CHECK_STATUS(EG_deleteObject);
          }

          //status = EG_getGeometry(surfaceOut, &oclass, &mtype, &ref, &ivec, &rvec);
          //CHECK_STATUS(EG_getGeometry);

          new_surface = surfaceOut;
        }

        FREE(ivec);
        FREE(rvec);
        FREE(Uknots);
        FREE(Vknots);

        /* Store the new surface */
        eNewSurfs[convert_faces[i].index] = new_surface;
      }
    }
  }

#if 0
  /* make bsurface from the ivec and rvec in builder */
  for (i = 0; i < n_inter_face; i++) {

    if (convert_faces[i].index == -1) continue;

    status = EG_getTopology(efaces[convert_faces[i].index], &eref, &oclass, &mtype, data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);

    for (iloop = 0; iloop < nloop; iloop++) {

      status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype, data, &nloop_edge, &eloop_edges, &senses);
      CHECK_STATUS(EG_getTopology);

      for (iedge = 0; iedge < nloop_edge; iedge++) {

        status = EG_attributeRet(eloop_edges[iedge], "ffdID", &atype, &alen, NULL, NULL, NULL);
        if (status == EGADS_NOTFOUND) continue;

        index = EG_indexBodyTopo(deform_body, eloop_edges[iedge]);
        status = EG_tolerance(eedges[iedge], &tol);
        CHECK_STATUS(EG_tolerance);

        if (eref == NULL) {
          status = EG_otherCurve(eNewSurfs[convert_faces[i].index], eloop_edges[iedge], 0., &epcurve);
          CHECK_STATUS(EG_otherCurve);

        } else {
          epcurve = eloop_edges[nloop_edge+iedge];
        }

        status = EG_getTopology(eloop_edges[iedge], &ecurve, &oclass, &mtype, data, &nnode, &enodes, &senses);
        CHECK_STATUS(EG_getTopology);

        if (epcurve->mtype != TRIMMED &&
            epcurve->mtype != BEZIER &&
            epcurve->mtype != BSPLINE) {
          status = EG_makeGeometry(context, PCURVE, TRIMMED, epcurve, NULL, data, &etmp);
          CHECK_STATUS(EG_getTopology);
          epcurve = etmp;
        }

        status = EG_otherCurve(eNewSurfs[convert_faces[i].index], epcurve, tol, &eNewCurvs[index-1]);
        CHECK_STATUS(EG_otherCurve);
      }
    }
  }
#endif

  /* build the body */
  status = buildBody(deform_body,
                     nNewNode, eNewNodes,
                     nNewCurv, eNewCurvs,
                     nNewSurf, eNewSurfs,
                     bspline_body,
                     LENMSG, message);
  CHECK_STATUS(buildBody);

  /* Check Face boundary conditions for parametric iso-lines */
  FREE(efaces);
  status = EG_getBodyTopos(*bspline_body, NULL, FACE, &nNewSurf, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  for (i = 0; i < n_inter_face; i++) {

    if (convert_faces[i].index == -1) continue;

    status = EG_getTopology(efaces[convert_faces[i].index], &esurf, &oclass, &mtype, data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_attributeRet(efaces[convert_faces[i].index], "ffdFixed", &atype, &alen, NULL, NULL, &pstr_face);
    if ((pstr_face != NULL) &&
        (strcasecmp(pstr_face,"G0") != 0 &&
         strcasecmp(pstr_face,"G1") != 0)) {
      snprintf(message, LENMSG, "'ffdFixed' Face attribute must be a string 'G0' or 'G1'\n");
      status = EGADS_ATTRERR;
      goto cleanup;
    }
    if (status == EGADS_NOTFOUND) continue;

    for (iloop = 0; iloop < nloop; iloop++) {

      status = EG_getTopology(eloops[iloop], &esurf, &oclass, &mtype, data, &nloop_edge, &eloop_edges, &senses);
      CHECK_STATUS(EG_getTopology);

      for (iedge = 0; iedge < nloop_edge; iedge++) {

        status = EG_attributeRet(eloop_edges[iedge], "ffdFixed", &atype, &alen, NULL, NULL, &pstr_edge);
        if (status == EGADS_NOTFOUND) {
          // Add the face attribute if not on the edge
          if (pstr_face != NULL) {
            status = EG_attributeAdd(eloop_edges[iedge], "ffdFixed", ATTRSTRING, 0, NULL, NULL, pstr_face);
            CHECK_STATUS(EG_attributeAdd);

            // Also add to the Nodes of the Edge
            status = EG_getTopology(eloop_edges[iedge], &ecurve, &oclass, &mtype, data, &nnode, &enodes, &senses);
            CHECK_STATUS(EG_getTopology);
            for (inode = 0; inode < nnode; inode++) {
              if (EGADS_NOTFOUND == EG_attributeRet(enodes[inode], "ffdFixed", &atype, &alen, NULL, NULL, NULL) ) {
                status = EG_attributeAdd(enodes[inode], "ffdFixed", ATTRSTRING, 0, NULL, NULL, pstr_face);
                CHECK_STATUS(EG_attributeAdd);
              }
            }
          }
          continue;
        }

        if ((atype != ATTRSTRING) ||
            (pstr_edge == NULL) ||
            (strcasecmp(pstr_edge,"G0") != 0 &&
             strcasecmp(pstr_edge,"G1") != 0)) {
          snprintf(message, LENMSG, "'ffdFixed' Edge attribute must be a string 'G0' or 'G1'\n");
          status = EGADS_ATTRERR;
          goto cleanup;
        }

        epcurve = eloop_edges[nloop_edge+iedge];

        if (EG_isIsoPCurve(epcurve,NULL,NULL,NULL) != EGADS_SUCCESS) {
          snprintf(message, LENMSG, "'ffdFixed' Edge '%d' must be an IsoCline on Face %d\n",
                   EG_indexBodyTopo(deform_body, eloop_edges[iedge]), convert_faces[i].index+1);
          status = EGADS_TOPOERR;
          goto cleanup;
        }
      }
    }
  }

  /* Check Face boundary conditions for parametric iso-lines */
  FREE(eedges);
  status = EG_getBodyTopos(*bspline_body, NULL, EDGE, &nNewCurv, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  /* Refine b-splines to the maximum degree */
  for (iedge = 0; iedge < nNewCurv; iedge++) {

    /* degenerate edges */
    if (eedges[iedge]->mtype == DEGENERATE) continue;

    status = EG_getTopology(eedges[iedge], &ecurve, &oclass, &edge_mtype, tdata, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_attributeRet(eedges[iedge], "ffdFixed", &atype, &alen, NULL, NULL, &pstr_edge);
    if ((pstr_edge != NULL) &&
        (strcasecmp(pstr_edge,"G0") != 0 &&
         strcasecmp(pstr_edge,"G1") != 0)) {
      snprintf(message, LENMSG, "'ffdFixed' Edge %d attribute must be a string 'G0' or 'G1'\n",
                                EG_indexBodyTopo(deform_body, eedges[convert_edges[i].index]));
      status = EGADS_ATTRERR;
      goto cleanup;
    }
    if (status == EGADS_NOTFOUND) continue;

    // Also add to the Nodes of the Edge
    status = EG_getTopology(eedges[iedge], &ecurve, &oclass, &mtype, data, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);
    for (inode = 0; inode < nnode; inode++) {
      if (EGADS_NOTFOUND == EG_attributeRet(enodes[inode], "ffdFixed", &atype, &alen, NULL, NULL, NULL) ) {
        status = EG_attributeAdd(enodes[inode], "ffdFixed", ATTRSTRING, 0, NULL, NULL, pstr_edge);
        CHECK_STATUS(EG_attributeAdd);
      }
    }
  }


  status = EGADS_SUCCESS;

cleanup:

  FREE(eNewNodes);
  FREE(eNewCurvs);
  FREE(eNewSurfs);

  FREE(eedges);
  FREE(efaces);

  FREE(Uknots);
  FREE(Vknots);

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
udpExecute(ego emodel,                  /* (in)  input model */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
  int    status = EGADS_SUCCESS;
  ego    context=NULL;
  ego    eref, *efaces = NULL, *eedges = NULL, *enodes = NULL, result[2]={NULL,NULL};
  ego    bspline_body, *ebodies = NULL, deform_body = NULL, ffd_body = NULL;
  ego    intersect_model=NULL, *edge_faces=NULL, *node_edges=NULL;
  ego    *inter_ebodies = NULL, *inter_faces = NULL, *inter_edges = NULL, *inter_nodes = NULL;
  int    n_inter_body, n_inter_face, n_inter_edge, n_inter_node;
  int    oclass, mtype, nface, nedge, nnode, periodic, *senses;
  int    i, nbody, atype, alen, index; //, n_node_edge; //n_edge_face
  int    identified_boxes = 0, rface_len = 0, nbspl_topo=0;
  double data[18];
  char   *message=NULL;
  const int     *pint;
  const double  *dint;
  //bool  fixed;

  convert_T *convert_faces = NULL, *convert_edges = NULL;

  udp_T   *udps = NULL;
  udp_T   *udp = NULL;

  ROUTINE(udpExecute);

  /* --------------------------------------------------------------- */

#ifdef DEBUG
  printf("udpExecute(emodel=%llx)\n", (long long)emodel);
#endif

  /* default return values */
  *ebody  = NULL;
  *nMesh  = 0;
  *string = NULL;

  MALLOC(message, char, LENMSG);
  message[0] = '\0';

  status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                          data, &nbody, &ebodies, &senses);
  CHECK_STATUS(EG_getTopology);

  if (oclass != MODEL) {
    snprintf(message, LENMSG, "Expecting a Model\n");
    status = EGADS_NOTMODEL;
    goto cleanup;
  } else if (nbody != 2) {
    snprintf(message, LENMSG, "Expecting Model to contain exactly two bodies, but %d bodies found.\n", nbody);
    status = EGADS_NOTBODY;
    goto cleanup;
  }

  /* There should be one body to deform (deform body), and one ffd box body. */
  for (i = 0; i < nbody; i++) {
    status = checkFFDbox(ebodies[i], message);
    if (status == EGADS_SUCCESS)
    {
      if (identified_boxes > 1) {
        snprintf(message, LENMSG, "Too many FFD boxes identified. \n");
        status  = EGADS_RANGERR;
        goto cleanup;
      } else {
        /* copy the FFD box so it can be returned */
        status = EG_copyObject(ebodies[i], NULL, &ffd_body);
        CHECK_STATUS(EG_copyObject);
        identified_boxes++;
      }
    }
    else if (status == EGADS_RANGERR)
    {
      goto cleanup;
    }
    else
    { /* The body is not an FFD body */
      if (deform_body != NULL) {
        snprintf(message, LENMSG, "More than one deform body identified. \n"
            "FFD box must have one each of ffd_uMin, ffd_vMin, and ffd_wMin attributes on faces. \n");
        status  = EGADS_NOTFOUND;
        goto cleanup;
      } else {
        /* copy the deformed body to place attributes on it */
        status = EG_copyObject(ebodies[i], NULL, &deform_body);
        CHECK_STATUS(EG_copyObject);
      }
    }
  }

  /* Cache copy of arguments for future use.
   * This also increments numUdp and copies udps[0] to udps[numUdp].
   * Caching should only be performed after checking for valid inputs.
   */
  status = cacheUdp(emodel);
  CHECK_STATUS(cacheUdp);
  udp = &udps[numUdp];

  /* Look at input bodies */
  status = EG_getContext(emodel, &context);
  CHECK_STATUS(cacheUdp);

  /* Attribute the deform body before intersection */
  status = EG_getBodyTopos(deform_body, NULL, FACE, &nface, &efaces);
  CHECK_STATUS(EG_getBodyTopos);

  /* attribute DEFORMBODY faces */
  for (i = 0; i < nface; i++) {
    status = EG_attributeAdd(efaces[i], DEFORMBODY, ATTRINT, 1, &i, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);
  }

  status = EG_getBodyTopos(deform_body, NULL, EDGE, &nedge, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  /* attribute DEFORMBODY edges */
  for (i = 0; i < nedge; i++) {
#if 0
    // delete attributes --> issue with replace faces and _body attributes for edges and nodes
    status = EG_attributeDel(eedges[i], NULL);
    CHECK_STATUS(EG_attributeDel);
#endif

    status = EG_attributeAdd(eedges[i], DEFORMBODY, ATTRINT, 1, &i, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);
  }

  status = EG_getBodyTopos(deform_body, NULL, NODE, &nnode, &enodes);
  CHECK_STATUS(EG_getBodyTopos);

  /* attribute DEFORMBODY nodes */
  for (i = 0; i < nnode; i++) {
#if 0
    // delete attributes --> issue with replace faces and _body attributes for edges and nodes
    status = EG_attributeDel(enodes[i], NULL);
    CHECK_STATUS(EG_attributeDel);
#endif

    status = EG_attributeAdd(enodes[i], DEFORMBODY, ATTRINT, 1, &i, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);
  }

  /* Intersect the FFD body and the deform body with full attribution enabled */
  status = EG_setFullAttrs(context, 1);
  CHECK_STATUS(EG_setFullAttrs);

  status = EG_generalBoolean(deform_body, ffd_body, INTERSECTION, 0, &intersect_model);
  CHECK_STATUS(EG_generalBoolean);

  status = EG_setFullAttrs(context, 0);
  CHECK_STATUS(EG_setFullAttrs);

//  remove("intersection.egads");
//  EG_saveModel(intersect_model, "intersection.egads");

  status = EG_getTopology(intersect_model, &eref, &oclass, &mtype, data, &n_inter_body, &inter_ebodies, &senses);
  CHECK_STATUS(EG_getTopology);

  if (n_inter_body == 0) {
    snprintf(message, LENMSG, "FFD box does not intersect the deform body. \n");
    status  = EGADS_RANGERR;
    goto cleanup;
  }

  if (n_inter_body > 1) {
    snprintf(message, LENMSG, "Expected one intersection body from FFD box and deform body, but intersection resulted in %d bodies. \n", n_inter_body);
    status  = EGADS_RANGERR;
    goto cleanup;
  }

  status = EG_getBodyTopos(inter_ebodies[0], NULL, FACE, &n_inter_face, &inter_faces);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(inter_ebodies[0], NULL, EDGE, &n_inter_edge, &inter_edges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(inter_ebodies[0], NULL, NODE, &n_inter_node, &inter_nodes);
  CHECK_STATUS(EG_getBodyTopos);


  /* Maps faces to convert to B-Splines to the index of efaces for the deform body
   * array has value -1 if the face or edge is not in intersection and does not need
   * to be converted to b-spline
   */
  MALLOC(convert_faces, convert_T, n_inter_face);
  MALLOC(convert_edges, convert_T, n_inter_edge);

  /* identify the faces to be deformed and attribute DEFORMBODY */
  for (i = 0; i < n_inter_face; i++) {
    convert_faces[i].index = -1;
    status = EG_attributeRet(inter_faces[i], DEFORMBODY, &atype, &alen, &pint, &dint, NULL);
    if (status == EGADS_SUCCESS) {
      index = pint[0];
      //if (EG_isEquivalent(inter_faces[i], efaces[index]) != EGADS_SUCCESS) continue;
      convert_faces[i].index = index;
      status = EG_getRange(inter_faces[i], convert_faces[i].range, &periodic);
      CHECK_STATUS(EG_getRange);
      status = EG_attributeAdd(efaces[index], FFDID, ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);
      rface_len++;
      nbspl_topo++;
    }
  }

  /* identify the edges to be deformed and attribute DEFORMBODY */
  for (i = 0; i < n_inter_edge; i++) {
    convert_edges[i].index = -1;
    status = EG_attributeRet(inter_edges[i], DEFORMBODY, &atype, &alen, &pint, &dint, NULL);
    if (status == EGADS_SUCCESS && inter_edges[i]->mtype != DEGENERATE) {
      index = pint[0];
      /*
      if (EG_isEquivalent(inter_edges[i], eedges[index]) != EGADS_SUCCESS) continue;

      status = EG_getBodyTopos(deform_body, eedges[index], FACE, &n_edge_face, &edge_faces);
      CHECK_STATUS(EG_getBodyTopos);

      fixed = false;
      for (j = 0; j < n_edge_face; j++) {
        status = EG_attributeRet(edge_faces[j], FFDID, &atype, &alen, &pint, NULL, NULL);
        if (status == EGADS_NOTFOUND){
          fixed = true;
          break;
        }
      }
      FREE(edge_faces);
      if (fixed) continue;
      */
      convert_edges[i].index = index;
      status = EG_getRange(inter_edges[i], convert_edges[i].range, &periodic);
      CHECK_STATUS(EG_getRange);
      status = EG_attributeAdd(eedges[index], FFDID, ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);
      nbspl_topo++;
    }
  }

  /* identify the nodes to be deformed and attribute DEFORMBODY */
  for (i = 0; i < n_inter_node; i++) {
    status = EG_attributeRet(inter_nodes[i], DEFORMBODY, &atype, &alen, &pint, &dint, NULL);
    if (status == EGADS_SUCCESS) {
      index = *pint;
/*
      status = EG_getBodyTopos(deform_body, enodes[index], EDGE, &n_node_edge, &node_edges);
      CHECK_STATUS(EG_getBodyTopos);

      fixed = false;
      for (j = 0; j < n_node_edge; j++) {
        status = EG_attributeRet(node_edges[j], FFDID, &atype, &alen, &pint, NULL, NULL);
        if (status == EGADS_NOTFOUND){
          fixed = true;
          break;
        }
      }
      FREE(node_edges);
      if (fixed) continue;
*/
      status = EG_attributeAdd(enodes[index], FFDID, ATTRINT, 1, &i, NULL, NULL);
      CHECK_STATUS(EG_attributeAdd);
      nbspl_topo++;
    }
  }

  if (nbspl_topo == 0) {
    snprintf(message, LENMSG, "No complete FACE/EDGE/NODE within the FFD box. \n");
    status  = EGADS_RANGERR;
    goto cleanup;
  }

  status = convertDeformToBSpline(udp, context, deform_body, n_inter_edge, convert_edges, n_inter_face,
                                  convert_faces, rface_len, &bspline_body, message);
  CHECK_STATUS(convertDeformToBSpline);

#if 0
  status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
  CHECK_STATUS(EG_getBodyTopos);

  for (i = 0; i < nedge; i++) {
    status = EG_attributeRet(eedges[i], FFDID, &atype, &alen, &pint, NULL, NULL);
    CHECK_STATUS(EG_attributeRet);
  }
#endif

  /* return FFD body and b-spline body for FFDdeform */
  result[0] = ffd_body;
  result[1] = bspline_body;
  status = EG_makeTopology(context, NULL, MODEL, 0, NULL, 2, result, NULL, ebody);
  CHECK_STATUS(EG_makeTopology);

  /* remember this model (body) */
  udps[numUdp].ebody = *ebody;
  status = EGADS_SUCCESS;

#ifdef DEBUG
  printf("udpExecute -> *ebody=%lx\n", (long)(*ebody));
#endif

cleanup:
  if (strlen(message) > 0) {
    *string = message;
    printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
    printf("ERROR in udpExecute: status = %d \n", status);
    FREE(message);
    *string = udpErrorStr(status);
  } else {
    FREE(message);
  }

  /* remove all the construction parts first */
  int outLevel = EG_setOutLevel(context, 0);
  EG_deleteObject(context);
  EG_setOutLevel(context, outLevel);

  /* delete temporary bodies */
  EG_deleteObject(deform_body);
  EG_deleteObject(intersect_model);

  FREE(enodes);
  FREE(eedges);
  FREE(efaces);

  FREE(edge_faces);
  FREE(node_edges);

  FREE(inter_nodes);
  FREE(inter_edges);
  FREE(inter_faces);

  FREE(convert_faces);
  FREE(convert_edges);

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
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int iudp, judp;

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

    /* this routine is not written yet */
    return EGADS_NOLOAD;
}
