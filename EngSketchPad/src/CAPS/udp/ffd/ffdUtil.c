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

#include "egads.h"
#include "egads_dot.h"
#include "common.h"
#include "OpenCSM.h"

/* Re-build a body given a set of possibly new geometric nodes/curves/surfaces.
 */

int buildBody(ego eOrigBody,                       /* (in) original body to rebuild */
 /*@unused@*/ const int nNewNode, ego *eNewNodes,  /* (in) mixture of orig and new Nodes    (Node BodyIndex ordered) */
              const int nNewCurv, ego *eNewCurvs,  /* (in) mixture of orig and new Curves   (Edge BodyIndex ordered) */
              const int nNewSurf, ego *eNewSurfs,  /* (in) mixture of orig and new Surfaces (Face BodyIndex ordered) */
              ego *eNewBody,                       /* (out) the new body with all geometry replaced */
              int lenmsg,                          /* (in) error message length */
              char *message)                       /* (out) error message */
{
  int             status = EGADS_SUCCESS;

  int i, replace, replace_loop, replace_face;
  int index, fullattr;
  int nnode, nedge, nloop;
  int iedge, iloop, iface;
  int oclass, mtype, edge_mtype, face_mtype, *senses, *loop_senses, *edge_senses;
  int nOrigFace, nOrigEdge, nOrigNode;
  int nReplaceFace, hasdot;
  ego context;
  ego *eOrigFaces=NULL, *eOrigLoops=NULL, *eOrigEdges=NULL, *eOrigNodes=NULL;
  ego eNewFace, *eNewLoops=NULL, *eNewEdges=NULL;
  ego eref, face_surf, *enodes, *loop_eedges, *eNewLoopEdges=NULL, edge_nodes[2];
  ego *replaceFaces=NULL;
  double edge_data[4], data[4], zero[4] = {0,0,0,0};
  double tolerance;

  ROUTINE(buildBody);

  /* --------------------------------------------------------------- */

  status = EG_getContext(eOrigBody, &context);
  CHECK_STATUS(EG_getBodyTopos);

  //status = EG_tolerance(eOrigBody, &tolerance);
  //CHECK_STATUS(EG_tolerance);
  // printf("Body tolerance before conversion to B-Spline = %le \n", tolerance);

  hasdot = EG_hasGeometry_dot(eOrigBody);

  status = EG_getBodyTopos(eOrigBody, NULL, NODE, &nOrigNode, &eOrigNodes);
  CHECK_STATUS(EG_getBodyTopos);
  SPLINT_CHECK_FOR_NULL(eOrigNodes);

  status = EG_getBodyTopos(eOrigBody, NULL, EDGE, &nOrigEdge, &eOrigEdges);
  CHECK_STATUS(EG_getBodyTopos);

  status = EG_getBodyTopos(eOrigBody, NULL, FACE, &nOrigFace, &eOrigFaces);
  CHECK_STATUS(EG_getBodyTopos);

  if (nOrigFace != nNewSurf) {
    snprintf(message, lenmsg, "Inconsistent Surface and Face counts\n");
    status = EGADS_RANGERR;
    goto cleanup;
  }

  if (nOrigEdge != nNewCurv) {
    snprintf(message, lenmsg, "Inconsistent Curve and Edge counts\n");
    status = EGADS_RANGERR;
    goto cleanup;
  }

  MALLOC(eNewEdges, ego, nOrigEdge);
  for (iedge = 0; iedge < nOrigEdge; iedge++) {
    SPLINT_CHECK_FOR_NULL(eOrigEdges);
    eNewEdges[iedge] = eOrigEdges[iedge];
  }

  for (iedge = 0; iedge < nOrigEdge; iedge++) {
    SPLINT_CHECK_FOR_NULL(eOrigEdges);
    status = EG_getTopology(eOrigEdges[iedge], &eref, &oclass, &edge_mtype, edge_data, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);

    /* Check if Node or Curve geometry is replaced */
    replace = 0;
    for (i = 0; i < nnode; i++) {
      edge_nodes[i] = enodes[i];
      index = EG_indexBodyTopo(eOrigBody, enodes[i]);
      if (eNewNodes[index-1] != eOrigNodes[index-1]) {
        edge_nodes[i] = eNewNodes[index-1];
        replace++;
      }
    }

    if (eref != eNewCurvs[iedge]) {
      eref = eNewCurvs[iedge];
      replace++;
    }

    /* Use original Edge if geometry is not changed */
    if (replace == 0) continue;

    /* Make the new Edge */
    status = EG_makeTopology(context, eref, EDGE, edge_mtype, edge_data, nnode, edge_nodes, NULL, &eNewEdges[iedge]);
    CHECK_STATUS(EG_makeTopology);

    if (hasdot == EGADS_SUCCESS) {
      status = EG_setRange_dot(eNewEdges[iedge], EDGE, edge_data, zero);
      CHECK_STATUS(EG_setRange_dot);
    }

    status = EG_tolerance(eNewEdges[iedge], &tolerance);
    CHECK_STATUS(EG_tolerance);

//    printf("New Edge %d tol %10.12e\n", iedge+1, tolerance);

    /* Preserve attributes */
    status = EG_attributeDup(eOrigEdges[iedge], eNewEdges[iedge]);
    CHECK_STATUS(EG_attributeDup);
  }


  // WireBody
  if (nOrigFace == 0) {

    /* Get the loop from the WIREBODY */
    status = EG_getTopology(eOrigBody, &eref, &oclass, &mtype, data, &nloop, &eOrigLoops, &loop_senses);
    CHECK_STATUS(EG_getTopology);
    SPLINT_CHECK_FOR_NULL(eOrigLoops);

    if (oclass != BODY || mtype != WIREBODY) {
      snprintf(message, lenmsg, "This should never happen...\n");
      status = EGADS_TOPOERR;
      goto cleanup;
    }

    MALLOC(eNewLoops, ego, nloop);

    /* Get the Edges from the loop */
    status = EG_getTopology(eOrigLoops[0], &eref, &oclass, &mtype, data, &nedge, &loop_eedges, &edge_senses);
    CHECK_STATUS(EG_getTopology);

    MALLOC(eNewLoopEdges, ego, nedge);

    /* Get the updated Edges */
    for (iedge = 0; iedge < nedge; iedge++) {
      index = EG_indexBodyTopo(eOrigBody, loop_eedges[iedge]);
      eNewLoopEdges[iedge] = eNewEdges[index-1];
    }

    /* Make the new loop */
    status = EG_makeTopology(context, NULL, LOOP, mtype, data, nedge, eNewLoopEdges, edge_senses, &eNewLoops[0]);
    CHECK_STATUS(EG_makeTopology);

    /* Make the new bspline WIREBODY */
    status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eNewLoops[0], loop_senses, eNewBody);
    CHECK_STATUS(EG_makeTopology);

    FREE(eNewLoopEdges);
    FREE(eNewLoops);

    goto cleanup;
  }

  nReplaceFace = 0;
  MALLOC(replaceFaces, ego, 2*nOrigFace);

  /* Recreate faces */
  for (iface = 0; iface < nOrigFace; iface++) {

    replace_face = 0;

    SPLINT_CHECK_FOR_NULL(eOrigFaces);
    status = EG_getTopology(eOrigFaces[iface], &face_surf, &oclass, &face_mtype, data, &nloop, &eOrigLoops, &loop_senses);
    CHECK_STATUS(EG_getTopology);
    SPLINT_CHECK_FOR_NULL(eOrigLoops);

    if (face_surf != eNewSurfs[iface]) {
      replace_face++;
    }

    MALLOC(eNewLoops, ego, nloop);

    for (iloop = 0; iloop < nloop; iloop++) {
      replace_loop = 0;

      eNewLoops[iloop] = eOrigLoops[iloop];

      status = EG_getTopology(eOrigLoops[iloop], &eref, &oclass, &mtype, data, &nedge, &loop_eedges, &edge_senses);
      CHECK_STATUS(EG_getTopology);

      MALLOC(eNewLoopEdges, ego, 2*nedge);

      /* find the edges in this loop in the eNewEdges array */
      for (iedge = 0; iedge < nedge; iedge++) {
        index = EG_indexBodyTopo(eOrigBody, loop_eedges[iedge]);
        eNewLoopEdges[iedge] = eNewEdges[index-1];

        if (loop_eedges[iedge] != eNewEdges[index-1])
          replace_loop++;
      }

      if (face_surf != eNewSurfs[iface]) {
        if (face_surf->mtype == PLANE) {
          /* construct the p-curves for a PLANE */
          for (iedge = 0; iedge < nedge; iedge++) {
            status = EG_otherCurve(eNewSurfs[iface], eNewLoopEdges[iedge], 0., &eNewLoopEdges[nedge+iedge]);
            CHECK_STATUS(EG_otherCurve);
          }
        } else {
          /* recycle the p-curves */
          for (iedge = 0; iedge < nedge; iedge++) {
            eNewLoopEdges[nedge+iedge] = loop_eedges[nedge+iedge];
          }
        }
      }

      if (face_surf != eNewSurfs[iface]) {
        eref = eNewSurfs[iface];
        replace_loop++;
      }

      replace_face += replace_loop;

      if (replace_loop > 0) {
        status = EG_makeTopology(context, eref, LOOP, mtype, data, nedge, eNewLoopEdges, edge_senses, &eNewLoops[iloop]);
        CHECK_STATUS(EG_makeTopology);
      }

      FREE(eNewLoopEdges);
    }

    if (replace_face > 0) {
      status = EG_makeTopology(context, eNewSurfs[iface], FACE, face_mtype, NULL, nloop, eNewLoops, loop_senses, &eNewFace);
      CHECK_STATUS(EG_makeTopology);

      status = EG_attributeDup(eOrigFaces[iface], eNewFace);
      CHECK_STATUS(EG_attributeDup);

      status = EG_tolerance(eNewFace, &tolerance);
      CHECK_STATUS(EG_tolerance);

//      printf("New Face %d tol %10.12e\n", iface+1, tolerance);

      replaceFaces[2*nReplaceFace  ] = eOrigFaces[iface]; // original face
      replaceFaces[2*nReplaceFace+1] = eNewFace; // new face
      nReplaceFace++;
    }

    FREE(eNewLoops);
  }


  status = fullattr = EG_setFullAttrs(context, 1);
  CHECK_STATUS(EG_setFullAttrs);

  if (nOrigFace == 1) {
    status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eNewFace, NULL, eNewBody);
    CHECK_STATUS(EG_makeTopology);
  } else {
    status = EG_replaceFaces(eOrigBody, nReplaceFace, replaceFaces, eNewBody);
    CHECK_STATUS(EG_replaceFaces);
  }

  status = EG_setFullAttrs(context, fullattr);
  CHECK_STATUS(EG_setFullAttrs);

  status = EG_tolerance(*eNewBody, &tolerance);
  CHECK_STATUS(EG_tolerance);

cleanup:
  FREE(eOrigFaces);
  FREE(eOrigEdges);
  FREE(eOrigNodes);

  FREE(eNewLoops);
  FREE(eNewEdges);

  FREE(eNewLoopEdges);
  FREE(replaceFaces);

  return status;
}
