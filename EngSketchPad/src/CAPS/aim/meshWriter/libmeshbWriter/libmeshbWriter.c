/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             libMeshb Mesh Writer Code
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "libmeshbWriter.h"
#include "libMeshb/sources/libmeshb7.h"

/* these values requested to emulate Feflo.a behavior */
#define EXPORT_MESHB_VERTEX_ID (1)
#define EXPORT_MESHB_2D_ID (1)
#define EXPORT_MESHB_3D_ID (0)
#define EXPORT_MESHB_VERTEX_3 (10000000)
#define EXPORT_MESHB_VERTEX_4 (200000000)

const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j, igroup, ielem, nPoint;
  int state, nglobal, id;
  int nNode, nEdge, nFace;
  int iedge, iface;
  int nNodeOffset, nEdgeOffset, nFaceOffset;
  int nNodeVerts, nEdgeVerts, nFaceVerts;
  int nLine, nTri;
  int ptype, pindex, plen, tlen, iglobal;
  int elem[3];
  const double *points, *uv, *t;
  const int *ptypes, *pindexs, *tris, *tric;
  double xyz[3], param[2];
  ego body, *faces=NULL,*edges=NULL;
  char filename[PATH_MAX];
  int xMeshConstant = (int)true, yMeshConstant = (int)true, zMeshConstant = (int)true;
  int64_t fileID=0;
  int meshVersion;
  aimMeshRef *meshRef = NULL;
  aimMeshData *meshData = NULL;

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  meshRef  = mesh->meshRef;
  meshData = mesh->meshData;

  if (meshData->dim != 2 && meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  meshVersion = 2;
  if (EXPORT_MESHB_VERTEX_3 < meshData->nVertex) meshVersion = 3;
  if (EXPORT_MESHB_VERTEX_4 < meshData->nVertex) meshVersion = 4;

  fileID = GmfOpenMesh(filename, GmfWrite, meshVersion, meshData->dim);

  if (fileID == 0) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }

  status = GmfSetKwd(fileID, GmfVertices, meshData->nVertex);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  if (meshData->dim == 2) {

    for (i = 0; i < meshData->nVertex; i++) {
      // Constant x?
      if (fabs(meshData->verts[i][0] - meshData->verts[0][0]) > 1E-7) {
        xMeshConstant = (int) false;
      }

      // Constant y?
      if (fabs(meshData->verts[i][1] - meshData->verts[0][1] ) > 1E-7) {
        yMeshConstant = (int) false;
      }

      // Constant z?
      if (fabs(meshData->verts[i][2] - meshData->verts[0][2]) > 1E-7) {
        zMeshConstant = (int) false;
      }
    }

    if (zMeshConstant == (int) false) {
      printf("libMeshb expects 2D meshes in the x-y plane... attempting to rotate mesh!");

      if (xMeshConstant == (int) true && yMeshConstant == (int) false) {
        printf("Swapping z and x coordinates!\n");
        for (i = 0; i < meshData->nVertex; i++) {
          status = GmfSetLin(fileID, GmfVertices, meshData->verts[i][2],
                                                  meshData->verts[i][1], EXPORT_MESHB_VERTEX_ID);
          if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }

      } else if(xMeshConstant == (int) false && yMeshConstant == (int) true) {

        printf("Swapping z and y coordinates!");
        for (i = 0; i < meshData->nVertex; i++) {
          status = GmfSetLin(fileID, GmfVertices, meshData->verts[i][0],
                                                  meshData->verts[i][2], EXPORT_MESHB_VERTEX_ID);
          if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }

      } else {
        AIM_ERROR(aimInfo, "Unable to rotate mesh!");
        status = CAPS_NOTFOUND;
        goto cleanup;
      }

    } else { // zMeshConstant == true
      // Write nodal coordinates as is
      for (i = 0; i < meshData->nVertex; i++) {
        status = GmfSetLin(fileID, GmfVertices, meshData->verts[i][0],
                                                meshData->verts[i][1], EXPORT_MESHB_VERTEX_ID);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

  } else {

    // Write nodal coordinates
    for (i = 0; i < meshData->nVertex; i++) {
      status = GmfSetLin(fileID, GmfVertices, meshData->verts[i][0],
                                              meshData->verts[i][1],
                                              meshData->verts[i][2], EXPORT_MESHB_VERTEX_ID);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }
  }


  // write out elements

  // count the number of EDGE/FACE elements
  nTri = nLine = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      if (edges[iedge]->mtype == DEGENERATE) continue;
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      AIM_STATUS(aimInfo, status);
      nLine += plen-1;
    }
    AIM_FREE(edges);

    // No triangles for an area mesh
    if (meshRef->type == aimAreaMesh) continue;

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);
      nTri += tlen;
    }
  }

  // Write out EDGE line elements
  status = GmfSetKwd(fileID, GmfEdges, nLine);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nEdgeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      if (status == EGADS_DEGEN) continue;
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < plen-1; j++) {

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 1, &elem[0]);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 2, &elem[1]);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        status = GmfSetLin(fileID, GmfEdges, meshRef->maps[i].map[elem[0]-1],
                                             meshRef->maps[i].map[elem[1]-1],
                                             nEdgeOffset + iedge + 1);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }
    nEdgeOffset += nEdge;
  }

  if (meshRef->type != aimAreaMesh) {
    // Write FACE triangle elements
    status = GmfSetKwd(fileID, GmfTriangles, nTri);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    nFaceOffset = 0;
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
      AIM_STATUS(aimInfo, status);

      for (iface = 0; iface < nFace; iface++) {
        status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        for (j = 0; j < tlen; j++) {
          /* triangle orientation flipped, per refine convention */
          status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 0], &elem[1/*0*/]);
          AIM_STATUS(aimInfo, status);
          status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 1], &elem[0/*1*/]);
          AIM_STATUS(aimInfo, status);
          status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 2], &elem[2]);
          AIM_STATUS(aimInfo, status);

          status = GmfSetLin(fileID, GmfTriangles, meshRef->maps[i].map[elem[0]-1],
                                                   meshRef->maps[i].map[elem[1]-1],
                                                   meshRef->maps[i].map[elem[2]-1],
                                                   nFaceOffset + iface + 1);
          if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }
      }

      nFaceOffset += nFace;
    }
  }

  // Write any remaining elements
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "libMeshb writer currently only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    if (meshData->elemGroups[igroup].elementTopo == aimLine) {
      // written with geometry above
    }

    else if (meshData->elemGroups[igroup].elementTopo == aimTri) {
      if (meshRef->type == aimAreaMesh) {
        status = GmfSetKwd(fileID, GmfTriangles, meshData->elemGroups[igroup].nElems);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

        // element connectivity 1-based
        nPoint = meshData->elemGroups[igroup].nPoint;
        for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
          /* triangle orientation flipped, per refine convention */
          status = GmfSetLin(fileID, GmfTriangles, meshData->elemGroups[igroup].elements[nPoint*ielem+1/*0*/],
                                                   meshData->elemGroups[igroup].elements[nPoint*ielem+0/*1*/],
                                                   meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                                                   igroup+1);
          if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }
      } else {
        // written with geometry above
      }
    }

    else if (meshData->elemGroups[igroup].elementTopo == aimQuad) {
      status = GmfSetKwd(fileID, GmfQuadrilaterals, meshData->elemGroups[igroup].nElems);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      // element connectivity 1-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
        status = GmfSetLin(fileID, GmfQuadrilaterals, meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                                                      meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                                                      meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                                                      meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                                                      igroup+1);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    else if (meshData->elemGroups[igroup].elementTopo == aimTet) {
      status = GmfSetKwd(fileID, GmfTetrahedra, meshData->elemGroups[igroup].nElems);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      // element connectivity 1-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
        status = GmfSetLin(fileID, GmfTetrahedra, meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                                                  meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                                                  meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                                                  meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                                                  0); // to be consistent with refine
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    } else {
      AIM_ERROR(aimInfo, "libMeshb writer element type currently not implemented! group %d type = %d",
                igroup, meshData->elemGroups[igroup].elementTopo);
      status = CAPS_IOERR;
      goto cleanup;
    }
  }

  nNodeVerts = nEdgeVerts = nFaceVerts = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, NODE, &nNode, NULL);
    AIM_STATUS(aimInfo, status);
    nNodeVerts += nNode;

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      if (edges[iedge]->mtype == DEGENERATE) continue;
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      AIM_STATUS(aimInfo, status);
      nEdgeVerts += plen;
    }
    AIM_FREE(edges);

    if (meshRef->type == aimAreaMesh) continue;

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);
      nFaceVerts += plen;
    }
  }

  // write out parametric coordinates

  // Write NODEs
  status = GmfSetKwd(fileID, GmfVerticesOnGeometricVertices, nNodeVerts);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nNodeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, NODE, &nNode, NULL);
    AIM_STATUS(aimInfo, status);

    for (j = 0; j < nglobal; j++) {
      status = EG_getGlobal(meshRef->maps[i].tess, j + 1, &ptype, &pindex, xyz);
      AIM_STATUS(aimInfo, status);
      if (ptype == 0) {
        status = GmfSetLin(fileID, GmfVerticesOnGeometricVertices,
                           meshRef->maps[i].map[j],
                           nNodeOffset + pindex);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    nNodeOffset += nNode;
  }

  // Write EDGEs
  status = GmfSetKwd(fileID, GmfVerticesOnGeometricEdges, nEdgeVerts);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nEdgeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      if (status == EGADS_DEGEN) continue;
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < plen; j++) {

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 1, &iglobal);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        id = nEdgeOffset + iedge + 1;
        status = GmfSetLin(fileID, GmfVerticesOnGeometricEdges,
                           meshRef->maps[i].map[iglobal-1],
                           id,
                           t[j],
                           (double)id);  // refine abuse of dist
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    nEdgeOffset += nEdge;
  }


  if (meshRef->type == aimAreaMesh) {

    // Write Triangles
    status = GmfSetKwd(fileID, GmfVerticesOnGeometricTriangles, meshData->nVertex);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    status = EG_statusTessBody(meshRef->maps[0].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
    AIM_STATUS(aimInfo, status);

    for (i = 0; i < meshData->nVertex; i++) {

      status = EG_invEvaluate(faces[0], meshData->verts[i], param, xyz);
      AIM_STATUS(aimInfo, status);

      id = 1;
      status = GmfSetLin(fileID, GmfVerticesOnGeometricTriangles,
                         i+1,
                         id,
                         param[0], param[1],
                         (double)id); // refine abuse of dist
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

    AIM_FREE(faces);

  } else {

    // Write FACEs
    status = GmfSetKwd(fileID, GmfVerticesOnGeometricTriangles, nFaceVerts);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    nFaceOffset = 0;
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
      AIM_STATUS(aimInfo, status);

      for (iface = 0; iface < nFace; iface++) {
        status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        for (j = 0; j < plen; j++) {

          status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, j + 1, &iglobal);
          AIM_STATUS(aimInfo, status);

          id = nFaceOffset + iface + 1;
          status = GmfSetLin(fileID, GmfVerticesOnGeometricTriangles,
                             meshRef->maps[i].map[iglobal-1],
                             id,
                             uv[2*j+0], uv[2*j+1],
                             (double)id); // refine abuse of dist
          if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }
      }

      nFaceOffset += nFace;
    }
  }

  //printf("Finished writing libMeshb file\n\n");

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(edges);
  AIM_FREE(faces);
  if (fileID != 0) GmfCloseMesh(fileID);
  return status;
}
