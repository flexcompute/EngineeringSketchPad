/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Exodus 3D Mesh Writer Code
 *
 * *      Copyright 2014-2023, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

// https://sandialabs.github.io/seacas-docs/html/index.html

#include <string.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "exodusWriter.h"

#include <exodusII.h>

#include <map>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <array>
#include <cmath>

const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}

template<int N>
static int getFaceIndex(const int *tElemConn,
                        const int nFace, const std::array<int,N> *tExodusFaceMap,
                        std::array<int,N> &aFaceConn)
{
  for(int i = 0; i < nFace; i++)
  {
    std::array<int,N> tCurFaceConn;
    for(int j = 0; j < N; j++)
      tCurFaceConn[j] = tElemConn[tExodusFaceMap[i][j]];
    std::sort(tCurFaceConn.begin(), tCurFaceConn.end());
    if (tCurFaceConn == aFaceConn)
      return i+1;
  }
  return -1;
}


static int getNodesetTopos(void *aimInfo,
                           const int numTopo,
                           ego *topos,
                           aimMeshTessMap *map,
                           std::map<std::string, std::set<int> > &nodesetGroups)
{
  int status = CAPS_SUCCESS;
  int i, j, itopo, iglobal;
  int atype, alen;
  const int *ints;
  const double *reals;
  const char *string;

  int          plen, tlen;
  const int    *tris, *tric, *ptype, *pindex;
  const double *points, *uv;

  char *rest=NULL, *token=NULL;

  const char *attributeKey = "exNodeset";

  for (itopo = 0; itopo < numTopo; itopo++) {

    status = EG_attributeRet(topos[itopo], attributeKey, &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_NOTFOUND) continue;
    AIM_STATUS(aimInfo, status);

    if (atype != ATTRSTRING) {
        AIM_ERROR(aimInfo, "Attribute '%s' should be a single string", attributeKey);
        status = EGADS_ATTRERR;
        goto cleanup;
    }

    char copy[1025];
    strncpy(copy, string, 1024);
    rest = copy;
    while((token = strtok_r(rest, ";", &rest))) {

      std::string key(token);

      if (topos[itopo]->oclass == FACE) {
        /* get the Face indexing and store it away */
        status = EG_getTessFace(map->tess, itopo+1, &plen, &points, &uv, &ptype, &pindex,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);
        for (i = 0; i < tlen; i++) {
          for (j = 0; j < 3; j++) {
            status = EG_localToGlobal(map->tess, itopo+1, tris[3*i+j], &iglobal);
            AIM_STATUS(aimInfo, status);

            nodesetGroups[key].insert(map->map[iglobal-1]);
          }
        }
      } else if (topos[itopo]->oclass == EDGE) {
        /* get the Edge indexing and store it away */
        status = EG_getTessEdge(map->tess, itopo+1, &plen, &points, &uv);
        AIM_STATUS(aimInfo, status);

        for (i = 0; i < plen; i++) {
          status = EG_localToGlobal(map->tess, -(itopo+1), i+1, &iglobal);
          AIM_STATUS(aimInfo, status);

          nodesetGroups[key].insert(map->map[iglobal-1]);
        }
      } else if (topos[itopo]->oclass == NODE) {
        /* get the Node indexing and store it away */
        status = EG_localToGlobal(map->tess, 0, (itopo+1), &iglobal);
        AIM_STATUS(aimInfo, status);

        nodesetGroups[key].insert(map->map[iglobal-1]);
      }
    }
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}


static int getSidesetBounds(void *aimInfo,
                            aimMeshRef *meshRef,
                            std::vector<std::map<std::string, std::vector<int>>> &sidesetBounds)
{
  int status = CAPS_SUCCESS;
  int imap, iface;
  int atype, alen, state, np;
  const int *ints;
  const double *reals;
  const char *string;
  int numFace;
  ego *faces, body;

  const char *attributeKey = "capsBound";

  sidesetBounds.resize(meshRef->nmap);

  for (imap = 0; imap < meshRef->nmap; imap++) {
    status = EG_statusTessBody(meshRef->maps[imap].tess, &body, &state, &np);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &faces);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < numFace; iface++) {

      status = EG_attributeRet(faces[iface], attributeKey, &atype, &alen, &ints, &reals, &string);
      if (status == EGADS_NOTFOUND) continue;
      AIM_STATUS(aimInfo, status);

      if (atype != ATTRSTRING) {
        AIM_ERROR(aimInfo, "Attribute '%s' should be a single string", attributeKey);
        status = EGADS_ATTRERR;
        goto cleanup;
      }

      sidesetBounds[imap][std::string(string)].push_back(iface);
    }

    AIM_FREE(faces);
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int  status; // Function return status

  int CPU_word_size = sizeof(double);
  int IO_word_size = sizeof(double);
  int exoid = 0;

  // Incoming bodies
  const char *intents;
  ego *bodies = NULL;
  int numBody = 0;

  int i, j, k, d, iface, igroup;
  int nPoint, nElems, nBnds=0, nCellElem=0, nCellGroup=0, ID = 0;
  int blk, bnd, cellID = 1;

  ex_init_params par;
  double *x=NULL, *y=NULL, *z=NULL;
  char **coord_names = NULL, **block_names = NULL, **side_set_names=NULL, **node_set_names=NULL;
  const char *name=NULL;

  char title[MAX_STR_LENGTH+1];
  int64_t num_attr = 0;

  char filename[PATH_MAX];

  ego body;
  int state, np;
  int imap;

  int atype, alen;
  const int *ints = NULL, *tessFaceQuadMap=NULL;
  const double *reals = NULL;
  const char *str = NULL;

  int plen, tlen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;


  int xMeshConstant = (int)true, yMeshConstant = (int)true, zMeshConstant = (int)true;

  int numTopo;
  ego *topos = NULL;
  int topoTypes[3] = {FACE, EDGE, NODE};

  aimMeshData *meshData = NULL;
  aimMeshRef *meshRef = NULL;

  std::array<int,2> tLineNodes;
  std::array<int,3> tTriNodes;
  std::array<int,4> tQuadNodes;

  std::array<int,2> tExodusTriFaceMap[3] = {{0,1},{1,2},{2,0}};
  std::array<int,2> tExodusQuadFaceMap[4] = {{0,1},{1,2},{2,3},{3,0}};
  std::array<int,3> tExodusTetFaceMap[4] = {{0,1,3},{1,2,3},{2,0,3},{0,2,1}};

  std::map<std::array<int,2>, std::pair<int,const int*> > mLineToTriMap;
  std::map<std::array<int,2>, std::pair<int,const int*> > mLineToQuadMap;
  std::map<std::array<int,3>, std::pair<int,const int*> > mFaceToTriMap;
  std::map<std::array<int,4>, std::pair<int,const int*> > mFaceToQuadMap;
  std::map<std::array<int,3>, std::pair<int,const int*> > mFaceToTetMap;
  std::map<std::string, std::vector<int> > sidesetGroups;
  std::vector<std::map<std::string, std::vector<int>>> sidesetBounds;
  std::map<std::string, std::set<int> > nodesetGroups;
  typedef std::map<std::array<int,2>, std::pair<int,const int*> >::const_iterator sidesetLine_iterator;
  typedef std::map<std::string, std::vector<int>>::const_iterator sideset_iterator;
  typedef std::map<std::string, std::set<int> >::const_iterator nodeset_iterator;

  printf("\nWriting exodus file ....\n");

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  if (mesh->meshData->dim != 2 && mesh->meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  exoid = ex_create(filename, EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC,
                    &CPU_word_size, &IO_word_size);
  if (exoid <= 0) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }

  snprintf(title,MAX_STR_LENGTH,"CAPS Generated");
  ex_copy_string(par.title,title,MAX_STR_LENGTH+1);

  meshRef = mesh->meshRef;
  meshData = mesh->meshData;

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "CAPS Exodus writer only supports linear mesh elements at the moment! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // count the number of trace element groups
    if (meshData->dim == 2 || meshRef->type == aimSurfaceMesh) {

      if (meshData->elemGroups[igroup].elementTopo != aimLine &&
          meshData->elemGroups[igroup].elementTopo != aimTri &&
          meshData->elemGroups[igroup].elementTopo != aimQuad) {
        AIM_ERROR(aimInfo, "CAPS 2D/Surface Exodus writer only supports triangle/quad meshes at the moment",
                  igroup, meshData->elemGroups[igroup].order);
        status = CAPS_IOERR;
        goto cleanup;
      }

      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad) {
        nCellElem += meshData->elemGroups[igroup].nElems;
        nCellGroup++;
      }

    } else {

      if (meshData->elemGroups[igroup].elementTopo != aimTri &&
          meshData->elemGroups[igroup].elementTopo != aimTet) {
        AIM_ERROR(aimInfo, "CAPS 3D Exodus writer only supports tetrahedral meshes at the moment",
                  igroup, meshData->elemGroups[igroup].order);
        status = CAPS_IOERR;
        goto cleanup;
      }

      if (
#ifdef EXODUS_SIDE_BLOCKS
          meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad ||
#endif
          meshData->elemGroups[igroup].elementTopo == aimTet ||
          meshData->elemGroups[igroup].elementTopo == aimPyramid ||
          meshData->elemGroups[igroup].elementTopo == aimPrism ||
          meshData->elemGroups[igroup].elementTopo == aimHex) {
        nCellElem += meshData->elemGroups[igroup].nElems;
        nCellGroup++;
      }

    }
  }

  /* agglomerate side set groups based on semi-colon separated groupName */
  for (igroup = 0, bnd = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->dim == 2) {
      if (meshData->elemGroups[igroup].elementTopo != aimLine)
        continue;
    } else {
      if (meshData->elemGroups[igroup].elementTopo != aimTri &&
          meshData->elemGroups[igroup].elementTopo != aimQuad)
        continue;
    }

    char copy[1025];
    if (meshData->elemGroups[igroup].groupName != NULL) {
      char *rest, *token;

      strncpy(copy, meshData->elemGroups[igroup].groupName, 1024);
      rest = copy;
      while((token = strtok_r(rest, ";", &rest))) {
        sidesetGroups[std::string(token)].push_back(igroup);
      }
    } else {
      snprintf(copy, 104, "BndGroup%d", bnd);
      sidesetGroups[std::string(copy)].push_back(igroup);
    }
    bnd++;
  }

  /* get sidesets for capsBound attributes */
  status = getSidesetBounds(aimInfo, meshRef, sidesetBounds);
  AIM_STATUS(aimInfo, status);

  /* agglomerate node set groups based on semi-colon separated exNodeset */
  for (imap = 0; imap < meshRef->nmap; imap++) {
    status = EG_statusTessBody(meshRef->maps[imap].tess, &body, &state, &np);
    AIM_STATUS(aimInfo, status);

    for (i = 0; i < 3; i++) {
      status = EG_getBodyTopos(body, NULL, topoTypes[i], &numTopo, &topos);
      AIM_STATUS(aimInfo, status);
      if (numTopo == 0) continue;
      AIM_NOTNULL(topos, aimInfo, status);

      status = getNodesetTopos(aimInfo, numTopo, topos, &meshRef->maps[imap], nodesetGroups);
      AIM_STATUS(aimInfo, status);
      AIM_FREE(topos);
    }

  }

  nBnds = sidesetGroups.size();
  for (std::size_t i = 0; i < sidesetBounds.size(); i++)
    nBnds += sidesetBounds[i].size();

  par.num_dim       = meshData->dim;
  par.num_nodes     = meshData->nVertex;
  par.num_edge      = 0;
  par.num_edge_blk  = 0;
  par.num_face      = 0;
  par.num_face_blk  = 0;
  par.num_elem      = nCellElem;
  par.num_elem_blk  = nCellGroup;
  par.num_node_sets = (int)nodesetGroups.size();
  par.num_edge_sets = 0;
  par.num_face_sets = 0;
  par.num_side_sets = nBnds;
  par.num_elem_sets = 0;
  par.num_node_maps = 0;
  par.num_edge_maps = 0;
  par.num_face_maps = 0;
  par.num_elem_maps = 0;
  par.num_assembly  = 0;
  par.num_blob      = 0;

  status = ex_put_init_ext(exoid, &par);
  AIM_STATUS(aimInfo, status);

  AIM_ALLOC(x, meshData->nVertex, double, aimInfo, status);
  AIM_ALLOC(y, meshData->nVertex, double, aimInfo, status);
  if (meshData->dim == 3)
    AIM_ALLOC(z, meshData->nVertex, double, aimInfo, status);


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
      printf("Exodus expects 2D meshes in the x-y plane... attempting to rotate mesh!\n");

      if (xMeshConstant == (int) true && yMeshConstant == (int) false) {
        printf("Swapping z and x coordinates!\n");
        for (i = 0; i < meshData->nVertex; i++) {
          x[i] = meshData->verts[i][2];
          y[i] = meshData->verts[i][1];
        }

      } else if(xMeshConstant == (int) false && yMeshConstant == (int) true) {

        printf("Swapping z and y coordinates!\n");
        for (i = 0; i < meshData->nVertex; i++) {
          x[i] = meshData->verts[i][0];
          y[i] = meshData->verts[i][2];
        }

      } else {
        AIM_ERROR(aimInfo, "Unable to rotate mesh!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
      }

    } else { // zMeshConstant == true
      // Write nodal coordinates as is
      for (i = 0; i < meshData->nVertex; i++) {
        x[i] = meshData->verts[i][0];
        y[i] = meshData->verts[i][1];
      }
    }

  } else {

    for (i = 0; i < meshData->nVertex; i++)
    {
      x[i] = meshData->verts[i][0];
      y[i] = meshData->verts[i][1];
      z[i] = meshData->verts[i][2];
    }
  }

  /* write all of the vertices */
  status = ex_put_coord(exoid, x, y, z);
  AIM_STATUS(aimInfo, status);

  AIM_FREE(x);
  AIM_FREE(y);
  AIM_FREE(z);

  /* write coordinate names */
  AIM_ALLOC(coord_names, meshData->dim, char*, aimInfo, status);
  for (d = 0; d < meshData->dim; d++)
    coord_names[d] = NULL;

  for (d = 0; d < meshData->dim; d++)
    AIM_ALLOC(coord_names[d], MAX_STR_LENGTH+1, char, aimInfo, status);

  snprintf(coord_names[0], MAX_STR_LENGTH,"x");
  snprintf(coord_names[1], MAX_STR_LENGTH,"y");
  if (meshData->dim == 3) snprintf(coord_names[2], MAX_STR_LENGTH,"z");

  status = ex_put_coord_names(exoid, coord_names);
  AIM_STATUS(aimInfo, status);

  for (d = 0; d < meshData->dim; d++)
    AIM_FREE(coord_names[d]);
  AIM_FREE(coord_names);


  AIM_ALLOC(block_names, nCellGroup, char*, aimInfo, status);
  for (blk = 0; blk < nCellGroup; blk++)
    block_names[blk] = NULL;

  // if there is only one body then look for the "_name" to indicate the block_name
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);
  if (numBody == 1 && (meshRef->type == aimAreaMesh || meshRef->type == aimVolumeMesh)) {
    AIM_NOTNULL(bodies, aimInfo, status);
    EG_attributeRet(bodies[0], "_name", &atype, &alen, &ints, &reals, &str);
  }

  for (blk = 0, igroup = 0; igroup < meshData->nElemGroup; igroup++)
  {
    if (meshData->dim == 2 || meshRef->type == aimSurfaceMesh) {
      if (meshData->elemGroups[igroup].elementTopo == aimLine)
        continue;
    } else {
#ifndef EXODUS_SIDE_BLOCKS
      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad)
        continue;
#endif
    }

         if (meshData->elemGroups[igroup].elementTopo == aimTri    ) name = "tri";
    else if (meshData->elemGroups[igroup].elementTopo == aimQuad   )
      if (meshData->dim == 2) name = "quad";
      else                    name = "shell";
    else if (meshData->elemGroups[igroup].elementTopo == aimTet    ) name = "tet";
    else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) name = "pyramid";
    else if (meshData->elemGroups[igroup].elementTopo == aimPrism  ) name = "prism";
    else if (meshData->elemGroups[igroup].elementTopo == aimHex    ) name = "hex";

    AIM_ALLOC(block_names[blk], MAX_STR_LENGTH+1, char, aimInfo, status);

    if (meshData->elemGroups[igroup].groupName != NULL)
      snprintf(block_names[blk], MAX_STR_LENGTH,"%s_%s", meshData->elemGroups[igroup].groupName, name);
    else if (str != NULL)
      snprintf(block_names[blk], MAX_STR_LENGTH,"%s_%s", str, name);
    else
      snprintf(block_names[blk], MAX_STR_LENGTH,"Block%d_%s", blk+1, name);

    blk++;
  }

  /* write the block names */
  status = ex_put_names(exoid, EX_ELEM_BLOCK, block_names);
  AIM_STATUS(aimInfo, status);


  /* write the sideset names */
  AIM_ALLOC(side_set_names, nBnds, char*, aimInfo, status);
  for (bnd = 0; bnd < nBnds; bnd++)
    side_set_names[bnd] = NULL;

  bnd = 0;
  for (sideset_iterator sideset = sidesetGroups.begin(); sideset != sidesetGroups.end(); ++sideset) {
    AIM_ALLOC(side_set_names[bnd], MAX_STR_LENGTH+1, char, aimInfo, status);
    snprintf(side_set_names[bnd], MAX_STR_LENGTH, "%s", sideset->first.c_str());
    bnd++;
  }

  for (std::size_t imap = 0; imap < sidesetBounds.size(); imap++) {
    for (sideset_iterator sideset = sidesetBounds[imap].begin(); sideset != sidesetBounds[imap].end(); ++sideset) {
      AIM_ALLOC(side_set_names[bnd], MAX_STR_LENGTH+1, char, aimInfo, status);
      snprintf(side_set_names[bnd], MAX_STR_LENGTH, "%s", sideset->first.c_str());
      bnd++;
    }
  }

  if (!sidesetGroups.empty() || !sidesetBounds.empty()) {
    status = ex_put_names(exoid, EX_SIDE_SET, side_set_names);
    AIM_STATUS(aimInfo, status);
  }

  /* write the nodeset names */
  if (!nodesetGroups.empty()) {
    AIM_ALLOC(node_set_names, nodesetGroups.size(), char*, aimInfo, status);
    for (i = 0; i < (int)nodesetGroups.size(); i++)
      node_set_names[i] = NULL;

    i = 0;
    for (nodeset_iterator nodeset = nodesetGroups.begin(); nodeset != nodesetGroups.end(); ++nodeset) {
      AIM_ALLOC(node_set_names[i], MAX_STR_LENGTH+1, char, aimInfo, status);
      snprintf(node_set_names[i], MAX_STR_LENGTH, "%s", nodeset->first.c_str());
      i++;
    }

    status = ex_put_names(exoid, EX_NODE_SET, node_set_names);
    AIM_STATUS(aimInfo, status);
  }


  /* write volume elements blocks */
  cellID = 1;
  blk = 0;
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->dim == 2 || meshRef->type == aimSurfaceMesh) {
      if (meshData->elemGroups[igroup].elementTopo == aimLine)
        continue;
    } else {
#ifndef EXODUS_SIDE_BLOCKS
      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad)
        continue;
#endif
    }

    nPoint = meshData->elemGroups[igroup].nPoint;
    nElems = meshData->elemGroups[igroup].nElems;

         if (meshData->elemGroups[igroup].elementTopo == aimTri    ) name = "tri";
    else if (meshData->elemGroups[igroup].elementTopo == aimQuad   )
      if (meshData->dim == 2) name = "quad";
      else                    name = "shell";
    else if (meshData->elemGroups[igroup].elementTopo == aimTet    ) name = "tet";
    else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) name = "pyramid";
    else if (meshData->elemGroups[igroup].elementTopo == aimPrism  ) name = "prism";
    else if (meshData->elemGroups[igroup].elementTopo == aimHex    ) name = "hex";

    printf("\tBlock: '%s', blk_id %d, num_entries %d\n", block_names[blk], blk+1, nElems);

    /* write the block info */
    status = ex_put_block(exoid, EX_ELEM_BLOCK, blk+1, name, nElems, nPoint, 0, 0, num_attr);
    AIM_STATUS(aimInfo, status);

    /* write the element connectivity */
    status = ex_put_conn(exoid, EX_ELEM_BLOCK, blk+1, meshData->elemGroups[igroup].elements, NULL, NULL);
    AIM_STATUS(aimInfo, status);
    blk++;

    /* construct the trace to cellID map */
    for(i = 0; i < nElems; i++, cellID++)
    {
      const int *tElemConn = meshData->elemGroups[igroup].elements + nPoint*i;
      if (meshData->dim == 2) {
        if (meshData->elemGroups[igroup].elementTopo == aimTri) {
          for(j = 0; j < 3; j++) // 3 faces on a tri
          {
            tLineNodes[0] = tElemConn[tExodusTriFaceMap[j][0]];
            tLineNodes[1] = tElemConn[tExodusTriFaceMap[j][1]];
            std::sort(tLineNodes.begin(), tLineNodes.end());
            mLineToTriMap[tLineNodes] = std::pair<int,const int*>(cellID, tElemConn);
          }
        } else {
          for(j = 0; j < 4; j++) // 4 faces on a quad
          {
            tLineNodes[0] = tElemConn[tExodusQuadFaceMap[j][0]];
            tLineNodes[1] = tElemConn[tExodusQuadFaceMap[j][1]];
            std::sort(tLineNodes.begin(), tLineNodes.end());
            mLineToQuadMap[tLineNodes] = std::pair<int,const int*>(cellID, tElemConn);
          }
        }
      } else if (meshRef->type == aimSurfaceMesh) {
        if (meshData->elemGroups[igroup].elementTopo == aimTri) {
          for(j = 0; j < 3; j++) // 3 faces on a tri
            tTriNodes[j] = tElemConn[j];
          std::sort(tTriNodes.begin(), tTriNodes.end());
          mFaceToTriMap[tTriNodes] = std::pair<int,const int*>(cellID, tElemConn);
        } else {
          for(j = 0; j < 4; j++) // 4 faces on a quad
            tQuadNodes[j] = tElemConn[j];
          std::sort(tQuadNodes.begin(), tQuadNodes.end());
          mFaceToQuadMap[tQuadNodes] = std::pair<int,const int*>(cellID, tElemConn);
        }
      } else {
        if (meshData->elemGroups[igroup].elementTopo == aimTet) {
          for(j = 0; j < 4; j++) // 4 faces on a tet
          {
            tTriNodes[0] = tElemConn[tExodusTetFaceMap[j][0]];
            tTriNodes[1] = tElemConn[tExodusTetFaceMap[j][1]];
            tTriNodes[2] = tElemConn[tExodusTetFaceMap[j][2]];
            std::sort(tTriNodes.begin(), tTriNodes.end());
            mFaceToTetMap[tTriNodes] = std::pair<int,const int*>(cellID, tElemConn);
          }
        }
      }
    }
  }

  /* write side sets */
  ID = 1;
  bnd = 0;
  for (sideset_iterator sideset = sidesetGroups.begin(); sideset != sidesetGroups.end(); ++sideset) {

    nElems = 0;
    for (i = 0; i < (int)sideset->second.size(); i++) {
      igroup = sideset->second[i];
      nElems += meshData->elemGroups[igroup].nElems;
    }

    std::vector<int> elem_list(nElems,-1);
    std::vector<int> side_list(nElems,-1);

    printf("\tSideset: '%s', set_id %d, num_entries %d\n", side_set_names[bnd++], ID, nElems);

    /* write the side set info */
    status = ex_put_set_param(exoid, EX_SIDE_SET, ID, nElems, 0);
    AIM_STATUS(aimInfo, status);

    int ielem = 0;
    for (i = 0; i < (int)sideset->second.size(); i++) {
      igroup = sideset->second[i];
      nElems = meshData->elemGroups[igroup].nElems;
      nPoint = meshData->elemGroups[igroup].nPoint;

      for(j = 0; j < nElems; j++, ielem++) {
        const int *tElemConn = meshData->elemGroups[igroup].elements + nPoint*j;
        std::pair<int,const int*> tAttachedElem;

        if (meshData->dim == 2) {
          tLineNodes[0] = tElemConn[0];
          tLineNodes[1] = tElemConn[1];

          std::sort(tLineNodes.begin(), tLineNodes.end());

          sidesetLine_iterator it = mLineToTriMap.find(tLineNodes);
          if (it != mLineToTriMap.end()) {
            tAttachedElem = it->second;
            elem_list[ielem] = tAttachedElem.first;
            side_list[ielem] = getFaceIndex<2>(tAttachedElem.second, 3, tExodusTriFaceMap, tLineNodes);
            continue;
          }

          it = mLineToQuadMap.find(tLineNodes);
          if (it != mLineToQuadMap.end()) {
            tAttachedElem = it->second;
            elem_list[ielem] = tAttachedElem.first;
            side_list[ielem] = getFaceIndex<2>(tAttachedElem.second, 4, tExodusQuadFaceMap, tLineNodes);
            continue;
          }

          AIM_ERROR(aimInfo, "Failed to located block element for segment %d %d", tLineNodes[0], tLineNodes[1]);
          status = CAPS_BADVALUE;
          goto cleanup;
        } else if (meshRef->type == aimSurfaceMesh) {

          if (meshData->elemGroups[igroup].elementTopo == aimTri) {
            for(k = 0; k < 3; k++) // 3 faces on a tri
              tTriNodes[k] = tElemConn[k];
            std::sort(tTriNodes.begin(), tTriNodes.end());

            tAttachedElem = mFaceToTriMap.at(tTriNodes);
          } else {
            for(k = 0; k < 4; k++) // 4 faces on a quad
              tQuadNodes[k] = tElemConn[k];
            std::sort(tQuadNodes.begin(), tQuadNodes.end());

            tAttachedElem = mFaceToQuadMap.at(tQuadNodes);
          }

          elem_list[ielem] = tAttachedElem.first;
          side_list[ielem] = 1;
        } else {
          tTriNodes[0] = tElemConn[0];
          tTriNodes[1] = tElemConn[1];
          tTriNodes[2] = tElemConn[2];

          std::sort(tTriNodes.begin(), tTriNodes.end());
          tAttachedElem = mFaceToTetMap.at(tTriNodes);

          elem_list[ielem] = tAttachedElem.first;
          side_list[ielem] = getFaceIndex<3>(tAttachedElem.second, 4, tExodusTetFaceMap, tTriNodes);
        }
      }
    }

    /* write the side set elem and side lists */
    status = ex_put_set(exoid, EX_SIDE_SET, ID, elem_list.data(), side_list.data());
    AIM_STATUS(aimInfo, status);

    ID++;
  }

  for (std::size_t imap = 0; imap < sidesetBounds.size(); imap++) {

    // check if the tessellation has a mixture of quad and tri
    status = EG_attributeRet(meshRef->maps[imap].tess, ".mixed",
                             &atype, &alen, &tessFaceQuadMap, &reals, &str);
    AIM_NOTFOUND(aimInfo, status);


    for (sideset_iterator sideset = sidesetBounds[imap].begin(); sideset != sidesetBounds[imap].end(); ++sideset) {

      nElems = 0;
      for (i = 0; i < (int)sideset->second.size(); i++) {

        iface = sideset->second[i];

        /* get the Face indexing and store it away */
        status = EG_getTessFace(meshRef->maps[imap].tess, iface+1, &plen, &points, &uv, &ptype, &pindex,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        if (tessFaceQuadMap != NULL) {
          tlen -= 2*tessFaceQuadMap[iface];
          nElems += tessFaceQuadMap[iface];
          }
        nElems += tlen;
      }

      std::vector<int> elem_list(nElems,-1);
      std::vector<int> side_list(nElems,-1);

      printf("\tSideset: '%s', set_id %d, num_entries %d\n", side_set_names[bnd++], ID, nElems);

      /* write the side set info */
      status = ex_put_set_param(exoid, EX_SIDE_SET, ID, nElems, 0);
      AIM_STATUS(aimInfo, status);

      int ielem = 0;
      for (i = 0; i < (int)sideset->second.size(); i++) {
        iface = sideset->second[i];

        /* get the Face indexing and store it away */
        status = EG_getTessFace(meshRef->maps[imap].tess, iface+1, &plen, &points, &uv, &ptype, &pindex,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        if (tessFaceQuadMap != NULL) {
          tlen -= 2*tessFaceQuadMap[iface];
        }

        int stride = 0;
        for(j = 0; j < tlen; j++, stride += 3, ielem++) {
          std::pair<int,const int*> tAttachedElem;

          status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 0], &tTriNodes[0]);
          AIM_STATUS(aimInfo, status);
          status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 1], &tTriNodes[1]);
          AIM_STATUS(aimInfo, status);
          status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 2], &tTriNodes[2]);
          AIM_STATUS(aimInfo, status);

          std::sort(tTriNodes.begin(), tTriNodes.end());

          if (meshRef->type == aimSurfaceMesh) {
            tAttachedElem = mFaceToTriMap.at(tTriNodes);

            elem_list[ielem] = tAttachedElem.first;
            side_list[ielem] = 1;

          } else {
            tAttachedElem = mFaceToTetMap.at(tTriNodes);

            elem_list[ielem] = tAttachedElem.first;
            side_list[ielem] = getFaceIndex<3>(tAttachedElem.second, 4, tExodusTetFaceMap, tTriNodes);
          }
        }

        if (tessFaceQuadMap != NULL && tessFaceQuadMap[iface] > 0) {
          for(j = 0; j < tessFaceQuadMap[iface]; j++, stride += 6, ielem++) {
            std::pair<int,const int*> tAttachedElem;

            status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 0], &tQuadNodes[0]);
            AIM_STATUS(aimInfo, status);
            status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 1], &tQuadNodes[1]);
            AIM_STATUS(aimInfo, status);
            status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 2], &tQuadNodes[2]);
            AIM_STATUS(aimInfo, status);
            status = EG_localToGlobal(meshRef->maps[imap].tess, iface + 1, tris[stride + 5], &tQuadNodes[3]);
            AIM_STATUS(aimInfo, status);

            std::sort(tQuadNodes.begin(), tQuadNodes.end());

            if (meshRef->type == aimSurfaceMesh) {
              tAttachedElem = mFaceToQuadMap.at(tQuadNodes);

              elem_list[ielem] = tAttachedElem.first;
              side_list[ielem] = 1;

            } else {
              AIM_ERROR(aimInfo, "Quad side-set for Hex elements not yet supported!");
              status = CAPS_NOTIMPLEMENT;
              goto cleanup;

//              tAttachedElem = mFaceToTetMap.at(tTriNodes);
//
//              elem_list[ielem] = tAttachedElem.first;
//              side_list[ielem] = getFaceIndex<3>(tAttachedElem.second, 4, tExodusTetFaceMap, tTriNodes);
            }
          }
        }
      }

      /* write the side set elem and side lists */
      status = ex_put_set(exoid, EX_SIDE_SET, ID, elem_list.data(), side_list.data());
      AIM_STATUS(aimInfo, status);

      ID++;
    }
  }


  /* write node sets */
  ID = 1;
  bnd = 0;
  for (nodeset_iterator nodeset = nodesetGroups.begin(); nodeset != nodesetGroups.end(); ++nodeset) {

    std::vector<int> node_list(nodeset->second.begin(),nodeset->second.end());
    nElems = node_list.size();

    AIM_NOTNULL(node_set_names, aimInfo, status);
    AIM_NOTNULL(node_set_names[bnd], aimInfo, status);
    printf("\tNodeset: '%s', set_id %d, num_entries %d\n", node_set_names[bnd++], ID, nElems);

    /* write the side set info */
    status = ex_put_set_param(exoid, EX_NODE_SET, ID, nElems, 0);
    AIM_STATUS(aimInfo, status);

    /* write the side set elem and side lists */
    status = ex_put_set(exoid, EX_NODE_SET, ID, node_list.data(), NULL);
    AIM_STATUS(aimInfo, status);

    ID++;
  }

  printf("Finished writing Exodus file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (exoid > 0) ex_close(exoid);

  AIM_FREE(topos);

  AIM_FREE(x);
  AIM_FREE(y);
  AIM_FREE(z);
  if (coord_names != NULL) {
    for (d = 0; d < meshData->dim; d++)
      AIM_FREE(coord_names[d]);
    AIM_FREE(coord_names);
  }
  if (block_names != NULL) {
    for (blk = 0; blk < nCellGroup; blk++)
      AIM_FREE(block_names[blk]);
    AIM_FREE(block_names);
  }
  if (side_set_names != NULL) {
    for (bnd = 0; bnd < nBnds; bnd++)
      AIM_FREE(side_set_names[bnd]);
    AIM_FREE(side_set_names);
  }
  if (node_set_names != NULL) {
    for (i = 0; i < (int)nodesetGroups.size(); i++)
      AIM_FREE(node_set_names[i]);
    AIM_FREE(node_set_names);
  }

  return status;
}
