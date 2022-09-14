/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Testing AIM 3D Mesh Writer Example Code
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "ugridWriter.h"


const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int  status; // Function return status
  int  nLine=0, nTri=0, nQuad=0;
  int  nTet=0, nPyramid=0, nPrism=0, nHex=0;
  int  igroup, ielem, nPoint, elementTopo=aimUnknownElem, ielemTopo, elemID = 0, nElems;
  char filename[PATH_MAX];
  FILE *fp=NULL;
  aimMeshData *meshData = NULL;

  printf("\nWriting ugrid file ....\n");

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  if (mesh->meshData->dim != 2 && mesh->meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  fp = fopen(filename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }

  meshData = mesh->meshData;

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "SU2 only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // count the number of element types
         if (meshData->elemGroups[igroup].elementTopo == aimLine   ) nLine    += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimTri    ) nTri     += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimQuad   ) nQuad    += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimTet    ) nTet     += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) nPyramid += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimPrism  ) nPrism   += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimHex    ) nHex     += meshData->elemGroups[igroup].nElems;
    else {
      AIM_ERROR(aimInfo, "Unknown group %d element topology: %d", igroup+1, meshData->elemGroups[igroup].elementTopo);
      status = CAPS_MISMATCH;
      goto cleanup;
    }
  }

  /* write a binary UGRID file */
  status = fwrite(&meshData->nVertex, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nTri,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nQuad,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nTet,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nPyramid, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nPrism,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nHex,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write all of the vertices */
  status = fwrite(meshData->verts, sizeof(aimMeshCoords), meshData->nVertex, fp);
  if (status != meshData->nVertex) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write triangles and quads */
  for (ielemTopo = 0; ielemTopo < 2; ielemTopo++) {
         if (ielemTopo == 0) elementTopo = aimTri;
    else if (ielemTopo == 1) elementTopo = aimQuad;

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
      if (meshData->elemGroups[igroup].elementTopo != elementTopo) continue;

      nPoint = meshData->elemGroups[igroup].nPoint;
      nElems = meshData->elemGroups[igroup].nElems;

      /* write the element connectivity */
      status = fwrite(meshData->elemGroups[igroup].elements, sizeof(int), nPoint*nElems, fp);
      if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }
  }

  /* write triangles and quads IDs*/
  for (ielemTopo = 0; ielemTopo < 2; ielemTopo++) {
         if (ielemTopo == 0) elementTopo = aimTri;
    else if (ielemTopo == 1) elementTopo = aimQuad;

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
      if (meshData->elemGroups[igroup].elementTopo != elementTopo) continue;

      nElems = meshData->elemGroups[igroup].nElems;
      elemID = meshData->elemGroups[igroup].ID;

      /* write the element ID */
      for (ielem = 0; ielem < nElems; ielem++) {
        status = fwrite(&elemID,  sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }
  }

  /* write volume elements*/
  for (ielemTopo = 0; ielemTopo < 4; ielemTopo++) {
         if (ielemTopo == 0) elementTopo = aimTet;
    else if (ielemTopo == 1) elementTopo = aimPyramid;
    else if (ielemTopo == 2) elementTopo = aimPrism;
    else if (ielemTopo == 3) elementTopo = aimHex;

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
      if (meshData->elemGroups[igroup].elementTopo != elementTopo) continue;

      nPoint = meshData->elemGroups[igroup].nPoint;
      nElems = meshData->elemGroups[igroup].nElems;

      /* write the element connectivity */
      status = fwrite(meshData->elemGroups[igroup].elements, sizeof(int), nPoint*nElems, fp);
      if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }
  }

  if (nTet+nPyramid+nPrism+nHex == 0) {

    status = fwrite(&nLine, sizeof(int), 1, fp);
    if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
      if (meshData->elemGroups[igroup].elementTopo != aimLine) continue;

      nPoint = meshData->elemGroups[igroup].nPoint;
      nElems = meshData->elemGroups[igroup].nElems;
      elemID = meshData->elemGroups[igroup].ID;

      /* write the element connectivity and ID */
      for (ielem = 0; ielem < nElems; ielem++) {
        status = fwrite(&meshData->elemGroups[igroup].elements[nPoint*ielem], sizeof(int), nPoint, fp);
        if (status != nPoint) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        status = fwrite(&elemID, sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }
  }

  printf("Finished writing SU2 file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  return status;
}
