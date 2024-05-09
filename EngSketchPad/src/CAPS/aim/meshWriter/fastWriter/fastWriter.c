/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             FAST 2D/3D Mesh Writer Code
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "fastWriter.h"


const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j, igroup, ielem, nPoint, nTri=0, nTet=0;
  char filename[PATH_MAX];
  FILE *fp=NULL;
  aimMeshData *meshData = NULL;

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);
  
  fp = fopen(filename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }
  
  meshData = mesh->meshData;

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "FAST only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    if (meshData->elemGroups[igroup].elementTopo == aimLine)
      continue;
    else if (meshData->elemGroups[igroup].elementTopo == aimTri)
      nTri += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimTet)
      nTet += meshData->elemGroups[igroup].nElems;
    else {
      AIM_ERROR(aimInfo, "FAST file format only supports Triangle and Tetrahedral elements!");
      status = CAPS_IOERR;
      goto cleanup;
    }
  }

  printf("\nWriting FAST file ....\n");

  // Number of nodes  Number of triangles Number of Tetra
  fprintf(fp, "%d %d %d\n", meshData->nVertex, nTri, nTet);

  // Write x coordinates for nodes
  for (i = 0; i < meshData->nVertex; i++) {
    fprintf(fp, "%.18e\n", meshData->verts[i][0]);
  }

  // Write y coordinates for nodes
  for (i = 0; i < meshData->nVertex; i++) {
    fprintf(fp, "%.18e\n", meshData->verts[i][1]);
  }

  // Write z coordinates for nodes
  for (i = 0; i < meshData->nVertex; i++) {
    fprintf(fp, "%.18e\n", meshData->verts[i][2]);
  }

  // Write connectivity for Triangles
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->elemGroups[igroup].elementTopo != aimTri) continue;

    for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {

      // element connectivity 1-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (j = 0; j < nPoint; j++ ) {
        fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]);
      }
      fprintf(fp, "\n");
    }
  }

  // Write out triangle group IDs
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->elemGroups[igroup].elementTopo != aimTri) continue;

    for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
      // element group ID
      fprintf(fp, "%d\n", meshData->elemGroups[igroup].ID);
    }
  }

  // Write connectivity for Tetrahedron
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->elemGroups[igroup].elementTopo != aimTet) continue;

    for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {

      // element connectivity 1-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (j = 0; j < nPoint; j++ ) {
        fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]);
      }
    }
    fprintf(fp, "\n");
  }

  printf("Finished writing FAST file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  return status;
}
