/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Wavefront OBJ Surface Mesh Writer Code
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "wavefrontWriter.h"

const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j;
  int  igroup, jgroup, kgroup, ielem, nPoint;
  int *nGrouping=NULL, **grouping=NULL, *nElemGroup=NULL;
  char filename[PATH_MAX];
  FILE *fp=NULL;
  aimMeshData *meshData = NULL;

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

  AIM_ALLOC(nGrouping, meshData->nElemGroup, int, aimInfo, status);
  AIM_ALLOC(nElemGroup, meshData->nElemGroup, int, aimInfo, status);
  AIM_ALLOC(grouping, meshData->nElemGroup, int*, aimInfo, status);
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    nGrouping[igroup] = 1;
    nElemGroup[igroup] = meshData->elemGroups[igroup].nElems;
    grouping[igroup] = NULL;
  }
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    AIM_ALLOC(grouping[igroup], meshData->nElemGroup, int, aimInfo, status);
    grouping[igroup][0] = igroup;
  }

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "Wavefront writer only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* Find all groups with the same ID and wrinte them together */
    for (jgroup = igroup+1; jgroup < meshData->nElemGroup; jgroup++) {

      if (meshData->elemGroups[igroup].ID == meshData->elemGroups[jgroup].ID &&
          aim_elemTopoDim( meshData->elemGroups[igroup].elementTopo ) ==
          aim_elemTopoDim( meshData->elemGroups[jgroup].elementTopo ) ) {
        if (nGrouping[jgroup] == 0) continue;
        grouping[igroup][nGrouping[igroup]++] = jgroup;
        nElemGroup[igroup] += meshData->elemGroups[jgroup].nElems;
        grouping[jgroup][0] = -1;
        nGrouping[jgroup] = 0;
        nElemGroup[jgroup] = 0;
      }
    }
  }


  printf("\nWriting Wavefront OBJ file ....\n");

  fprintf(fp,"# Mesh file written by CAPS\n");

  // Write vertex coordinates
  for (i = 0; i < meshData->nVertex; i++)
    fprintf(fp,"v %.18e %.18e %.18e\n", meshData->verts[i][0],
                                        meshData->verts[i][1],
                                        meshData->verts[i][2]);

  /* Write elements by groupName */
  for (jgroup = 0; jgroup < meshData->nElemGroup; jgroup++) {

    if (grouping[jgroup][0] == -1) continue;

    igroup = grouping[jgroup][0];

    fprintf(fp, "\ng %s group%d\n", meshData->elemGroups[igroup].groupName,
            meshData->elemGroups[igroup].ID);

    for (kgroup = 0; kgroup < nGrouping[jgroup]; kgroup++) {
      igroup = grouping[jgroup][kgroup];

      nPoint = meshData->elemGroups[igroup].nPoint;

      // element connectivity 1-based
      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad) {

        for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
          fprintf(fp, "f ");
          for (j = 0; j < nPoint; j++ ) {
            fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]);
          }
          fprintf(fp, "\n");
        }
      }
    }
  }

  printf("Finished writing Wavefront OBJ file\n\n");

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(nGrouping);
  AIM_FREE(nElemGroup);
  if (grouping != NULL) {
    for (igroup = 0; igroup < meshData->nElemGroup; igroup++)
      AIM_FREE(grouping[igroup]);
  }
  AIM_FREE(grouping);

  if (fp != NULL) fclose(fp);
  return status;
}
