/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Tecplot Mesh Writer Code
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "tecplotWriter.h"


static void writeZone(FILE *fp, int dim, int nnode, int nElems, const aimMeshElemGroup* elemGroup)
{
  const char *elemType = NULL;
  char title[128];
  int d;

  switch (elemGroup->elementTopo) {
  case aimLine:
    elemType = "FELINESEG";
    break;
  case aimTri:
  case aimQuad:
    elemType = "FEQUADRILATERAL";
    break;
  case aimTet:
  case aimPyramid:
  case aimPrism:
  case aimHex:
    elemType = "FEBRICK";
    break;
  case aimUnknownElem:
    elemType = "UNKNOWN";
    break;
  }

  if (elemGroup->groupName != NULL) {
    snprintf(title, 128, "%s", elemGroup->groupName);
  } else {
    snprintf(title, 128, "ID_%d", elemGroup->ID);
  }

  fprintf( fp, "ZONE T=\"%s\", DATAPACKING=BLOCK, NODES=%d, ELEMENTS=%d, ZONETYPE=%s", title, nnode, nElems, elemType );

  fprintf( fp, ", VARLOCATION=([");
  for (d = 0; d < dim-1; d++)
    fprintf( fp, "%d,", d+1);
  fprintf( fp, "%d]=NODAL)", dim);

  fprintf( fp, ", DT=(");
  for (d = 0; d < dim; d++)
    fprintf( fp, "DOUBLE ");
  fprintf( fp, ")\n");
}


const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j, d, igroup, jgroup, kgroup, ielem, nPoint;
  char filename[PATH_MAX];
  int xMeshConstant = (int)true, yMeshConstant = (int)true, zMeshConstant = (int)true;
  int *nGrouping=NULL, **grouping=NULL, *nElemGroup=NULL;
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
      AIM_ERROR(aimInfo, "Tecplot writer only supports linear mesh elements! group %d order = %d",
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

  printf("\nWriting Tecplot file ....\n");

  fprintf(fp,"VARIABLES = \"X\", \"Y\"");
  if (meshData->dim == 2)
    fprintf(fp,"\n");
  else
    fprintf(fp,", \"Z\"\n");


  for (jgroup = 0; jgroup < meshData->nElemGroup; jgroup++) {

    if (grouping[jgroup][0] == -1) continue;

    igroup = grouping[jgroup][0];
    writeZone(fp, meshData->dim, meshData->nVertex, nElemGroup[igroup], &meshData->elemGroups[igroup]);

    /* Write coordinates for first zone, and share for remaining */
    if (jgroup == 0) {

      if (meshData->dim == 2) {

        for (i = 0; i < meshData->nVertex; i++) {
          // Constant x?
          if ((meshData->verts[i][0] - meshData->verts[0][0]) > 1E-7) {
            xMeshConstant = (int) false;
          }

          // Constant y?
          if ((meshData->verts[i][1] - meshData->verts[0][1] ) > 1E-7) {
            yMeshConstant = (int) false;
          }

          // Constant z?
          if ((meshData->verts[i][2] - meshData->verts[0][2]) > 1E-7) {
            zMeshConstant = (int) false;
          }
        }

        if (zMeshConstant == (int) false) {
          printf("Tecplot expects 2D meshes in the x-y plane... attempting to rotate mesh!\n");

          if (xMeshConstant == (int) true && yMeshConstant == (int) false) {
            printf("Swapping z and x coordinates!\n");
            for (i = 0; i < meshData->nVertex; i++)
              fprintf(fp,"%.18e\n", meshData->verts[i][2]);

            for (i = 0; i < meshData->nVertex; i++)
              fprintf(fp,"%.18e\n", meshData->verts[i][1]);

          } else if(xMeshConstant == (int) false && yMeshConstant == (int) true) {

            printf("Swapping z and y coordinates!\n");
            for (i = 0; i < meshData->nVertex; i++)
              fprintf(fp,"%.18e\n", meshData->verts[i][0]);

            for (i = 0; i < meshData->nVertex; i++)
              fprintf(fp,"%.18e\n", meshData->verts[i][2]);

          } else {
            AIM_ERROR(aimInfo, "Unable to rotate mesh!\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
          }

        } else { // zMeshConstant == true
          // Write nodal coordinates as is
          for (i = 0; i < meshData->nVertex; i++)
            fprintf(fp,"%.18e\n", meshData->verts[i][0]);
          for (i = 0; i < meshData->nVertex; i++)
            fprintf(fp,"%.18e\n", meshData->verts[i][1]);
        }

      } else {

        // Write nodal coordinates
        for (i = 0; i < meshData->nVertex; i++)
          fprintf(fp,"%.18e\n", meshData->verts[i][0]);
        for (i = 0; i < meshData->nVertex; i++)
          fprintf(fp,"%.18e\n", meshData->verts[i][1]);
        for (i = 0; i < meshData->nVertex; i++)
          fprintf(fp,"%.18e\n", meshData->verts[i][2]);
      }

    } else {
      // Share nodal coordinates
      fprintf( fp, ", VARSHARELIST=([");
      for (d = 0; d < meshData->dim-1; d++)
        fprintf( fp, "%d,", d+1);
      fprintf( fp, "%d]=1)\n", meshData->dim);
    }

    for (kgroup = 0; kgroup < nGrouping[jgroup]; kgroup++) {
      igroup = grouping[jgroup][kgroup];

      nPoint = meshData->elemGroups[igroup].nPoint;

      // element connectivity 1-based
      if (meshData->elemGroups[igroup].elementTopo == aimTri) {

         for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
           fprintf(fp,"%d %d %d %d\n",
                   meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+2]);
         }

      } else if (meshData->elemGroups[igroup].elementTopo == aimTet) {

         for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
           fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                   meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                   meshData->elemGroups[igroup].elements[nPoint*ielem+3]);
         }

      } else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) {

        for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
          fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                  meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+4],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+4],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+4],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+4]);
        }

      } else if (meshData->elemGroups[igroup].elementTopo == aimPrism) {

        for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
          fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                  meshData->elemGroups[igroup].elements[nPoint*ielem+0],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+1],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+2],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+3],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+4],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+5],
                  meshData->elemGroups[igroup].elements[nPoint*ielem+5]);
        }

      } else { // Line, Quad, Hex

        for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
          for (j = 0; j < nPoint; j++ ) {
            fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]);
          }
          fprintf(fp, "\n");
        }
      }
    }
  }

  printf("Finished writing Tecplot file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  AIM_FREE(nGrouping);
  AIM_FREE(nElemGroup);
  if (grouping != NULL) {
    for (igroup = 0; igroup < meshData->nElemGroup; igroup++)
      AIM_FREE(grouping[igroup]);
  }
  AIM_FREE(grouping);

  return status;
}
