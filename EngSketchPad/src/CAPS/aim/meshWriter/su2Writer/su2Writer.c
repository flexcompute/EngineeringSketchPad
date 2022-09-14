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

#include "su2Writer.h"


const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}


int meshWrite(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j, igroup, ielem, nPoint, nBnds=0, elementType, nCellElem=0, elemID = 0;
  char filename[PATH_MAX];
  int xMeshConstant = (int)true, yMeshConstant = (int)true, zMeshConstant = (int)true;
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

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "SU2 only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // count the number of trace element groups
    if (meshData->dim == 2) {

      if (meshData->elemGroups[igroup].elementTopo == aimLine)
        nBnds += 1;

      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad)
        nCellElem += meshData->elemGroups[igroup].nElems;

    } else {

      if (meshData->elemGroups[igroup].elementTopo == aimTri ||
          meshData->elemGroups[igroup].elementTopo == aimQuad)
        nBnds += 1;

      if (meshData->elemGroups[igroup].elementTopo == aimTet ||
          meshData->elemGroups[igroup].elementTopo == aimPyramid ||
          meshData->elemGroups[igroup].elementTopo == aimPrism ||
          meshData->elemGroups[igroup].elementTopo == aimHex)
        nCellElem += meshData->elemGroups[igroup].nElems;

    }
  }

  // Important: SU2 wants elements/index to start at 0 we assume everything coming in starts at 1

  printf("\nWriting SU2 file ....\n");

  // Dimensionality
  fprintf(fp,"NDIME= %d\n", meshData->dim);

  // Number of elements
  fprintf(fp,"NELEM= %d\n", nCellElem);

  // SU2 wants elements/index to start at 0 - assume everything starts at 1
  for (i = 0; i < meshData->nTotalElems; i++) {

      igroup = meshData->elemMap[i][0];
      ielem  = meshData->elemMap[i][1];

      if (meshData->dim == 2) {

        if (meshData->elemGroups[igroup].elementTopo == aimLine)
          continue;

        if      ( meshData->elemGroups[igroup].elementTopo == aimTri)  elementType = 5;
        else if ( meshData->elemGroups[igroup].elementTopo == aimQuad) elementType = 9;
        else {
            AIM_ERROR(aimInfo, "Unrecognized elementTopo %d for SU2!",
                      meshData->elemGroups[igroup].elementTopo);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

      } else {
        if (meshData->elemGroups[igroup].elementTopo == aimTri ||
            meshData->elemGroups[igroup].elementTopo == aimQuad)
          continue;

        if      ( meshData->elemGroups[igroup].elementTopo == aimTet)     elementType = 10;
        else if ( meshData->elemGroups[igroup].elementTopo == aimPyramid) elementType = 14;
        else if ( meshData->elemGroups[igroup].elementTopo == aimPrism)   elementType = 13;
        else if ( meshData->elemGroups[igroup].elementTopo == aimHex)     elementType = 12;
        else {
            AIM_ERROR(aimInfo, "Unrecognized elementTopo %d for SU2!",
                      meshData->elemGroups[igroup].elementTopo);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
      }

      // element type ID
      fprintf(fp, "%d ", elementType);

      // element connectivity 0-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (j = 0; j < nPoint; j++ ) {
        fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]-1);
      }

      // element ID
      fprintf(fp, "%d\n", elemID++);
  }

  // Number of points
  fprintf(fp,"NPOIN= %d\n", meshData->nVertex);

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
      printf("SU2 expects 2D meshes be in the x-y plane... attempting to rotate mesh!\n");

      if (xMeshConstant == (int) true && yMeshConstant == (int) false) {
        printf("Swapping z and x coordinates!\n");
        for (i = 0; i < meshData->nVertex; i++)
          fprintf(fp,"%.18e %.18e %d\n", meshData->verts[i][2],
                                         meshData->verts[i][1],
                                         i);  // Connectivity starts at 0

      } else if(xMeshConstant == (int) false && yMeshConstant == (int) true) {

        printf("Swapping z and y coordinates!\n");
        for (i = 0; i < meshData->nVertex; i++)
          fprintf(fp,"%.18e %.18e %d\n", meshData->verts[i][0],
                                         meshData->verts[i][2],
                                         i);  // Connectivity starts at 0

      } else {
        AIM_ERROR(aimInfo, "Unable to rotate mesh!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
      }

    } else { // zMeshConstant == true
      // Write nodal coordinates as is
      for (i = 0; i < meshData->nVertex; i++)
        fprintf(fp,"%.18e %.18e %d\n", meshData->verts[i][0],
                                       meshData->verts[i][1],
                                       i);  // Connectivity starts at 0
    }

  } else {

    // Write nodal coordinates
    for (i = 0; i < meshData->nVertex; i++)
      fprintf(fp,"%.18e %.18e %.18e %d\n", meshData->verts[i][0],
                                           meshData->verts[i][1],
                                           meshData->verts[i][2],
                                           i);  // Connectivity starts at 0
  }

  // Number of boundary groups
  fprintf(fp,"NMARK= %d\n", nBnds);

  // Finally lets write the connectivity for the boundaries
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->dim == 2) {

      if ( meshData->elemGroups[igroup].elementTopo == aimTri ||
           meshData->elemGroups[igroup].elementTopo == aimQuad) continue;

      if ( meshData->elemGroups[igroup].elementTopo == aimLine) elementType = 3;
      else {
          AIM_ERROR(aimInfo, "Unrecognized trace elementTopo %d for SU2!",
                    meshData->elemGroups[igroup].elementTopo);
          status = CAPS_BADVALUE;
          goto cleanup;
      }

    } else {
      if ( meshData->elemGroups[igroup].elementTopo == aimTet ||
           meshData->elemGroups[igroup].elementTopo == aimPyramid ||
           meshData->elemGroups[igroup].elementTopo == aimPrism ||
           meshData->elemGroups[igroup].elementTopo == aimHex) continue;

      if      ( meshData->elemGroups[igroup].elementTopo == aimTri)  elementType = 5;
      else if ( meshData->elemGroups[igroup].elementTopo == aimQuad) elementType = 9;
      else {
          AIM_ERROR(aimInfo, "Unrecognized trace elementTopo %d for SU2!",
                    meshData->elemGroups[igroup].elementTopo);
          status = CAPS_BADVALUE;
          goto cleanup;
      }
    }

    // Probably eventually want to change this to a string tag
    fprintf(fp,"MARKER_TAG= BC_%d\n", meshData->elemGroups[igroup].ID);
    // see note at the beginning of function

    //Number of elements with a particular ID
    fprintf(fp,"MARKER_ELEMS= %d\n", meshData->elemGroups[igroup].nElems);

    for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {

      // element type ID
      fprintf(fp, "%d ", elementType);

      // element connectivity 0-based
      nPoint = meshData->elemGroups[igroup].nPoint;
      for (j = 0; j < nPoint; j++ ) {
        fprintf(fp, "%d ", meshData->elemGroups[igroup].elements[nPoint*ielem+j]-1);
      }
      fprintf(fp, "\n");
    }
  }

  printf("Finished writing SU2 file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  return status;
}
