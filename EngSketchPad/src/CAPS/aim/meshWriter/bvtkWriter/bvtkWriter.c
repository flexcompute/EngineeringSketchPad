/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Binary VTK 3D Mesh Writer Code
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "aimUtil.h"
#include "aimMesh.h"

#include "bvtkWriter.h"

/* Legacy VTK binary file assumes data stored as BigEndian */
static
void byteswap(void *x, void *y, size_t size)
{
  size_t i;
  char *xp = (char *)x;
  char *yp = (char *)y;

  for (i = 0; i < size; i++)
    *(yp + size-1-i) = *(xp + i);
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
  int  igroup, i, d, ipnt, ielem, nPoint, nElems;
  int  elemType, iswap;
  int  nCell, length, iconn;
  char filename[PATH_MAX];
  double dswap;
  FILE *fp=NULL;
  aimMeshData *meshData = NULL;

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  if (mesh->meshData->dim != 2 && mesh->meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }
  
  meshData = mesh->meshData;

  nCell = 0;
  length = 0;
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "VTK only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // count the number of elements and connectivity table size
    nPoint = meshData->elemGroups[igroup].nPoint;
    nElems = meshData->elemGroups[igroup].nElems;

    nCell += nElems;
    length += nElems*(nPoint+1);
  }

  printf("\nWriting Binary VTK file ....\n");

  /* open the file */
  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  fp = fopen(filename, "wb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }

  /* header information */
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"CAPS Unstructured Grid\n");
  fprintf(fp,"BINARY\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n", meshData->nVertex);

  /* write all of the vertices */
  for (i = 0; i < meshData->nVertex; i++) {
    for (d = 0; d < 3; d++) {
      byteswap(&meshData->verts[i][d], &dswap, sizeof(double));
      status = fwrite(&dswap, sizeof(double), 1, fp);
      if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }
  }

  /*----------------------------*/
  /* write element connectivity */
  /*----------------------------*/
  fprintf(fp,"CELLS %d %d\n", nCell, length);

  if (meshData->elemMap != NULL) {

    /* Write elements in order from the mesh generator */
    for (i = 0; i < meshData->nTotalElems; i++) {

      igroup = meshData->elemMap[i][0];
      ielem  = meshData->elemMap[i][1];

      nPoint = meshData->elemGroups[igroup].nPoint;

      /* number of points in the element */
      byteswap(&nPoint, &iswap, sizeof(int));
      status = fwrite(&iswap, sizeof(int), 1, fp);
      if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      /* VTK is 0-based indexed */
      for (ipnt = 0; ipnt < nPoint; ipnt++) {
        iconn = meshData->elemGroups[igroup].elements[nPoint*ielem+ipnt]-1;
        byteswap(&iconn, &iswap, sizeof(int));
        status = fwrite(&iswap, sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

  } else {

    /* Write elements by grouping */
    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

      nPoint = meshData->elemGroups[igroup].nPoint;
      nElems = meshData->elemGroups[igroup].nElems;

      for (ielem = 0; ielem < nElems; ielem++) {

        /* number of points in the element */
        byteswap(&nPoint, &iswap, sizeof(int));
        status = fwrite(&iswap, sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

        /* VTK is 0-based indexed */
        for (ipnt = 0; ipnt < nPoint; ipnt++) {
          iconn = meshData->elemGroups[igroup].elements[nPoint*ielem+ipnt]-1;
          byteswap(&iconn, &iswap, sizeof(int));
          status = fwrite(&iswap, sizeof(int), 1, fp);
          if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        }
      }
    }

  }

  /*----------------------------*/
  /* write cell types           */
  /*----------------------------*/
  fprintf(fp,"CELL_TYPES %d\n", nCell);

  if (meshData->elemMap != NULL) {

    /* Write elements in order from the mesh generator */
    for (i = 0; i < meshData->nTotalElems; i++) {

      igroup = meshData->elemMap[i][0];

           if (meshData->elemGroups[igroup].elementTopo == aimLine)    elemType = 3;
      else if (meshData->elemGroups[igroup].elementTopo == aimTri)     elemType = 5;
      else if (meshData->elemGroups[igroup].elementTopo == aimQuad)    elemType = 9;
      else if (meshData->elemGroups[igroup].elementTopo == aimTet)     elemType = 10;
      else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) elemType = 14;
      else if (meshData->elemGroups[igroup].elementTopo == aimPrism)   elemType = 13;
      else if (meshData->elemGroups[igroup].elementTopo == aimHex)     elemType = 12;
      else { AIM_ERROR(aimInfo, "Unknown element type"); status = CAPS_NOTIMPLEMENT; goto cleanup; }

      byteswap(&elemType, &iswap, sizeof(int));

      /* write the element type */
      status = fwrite(&iswap, sizeof(int), 1, fp);
      if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

  } else {

    /* Write elements by grouping */
    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
           if (meshData->elemGroups[igroup].elementTopo == aimLine)    elemType = 3;
      else if (meshData->elemGroups[igroup].elementTopo == aimTri)     elemType = 5;
      else if (meshData->elemGroups[igroup].elementTopo == aimQuad)    elemType = 9;
      else if (meshData->elemGroups[igroup].elementTopo == aimTet)     elemType = 10;
      else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) elemType = 14;
      else if (meshData->elemGroups[igroup].elementTopo == aimPrism)   elemType = 13;
      else if (meshData->elemGroups[igroup].elementTopo == aimHex)     elemType = 12;
      else { AIM_ERROR(aimInfo, "Unknown element type"); status = CAPS_NOTIMPLEMENT; goto cleanup; }

      nElems = meshData->elemGroups[igroup].nElems;

      byteswap(&elemType, &iswap, sizeof(int));

      /* write the element type */
      for (ielem = 0; ielem < nElems; ielem++) {
        status = fwrite(&iswap, sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

  }

  /*----------------------------*/
  /* write cell IDs             */
  /*----------------------------*/
  fprintf(fp, "CELL_DATA %d\n", nCell);
  fprintf(fp, "SCALARS ID int 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  if (meshData->elemMap != NULL) {

    /* Write elements in order from the mesh generator */
    for (i = 0; i < meshData->nTotalElems; i++) {

      igroup = meshData->elemMap[i][0];

      /* write the element ID */
      byteswap(&meshData->elemGroups[igroup].ID, &iswap, sizeof(int));

      status = fwrite(&iswap, sizeof(int), 1, fp);
      if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

  } else {

    /* Write elements by grouping */
    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

      nElems = meshData->elemGroups[igroup].nElems;

      /* write the element ID */
      for (ielem = 0; ielem < nElems; ielem++) {
        byteswap(&meshData->elemGroups[igroup].ID, &iswap, sizeof(int));

        status = fwrite(&iswap, sizeof(int), 1, fp);
        if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

  }

  printf("Finished writing Binary VTK file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  return status;
}
