/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Binary STL 3D Surface Mesh Writer Code
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

#include "bstlWriter.h"

static
void triNormal(double *p1,
               double *p2,
               double *p3,
               double norm[])
{
    double a[3]={0.0,0.0,0.0};
    double b[3]={0.0,0.0,0.0};
    double mag = 0.0;

    // a x b

    // Create two vectors from points
    a[0]= p2[0]-p1[0];
    a[1]= p2[1]-p1[1];
    a[2]= p2[2]-p1[2];

    b[0]= p3[0]-p1[0];
    b[1]= p3[1]-p1[1];
    b[2]= p3[2]-p1[2];

    // Take the cross product
    norm[0] = a[1]*b[2]-a[2]*b[1];

    norm[1] = a[2]*b[0]-a[0]*b[2];

    norm[2] = a[0]*b[1]-a[1]*b[0];

    // Normalize vector
    mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

    norm[0] = norm[0]/mag;
    norm[1] = norm[1]/mag;
    norm[2] = norm[2]/mag;
}

static
void writeFloat(double *p, FILE *fp)
{
  float f[3] = {p[0], p[1], p[2]};
  fwrite(f, sizeof(float), 3, fp);
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
  int  igroup, ielem, nPoint, itri, nTri=0;
  short attr = 0;
  char filename[PATH_MAX];
  double norm[3] = {0,0,0};
  double *p1, *p2, *p3;
  FILE *fp=NULL;
  aimMeshData *meshData = NULL;
  char header[80] = "CAPS_STL";

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  if (mesh->meshData->dim != 2 && mesh->meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);
  
  fp = fopen(filename, "wb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }
  
  meshData = mesh->meshData;

  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
    if (meshData->elemGroups[igroup].order != 1) {
      AIM_ERROR(aimInfo, "STL only supports linear mesh elements! group %d order = %d",
                igroup, meshData->elemGroups[igroup].order);
      status = CAPS_IOERR;
      goto cleanup;
    }

    if (meshData->elemGroups[igroup].elementTopo == aimTri)
      nTri += meshData->elemGroups[igroup].nElems;
    else if (meshData->elemGroups[igroup].elementTopo == aimQuad)
      nTri += 2*meshData->elemGroups[igroup].nElems;
  }

  printf("\nWriting Binary STL file ....\n");

  fwrite(header, sizeof(char), 80, fp);
  fwrite(&nTri, sizeof(int), 1, fp);

  // Write Normals and Points (very redundant file format...)
  for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

    if (meshData->elemGroups[igroup].elementTopo == aimTri) {

      nPoint = meshData->elemGroups[igroup].nPoint;

      for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {

        p1 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+0]-1][0];
        p2 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+1]-1][0];
        p3 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+2]-1][0];

        triNormal(p1, p2, p3, norm);

        writeFloat(norm, fp);
        writeFloat(p1, fp);
        writeFloat(p2, fp);
        writeFloat(p3, fp);
        fwrite(&attr, sizeof(short), 1, fp);
      }

    } else if (meshData->elemGroups[igroup].elementTopo == aimQuad) {

      nPoint = meshData->elemGroups[igroup].nPoint;

      for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {

        for (itri = 0; itri < 2; itri++) {
          if (itri == 0) {
            p1 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+0]-1][0];
            p2 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+1]-1][0];
            p3 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+2]-1][0];
          } else {
            p1 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+0]-1][0];
            p2 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+2]-1][0];
            p3 = &meshData->verts[meshData->elemGroups[igroup].elements[nPoint*ielem+3]-1][0];
          }

          triNormal(p1, p2, p3, norm);

          writeFloat(norm, fp);
          writeFloat(p1, fp);
          writeFloat(p2, fp);
          writeFloat(p3, fp);
          fwrite(&attr, sizeof(short), 1, fp);
        }
      }
    }
  }

  printf("Finished writing Binary STL file\n\n");

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  return status;
}
