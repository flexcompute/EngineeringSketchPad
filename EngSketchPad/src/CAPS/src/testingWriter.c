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


const char *meshExtension()
{
/*@-observertrans@*/
  return ".txt";
/*@+observertrans@*/
}


int meshWrite(/*@unused@*/ void *aimInfo, aimMesh *mesh)
{
  int  i, len;
  char *full;
  FILE *fp;
  
  len = strlen(mesh->meshRef->fileName) + 5;
  
  full = (char *) malloc(len*sizeof(char));
  if (full == NULL) return EGADS_MALLOC;
  for (i = 0; i < len-5; i++) full[i] = mesh->meshRef->fileName[i];
  full[len-5] = '.';
  full[len-4] = 't';
  full[len-3] = 'x';
  full[len-2] = 't';
  full[len-1] =  0;
  
  fp = fopen(full, "w");
  if (fp == NULL) {
    printf(" testingWriter Error: Cannot open file %s\n", full);
    free(full);
    return CAPS_IOERR;
  }
  free(full);
  
  fprintf(fp, "Output by testingWriter so/DLL\n");
  fclose(fp);

  return CAPS_SUCCESS;
}
