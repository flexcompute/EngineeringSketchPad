/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Volume Mesh Functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>

#ifdef WIN32
#include <Windows.h>
#include <direct.h>
#define strcasecmp stricmp
#define F_OK       0
#else
#include <dirent.h>
#include <dlfcn.h>
#include <unistd.h>
#endif

#include "aimUtil.h"
#include "aimMesh.h"


static /*@null@*/ DLL writerDLopen(const char *name)
{
  int             i, len;
  DLL             dll;
  char            *full, *env;
#ifdef WIN32
  WIN32_FIND_DATA ffd;
  char            dir[MAX_PATH];
  size_t          length_of_arg;
  HANDLE hFind =  INVALID_HANDLE_VALUE;
#else
/*@-unrecog@*/
  char            dir[PATH_MAX];
/*@+unrecog@*/
  struct dirent   *de;
  DIR             *dr;
#endif

  env = getenv("ESP_ROOT");
  if (env == NULL) {
    printf(" Information: Could not find $ESP_ROOT\n");
    return NULL;
  }

  if (name == NULL) {
    printf(" Information: Dynamic Loader invoked with NULL name!\n");
    return NULL;
  }
  len  = strlen(name);
  full = (char *) malloc((len+5)*sizeof(char));
  if (full == NULL) {
    printf(" Information: Dynamic Loader MALLOC Error!\n");
    return NULL;
  }
  for (i = 0; i < len; i++) full[i] = name[i];
  full[len  ] = '.';
#ifdef WIN32
  full[len+1] = 'D';
  full[len+2] = 'L';
  full[len+3] = 'L';
  full[len+4] =  0;
  _snprintf(dir, MAX_PATH, "%s\\lib\\*", env);
  hFind = FindFirstFile(dir, &ffd);
  if (INVALID_HANDLE_VALUE == hFind) {
    printf(" Information: Dynamic Loader could not open %s\n", dir);
    free(full);
    return NULL;
  }
  i = 0;
  do {
    if (strcasecmp(full, ffd.cFileName) == 0) i++;
  } while (FindNextFile(hFind, &ffd) != 0);
  FindClose(hFind);
  if (i > 1) {
    printf(" Information: Dynamic Loader more than 1 file: %s!\n", full);
    free(full);
    return NULL;
  }

  if (i == 1) {
    hFind = FindFirstFile(dir, &ffd);
    do {
      if (strcasecmp(full, ffd.cFileName) == 0) break;
    } while (FindNextFile(hFind, &ffd) != 0);
    dll = LoadLibrary(ffd.cFileName);
    FindClose(hFind);
  } else {
    dll = LoadLibrary(full);
  }

#else

  full[len+1] = 's';
  full[len+2] = 'o';
  full[len+3] =  0;
  snprintf(dir, PATH_MAX, "%s/lib", env);
  dr = opendir(dir);
  if (dr == NULL) {
    printf(" Information: Dynamic Loader could not open %s\n", dir);
    free(full);
    return NULL;
  }
  i = 0;
  while ((de = readdir(dr)) != NULL)
/*@-unrecog@*/
    if (strcasecmp(full, de->d_name) == 0) i++;
/*@+unrecog@*/
  closedir(dr);
  if (i > 1) {
    printf(" Information: Dynamic Loader more than 1 file: %s!\n", full);
    free(full);
    return NULL;
  }

  if (i == 1) {
    dr = opendir(dir);
    while ((de = readdir(dr)) != NULL)
      if (strcasecmp(full, de->d_name) == 0) break;
/*@-nullderef@*/
    dll = dlopen(de->d_name, RTLD_NOW | RTLD_NODELETE);
/*@-nullderef@*/
    closedir(dr);
  } else {
    dll = dlopen(full, RTLD_NOW | RTLD_NODELETE);
  }
#endif

  if (!dll) {
    printf(" Information: Dynamic Loader Error for %s\n", full);
#ifndef WIN32
/*@-mustfreefresh@*/
    printf("              %s\n", dlerror());
/*@+mustfreefresh@*/
#endif
    free(full);
    return NULL;
  }
  free(full);

  return dll;
}


static void writerDLclose(/*@only@*/ DLL dll)
{
#ifdef WIN32
  FreeLibrary(dll);
#else
  dlclose(dll);
#endif
}


static wrDLLFunc writerDLget(DLL dll, const char *symname)
{
  wrDLLFunc data;

#ifdef WIN32
  data = (wrDLLFunc) GetProcAddress(dll, symname);
#else
/*@-castfcnptr@*/
  data = (wrDLLFunc) dlsym(dll, symname);
/*@+castfcnptr@*/
#endif

  return data;
}


static int writerDLoaded(writerContext cntxt, const char *name)
{
  int i;

  for (i = 0; i < cntxt.aimWriterNum; i++)
    if (strcasecmp(name, cntxt.aimWriterName[i]) == 0) return i;

  return -1;
}


static int writerDYNload(writerContext *cntxt, const char *name)
{
  int i, len, ret;
  DLL dll;

  if (cntxt->aimWriterNum >= MAXWRITER) {
    printf(" Information: Number of Writers > %d!\n", MAXWRITER);
    return EGADS_INDEXERR;
  }
  dll = writerDLopen(name);
  if (dll == NULL) return EGADS_NULLOBJ;

  ret                      = cntxt->aimWriterNum;
  cntxt->aimExtension[ret] = (AIMext)    writerDLget(dll, "meshExtension");
  cntxt->aimWriter[ret]    = (AIMwriter) writerDLget(dll, "meshWrite");
  if (cntxt->aimExtension[ret] == NULL) {
    writerDLclose(dll);
    printf(" Error: Missing symbol 'meshExtension' in %s\n", name);
    return EGADS_EMPTY;
  }
  if (cntxt->aimWriter[ret] == NULL) {
    writerDLclose(dll);
    printf(" Error: Missing symbol 'meshWrite' in %s\n", name);
    return EGADS_EMPTY;
  }

  len = strlen(name) + 1;
  cntxt->aimWriterName[ret] = (char *) EG_alloc(len*sizeof(char));
  if (cntxt->aimWriterName[ret] == NULL) {
    writerDLclose(dll);
    return EGADS_MALLOC;
  }
  for (i = 0; i < len; i++) cntxt->aimWriterName[ret][i] = name[i];
  cntxt->aimWriterDLL[ret] = dll;
  cntxt->aimWriterNum++;

  return ret;
}


/*@null@*/ /*@dependent@*/ static const char *
aim_writerExtension(void *aimStruc, const char *writerName)
{
  int     i;
  aimInfo *aInfo;

  aInfo = (aimInfo *) aimStruc;
  i     = writerDLoaded(aInfo->wCntxt, writerName);
  if (i == -1) {
    i  = writerDYNload(&aInfo->wCntxt, writerName);
    if (i < 0) return NULL;
  }

  return aInfo->wCntxt.aimExtension[i]();
}


static int
aim_writeMesh(void *aimStruc, const char *writerName, aimMesh *mesh)
{
  int     i;
  aimInfo *aInfo;

  aInfo = (aimInfo *) aimStruc;
  i     = writerDLoaded(aInfo->wCntxt, writerName);
  if (i == -1) {
    i  = writerDYNload(&aInfo->wCntxt, writerName);
    if (i < 0) return i;
  }

  return aInfo->wCntxt.aimWriter[i](aimStruc, mesh);
}


/* ************************* Exposed Functions **************************** */

int
aim_deleteMeshes(void *aimStruc, aimMeshRef *meshRef)
{
  int           i;
#ifdef WIN32
  char          file[MAX_PATH];
#else
  char          file[PATH_MAX];
#endif
  const char    *ext;
  aimInfo       *aInfo;
  writerContext cntxt;

  aInfo = (aimInfo *) aimStruc;
  if (aInfo == NULL)                    return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  cntxt = aInfo->wCntxt;

  for (i = 0; i < cntxt.aimWriterNum; i++) {
    ext = aim_writerExtension(aimStruc, cntxt.aimWriterName[i]);
    if (ext == NULL) continue;
#ifdef WIN32
    _snprintf(file, MAX_PATH, "%s%s", meshRef->fileName, ext);
    _unlink(file);
#else
    snprintf(file, PATH_MAX, "%s%s", meshRef->fileName, ext);
    unlink(file);
#endif
  }

  return CAPS_SUCCESS;
}


int
aim_queryMeshes(void *aimStruc, int index, aimMeshRef *meshRef)
{
  int          i, j, k, ret, nWrite = 0;
  char         *writerName[MAXWRITER];
#ifdef WIN32
  char         file[MAX_PATH];
#else
  char         file[PATH_MAX];
#endif
  const char   *ext;
  aimInfo      *aInfo;
  capsObject   *vobject, *source, *last;
  capsProblem  *problem;
  capsAnalysis *analysis, *another;
  capsValue    *value;

  aInfo    = (aimInfo *) aimStruc;
  if (aInfo    == NULL)                 return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  if ((index < 1) || (index > analysis->nAnalysisOut)) return CAPS_BADINDEX;
  vobject = analysis->analysisOut[index-1];
  if (vobject == NULL)            return CAPS_NULLOBJ;
  value = (capsValue *) vobject->blind;
  if (value == NULL)              return CAPS_NULLBLIND;
  if (value->type != PointerMesh) return CAPS_BADTYPE;

  /* find all of our linked writers */
  for (i = 0; i < problem->nAnalysis; i++) {
    if (problem->analysis[i]        == NULL) continue;
    if (problem->analysis[i]->blind == NULL) continue;
    another = problem->analysis[i]->blind;
    if (analysis == another)                 continue;
    for (j = 0; j < another->nAnalysisIn; j++) {
      source = another->analysisIn[j];
      last   = NULL;
      do {
        if (source->magicnumber != CAPSMAGIC)      break;
        if (source->type        != VALUE)          break;
        if (source->blind       == NULL)           break;
        value  = (capsValue *) source->blind;
        if (value->link == another->analysisIn[j]) break;
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      if (last != vobject) continue;
      /* we hit our object from a AnalysisIn link */
      value = (capsValue *) another->analysisIn[j]->blind;
      if (value->meshWriter == NULL) {
        printf(" CAPS Warning: Link found but NULL writer!\n");
        continue;
      }
      for (k = 0; k < nWrite; k++)
        if (strcmp(writerName[k], value->meshWriter) == 0) break;
      if (k != nWrite) continue;
      writerName[nWrite] = value->meshWriter;
      nWrite++;
    }
  }
  if (nWrite == 0) return CAPS_NOTFOUND;

  /* look at the extensions for our files -- do they exist? */
  for (ret = i = 0; i < nWrite; i++) {
    ext = aim_writerExtension(aimStruc, writerName[i]);
    if (ext == NULL) continue;
#ifdef WIN32
    _snprintf(file, MAX_PATH, "%s%s", meshRef->fileName, ext);
    if (_access(file, F_OK) != 0) ret++;
#else
    snprintf(file,  PATH_MAX, "%s%s", meshRef->fileName, ext);
    if (access(file,  F_OK) != 0) ret++;
#endif
  }

  return ret;
}


int
aim_writeMeshes(void *aimStruc, int index, aimMesh *mesh)
{
  int          i, j, k, stat, nWrite = 0;
  char         *writerName[MAXWRITER];
#ifdef WIN32
  char         file[MAX_PATH];
#else
  char         file[PATH_MAX];
#endif
  const char   *ext;
  aimInfo      *aInfo;
  capsObject   *vobject, *source, *last;
  capsProblem  *problem;
  capsAnalysis *analysis, *another;
  capsValue    *value;

  aInfo    = (aimInfo *) aimStruc;
  if (aInfo == NULL)                    return CAPS_NULLOBJ;
  if (aInfo->magicnumber != CAPSMAGIC)  return CAPS_BADOBJECT;
  problem  = aInfo->problem;
  analysis = (capsAnalysis *) aInfo->analysis;

  if ((index < 1) || (index > analysis->nAnalysisOut)) return CAPS_BADINDEX;
  vobject = analysis->analysisOut[index-1];
  if (vobject == NULL)            return CAPS_NULLOBJ;
  value = (capsValue *) vobject->blind;
  if (value == NULL)              return CAPS_NULLBLIND;
  if (value->type != PointerMesh) return CAPS_BADTYPE;

  /* find all of our linked writers */
  for (i = 0; i < problem->nAnalysis; i++) {
    if (problem->analysis[i]        == NULL) continue;
    if (problem->analysis[i]->blind == NULL) continue;
    another = problem->analysis[i]->blind;
    if (analysis == another)                 continue;
    for (j = 0; j < another->nAnalysisIn; j++) {
      source = another->analysisIn[j];
      last   = NULL;
      do {
        if (source->magicnumber != CAPSMAGIC)      break;
        if (source->type        != VALUE)          break;
        if (source->blind       == NULL)           break;
        value  = (capsValue *) source->blind;
        if (value->link == another->analysisIn[j]) break;
        last   = source;
        source = value->link;
      } while (value->link != NULL);
      if (last != vobject) continue;
      /* we hit our object from a AnalysisIn link */
      value = (capsValue *) another->analysisIn[j]->blind;
#ifdef DEBUG
      printf(" *** Found Writer = %s ***\n", value->meshWriter);
#endif
      if (value->meshWriter == NULL) {
        printf(" CAPS Error: Link found but NULL writer!\n");
        return CAPS_NOTFOUND;
      }
      for (k = 0; k < nWrite; k++)
        if (strcmp(writerName[k], value->meshWriter) == 0) break;
      if (k != nWrite) continue;
      writerName[nWrite] = value->meshWriter;
      nWrite++;
    }
  }
  if (nWrite == 0) return CAPS_NOTFOUND;

  /* write the files */
  for (i = 0; i < nWrite; i++) {
    /* does the file exist? */
    ext = aim_writerExtension(aimStruc, writerName[i]);
    if (ext == NULL) continue;
#ifdef WIN32
    _snprintf(file, MAX_PATH, "%s%s", mesh->meshRef->fileName, ext);
    if (_access(file, F_OK) == 0) continue;
#else
    snprintf(file,  PATH_MAX, "%s%s", mesh->meshRef->fileName, ext);
    if (access(file,  F_OK) == 0) continue;
#endif
    /* no -- write it */
    stat = aim_writeMesh(aimStruc, writerName[i], mesh);
    if (stat != CAPS_SUCCESS) {
      printf(" CAPS Error: aim_writeMesh = %d for Writer = %s!\n",
             stat, writerName[i]);
      return stat;
    }
  }

  return CAPS_SUCCESS;
}


int
aim_initMeshRef(aimMeshRef *meshRef)
{
  if (meshRef == NULL) return CAPS_NULLOBJ;

  meshRef->nmap    = 0;
  meshRef->maps     = NULL;
  meshRef->fileName = NULL;

  return CAPS_SUCCESS;
}


int
aim_freeMeshRef(/*@null@*/ aimMeshRef *meshRef)
{
  int i;
  if (meshRef == NULL) return CAPS_NULLOBJ;

  if (meshRef->maps != NULL)
    for (i = 0; i < meshRef->nmap; i++)
      AIM_FREE(meshRef->maps[i].map);

  AIM_FREE(meshRef->maps);
  AIM_FREE(meshRef->fileName);

  aim_initMeshRef(meshRef);

  return CAPS_SUCCESS;
}


int
aim_initMeshData(aimMeshData *meshData)
{
  if (meshData == NULL) return CAPS_NULLOBJ;

  meshData->dim = 0;

  meshData->nVertex = 0;
  meshData->verts   = NULL;

  meshData->nElemGroup = 0;
  meshData->elemGroups = NULL;

  meshData->nTotalElems = 0;
  meshData->elemMap     = NULL;

  return CAPS_SUCCESS;
}


int
aim_freeMeshData(/*@null@*/ aimMeshData *meshData)
{
  int i;
  if (meshData == NULL) return CAPS_SUCCESS;

  AIM_FREE(meshData->verts);

  for (i = 0; i < meshData->nElemGroup; i++) {
    AIM_FREE(meshData->elemGroups[i].groupName);
    AIM_FREE(meshData->elemGroups[i].elements);
  }
  AIM_FREE(meshData->elemGroups);

  AIM_FREE(meshData->elemMap);

  aim_initMeshData(meshData);

  return CAPS_SUCCESS;
}


int
aim_addMeshElemGroup(void *aimStruc,
                     /*@null@*/ const char* groupName,
                     int ID,
                     enum aimMeshElem elementTopo,
                     int order,
                     int nPoint,
                     aimMeshData *meshData)
{
  int status = CAPS_SUCCESS;

  /* resize the element group */
  AIM_REALL(meshData->elemGroups, meshData->nElemGroup+1, aimMeshElemGroup, aimStruc, status);

  meshData->elemGroups[meshData->nElemGroup].groupName = NULL;
  meshData->elemGroups[meshData->nElemGroup].ID = ID;
  meshData->elemGroups[meshData->nElemGroup].elementTopo = elementTopo;
  meshData->elemGroups[meshData->nElemGroup].order = order;
  meshData->elemGroups[meshData->nElemGroup].nPoint = nPoint;
  meshData->elemGroups[meshData->nElemGroup].nElems = 0;
  meshData->elemGroups[meshData->nElemGroup].elements = NULL;

  meshData->nElemGroup++;

  if (groupName != NULL)
    AIM_STRDUP(meshData->elemGroups[meshData->nElemGroup].groupName, groupName, aimStruc, status);

cleanup:
  return status;
}


int
aim_addMeshElem(void *aimStruc,
                int nElems,
                aimMeshElemGroup *elemGroup)
{
  int status = CAPS_SUCCESS;

  /* resize the element group */
  AIM_REALL(elemGroup->elements, elemGroup->nPoint*(elemGroup->nElems + nElems), int, aimStruc, status);

  elemGroup->nElems += nElems;

cleanup:
  return status;
}


int
aim_readBinaryUgridHeader(void *aimStruc, aimMeshRef *meshRef,
                          int *nVertex, int *nTri, int *nQuad,
                          int *nTet, int *nPyramid, int *nPrism, int *nHex)
{
  int    status = CAPS_SUCCESS;

  char filename[PATH_MAX];
  FILE *fp = NULL;

  if (meshRef  == NULL) return CAPS_NULLOBJ;
  if (meshRef->fileName  == NULL) return CAPS_NULLOBJ;

  snprintf(filename, PATH_MAX, "%s%s", meshRef->fileName, ".lb8.ugrid");

  fp = fopen(filename, "rb");
  if (fp == NULL) {
    AIM_ERROR(aimStruc, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  /* read a binary UGRID file */
  status = fread(nVertex, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nTri,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nQuad,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nTet,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nPyramid, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nPrism,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(nHex,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

  status = CAPS_SUCCESS;

cleanup:
/*@-dependenttrans@*/
  if (fp != NULL) fclose(fp);
/*@+dependenttrans@*/
  fp = NULL;
  return status;
}


int
aim_readBinaryUgrid(void *aimStruc, aimMesh *mesh)
{
  int    status = CAPS_SUCCESS;

  int    nLine, nTri, nQuad;
  int    nTet, nPyramid, nPrism, nHex;
  int    i, j, elementIndex, nPoint, nElems, bcID, igroup;
  int    line[2], *mapGroupID=NULL, nMapGroupID=0;
  enum aimMeshElem elementTopo;
  char filename[PATH_MAX];
  FILE *fp = NULL, *fpID = NULL;
  aimMeshData *meshData = NULL;

  if (mesh           == NULL) return CAPS_NULLOBJ;
  if (mesh->meshRef  == NULL) return CAPS_NULLOBJ;
  if (mesh->meshRef->fileName  == NULL) return CAPS_NULLOBJ;

  status = aim_freeMeshData(mesh->meshData);
  AIM_STATUS(aimStruc, status);
  AIM_FREE(mesh->meshData);

  AIM_ALLOC(meshData, 1, aimMeshData, aimStruc, status);
  status = aim_initMeshData(meshData);
  AIM_STATUS(aimStruc, status);

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, ".lb8.ugrid");

  fp = fopen(filename, "rb");
  if (fp == NULL) {
    AIM_ERROR(aimStruc, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // Open a 2nd copy to read in the BC ID's
  fpID = fopen(filename, "rb");
  if (fpID == NULL) {
    AIM_ERROR(aimStruc, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  /* read a binary UGRID file */
  status = fread(&meshData->nVertex, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nTri,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nQuad,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nTet,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nPyramid, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nPrism,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
  status = fread(&nHex,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

  /* skip the header */
  status = fseek(fpID, 7*sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

  /*
  printf("\n Header from UGRID file: %d  %d %d  %d %d %d %d\n", numNode,
         numTriangle, numQuadrilateral, numTetrahedral, numPyramid, numPrism,
         numHexahedral);
   */

  AIM_ALLOC(meshData->verts, meshData->nVertex, aimMeshCoords, aimStruc, status);

  /* read all of the vertices */
  status = fread(meshData->verts, sizeof(aimMeshCoords), meshData->nVertex, fp);
  if (status != meshData->nVertex) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

  /* skip the verts */
  status = fseek(fpID, meshData->nVertex*sizeof(aimMeshCoords), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }


  // Numbers
  meshData->nTotalElems =  nTri     +
                           nQuad    +
                           nTet     +
                           nPyramid +
                           nPrism   +
                           nHex;

  // allocate the element map that maps back to the original element numbering
  AIM_ALLOC(meshData->elemMap, meshData->nTotalElems, aimMeshIndices, aimStruc, status);

  /*
  printf("Volume mesh:\n");
  printf("\tNumber of nodes          = %d\n", meshData->nVertex);
  printf("\tNumber of elements       = %d\n", meshData->nTotalElems);
  printf("\tNumber of triangles      = %d\n", numTriangle);
  printf("\tNumber of quadrilatarals = %d\n", numQuadrilateral);
  printf("\tNumber of tetrahedrals   = %d\n", numTetrahedral);
  printf("\tNumber of pyramids       = %d\n", numPyramid);
  printf("\tNumber of prisms         = %d\n", numPrism);
  printf("\tNumber of hexahedrals    = %d\n", numHexahedral);
*/

  /* skip the Tri+Quad elements for the ID */
  status = fseek(fpID, (3*nTri+4*nQuad)*sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; goto cleanup; }

  /* Start of element index */
  elementIndex = 0;

  /* Elements triangles */
  nPoint = 3;
  elementTopo = aimTri;
  for (i = 0; i < nTri; i++) {

    /* read the BC ID */
    status = fread(&bcID, sizeof(int), 1, fpID);
    if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
    if (bcID <= 0) {
      AIM_ERROR(aimStruc, "BC ID must be a positive number: %d!", bcID);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* add the BC group if necessary */
    if (bcID > nMapGroupID) {
      AIM_REALL(mapGroupID, bcID, int, aimStruc, status);
      for (j = nMapGroupID; j < bcID; j++) mapGroupID[j] = -1;
      nMapGroupID = bcID;
    }
    if (mapGroupID[bcID-1] == -1) {
      status = aim_addMeshElemGroup(aimStruc, NULL, bcID, elementTopo, 1, nPoint, meshData);
      AIM_STATUS(aimStruc, status);
      mapGroupID[bcID-1] = meshData->nElemGroup-1;
    }

    igroup = mapGroupID[bcID-1];

    /* add the element to the group */
    status = aim_addMeshElem(aimStruc, 1, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(&meshData->elemGroups[igroup].elements[nPoint*(meshData->elemGroups[igroup].nElems-1)],
                   sizeof(int), nPoint, fp);
    if (status != nPoint) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    meshData->elemMap[elementIndex][0] = igroup;
    meshData->elemMap[elementIndex][1] = meshData->elemGroups[igroup].nElems-1;

    elementIndex += 1;
  }

  /* Elements quadrilateral */
  nPoint = 4;
  elementTopo = aimQuad;
  for (i = 0; i < nQuad; i++) {

    /* read the BC ID */
    status = fread(&bcID, sizeof(int), 1, fpID);
    if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
    if (bcID <= 0) {
      AIM_ERROR(aimStruc, "BC ID must be a positive number: %d!", bcID);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* add the BC group if necessary */
    if (bcID > nMapGroupID) {
      AIM_REALL(mapGroupID, bcID, int, aimStruc, status);
      for (j = nMapGroupID; j < bcID; j++) mapGroupID[j] = -1;
      nMapGroupID = bcID;
    }
    if (mapGroupID[bcID-1] == -1) {
      status = aim_addMeshElemGroup(aimStruc, NULL, bcID, elementTopo, 1, nPoint, meshData);
      AIM_STATUS(aimStruc, status);
      mapGroupID[bcID-1] = meshData->nElemGroup-1;
    }

    igroup = mapGroupID[bcID-1];

    /* add the element to the group */
    status = aim_addMeshElem(aimStruc, 1, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(&meshData->elemGroups[igroup].elements[nPoint*(meshData->elemGroups[igroup].nElems-1)],
                   sizeof(int), nPoint, fp);
    if (status != nPoint) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    meshData->elemMap[elementIndex][0] = igroup;
    meshData->elemMap[elementIndex][1] = meshData->elemGroups[igroup].nElems-1;

    elementIndex += 1;
  }

  /*@-dependenttrans@*/
  fclose(fpID); fpID = NULL;
  /*@+dependenttrans@*/

  // skip face ID section of the file
  status = fseek(fp, (nTri + nQuad)*sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; goto cleanup; }

  // Elements Tetrahedral
  nPoint = 4;
  elementTopo = aimTet;
  nElems = nTet;

  if (nElems > 0) {
    status = aim_addMeshElemGroup(aimStruc, NULL, 0, elementTopo, 1, nPoint, meshData);
    AIM_STATUS(aimStruc, status);

    igroup = meshData->nElemGroup-1;

    /* add the elements to the group */
    status = aim_addMeshElem(aimStruc, nElems, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(meshData->elemGroups[igroup].elements,
                   sizeof(int), nPoint*nElems, fp);
    if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    for (i = 0; i < nElems; i++) {
      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = i;

      elementIndex += 1;
    }
  }

  // Elements Pyramid
  nPoint = 5;
  elementTopo = aimPyramid;
  nElems = nPyramid;

  if (nElems > 0) {
    status = aim_addMeshElemGroup(aimStruc, NULL, 0, elementTopo, 1, nPoint, meshData);
    AIM_STATUS(aimStruc, status);

    igroup = meshData->nElemGroup-1;

    /* add the elements to the group */
    status = aim_addMeshElem(aimStruc, nElems, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(meshData->elemGroups[igroup].elements,
                   sizeof(int), nPoint*nElems, fp);
    if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    for (i = 0; i < nElems; i++) {
      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = i;

      elementIndex += 1;
    }
  }


  // Elements Prism
  nPoint = 6;
  elementTopo = aimPrism;
  nElems = nPrism;

  if (nElems > 0) {
    status = aim_addMeshElemGroup(aimStruc, NULL, 0, elementTopo, 1, nPoint, meshData);
    AIM_STATUS(aimStruc, status);

    igroup = meshData->nElemGroup-1;

    /* add the elements to the group */
    status = aim_addMeshElem(aimStruc, nElems, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(meshData->elemGroups[igroup].elements,
                   sizeof(int), nPoint*nElems, fp);
    if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    for (i = 0; i < nElems; i++) {
      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = i;

      elementIndex += 1;
    }
  }


  // Elements Hex
  nPoint = 8;
  elementTopo = aimHex;
  nElems = nHex;

  if (nElems > 0) {
    status = aim_addMeshElemGroup(aimStruc, NULL, 0, elementTopo, 1, nPoint, meshData);
    AIM_STATUS(aimStruc, status);

    igroup = meshData->nElemGroup-1;

    /* add the elements to the group */
    status = aim_addMeshElem(aimStruc, nElems, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimStruc, status);

    /* read the element connectivity */
    status = fread(meshData->elemGroups[igroup].elements,
                   sizeof(int), nPoint*nElems, fp);
    if (status != nPoint*nElems) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    for (i = 0; i < nElems; i++) {
      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = i;

      elementIndex += 1;
    }
  }


  if (nTet+nPyramid+nPrism+nHex == 0) {
    meshData->dim = 2;

    status = fread(&nLine, sizeof(int), 1, fp);
    if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }

    meshData->nTotalElems += nLine;
    AIM_REALL(meshData->elemMap, meshData->nTotalElems, aimMeshIndices, aimStruc, status);

    // Elements Hex
    nPoint = 2;
    elementTopo = aimLine;
    nElems = nLine;

    for (i = 0; i < nElems; i++) {
      /* read the element connectivity */
      status = fread(line, sizeof(int), nPoint, fp);
      if (status != nPoint) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
      status = fread(&bcID, sizeof(int), 1, fp);
      if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimStruc, status); }
      if (bcID <= 0) {
        AIM_ERROR(aimStruc, "BC ID must be a positive number: %d!", bcID);
        status = CAPS_IOERR;
        goto cleanup;
      }

      /* add the BC group if necessary */
      if (bcID > nMapGroupID) {
        AIM_REALL(mapGroupID, bcID, int, aimStruc, status);
        for (j = nMapGroupID; j < bcID; j++) mapGroupID[j] = -1;
        nMapGroupID = bcID;
      }
      if (mapGroupID[bcID-1] == -1) {
        status = aim_addMeshElemGroup(aimStruc, NULL, bcID, elementTopo, 1, nPoint, meshData);
        AIM_STATUS(aimStruc, status);
        mapGroupID[bcID-1] = meshData->nElemGroup-1;
      }

      igroup = mapGroupID[bcID-1];

      /* add the elements to the group */
      status = aim_addMeshElem(aimStruc, 1, &meshData->elemGroups[igroup]);
      AIM_STATUS(aimStruc, status);

      meshData->elemGroups[igroup].elements[nPoint*meshData->elemGroups[igroup].nElems-2] = line[0];
      meshData->elemGroups[igroup].elements[nPoint*meshData->elemGroups[igroup].nElems-1] = line[1];

      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = i;

      elementIndex += 1;
    }

  } else {
    meshData->dim = 3;
  }

  mesh->meshData = meshData;
  meshData = NULL;

  status = CAPS_SUCCESS;

cleanup:
  if (status != CAPS_SUCCESS) {
    aim_freeMeshData(meshData);
    AIM_FREE(meshData);
  }
  AIM_FREE(mapGroupID);

/*@-dependenttrans@*/
  if (fp   != NULL) fclose(fp);
  if (fpID != NULL) fclose(fpID);
/*@+dependenttrans@*/

  return status;
}
