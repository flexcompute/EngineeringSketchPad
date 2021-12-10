// structure and routines to make element vertexes to an element array index

#include "capsTypes.h"

#include "hashElement.h"

static int initiate_hashElem(hashElem* elem)
{
  int i;

  for (i = 0; i < NVERT; i++)
    elem->verts[i] = -1;

  return CAPS_SUCCESS;
}

static int initiate_hashElemLink(hashElemLink* link)
{

  link->elemIndex = NULL;
  link->nelem = 0;

  return CAPS_SUCCESS;
}

int initiate_hashTable(hashElemTable *table)
{

  table->nvertex = 0;
  table->lookup = NULL;
  table->elements = NULL;

  return CAPS_SUCCESS;
}

int destroy_hashTable(hashElemTable *table)
{
  int i;

  for (i = 0; i < table->nvertex; i++)
    EG_free(table->lookup[i].elemIndex);

  EG_free(table->lookup);
  EG_free(table->elements);

  table->nvertex = 0;
  table->lookup = NULL;
  table->elements = NULL;

  return CAPS_SUCCESS;
}

int allocate_hashTable(int nvertex, int nelem, hashElemTable *table)
{
  int status = CAPS_SUCCESS;
  int i;

  destroy_hashTable(table);

  table->lookup = (hashElemLink *)EG_alloc((nvertex+1)*sizeof(hashElemLink));
  if (table->lookup == NULL) { status = EGADS_MALLOC; goto cleanup; }

  table->nvertex = nvertex+1;
  for (i = 0; i <= nvertex; i++)
    initiate_hashElemLink(table->lookup+i);

  table->elements = (hashElem *)EG_alloc(nelem*sizeof(hashElem));
  if (table->elements == NULL) { status = EGADS_MALLOC; goto cleanup; }

  for (i = 0; i < nelem; i++)
    initiate_hashElem(table->elements+i);

  status = CAPS_SUCCESS;

  cleanup:
    if (status != CAPS_SUCCESS)
      printf("Premature exit in initiate_hashElem status = %d\n", status);

    if (status != CAPS_SUCCESS) {
      EG_free(table->lookup); table->lookup = NULL;
      EG_free(table->elements); table->elements = NULL;
    }
    return status;
}

static void swapi(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// A function to implement bubble sort
static void bubbleSort(int n, int vert[])
{
  int i, j;
  
  for (i = 0; i < n-1; i++)
    // Last i elements are already in place
    for (j = 0; j < n-i-1; j++)
      if (vert[j] > vert[j+1]) {
        swapi(&vert[j], &vert[j+1]);
      }
}

// vertex indexing is assumed to be 1-based
int hash_addElement(int nvertex, int vertex[], int elemIndex,
                    hashElemTable *table)
{

  int status = CAPS_SUCCESS;
  int i;
  int verts[NVERT] = {0};
  hashElemLink* hlink;

  if (nvertex > NVERT+1) return CAPS_BADRANK;

  for (i = 0; i < nvertex; i++)
    verts[i] = vertex[i];

  // sort the vertexes
  bubbleSort(nvertex, verts);

  // add the element to the lookup table
  hlink = table->lookup + verts[0];

  // increase the number of elements
  hlink->nelem++;
  hlink->elemIndex = (int *)EG_reall(hlink->elemIndex, hlink->nelem*sizeof(int));
  if (hlink->elemIndex == NULL) { status = EGADS_MALLOC; goto cleanup; }

  // save the element index and vertexes
  hlink->elemIndex[hlink->nelem-1] = elemIndex;
  for (i = 0; i < nvertex; i++)
    table->elements[elemIndex].verts[i] = verts[i];

  status = CAPS_SUCCESS;
  cleanup:
    if (status != CAPS_SUCCESS) printf("Premature exit in addElement_hashTable status = %d\n", status);

    return status;
}

// get an element index based on the vertexes of the element
int hash_getIndex(int nvertex, int vertex[], hashElemTable *table,
                  int *elemIndex)
{

  int status = CAPS_SUCCESS;
  int i;
  int found = (int)false;
  int elemInd, ielem;
  int verts[NVERT] = {0};
  hashElemLink* hlink;

  if (nvertex > NVERT+1) return CAPS_BADRANK;

  for (i = 0; i < nvertex; i++)
    verts[i] = vertex[i];

  // sort the vertexes
  bubbleSort(nvertex, verts);

  // add the element to the lookup table
  if (table->lookup[verts[0]].nelem == 0) {
    printf("Error: Vertex %d is not in the hash table\n", verts[0]);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // look for the element in the pool
  hlink = table->lookup+verts[0];
  for (ielem = 0; ielem < hlink->nelem; ielem++) {

    elemInd = hlink->elemIndex[ielem];

    found = (int)true;
    for (i = 0; i < nvertex; i++)
      found = (int)(found && (table->elements[elemInd].verts[i] == verts[i]));

    if (found == (int)true) break;
  }

  if (found == (int)false) {
    printf("Could not find element in hash table!!\n");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // found the element index
  *elemIndex = elemInd;

  status = CAPS_SUCCESS;

  cleanup:
    if (status != CAPS_SUCCESS)
      printf("Premature exit in getIndex_hashTable status = %d\n", status);

    return status;
}

