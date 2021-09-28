#ifndef HASH_ELEMENT_H
#define HASH_ELEMENT_H

#define NVERT 4

typedef struct {
  int verts[NVERT];
} hashElem;

typedef struct {
  int* elemIndex;
  int nelem;
} hashElemLink;

typedef struct {

  // table of nvertex long that points into elements array
  int nvertex;
  hashElemLink* lookup;

  // array of elements with sorted vertexes
  hashElem* elements;

} hashElemTable;

int initiate_hashTable(hashElemTable *table);
int allocate_hashTable(int nvertex, int nelem, hashElemTable *table);
int destroy_hashTable(hashElemTable *table);
int hash_addElement(int nvertex, int vertex[], int elemIndex, hashElemTable *table);
int hash_getIndex(int nvertex, int vertex[], hashElemTable *table, int *elemIndex);

#endif //HASH_ELEMENT_H
