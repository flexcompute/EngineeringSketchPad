#include <stdio.h>
#include <stdlib.h>


#define CHUNK    1024
#define MAGIC 9876543


  typedef struct {
    int vert2;            /* second vert */
    int index;            /* the stored index */
    int next;             /* next link */
  } vhEnt;


  typedef struct {
    int   magic;
    int   bias;
    int   nVert;
    int   *lookup;
    int   current;
    int   nentities;
    int   mentities;
    vhEnt *entities;
  } vhStruct;


int VH_init(int bias, int nVert, int startIndex, void **vhash)
{
  int      i;
  vhStruct *vh;
  
  *vhash = NULL;
  if ((bias != 0) && (bias != 1)) return 4;
  vh     = (vhStruct *) malloc(sizeof(vhStruct));
  if (vh == NULL) return 2;
  
  vh->magic     = MAGIC;
  vh->bias      = bias;
  vh->nVert     = nVert;
  vh->current   = startIndex;
  vh->nentities = 0;
  vh->mentities = CHUNK;
  vh->lookup    = (int *)   malloc(nVert*sizeof(int));
  vh->entities  = (vhEnt *) malloc(CHUNK*sizeof(vhEnt));
  if ((vh->lookup == NULL) || (vh->entities == NULL)) {
    if (vh->lookup   != NULL) free(vh->lookup);
    if (vh->entities != NULL) free(vh->entities);
    free(vh);
    return 2;
  }
  for (i = 0; i < nVert; i++) vh->lookup[i] = -1;
  
  *vhash = vh;
  return 0;
}


int VH_index(void *vhash, int index1, int index2, int *newIndex)
{
  int      i1, i2, n, link;
  vhEnt    *tmp;
  vhStruct *vh;
  
  *newIndex = -1;
  if (vhash == NULL) return 3;
  
  vh = (vhStruct *) vhash;
  if (vh->magic != MAGIC) return 3;
  if (index1 == index2)   return 4;
  if (index1 <  index2) {
    i1 = index1 - vh->bias;
    i2 = index2 - vh->bias;
  } else {
    i1 = index2 - vh->bias;
    i2 = index1 - vh->bias;
  }
  if ((i1 < 0) || (i1 >= vh->nVert) || (i2 < 0) || (i2 >= vh->nVert)) return 3;
  link = vh->lookup[i1];
  
  /* nothing in lookup table -- must be new */
  if (link == -1) {
    if (vh->nentities == vh->mentities) {
      n   = vh->mentities + CHUNK;
      tmp = (vhEnt *) realloc(vh->entities, n*sizeof(vhEnt));
      if (tmp == NULL) return 2;
      vh->entities  = tmp;
      vh->mentities = n;
    }
    vh->entities[vh->nentities].vert2 = i2;
    vh->entities[vh->nentities].index = vh->current;
    vh->entities[vh->nentities].next  = -1;
    vh->lookup[i1]                    = vh->nentities;
    vh->nentities++;
    
    *newIndex = vh->current;
    vh->current++;
    return 1;
  }
  
  /* follow the links */
  for (;;) {
    if (vh->entities[link].vert2 == i2) {
      /* found our entity */
      *newIndex = vh->entities[link].index;
      return 0;
    }
    if (vh->entities[link].next == -1) break;
    link = vh->entities[link].next;
  }
  
  /* not found -- add */
  if (vh->nentities == vh->mentities) {
    n   = vh->mentities + CHUNK;
    tmp = (vhEnt *) realloc(vh->entities, n*sizeof(vhEnt));
    if (tmp == NULL) return 2;
    vh->entities  = tmp;
    vh->mentities = n;
  }
  vh->entities[vh->nentities].vert2 = i2;
  vh->entities[vh->nentities].index = vh->current;
  vh->entities[vh->nentities].next  = -1;
  vh->entities[link].next           = vh->nentities;
  vh->nentities++;
  
  *newIndex = vh->current;
  vh->current++;
  return 1;
}


int VH_populate(void *vhash, int index1, int index2, int index3)
{
  int      i1, i2, i3, n, link;
  vhEnt    *tmp;
  vhStruct *vh;

//  *newIndex = -1;
  if (vhash == NULL) return 3;

  vh = (vhStruct *) vhash;
  if (vh->magic != MAGIC) return 3;
  if (index1 == index2)   return 4;
  if (index1 <  index2) {
    i1 = index1 - vh->bias;
    i2 = index2 - vh->bias;
  } else {
    i1 = index2 - vh->bias;
    i2 = index1 - vh->bias;
  }
  i3 = index3; // - vh->bias;
  if ((i1 < 0) || (i1 >= vh->nVert) || (i2 < 0) || (i2 >= vh->nVert)) return 3;
  link = vh->lookup[i1];

  /* nothing in lookup table -- must be new */
  if (link == -1) {
    if (vh->nentities == vh->mentities) {
      n   = vh->mentities + CHUNK;
      tmp = (vhEnt *) realloc(vh->entities, n*sizeof(vhEnt));
      if (tmp == NULL) return 2;
      vh->entities  = tmp;
      vh->mentities = n;
    }
    vh->entities[vh->nentities].vert2 = i2;
    vh->entities[vh->nentities].index = i3; //vh->current;
    vh->entities[vh->nentities].next  = -1;
    vh->lookup[i1]                    = vh->nentities;
    vh->nentities++;

//    *newIndex = vh->current;
//    vh->current++;
    return 1;
  }

  /* follow the links */
  for (;;) {
    if (vh->entities[link].vert2 == i2) {
      /* found our entity */
//      *newIndex = vh->entities[link].index;
      return 0;
    }
    if (vh->entities[link].next == -1) break;
    link = vh->entities[link].next;
  }

  /* not found -- add */
  if (vh->nentities == vh->mentities) {
    n   = vh->mentities + CHUNK;
    tmp = (vhEnt *) realloc(vh->entities, n*sizeof(vhEnt));
    if (tmp == NULL) return 2;
    vh->entities  = tmp;
    vh->mentities = n;
  }
  vh->entities[vh->nentities].vert2 = i2;
  vh->entities[vh->nentities].index = i3; //vh->current;
  vh->entities[vh->nentities].next  = -1;
  vh->entities[link].next           = vh->nentities;
  vh->nentities++;

//  *newIndex = vh->current;
//  vh->current++;
  return 1;
}


void VH_free(void **vhash)
{
  vhStruct *vh;
  
  if (*vhash == NULL) return;

  vh = (vhStruct *) *vhash;
  if (vh->magic != MAGIC) return;
#ifdef STANDALONE
  printf(" VH_free with %d/%d entities!\n", vh->nentities, vh->mentities);
#endif
  free(vh->lookup);
  free(vh->entities);
  free(vh);
  
  *vhash = NULL;
}


#ifdef STANDALONE
int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int  i, j, k, index, stat;
  void *ptr;
  
  /* initialize */
  stat = VH_init(1, 256, 0, &ptr);
  if (stat != 0) {
    printf(" Error: VH_init = %d\n", stat);
    return 1;
  }
  
  /* fill'er up */
  for (i = 0; i < 255; i++)
    for (j = i+1; j < 256; j++) {
      stat = VH_index(ptr, i+1, j+1, &index);
      if (stat != 1) {
        printf(" Error: 1 VH_index%d %d = %d\n", i+1, j+1, stat);
        return 1;
      }
    }
  
  /* check that we can retrieve what we have filled */
  for (k = i = 0; i < 255; i++)
    for (j = i+1; j < 256; j++, k++) {
      stat = VH_index(ptr, j+1, i+1, &index);
      if (stat != 0) {
        printf(" Error: 2 VH_index %d %d = %d\n", i+1, j+1, stat);
        return 1;
      }
      if (index != k) printf(" Warning: %d != %d\n", index, k);
    }
  
  /* cleanup */
  VH_free(&ptr);
  
  return 0;
}
#endif
