#include <stdlib.h>
#ifdef STANDALONE
#include <stdio.h>
#endif


#define VBLOCK 2048             /* allocation chunk */

  typedef struct {
    int vert;                   /* other vert index */
    int thread;                 /* next vert in list */
  } valence;



static void
fillValence(int i0, int i1, int *vti, int *vtable, valence *vthread)
{
  int thread, last;
  
  /* have we seen this vert? */
  thread = vtable[i0];
  last   = -1;
  while (thread != -1) {
    if (vthread[thread].vert == i1) return;
    last   = thread;
    thread = vthread[thread].thread;
  }
  
  /* new entry */
  if (last == -1) {
    vtable[i0] = *vti;
  } else {
    vthread[last].thread = *vti;
  }
  vthread[*vti].vert   = i1;
  vthread[*vti].thread = -1;
  *vti += 1;
}


int
maxValence(int nvert, int ntri, int *tris)
{
  int     i, j, i0, i1, i2, vtlen, vti, thread, cnt, maxv = 0;
  int     *vtable;
  valence *vthread, *vtmp;
  static int combo[3][3] = {{0,1,2}, {1,0,2}, {2,0,1}};

  vtable = (int *) malloc(nvert*sizeof(int));
  if (vtable == NULL) return -1;
  for (i = 0; i < nvert; i++) vtable[i] = -1;

  vtlen   = VBLOCK;
  vthread = (valence *) malloc(vtlen*sizeof(valence));
  if (vthread == NULL) {
    free(vtable);
    return -2;
  }

  /* fill the valence table */
  vti = 0;
  for (i = 0; i < ntri; i++)
    for (j = 0; j < 3; j++) {
      if (vtlen < vti+2) {
        vtmp = (valence *) realloc(vthread,(vtlen+VBLOCK)*sizeof(valence));
        if (vtmp == NULL) {
          free(vthread);
          free(vtable);
          return -3;
        }
        vthread  = vtmp;
        vtlen   += VBLOCK;
      }
      i0 = tris[3*i+combo[j][0]]-1;
      i1 = tris[3*i+combo[j][1]]-1;
      i2 = tris[3*i+combo[j][2]]-1;
      fillValence(i0, i1, &vti, vtable, vthread);
      fillValence(i0, i2, &vti, vtable, vthread);
    }
  
  /* count the max number */
  for (i = 0; i < nvert; i++) {
    thread = vtable[i];
    cnt    = 0;
    while (thread != -1) {
      cnt++;
      thread = vthread[thread].thread;
    }
    if (cnt > maxv) maxv = cnt;
  }
#ifdef STANDALONE
  printf(" table length = %d\n", vti);
#endif
  
  free(vthread);
  free(vtable);
  return maxv;
}


#ifdef STANDALONE
int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  static int triangles[3*32] = {
     1,  2,  6,
     7,  6,  2,
     2,  3,  7,
     8,  7,  3,
     3,  4,  8,
     9,  8,  4,
     4,  5,  9,
    10,  9,  5,
     6,  7, 11,
    12, 11,  7,
     7,  8, 12,
    13, 12,  8,
     8,  9, 13,
    14, 13,  9,
     9, 10, 14,
    15, 14, 10,
    11, 12, 16,
    17, 16, 12,
    12, 13, 17,
    18, 17, 13,
    13, 14, 18,
    19, 18, 14,
    14, 15, 19,
    20, 19, 15,
    16, 17, 21,
    22, 21, 17,
    17, 18, 22,
    23, 22, 18,
    18, 19, 23,
    24, 23, 19,
    19, 20, 24,
    25, 24, 20 };
  
  printf(" max  Valence = %d\n", maxValence(25, 32, triangles));
  return 0;
}
#endif
