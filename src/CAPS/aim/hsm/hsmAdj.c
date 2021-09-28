#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "caps.h"

#include "meshUtils.h"    // Meshing utilities

#include "hsmTypes.h" // Bring in hsm structures

#define MAX(A,B)  (((A) < (B)) ? (B) : (A))

static void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// A function to implement bubble sort
static void bubbleSort(int arr[], int n)
{
   int i, j;
   for (i = 0; i < n-1; i++)

       // Last i elements are already in place
       for (j = 0; j < n-i-1; j++)
           if (arr[j] > arr[j+1])
              swap(&arr[j], &arr[j+1]);
}


// Generates the adjacency structure for HSM
int hsm_Adjacency(meshStruct *feaMesh,
                  const int numJoint, // Number of joints
                  const int *kjoint,  // Joint connectivity
                  int *maxAdjacency,    // Max valence
                  int **xadj_out,     // Pointers (indices) into adj (freeable)
                  int **adj_out) {    // The adjacency lists (freeable)

    int status; // Function return status

    int i, j, k, m, n;
    int *xadj = NULL, *adj = NULL;
    int **columns = NULL, *ncol = NULL;
    int numConnectivity;
    int found;

    if (xadj_out == NULL) {
        printf("Error xadj_out == NULL!\n");
        return CAPS_NULLVALUE;
    }

    if (adj_out == NULL) {
        printf("Error adj_out == NULL!\n");
        return CAPS_NULLVALUE;
    }

    *xadj_out = NULL;
    *adj_out  = NULL;

    ncol    = (int * ) EG_alloc(feaMesh->numNode*sizeof(int ));
    columns = (int **) EG_alloc(feaMesh->numNode*sizeof(int*));
    if ((ncol == NULL) || (columns == NULL)) {
        if (ncol    != NULL) EG_free(ncol);
        if (columns != NULL) EG_free(columns);
        return EGADS_MALLOC;
    }
    memset(ncol,    0, feaMesh->numNode*sizeof(int));
    memset(columns, 0, feaMesh->numNode*sizeof(int*));

    // Loop through the elements
    for (i = 0; i < feaMesh->numElement; i++) {

        // ignore line elements
        if (feaMesh->element[i].elementType == Line) continue;

        numConnectivity = mesh_numMeshElementConnectivity(&feaMesh->element[i]);
        if (numConnectivity < 0) { status = numConnectivity; goto cleanup; }

        // Loop through element connectivity
        for (j = 0; j < numConnectivity; j++) {

            // The row (0-based) for the current node
            n = feaMesh->element[i].connectivity[j] - 1;

            // check to see of the other nodes in the element should be added to the list
            for (k = 0; k < numConnectivity; k++) {
                found = (int) false;
                for (m = 0; m < ncol[n]; m++) {
                    if (columns[n][m] == feaMesh->element[i].connectivity[k] ) {
                        found = (int) true;
                        break;
                    }
                }
                if (found == (int) true) continue;

                columns[n] = (int *) EG_reall(columns[n], (ncol[n]+1)*sizeof(int));
                if (columns[n] == NULL) { status = EGADS_MALLOC; goto cleanup; }

                // add column with 1-based indexing
                columns[n][ncol[n]] = feaMesh->element[i].connectivity[k];
                ncol[n]++;
            }
        }
    }

    // collapse and wipe joint equations
    for (j = 0; j < numJoint; j++) {

      m = kjoint[2*j + 0];
      n = kjoint[2*j + 1];

      columns[n] = (int *) EG_reall(columns[n], (ncol[n]+ncol[m])*sizeof(int));
      if (columns[n] == NULL) { status = EGADS_MALLOC; goto cleanup; }

      // copy over the columns from m to n
      for (i = 0; i < ncol[m]; i++) {
          columns[n][ncol[n] + i] = columns[m][i];
      }
      //columns[n][ncol[n] + ncol[m]] = m; // add the node m to the list
      ncol[n] = ncol[n] + ncol[m];

      // wipe the other column as just a single connection to n
      columns[m] = (int *) EG_reall(columns[m], 2*sizeof(int));
      if (columns[m] == NULL) { status = EGADS_MALLOC; goto cleanup; }

      ncol[m] = 2;
      columns[m][0] = n+1;
      columns[m][1] = m+1;
    }

    // populate the pointer array so the adjacency can be allocated
    xadj = (int *) EG_alloc((feaMesh->numNode+1)*sizeof(int));
    if (xadj == NULL) { status = EGADS_MALLOC; goto cleanup; }
    memset(xadj, 0, (feaMesh->numNode+1)*sizeof(int));
    *maxAdjacency = 0;

    xadj[0] = 1; // 1-based indexing
    for (i = 0; i < feaMesh->numNode; i++) {
        xadj[i+1] = xadj[i] + ncol[i];
        bubbleSort(columns[i], ncol[i]);
        *maxAdjacency = MAX(*maxAdjacency, ncol[i]);
    }

    adj = (int *) EG_alloc(xadj[feaMesh->numNode]*sizeof(int));
    if (adj == NULL) { status = EGADS_MALLOC; goto cleanup; }
    memset(adj, 0, (xadj[feaMesh->numNode]-1)*sizeof(int));

    // populate the adjacency
    for (i = 0; i < feaMesh->numNode; i++) {
        for (j = 0; j < ncol[i]; j++) {
            adj[xadj[i]-1 + j] = columns[i][j];
        }
    }

//#define WRITE_MATRIX_MARKET
#ifdef WRITE_MATRIX_MARKET
    FILE *file = aim_fopen(aimInfo, "A.mtx", "w");

    //Write the banner
    fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(file, "%d %d %d\n", feaMesh->numNode, feaMesh->numNode, xadj[feaMesh->numNode]-1);

    //Write out the matrix data
    //Add one to the column index as the file format is 1-based
    for ( int row = 0; row < feaMesh->numNode; row++ )
      for ( int k = xadj[row]-1; k < xadj[row+1]-1; k++ )
        fprintf(file, "%d %d 100\n", row+1, adj[k]);

    fclose(file); file = NULL;
#endif

    status = CAPS_SUCCESS;

cleanup:
    for (i = 0; i < feaMesh->numNode; i++) {
        EG_free(columns[i]);
    }
    EG_free(columns);
    EG_free(ncol);

    if (status != CAPS_SUCCESS) {
        printf("Error: Premature exit in hsmAdj, status %d\n", status);
        EG_free(xadj);
        EG_free(adj);
    } else {
      *xadj_out = xadj;
      *adj_out = adj;
    }

    return status;
}

