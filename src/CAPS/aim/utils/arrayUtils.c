// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
#include "capsTypes.h"  // Bring in CAPS types

// Allocate an integer matrix
int array_allocIntegerMatrix(int numRow, int numCol, int defaultValue, int ***matOut) {
    int status;
    int i, j;
    int **mat=NULL;

    *matOut = NULL;

    mat = (int **) EG_alloc(numRow*sizeof(int *));
    if (mat == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < numRow; i++) {
        mat[i] = (int *) EG_alloc(numCol*sizeof(int));

        if (mat[i] == NULL) {

            for (j = 0; j < i; j++) EG_free(mat[j]);

            status = EGADS_MALLOC;
            goto cleanup;
        }
    }

    for (i = 0; i < numRow; i++) {
        for(j = 0; j < numCol; j++) {
            mat[i][j] = defaultValue;
        }
    }

    *matOut = mat;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_allocIntegerMatrix, status = %d\n", status);
        return status;
}

// Free an integer matrix
int array_freeIntegerMatrix(int numRow, /*@unused@*/int numCol, int ***matOut) {
    int status;
    int i;
    int **mat=NULL;

    mat = *matOut;

    if (mat == NULL) return CAPS_SUCCESS;

    for (i = 0; i < numRow; i++) {
        if (mat[i] != NULL) EG_free(mat[i]);
    }

    if (mat != NULL) EG_free(mat);

    *matOut = NULL;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_freeIntegerMatrix, status = %d\n", status);
        return status;
}

// Allocate an double matrix
int array_allocDoubleMatrix(int numRow, int numCol, double defaultValue, double ***matOut) {
    int status;
    int i, j;
    double **mat=NULL;

    *matOut = NULL;

    mat = (double **) EG_alloc(numRow*sizeof(double *));
    if (mat == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < numRow; i++) {
        mat[i] = (double *) EG_alloc(numCol*sizeof(double));

        if (mat[i] == NULL) {

            for (j = 0; j < i; j++) EG_free(mat[j]);

            status = EGADS_MALLOC;
            goto cleanup;
        }
    }

    for (i = 0; i < numRow; i++) {
        for(j = 0; j < numCol; j++) {
            mat[i][j] = defaultValue;
        }
    }


    *matOut = mat;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_allocDoubleMatrix, status = %d\n", status);
        return status;
}

// Free an double matrix
int array_freeDoubleMatrix(int numRow, /*@unused@*/int numCol, double ***matOut) {
    int status;
    int i;
    double **mat=NULL;

    mat = *matOut;

    for (i = 0; i < numRow; i++) {
        if (mat[i] != NULL) EG_free(mat[i]);
    }

    if (mat != NULL) EG_free(mat);

    *matOut = NULL;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_freeDoubleMatrix, status = %d\n", status);
        return status;
}


// Set default integer array value
int array_setIntegerVectorValue(int numRow, int defaultValue, int **arrOut) {
    int status;
    int i;
    int *arr=NULL;

    if (*arrOut == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    arr = *arrOut;

    for (i = 0; i < numRow; i++) {
        arr[i] = defaultValue;
    }

    *arrOut = arr;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_setIntegerVectorValue, status = %d\n", status);
        return status;
}

// Set default double array value
int array_setDoubleVectorValue(int numRow, double defaultValue, double **arrOut) {
    int status;
    int i;
    double *arr=NULL;

    if (*arrOut == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    arr = *arrOut;

    for (i = 0; i < numRow; i++) {
        arr[i] = defaultValue;
    }

    *arrOut = arr;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_setDoubleVectorValue, status = %d\n", status);
        return status;
}

// Allocate an integer array
int array_allocIntegerVector(int numRow, int defaultValue, int **arrOut) {
    int status;
    int *arr=NULL;

    *arrOut = NULL;

    arr = (int *) EG_alloc(numRow*sizeof(int));
    if (arr == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    status = array_setIntegerVectorValue(numRow, defaultValue, &arr);
    if (status != CAPS_SUCCESS) goto cleanup;

    *arrOut = arr;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_allocIntegerArray, status = %d\n", status);
        return status;
}

// Allocate an double array
int array_allocDoubleVector(int numRow, double defaultValue, double **arrOut) {
    int status;
    double *arr=NULL;

    *arrOut = NULL;

    arr = (double *) EG_alloc(numRow*sizeof(double));
    if (arr == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    status = array_setDoubleVectorValue(numRow, defaultValue, &arr);
    if (status != CAPS_SUCCESS) goto cleanup;

    *arrOut = arr;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_allocDoubleArray, status = %d\n", status);
        return status;
}

// Max value in an double array
int array_maxDoubleValue(int numRow, double *arr, int *index, double *value) {
    int status;

    int i;

    int maxIndex;
    int maxValue;

    if (arr == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    if (numRow <= 0) {
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    maxIndex = 0;
    maxValue = arr[0];

    for (i = 1; i < numRow; i++) {
        if (maxValue >= arr[i]) continue;

        maxIndex = i;
        maxValue = arr[i];
    }

    *index = maxIndex;
    *value = maxValue;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in array_maxDoubleValue, status = %d\n", status);
        return status;
}
