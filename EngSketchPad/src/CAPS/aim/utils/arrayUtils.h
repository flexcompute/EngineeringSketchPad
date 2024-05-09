// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_ARRAYUTILS_H_
#define _AIM_UTILS_ARRAYUTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

// Allocate an integer matrix
int array_allocIntegerMatrix(int numRow, int numCol, int defaultValue, int ***matOut);

// Free an integer matrix
int array_freeIntegerMatrix(int numRow, int numCol, int ***matOut);

// Allocate an double matrix
int array_allocDoubleMatrix(int numRow, int numCol, double defaultValue, double ***matOut);

// Free an double matrix
int array_freeDoubleMatrix(int numRow, int numCol, double ***matOut);

// Set default integer array value
int array_setIntegerVectorValue(int numRow, int defaultValue, int **arrOut);

// Set default double array value
int array_setDoubleVectorValue(int numRow, double defaultValue, double **arrOut);

// Allocate an integer array
int array_allocIntegerVector(int numRow, int defaultValue, int **arrOut);

// Allocate an double array
int array_allocDoubleVector(int numRow, double defaultValue, double **arrOut);

// Max value in an double array
int array_maxDoubleValue(int numRow, double *arr, int *index, double *value);

// Remove duplicates in a integer array - if in2 == NULL, in1 is simply copied
int array_removeIntegerDuplicate(int numIn1, int *in1, int numIn2, /*@null@*/ int *in2, int *numOut, int **out);

#ifdef __cplusplus
}
#endif

#endif //_AIM_UTILS_ARRAYUTILS_H_
