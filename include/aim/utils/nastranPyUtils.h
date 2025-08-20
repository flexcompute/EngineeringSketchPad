// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_NASTRANPYUTILS_H_
#define _AIM_UTILS_NASTRANPYUTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

// General container for to map attribute names to an assigned index
typedef struct {

  // Flag to see if Python was already initialized - i.e. the AIM was called from within Python
  int initPy;

  // The PyObject op2 model
  void *pymodel;

} PyOP2model;

// Open and read nastran OP2 file
int nastran_openOP2(void *aimInfo, const char *filename, PyOP2model *model);

// Close the model
int nastran_closeOP2(PyOP2model *model);

// Read displacement values for a Nastran OP2 file and load it into a dataMatrix
int nastran_readOP2Displacement(void *aimInfo, PyOP2model *model,
                                int itime, int subcaseId, int *numGridPoint, double ***dataMatrix);

// Read eigen values for a Nastran OP2 file and load it into a dataMatrix
int nastran_readOP2EigenValue(void *aimInfo, PyOP2model *model,
                              int subcaseId, int *numGridPoint, double ***dataMatrix);

// Read weight value from a Nastran OP2 file
int nastran_readOP2GridPointWeightGeneratorOutput(void *aimInfo, PyOP2model *model,
                                                  double *mass, double cg[3],
                                                  double is[6], double iq[3], double q[9]);

// Read objective values for a Nastran OP2 file  and liad it into a dataMatrix[numPoint]
int nastran_readOP2Objective(void *aimInfo, PyOP2model *model, int *numPoint,  double **dataMatrix);

#ifdef __cplusplus
}
#endif

#endif // _AIM_UTILS_NASTRANUTILS_H_
