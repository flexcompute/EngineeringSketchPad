#include "feaTypes.h"  // Bring in FEA structures

#ifdef __cplusplus
extern "C" {
#endif

typedef double DOUBLE_2[2];
typedef double DOUBLE_6[6];

// Write Abaqus analysis card from a feaAnalysis structure
int abaqus_writeAnalysisCard(void *aimInfo,
                             FILE *fp,
                             const int numLoad,
                  /*@null@*/ const feaLoadStruct feaLoad[],
                             const feaAnalysisStruct *feaAnalysis,
                             const meshStruct *mesh);

// Write Abaqus constraint card from a feaConstraint structure
int abaqus_writeConstraintCard(FILE *fp, const feaConstraintStruct *feaConstraint);

// Write a Abaqus Material card from a feaMaterial structure
int abaqus_writeMaterialCard(FILE *fp, const feaMaterialStruct *feaMaterial);

// Write a Abaqus Property card from a feaProperty structure
int abaqus_writePropertyCard(FILE *fp, const feaPropertyStruct *feaProperty);

// Write Abaqus load card from a feaLoad structure
int abaqus_writeLoadCard(FILE *fp, const feaLoadStruct *feaLoad, const meshStruct *mesh);

// Write a Abaqus Element Set card
int abaqus_writeElementSet(void *aimInfo, FILE *fp, char *setName, int numElement, int elementSet[], int numName, /*@null@*/ char *name[]);

// Write a Abaqus Node Set card
int abaqus_writeNodeSet(void* aimInfo, FILE *fp, char *setName, int numNode, int nodeSet[], int numName, /*@null@*/ char *name[]);

// Write element and node sets for various components of the problem
int abaqus_writeAllSets(void* aimInfo, FILE *fp, const feaProblemStruct *feaProblem);

// Read displacement information from .dat file
int abaqus_readDATDisplacement(void *aimInfo, const char *filename,
                               int numGridPoint, int **nodeID, DOUBLE_6 **dataMatrix);

// Read mises information from .dat file
int abaqus_readDATMises(void *aimInfo, const char *filename,
                        int numElement, int **elemID, double **elemData);

int abaqus_readFIL(void *aimInfo, char *filename, int field, int *numData,  double ***dataMatrix);

#ifdef __cplusplus
}
#endif
