// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "capsTypes.h" // Bring in CAPS types

#ifdef __cplusplus
extern "C" {
#endif

// Read data from Mystran OUTPUT4 file and load it into a capsValue
int mystran_readOutput4Data(FILE *fp, const char *keyword, capsValue *val);

// Read data from a Mystran F06 file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int mystran_readF06EigenVector(FILE *fp, int *numEigenVector, int *numGridPoint, double ***dataMatrix);

// Read data from a Mystran F06 file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int mystran_readF06Displacement(FILE *fp, int subcaseId, int *numGridPoint, double ***dataMatrix);

#ifdef __cplusplus
}
#endif
