#ifndef __ZAERO_UTILS_H__
#define __ZAERO_UTILS_H__

#include "zaeroTypes.h"

// general utility

int zaero_setBlankIntegerArray(int size, int **array);

int zaero_setBlankRealArray(int size, double **array);

int zaero_setBlankStringArray(int size, char ***array);

int zaero_allEqual(int type, void *a, void *b, ...);

// file utility

void zaero_readFile(FILE *fp, char **buffer);

// int  zaero_fileExists(char *filename);

int zaero_getFileLastModified(char *filePath, int *timestamp);

/* finders */

int zaero_findUAICByName(void *aimInfo,
                         const char *matchName,
                         const zaeroProblemStruct *zaeroProblem,
                         zaeroUAICStruct **matchUAIC);

int zaero_findSubcaseByName(void *aimInfo,
                            const char *matchName,
                            const zaeroProblemStruct *zaeroProblem,
                            zaeroSubcaseStruct **matchSubcase);

#endif // __ZAERO_UTILS_H__
