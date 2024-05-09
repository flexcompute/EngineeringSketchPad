#ifndef __ZAERO_GRAPHICAL_H__
#define __ZAERO_GRAPHICAL_H__

#include "zaeroTypes.h"

int zaero_textFileGenerationAero(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                                 const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID);

int zaero_textFileGenerationFlutter(void *aimInfo,
                                    FILE *fp,
                                    const char *jsonDict,
                                    const zaeroProblemStruct *zaeroProblem,
                                    const feaFileTypeEnum formatType,
                                    int *setID);

int zaero_textFileGenerationMode(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                                 const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID);

int zaero_textFileGenerationTrim(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                                 const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID);

int zaero_textFileGenerationVG(void *aimInfo,
                               FILE *fp,
                               const char *jsonDict,
                               const zaeroProblemStruct *zaeroProblem,
                               const feaFileTypeEnum formatType,
                               int *setID);

#endif // __ZAERO_GRAPHICAL_H__
