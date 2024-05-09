#ifndef __ZAERO_DISCIPLINE_H__
#define __ZAERO_DISCIPLINE_H__

#include "zaeroTypes.h"

// LinearFlutter ///////////////////////////////////////////////////////////////

/* define subcase discipline with LinearFlutter analysis input */
int zaero_getLinearFlutterDiscipline(void *aimInfo,
                                     char *jsonDict,
                                     zaeroSubcaseStruct *subcase);

/* write LinearFlutter analysis cards from defined subcase */
int zaero_flutterAnalysis(void *aimInfo,
                          FILE *fp,
                          const zaeroSubcaseStruct *subcase,
                          const zaeroProblemStruct *zaeroProblem,
                          const cfdUnitsStruct *units,
                          const feaFileTypeEnum formatType);

// StaticAeroelastic / Trim ////////////////////////////////////////////////////

/* define subcase discipline with Trim analysis input */
int zaero_getTrimDiscipline(void *aimInfo,
                            const char *jsonDict,
                            const cfdUnitsStruct *units,
                            zaeroSubcaseStruct *subcase);

/* define Trim analysis variables from input tuples*/
int zaero_getTrimVariables(void *aimInfo,
                           int numtrimVarTuple,
                           capsTuple trimVarTuple[],
                           zaeroProblemStruct *zaeroProblem);

/* write StaticAeroelastic / Trim analysis cards from defined subcase */
int zaero_trimAnalysis(void *aimInfo,
                       FILE *fp,
                       const zaeroSubcaseStruct *subcase,
                       const zaeroProblemStruct *zaeroProblem,
                       const feaFileTypeEnum formatType);

#endif // __ZAERO_DISCIPLINE_H__
