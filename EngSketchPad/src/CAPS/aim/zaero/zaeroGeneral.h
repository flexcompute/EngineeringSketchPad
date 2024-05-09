#include "zaeroTypes.h"
#include "cfdTypes.h"

// create (write cards for) aerodynamic wing component from feaAeroStruct
int zaero_aerodynamicWingComponent(void *aimInfo,
                                   FILE *fp,
                                   const feaAeroStruct *feaAero,
                                   const int compIndex,
                                   const feaFileTypeEnum formatType);

// create (write cards for) surface spline method from a feaAeroStruct
int zaero_splineMethod(void *aimInfo,
                       FILE *fp,
                       const feaAeroStruct *feaAero,
                       const zaeroSplineStruct *spline,
                       const feaFileTypeEnum formatType);

// create (write cards for) modeling physical data from a feaAeroRefStruct
int zaero_modelPhysicalData(void *aimInfo,
                            FILE *fp,
                            const zaeroHFGStruct *hfg,
                            const cfdUnitsStruct *units,
                            const feaFileTypeEnum formatType);

// create (write cards for) unsteady aerodynamics generator from zaeroUAICStruct
int zaero_unsteadyAerodynamicsGenerator(void *aimInfo,
                                        FILE *fp,
                                        const zaeroUAICStruct *uaic,
                                        const feaFileTypeEnum formatType);

/* finders */

// int zaero_findUAICByName(char *matchName, zaeroProblemStruct *zaeroProblem,
//                          zaeroUAICStruct **matchUAIC);
