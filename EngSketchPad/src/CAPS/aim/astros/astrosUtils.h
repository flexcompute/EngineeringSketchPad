// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "feaTypes.h"  // Bring in FEA structures
#include "vlmTypes.h"  // Bring in VLM structures

#ifdef __cplusplus
extern "C" {
#endif

// Write a Astros connections card from a feaConnection structure
int astros_writeConnectionCard(FILE *fp, feaConnectionStruct *feaConnect,
                               feaFileFormatStruct *feaFileFormat);

// Write out PLYLIST Card.
int astros_writePlyListCard(FILE *fp, feaFileFormatStruct *feaFileFormat,
                            int id, int numVal, int values[]);

// Write Nastran load card from a feaLoad structure
int astros_writeLoadCard(FILE *fp, meshStruct *mesh, feaLoadStruct *feaLoad,
                         feaFileFormatStruct *feaFileFormat);

// Write a Astros AEROS card from a feaAeroRef structure
int astros_writeAEROSCard(FILE *fp, feaAeroRefStruct *feaAeroRef,
                          feaFileFormatStruct *feaFileFormat);

// Write a Astros AEROS card from a feaAeroRef structure
int astros_writeAEROCard(FILE *fp, feaAeroRefStruct *feaAeroRef,
                         feaFileFormatStruct *feaFileFormat);

// Write out FLFact Card.
int astros_writeFLFactCard(FILE *fp, feaFileFormatStruct *feaFileFormat,
		                   int id, int numVal, double values[]) ;

// Write out AEFact Card.
int astros_writeAEFactCard(FILE *fp, feaFileFormatStruct *feaFileFormat,
		                   int id, int numVal, double values[]) ;

// Check to make for the bodies' topology are acceptable for airfoil shapes
// Return: CAPS_SUCCESS if everything is ok
// 		   CAPS_SOURCEERR if not acceptable
// 		   CAPS_* if something else went wrong
int astros_checkAirfoil(void *aimInfo,
						feaAeroStruct *feaAero);

// Write out all the Aero cards necessary to define the VLM surface
int astros_writeAeroData(void *aimInfo, FILE *fp,
			 int useAirfoilShape, // = true use the airfoils shape, = false panel
                         feaAeroStruct *feaAero,
                         feaFileFormatStruct *feaFileFormat);

// Write Astros CAERO6 cards from a feaAeroStruct
int astros_writeCAeroCard(FILE *fp, feaAeroStruct *feaAero,
                          feaFileFormatStruct *feaFileFormat);

// Write out all the Airfoil cards for each each of a surface
int astros_writeAirfoilCard(FILE *fp,
                            int useAirfoilShape, // = true use the airfoils shape, = false panel
                            feaAeroStruct *feaAero,
                            feaFileFormatStruct *feaFileFormat);

// Write Astros Spline1 cards from a feaAeroStruct
int astros_writeAeroSplineCard(FILE *fp, feaAeroStruct *feaAero,
                               feaFileFormatStruct *feaFileFormat);

// Write Astros constraint card from a feaConstraint structure
int astros_writeConstraintCard(FILE *fp, int feaConstraintSetID,
                               feaConstraintStruct *feaConstraint,
                               feaFileFormatStruct *feaFileFormat);

// Write Astros support card from a feaSupport structure
int astros_writeSupportCard(FILE *fp, feaSupportStruct *feaSupport,
                            feaFileFormatStruct *feaFileFormat);

// Write a Astros Property card from a feaProperty structure w/ design parameters
int astros_writePropertyCard(FILE *fp, feaPropertyStruct *feaProperty,
                             feaFileFormatStruct *feaFileFormat,
                             int numDesignVariable,
                             feaDesignVariableStruct feaDesignVariable[]);

// Write Astros element cards not supported by mesh_writeNastran in meshUtils.c
int astros_writeSubElementCard(FILE *fp, meshStruct *feaMesh, int numProperty,
                               feaPropertyStruct *feaProperty,
                               feaFileFormatStruct *feaFileFormat);

// Write a Astros Analysis card from a feaAnalysis structure
int astros_writeAnalysisCard(FILE *fp, feaAnalysisStruct *feaAnalysis,
                             feaFileFormatStruct *feaFileFormat);

// Write design variable/optimization information from a feaDesignVariable structure
int astros_writeDesignVariableCard(FILE *fp, feaDesignVariableStruct *feaDesignVariable,
                                   int numProperty, feaPropertyStruct feaProperty[],
                                   feaFileFormatStruct *feaFileFormat);

// Write design constraint/optimization information from a feaDesignConstraint structure
int astros_writeDesignConstraintCard(FILE *fp, int feaDesignConstraintSetID,
                                     feaDesignConstraintStruct *feaDesignConstraint,
                                     int numMaterial, feaMaterialStruct feaMaterial[],
                                     int numProperty, feaPropertyStruct feaProperty[],
                                     feaFileFormatStruct *feaFileFormat);

// Read data from a Astros OUT file to determine the number of eignevalues
int astros_readOUTNumEigenValue(FILE *fp, int *numEigenVector);

// Read data from a Astros OUT file to determine the number of grid points
int astros_readOUTNumGridPoint(FILE *fp, int *numGridPoint);

// Read data from a Astros OUT file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int astros_readOUTEigenVector(FILE *fp, int *numEigenVector, int *numGridPoint,
                              double ***dataMatrix);

// Read data from a Astros OUT file and load it into a dataMatrix[numEigenVector][5]
// where variables are eigenValue, eigenValue(radians), eigenValue(cycles),
// generalized mass, and generalized stiffness.      MASS              STIFFNESS
int astros_readOUTEigenValue(FILE *fp, int *numEigenVector, double ***dataMatrix);

// Read data from a Astros OUT file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int astros_readOUTDisplacement(FILE *fp, int subcaseId, int *numGridPoint,
                               double ***dataMatrix);

// Write geometric parametrization - only valid for modifications made by Bob Canfield to Astros
int astros_writeGeomParametrization(FILE *fp, void *aimInfo,
                                    int numDesignVariable,
                                    feaDesignVariableStruct *feaDesignVariable,
                                    int numGeomIn, capsValue *geomInVal,
                                    meshStruct *feaMesh,
                                    feaFileFormatStruct *feaFileFormat);

// Write a mesh contained in the mesh structure in Astros format (*.bdf)
int astros_writeMesh(void *aimInfo,
                     char *fname,
                     int asciiFlag, // 0 for binary, anything else for ascii
                     meshStruct *mesh,
                     feaFileTypeEnum gridFileType,
                     int numDesignVariable,
                     feaDesignVariableStruct feaDesignVariable[],
                     double scaleFactor);

#ifdef __cplusplus
}
#endif
