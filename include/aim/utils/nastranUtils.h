// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_NASTRANUTILS_H_
#define _AIM_UTILS_NASTRANUTILS_H_

#include "feaTypes.h"  // Bring in FEA structures
#include "vlmTypes.h"  // Bring in VLM structures

#ifdef __cplusplus
extern "C" {
#endif

// Scaling factors for Nastran VLM IDs
#define NASTRAN_SURFID 10000

// Write SET case control card
int nastran_writeSetCard(FILE *fp, int n, int numSetID, int *setID);

// Write a Nastran element cards not supported by mesh_writeNastran in meshUtils.c
int nastran_writeSubElementCard(void *aimInfo, FILE *fp, const meshStruct *feaMesh, int numProperty, const feaPropertyStruct *feaProperty, const feaFileFormatStruct *feaFileFormat);

// Write a Nastran connections card from a feaConnection structure
int nastran_writeConnectionCard(FILE *fp, const feaConnectionStruct *feaConnect, const feaFileFormatStruct *feaFileFormat);

// Write a Nastran AERO card from a feaAeroRef structure
int nastran_writeAEROCard(FILE *fp, const feaAeroRefStruct *feaAeroRef, const feaFileFormatStruct *feaFileFormat);

// Write a Nastran AEROS card from a feaAeroRef structure
int nastran_writeAEROSCard(FILE *fp, const feaAeroRefStruct *feaAeroRef, const feaFileFormatStruct *feaFileFormat);

// Write Nastran Spline1 cards from a feaAeroStruct
int nastran_writeAeroSplineCard(void *aimInfo, FILE *fp, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat);

// Write Nastran AESURF cards for controls from feaAeroStructs
int nastran_writeControls(void *aimInfo, FILE *fp, const int numAero, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat);

// Write Nastran CAERO1 cards from a feaAeroStruct
int nastran_writeCAeroCard(void *aimInfo, FILE *fp, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat);

// Write Nastran coordinate system card from a feaCoordSystemStruct structure
int nastran_writeCoordinateSystemCard(FILE *fp, const feaCoordSystemStruct *feaCoordSystem, const feaFileFormatStruct *feaFileFormat);

// Write combined Nastran constraint card from a set of constraint IDs.
// 	The combined constraint ID is set through the constraintID variable.
int nastran_writeConstraintADDCard(FILE *fp, int constraintID, int numSetID, int constraintSetID[], const feaFileFormatStruct *feaFileFormat);

// Write Nastran constraint card from a feaConstraint structure
int nastran_writeConstraintCard(FILE *fp, const feaConstraintStruct *feaConstraint, const feaFileFormatStruct *feaFileFormat);

// Write Nastran support card from a feaSupport structure - withID = NULL or false SUPORT, if withID = true SUPORT1
int nastran_writeSupportCard(FILE *fp, const feaSupportStruct *feaSupport, const feaFileFormatStruct *feaFileFormat, const int *withID);

// Write a Nastran Material card from a feaMaterial structure
int nastran_writeMaterialCard(FILE *fp, const feaMaterialStruct *feaMaterial, const feaFileFormatStruct *feaFileFormat);

// Write a Nastran Property card from a feaProperty structure
int nastran_writePropertyCard(FILE *fp, const feaPropertyStruct *feaProperty, const feaFileFormatStruct *feaFileFormat);

// Write a combined Nastran load card from a set of load IDs. Uses the feaLoad structure to get the local scale factor
// 	for the load. The overall load scale factor is set to 1. The combined load ID is set through the loadID variable.
int nastran_writeLoadADDCard(FILE *fp, int loadID, int numSetID, int loadSetID[], feaLoadStruct feaLoad[], const feaFileFormatStruct *feaFileFormat);

// Write Nastran load card from a feaLoad structure
int nastran_writeLoadCard(FILE *fp, const feaLoadStruct *feaLoad, const feaFileFormatStruct *feaFileFormat);

// Write Nastran analysis card from a feaAnalysis structure
int nastran_writeAnalysisCard(FILE *fp, feaAnalysisStruct *feaAnalysis, const feaFileFormatStruct *feaFileFormat);

// Write a combined Nastran design constraint card from a set of constraint IDs.
//  The combined design constraint ID is set through the constraint ID variable.
int nastran_writeDesignConstraintADDCard(FILE *fp, int constraintID, int numSetID, const int designConstraintSetID[], const feaFileFormatStruct *feaFileFormat);

// Write design constraint/optimization information from a feaDesignConstraint structure
int nastran_writeDesignConstraintCard(FILE *fp, const feaDesignConstraintStruct *feaDesignConstraint, const feaProblemStruct *feaProblem, const feaFileFormatStruct *feaFileFormat);

// Write design variable/optimization information from a feaDesignVariable structure
int nastran_writeDesignVariableCard(void *aimInfo, FILE *fp, const feaDesignVariableStruct *feaDesignVariable, const feaFileFormatStruct *feaFileFormat);

// Write design variable relation information from a feaDesignVariableRelation structure
int nastran_writeDesignVariableRelationCard(void *aimInfo, FILE *fp, const feaDesignVariableRelationStruct *feaDesignVariableRelation, const feaProblemStruct *feaProblem, const feaFileFormatStruct *feaFileFormat);

// Write equation information from a feaDesignEquation structure
int nastran_writeDesignEquationCard(FILE *fp, const feaDesignEquationStruct *feaEquation, const feaFileFormatStruct *fileFormat);

// Write design table constants information from a feaDesignTable structure
int nastran_writeDesignTableCard(FILE *fp, const feaDesignTableStruct *feaDesignTable, const feaFileFormatStruct *fileFormat);

// Write design response information from a feaDesignResponse structure
int nastran_writeDesignResponseCard(FILE *fp, const feaDesignResponseStruct *feaDesignResponse, const feaFileFormatStruct *fileFormat);

// Write design equation response information from a feaDesignEquationResponse structure
int nastran_writeDesignEquationResponseCard(FILE *fp, const feaDesignEquationResponseStruct *feaDesignEquationResponse, const feaProblemStruct *feaProblem, const feaFileFormatStruct *fileFormat);

// Write design optimization parameter information from a feaDesignOptParam struct
int nastran_writeDesignOptParamCard(FILE *fp, const feaDesignOptParamStruct *feaDesignOptParam, const feaFileFormatStruct *fileFormat);

// Write mass increment information from a feaMassIncrement structure
int nastran_writeMasssetCard(FILE *fp, int incrementID, const feaFileFormatStruct *feaFileFormat);

// Write mass increment information from a feaMassIncrement structure
int nastran_writeMassIncrementSet(FILE *fp, int incrementID, int numIncrementID, const int massIncrementSetID[], const feaMassIncrementStruct *feaMassIncrement, const feaFileFormatStruct *feaFileFormat);

// Read data from a Nastran F06 file to determine the number of eignevalues
int nastran_readF06NumEigenValue(FILE *fp, int *numEigenVector);

// Read data from a Nastran F06 file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int nastran_readF06EigenVector(FILE *fp, int *numEigenVector, int *numGridPoint, double ***dataMatrix);

// Read data from a Nastran F06 file and load it into a dataMatrix[numEigenVector][5]
// where variables are eigenValue, eigenValue(radians), eigenValue(cycles), generalized mass, and generalized stiffness.                                                                   MASS              STIFFNESS
int nastran_readF06EigenValue(FILE *fp, int *numEigenVector, double ***dataMatrix);

// Read data from a Nastran F06 file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int nastran_readF06Displacement(void *aimInfo, FILE *fp, int subcaseId, int *numGridPoint, double ***dataMatrix);

// Read weight value from a Nastran F06 file
int nastran_readF06GridPointWeightGeneratorOutput(FILE *fp, double *mass, double cg[3],
                                                  double is[6], double iq[3], double q[9]);

int nastran_writeAeroCamberTwist(void *aimInfo, FILE *fp, int numAero, feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat);

#ifdef __cplusplus
}
#endif

#endif // _AIM_UTILS_NASTRANUTILS_H_
