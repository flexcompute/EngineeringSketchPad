// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "feaTypes.h"  // Bring in FEA structures
#include "vlmTypes.h"  // Bring in VLM structures

#ifdef __cplusplus
extern "C" {
#endif


// Write SET case control card
int nastran_writeSetCard(FILE *fp, int n, int numSetID, int *setID);

// Write a Nastran element cards not supported by mesh_writeNastran in meshUtils.c
int nastran_writeSubElementCard(FILE *fp, meshStruct *feaMesh, int numProperty, feaPropertyStruct *feaProperty, feaFileFormatStruct *feaFileFormat);

// Write a Nastran connections card from a feaConnection structure
int nastran_writeConnectionCard(FILE *fp, feaConnectionStruct *feaConnect, feaFileFormatStruct *feaFileFormat);

// Write a Nastran AERO card from a feaAeroRef structure
int nastran_writeAEROCard(FILE *fp, feaAeroRefStruct *feaAeroRef, feaFileFormatStruct *feaFileFormat);

// Write a Nastran AEROS card from a feaAeroRef structure
int nastran_writeAEROSCard(FILE *fp, feaAeroRefStruct *feaAeroRef, feaFileFormatStruct *feaFileFormat);

// Write Nastran SET1 card from a feaAeroStruct
int nastran_writeSet1Card(FILE *fp, feaAeroStruct *feaAero, feaFileFormatStruct *feaFileFormat);

// Write Nastran Spline1 cards from a feaAeroStruct
int nastran_writeAeroSplineCard(FILE *fp, feaAeroStruct *feaAero, feaFileFormatStruct *feaFileFormat);

// Write Nastran CAERO1 cards from a feaAeroStruct
int nastran_writeCAeroCard(FILE *fp, feaAeroStruct *feaAero, feaFileFormatStruct *feaFileFormat);

// Write Nastran coordinate system card from a feaCoordSystemStruct structure
int nastran_writeCoordinateSystemCard(FILE *fp, feaCoordSystemStruct *feaCoordSystem, feaFileFormatStruct *feaFileFormat);

// Write combined Nastran constraint card from a set of constraint IDs.
// 	The combined constraint ID is set through the constraintID variable.
int nastran_writeConstraintADDCard(FILE *fp, int constraintID, int numSetID, int constraintSetID[], feaFileFormatStruct *feaFileFormat);

// Write Nastran constraint card from a feaConstraint structure
int nastran_writeConstraintCard(FILE *fp, feaConstraintStruct *feaConstraint, feaFileFormatStruct *feaFileFormat);

// Write Nastran support card from a feaSupport structure - withID = NULL or false SUPORT, if withID = true SUPORT1
int nastran_writeSupportCard(FILE *fp, feaSupportStruct *feaSupport, feaFileFormatStruct *feaFileFormat, int *withID);

// Write a Nastran Material card from a feaMaterial structure
int nastran_writeMaterialCard(FILE *fp, feaMaterialStruct *feaMaterial, feaFileFormatStruct *feaFileFormat);

// Write a Nastran Property card from a feaProperty structure
int nastran_writePropertyCard(FILE *fp, feaPropertyStruct *feaProperty, feaFileFormatStruct *feaFileFormat);

// Write a combined Nastran load card from a set of load IDs. Uses the feaLoad structure to get the local scale factor
// 	for the load. The overall load scale factor is set to 1. The combined load ID is set through the loadID variable.
int nastran_writeLoadADDCard(FILE *fp, int loadID, int numSetID, int loadSetID[], feaLoadStruct feaLoad[], feaFileFormatStruct *feaFileFormat);

// Write Nastran load card from a feaLoad structure
int nastran_writeLoadCard(FILE *fp, feaLoadStruct *feaLoad, feaFileFormatStruct *feaFileFormat);

// Write Nastran analysis card from a feaAnalysis structure
int nastran_writeAnalysisCard(FILE *fp, feaAnalysisStruct *feaAnalysis, feaFileFormatStruct *feaFileFormat);

// Write a combined Nastran design constraint card from a set of constraint IDs.
//  The combined design constraint ID is set through the constraint ID variable.
int nastran_writeDesignConstraintADDCard(FILE *fp, int constraintID, int numSetID, int designConstraintSetID[], feaFileFormatStruct *feaFileFormat);

// Write design constraint/optimization information from a feaDesignConstraint structure
int nastran_writeDesignConstraintCard(FILE *fp, feaDesignConstraintStruct *feaDesignConstraint, feaFileFormatStruct *feaFileFormat);

// Write design variable/optimization information from a feaDesignVariable structure
int nastran_writeDesignVariableCard(FILE *fp, feaDesignVariableStruct *feaDesignVariable, feaFileFormatStruct *feaFileFormat);

// Write design variable relation information from a feaDesignVariableRelation structure
int nastran_writeDesignVariableRelationCard(void *aimInfo, FILE *fp, feaDesignVariableRelationStruct *feaDesignVariableRelation, feaProblemStruct *feaProblem, feaFileFormatStruct *feaFileFormat);

// Write equation information from a feaDesignEquation structure
int nastran_writeDesignEquationCard(FILE *fp, feaDesignEquationStruct *feaEquation, feaFileFormatStruct *fileFormat);

// Write design table constants information from a feaDesignTable structure 
int nastran_writeDesignTableCard(FILE *fp, feaDesignTableStruct *feaDesignTable, feaFileFormatStruct *fileFormat);

// Write design response information from a feaDesignResponse structure 
int nastran_writeDesignResponseCard(FILE *fp, feaDesignResponseStruct *feaDesignResponse, feaFileFormatStruct *fileFormat);

// Write design equation response information from a feaDesignEquationResponse structure 
int nastran_writeDesignEquationResponseCard(FILE *fp, feaDesignEquationResponseStruct *feaDesignEquationResponse, feaProblemStruct *feaProblem, feaFileFormatStruct *fileFormat);

// Write design optimization parameter information from a feaDesignOptParam struct
int nastran_writeDesignOptParamCard(FILE *fp, feaDesignOptParamStruct *feaDesignOptParam, feaFileFormatStruct *fileFormat);

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
int nastran_readF06Displacement(FILE *fp, int subcaseId, int *numGridPoint, double ***dataMatrix);

// Read objective values for a Nastran OP2 file  and liad it into a dataMatrix[numPoint]
int nastran_readOP2Objective(char *filename, int *numPoint,  double **dataMatrix);

int nastran_writeAeroCamberTwist(FILE *fp, int numAero, feaAeroStruct *feaAero, feaFileFormatStruct *feaFileFormat);

#ifdef __cplusplus
}
#endif
