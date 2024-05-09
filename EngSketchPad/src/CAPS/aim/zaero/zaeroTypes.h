#ifndef __ZAERO_TYPES_H__
#define __ZAERO_TYPES_H__

#include "capsTypes.h"
#include "feaTypes.h"

typedef enum {
    UnknownFlutterFunction,
    FixMach,
    FixMachDensity,
    FixMachAtmos,
    FixAltAtmos
} zaeroFlutterFuncEnum;

typedef enum {
    UnknownDiscipline,
    LinearFlutter,
    ParamFlutter,
    Aeroservoelastic,
    StaticAeroelastic,
    ManeuverLoads,
    EjectionLoads,
    GustLoads,
    MFTGustLoads,
    NonLinearFlutter
} zaeroDisciplineTypeEnum;

// HFG Module - Aerodynamic Model Input
typedef struct {

    char *XZSymmetry; // TODO: deduct from inherited feaAeroRef if possible ?

    char *flip;

    double refChord;
    double refSpan;
    double refArea;

    double refCenter[3];

} zaeroHFGStruct;

// FEM Module - Modal Data Importer
typedef struct {

    /* Filename of the external file that contains the free vibration solution
     * of the structural finite model
     */
    char *filename;

    /* Name of the structural finite model code that generates the output file
     * containing the free vibration solution of the structural finite model
     * Options: ['MSC', 'UAI', 'CSA', 'NE', 'ASTROS', 'IDEAS', 'ELFINI',
     *           'GENESIS', 'ABAQUS', 'FREE']
     * Refer to the manual for more information on those options
     */
    char *form;

    /* Boundary condition of the structural finite element model
     */
    char *boundary;

    /* Specifies the degrees of freedom of the rigid bodymodes of the structural
     * finite element model
     */
    int suport[2];

    /* Print options
     */
    int printFlag;

    /* Whether to use the A-set (Analysis Set) degrees of freedom of the
     * finite element model
     */
    int asetFlag;

} zaeroFEMStruct;

// UAIC Module - Unsteady Aerodynamics Data Generation (AIC) Matrices
typedef struct {

    // Name identifying this UAIC configuration
    char *name;

    // ID identifying this UAIC configuration
    int id;

    // Mach number
    double machNumber;

    /* Aerodynamic method flag
     * METHOD = 0    for the ZONA6/ZONA7/ZSAP method
     * METHOD = 1    for the ZTAIC method
     * METHOD +/- 2  for the ZONA7U method
     * METHOD = 3    for the ZTRAN method
     */
    int methodFlag;

    // Filename to store/load generated AIC matrices (max 16 char)
    char aicFilename[17];

    // Reduced frequecies
    int numReducedFreq;
    double *reducedFreq;

    /* Whether to SAVE or ACQUIRE AIC matrices
     * - 0: SAVE
     * - 1: ACQUIRE
     */
    int saveFlag;

    // Print flag
    int printFlag;

} zaeroUAICStruct;

// SPLINE Module - Aerodynamic & FEM Model Interconnection
typedef struct {

    /* spline method, can be one of:
     *   - 0: zero displacement
     *   - 1: surface spline
     *   - 2: beam spline
     *   - 3: 3-D spline
     */
    int method;

    /* multiplication factor to obtain a small tolerance to detect any
     * duplicated location of structural grid points. Tolerance is computed
     * by eps x refChord
     */
    double eps;

    /* linear attachment flexibility */
    double attachFlex;

} zaeroSplineStruct;

typedef struct {

    /* flutter analysis function, can be one of:
     *   - FixMach (1)       : non-matched point flutter analysis at a fixed
     *                         Mach number for various velocity and density
     *                         pairs
     *   - FixMachDensity (2): non-matched point flutter analysis at a fixed
     *                         Mach number and density pair for various
     *                         velocities
     *   - FixMachAtmos (3)  : matched point flutter analysis at a fixed Mach
     *                         number and an atmosphere table for various
     *                         altitudes
     *   - FixAltAtmos (4)   : matched point flutter analysis at a fixed
     *                         altitude and an atmosphere table for various
     *                         Mach numbers
     */
    zaeroFlutterFuncEnum function;

    int numModes;

    char *boundaryCondition;

    // FixMachDensity variables
    double density;

    double refVelocity;

    int numVelocities;
    double *velocities;

    // modal damping variables
    int numDampingFreq;
    double *dampingFreq;

    int numDampingValues;
    double *dampingValues;

    char *dampingUnits;

    // print flags
    int printFlag;

} zaeroLinearFlutterStruct;

typedef struct {

    // trim variable unique name
    char *name;

    // label
    char *label;

    /* value
     *   - "none"        : (only applicable to degrees of freedom associated
                            with translational and angular acceleration)
                            the variable is eliminated from the trim system
     *   - "free"        : the variable is a free degree of freedom. The value
     *                     is unknown and to be solved by the trim system
     *   - "<Real Value>": the variable (constraint) is fixed and given by the
     *                     real value
     */
    char *value;

    // limits
    double lowerLimit;
    double upperLimit;

    // zaeroTrimLink link;

    char *dmi;

    // char *sym;

    // initial guess of the trim variable for the minimization computation
    //   of an over-determined trim
    double initial;

    // // user input aerodynamic stability derivatives
    // char *aeroStabilityDerivatives;

} zaeroTrimVarStruct;

typedef struct {

    double dynamicPressure;

    // zaeroTrimObjectiveStruct objective;

    // zaeroTrimConstraintStruct constraint;

    double vectorToCG[3];

    double gravityAcceleration;

    double weight;

    double weightMomentOfInertia[6];

    char *accUnit;

    // trim variables
    int numVariables;
    char **variables;

} zaeroTrimStruct;

typedef struct {

    // analysis unique name
    char *name;

    // Defines a subtitle of each subcase section by a character string.
    // (up to 72 characters)
    char *subtitle;

    // Defines subcase id
    int subcaseID;

    // Provides additional description of the subcase by a character string.
    // (up to 72 characters)
    char *label;

    // the unique id for this analysis
    int analysisID;

    // type of analysis (currently NOT USED)
    // TODO: merge zaeroDisciplineTypes to analysisType?
    analysisTypeEnum analysisType;

    // discipline
    zaeroDisciplineTypeEnum disciplineType;

    void *discipline;

    // name of UAIC used by this analysis subcase
    char *uaicName;

} zaeroSubcaseStruct;

typedef struct {

    // zaero input file
    char *input;

    // zaero output file
    char *output;

    // // fem input files
    // int numFEM;
    // capsTuple *fem;

    // aic filenames
    int numAIC;
    capsTuple *aic;

} zaeroArtifactsStruct;

typedef struct {

    // Describes the job by a character string. (up to 72 characters)
    char *title;

    // Controls echo (printout) of the Bulk Data Section
    char *echo;

    // HFG module data
    zaeroHFGStruct hfg;

    // FEM module data
    // Array of FEMs for each free vibration solution to import (at most 2)
    int numFEMs;
    zaeroFEMStruct FEMs[2];

    // UAIC module data
    // Array of UAICs for each Mach Number - Reduced Frequency pair defined
    int numUAICs;
    zaeroUAICStruct *UAICs;

    // FEA Aerodynamic information
    int numAero;
    feaAeroStruct *feaAero; // size = [numAero]

    // SPLINE module data
    zaeroSplineStruct spline;

    // Array of defined Trim variables
    int numTrimVariables;
    zaeroTrimVarStruct *trimVariables;

    // Array of defined subcases
    int numSubcases;
    zaeroSubcaseStruct *subcases;


} zaeroProblemStruct;


// "constructors" / "destructors"

int initiate_zaeroLinearFlutterStruct(zaeroLinearFlutterStruct *flutter);

int destroy_zaeroLinearFlutterStruct(zaeroLinearFlutterStruct *flutter);

int initiate_zaeroTrimVarStruct(zaeroTrimVarStruct *variable);

int destroy_zaeroTrimVarStruct(zaeroTrimVarStruct *variable);

int initiate_zaeroTrimStruct(zaeroTrimStruct *trim);

int destroy_zaeroTrimStruct(zaeroTrimStruct *trim);

int initiate_zaeroHFGStruct(zaeroHFGStruct *hfg);

int destroy_zaeroHFGStruct(zaeroHFGStruct *hfg);

int initiate_zaeroFEMStruct(zaeroFEMStruct *fem);

int destroy_zaeroFEMStruct(zaeroFEMStruct *fem);

int initiate_zaeroUAICStruct(zaeroUAICStruct *uaic);

int destroy_zaeroUAICStruct(zaeroUAICStruct *uaic);

int initiate_zaeroSplineStruct(zaeroSplineStruct *spline);

int destroy_zaeroSplineStruct(zaeroSplineStruct *spline);

int initiate_zaeroSubcaseStruct(zaeroSubcaseStruct *zaeroSubcase);

int destroy_zaeroSubcaseStruct(zaeroSubcaseStruct *zaeroSubcase);

int initiate_zaeroProblemStruct(zaeroProblemStruct *zaeroProblem);

int destroy_zaeroProblemStruct(zaeroProblemStruct *zaeroProblem);

int initiate_zaeroArtifactsStruct(zaeroArtifactsStruct *artifacts);

int destroy_zaeroArtifactsStruct(zaeroArtifactsStruct *artifacts);


#endif // __ZAERO_TYPES_H__
