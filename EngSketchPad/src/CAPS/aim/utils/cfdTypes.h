// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Structures for general CFD analysis - Written by Dr. Ryan Durscher AFRL/RQVC

#ifndef CFDTYPES_H
#define CFDTYPES_H

typedef enum {UnknownBoundary, Inviscid, Viscous, Farfield, Extrapolate, Freestream,
              BackPressure, Symmetry, SubsonicInflow, SubsonicOutflow,
              MassflowIn, MassflowOut, FixedInflow, FixedOutflow, MachOutflow} cfdSurfaceTypeEnum;

typedef enum {DesignVariableUnknown, DesignVariableGeometry, DesignVariableAnalysis} cfdDesignVariableTypeEnum;

// Structure to hold CFD surface information
typedef struct {
    char   *name;

    cfdSurfaceTypeEnum surfaceType; // "Global" boundary condition types

    int    bcID;                 // ID of boundary

    // Wall specific properties
    int    wallTemperatureFlag;  // Temperature flag
    double wallTemperature;      // Temperature value -1 = adiabatic ; >0 = isothermal
    double wallHeatFlux;         // Wall heat flux. to use Temperature flag should be true and wallTemperature < 0

    // Symmetry plane
    int symmetryPlane;        // Symmetry flag / plane

    // Stagnation quantities
    double totalPressure;    // Total pressure
    double totalTemperature; // Total temperature
    double totalDensity;     // Total density

    // Static quantities
    double staticPressure;   // Static pressure
    double staticTemperature;// Static temperature
    double staticDensity;    // Static temperature

    // Velocity components
    double uVelocity;  // x-component of velocity
    double vVelocity;  // y-component of velocity
    double wVelocity;  // z-component of velocity
    double machNumber; // Mach number

    // Massflow
    double massflow; // Mass flow through a boundary

} cfdSurfaceStruct;

// Structure to hold CFD boundary condition information
typedef struct {
    char             *name;   // Name of BCsStruct
    int              numSurfaceProp;       // Number of unique BC ids
    cfdSurfaceStruct *surfaceProp;  // Surface properties for each bc - length of numSurfaceProp

} cfdBoundaryConditionStruct;

// Structure to hold EigenValue information as used in CFD solvers
typedef struct {

    char *name;

    int modeNumber;

    double frequency;
    double damping;

    double generalMass;
    double generalDisplacement;
    double generalVelocity;
    double generalForce;

} cfdEigenValueStruct;

// Structure to hold a collection of EigenValue information as used in CFD solvers
typedef struct {
    int surfaceID;

    int numEigenValue;
    cfdEigenValueStruct *eigenValue;

    double freestreamVelocity;
    double freestreamDynamicPressure;
    double lengthScaling;

} cfdModalAeroelasticStruct;

// Structure to hold design variable information as used in CFD solvers
typedef struct {
    char *name;

    cfdDesignVariableTypeEnum type; // variable type

    const capsValue *var;

    double *value;        // value of the variable [length]
    double *lowerBound;   // Lower bounds of variable [length]
    double *upperBound;   // Upper bounds of variable [length]
    double *typicalSize;  // typical size of variable [length]

} cfdDesignVariableStruct;

// Structure to hold a component of an output functional used in CFD solvers
typedef struct {

    char *name;

    double target;
    double weight;
    double power;
    double bias;

    int frame;
    int form;

    int bcID;
    char *boundaryName;

} cfdDesignFunctionalCompStruct;

// Structure to hold output functional information as used in CFD solvers
typedef struct {

    char *name;

    int numComponent;
    cfdDesignFunctionalCompStruct *component;

    double value;                  // computed objective function value

    int numDesignVariable;
    cfdDesignVariableStruct *dvar; // dObjective/dDesignVariable: numDesignVariable in length

} cfdDesignFunctionalStruct;

// Structure to hold a collection of optimization information as used in CFD solvers
typedef struct {

    int numDesignFunctional;
    cfdDesignFunctionalStruct *designFunctional; // [numObjective]

    int numDesignVariable;
    cfdDesignVariableStruct *designVariable; // [numDesignVariable]

} cfdDesignStruct;


// Structure to hold CFD unit system
typedef struct {

    // base units
    char *length;      // Length unit
    char *mass;        // mass unit
    char *time;        // time unit
    char *temperature; // temperature unit

    // derived units
    char *density;      // density unit
    char *pressure;     // pressure unit
    char *speed;        // speed unit
    char *acceleration; // acceleration unit
    char *force;        // force unit
    char *viscosity;    // viscosity unit
    char *area;         // area unit

    // coefficient units
    char *Cpressure; // Pressure Coefficient unit
    char *Cforce;    // Force Coefficient unit
    char *Cmoment;   // Moment Coefficient unit

} cfdUnitsStruct;

#endif
