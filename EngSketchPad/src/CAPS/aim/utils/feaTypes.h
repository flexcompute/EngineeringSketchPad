// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Structures for general FEA analysis - Written by Dr. Ryan Durscher AFRL/RQVC

#ifndef FEATYPES_H
#define FEATYPES_H

#include "meshTypes.h" // Bring in mesh structures
#include "vlmTypes.h"  // Bring in vortex lattice methods structures

                  // Linear ----------------------------------- Linear
typedef enum {UnknownMaterial, Isotropic, Anisothotropic, Orthotropic, Anisotropic} materialTypeEnum;

typedef enum {UnknownProperty, ConcentratedMass, Rod, Bar, Beam, Shear, Shell, Composite, Solid, Membrane} propertyTypeEnum;

typedef enum {UnknownConstraint, Displacement, ZeroDisplacement} constraintTypeEnum;

typedef enum {UnknownLoad, GridForce, GridMoment, LineForce, LineMoment, Gravity, Pressure, PressureDistribute, Rotational, Thermal, PressureExternal} loadTypeEnum;

typedef enum {UnknownAnalysis, Modal, Static, Optimization, AeroelasticTrim, AeroelasticFlutter} analysisTypeEnum;

typedef enum {UnknownFileType, SmallField, LargeField, FreeField} feaFileTypeEnum;

typedef enum {UnknownDesignVar, MaterialDesignVar, PropertyDesignVar, ElementDesignVar} feaDesignVariableTypeEnum;

typedef enum {UnknownCoordSystem, RectangularCoordSystem, SphericalCoordSystem, CylindricalCoordSystem } feaCoordSystemTypeEnum;

typedef enum {UnknownConnection, Mass, Spring, Damper, RigidBody, RigidBodyInterpolate} feaConnectionTypeEnum;

// Indexing of massIntertia in feaPropertyStruct
typedef enum { I11, I21, I22, I31, I32, I33 } feaMassInertia;


// Structure to hold aerodynamic reference information
typedef struct {

    int coordSystemID; // Aerodynamic coordinate sytem id
    int rigidMotionCoordSystemID; // Reference coordinate system identification for rigid body motions.

    double refChord; // Reference chord length.      Reference span.  (Real > 0.0)
    double refSpan; // Reference span
    double refArea; // Reference area

    int symmetryXZ; // Symmetry key for the aero coordinate x-z plane.  (Integer = +1 for symmetry, 0 for no symmetry,
                      // and -1 for antisymmetry; Default = 0)
    int symmetryXY; // The symmetry key for the aero coordinate x-y plane can be used to simulate ground effects.
                    // (Integer = +1 for antisymmetry, 0 for no symmetry, and -1 for symmetry; Default = 0)
} feaAeroRefStruct;

// Structure to hold aerodynamic information
typedef struct {

    char *name; // Aero name

    int surfaceID; // Surface ID
    int coordSystemID; // Coordinate system ID

    int numGridID; // Number of grid IDs in grid ID set for the spline
    int *gridIDSet; // List of grid IDs to apply spline to. size = [numGridID]

    vlmSurfaceStruct vlmSurface; // VLM surface structure

} feaAeroStruct;

// Structure to hold connection information
typedef struct {
    char *name; // Connection name

    int connectionID; // Connection ID

    feaConnectionTypeEnum connectionType; // Connection type

    int elementID;

    // RBE2 - dependent degrees of freedom
    int connectivity[2]; // Grid IDs - 0 index = Independent grid ID, 1 index = Dependent grid ID
    int dofDependent;

    // Spring (scalar)
    double stiffnessConst;
    double dampingConst;
    double stressCoeff;
    int    componentNumberStart;
    int    componentNumberEnd;

    // Damper (scalar) - see spring for additional entries

    // Mass (scalar) - see spring for additional entries
    double mass;

    // RBE3 - master
    // int slaveID; // Dependent - use connectivity[1]
    //int slaveComponent;  - use dofDependent
    int numMaster; // Independent
    int *masterIDSet; // [numMaster]
    double *masterWeighting; // [numMaster]
    int *masterComponent; // [numMaster]

} feaConnectionStruct;

// Structure to hold coordinate system information
typedef struct {

    char *name; // Coordinate system name

    feaCoordSystemTypeEnum coordSystemType;  // Coordinate system type

    int coordSystemID; // ID number of coordinate system

    int refCoordSystemID; // ID of reference coordinate system

    double origin[3]; // x, y, and z coordinates for the origin

    double normal1[3]; // First normal direction
    double normal2[3]; // Second normal direction
    double normal3[3]; // Third normal direction - found from normal1 x normal2

} feaCoordSystemStruct;

// Structure to hold design variable relation information
typedef struct {

    char *name;

    feaDesignVariableTypeEnum relationType;

    int relationID;

    int numDesignVariable;
    char **designVariableNameSet; // Design variable borrowed references, size = [numDesignVariable]
 
    int fieldPosition; //  Position in card to apply design variable to
    char *fieldName; // Name of property/material to apply design variable to

    double constantRelationCoeff; // Constant terms of relation
    double *linearRelationCoeff; // Coefficients of linear relation, size = [numDesignVariable]

} feaDesignVariableRelationStruct;

// Structure to hold design variable information
typedef struct {

    char *name;

    feaDesignVariableTypeEnum designVariableType;

    int designVariableID; //  ID number of design variable

    double initialValue; // Initial value of design variable
    double lowerBound;   // Lower bounds of variable
    double upperBound;   // Upper bounds of variable
    double maxDelta;     // Move limit for design variable

    int numDiscreteValue; // Number of discrete values that a design variable can assume
    double *discreteValue; // List of discrete values that a design variable can assume

    int numMaterialID; // Number of materials to apply the design variable to
    int *materialSetID; // List of materials IDs
    int *materialSetType; // List of materials types corresponding to the materialSetID

    int numPropertyID;   // Number of property ID to apply the design variable to
    int *propertySetID; // List of property IDs
    int *propertySetType; // List of property types corresponding to the propertySetID

    int numElementID; // Number of element ID to apply the design variable to
    int *elementSetID; // List of element IDs
    int *elementSetType; // List of element types corresponding to the elementSetID
    int *elementSetSubType; // List of element subtypes correspoding to the elementSetID

    int fieldPosition; //  Position in card to apply design variable to
    char *fieldName; // Name of property/material to apply design variable to

    int numIndependVariable; // Number of independent variables this variables depends on
    char **independVariable; // List of independent variable names, size[numIndependVariable]
    int *independVariableID; // List of independent variable designVariableIDs, size[numIndependVariable]
    double *independVariableWeight; // List of independent variable weights, size[numIndependVariable]

    double variableWeight[2]; // Weight to apply to if variable is dependent

} feaDesignVariableStruct;


// Structure to hold design constraint information
typedef struct {
    char *name;

    int designConstraintID;  // ID number of design constraint

    char *responseType;  // Response type options for DRESP1 Entry

    double lowerBound;   // Lower bounds of design response
    double upperBound;   // Upper bounds of design response

    int numPropertyID;   // Number of property ID to apply the design variable to
    int *propertySetID; // List of property IDs
    int *propertySetType; // List of property types corresponding to the propertySetID

    int fieldPosition; //  Position in card to apply design variable to
    char *fieldName; // Name of property/material to apply design variable to


} feaDesignConstraintStruct;


typedef struct {

    int equationID;

    char *name; // Unique name

    int equationArraySize;
    char **equationArray;

} feaDesignEquationStruct;

typedef struct {

    int responseID;

    char *name; // Unique name

    char *responseType;
    char *propertyType;

    int region;

    int component; // Component number
    int itemCode; // Item code
    int modeNumber; // Mode number

    int lamina; // Lamina number
    double frequency; // Frequency value
    double time; // Time value
    int restraintFlag; // Restraint flag

    int gridID; // Grid ID
    int propertyID; // Property entry ID

} feaDesignResponseStruct;

typedef struct {

    int equationResponseID; 

    char *name; // Unique name

    char *equationName; // Name of the equation

    int region; // Region identifier for constant screening

    // DESVAR section
    int numDesignVariable;
    char **designVariableNameSet; // Design variable names, size = [numDesignVariable]

    // DTABLE section
    int numConstant;
    char **constantLabelSet; // Labels of the table constants, size = [numConstant]

    // DRESP1 section
    int numResponse;
    char **responseNameSet; // Names of design sensitivity response quantities, size = [numResponse]

    // DNODE section
    int numGrid;
    int *gridIDSet; // Grid IDs, size = [numGrid]
    int *dofNumberSet; // Degree of freedom numbers, size = [numGrid]

    // DRESP2 section
    int numEquationResponse;
    char **equationResponseNameSet; // Names of design sensitivity equation response quantities, size = [numEquationResponse]

} feaDesignEquationResponseStruct;

typedef struct {

    int numConstant;
    char **constantLabel; // size = [numConstant]
    double *constantValue; // size = [numConstant]

} feaDesignTableStruct;

typedef struct {

    int numParam;
    char **paramLabel; // size = [numParam]
    void **paramValue; // size = [numParam]
    int *paramType; // size = [numParam]

} feaDesignOptParamStruct;

// Structure to hold formatting information relevant to FEA file output
typedef struct {

    feaFileTypeEnum fileType;

    feaFileTypeEnum gridFileType; // GRID* 8 space + 4 x 16 space field entry instead of 10 x 8 space

} feaFileFormatStruct;

// Structure to hold FEA property information
typedef struct {
    char   *name;

    propertyTypeEnum propertyType;

    int propertyID; // ID number of property

    int materialID; // ID number of material
    char *materialName; // Name of material associated with material ID

    // Rods
    double crossSecArea; // Bar cross-sectional area
    double torsionalConst; // Torsional constant
    double torsionalStressReCoeff; // Torsional stress recovery coefficient
    double massPerLength; // Mass per unit length or Non-structural mass per unit length

    // Bar - see rod for additional variables
    double zAxisInertia; // Section moment of inertia about the z axis
    double yAxisInertia; // Section moment of inertia about the y axis
    double yCoords[4]; // Element y, z coordinates, in the bar cross-section, of
    double zCoords[4]; //    of four points at which to recover stresses
    double areaShearFactors[2]; // Area factors for shear
    double crossProductInertia; // Section cross-product of inertia

    char *crossSecType; //Type of cross section
    double crossSecDimension[10]; // Dimensions

    double orientationVec[3]; // Orientation vector

    // Shear

    // Shell
    double membraneThickness;  // Membrane thickness
    int materialBendingID;     // ID number of material for bending - if not specified and bendingInertiaRatio > 0 this value defaults to materialID
    double bendingInertiaRatio;// Ratio of actual bending moment inertia (I) to bending inertia of a solid
                               //   plate of thickness - default 1.0
    int materialShearID;       // ID number of material for shear - if not specified and shearMembraneRatio > 0 this value defaults to materialID
    double shearMembraneRatio; // Ratio of shear to membrane thickness  - default 5/6
    double massPerArea;        // Mass per unit area or Non-structural mass per unit area
    //double neutralPlaneDist[2];// Distances from the neutral plane of the plate to locations where
                               //   stress is calculate

    double zOffsetRel; // Offset from the surface of grid points to the element reference plane

    // Composite - more to be added
    //double distanceFromRefPlane; // Distance from reference plane to bottom surface of the element
    double compositeShearBondAllowable;  // Shear stress limit for bonding between plies
    char *compositeFailureTheory;    // HILL, HOFF, TSAI, STRN
    int compositeSymmetricLaminate;      // 1- SYM only half the plies are specified, for odd number plies 1/2 thickness of center ply is specified
                                // the first ply is the bottom ply in the stack, default (0) all plies specified
    int numPly;               // number of plies, size o
    int *compositeMaterialID;      // Vector of material ID's from bottom to top for all plies .size=[numPly]
    double *compositeThickness;        // Vector of thicknesses from bottom to top for all plies .size=[numPly]
    double *compositeOrientation;      // Vector of orientations from bottom to top for all plies .size=[numPly]

    // Solid

    // Concentrated Mass
    double mass; // Mass value
    double massOffset[3]; // Offset distance from the grid point to the center of gravity
    double massInertia[6]; // Mass moment of inertia measured at the mass center of gravity

} feaPropertyStruct;

// Structure to hold FEA material information
typedef struct {
    char   *name;

    materialTypeEnum materialType;

    int materialID; // ID number of material

    double youngModulus;    // E - Young's Modulus [Longitudinal if distinction is made]
    double shearModulus;    // G - Shear Modulus
    double poissonRatio;    // Poisson's Ratio
    double density;         // Rho - material mass density
    double thermalExpCoeff; // Coefficient of thermal expansion
    double temperatureRef;  // Reference temperature
    double dampingCoeff;    // Damping coefficient
    double yieldAllow;      // yield allowable for the isotropic material, populates tension
    double tensionAllow;    // Tension allowable for the material
    double compressAllow;   // Compression allowable for the material
    double shearAllow;      // Shear allowable for the material

    double youngModulusLateral; // Young's Modulus in the lateral direction
    double shearModulusTrans1Z; // Transverse shear modulus in the 1-Z plane
    double shearModulusTrans2Z; // Transverse shear modulus in the 2-Z plane

    double tensionAllowLateral;    // Lateral Tension allowable for the material
    double compressAllowLateral;   // Lateral Compression allowable for the material
    double thermalExpCoeffLateral; // Lateral Coefficient of thermal expansion
    int allowType; // 0 for stress, 1 for strain

} feaMaterialStruct;

// Structure to hold FEA unit system
typedef struct {

    char *length;      // Length unit
    char *mass;        // Mass unit

    char *pressure;    // Pressure unit
    char *densityVol;  // Density unit (mass/volume)
    char *densityArea; // Density unit (mass/area)
    char *temperature; // Temperature unit

    char *momentOfInertia; // moment of inertia (mass*length^2)

} feaUnitsStruct;


// Structure to hold FEA constraints
typedef struct {

    char *name;

    constraintTypeEnum constraintType;

    int constraintID; // ID number of constraint

    int numGridID; // Number of grid IDs in grid ID set
    int *gridIDSet; // List of grid IDs to apply constraint to. size = [numGridID]

    int dofConstraint; // Number to indicate DOF constraints
                       // For example 123 means x,y, and displacement are constrained,
                       //     123456 mean x,y,z, and all rotations are constrained,
                       //     23 means just y and z displacement are constrained,
                       //      etc.....

    double gridDisplacement; // The value for the displacement at gridID with constraint dofConstraint

} feaConstraintStruct;

// Structure to hold FEA supports
typedef struct {

    char *name;

    int supportID; // ID number of support

    int numGridID; // Number of grid IDs in grid ID set
    int *gridIDSet; // List of grid IDs to apply support to. size = [numGridID]

    int dofSupport; // Number to indicate DOF supports
                    // For example 123 means x,y, and displacement are supported,
                    //     123456 mean x,y,z, and all rotations are supported,
                    //     23 means just y and z displacement are supported,
                    //      etc.....

} feaSupportStruct;


// Structure to hold FEA loads
typedef struct {

    char *name;

    loadTypeEnum loadType;

    int loadID; // ID number of load

    double loadScaleFactor; // Scale factor for when combining loads

    // Concentrated force at a grid point
    int numGridID; // Number of grid IDs in grid ID set
    int *gridIDSet; // List of grid IDs to apply the constraint to. size = [numGridID]
    int coordSystemID; // Component number of coordinate system in which force vector is specified
    double forceScaleFactor; // Overall scale factor for the force
    double directionVector[3];   // [0]-x, [1]-y, [2]-z components of the force vector

    // Concentrated moment at a grid point (also uses coordSystemID and directionVector)
    double momentScaleFactor; // Overall scale factor for the moment

    // Gravitational load (also uses coordSystemID and directionVector)
    double gravityAcceleration; // Gravitational acceleration

    // Pressure load
    double pressureForce; // Pressure value
    double pressureDistributeForce[4]; // Pressure load at a specified grid location in the element

    double *pressureMultiDistributeForce; // Unique pressure load at a specified grid location for
                                          // each element in elementIDSet size = [4*numElementID] - used in type PressureExternal
                                          // where the pressure force is being provided by an external source (i.e. data transfer)

    int numElementID; // Number of elements IDs in element ID set
    int *elementIDSet; // List element IDs in which to apply the load. size = [numElementID]


    // Rotational velocity (also uses coordSystemID and directionVector)
    double angularVelScaleFactor; // Overall scale factor for the angular velocity
    double angularAccScaleFactor; // Overall scale factor for the angular acceleration

    // Thermal load - currently the temperature at a grid point - use gridIDSet
    double temperature; // Temperature value at given grid point
    double temperatureDefault; // Default temperature of grid point explicitly not used

} feaLoadStruct;

// Structure to hold FEA analysis
typedef struct {

    char *name; // Analysis name

    analysisTypeEnum analysisType; // Type of analysis

    int analysisID; // ID number of analysis

    // Loads for the analysis
    int numLoad;     // Number of loads in the analysis
    int *loadSetID; // List of the load IDs

    // Constraints for the analysis
    int numConstraint;   // Number of constraints in the analysis
    int *constraintSetID; // List of constraint IDs

    // Supports for the analysis
    int numSupport;   // Number of supports in the analysis
    int *supportSetID; // List of support IDs

    // Optimization - constraints
    int numDesignConstraint; // Number of design constraints
    int *designConstraintSetID; // List of design constraint IDs

    // Optimization - design response spanning
    int numDesignResponse; // Number of design responses
    int *designResponseSetID; // List of design response IDs

    // Eigenvalue
    char *extractionMethod;
    double frequencyRange[2];
    int numEstEigenvalue;
    int numDesiredEigenvalue;
    char *eigenNormaliztion;
    int gridNormaliztion;
    int componentNormaliztion;

    int lanczosMode; //Lanczos mode for calculating eigenvalues
    char *lanczosType; //Lanczos matrix type (DPB, DGB)

    // Trim
    int numMachNumber; // number of Mach numbers - flutter can have more than 1 per case
    double *machNumber; // Mach number
    double dynamicPressure; // Dynamic pressure
    double density; // Density of Air
    char *aeroSymmetryXY;
    char *aeroSymmetryXZ;

    int numRigidVariable; // Number of rigid trim variables
    char **rigidVariable; // List of character labels identifying rigid trim variables, size=[numRigidVariables]

    int numRigidConstraint; // Number of rigid trim constrained variables
    char **rigidConstraint; // List of character labels identifying rigid constrained trim variables, size=[numRigidConstraint]
    double *magRigidConstraint; // Magnitude of rigid constrained trim variables, size=[numRigidConstraint]

    int numControlConstraint; // Number of control surface constrained variables
    char **controlConstraint; // List of character labels identifying control surfaces to be constrained trim variables, size=[numControlConstraint]
    double *magControlConstraint; // Magnitude of control surface constrained variables, size=[numControlConstraint]

    // Flutter
    int numReducedFreq; // number of reduced frequencies numbers - flutter can have up to 8
    double *reducedFreq; // reduced frequency

} feaAnalysisStruct;

// Structure to hold FEA problem information
typedef struct {

    // Note: Setting order is important.
    // 1. Materials should be set before properties.
    // 2. Mesh should be set before loads and constraints
    // 3. Constraints and loads should be set before analysis
    // 4. Optimization should be set after properties, but before analysis

    // Analysis information
    int numAnalysis;
    feaAnalysisStruct *feaAnalysis; // size - [numAnalysis]

    // Property information
    int numProperty;
    feaPropertyStruct *feaProperty; // size = [numProperty]

    // Material information
    int numMaterial;
    feaMaterialStruct *feaMaterial; // size = [numMaterial]

    // Constraint information
    int numConstraint;
    feaConstraintStruct *feaConstraint; // size = [numConstraint]

    // Support information
    int numSupport;
    feaSupportStruct *feaSupport; // size = [numSupport]

    // Load information
    int numLoad;
    feaLoadStruct *feaLoad;  // size = [numLoad]

    // Mesh
    meshStruct feaMesh;

    // Connection information
    int numConnect;
    feaConnectionStruct *feaConnect;

    // Output formating
    feaFileFormatStruct feaFileFormat;

    // Optimization - Design variables
    int numDesignVariable;
    feaDesignVariableStruct *feaDesignVariable; // size = [numDesignVariable]

    // Optimization - Design variables relations
    int numDesignVariableRelation;
    feaDesignVariableRelationStruct *feaDesignVariableRelation; // size = [numDesignVariableRelation]

    // Optimization - Design constraints
    int numDesignConstraint;
    feaDesignConstraintStruct *feaDesignConstraint; // size = [numDesignConstraint]

    // Optimization - Equations
    int numEquation;
    feaDesignEquationStruct *feaEquation;

    // Optimization - Table Constants
    feaDesignTableStruct feaDesignTable;

    // Optimization - Table Constants
    feaDesignOptParamStruct feaDesignOptParam;

    // Optimization - Design Sensitivity Response Quantities
    int numDesignResponse;
    feaDesignResponseStruct *feaDesignResponse;

    // Optimization - Design Sensitivity Equation Response Quantities
    int numEquationResponse;
    feaDesignEquationResponseStruct *feaEquationResponse;

    // Coordinate systems
    int numCoordSystem;
    feaCoordSystemStruct *feaCoordSystem; // size = [numCoordSystem]

    // Aerodynamic information
    int numAero;
    feaAeroStruct *feaAero; // size = [numAero]

    feaAeroRefStruct feaAeroRef;

} feaProblemStruct;

#endif
