// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "capsTypes.h" // Bring in CAPS types
#include "cfdTypes.h"  // Bring in CFD structures
#include "miscTypes.h"  // Bring in miscellaneous types

#ifdef __cplusplus
extern "C" {
#endif

// Fill bcProps in a cfdBoundaryConditionStruct format with boundary condition information from incoming BC Tuple
int cfd_getBoundaryCondition(void *aimInfo,
                             int numTuple,
                             capsTuple bcTuple[],
                             mapAttrToIndexStruct *attrMap,
                             cfdBoundaryConditionStruct *bcProps);

// Fill modalAeroelastic in a cfdModalAeroelasticStruct format with modal aeroelastic data from Modal Tuple
int cfd_getModalAeroelastic(int numTuple,
                            capsTuple modalTuple[],
                            cfdModalAeroelasticStruct *modalAeroelastic);

// Get the design variables from a capsTuple
int cfd_getDesignVariable(void *aimInfo,
                          int numDesignVariableTuple,
                          capsTuple designVariableTuple[],
                          int *numDesignVariable,
                          cfdDesignVariableStruct *variable[]);

// Fill objective in a cfdDesignFunctionalStruct format with objective data from Objective Tuple
int cfd_getDesignFunctional(void *aimInfo,
                           int numObjectiveTuple,
                           capsTuple objectiveTuple[],
                           cfdBoundaryConditionStruct *bcProps,
                           int numDesignVariable,
               /*@null@*/  cfdDesignVariableStruct variables[],
                           int *numObjective,
                           cfdDesignFunctionalStruct *objective[]);


// Initiate (0 out all values and NULL all pointers) of surfaceProps in the cfdSurfaceStruct structure format
int initiate_cfdSurfaceStruct(cfdSurfaceStruct *surfaceProps);

// Destroy (0 out all values and NULL all pointers) of bcProps in the cfdBoundaryConditionStruct structure format
int destroy_cfdSurfaceStruct(cfdSurfaceStruct *surfaceProps);

// Initiate (0 out all values and NULL all pointers) of bcProps in the cfdBoundaryConditionStruct structure format
int initiate_cfdBoundaryConditionStruct(cfdBoundaryConditionStruct *bcProps);

// Destroy (0 out all values and NULL all pointers) of surfaceProps in the cfdSurfaceStruct structure format
int destroy_cfdBoundaryConditionStruct(cfdBoundaryConditionStruct *bcProps);

// Initiate (0 out all values and NULL all pointers) of eigenValue in the cfdEigenValueStruct structure format
int initiate_cfdEigenValueStruct(cfdEigenValueStruct *eigenValue);

// Destroy (0 out all values and NULL all pointers) of eigenValue in the cfdEigenValueStruct structure format
int destroy_cfdEigenValueStruct(cfdEigenValueStruct *eigenValue);

// Initiate (0 out all values and NULL all pointers) of modalAeroelastic in the cfdModalAeroelasticStruct structure format
int initiate_cfdModalAeroelasticStruct(cfdModalAeroelasticStruct *modalAeroelastic);

// Destroy (0 out all values and NULL all pointers) of modalAeroelastic in the cfdModalAeroelasticStruct structure format
int destroy_cfdModalAeroelasticStruct(cfdModalAeroelasticStruct *modalAeroelastic);

// Initiate (0 out all values and NULL all pointers) of designVariable in the cfdDesignVariableStruct structure format
int initiate_cfdDesignVariableStruct(cfdDesignVariableStruct *designVariable);

// Destroy (0 out all values and NULL all pointers) of designVariable in the cfdDesignVariableStruct structure format
int destroy_cfdDesignVariableStruct(cfdDesignVariableStruct *designVariable);

// Copy cfdDesignVariableStruct structure
int copy_cfdDesignVariableStruct(void *aimInfo, cfdDesignVariableStruct *designVariable, cfdDesignVariableStruct *copy);

// Allocate cfdDesignVariableStruct structure
int allocate_cfdDesignVariableStruct(void *aimInfo, const char *name, const capsValue *var, cfdDesignVariableStruct *designVariable);

// Initiate (0 out all values and NULL all pointers) of objective in the cfdDesignFunctionalCompStruct structure format
int initiate_cfdDesignFunctionalCompStruct(cfdDesignFunctionalCompStruct *comp);

// Destroy (0 out all values and NULL all pointers) of objective in the cfdDesignFunctionalCompStruct structure format
int destroy_cfdDesignFunctionalCompStruct(cfdDesignFunctionalCompStruct *comp);

// Initiate (0 out all values and NULL all pointers) of objective in the cfdDesignFunctionalStruct structure format
int initiate_cfdDesignFunctionalStruct(cfdDesignFunctionalStruct *objective);

// Destroy (0 out all values and NULL all pointers) of objective in the cfdDesignFunctionalStruct structure format
int destroy_cfdDesignFunctionalStruct(cfdDesignFunctionalStruct *objective);

// Initiate (0 out all values and NULL all pointers) of objective in the cfdDesignStruct structure format
int initiate_cfdDesignStruct(cfdDesignStruct *design);

// Destroy (0 out all values and NULL all pointers) of objective in the cfdDesignStruct structure format
int destroy_cfdDesignStruct(cfdDesignStruct *design);

// Initiate (0 out all values and NULL all pointers) of objective in the cfdUnitsStruct structure format
int initiate_cfdUnitsStruct(cfdUnitsStruct *unit);

// Destroy (0 out all values and NULL all pointers) of objective in the cfdUnitsStruct structure format
int destroy_cfdUnitsStruct(cfdUnitsStruct *unit);

// Compute derived units from bas units
int cfd_cfdDerivedUnits(void *aimInfo, cfdUnitsStruct *units);

// Compute coefficient units from reference quantities
int cfd_cfdCoefficientUnits(void *aimInfo,
                            double length  , const char *lengthUnit,
                            double area    , const char *areaUnit,
                            double density , const char *densityUnit,
                            double speed   , const char *speedUnit,
                            double pressure, const char *pressureUnit,
                            cfdUnitsStruct *units);

#ifdef __cplusplus
}
#endif
