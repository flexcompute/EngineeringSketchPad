// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifdef __cplusplus
extern "C" {
#endif

// Extract the FEPOINT Tecoplot data from a FUN3D Aero-Loads file (connectivity is ignored) - dataMatrix = [numVariable][numDataPoint]
int fun3d_readAeroLoad(void *aimInfo, char *filename, int *numVariable, char **variableName[],
                       int *numDataPoint, double ***dataMatrix);

// Create a 3D mesh for FUN3D from a 2D mesh
int fun3d_2DMesh(void *aimInfo,
                 aimMeshRef *meshRef,
                 const char *projectName,
                 mapAttrToIndexStruct *attrMap,
                 cfdBoundaryConditionStruct *bcProps);

// Write FUN3D data transfer files
int fun3d_dataTransfer(void *aimInfo,
                       const char *projectName,
                       mapAttrToIndexStruct *groupMap,
                       cfdBoundaryConditionStruct bcProps,
                       aimMeshRef *meshRef,
                       /*@null@*/ cfdModalAeroelasticStruct *eigenVector);

// Write FUN3D fun3d.nml file
int fun3d_writeNML(void *aimInfo, capsValue *aimInputs,
                   cfdBoundaryConditionStruct bcProps);

// Write FUN3D movingbody.input file
int fun3d_writeMovingBody(void *aimInfo, double fun3dVersion,
                          cfdBoundaryConditionStruct bcProps,
                          cfdModalAeroelasticStruct *modalAeroelastic);

// Write FUN3D parametrization/sensitivity file
int fun3d_writeParameterization(void *aimInfo,
                                int numDesignVariable,
                                cfdDesignVariableStruct designVariable[],
                                int writeSensitivity,
                                aimMeshRef *meshRef);

// Write FUN3D rubber.data file
int fun3d_writeRubber(void *aimInfo,
                      cfdDesignStruct design,
                      double fun3dVersion,
                      aimMeshRef *meshRef);

// Read FUN3D rubber.data file
int fun3d_readRubber(void *aimInfo,
                     cfdDesignStruct design,
                     double fun3dVersion);

// Make FUN3D directory structure/tree
int fun3d_makeDirectory(void *aimInfo);

#ifdef __cplusplus
}
#endif
