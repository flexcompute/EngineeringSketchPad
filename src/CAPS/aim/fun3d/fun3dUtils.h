// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Mach,
  Re,
  Viscoux,
  Equation_Type,
  Alpha,
  Beta,
  Overwrite_NML,
  Mesh_Format,
  Mesh_ASCII_Flag,
  Num_Iter,
  CFL_Schedule,
  CFL_Schedule_Iter,
  Restart_Read,
  Boundary_Condition,
  Use_Python_NML,
  Pressure_Scale_Factor,
  Pressure_Scale_Offset,
  NonInertial_Rotation_Rate,
  NonInertial_Rotation_Center,
  Two_Dimensional,
  Modal_Aeroelastic,
  Modal_Ref_Velocity,
  Modal_Ref_Length,
  Modal_Ref_Dynamic_Pressure,
  Time_Accuracy,
  Time_Step,
  Num_Subiter,
  Temporal_Error,
  Reference_Area,
  Moment_Length,
  Moment_Center,
  FUN3D_Version,
  Design_Variable,
  Design_Objective,
  Mesh,
  NUMINPUT = Mesh                /* Total number of inputs */
};

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
                                aimMeshRef *meshRef);

// Write FUN3D rubber.data file
int fun3d_writeRubber(void *aimInfo,
                      cfdDesignStruct design,
                      double fun3dVersion,
                      aimMeshRef *meshRef);

// Make FUN3D directory structure/tree
int fun3d_makeDirectory(void *aimInfo);

#ifdef __cplusplus
}
#endif
