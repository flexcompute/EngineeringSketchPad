# This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

# Cython declarations
# cython: language_level=3

# Import Cython *.pxds 
cimport cCAPS
cimport cfun3dNamelist
cimport cAIMUtil

# Python version
try:
    from cpython.version cimport PY_MAJOR_VERSION
except: 
    print("Error: Unable to import cpython.version")
#    raise ImportError

# Import Python modules
import sys
import os

try:
    import f90nml
except:
    print("\tUnable to import f90nml module - Try ' pip install f90nml '")
    raise

cdef extern from *:
    """
    #if defined(WIN32)
      #define PATH_MAX _MAX_PATH
    #else
      #include <limits.h>
    #endif
    """
    const int PATH_MAX

# Handles all byte to uni-code and and byte conversion. 
def _str(data,ignore_dicts = False, lower = True):
        
    if PY_MAJOR_VERSION >= 3:
        if isinstance(data, bytes):
            return data.decode("utf-8")
   
        if isinstance(data,list):
            return [_str(item,ignore_dicts=True) for item in data]
        if isinstance(data, dict) and not ignore_dicts:
            return {
                    _str(key,ignore_dicts=True): _str(value, ignore_dicts=True)
                    for key, value in data.iteritems()
                    }
    return data.lower() if lower else data


def read_fun3d_nml(fname):
#############################################################
### Read and evaluate fun3d namelist ########################
#############################################################    

    if os.path.isfile(fname):
        print('\n\tReading ' + fname + ' .....')
        # Read namelist
        nml = f90nml.read(fname)
        nml.row_major=False
        return nml
    else:
        print("\t" + fname + " not found!")

    return (None)

def write_fun3d_nml(nml,fname):
#############################################################
### Write an arbitary namelist to file ######################
#############################################################

    nml.write(fname, force=True)
        
cdef public int fun3d_writeNMLPython(void *aimInfo, 
                                     cCAPS.capsValue *aimInputs, 
                                     cfun3dNamelist.cfdBoundaryConditionStruct bcProps) except -1:
    
    cdef:
        char aimFile[PATH_MAX]

    index  = cAIMUtil.aim_getIndex(aimInfo, "Design_Functional", cCAPS.ANALYSISIN)-1
    index2 = cAIMUtil.aim_getIndex(aimInfo, "Design_SensFile", cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull or aimInputs[index2].vals.integer == cCAPS.true:
        Flow = b"Flow" + os.sep.encode()
    else:
        Flow = b""

    cAIMUtil.aim_file(aimInfo, Flow + b'fun3d.nml', aimFile)
    fname = _str(aimFile, lower=False)

    nml = None
    
    # Load fun3d namelis - if Overwrite_NML is true create a new namelist 
    index = cAIMUtil.aim_getIndex(aimInfo, "Overwrite_NML",cCAPS.ANALYSISIN)-1
    if aimInputs[index].vals.integer == cCAPS.false:
        nml = read_fun3d_nml(fname)
    else:
        print("\tOverwrite_NML is set to 'True' - a new namelist will be created")
        
    # Create a new namelist if no namelist is found  
    if nml is None:
        nml = f90nml.Namelist()
        print("\tCreating namelist")
    else:    
        print("\tAppending namelist")
    
    # Go through varible options/namelists
    
    # &project
    index = cAIMUtil.aim_getIndex(aimInfo, "Proj_Name", cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'project' not in nml:
            nml['project'] = f90nml.Namelist()
        
        nml['project']['project_rootname'] = _str(aimInputs[index].vals.string)

    # &raw_grid
    #index = cAIMUtil.aim_getIndex(aimInfo, "Mesh_Format",  cCAPS.ANALYSISIN)-1
    #if aimInputs[index].nullVal != cCAPS.IsNull:
    #    if 'raw_grid' not in nml:
    #        nml['raw_grid'] = f90nml.Namelist()
    #
    #    nml['raw_grid']['grid_format'] = _str(aimInputs[index].vals.string)
    #    
    #    index = cAIMUtil.aim_getIndex(aimInfo, "Mesh_ASCII_Flag",  cCAPS.ANALYSISIN)-1
    #    if aimInputs[index].vals.integer == <int> cCAPS.true:
    #       nml['raw_grid']['data_format'] = 'ascii'
    #    else:
    #       nml['raw_grid']['data_format'] = 'stream'
    if 'raw_grid' not in nml:
        nml['raw_grid'] = f90nml.Namelist()
    
    nml['raw_grid']['grid_format'] = "AFLR3"
    nml['raw_grid']['data_format'] = 'stream'

    index = cAIMUtil.aim_getIndex(aimInfo, "Two_Dimensional",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].vals.integer == <int> cCAPS.true:
        if 'raw_grid' not in nml:
            nml['raw_grid'] = f90nml.Namelist()

        nml['raw_grid']['twod_mode'] = True
    
    # &reference_physical_properties
    index = cAIMUtil.aim_getIndex(aimInfo, "Mach",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'reference_physical_properties' not in nml:
            nml['reference_physical_properties'] = f90nml.Namelist()
             
        nml['reference_physical_properties']['mach_number'] = aimInputs[index].vals.real
 
    index = cAIMUtil.aim_getIndex(aimInfo, "Re",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'reference_physical_properties' not in nml:
            nml['reference_physical_properties'] = f90nml.Namelist()
             
        nml['reference_physical_properties']['reynolds_number'] = aimInputs[index].vals.real

    index = cAIMUtil.aim_getIndex(aimInfo, "Alpha",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'reference_physical_properties' not in nml:
            nml['reference_physical_properties'] = f90nml.Namelist()
             
        nml['reference_physical_properties']['angle_of_attack'] = aimInputs[index].vals.real

    index = cAIMUtil.aim_getIndex(aimInfo, "Beta",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'reference_physical_properties' not in nml:
            nml['reference_physical_properties'] = f90nml.Namelist()
             
        nml['reference_physical_properties']['angle_of_yaw'] = aimInputs[index].vals.real
        
    index = cAIMUtil.aim_getIndex(aimInfo, "Temperature",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'reference_physical_properties' not in nml:
            nml['reference_physical_properties'] = f90nml.Namelist()
             
        nml['reference_physical_properties']['temperature'] = aimInputs[index].vals.real
        
        if aimInputs[index].units: 
            nml['reference_physical_properties']['temperature_units'] = _str(aimInputs[index].units)
    
    # &governing_equations
    index = cAIMUtil.aim_getIndex(aimInfo, "Viscous",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'governing_equations' not in nml:
            nml['governing_equations'] = f90nml.Namelist()
             
        nml['governing_equations']['viscous_terms'] = _str(aimInputs[index].vals.string)

    index = cAIMUtil.aim_getIndex(aimInfo, "Equation_Type",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'governing_equations' not in nml:
            nml['governing_equations'] = f90nml.Namelist()
             
        nml['governing_equations']['eqn_type'] = _str(aimInputs[index].vals.string)

    # &nonlinear_solver_parameters
    index = cAIMUtil.aim_getIndex(aimInfo, "Time_Accuracy",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
             
        nml['nonlinear_solver_parameters']['time_accuracy'] = _str(aimInputs[index].vals.string)

    index = cAIMUtil.aim_getIndex(aimInfo, "Time_Step",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
             
        nml['nonlinear_solver_parameters']['time_step_nondim'] = aimInputs[index].vals.real

    index = cAIMUtil.aim_getIndex(aimInfo, "Num_Subiter",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
             
        nml['nonlinear_solver_parameters']['subiterations'] = aimInputs[index].vals.integer

    index = cAIMUtil.aim_getIndex(aimInfo, "Temporal_Error",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
    
        nml['nonlinear_solver_parameters']['temporal_err_control'] = True
        nml['nonlinear_solver_parameters']['temporal_err_floor'] = aimInputs[index].vals.real

    index = cAIMUtil.aim_getIndex(aimInfo, "CFL_Schedule",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
    
        if 'schedule_cfl' not in nml['nonlinear_solver_parameters']:
            nml['nonlinear_solver_parameters']['schedule_cfl'] = [None]*2

        nml['nonlinear_solver_parameters']['schedule_cfl'][0] = aimInputs[index].vals.reals[0]
        nml['nonlinear_solver_parameters']['schedule_cfl'][1] = aimInputs[index].vals.reals[1]

    index = cAIMUtil.aim_getIndex(aimInfo, "CFL_Schedule_Iter",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'nonlinear_solver_parameters' not in nml:
            nml['nonlinear_solver_parameters'] = f90nml.Namelist()
    
        if 'schedule_iteration' not in nml['nonlinear_solver_parameters']:
            nml['nonlinear_solver_parameters']['schedule_iteration'] = [None]*2

        nml['nonlinear_solver_parameters']['schedule_iteration'][0] = aimInputs[index].vals.integers[0]
        nml['nonlinear_solver_parameters']['schedule_iteration'][1] = aimInputs[index].vals.integers[1]

    # &code_run_control
    index = cAIMUtil.aim_getIndex(aimInfo, "Num_Iter",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'code_run_control' not in nml:
            nml['code_run_control'] = f90nml.Namelist()
             
        nml['code_run_control']['steps'] = aimInputs[index].vals.integer

    index = cAIMUtil.aim_getIndex(aimInfo, "Restart_Read",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'code_run_control' not in nml:
            nml['code_run_control'] = f90nml.Namelist()
             
        nml['code_run_control']['restart_read'] = _str(aimInputs[index].vals.string)

    # &force_moment_integ_properties
    index = cAIMUtil.aim_getIndex(aimInfo, "Reference_Area",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'force_moment_integ_properties' not in nml:
            nml['force_moment_integ_properties'] = f90nml.Namelist()
  
        nml['force_moment_integ_properties']['area_reference']  = aimInputs[index].vals.real
    
    index = cAIMUtil.aim_getIndex(aimInfo, "Moment_Length",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'force_moment_integ_properties' not in nml:
            nml['force_moment_integ_properties'] = f90nml.Namelist()
        
        nml['force_moment_integ_properties']['x_moment_length'] = aimInputs[index].vals.reals[0]
        nml['force_moment_integ_properties']['y_moment_length'] = aimInputs[index].vals.reals[1]
  
    index = cAIMUtil.aim_getIndex(aimInfo, "Moment_Center",  cCAPS.ANALYSISIN)-1
    if aimInputs[index].nullVal != cCAPS.IsNull:
        if 'force_moment_integ_properties' not in nml:
            nml['force_moment_integ_properties'] = f90nml.Namelist()
     
        nml['force_moment_integ_properties']['x_moment_center'] = aimInputs[index].vals.reals[0]
        nml['force_moment_integ_properties']['y_moment_center'] = aimInputs[index].vals.reals[1]
        nml['force_moment_integ_properties']['z_moment_center'] = aimInputs[index].vals.reals[2]

    # &boundary_conditions
    for i in range(bcProps.numSurfaceProp):
        if 'boundary_conditions' not in nml:
            nml['boundary_conditions'] = f90nml.Namelist()
        
        index = bcProps.surfaceProp[i].bcID-1
        
        length = bcProps.surfaceProp[i].bcID
        
        # Wall Temperature 
        if bcProps.surfaceProp[i].wallTemperatureFlag == cCAPS.true:
            
            if 'wall_temperature' not in nml['boundary_conditions']:
                nml['boundary_conditions']['wall_temperature'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['wall_temperature']):
                    for _ in range(len(nml['boundary_conditions']['wall_temperature']), length):
                        nml['boundary_conditions']['wall_temperature'].append(None)
                             
            if 'wall_temp_flag' not in nml['boundary_conditions']:
                nml['boundary_conditions']['wall_temp_flag'] = [None]*length
            else:
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['wall_temp_flag']):
                    for _ in range(len(nml['boundary_conditions']['wall_temp_flag']), length):
                        nml['boundary_conditions']['wall_temp_flag'].append(None)
                        
            nml['boundary_conditions']['wall_temperature'][index] = bcProps.surfaceProp[i].wallTemperature
            nml['boundary_conditions']['wall_temp_flag'][index] = True
            
        # Total pressure and temperature 
        if bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.SubsonicInflow:

            if 'total_pressure_ratio' not in nml['boundary_conditions']:
                nml['boundary_conditions']['total_pressure_ratio'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['total_pressure_ratio']):
                    for _ in range(len(nml['boundary_conditions']['total_pressure_ratio']), length):
                        nml['boundary_conditions']['total_pressure_ratio'].append(None)
                             
            nml['boundary_conditions']['total_pressure_ratio'][index] = bcProps.surfaceProp[i].totalPressure
        
            if 'total_temperature_ratio' not in nml['boundary_conditions']:
                nml['boundary_conditions']['total_temperature_ratio'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['total_temperature_ratio']):
                    for _ in range(len(nml['boundary_conditions']['total_temperature_ratio']), length):
                        nml['boundary_conditions']['total_temperature_ratio'].append(None)
                             
            nml['boundary_conditions']['total_temperature_ratio'][index] = bcProps.surfaceProp[i].totalTemperature
        

        # Static pressure 
        if (bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.BackPressure or
            bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.SubsonicOutflow):

            if 'static_pressure_ratio' not in nml['boundary_conditions']:
                nml['boundary_conditions']['static_pressure_ratio'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['static_pressure_ratio']):
                    for _ in range(len(nml['boundary_conditions']['static_pressure_ratio']), length):
                        nml['boundary_conditions']['static_pressure_ratio'].append(None)
                             
            nml['boundary_conditions']['static_pressure_ratio'][index] = bcProps.surfaceProp[i].staticPressure

        # Mach number 
        if (bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.MachOutflow or
            bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.MassflowOut):

            if 'mach_bc' not in nml['boundary_conditions']:
                nml['boundary_conditions']['mach_bc'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['mach_bc']):
                    for _ in range(len(nml['boundary_conditions']['mach_bc']), length):
                        nml['boundary_conditions']['mach_bc'].append(None)
                             
            nml['boundary_conditions']['mach_bc'][index] = bcProps.surfaceProp[i].machNumber
        
        # Massflow
        if (bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.MassflowIn or
            bcProps.surfaceProp[i].surfaceType == cfun3dNamelist.MassflowOut):

            if 'massflow' not in nml['boundary_conditions']:
                nml['boundary_conditions']['massflow'] = [None]*length
            else: 
                if bcProps.surfaceProp[i].bcID > len(nml['boundary_conditions']['massflow']):
                    for _ in range(len(nml['boundary_conditions']['massflow']), length):
                        nml['boundary_conditions']['massflow'].append(None)
                             
            nml['boundary_conditions']['massflow'][index] = bcProps.surfaceProp[i].massflow

    # &noninertial_reference_frame
    if (aimInputs[cAIMUtil.aim_getIndex(aimInfo, "NonInertial_Rotation_Rate", cCAPS.ANALYSISIN)-1].nullVal != cCAPS.IsNull or
        aimInputs[cAIMUtil.aim_getIndex(aimInfo, "NonInertial_Rotation_Center",cCAPS.ANALYSISIN)-1].nullVal != cCAPS.IsNull):
                                   
        if 'noninertial_reference_frame' not in nml:
            nml['noninertial_reference_frame'] = f90nml.Namelist()
       
        nml['noninertial_reference_frame']['noninertial'] = True
    
        index = cAIMUtil.aim_getIndex(aimInfo, "NonInertial_Rotation_Center", cCAPS.ANALYSISIN)-1
        
        nml['noninertial_reference_frame']['rotation_center_x'] = aimInputs[index].vals.reals[0]
        nml['noninertial_reference_frame']['rotation_center_y'] = aimInputs[index].vals.reals[1] 
        nml['noninertial_reference_frame']['rotation_center_z'] = aimInputs[index].vals.reals[2]
    
        index = cAIMUtil.aim_getIndex(aimInfo, "NonInertial_Rotation_Rate", cCAPS.ANALYSISIN)-1
        
        nml['noninertial_reference_frame']['rotation_rate_x'] = aimInputs[index].vals.reals[0]
        nml['noninertial_reference_frame']['rotation_rate_y'] = aimInputs[index].vals.reals[1]
        nml['noninertial_reference_frame']['rotation_rate_z'] = aimInputs[index].vals.reals[2]
          
    # Write new namelist
    write_fun3d_nml(nml,fname)