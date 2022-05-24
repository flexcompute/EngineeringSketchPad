#include <string.h>
#include <stdio.h>

#include "egads.h"     // Bring in egads utilss
#include "capsTypes.h" // Bring in CAPS types
#include "aimUtil.h"   // Bring in AIM utils
#include "aimMesh.h"// Bring in AIM meshing utils

#include "miscUtils.h" // Bring in misc. utility functions
#include "meshUtils.h" // Bring in meshing utility functions
#include "cfdTypes.h"  // Bring in cfd specific types
#include "su2Utils.h"  // Bring in su2 utility header

// Write SU2 configuration file for version Cardinal (4.0)
int su2_writeCongfig_Cardinal(void *aimInfo,  capsValue *aimInputs,
                              const char *meshfilename,
                              cfdBoundaryConditionStruct bcProps)
{

    int status; // Function return status

    int i, slen; // Indexing

    int stringLength;

    // For SU2 boundary tagging
    int counter;

    char *filename = NULL;
    FILE *fp = NULL;
    char fileExt[] = ".cfg";

    printf("Write SU2 configuration file for version \"Cardinal\"\n");

    stringLength = 1
                   + strlen(aimInputs[Proj_Name-1].vals.string)
                   + strlen(fileExt);

    filename = (char *) EG_alloc((stringLength +1)*sizeof(char));
    if (filename == NULL) {
        status =  EGADS_MALLOC;
        goto cleanup;
    }

    strcpy(filename, aimInputs[Proj_Name-1].vals.string);
    strcat(filename, fileExt);

    fp = aim_fopen(aimInfo, filename,"w");
    if (fp == NULL) {
        status =  CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    fprintf(fp,"%%                                                                              %%\n");
    fprintf(fp,"%% SU2 configuration file                                                       %%\n");
    fprintf(fp,"%% Created by SU2AIM for Project: \"%s\"\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%% File Version 4.1.1/2 \"Cardinal\"                                              %%\n");
    fprintf(fp,"%%                                                                              %%\n");
    fprintf(fp,"%% Please report bugs/comments/suggestions to NBhagat1@UDayton.edu              %%\n");
    fprintf(fp,"%%                                                                              %%\n");
    fprintf(fp,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Physical governing equations (EULER, NAVIER_STOKES,\n");
    fprintf(fp,"%%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,\n");
    fprintf(fp,"%%                               POISSON_EQUATION)\n");
    string_toUpperCase(aimInputs[Physical_Problem-1].vals.string);
    fprintf(fp,"PHYSICAL_PROBLEM= %s\n", aimInputs[Physical_Problem-1].vals.string);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Specify turbulence model (NONE, SA, SA_NEG, SST)\n");
    string_toUpperCase(aimInputs[Turbulence_Model-1].vals.string);
    fprintf(fp,"KIND_TURB_MODEL = %s\n", aimInputs[Turbulence_Model-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT)\n");
    fprintf(fp,"MATH_PROBLEM= DIRECT\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Restart solution (NO, YES)\n");
    fprintf(fp,"RESTART_SOL= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Regime type (COMPRESSIBLE, INCOMPRESSIBLE)\n");
    string_toUpperCase(aimInputs[Equation_Type-1].vals.string);
    fprintf(fp,"REGIME_TYPE= %s\n", aimInputs[Equation_Type-1].vals.string);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% System of measurements (SI, US)\n");
    fprintf(fp,"%% International system of units (SI): ( meters, kilograms, Kelvins,\n");
    fprintf(fp,"%%                                       Newtons = kg m/s^2, Pascals = N/m^2, \n");
    fprintf(fp,"%%                                       Density = kg/m^3, Speed = m/s,\n");
    fprintf(fp,"%%                                       Equiv. Area = m^2 )\n");
    fprintf(fp,"%% United States customary units (US): ( inches, slug, Rankines, lbf = slug ft/s^2, \n");
    fprintf(fp,"%%                                       psf = lbf/ft^2, Density = slug/ft^3, \n");
    fprintf(fp,"%%                                       Speed = ft/s, Equiv. Area = ft^2 )\n");
    string_toUpperCase(aimInputs[Unit_System-1].vals.string);
    fprintf(fp,"SYSTEM_MEASUREMENTS= %s\n", aimInputs[Unit_System-1].vals.string);
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mach number (non-dimensional, based on the free-stream values)\n");

    fprintf(fp,"MACH_NUMBER= %f\n", aimInputs[Mach-1].vals.real);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Angle of attack (degrees, only for compressible flows)\n");
    fprintf(fp,"AoA= %f\n", aimInputs[Alpha-1].vals.real);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Activate fixed lift mode (specify a CL instead of AoA, NO/YES)\n");
    fprintf(fp,"FIXED_CL_MODE= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Target coefficient of lift for fixed lift mode (0.0 by default)\n");
    fprintf(fp,"TARGET_CL= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Iterations to re-evaluate the angle of attack (100 by default)\n");
    fprintf(fp,"%%ITER_FIXED_CL= 100\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for fixed CL mode (0.1 by default)\n");
    //fprintf(fp,"DAMP_FIXED_CL= 0.2\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Side-slip angle (degrees, only for compressible flows)\n");
    fprintf(fp,"SIDESLIP_ANGLE= %f\n", aimInputs[Beta-1].vals.real);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Init option to choose between Reynolds (default) or thermodynamics quantities\n");
    fprintf(fp,"%% for initializing the solution (REYNOLDS, TD_CONDITIONS)\n");
    if (aimInputs[Init_Option-1].nullVal == NotNull) {
      string_toUpperCase(aimInputs[Init_Option-1].vals.string);
      fprintf(fp,"INIT_OPTION= %s\n", aimInputs[Init_Option-1].vals.string);
    }
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream option to choose between density and temperature (default) for\n");
    fprintf(fp,"%% initializing the solution (TEMPERATURE_FS, DENSITY_FS)\n");
    fprintf(fp,"FREESTREAM_OPTION= TEMPERATURE_FS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)\n");
    if (aimInputs[Freestream_Pressure-1].nullVal == NotNull) {
        fprintf(fp,"FREESTREAM_PRESSURE= %f\n",
                aimInputs[Freestream_Pressure-1].vals.real);
    } else {
        fprintf(fp,"FREESTREAM_PRESSURE= 101325.0\n");
    }
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream temperature (288.15 K, 518.67 R by default)\n");
    if (aimInputs[Freestream_Temperature-1].nullVal == NotNull) {
        fprintf(fp,"FREESTREAM_TEMPERATURE= %f\n",
                aimInputs[Freestream_Temperature-1].vals.real);
    }

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reynolds number (non-dimensional, based on the free-stream values)\n");
    fprintf(fp,"REYNOLDS_NUMBER= %e\n", aimInputs[Re-1].vals.real);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reynolds length (1 m, 1 inch by default)\n");
    fprintf(fp,"REYNOLDS_LENGTH= 1.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- INCOMPRESSIBLE FREE-STREAM DEFINITION ------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream density (1.2886 Kg/m^3, 0.0025 slug/ft^3 by default)\n");
    if (aimInputs[Freestream_Density-1].nullVal == NotNull) {
        fprintf(fp,"FREESTREAM_DENSITY= %f\n",
                aimInputs[Freestream_Density-1].vals.real);
    }

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream velocity (1.0 m/s, 1.0 ft/s by default)\n");
    if (aimInputs[Freestream_Velocity-1].nullVal == NotNull) {
        fprintf(fp,"FREESTREAM_VELOCITY= (%f, 0.0, 0.0) \n",
                aimInputs[Freestream_Velocity-1].vals.real);
    } else {
        fprintf(fp,"FREESTREAM_VELOCITY= (1.0, 0.0, 0.0)\n");
    }
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Free-stream viscosity (1.853E-5 N s/m^2, 3.87E-7 lbf s/ft^2 by default)\n");
    if (aimInputs[Freestream_Viscosity-1].nullVal == NotNull) {
        fprintf(fp,"FREESTREAM_VISCOSITY= %e\n",
                aimInputs[Freestream_Viscosity-1].vals.real);
    }

    fprintf(fp,"\n");
    fprintf(fp,"%% ---------------------- REFERENCE VALUE DEFINITION ---------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reference origin for moment computation (m or in)\n");
    if (aimInputs[Moment_Center-1].nullVal == NotNull) {
        fprintf(fp,"REF_ORIGIN_MOMENT_X= %f\n",
                aimInputs[Moment_Center-1].vals.reals[0]);
        fprintf(fp,"REF_ORIGIN_MOMENT_Y= %f\n",
                aimInputs[Moment_Center-1].vals.reals[1]);
        fprintf(fp,"REF_ORIGIN_MOMENT_Z= %f\n",
                aimInputs[Moment_Center-1].vals.reals[2]);

    } else {

        fprintf(fp,"REF_ORIGIN_MOMENT_X= 0.00\n");
        fprintf(fp,"REF_ORIGIN_MOMENT_Y= 0.00\n");
        fprintf(fp,"REF_ORIGIN_MOMENT_Z= 0.00\n");
    }

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reference length for pitching, rolling, and yawing non-dimensional\n");
    fprintf(fp,"%% moment (m or in)\n");
    if (aimInputs[Moment_Length-1].nullVal == NotNull) {
        fprintf(fp,"REF_LENGTH_MOMENT= %f\n", aimInputs[Moment_Length-1].vals.real);
    } else {
        fprintf(fp,"REF_LENGTH_MOMENT= 1.00\n");
    }
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reference area for force coefficients (0 implies automatic\n");
    fprintf(fp,"%% calculation) (m^2 or in^2)\n");
    if (aimInputs[Reference_Area-1].nullVal == NotNull) {
        fprintf(fp,"REF_AREA= %f\n", aimInputs[Reference_Area-1].vals.real);
    } else {
        fprintf(fp,"REF_AREA= 1.00\n");
    }
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,\n");
    fprintf(fp,"%%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)\n");
    string_toUpperCase(aimInputs[Reference_Dimensionalization-1].vals.string);
    fprintf(fp,"REF_DIMENSIONALIZATION= %s\n",
            aimInputs[Reference_Dimensionalization-1].vals.string);
    fprintf(fp,"\n");

    fprintf(fp,"%% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Different gas model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS)\n");
    fprintf(fp,"FLUID_MODEL= STANDARD_AIR\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Ratio of specific heats (1.4 default and the value is hardcoded\n");
    fprintf(fp,"%%                          for the model STANDARD_AIR)\n");
    fprintf(fp,"GAMMA_VALUE= 1.4\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Specific gas constant (287.058 J/kg*K default and this value is hardcoded \n");
    fprintf(fp,"%%                        for the model STANDARD_AIR)\n");
    fprintf(fp,"GAS_CONSTANT= 287.058\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Critical Temperature (131.00 K by default)\n");
    fprintf(fp,"CRITICAL_TEMPERATURE= 131.00\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Critical Pressure (3588550.0 N/m^2 by default)\n");
    fprintf(fp,"CRITICAL_PRESSURE= 3588550.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Critical Density (263.0 Kg/m3 by default)\n");
    fprintf(fp,"CRITICAL_DENSITY= 263.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Acentri factor (0.035 (air))\n");
    fprintf(fp,"ACENTRIC_FACTOR= 0.035\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------------- VISCOSITY MODEL ---------------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY).\n");
    fprintf(fp,"VISCOSITY_MODEL= SUTHERLAND\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Molecular Viscosity that would be constant (1.716E-5 by default)\n");
    fprintf(fp,"MU_CONSTANT= 1.716E-5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Sutherland Viscosity Ref (1.716E-5 default value for AIR SI)\n");
    fprintf(fp,"MU_REF= 1.716E-5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Sutherland Temperature Ref (273.15 K default value for AIR SI)\n");
    fprintf(fp,"MU_T_REF= 273.15\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Sutherland constant (110.4 default value for AIR SI)\n");
    fprintf(fp,"SUTHERLAND_CONSTANT= 110.4\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL).\n");
    fprintf(fp,"CONDUCTIVITY_MODEL= CONSTANT_PRANDTL\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Molecular Thermal Conductivity that would be constant (0.0257 by default)\n");
    fprintf(fp,"KT_CONSTANT= 0.0257\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------------------- UNSTEADY SIMULATION -------------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Unsteady simulation (NO, TIME_STEPPING, DUAL_TIME_STEPPING-1ST_ORDER, \n");
    fprintf(fp,"%%                      DUAL_TIME_STEPPING-2ND_ORDER, TIME_SPECTRAL)\n");
    fprintf(fp,"UNSTEADY_SIMULATION= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time Step for dual time stepping simulations (s) -- Only used when UNST_CFL_NUMBER = 0.0\n");
    fprintf(fp,"UNST_TIMESTEP= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Total Physical Time for dual time stepping simulations (s)\n");
    fprintf(fp,"UNST_TIME= 50.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Unsteady Courant-Friedrichs-Lewy number of the finest grid\n");
    fprintf(fp,"UNST_CFL_NUMBER= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of internal iterations (dual time method)\n");
    fprintf(fp,"UNST_INT_ITER= 200\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Integer number of periodic time instances for Time Spectral\n");
    fprintf(fp,"TIME_INSTANCES= 1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Iteration number to begin unsteady restarts\n");
    fprintf(fp,"UNST_RESTART_ITER= 0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ----------------------- DYNAMIC MESH DEFINITION -----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Dynamic mesh simulation (NO, YES)\n");
    fprintf(fp,"GRID_MOVEMENT= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Type of dynamic mesh (NONE, RIGID_MOTION, DEFORMING, ROTATING_FRAME,\n");
    fprintf(fp,"%%                       MOVING_WALL, STEADY_TRANSLATION, FLUID_STRUCTURE,\n");
    fprintf(fp,"%%                       AEROELASTIC, ELASTICITY, EXTERNAL,\n");
    fprintf(fp,"%%                       AEROELASTIC_RIGID_MOTION, GUST)\n");
    fprintf(fp,"GRID_MOVEMENT_KIND= DEFORMING\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Motion mach number (non-dimensional). Used for initializing a viscous flow\n");
    fprintf(fp,"%% with the Reynolds number and for computing force coeffs. with dynamic meshes.\n");
    fprintf(fp,"MACH_MOTION= 0.8\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Moving wall boundary marker(s) (NONE = no marker, ignored for RIGID_MOTION)\n");

    fprintf(fp,"MARKER_MOVING= (" );
    counter = 0;
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
                bcProps.surfaceProp[i].surfaceType == Viscous) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");
    //fprintf(fp,"MARKER_MOVING= ( NONE )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Coordinates of the motion origin\n");
    fprintf(fp,"MOTION_ORIGIN_X= 0.0\n");
    fprintf(fp,"MOTION_ORIGIN_Y= 0.0\n");
    fprintf(fp,"MOTION_ORIGIN_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Angular velocity vector (rad/s) about the motion origin\n");
    fprintf(fp,"ROTATION_RATE_X = 0.0\n");
    fprintf(fp,"ROTATION_RATE_Y = 0.0\n");
    fprintf(fp,"ROTATION_RATE_Z = 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Pitching angular freq. (rad/s) about the motion origin\n");
    fprintf(fp,"PITCHING_OMEGA_X= 0.0\n");
    fprintf(fp,"PITCHING_OMEGA_Y= 0.0\n");
    fprintf(fp,"PITCHING_OMEGA_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Pitching amplitude (degrees) about the motion origin\n");
    fprintf(fp,"PITCHING_AMPL_X= 0.0\n");
    fprintf(fp,"PITCHING_AMPL_Y= 0.0\n");
    fprintf(fp,"PITCHING_AMPL_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Pitching phase offset (degrees) about the motion origin\n");
    fprintf(fp,"PITCHING_PHASE_X= 0.0\n");
    fprintf(fp,"PITCHING_PHASE_Y= 0.0\n");
    fprintf(fp,"PITCHING_PHASE_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Translational velocity (m/s) in the x, y, & z directions\n");
    fprintf(fp,"TRANSLATION_RATE_X = 0.0\n");
    fprintf(fp,"TRANSLATION_RATE_Y = 0.0\n");
    fprintf(fp,"TRANSLATION_RATE_Z = 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Plunging angular freq. (rad/s) in x, y, & z directions\n");
    fprintf(fp,"PLUNGING_OMEGA_X= 0.0\n");
    fprintf(fp,"PLUNGING_OMEGA_Y= 0.0\n");
    fprintf(fp,"PLUNGING_OMEGA_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Plunging amplitude (m) in x, y, & z directions\n");
    fprintf(fp,"PLUNGING_AMPL_X= 0.0\n");
    fprintf(fp,"PLUNGING_AMPL_Y= 0.0\n");
    fprintf(fp,"PLUNGING_AMPL_Z= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Move Motion Origin for marker moving (1 or 0)\n");
    fprintf(fp,"MOVE_MOTION_ORIGIN = 0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------- AEROELASTIC SIMULATION (Typical Section Model) ---------------%%\n");
    fprintf(fp,"%% Activated by GRID_MOVEMENT_KIND option\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% The flutter speed index (modifies the freestream condition in the solver)\n");
    fprintf(fp,"FLUTTER_SPEED_INDEX = 0.6\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Natural frequency of the spring in the plunging direction (rad/s)\n");
    fprintf(fp,"PLUNGE_NATURAL_FREQUENCY = 100\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Natural frequency of the spring in the pitching direction (rad/s)\n");
    fprintf(fp,"PITCH_NATURAL_FREQUENCY = 100\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% The airfoil mass ratio\n");
    fprintf(fp,"AIRFOIL_MASS_RATIO = 60\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Distance in semichords by which the center of gravity lies behind\n");
    fprintf(fp,"%% the elastic axis\n");
    fprintf(fp,"CG_LOCATION = 1.8\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% The radius of gyration squared (expressed in semichords)\n");
    fprintf(fp,"%% of the typical section about the elastic axis\n");
    fprintf(fp,"RADIUS_GYRATION_SQUARED = 3.48\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Solve the aeroelastic equations every given number of internal iterations\n");
    fprintf(fp,"AEROELASTIC_ITER = 3\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------------- GUST SIMULATION ---------------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Apply a wind gust (NO, YES)\n");
    fprintf(fp,"WIND_GUST = NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Type of gust (NONE, TOP_HAT, SINE, ONE_M_COSINE, VORTEX, EOG)\n");
    fprintf(fp,"GUST_TYPE = NONE\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Direction of the gust (X_DIR or Y_DIR)\n");
    fprintf(fp,"GUST_DIR = Y_DIR\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Gust wavelenght (meters)\n");
    fprintf(fp,"GUST_WAVELENGTH= 10.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of gust periods\n");
    fprintf(fp,"GUST_PERIODS= 1.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Gust amplitude (m/s)\n");
    fprintf(fp,"GUST_AMPL= 10.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time at which to begin the gust (sec)\n");
    fprintf(fp,"GUST_BEGIN_TIME= 0.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Location at which the gust begins (meters) */\n");
    fprintf(fp,"GUST_BEGIN_LOC= 0.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------------------ SUPERSONIC SIMULATION ------------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Evaluate equivalent area on the Near-Field (NO, YES)\n");
    fprintf(fp,"EQUIV_AREA= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Integration limits of the equivalent area ( xmin, xmax, Dist_NearField )\n");
    fprintf(fp,"EA_INT_LIMIT= ( 1.6, 2.9, 1.0 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Equivalent area scale factor ( EA should be ~ force based objective functions )\n");
    fprintf(fp,"EA_SCALE_FACTOR= 1.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Fix an azimuthal line due to misalignments of the near-field\n");
    fprintf(fp,"FIX_AZIMUTHAL_LINE= 90.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Drag weight in sonic boom Objective Function (from 0.0 to 1.0)\n");
    fprintf(fp,"DRAG_IN_SONICBOOM= 0.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------------- ENGINE SIMULATION --------------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for the engine inflow.\n");
    fprintf(fp,"DAMP_ENGINE_INFLOW= 0.95\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for the engine exhaust.\n");
    fprintf(fp,"DAMP_ENGINE_EXHAUST= 0.95\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for the engine bleed.\n");
    fprintf(fp,"%%DAMP_ENGINE_BLEED= 0.95\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine nu factor (SA model).\n");
    fprintf(fp,"ENGINE_NU_FACTOR= 30.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Initialization with a subsonic flow around the engine.\n");
    fprintf(fp,"SUBSONIC_ENGINE= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Coordinates of the box that defines the subsonic region (Xmin, Ymin, Zmin,\n");
    fprintf(fp,"%%                                                          Xmax, Ymax, Zmax)\n");
    fprintf(fp,"%%SUBSONIC_ENGINE_BOX= ( -0.5, -0.49, 0.0, 2.5, 0.49, 0.0 )\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------- INVERSE DESIGN SIMULATION -----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Evaluate an inverse design problem using Cp (NO, YES)\n");
    fprintf(fp,"INV_DESIGN_CP= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Evaluate an inverse design problem using heat flux (NO, YES)\n");
    fprintf(fp,"INV_DESIGN_HEATFLUX= NO\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Euler wall boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_EULER= (" );

    counter = 0; // Euler boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {

        if (bcProps.surfaceProp[i].surfaceType == Inviscid) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( marker name, constant heat flux (J/m^2), ... )\n");
    fprintf(fp,"MARKER_HEATFLUX= (");

    counter = 0; // Viscous boundary w/ heat flux
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Viscous &&
                bcProps.surfaceProp[i].wallTemperatureFlag == (int) true &&
                bcProps.surfaceProp[i].wallTemperature < 0) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID, bcProps.surfaceProp[i].wallHeatFlux);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Navier-Stokes (no-slip), isothermal wall marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( marker name, constant wall temperature (K), ... )\n");
    fprintf(fp,"MARKER_ISOTHERMAL= (");

    counter = 0; // Viscous boundary w/ isothermal wall
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Viscous &&
            bcProps.surfaceProp[i].wallTemperatureFlag == (int) true &&
            bcProps.surfaceProp[i].wallTemperature >= 0) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID, bcProps.surfaceProp[i].wallTemperature);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Far-field boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_FAR= (" );

    counter = 0; // Farfield boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Farfield) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Symmetry boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_SYM= (" );

    counter = 0; // Symmetry boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Symmetry) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Near-Field boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_NEARFIELD= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Zone interface boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_INTERFACE= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Actuator disk boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%% rotation_angle_x-axis, rotation_angle_y-axis, rotation_angle_z-axis,\n");
    fprintf(fp,"%% root radius, tip radius, pressure jump, temperature jump, rev/min,\n");
    fprintf(fp,"%% uniform(0)/linear(1) distribution, ... )\n");
    fprintf(fp,"MARKER_ACTDISK= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Inlet boundary type (TOTAL_CONDITIONS, MASS_FLOW)\n");
    fprintf(fp,"INLET_TYPE= TOTAL_CONDITIONS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Inlet boundary marker(s) with the following formats (NONE = no marker) \n");
    fprintf(fp,"%% Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x, \n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"%% Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x, \n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"%% Incompressible: (inlet marker, NULL, velocity magnitude, flow_direction_x,\n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"MARKER_INLET= ( ");

    counter = 0; // Subsonic Inflow
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == SubsonicInflow) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f, %f, %f, %f, %f", bcProps.surfaceProp[i].bcID,
                                                  bcProps.surfaceProp[i].totalTemperature,
                                                  bcProps.surfaceProp[i].totalPressure,
                                                  bcProps.surfaceProp[i].uVelocity,
                                                  bcProps.surfaceProp[i].vVelocity,
                                                  bcProps.surfaceProp[i].wVelocity);
            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Supersonic inlet boundary marker(s) (NONE = no marker) \n");
    fprintf(fp,"%% Format: (inlet marker, temperature, static pressure, velocity_x, \n");
    fprintf(fp,"%%           velocity_y, velocity_z, ... ), i.e. primitive variables specified.\n");
    fprintf(fp,"MARKER_SUPERSONIC_INLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Outlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( outlet marker, back pressure (static), ... )\n");
    fprintf(fp,"MARKER_OUTLET= ( ");

    counter = 0; // Outlet boundary
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == BackPressure ||
            bcProps.surfaceProp[i].surfaceType == SubsonicOutflow) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID,
                                  bcProps.surfaceProp[i].staticPressure);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Supersonic outlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_SUPERSONIC_OUTLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Periodic boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y, \n");
    fprintf(fp,"%% rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis, \n");
    fprintf(fp,"%% rotation_angle_z-axis, translation_x, translation_y, translation_z, ... )\n");
    fprintf(fp,"MARKER_PERIODIC= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine inflow boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (engine inflow marker, fan face Mach, ... )\n");
    fprintf(fp,"MARKER_ENGINE_INFLOW= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine bleed boundary marker(s) with the following formats (NONE = no marker)\n");
    fprintf(fp,"%% Format: (engine bleed marker, mass flow rate, total temp, ... )\n");
    fprintf(fp,"%%MARKER_ENGINE_BLEED= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine exhaust boundary marker(s) with the following formats (NONE = no marker) \n");
    fprintf(fp,"%% Format: (engine exhaust marker, total nozzle temp, total nozzle pressure, ... )\n");
    fprintf(fp,"MARKER_ENGINE_EXHAUST= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Displacement boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( displacement marker, displacement value normal to the surface, ... )\n");
    fprintf(fp,"MARKER_NORMAL_DISPL= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Load boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( load marker, force value normal to the surface, ... )\n");
    fprintf(fp,"MARKER_NORMAL_LOAD= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Pressure boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( pressure marker )\n");
    fprintf(fp,"MARKER_PRESSURE= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Neumann bounday marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_NEUMANN= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Dirichlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_DIRICHLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Riemann boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (marker, data kind flag, list of data)\n");
    fprintf(fp,"MARKER_RIEMANN= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Non Reflecting boundary conditions marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (marker, data kind flag, list of data)\n");
    fprintf(fp,"MARKER_NRBC= ( NONE )\n");
    fprintf(fp,"%% ------------------------ SURFACES IDENTIFICATION ----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker(s) of the surface in the surface flow solution file\n");
    fprintf(fp,"MARKER_PLOTTING= (" );

    counter = 0; // Surface marker
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
            bcProps.surfaceProp[i].surfaceType == Viscous) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }
    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    // write monitoring information
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker(s) of the surface where the non-dimensional coefficients are evaluated.\n");
    fprintf(fp,"MARKER_MONITORING= (" );
    status = su2_marker(aimInfo, "Surface_Monitor", aimInputs, fp, bcProps);
    if (status != CAPS_SUCCESS) goto cleanup;

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker(s) of the surface where obj. func. (design problem) will be evaluated\n");
    fprintf(fp,"MARKER_DESIGNING = ( NONE )\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)\n");
    fprintf(fp,"NUM_METHOD_GRAD= GREEN_GAUSS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% CFL number (stating value for the adaptive CFL number)\n");
    fprintf(fp,"CFL_NUMBER= %f\n", aimInputs[CFL_Number-1].vals.real);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Adaptive CFL number (NO, YES)\n");
    fprintf(fp,"CFL_ADAPT= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,\n");
    fprintf(fp,"%%                                        CFL max value )\n");
    fprintf(fp,"CFL_ADAPT_PARAM= ( 1.5, 0.5, 1.25, 50.0 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Maximum Delta Time in local time stepping simulations\n");
    fprintf(fp,"MAX_DELTA_TIME= 1E6\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Runge-Kutta alpha coefficients\n");
    fprintf(fp,"RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Objective function in optimization problem (DRAG, LIFT, SIDEFORCE, MOMENT_X,\n");
    fprintf(fp,"%%                                             MOMENT_Y, MOMENT_Z, EFFICIENCY,\n");
    fprintf(fp,"%%                                             EQUIVALENT_AREA, NEARFIELD_PRESSURE,\n");
    fprintf(fp,"%%                                             FORCE_X, FORCE_Y, FORCE_Z, THRUST,\n");
    fprintf(fp,"%%                                             TORQUE, FREE_SURFACE, TOTAL_HEATFLUX,\n");
    fprintf(fp,"%%                                             MAXIMUM_HEATFLUX, INVERSE_DESIGN_PRESSURE,\n");
    fprintf(fp,"%%                                             INVERSE_DESIGN_HEATFLUX, AVG_TOTAL_PRESSURE, \n");
    fprintf(fp,"%%                                             MASS_FLOW_RATE)\n");
    fprintf(fp,"OBJECTIVE_FUNCTION= DRAG\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ----------------------- SLOPE LIMITER DEFINITION ----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reference element length for computing the slope and sharp edges \n");
    fprintf(fp,"%%                              limiters (0.1 m, 5.0 in by default)\n");
    fprintf(fp,"REF_ELEM_LENGTH= 0.1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Coefficient for the limiter\n");
    fprintf(fp,"LIMITER_COEFF= 0.3\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Freeze the value of the limiter after a number of iterations\n");
    fprintf(fp,"LIMITER_ITER= 999999\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Coefficient for the sharp edges limiter\n");
    fprintf(fp,"SHARP_EDGES_COEFF= 3.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reference coefficient (sensitivity) for detecting sharp edges.\n");
    fprintf(fp,"REF_SHARP_EDGES= 3.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Remove sharp edges from the sensitivity evaluation (NO, YES)\n");
    fprintf(fp,"SENS_REMOVE_SHARP= NO\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI, \n");
    fprintf(fp,"%%                                                      SMOOTHER_ILU0, SMOOTHER_LUSGS, \n");
    fprintf(fp,"%%                                                      SMOOTHER_LINELET)\n");
    fprintf(fp,"LINEAR_SOLVER= FGMRES\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Preconditioner of the Krylov linear solver (ILU0, LU_SGS, LINELET, JACOBI)\n");
    fprintf(fp,"LINEAR_SOLVER_PREC= LU_SGS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Minimum error of the linear solver for implicit formulations\n");
    fprintf(fp,"LINEAR_SOLVER_ERROR= 1E-4\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Max number of iterations of the linear solver for the implicit formulation\n");
    fprintf(fp,"LINEAR_SOLVER_ITER= 5\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------------- MULTIGRID PARAMETERS -----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Multi-grid levels (0 = no multi-grid)\n");
    fprintf(fp,"MGLEVEL= %d\n", aimInputs[MultiGrid_Level-1].vals.integer);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Multi-grid cycle (V_CYCLE, W_CYCLE, FULLMG_CYCLE)\n");
    fprintf(fp,"MGCYCLE= V_CYCLE\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Multi-grid pre-smoothing level\n");
    fprintf(fp,"MG_PRE_SMOOTH= ( 1, 2, 3, 3 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Multi-grid post-smoothing level\n");
    fprintf(fp,"MG_POST_SMOOTH= ( 0, 0, 0, 0 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Jacobi implicit smoothing of the correction\n");
    fprintf(fp,"MG_CORRECTION_SMOOTH= ( 0, 0, 0, 0 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for the residual restriction\n");
    fprintf(fp,"MG_DAMP_RESTRICTION= 0.75\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Damping factor for the correction prolongation\n");
    fprintf(fp,"MG_DAMP_PROLONGATION= 0.75\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, HLLC,\n");
    fprintf(fp,"%%                              TURKEL_PREC, MSW)\n");
    string_toUpperCase(aimInputs[Convective_Flux-1].vals.string);
    fprintf(fp,"CONV_NUM_METHOD_FLOW= %s\n", aimInputs[Convective_Flux-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)\n");
    fprintf(fp,"SPATIAL_ORDER_FLOW= 2ND_ORDER_LIMITER\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Slope limiter (VENKATAKRISHNAN, BARTH_JESPERSEN)\n");
    fprintf(fp,"SLOPE_LIMITER_FLOW= VENKATAKRISHNAN\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Entropy fix coefficient (0.0 implies no entropy fixing, 1.0 implies scalar\n");
    fprintf(fp,"%%                          artificial dissipation)\n");
    fprintf(fp,"ENTROPY_FIX_COEFF= 0.001\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% 1st, 2nd and 4th order artificial dissipation coefficients\n");
    fprintf(fp,"AD_COEFF_FLOW= ( 0.15, 0.5, 0.02 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Viscous limiter (NO, YES)\n");
    fprintf(fp,"VISCOUS_LIMITER_FLOW= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)\n");
    fprintf(fp,"TIME_DISCRE_FLOW= EULER_IMPLICIT\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Relaxation coefficient\n");
    fprintf(fp,"RELAXATION_FACTOR_FLOW= 1.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Convective numerical method (SCALAR_UPWIND)\n");
    fprintf(fp,"CONV_NUM_METHOD_TURB= SCALAR_UPWIND\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)\n");
    fprintf(fp,"SPATIAL_ORDER_TURB= 1ST_ORDER\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Slope limiter (VENKATAKRISHNAN)\n");
    fprintf(fp,"SLOPE_LIMITER_TURB= VENKATAKRISHNAN\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Viscous limiter (NO, YES)\n");
    fprintf(fp,"VISCOUS_LIMITER_TURB= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time discretization (EULER_IMPLICIT)\n");
    fprintf(fp,"TIME_DISCRE_TURB= EULER_IMPLICIT\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reduction factor of the CFL coefficient in the turbulence problem\n");
    fprintf(fp,"CFL_REDUCTION_TURB= 1.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Relaxation coefficient\n");
    fprintf(fp,"RELAXATION_FACTOR_TURB= 1.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------- HEAT NUMERICAL METHOD DEFINITION ----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Value of the thermal diffusivity\n");
    fprintf(fp,"THERMAL_DIFFUSIVITY= 1.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ---------------- ADJOINT-FLOW NUMERICAL METHOD DEFINITION -------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Convective numerical method (JST, LAX-FRIEDRICH, ROE)\n");
    fprintf(fp,"CONV_NUM_METHOD_ADJFLOW= JST\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)\n");
    fprintf(fp,"SPATIAL_ORDER_ADJFLOW= 2ND_ORDER\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Slope limiter (VENKATAKRISHNAN, SHARP_EDGES, WALL_DISTANCE)\n");
    fprintf(fp,"SLOPE_LIMITER_ADJFLOW= VENKATAKRISHNAN\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% 1st, 2nd, and 4th order artificial dissipation coefficients\n");
    fprintf(fp,"AD_COEFF_ADJFLOW= ( 0.15, 0.5, 0.02 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT)\n");
    fprintf(fp,"TIME_DISCRE_ADJFLOW= EULER_IMPLICIT\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Relaxation coefficient\n");
    fprintf(fp,"RELAXATION_FACTOR_ADJFLOW= 1.0\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reduction factor of the CFL coefficient in the adjoint problem\n");
    fprintf(fp,"CFL_REDUCTION_ADJFLOW= 0.8\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Limit value for the adjoint variable\n");
    fprintf(fp,"LIMIT_ADJFLOW= 1E6\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Multigrid adjoint problem (NO, YES)\n");
    fprintf(fp,"MG_ADJFLOW= YES\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ---------------- ADJOINT-TURBULENT NUMERICAL METHOD DEFINITION --------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Convective numerical method (SCALAR_UPWIND)\n");
    fprintf(fp,"CONV_NUM_METHOD_ADJTURB= SCALAR_UPWIND\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"SPATIAL_ORDER_ADJTURB= 1ST_ORDER\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Slope limiter (VENKATAKRISHNAN)\n");
    fprintf(fp,"SLOPE_LIMITER_ADJTURB= VENKATAKRISHNAN\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Time discretization (EULER_IMPLICIT)\n");
    fprintf(fp,"TIME_DISCRE_ADJTURB= EULER_IMPLICIT\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Reduction factor of the CFL coefficient in the adjoint turbulent problem\n");
    fprintf(fp,"CFL_REDUCTION_ADJTURB= 0.01\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ----------------------- GEOMETRY EVALUATION PARAMETERS ----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Geometrical evaluation mode (FUNCTION, GRADIENT)\n");
    fprintf(fp,"%%GEO_MODE= FUNCTION\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker(s) of the surface where geometrical based func. will be evaluated\n");
    fprintf(fp,"%%GEO_MARKER= ( airfoil )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of airfoil sections\n");
    fprintf(fp,"%%GEO_NUMBER_SECTIONS= 5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Orientation of airfoil sections (X_AXIS, Y_AXIS, Z_AXIS)\n");
    fprintf(fp,"%%GEO_ORIENTATION_SECTIONS= Y_AXIS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Location (coordinate) of the airfoil sections (MinValue, MaxValue)\n");
    fprintf(fp,"%%GEO_LOCATION_SECTIONS= (1.5, 3.5)\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Plot loads and Cp distributions on each airfoil section\n");
    fprintf(fp,"%%GEO_PLOT_SECTIONS= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of section cuts to make when calculating internal volume\n");
    fprintf(fp,"%%GEO_VOLUME_SECTIONS= 101\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ------------------------- GRID ADAPTATION STRATEGY --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Kind of grid adaptation (NONE, PERIODIC, FULL, FULL_FLOW, GRAD_FLOW,\n");
    fprintf(fp,"%%                          FULL_ADJOINT, GRAD_ADJOINT, GRAD_FLOW_ADJ, ROBUST,\n");
    fprintf(fp,"%%                          FULL_LINEAR, COMPUTABLE, COMPUTABLE_ROBUST,\n");
    fprintf(fp,"%%                          REMAINING, WAKE, SMOOTHING, SUPERSONIC_SHOCK)\n");
    fprintf(fp,"KIND_ADAPT= FULL_FLOW\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Percentage of new elements (%% of the original number of elements)\n");
    fprintf(fp,"NEW_ELEMS= 5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Scale factor for the dual volume\n");
    fprintf(fp,"DUALVOL_POWER= 0.5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Adapt the boundary elements (NO, YES)\n");
    fprintf(fp,"ADAPT_BOUNDARY= YES\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% ----------------------- DESIGN VARIABLE PARAMETERS --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Kind of deformation (TRANSLATION, ROTATION, SCALE,\n");
    fprintf(fp,"%%                      FFD_SETTING,\n");
    fprintf(fp,"%%                      FFD_CONTROL_POINT, FFD_CAMBER, FFD_THICKNESS\n");
    fprintf(fp,"%%                      FFD_DIHEDRAL_ANGLE, FFD_TWIST_ANGLE, FFD_ROTATION,\n");
    fprintf(fp,"%%                      FFD_CONTROL_POINT_2D, FFD_CAMBER_2D, FFD_THICKNESS_2D,\n");
    fprintf(fp,"%%                      HICKS_HENNE, PARABOLIC, NACA_4DIGITS, AIRFOIL)\n");
    fprintf(fp,"DV_KIND= SURFACE_FILE\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker of the surface in which we are going apply the shape deformation\n");
    fprintf(fp,"DV_MARKER= (" );

    // default to all inviscid and viscous surfaces if Surface_Deform is not set
    if (aimInputs[Surface_Deform-1].nullVal == IsNull) {

        counter = 0;
        for (i = 0; i < bcProps.numSurfaceProp ; i++) {
            if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
                bcProps.surfaceProp[i].surfaceType == Viscous) {

                if (counter > 0) fprintf(fp, ",");
                fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

                counter += 1;
            }
        }

        if(counter == 0) fprintf(fp," NONE");
        fprintf(fp," )\n");
    } else {
        status = su2_marker(aimInfo, "Surface_Deform", aimInputs, fp, bcProps);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Parameters of the shape deformation\n");
    fprintf(fp,"%% - TRANSLATION ( x_Disp, y_Disp, z_Disp ), as a unit vector\n");
    fprintf(fp,"%% - ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%% - SCALE ( 1.0 )\n");
    fprintf(fp,"%% - FFD_SETTING ( 1.0 )\n");
    fprintf(fp,"%% - FFD_CONTROL_POINT ( FFD_BoxTag, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )\n");
    fprintf(fp,"%% - FFD_CAMBER ( FFD_BoxTag, i_Ind, j_Ind )\n");
    fprintf(fp,"%% - FFD_THICKNESS ( FFD_BoxTag, i_Ind, j_Ind )\n");
    fprintf(fp,"%% - FFD_DIHEDRAL_ANGLE ( FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%% - FFD_TWIST_ANGLE ( FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%% - FFD_ROTATION ( FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%% - FFD_CONTROL_POINT_2D ( FFD_BoxTag, i_Ind, j_Ind, x_Disp, y_Disp )\n");
    fprintf(fp,"%% - FFD_CAMBER_2D ( FFD_BoxTag, i_Ind )\n");
    fprintf(fp,"%% - FFD_THICKNESS_2D ( FFD_BoxTag, i_Ind )\n");
    fprintf(fp,"%% - HICKS_HENNE ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc )\n");
    fprintf(fp,"%% - PARABOLIC ( Center, Thickness )\n");
    fprintf(fp,"%% - NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )\n");
    fprintf(fp,"%% - AIRFOIL ( 1.0 )\n");
    fprintf(fp,"DV_PARAM= ( 1, 0.5 )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Value of the shape deformation\n");
    fprintf(fp,"DV_VALUE= 0.01\n");
    fprintf(fp,"\n");

    fprintf(fp,"MOTION_FILENAME=%s_motion.dat\n", aimInputs[Proj_Name-1].vals.string);

    fprintf(fp,"\n");
    fprintf(fp,"%% ------------------------ GRID DEFORMATION PARAMETERS ------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Linear solver or smoother for implicit formulations (FGMRES, RESTARTED_FGMRES, BCGSTAB)\n");
    fprintf(fp,"DEFORM_LINEAR_SOLVER= FGMRES\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of smoothing iterations for mesh deformation\n");
    fprintf(fp,"DEFORM_LINEAR_ITER= 500\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of nonlinear deformation iterations (surface deformation increments)\n");
    fprintf(fp,"DEFORM_NONLINEAR_ITER= 3\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Print the residuals during mesh deformation to the console (YES, NO)\n");
    fprintf(fp,"DEFORM_CONSOLE_OUTPUT= YES\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Factor to multiply smallest cell volume for deform tolerance (0.001 default)\n");
    fprintf(fp,"DEFORM_TOL_FACTOR = 0.001\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Type of element stiffness imposed for FEA mesh deformation (INVERSE_VOLUME, \n");
    fprintf(fp,"%%                                          WALL_DISTANCE, CONSTANT_STIFFNESS)\n");
    fprintf(fp,"DEFORM_STIFFNESS_TYPE= INVERSE_VOLUME\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Visualize the deformation (NO, YES)\n");
    fprintf(fp,"VISUALIZE_DEFORMATION= YES\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% -------------------- FREE-FORM DEFORMATION PARAMETERS -----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Tolerance of the Free-Form Deformation point inversion\n");
    fprintf(fp,"FFD_TOLERANCE= 1E-10\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Maximum number of iterations in the Free-Form Deformation point inversion\n");
    fprintf(fp,"FFD_ITERATIONS= 500\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% FFD box definition: 3D case (FFD_BoxTag, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,\n");
    fprintf(fp,"%%                              X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7, X8, Y8, Z8)\n");
    fprintf(fp,"%%                     2D case (FFD_BoxTag, X1, Y1, 0.0, X2, Y2, 0.0, X3, Y3, 0.0, X4, Y4, 0.0,\n");
    fprintf(fp,"%%                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)\n");
    fprintf(fp,"FFD_DEFINITION= (MAIN_BOX, 0.5, 0.25, -0.25, 1.5, 0.25, -0.25, 1.5, 0.75, -0.25, 0.5, 0.75, -0.25, 0.5, 0.25, 0.25, 1.5, 0.25, 0.25, 1.5, 0.75, 0.25, 0.5, 0.75, 0.25)\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% FFD box degree: 3D case (x_degree, y_degree, z_degree)\n");
    fprintf(fp,"%%                 2D case (x_degree, y_degree, 0)\n");
    fprintf(fp,"FFD_DEGREE= (10, 10, 1)\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Surface continuity at the intersection with the FFD (1ST_DERIVATIVE, 2ND_DERIVATIVE)\n");
    fprintf(fp,"FFD_CONTINUITY= 2ND_DERIVATIVE\n");
    fprintf(fp,"\n");
    fprintf(fp,"%% --------------------------- CONVERGENCE PARAMETERS --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of total iterations\n");
    fprintf(fp,"EXT_ITER= %d\n", aimInputs[Num_Iter-1].vals.integer);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Convergence criteria (CAUCHY, RESIDUAL)\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"CONV_CRITERIA= RESIDUAL\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Residual reduction (order of magnitude with respect to the initial value)\n");
    fprintf(fp,"RESIDUAL_REDUCTION= %d\n", aimInputs[Residual_Reduction-1].vals.integer);

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Min value of the residual (log10 of the residual)\n");
    fprintf(fp,"RESIDUAL_MINVAL= -8\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Start convergence criteria at iteration number\n");
    fprintf(fp,"STARTCONV_ITER= 10\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Number of elements to apply the criteria\n");
    fprintf(fp,"CAUCHY_ELEMS= 100\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Epsilon to control the series convergence\n");
    fprintf(fp,"CAUCHY_EPS= 1E-10\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Direct function to apply the convergence criteria (LIFT, DRAG, NEARFIELD_PRESS)\n");
    fprintf(fp,"CAUCHY_FUNC_FLOW= DRAG\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Adjoint function to apply the convergence criteria (SENS_GEOMETRY, SENS_MACH)\n");
    fprintf(fp,"CAUCHY_FUNC_ADJFLOW= SENS_GEOMETRY\n");
    fprintf(fp,"\n");

    fprintf(fp,"%% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh input file\n");
    fprintf(fp,"MESH_FILENAME= %s\n", meshfilename);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh input file format (SU2, CGNS)\n");
    fprintf(fp,"MESH_FORMAT= SU2\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh output file\n");
    fprintf(fp,"MESH_OUT_FILENAME= %s.su2\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Restart flow input file\n");
    fprintf(fp,"SOLUTION_FLOW_FILENAME= solution_flow.dat\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Restart adjoint input file\n");
    fprintf(fp,"SOLUTION_ADJ_FILENAME= solution_adj.dat\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW,\n");
    fprintf(fp,"%%                     FIELDVIEW, FIELDVIEW_BINARY)\n");
    string_toUpperCase(aimInputs[Output_Format-1].vals.string);
    fprintf(fp,"OUTPUT_FORMAT= %s\n", aimInputs[Output_Format-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file convergence history (w/o extension) \n");
    fprintf(fp,"CONV_FILENAME= history_%s\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file with the forces breakdown\n");
    fprintf(fp,"BREAKDOWN_FILENAME= forces_breakdown_%s.dat\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file restart flow\n");
    fprintf(fp,"RESTART_FLOW_FILENAME= restart_flow_%s.dat\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file restart adjoint\n");
    fprintf(fp,"RESTART_ADJ_FILENAME= restart_adj.dat\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file flow (w/o extension) variables\n");
    fprintf(fp,"VOLUME_FLOW_FILENAME= flow_%s\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file adjoint (w/o extension) variables\n");
    fprintf(fp,"VOLUME_ADJ_FILENAME= adjoint\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output Objective function\n");
    fprintf(fp,"VALUE_OBJFUNC_FILENAME= of_eval.dat\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output objective function gradient (using continuous adjoint)\n");
    fprintf(fp,"GRAD_OBJFUNC_FILENAME= of_grad.dat\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file surface flow coefficient (w/o extension)\n");
    fprintf(fp,"SURFACE_FLOW_FILENAME= surface_flow_%s\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file surface adjoint coefficient (w/o extension)\n");
    fprintf(fp,"SURFACE_ADJ_FILENAME= surface_adjoint\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Writing solution file frequency\n");
    fprintf(fp,"WRT_SOL_FREQ= 1000\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Writing solution file frequency for physical time steps (dual time)\n");
    fprintf(fp,"WRT_SOL_FREQ_DUALTIME= 1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Writing convergence history frequency\n");
    fprintf(fp,"WRT_CON_FREQ= 1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Writing convergence history frequency (dual time, only written to screen)\n");
    fprintf(fp,"WRT_CON_FREQ_DUALTIME= 10\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output residual values in the solution files\n");
    fprintf(fp,"WRT_RESIDUALS= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output limiters values in the solution files\n");
    fprintf(fp,"WRT_LIMITERS= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output the sharp edges detector\n");
    fprintf(fp,"WRT_SHARPEDGES= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Minimize the required output memory\n");
    fprintf(fp,"LOW_MEMORY_OUTPUT= NO\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Verbosity of console output: NONE removes minor MPI overhead (NONE, HIGH)\n");
    fprintf(fp,"CONSOLE_OUTPUT_VERBOSITY= HIGH\n");
    fprintf(fp,"\n");

    fprintf(fp,"%% --------------------- OPTIMAL SHAPE DESIGN DEFINITION -----------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Available flow based objective functions or constraint functions\n");
    fprintf(fp,"%%    DRAG, LIFT, SIDEFORCE, EFFICIENCY,\n");
    fprintf(fp,"%%    FORCE_X, FORCE_Y, FORCE_Z,\n");
    fprintf(fp,"%%    MOMENT_X, MOMENT_Y, MOMENT_Z,\n");
    fprintf(fp,"%%    THRUST, TORQUE, FIGURE_OF_MERIT,\n");
    fprintf(fp,"%%    EQUIVALENT_AREA, NEARFIELD_PRESSURE, \n");
    fprintf(fp,"%%    TOTAL_HEATFLUX, MAXIMUM_HEATFLUX,\n");
    fprintf(fp,"%%    INVERSE_DESIGN_PRESSURE, INVERSE_DESIGN_HEATFLUX,\n");
    fprintf(fp,"%%    FREE_SURFACE, AVG_TOTAL_PRESSURE, MASS_FLOW_RATE\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Available geometrical based objective functions or constraint functions\n");
    fprintf(fp,"%%    MAX_THICKNESS, 1/4_THICKNESS, 1/2_THICKNESS, 3/4_THICKNESS, AREA, AOA, CHORD, \n");
    fprintf(fp,"%%    MAX_THICKNESS_SEC1, MAX_THICKNESS_SEC2, MAX_THICKNESS_SEC3, MAX_THICKNESS_SEC4, MAX_THICKNESS_SEC5, \n");
    fprintf(fp,"%%    1/4_THICKNESS_SEC1, 1/4_THICKNESS_SEC2, 1/4_THICKNESS_SEC3, 1/4_THICKNESS_SEC4, 1/4_THICKNESS_SEC5, \n");
    fprintf(fp,"%%    1/2_THICKNESS_SEC1, 1/2_THICKNESS_SEC2, 1/2_THICKNESS_SEC3, 1/2_THICKNESS_SEC4, 1/2_THICKNESS_SEC5, \n");
    fprintf(fp,"%%    3/4_THICKNESS_SEC1, 3/4_THICKNESS_SEC2, 3/4_THICKNESS_SEC3, 3/4_THICKNESS_SEC4, 3/4_THICKNESS_SEC5, \n");
    fprintf(fp,"%%    AREA_SEC1, AREA_SEC2, AREA_SEC3, AREA_SEC4, AREA_SEC5, \n");
    fprintf(fp,"%%    AOA_SEC1, AOA_SEC2, AOA_SEC3, AOA_SEC4, AOA_SEC5, \n");
    fprintf(fp,"%%    CHORD_SEC1, CHORD_SEC2, CHORD_SEC3, CHORD_SEC4, CHORD_SEC5\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Available design variables\n");
    fprintf(fp,"%%    HICKS_HENNE   (  1, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc )\n");
    fprintf(fp,"%%    COSINE_BUMP (  2, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc, x_Size )\n");
    fprintf(fp,"%%    SPHERICAL   (  3, Scale | Mark. List | ControlPoint_Index, Theta_Disp, R_Disp )\n");
    fprintf(fp,"%%    NACA_4DIGITS  (  4, Scale | Mark. List | 1st digit, 2nd digit, 3rd and 4th digit )\n");
    fprintf(fp,"%%    DISPLACEMENT  (  5, Scale | Mark. List | x_Disp, y_Disp, z_Disp )\n");
    fprintf(fp,"%%    ROTATION    (  6, Scale | Mark. List | x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )\n");
    fprintf(fp,"%%    FFD_CONTROL_POINT (  7, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind, k_Ind, x_Mov, y_Mov, z_Mov )\n");
    fprintf(fp,"%%    FFD_DIHEDRAL_ANGLE  (  8, Scale | Mark. List | FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%%    FFD_TWIST_ANGLE   (  9, Scale | Mark. List | FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%%    FFD_ROTATION  ( 10, Scale | Mark. List | FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%%    FFD_CAMBER  ( 11, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind )\n");
    fprintf(fp,"%%    FFD_THICKNESS   ( 12, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind )\n");
    fprintf(fp,"%%    FOURIER     ( 14, Scale | Mark. List | Lower(0)/Upper(1) side, index, cos(0)/sin(1) )\n");
    fprintf(fp,"%%    FFD_CONTROL_POINT_2D (  15, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind, x_Mov, y_Mov )\n");
    fprintf(fp,"%%    FFD_CAMBER_2D   ( 16, Scale | Mark. List | FFD_BoxTag, i_Ind )\n");
    fprintf(fp,"%%    FFD_THICKNESS_2D  ( 17, Scale | Mark. List | FFD_BoxTag, i_Ind )\n");
    fprintf(fp,"%%    FFD_CONTROL_SURFACE   ( 18, Scale | Mark. List | FFD_BoxTag, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Optimization objective function with scaling factor\n");
    fprintf(fp,"%% ex= Objective * Scale\n");
    fprintf(fp,"OPT_OBJECTIVE= DRAG * 0.001\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Optimization constraint functions with scaling factors, separated by semicolons\n");
    fprintf(fp,"%% ex= (Objective = Value ) * Scale, use '>','<','='\n");
    fprintf(fp,"OPT_CONSTRAINT= ( LIFT > 0.328188 ) * 0.001; ( MOMENT_Z > 0.034068 ) * 0.001; ( MAX_THICKNESS > 0.11 ) * 0.001\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Maximum number of iterations\n");
    fprintf(fp,"OPT_ITERATIONS= 100\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Requested accuracy\n");
    fprintf(fp,"OPT_ACCURACY= 1E-6\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Upper bound for each design variable\n");
    fprintf(fp,"OPT_BOUND_UPPER= 0.1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Lower bound for each design variable\n");
    fprintf(fp,"OPT_BOUND_LOWER= -0.1\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Optimization design variables, separated by semicolons\n");
    fprintf(fp,"DEFINITION_DV= ( 1, 1.0 | airfoil | 0, 0.05 ); ( 1, 1.0 | airfoil | 0, 0.10 ); ( 1, 1.0 | airfoil | 0, 0.15 ); ( 1, 1.0 | airfoil | 0, 0.20 ); ( 1, 1.0 | airfoil | 0, 0.25 ); ( 1, 1.0 | airfoil | 0, 0.30 ); ( 1, 1.0 | airfoil | 0, 0.35 ); ( 1, 1.0 | airfoil | 0, 0.40 ); ( 1, 1.0 | airfoil | 0, 0.45 ); ( 1, 1.0 | airfoil | 0, 0.50 ); ( 1, 1.0 | airfoil | 0, 0.55 ); ( 1, 1.0 | airfoil | 0, 0.60 ); ( 1, 1.0 | airfoil | 0, 0.65 ); ( 1, 1.0 | airfoil | 0, 0.70 ); ( 1, 1.0 | airfoil | 0, 0.75 ); ( 1, 1.0 | airfoil | 0, 0.80 ); ( 1, 1.0 | airfoil | 0, 0.85 ); ( 1, 1.0 | airfoil | 0, 0.90 ); ( 1, 1.0 | airfoil | 0, 0.95 ); ( 1, 1.0 | airfoil | 1, 0.05 ); ( 1, 1.0 | airfoil | 1, 0.10 ); ( 1, 1.0 | airfoil | 1, 0.15 ); ( 1, 1.0 | airfoil | 1, 0.20 ); ( 1, 1.0 | airfoil | 1, 0.25 ); ( 1, 1.0 | airfoil | 1, 0.30 ); ( 1, 1.0 | airfoil | 1, 0.35 ); ( 1, 1.0 | airfoil | 1, 0.40 ); ( 1, 1.0 | airfoil | 1, 0.45 ); ( 1, 1.0 | airfoil | 1, 0.50 ); ( 1, 1.0 | airfoil | 1, 0.55 ); ( 1, 1.0 | airfoil | 1, 0.60 ); ( 1, 1.0 | airfoil | 1, 0.65 ); ( 1, 1.0 | airfoil | 1, 0.70 ); ( 1, 1.0 | airfoil | 1, 0.75 ); ( 1, 1.0 | airfoil | 1, 0.80 ); ( 1, 1.0 | airfoil | 1, 0.85 ); ( 1, 1.0 | airfoil | 1, 0.90 ); ( 1, 1.0 | airfoil | 1, 0.95 )\n");
    fprintf(fp,"%%\n");
    if (aimInputs[Input_String-1].nullVal != IsNull) {
        fprintf(fp,"%% CAPS Input_String\n");
        for (slen = i = 0; i < aimInputs[Input_String-1].length; i++) {
            string_toUpperCase(aimInputs[Input_String-1].vals.string + slen);
            fprintf(fp,"%s\n", aimInputs[Input_String-1].vals.string + slen);
            slen += strlen(aimInputs[Input_String-1].vals.string + slen) + 1;
        }
    }
    fprintf(fp,"\n");

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:

        if (filename != NULL) EG_free(filename);
        if (fp != NULL) fclose(fp);

        return status;
}
