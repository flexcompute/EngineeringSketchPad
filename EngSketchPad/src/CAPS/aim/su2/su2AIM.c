/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             SU2 AIM
 *
 *     Modified from code written by Dr. Ryan Durscher AFRL/RQVC
 *
 */


/*!\mainpage Introduction
 *
 * \section overviewSU2 SU2 AIM Overview
 * This module can be used to interface with the open-source CFD code, SU2 \cite Palacios2013 \cite Palacios2014
 *  with geometry in the CAPS system. For SU2 capabilities and related documentation, please refer to
 * http://su2.stanford.edu/. SU2 expects a volume mesh file and a corresponding
 * configuration file to perform the analysis.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsSU2 and \ref aimOutputsSU2, respectively.
 *
 * Details on the use of units are outlined in \ref aimUnitsSU2.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferSU2
 *
 * \subsection meshSU2 Automatic generation of SU2 Mesh file
 * The volume mesh file from SU2 AIM is written in native SU2
 * format ("filename.su2"). The description of the native SU2 mesh can be
 * found SU2 website. For the automatic generation of mesh file, SU2 AIM
 * depends on Mesh AIMs, for example, TetGen or AFLR4/3 AIM.
 *
 * \subsection configSU2 Automatic generation of SU2 Configuration file
 * The configuration file ("filename.cfg") from SU2 AIM is automatically
 * created by using the flow features and boundary conditions that were set in
 * the driver program as a user input. For the rest of the configuration
 * variables, default set of values are provided for a general execution.
 * If desired, a user has freedom to manually (a) change these variables based
 * on  personal preference, or (b) override the configuration file with unique
 * configuration variables.
 *
 * \section su2examples SU2 Examples
 *  Here is an example that illustrated use of SU2 AIM \ref su2Example. Note
 *  this AIM uses TetGen AIM for volume mesh generation.
 *
 */

#ifdef HAVE_PYTHON
#include "Python.h" // Bring in Python API
#endif

#include <string.h>
#include <ctype.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"
#include "su2Writer.h"

#include "meshUtils.h"
#include "cfdUtils.h"
#include "miscUtils.h"
#include "su2Utils.h"


#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#else
#include <unistd.h>
#endif

#define MXCHAR  255

//#define DEBUG

#define NUMOUTPUT  3*9

typedef struct {

    // SU2 project name
    char *projectName;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    // Point to caps input value for version of Su2
    capsValue *su2Version;

    // Pointer to caps input value for scaling pressure during data transfer
    capsValue *pressureScaleFactor;

    // Pointer to caps input value for offset pressure during data transfer
    capsValue *pressureScaleOffset;

    // Units structure
    cfdUnitsStruct units;

} aimStorage;

//#include "su2DataExchange.c"


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@null@*/ /*@unused@*/ const char *unitSys, void *aimInfo,
                  void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS; // Function return
    int *ints=NULL, i;
    char **strs=NULL;
    const char *keyWord;
    char *keyValue = NULL;
    double real = 1;
    cfdUnitsStruct *units=NULL;

    aimStorage *su2Instance=NULL;

    #ifdef DEBUG
        printf("\n su2AIM/aimInitialize   inst = %d!\n", inst);
    #endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 5;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
    strs[0]  = EG_strdup("Pressure");
    strs[1]  = EG_strdup("P");
    strs[2]  = EG_strdup("Cp");
    strs[3]  = EG_strdup("CoefficientOfPressure");
    strs[4]  = EG_strdup("Displacement");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
    *fnames  = strs;

    /* specify the dimension of each field variable */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = 1;
    ints[1]  = 1;
    ints[2]  = 1;
    ints[3]  = 1;
    ints[4]  = 3;
    *franks   = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldOut;
    ints[1]  = FieldOut;
    ints[2]  = FieldOut;
    ints[3]  = FieldOut;
    ints[4]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate su2Instance
    AIM_ALLOC(su2Instance, 1, aimStorage, aimInfo, status);
    *instStore = su2Instance;

    // Set initial values for su2Instance
    su2Instance->projectName = NULL;

    // Version
    su2Instance->su2Version = NULL;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&su2Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    // Pointer to caps input value for scaling pressure during data transfer
    su2Instance->pressureScaleFactor = NULL;

    // Pointer to caps input value for off setting pressure during data transfer
    su2Instance->pressureScaleOffset = NULL;

    initiate_cfdUnitsStruct(&su2Instance->units);

    /*! \page aimUnitsSU2 AIM Units
     *  A unit system may be optionally specified during AIM instance initiation. If
     *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
     *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
     *  unit consistent. A unit system may be specified via a JSON string dictionary for example:
     *  unitSys = "{"mass": "kg", "length": "m", "time":"seconds", "temperature": "K"}"
     */
    if (unitSys != NULL) {
      units = &su2Instance->units;

      // Do we have a json string?
      if (strncmp( unitSys, "{", 1) != 0) {
        AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
        return CAPS_BADVALUE;
      }

      /*! \page aimUnitsSU2
       *  \section jsonStringSU2 JSON String Dictionary
       *  The key arguments of the dictionary are described in the following:
       *
       *  <ul>
       *  <li> <B>mass = "None"</B> </li> <br>
       *  Mass units - e.g. "kilogram", "k", "slug", ...
       *  </ul>
       */
      keyWord = "mass";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->mass = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->mass, &real, "kg", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->mass, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsSU2
       *  <ul>
       *  <li> <B>length = "None"</B> </li> <br>
       *  Length units - e.g. "meter", "m", "inch", "in", "mile", ...
       *  </ul>
       */
      keyWord = "length";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->length = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->length, &real, "m", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->length, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsSU2
       *  <ul>
       *  <li> <B>time = "None"</B> </li> <br>
       *  Time units - e.g. "second", "s", "minute", ...
       *  </ul>
       */
      keyWord = "time";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->time = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->time, &real, "s", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->time, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsSU2
       *  <ul>
       *  <li> <B>temperature = "None"</B> </li> <br>
       *  Temperature units - e.g. "Kelvin", "K", "degC", ...
       *  </ul>
       */
      keyWord = "temperature";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->temperature = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->temperature, &real, "K", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->temperature, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      status = cfd_cfdDerivedUnits(aimInfo, units);
      AIM_STATUS(aimInfo, status);
    }

cleanup:
    if (status != CAPS_SUCCESS) {
        /* release all possibly allocated memory on error */
        if (*fnames != NULL)
          for (i = 0; i < *nFields; i++) AIM_FREE((*fnames)[i]);
        AIM_FREE(*franks);
        AIM_FREE(*fInOut);
        AIM_FREE(*fnames);
        AIM_FREE(*instStore);
        *nFields = 0;
    }

    return status;
}


int aimInputs(void *instStore, /*@unused@*/ void *aimInfo, int index,
              char **ainame, capsValue *defval)
{
    int status = CAPS_SUCCESS;
    aimStorage *su2Instance;
    cfdUnitsStruct *units=NULL;

#ifdef DEBUG
    printf(" su2AIM/aimInputs  index = %d!\n", index);
#endif

    *ainame = NULL;

    // SU2 Inputs
    /*! \page aimInputsSU2 AIM Inputs
     * For the description of the configuration variables, associated values,
     * and available options refer to the template configuration file that is
     * distributed with SU2.
     * Note: The configuration file is dependent on the version of SU2 used.
     * This configuration file that will be auto generated is compatible with
     * SU2 4.1.1. (Cardinal), 5.0.0 (Raven), 6.2.0 (Falcon) or 7.2.0 (Blackbird - Default)
     */

    su2Instance = (aimStorage *) instStore;

    if (su2Instance != NULL) units = &su2Instance->units;

    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("su2_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsSU2
         * - <B> Proj_Name = "su2_CAPS"</B> <br>
         * This corresponds to the project "root" name.
         */

    } else if (index == Mach) {
        *ainame              = EG_strdup("Mach"); // Mach number
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsSU2
         * - <B> Mach = NULL</B> <br>
         * Mach number; this corresponds to the MACH_NUMBER keyword in the configuration file.
         */

    } else if (index == Re) {
        *ainame              = EG_strdup("Re"); // Reynolds number
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsSU2
         * - <B> Re = NULL</B> <br>
         * Reynolds number; this corresponds to the REYNOLDS_NUMBER keyword in the configuration file.
         */
    } else if (index == Math_Problem) {
        *ainame              = EG_strdup("Math_Problem");
        defval->type         = String;
        defval->vals.string  = NULL;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Change;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Direct");

        /*! \page aimInputsSU2
         * - <B> Math_Problem = "Direct"</B> <br>
         * Math problem type; this corresponds to the MATH_PROBLEM keyword in the configuration file.
         * Options: DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, ... see SU2 template for additional options.
         */

    } else if (index == Physical_Problem) {
        *ainame              = EG_strdup("Physical_Problem");
        defval->type         = String;
        defval->vals.string  = NULL;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Change;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Euler");

        /*! \page aimInputsSU2
         * - <B> Physical_Problem = "Euler"</B> <br>
         * Physical problem type; this corresponds to the PHYSICAL_PROBLEM keyword in the configuration file.
         * Options: Euler, Navier_Stokes, Wave_Equation, ... see SU2 template for additional options.
         */

    } else if (index == Equation_Type) {
        *ainame              = EG_strdup("Equation_Type"); // Equation type
        defval->type         = String;
        defval->vals.string  = NULL;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.string  = EG_strdup("Compressible");

        /*! \page aimInputsSU2
         * - <B> Equation_Type = "Compressible"</B> <br>
         * Equation regime type; this corresponds to the REGIME_TYPE keyword in the configuration file.
         * Options: Compressible or Incompressible.
         */

    } else if (index == Turbulence_Model) {
        *ainame              = EG_strdup("Turbulence_Model");
        defval->type         = String;
        defval->vals.string  = NULL;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.string  = EG_strdup("SA_NEG");

        /*! \page aimInputsSU2
         * - <B> Turbulence_Model = "SA_NEG"</B> <br>
         * RANS turbulence model; this corresponds to the KIND_TURB_MODEL keyword in the configuration file.
         * Options: NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP, SST_SUST
         */

    } else if (index == Alpha) {
        *ainame              = EG_strdup("Alpha");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, "degree", aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Alpha = 0.0</B> <br>
         * Angle of attack [degree]; this corresponds to the AoA keyword in the configuration file.
         */
    } else if (index == Beta) {
        *ainame              = EG_strdup("Beta");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, "degree", aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Beta = 0.0</B> <br>
         * Side slip angle [degree]; this corresponds to the SIDESLIP_ANGLE keyword in the configuration file.
         */

    } else if (index == Init_Option) {
        *ainame              = EG_strdup("Init_Option");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Reynolds");

        /*! \page aimInputsSU2
         * - <B> Init_Option = "Reynolds"</B> <br>
         * Init option to choose between Reynolds (default) or thermodynamics quantities for
         * initializing the solution (REYNOLDS, TD_CONDITIONS); this corresponds to the INIT_OPTION keyword in the configuration file.
         */

    } else if (index == Overwrite_CFG) {
        *ainame              = EG_strdup("Overwrite_CFG");
        defval->type         = Boolean;
        defval->vals.integer = (int) true;
        defval->nullVal      = NotNull;

        /*! \page aimInputsSU2
         * - <B> Overwrite_CFG = True</B> <br>
         * Provides permission to overwrite configuration file. If set to False a new configuration file won't be generated.
         */

    } else if (index == Num_Iter) {
        *ainame              = EG_strdup("Num_Iter");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 9999;

        /*! \page aimInputsSU2
         * - <B> Num_Iter = 9999</B> <br>
         * Number of total iterations; this corresponds to the EXT_ITER keyword in the configuration file.
         */

    } else if (index == CFL_Number) {
        *ainame              = EG_strdup("CFL_Number");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 10.0;

        /*! \page aimInputsSU2
         * - <B> CFL_Number = 10.0</B> <br>
         *  Courant–Friedrichs–Lewy number; this corresponds to the CFL_NUMBER keyword in the configuration file.
         */

    } else if (index == Boundary_Condition) {
        *ainame              = EG_strdup("Boundary_Condition");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsSU2
         * - <B>Boundary_Condition = NULL </B> <br>
         * See \ref cfdBoundaryConditions for additional details.
         */
    } else if (index == MultiGrid_Level) {
        *ainame              = EG_strdup("MultiGrid_Level");
        defval->type         = Integer;
        defval->vals.integer = 2;
        defval->units        = NULL;
        defval->dim          = Scalar;

        /*! \page aimInputsSU2
         * - <B> MultiGrid_Level = 2</B> <br>
         *  Number of multi-grid levels; this corresponds to the MGLEVEL keyword in the configuration file.
         */

    } else if (index == Residual_Reduction) {
        *ainame              = EG_strdup("Residual_Reduction");
        defval->type         = Integer;
        defval->vals.integer = 6;
        defval->units        = NULL;
        defval->dim          = Scalar;

        /*! \page aimInputsSU2
         * - <B> Residual_Reduction = 6</B> <br>
         *  Residual reduction (order of magnitude with respect to the initial value);
         *  this corresponds to the RESIDUAL_REDUCTION keyword in the configuration file.
         */

    } else if (index == Unit_System) {
        *ainame              = EG_strdup("Unit_System");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("SI");
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;

        /*! \page aimInputsSU2
         * - <B> Unit_System = "SI"</B> <br>
         * Measurement unit system; this corresponds to the SYSTEM_MEASUREMENTS keyword in the configuration file. See SU2
         * template for additional details.
         */

    } else if (index == Reference_Dimensionalization) {
        *ainame              = EG_strdup("Reference_Dimensionalization");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Dimensional");
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;

        /*! \page aimInputsSU2
         * - <B> Reference_Dimensionalization = NULL</B> <br>
         * Reference dimensionalization; this corresponds to the REF_DIMENSIONALIZATION keyword in the configuration file. See SU2
         * template for additional details.
         */

    } else if (index == Freestream_Pressure) {
        *ainame              = EG_strdup("Freestream_Pressure");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        if (units != NULL && units->pressure != NULL) {
            AIM_STRDUP(defval->units, units->pressure, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Freestream_Pressure = NULL</B> <br>
         * Freestream reference pressure; this corresponds to the FREESTREAM_PRESSURE keyword in the configuration file. See SU2
         * template for additional details.
         */

    } else if (index == Freestream_Temperature) {
        *ainame              = EG_strdup("Freestream_Temperature");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        if (units != NULL && units->temperature != NULL) {
            AIM_STRDUP(defval->units, units->temperature, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Freestream_Temperature = NULL</B> <br>
         * Freestream reference temperature; this corresponds to the FREESTREAM_TEMPERATURE keyword in the configuration file. See SU2
         * template for additional details.
         */

    } else if (index == Freestream_Density) {
        *ainame              = EG_strdup("Freestream_Density");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        if (units != NULL && units->density != NULL) {
            AIM_STRDUP(defval->units, units->density, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Freestream_Density = NULL</B> <br>
         * Freestream reference density; this corresponds to the FREESTREAM_DENSITY keyword in the configuration file. See SU2
         * template for additional details.
         */
    } else if (index == Freestream_Velocity) {
        *ainame              = EG_strdup("Freestream_Velocity");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        if (units != NULL && units->speed != NULL) {
            AIM_STRDUP(defval->units, units->speed, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Freestream_Velocity = NULL</B> <br>
         * Freestream reference velocity; this corresponds to the FREESTREAM_VELOCITY keyword in the configuration file. See SU2
         * template for additional details.
         */
    } else if (index == Freestream_Viscosity) {
        *ainame              = EG_strdup("Freestream_Viscosity");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        if (units != NULL && units->viscosity != NULL) {
            AIM_STRDUP(defval->units, units->viscosity, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B> Freestream_Viscosity = NULL</B> <br>
         * Freestream reference viscosity; this corresponds to the FREESTREAM_VISCOSITY keyword in the configuration file. See SU2
         * template for additional details.
         */
    } else if (index == Moment_Center) {
        *ainame               = EG_strdup("Moment_Center");
        defval->type          = Double;
        defval->dim           = Vector;
        defval->nrow          = 3;
        defval->ncol          = 1;
        AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
        defval->vals.reals[0] = 0.0;
        defval->vals.reals[1] = 0.0;
        defval->vals.reals[2] = 0.0;
        defval->nullVal       = IsNull;
        defval->lfixed        = Fixed;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, units->length, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B>Moment_Center = NULL, [0.0, 0.0, 0.0]</B> <br>
         * Array values correspond to the x_moment_center, y_moment_center, and z_moment_center variables; which correspond
         * to the REF_ORIGIN_MOMENT_X, REF_ORIGIN_MOMENT_Y, and REF_ORIGIN_MOMENT_Z variables respectively in the SU2
         * configuration script.
         * Alternatively, the geometry (body) attributes "capsReferenceX", "capsReferenceY",
         * and "capsReferenceZ" may be used to specify the x-, y-, and z- moment centers, respectively
         * (note: values set through the AIM input will supersede the attribution values).
         */
    } else if (index == Moment_Length) {
        *ainame               = EG_strdup("Moment_Length");
        defval->type          = Double;
        defval->dim           = Scalar;
        defval->nullVal       = IsNull;
        defval->lfixed        = Fixed;
        defval->vals.real     = 1.0;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, units->length, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B>Moment_Length = NULL, 1.0</B> <br>
         * Reference length for pitching, rolling, and yawing non-dimensional; which correspond
         * to the REF_LENGTH_MOMENT. Alternatively, the geometry (body) attribute "capsReferenceSpan" may be
         * used to specify the x-, y-, and z- moment lengths, respectively (note: values set through
         * the AIM input will supersede the attribution values).
         */
    } else if (index == Reference_Area) {
        *ainame              = EG_strdup("Reference_Area");
        defval->type         = Double;
        defval->nullVal      = IsNull;
        defval->lfixed       = Change;
        defval->dim          = Scalar;
        defval->vals.real    = 1.0;
        if (units != NULL && units->area != NULL) {
            AIM_STRDUP(defval->units, units->area, aimInfo, status);
        }

        /*! \page aimInputsSU2
         * - <B>Reference_Area = NULL </B> <br>
         * This sets the reference area for used in force and moment calculations;
         * this corresponds to the REF_AREA keyword in the configuration file.
         * Alternatively, the geometry (body) attribute "capsReferenceArea" maybe used to specify this variable
         * (note: values set through the AIM input will supersede the attribution value).
         */

    } else if (index == Pressure_Scale_Factor) {
        *ainame              = EG_strdup("Pressure_Scale_Factor");
        defval->type         = Double;
        defval->vals.real    = 1.0;
        if (units != NULL && units->pressure != NULL) {
            AIM_STRDUP(defval->units, units->pressure, aimInfo, status);
        }

        if (su2Instance != NULL) su2Instance->pressureScaleFactor = defval;

        /*! \page aimInputsSU2
         * - <B>Pressure_Scale_Factor = 1.0</B> <br>
         * Value to scale Cp or Pressure data when transferring data. Data is scaled based on Pressure = Pressure_Scale_Factor*Cp + Pressure_Scale_Offset.
         */

    } else if (index == Pressure_Scale_Offset) {
        *ainame              = EG_strdup("Pressure_Scale_Offset");
        defval->type         = Double;
        defval->vals.real    = 0.0;
        if (units != NULL && units->pressure != NULL) {
            AIM_STRDUP(defval->units, units->pressure, aimInfo, status);
        }

        if (su2Instance != NULL) su2Instance->pressureScaleOffset = defval;

        /*! \page aimInputsSU2
         * - <B>Pressure_Scale_Offset = 0.0</B> <br>
         * Value to offset Cp or Pressure data when transferring data.
         * Data is scaled based on Pressure = Pressure_Scale_Factor*Cp + Pressure_Scale_Offset.
         */

    } else if (index == Output_Format) {
        *ainame              = EG_strdup("Output_Format");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Paraview");
        defval->lfixed       = Change;

        /*! \page aimInputsSU2
         * - <B> Output_Format = "Paraview"</B> <br>
         * Output file format; this corresponds to the OUTPUT_FORMAT keyword in the configuration file. See SU2
         * template for additional details.
         */

    }  else if (index == Two_Dimensional) {
        *ainame              = EG_strdup("Two_Dimensional");
        defval->type         = Boolean;
        defval->vals.integer = (int) false;

        /*! \page aimInputsSU2
         * - <B>Two_Dimensional = False</B> <br>
         * Run SU2 in 2D mode.
         */

    } else if (index == Convective_Flux) {
        *ainame              = EG_strdup("Convective_Flux");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Roe");
        defval->lfixed       = Change;

        /*! \page aimInputsSU2
         * - <B> Convective_Flux = "Roe"</B> <br>
         * Numerical method for convective (inviscid) flux construction; this corresponds to the CONV_NUM_METHOD_FLOW keyword in the configuration file. See SU2
         * template for additional details.
         */

    } else if (index == SU2_Version) {
        *ainame              = EG_strdup("SU2_Version");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Blackbird");
        defval->lfixed       = Change;

        /*! \page aimInputsSU2
         * - <B>SU2_Version = "Blackbird"</B> <br>
         * SU2 version to generate specific configuration file. Options: "Cardinal(4.0)", "Raven(5.0)", "Falcon(6.2)" or "Blackbird(7.3.0)".
         */

      if (su2Instance != NULL) su2Instance->su2Version = defval;

    } else if (index == Surface_Monitor) {
      *ainame             = EG_strdup("Surface_Monitor");
      defval->type        = String;
      defval->dim         = Vector;
      defval->vals.string = NULL;
      defval->nullVal     = IsNull;
      defval->lfixed      = Change;

      /*! \page aimInputsSU2
       * - <B>Surface_Monitor = NULL</B> <br>
       * Array of surface names where the non-dimensional coefficients are evaluated
       */
    } else if (index == Surface_Deform) {
      *ainame             = EG_strdup("Surface_Deform");
      defval->type        = String;
      defval->dim         = Vector;
      defval->vals.string = NULL;
      defval->nullVal     = IsNull;
      defval->lfixed      = Change;

      /*! \page aimInputsSU2
       * - <B>Surface_Deform = NULL</B> <br>
       * Array of surface names that should be deformed. Defaults to all inviscid and viscous surfaces.
       */

    } else if (index == Input_String) {
      *ainame             = EG_strdup("Input_String");
      defval->type        = String;
      defval->vals.string = NULL;
      defval->nullVal     = IsNull;
      defval->lfixed      = Change;
      defval->dim         = Vector;

      /*! \page aimInputsSU2
       * - <B>Input_String = NULL</B> <br>
       * Array of input strings that will be written as is to the end of the SU2 cfg file.
       */

    } else if (index == Mesh) {
        *ainame             = AIM_NAME(Mesh);
        defval->type        = PointerMesh;
        defval->nrow        = 1;
        defval->lfixed      = Fixed;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;
        AIM_STRDUP(defval->meshWriter, MESHWRITER, aimInfo, status);

        /*! \page aimInputsSU2
         * - <B>Mesh = NULL</B> <br>
         * A Area_Mesh or Volume_Mesh link for 2D and 3D calculations respectively.
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
    // Function return flag
    int status;

    // Python linking
    int pythonLinked = (int) false;
    int usePython = (int) false;

    // Indexing
    int i, body;

    // EGADS return values
    int          atype, alen;
    const int    *ints;
    const char   *string, *intents;
    const double *reals;

    // Output filename
    char meshfilename[PATH_MAX];
    char filepath[PATH_MAX];

    // File pointer if we are going to internally write the namelist
    // FILE *fnml;

    // AIM input bodies
    int  numBody;
    ego *bodies = NULL;

    const char *bodyLunits=NULL;
    cfdUnitsStruct *units=NULL;

    // Boundary/surface properties
    cfdBoundaryConditionStruct bcProps;

    // Boundary conditions container - for writing .mapbc file
    bndCondStruct bndConds;

    // Volume Mesh obtained from meshing AIM
    aimMeshRef *meshRef = NULL;

    // Discrete data transfer variables
    char **boundNames = NULL;
    int numBoundName;

    int attrLevel = 0;
    int withMotion = (int) false;

    aimStorage *su2Instance;

    su2Instance = (aimStorage *) instStore;
    units = &su2Instance->units;

    status = initiate_bndCondStruct(&bndConds);
    if (status != CAPS_SUCCESS) return status;

    // Get boundary conditions
    status = initiate_cfdBoundaryConditionStruct(&bcProps);
    if (status != CAPS_SUCCESS) return status;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" su2AIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif

    if ((numBody <= 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }

    if (aimInputs == NULL) return CAPS_BADVALUE;

    if (units->length != NULL) {
      // Get length units
      status = check_CAPSLength(numBody, bodies, &bodyLunits);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No units assigned *** capsLength is not set in *.csm file!");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    // Get reference quantities from the bodies
    for (body=0; body < numBody; body++) {

        if (aimInputs[Reference_Area-1].nullVal == IsNull) {
            status = EG_attributeRet(bodies[body], "capsReferenceArea", &atype,
                                     &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL) {
                    aimInputs[Reference_Area-1].vals.real = reals[0];
                    aimInputs[Reference_Area-1].nullVal   = NotNull;
                    if (bodyLunits != NULL) {
                        AIM_FREE(aimInputs[Reference_Area-1].units);
                        status = aim_unitMultiply(aimInfo, bodyLunits, bodyLunits, &aimInputs[Reference_Area-1].units);
                        AIM_STATUS(aimInfo, status);
                    }
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceArea should be followed by a single real value!");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }
        }

        if (aimInputs[Moment_Length-1].nullVal == IsNull) {
            status = EG_attributeRet(bodies[body], "capsReferenceSpan", &atype,
                                     &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS){

                if (atype == ATTRREAL) {
                    aimInputs[Moment_Length-1].vals.real = reals[0];
                    aimInputs[Moment_Length-1].nullVal   = NotNull;
                    if (bodyLunits != NULL) {
                        AIM_FREE(aimInputs[Moment_Length-1].units);
                        AIM_STRDUP( aimInputs[Moment_Length-1].units, bodyLunits, aimInfo, status );
                    }
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceSpan should be followed by a single real value!");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }
        }

        if (aimInputs[Moment_Center-1].nullVal == IsNull) {

            status = EG_attributeRet(bodies[body], "capsReferenceX", &atype,
                                     &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL) {
                    aimInputs[Moment_Center-1].vals.reals[0] = reals[0];
                    aimInputs[Moment_Center-1].nullVal       = NotNull;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceX should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[body], "capsReferenceY", &atype,
                                     &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL) {
                    aimInputs[Moment_Center-1].vals.reals[1] = reals[0];
                    aimInputs[Moment_Center-1].nullVal       = NotNull;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceY should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[body], "capsReferenceZ", &atype,
                                     &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL) {
                    aimInputs[Moment_Center-1].vals.reals[2] = reals[0];
                    aimInputs[Moment_Center-1].nullVal       = NotNull;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            if (bodyLunits != NULL) {
                AIM_FREE(aimInputs[Moment_Center-1].units);
                AIM_STRDUP( aimInputs[Moment_Center-1].units, bodyLunits, aimInfo, status );
            }
        }

    }

    if (units->length != NULL) {
      if (aimInputs[Moment_Length-1      ].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Moment_Length      , "Cannot be NULL with unitSys != NULL"); status = CAPS_BADVALUE; goto cleanup;
      }
      if (aimInputs[Reference_Area-1     ].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Reference_Area     , "Cannot be NULL with unitSys != NULL"); status = CAPS_BADVALUE; goto cleanup;
      }
      if (aimInputs[Freestream_Density-1 ].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Freestream_Density , "Cannot be NULL with unitSys != NULL"); status = CAPS_BADVALUE; goto cleanup;
      }
      if (aimInputs[Freestream_Velocity-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Freestream_Velocity, "Cannot be NULL with unitSys != NULL"); status = CAPS_BADVALUE; goto cleanup;
      }
      if (aimInputs[Freestream_Pressure-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Freestream_Pressure, "Cannot be NULL with unitSys != NULL"); status = CAPS_BADVALUE; goto cleanup;
      }

      status = cfd_cfdCoefficientUnits(aimInfo,
                                       aimInputs[Moment_Length-1      ].vals.real, aimInputs[Moment_Length-1      ].units,
                                       aimInputs[Reference_Area-1     ].vals.real, aimInputs[Reference_Area-1     ].units,
                                       aimInputs[Freestream_Density-1 ].vals.real, aimInputs[Freestream_Density-1 ].units,
                                       aimInputs[Freestream_Velocity-1].vals.real, aimInputs[Freestream_Velocity-1].units,
                                       aimInputs[Freestream_Pressure-1].vals.real, aimInputs[Freestream_Pressure-1].units,
                                       units);
      AIM_STATUS(aimInfo, status);
    }

    // Check to see if python was linked
#ifdef HAVE_PYTHON
    pythonLinked = (int) true;
#else
    pythonLinked = (int) false;
#endif

    // Should we use python even if it was linked? - Not implemented
    /*usePython = aimInputs[Use_Python_NML-1].vals.integer;
    if (usePython == (int) true && pythonLinked == (int) false) {

        printf("Use of Python library requested but not linked!\n");
        usePython = (int) false;

    } else if (usePython == (int) false && pythonLinked == (int) true) {

        printf("Python library was linked, but will not be used!\n");
    }
    */

    // Get project name
    su2Instance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Get attribute to index mapping
    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
        su2Instance->groupMap.mapName == NULL) {
        if (aimInputs[Two_Dimensional-1].vals.integer == (int) true) {
          attrLevel = 2; // Only search down to the edge level of the EGADS body
        } else {
          attrLevel = 1; // Only search down to the face level of the EGADS body
        }

        // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
        status = create_CAPSGroupAttrToIndexMap(numBody,
                                                bodies,
                                                attrLevel,
                                                &su2Instance->groupMap);
        if (status != CAPS_SUCCESS) return status;
    }

    if (aimInputs[Boundary_Condition-1].nullVal ==  IsNull) {
      AIM_ANALYSISIN_ERROR(aimInfo, Boundary_Condition, "No boundary conditions provided!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get boundary conditions
    status = cfd_getBoundaryCondition( aimInfo,
                                       aimInputs[Boundary_Condition-1].length,
                                       aimInputs[Boundary_Condition-1].vals.tuple,
                                       &su2Instance->groupMap,
                                       &bcProps);
    AIM_STATUS(aimInfo, status);

    if (aimInputs[Mesh-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Mesh, "'Mesh' input must be linked to an output 'Area_Mesh' or 'Volume_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;
    AIM_NOTNULL(meshRef, aimInfo, status);

    if (status == CAPS_SUCCESS) {

        status = populate_bndCondStruct_from_bcPropsStruct(&bcProps, &bndConds);
        if (status != CAPS_SUCCESS) return status;

        // Replace dummy values in bcVal with SU2 specific values (more to add here)
        printf("Writing boundary flags\n");
        for (i = 0; i < bcProps.numSurfaceProp ; i++) {
            printf(" - bcProps.surfaceProp[%d].surfaceType = %u\n",
                   i, bcProps.surfaceProp[i].surfaceType);

            // {UnknownBoundary, Inviscid, Viscous, Farfield, Freestream,
            //  BackPressure, Symmetry, SubsonicInflow, SubsonicOutflow}


            if      (bcProps.surfaceProp[i].surfaceType == Inviscid       ) bndConds.bcVal[i] = 3000;
            else if (bcProps.surfaceProp[i].surfaceType == Viscous        ) bndConds.bcVal[i] = 4000;
            else if (bcProps.surfaceProp[i].surfaceType == Farfield       ) bndConds.bcVal[i] = 5000;
            else if (bcProps.surfaceProp[i].surfaceType == Extrapolate    ) bndConds.bcVal[i] = 5026;
            else if (bcProps.surfaceProp[i].surfaceType == Freestream     ) bndConds.bcVal[i] = 5050;
            else if (bcProps.surfaceProp[i].surfaceType == BackPressure   ) bndConds.bcVal[i] = 5051;
            else if (bcProps.surfaceProp[i].surfaceType == SubsonicInflow ) bndConds.bcVal[i] = 7011;
            else if (bcProps.surfaceProp[i].surfaceType == SubsonicOutflow) bndConds.bcVal[i] = 7012;
            else if (bcProps.surfaceProp[i].surfaceType == MassflowIn     ) bndConds.bcVal[i] = 7036;
            else if (bcProps.surfaceProp[i].surfaceType == MassflowOut    ) bndConds.bcVal[i] = 7031;
            else if (bcProps.surfaceProp[i].surfaceType == FixedInflow    ) bndConds.bcVal[i] = 7100;
            else if (bcProps.surfaceProp[i].surfaceType == FixedOutflow   ) bndConds.bcVal[i] = 7105;
            else if (bcProps.surfaceProp[i].surfaceType == Symmetry       ) {

                if      (bcProps.surfaceProp[i].symmetryPlane == 1) bndConds.bcVal[i] = 6021;
                else if (bcProps.surfaceProp[i].symmetryPlane == 2) bndConds.bcVal[i] = 6022;
                else if (bcProps.surfaceProp[i].symmetryPlane == 3) bndConds.bcVal[i] = 6023;
            }
        }
        printf("Done boundary flags\n");

        status = destroy_bndCondStruct(&bndConds);
        if (status != CAPS_SUCCESS) return status;
    }

    // If data transfer is ok ....
    withMotion = (int) false;
    if (meshRef->nmap > 0) {
        //See if we have data transfer information
        status = aim_getBounds(aimInfo, &numBoundName, &boundNames);
        if (status == CAPS_SUCCESS && numBoundName > 0) {
            status = su2_dataTransfer(aimInfo,
                                      su2Instance->projectName,
                                      meshRef);
            if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND) goto cleanup;

            withMotion = (int) true;
        }
    } // End if data transfer ok

    status = aim_relPath(aimInfo, meshRef->fileName, ".", filepath);
    AIM_STATUS(aimInfo, status);
    snprintf(meshfilename, PATH_MAX, "%s%s", filepath, MESHEXTENSION);

    if (usePython == (int) true && pythonLinked == (int) true) {
        ///////////////////////////////////////////////////////////////////////////////////
        // Place holder for writing SU2 CFG file using python. Currently NOT implemented //
        ///////////////////////////////////////////////////////////////////////////////////
    } else {

        if (aimInputs[Overwrite_CFG-1].vals.integer == (int) false) {

            printf("Since Python was not linked and/or being used, the \"Overwrite_CFG\" input needs to be set to \"true\" to give");
            printf(" permission to create a new SU2 cfg. SU2 CFG will NOT be updated!!\n");

        } else {

            printf("Warning: The su2 cfg file will be overwritten!\n");

            if (strcasecmp(aimInputs[SU2_Version-1].vals.string, "Cardinal") == 0) {

                status = su2_writeCongfig_Cardinal(aimInfo, aimInputs, meshfilename, bcProps);

            } else if (strcasecmp(aimInputs[SU2_Version-1].vals.string, "Raven") == 0) {

                status = su2_writeCongfig_Raven(aimInfo, aimInputs, meshfilename, bcProps);

            } else if (strcasecmp(aimInputs[SU2_Version-1].vals.string, "Falcon") == 0) {

                status = su2_writeCongfig_Falcon(aimInfo, aimInputs, meshfilename, bcProps, withMotion);

            } else if (strcasecmp(aimInputs[SU2_Version-1].vals.string, "Blackbird") == 0) {

                status = su2_writeCongfig_Blackbird(aimInfo, aimInputs, meshfilename, bcProps, withMotion);

            } else {

                AIM_ERROR(aimInfo, "Unrecognized 'SU2_Version' = %s! Valid choices are Cardinal, Raven, Falcon, or Blackbird.\n",
                          aimInputs[SU2_Version-1].vals.string);
                status = CAPS_BADVALUE;
            }

            if (status != CAPS_SUCCESS){
                goto cleanup;
            }

        }
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(boundNames);

    (void) destroy_cfdBoundaryConditionStruct(&bcProps);
    // (void) destroy_cfdModalAeroelasticStruct(&modalAeroelastic);

    (void) destroy_bndCondStruct(&bndConds);

    return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  return CAPS_SUCCESS;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int index, char **aoname, capsValue *form)
{
	// SU2 Outputs
    /*! \page aimOutputsSU2 AIM Outputs
     * After successful completion, SU2 writes results in various files.
     * The data from these files can be directly viewed, visualized, and or
     * used for further postprocessing.
     *
     * One of the files is ("forces_breakdown.dat") which summarizes
     * convergence including flow properties, numerical parameters, and
     * resulting force and moment values. As an AIM output, this file is
     * parsed for force and moment coefficients, and printed as closing
     * remarks.
     */

    int numOutVars = 9; // Grouped

    #ifdef DEBUG
        printf(" su2AIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
    #endif

    // Total Forces - Pressure + Viscous
    if      (index == 1) *aoname = EG_strdup("CLtot");
    else if (index == 2) *aoname = EG_strdup("CDtot");
    else if (index == 3) *aoname = EG_strdup("CSFtot");
    else if (index == 4) *aoname = EG_strdup("CMXtot");
    else if (index == 5) *aoname = EG_strdup("CMYtot");
    else if (index == 6) *aoname = EG_strdup("CMZtot");
    else if (index == 7) *aoname = EG_strdup("CXtot");
    else if (index == 8) *aoname = EG_strdup("CYtot");
    else if (index == 9) *aoname = EG_strdup("CZtot");

    /*! \page aimOutputsSU2
     * Net Forces - Pressure + Viscous:
     * - <B>CLtot</B> = The lift coefficient.
     * - <B>CDtot</B> = The drag coefficient.
     * - <B>CSFtot</B> = The skin friction coefficient.
     * - <B>CMXtot</B> = The moment coefficient about the x-axis.
     * - <B>CMYtot</B> = The moment coefficient about the y-axis.
     * - <B>CMZtot</B> = The moment coefficient about the z-axis.
     * - <B>CXtot</B> = The force coefficient about the x-axis.
     * - <B>CYtot</B> = The force coefficient about the y-axis.
     * - <B>CZtot</B> = The force coefficient about the z-axis.
     * .
     *
     */

    // Pressure Forces
    else if (index == 1 + numOutVars) *aoname = EG_strdup("CLtot_p");
    else if (index == 2 + numOutVars) *aoname = EG_strdup("CDtot_p");
    else if (index == 3 + numOutVars) *aoname = EG_strdup("CSFtot_p");
    else if (index == 4 + numOutVars) *aoname = EG_strdup("CMXtot_p");
    else if (index == 5 + numOutVars) *aoname = EG_strdup("CMYtot_p");
    else if (index == 6 + numOutVars) *aoname = EG_strdup("CMZtot_p");
    else if (index == 7 + numOutVars) *aoname = EG_strdup("CXtot_p");
    else if (index == 8 + numOutVars) *aoname = EG_strdup("CYtot_p");
    else if (index == 9 + numOutVars) *aoname = EG_strdup("CZtot_p");

    /*! \page aimOutputsSU2
     * Pressure Forces:
     * - <B>CLtot_p</B> = The lift coefficient - pressure contribution only.
     * - <B>CDtot_p</B> = The drag coefficient - pressure contribution only.
     * - <B>CSFtot_p</B> = The skin friction coefficient - pressure contribution only.
     * - <B>CMXtot_p</B> = The moment coefficient about the x-axis - pressure contribution only.
     * - <B>CMYtot_p</B> = The moment coefficient about the y-axis - pressure contribution only.
     * - <B>CMZtot_p</B> = The moment coefficient about the z-axis - pressure contribution only.
     * - <B>CXtot_p</B> = The force coefficient about the x-axis - pressure contribution only.
     * - <B>CYtot_p</B> = The force coefficient about the y-axis - pressure contribution only.
     * - <B>CZtot_p</B> = The force coefficient about the z-axis - pressure contribution only.
     * .
     *
     */

    // Viscous Forces
    else if (index == 1 + 2*numOutVars) *aoname = EG_strdup("CLtot_v");
    else if (index == 2 + 2*numOutVars) *aoname = EG_strdup("CDtot_v");
    else if (index == 3 + 2*numOutVars) *aoname = EG_strdup("CSFtot_v");
    else if (index == 4 + 2*numOutVars) *aoname = EG_strdup("CMXtot_v");
    else if (index == 5 + 2*numOutVars) *aoname = EG_strdup("CMYtot_v");
    else if (index == 6 + 2*numOutVars) *aoname = EG_strdup("CMZtot_v");
    else if (index == 7 + 2*numOutVars) *aoname = EG_strdup("CXtot_v");
    else if (index == 8 + 2*numOutVars) *aoname = EG_strdup("CYtot_v");
    else if (index == 9 + 2*numOutVars) *aoname = EG_strdup("CZtot_v");

    /*! \page aimOutputsSU2
     * Viscous Forces:
     * - <B>CLtot_p</B> = The lift coefficient - viscous contribution only.
     * - <B>CDtot_p</B> = The drag coefficient - viscous contribution only.
     * - <B>CSFtot_p</B> = The skin friction coefficient - viscous contribution only.
     * - <B>CMXtot_p</B> = The moment coefficient about the x-axis - viscous contribution only.
     * - <B>CMYtot_p</B> = The moment coefficient about the y-axis - viscous contribution only.
     * - <B>CMZtot_p</B> = The moment coefficient about the z-axis - viscous contribution only.
     * - <B>CXtot_p</B> = The force coefficient about the x-axis - viscous contribution only.
     * - <B>CYtot_p</B> = The force coefficient about the y-axis - viscous contribution only.
     * - <B>CZtot_p</B> = The force coefficient about the z-axis - viscous contribution only.
     * .
     *
     */

    //if (index <= 3*numOutVars) {
    form->type    = Double;
    form->dim     = Scalar;
    //}

    return CAPS_SUCCESS;
}


// Calculate SU2 output
int aimCalcOutput(void *instStore, void *aimInfo, int index, capsValue *val)
{
    int status = CAPS_SUCCESS;

    char *strKeyword = NULL; //, *bndSectionKeyword = NULL; // Keyword strings

    char *strValue;    // Holder string for value of keyword

    char *filename = NULL; // File to open
    char fileExtension[] = ".dat";
    char fileSuffix[] = "forces_breakdown_";

    size_t linecap = 0;
    char *line = NULL; // Temporary line holder

    FILE *fp; // File pointer

    int numOutVars = 9; // Grouped

    cfdUnitsStruct *units=NULL;

    capsValue *SurfaceMonitor = NULL;
    aimStorage *su2Instance;

#ifdef DEBUG
    printf(" su2AIM/aimCalcOutput  index = %d!\n", index);
#endif
    su2Instance = (aimStorage *) instStore;
    units = &su2Instance->units;

    val->vals.real = 0.0; // Set default value

    //printf("index = %d\n",index);

    // Set the "search" string(s)
//    if (index <= numOutVars) {
//        bndSectionKeyword  = (char *) "Forces";
//    }

    if      (index == 1 || index == 1 + numOutVars || index == 1 + 2*numOutVars) strKeyword = (char *) "CL:";
    else if (index == 2 || index == 2 + numOutVars || index == 2 + 2*numOutVars) strKeyword = (char *) "CD:";
    else if (index == 3 || index == 3 + numOutVars || index == 3 + 2*numOutVars) strKeyword = (char *) "CSF:";
    else if (index == 4 || index == 4 + numOutVars || index == 4 + 2*numOutVars) strKeyword = (char *) "CMx:";
    else if (index == 5 || index == 5 + numOutVars || index == 5 + 2*numOutVars) strKeyword = (char *) "CMy:";
    else if (index == 6 || index == 6 + numOutVars || index == 6 + 2*numOutVars) strKeyword = (char *) "CMz:";
    else if (index == 7 || index == 7 + numOutVars || index == 7 + 2*numOutVars) strKeyword = (char *) "CFx:";
    else if (index == 8 || index == 8 + numOutVars || index == 8 + 2*numOutVars) strKeyword = (char *) "CFy:";
    else if (index == 9 || index == 9 + numOutVars || index == 9 + 2*numOutVars) strKeyword = (char *) "CFz:";
    else {
        printf("Unrecognized output variable index - %d\n", index);
        return CAPS_BADINDEX;
    }

    status = aim_getValue(aimInfo, Surface_Monitor, ANALYSISIN, &SurfaceMonitor);
    AIM_STATUS(aimInfo, status);

    if (SurfaceMonitor == NULL) {
        AIM_ERROR(aimInfo, "Forces are not available because 'Surface_Monitor' is NULL.");
        return CAPS_BADVALUE;
    }

    if (SurfaceMonitor->nullVal == IsNull) {
        AIM_ERROR(aimInfo, "Forces are not available because 'Surface_Monitor' is not specified.\n");
        return CAPS_BADVALUE;
    }

    // Open SU2 force file
    filename = (char *) EG_alloc((strlen(su2Instance->projectName) +
                                  strlen(fileSuffix) + strlen(fileExtension) +1)*sizeof(char));
    if (filename == NULL) return EGADS_MALLOC;

/*@-bufferoverflowhigh@*/
    sprintf(filename, "%s%s%s", fileSuffix, su2Instance->projectName, fileExtension);
/*@+bufferoverflowhigh@*/

    fp = aim_fopen(aimInfo, filename, "r");
    AIM_FREE(filename); // Free filename allocation

    if (fp == NULL) {
#ifdef DEBUG
        printf(" su2AIM/aimCalcOutput Cannot open Output file!\n");
#endif
        return CAPS_IOERR;
    } else {
#ifdef DEBUG
        printf(" Found force file: %s%s%s",
               fileSuffix, su2Instance->projectName, fileExtension);
#endif
    }

    // Scan the file for the string
    strValue = NULL;
    while( !feof(fp) ) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

/*@-nullpass@*/
        strValue = strstr(line, strKeyword); // Totals
/*@+nullpass@*/
        if( strValue != NULL) {

            if ( index <= numOutVars ) { // Totals
              strValue += strlen(strKeyword);
            } else if ( index <= 2*numOutVars ) { // Pressure
                strValue = strstr(line, "Pressure (");
                strValue = strstr(strValue, "):");
                strValue += 2;
            } else if ( index <= 3*numOutVars ) { // Friction
                strValue = strstr(line, "Friction (");
                strValue = strstr(strValue, "):");
                strValue += 2;
            }

            sscanf(strValue, "%lf", &val->vals.real);
            val->nullVal = NotNull;
            break;
        }
    }

    EG_free(line);
    fclose(fp);

#ifdef DEBUG
    printf(" Done\n Closing force file.\n");
#endif

    if (strValue == NULL) {
        AIM_ERROR(aimInfo, " su2AIM/aimCalcOutput Cannot find %s in Output file!\n",
               strKeyword);
        return CAPS_NOTFOUND;
    }

    // assign units now because they are not know until after preAnalysis
    if (units->length != NULL) {
      AIM_FREE(val->units);
      if (strncmp(strKeyword, "CM", 2) == 0) {
        AIM_STRDUP(val->units, units->Cmoment, aimInfo, status);
      } else {
        AIM_STRDUP(val->units, units->Cforce, aimInfo, status);
      }
    }

    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(filename);

    return status;
}


void aimCleanup(void *instStore)
{
    aimStorage *su2Instance;

#ifdef DEBUG
    printf(" su2AIM/aimCleanup!\n");
#endif
    su2Instance = (aimStorage *) instStore;

    // Clean up su2Instance data
    destroy_mapAttrToIndexStruct(&su2Instance->groupMap);

    // SU2 project name
    su2Instance->projectName = NULL;

    // Pointer to caps input value for scaling pressure during data transfer
    su2Instance->pressureScaleFactor = NULL;

    // Pointer to caps input value for offset pressure during data transfer
    su2Instance->pressureScaleOffset = NULL;

    // Cleanup units
    destroy_cfdUnitsStruct(&su2Instance->units);

    EG_free(su2Instance);

}


/************************************************************************/
// CAPS transferring functions
void aimFreeDiscrPtr(void *ptr)
{
    AIM_FREE(ptr);
}


/************************************************************************/
int aimDiscr(char *tname, capsDiscr *discr)
{
    int i; // Indexing

    int status; // Function return status

    int numBody;
    aimStorage *su2Instance;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;
    capsValue *meshVal;

    // Volume Mesh obtained from meshing AIM
    aimMeshRef *meshRef;

#ifdef DEBUG
    printf(" su2AIM/aimDiscr: tname = %s!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;
    su2Instance = (aimStorage *) discr->instStore;

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    AIM_STATUS(discr->aInfo, status);

    if (bodies == NULL) {
        AIM_ERROR(discr->aInfo, " su2AIM/aimDiscr: No Bodies!\n");
        return CAPS_NOBODIES;
    }

    // Get the mesh input Value
    status = aim_getValue(discr->aInfo, Mesh, ANALYSISIN, &meshVal);
    AIM_STATUS(discr->aInfo, status);

    if (meshVal->nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(discr->aInfo, Mesh, "'Mesh' input must be linked to an output 'Area_Mesh' or 'Volume_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    meshRef = (aimMeshRef *)meshVal->vals.AIMptr;
    AIM_NOTNULL(meshRef, discr->aInfo, status);

    if (meshRef->nmap == 0 || meshRef->maps == NULL) {
        AIM_ERROR(discr->aInfo, "No surface mesh map in volume mesh - data transfer isn't possible.\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aim_newGeometry(discr->aInfo) == CAPS_SUCCESS &&
        su2Instance->groupMap.numAttribute == 0) {
        // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
        status = create_CAPSGroupAttrToIndexMap(numBody,
                                                bodies,
                                                1, // Only search down to the face level of the EGADS body
                                                &su2Instance->groupMap);
        AIM_STATUS(discr->aInfo, status);
    }

    // Do we have an individual surface mesh for each body
    if (meshRef->nmap != numBody) {
        AIM_ERROR(  discr->aInfo, "Number of surface mesh in the linked volume mesh (%d) does not match the number");
        AIM_ADDLINE(discr->aInfo,"of bodies (%d) - data transfer is NOT possible.", meshRef->nmap,numBody);
        status = CAPS_MISMATCH;
        goto cleanup;
    }

    // To this point is doesn't appear that the volume mesh has done anything bad to our surface mesh(es)

    // Store away the tessellation now
    AIM_ALLOC(tess, numBody, ego, discr->aInfo, status);
    for (i = 0; i < numBody; i++) {
        tess[i] = meshRef->maps[i].tess;
    }

    status = mesh_fillDiscr(tname, &su2Instance->groupMap, numBody, tess, discr);
    AIM_STATUS(discr->aInfo, status);

#ifdef DEBUG
    printf(" su2AIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
#endif

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(tess);

    return status;
}


int aimLocateElement(capsDiscr *discr, double *params, double *param,
                     int *bIndex, int *eIndex, double *bary)
{
    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


int aimTransfer(capsDiscr *discr, const char *dataName, int numPoint,
                int dataRank, double *dataVal, /*@unused@*/ char **units)
{
    /*! \page dataTransferSU2 AIM Data Transfer
     *
     * The SU2 AIM has the ability to transfer surface data (e.g. pressure distributions) to and from the AIM
     * using the conservative and interpolative data transfer schemes in CAPS. Currently these transfers may only
     * take place on triangular meshes.
     *
     * \section dataFromSU2 Data transfer from SU2 (FieldOut)
     *
     * <ul>
     *  <li> <B>"Cp", or "CoefficientOfPressure"</B> </li> <br>
     *  Loads the coefficient of pressure distribution from surface_flow_[project_name].cvs file.
     *  This distribution may be scaled based on
     *  Pressure = Pressure_Scale_Factor*Cp + Pressure_Scale_Offset, where "Pressure_Scale_Factor"
     *  and "Pressure_Scale_Offset" are AIM inputs (\ref aimInputsSU2)
     *
     * <li> <B>"Pressure" or "P" </B> </li> <br>
     *  Loads the pressure distribution from surface_flow_[project_name].cvs file.
     *  This distribution may be scaled based on
     *  Pressure = Pressure_Scale_Factor*Pressure + Pressure_Scale_Offset, where "Pressure_Scale_Factor"
     *  and "Pressure_Scale_Offset" are AIM inputs (\ref aimInputsSU2)
     * </ul>
     *
     */ // Reset of this block comes from su2Util.c

    int status, status2; // Function return status
    int i, j, dataPoint, bIndex; // Indexing
    aimStorage *su2Instance;

    // Aero-Load data variables
    int numVariable=0;
    int numDataPoint=0;
    char **variableName = NULL;
    double **dataMatrix = NULL;
    double dataScaleFactor = 1.0;
    double dataScaleOffset = 0.0;
    //char *dataUnits;

    // Indexing in data variables
    int globalIDIndex = -99;
    int variableIndex = -99;

    // Variables used in global node mapping
    //int *storage;
    int globalNodeID;
    int found = (int) false;

    // Filename stuff
    char *filename = NULL; //"pyCAPS_SU2_Tetgen_ddfdrive_bndry1.dat";

#ifdef DEBUG
    printf(" su2AIM/aimTransfer name = %s   npts = %d/%d!\n",
           dataName, numPoint, dataRank);
#endif

    if (strcasecmp(dataName, "Pressure") != 0 &&
        strcasecmp(dataName, "P")        != 0 &&
        strcasecmp(dataName, "Cp")       != 0 &&
        strcasecmp(dataName, "CoefficientOfPressure") != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    su2Instance = (aimStorage *) discr->instStore;

    //Get the appropriate parts of the tessellation to data
    //storage = (int *) discr->ptrm;

    // Zero out data
    for (i = 0; i < numPoint; i++) {
        for (j = 0; j < dataRank; j++ ) {
            dataVal[dataRank*i+j] = 0;
        }
    }

    AIM_ALLOC(filename, (strlen(su2Instance->projectName) +
                         strlen("surface_flow_.csv") + 1), char, discr->aInfo, status);

/*@-bufferoverflowhigh@*/
    sprintf(filename,"%s%s%s", "surface_flow_", su2Instance->projectName, ".csv");
/*@+bufferoverflowhigh@*/

    status = su2_readAeroLoad(discr->aInfo,
                              filename,
                              &numVariable,
                              &variableName,
                              &numDataPoint,
                              &dataMatrix);
    AIM_STATUS(discr->aInfo, status);
    AIM_FREE(filename);
    if (variableName == NULL) {
        AIM_ERROR(discr->aInfo, "NULL variableName!");
        return CAPS_NULLNAME;
    }

    //printf("Number of variables %d\n", numVariable);
    // Output some of the first row of the data
    //for (i = 0; i < numVariable; i++) printf("Variable %d - %.6f\n", i, dataMatrix[i][0]);

    // Loop through the variable list to find which one is the global node ID variable
    for (i = 0; i < numVariable; i++) {
        if (strcasecmp("Global_Index", variableName[i]) == 0 ||
            strcasecmp("PointID"     , variableName[i]) == 0) {
            globalIDIndex = i;
            break;
        }
    }

    if (globalIDIndex == -99) {
        AIM_ERROR(discr->aInfo, "Global node number variable not found in data file\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Loop through the variable list to see if we can find the transfer data name
    for (i = 0; i < numVariable; i++) {

        if (strcasecmp(dataName, "Pressure") == 0 ||
            strcasecmp(dataName, "P")        == 0) {

            if (dataRank != 1) {
                AIM_ERROR(discr->aInfo, "Data transfer rank should be 1 not %d\n", dataRank);
                status = CAPS_BADRANK;
                goto cleanup;
            }

            dataScaleFactor = su2Instance->pressureScaleFactor->vals.real;
            dataScaleOffset = su2Instance->pressureScaleOffset->vals.real;

            if (strcasecmp("Pressure", variableName[i]) == 0) {
                variableIndex = i;
                if (su2Instance->units.pressure != NULL)
                    AIM_STRDUP(*units, su2Instance->units.pressure, discr->aInfo, status);
                break;
            }
        }

        if (strcasecmp(dataName, "Cp")       == 0 ||
            strcasecmp(dataName, "CoefficientOfPressure") == 0) {

            if (dataRank != 1) {
                AIM_ERROR(discr->aInfo, "Data transfer rank should be 1 not %d\n", dataRank);
                status = CAPS_BADRANK;
                goto cleanup;
            }

            dataScaleFactor = su2Instance->pressureScaleFactor->vals.real;
            dataScaleOffset = su2Instance->pressureScaleOffset->vals.real;

            if (strcasecmp(su2Instance->su2Version->vals.string, "Cardinal") == 0) {
                if (strcasecmp("C<sub>p</sub>", variableName[i]) == 0) {
                    variableIndex = i;
                    if (su2Instance->units.Cpressure != NULL)
                        AIM_STRDUP(*units, su2Instance->units.Cpressure, discr->aInfo, status);
                    break;
                }
            } else {
                if (strcasecmp("Pressure_Coefficient", variableName[i]) == 0) {
                    variableIndex = i;
                    if (su2Instance->units.Cpressure != NULL)
                        AIM_STRDUP(*units, su2Instance->units.Cpressure, discr->aInfo, status);
                    break;
                }
            }
        }
    }

    if (variableIndex == -99) {
        AIM_ERROR(discr->aInfo, "Variable %s not found in data file", dataName);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }
    if (dataMatrix == NULL) {
        AIM_ERROR(discr->aInfo, " dataMatrix is NULL!\n");
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    for (i = 0; i < numPoint; i++) {

        bIndex       = discr->tessGlobal[2*i  ];
        globalNodeID = discr->tessGlobal[2*i+1] +
                       discr->bodys[bIndex-1].globalOffset;

        found = (int) false;
        for (dataPoint = 0; dataPoint < numDataPoint; dataPoint++) {
            if ((int) dataMatrix[globalIDIndex][dataPoint] +1 ==  globalNodeID) { // SU2 meshes are 0-index
                found = (int) true;
                break;
            }
        }

        if (found == (int) true) {
            for (j = 0; j < dataRank; j++) {

                // Add something for units - aim_covert()

                dataVal[dataRank*i+j] = dataMatrix[variableIndex][dataPoint]*dataScaleFactor + dataScaleOffset;
                //dataVal[dataRank*i+j] = 99;

            }
        } else {
            AIM_ERROR(discr->aInfo, "Error: Unable to find node %d!\n", globalNodeID);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    // Free data matrix
    if (dataMatrix != NULL) {
        for (i = 0; i < numVariable; i++) {
            AIM_FREE(dataMatrix[i]);
        }

        AIM_FREE(dataMatrix);
    }

    // Free variable list
    status2 = string_freeArray(numVariable, &variableName);
    if (status2 != CAPS_SUCCESS) return status2;

    AIM_FREE(filename);

    return status;
}


int aimInterpolation(capsDiscr *discr, /*@unused@*/ const char *name,
                     int bIndex, int eIndex, double *bary, int rank,
                     double *data, double *result)
{
#ifdef DEBUG
    printf(" su2AIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


int aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                      int bIndex, int eIndex, double *bary, int rank,
                      double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" su2AIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


int aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                   int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" su2AIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                    int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" su2AIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
