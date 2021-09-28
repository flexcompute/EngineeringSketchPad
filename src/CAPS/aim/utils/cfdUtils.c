// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// CFD analysis related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include <string.h>
#include "capsTypes.h"  // Bring in CAPS types
#include "cfdTypes.h"   // Bring in cfd structures
#include "cfdUtils.h"
#include "miscUtils.h"  // Bring in miscellaneous utilities
#include "aimUtil.h"      // Bring in AIM utils

#include <math.h>

#ifdef WIN32
#define strcasecmp  stricmp
#endif

// Fill bcProps in a cfdBCStruct format with boundary condition information from incoming BC Tuple
int cfd_getBoundaryCondition(void *aimInfo,
                             int numTuple,
                             capsTuple bcTuple[],
                             mapAttrToIndexStruct *attrMap,
                             cfdBoundaryConditionStruct *bcProps) {

    /*! \page cfdBoundaryConditions CFD Boundary Conditions
     * Structure for the boundary condition tuple  = ("CAPS Group Name", "Value").
     * "CAPS Group Name" defines the capsGroup on which the boundary condition should be applied.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringCFDBoundary) or a single string keyword string
     *	(see Section \ref keyStringCFDBoundary)
     */

    int status; //Function return

    int i, j; // Indexing

    int found = (int) false;

    int bcIndex;

    char *keyValue = NULL;
    char *keyWord = NULL;

    // Destroy our bcProps structures coming in if aren't 0 and NULL already
    status = destroy_cfdBoundaryConditionStruct(bcProps);
    AIM_STATUS(aimInfo, status);

    printf("Getting CFD boundary conditions\n");

    //bcProps->numSurfaceProp = attrMap->numAttribute;
    bcProps->numSurfaceProp = numTuple;

    if (bcProps->numSurfaceProp > 0) {
        bcProps->surfaceProp = (cfdSurfaceStruct *) EG_alloc(bcProps->numSurfaceProp * sizeof(cfdSurfaceStruct));
        if (bcProps->surfaceProp == NULL) return EGADS_MALLOC;

    } else {
        printf("\tWarning: Number of Boundary Conditions is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < bcProps->numSurfaceProp; i++) {
        status = initiate_cfdSurfaceStruct(&bcProps->surfaceProp[i]);
        AIM_STATUS(aimInfo, status);
    }

    //printf("Number of tuple pairs = %d\n", numTuple);

    for (i = 0; i < numTuple; i++) {

        printf("\tBoundary condition name - %s\n", bcTuple[i].name);

        status = get_mapAttrToIndexIndex(attrMap, (const char *) bcTuple[i].name, &bcIndex);
        if (status == CAPS_NOTFOUND) {
            bcIndex = aim_getIndex(aimInfo, "Boundary_Condition", ANALYSISIN);
            AIM_ANALYSISIN_ERROR(aimInfo, bcIndex, "BC name \"%s\" not found in capsGroup attributes", bcTuple[i].name);
            AIM_ADDLINE(aimInfo, "Available capsGroup attributes:");
            for (j = 0; j < attrMap->numAttribute; j++) {
                AIM_ADDLINE(aimInfo, "\t%s", attrMap->attributeName[j]);
            }
            return status;
        }

        // Replaced bcIndex with i
        // bcIndex is 1 bias
        bcProps->surfaceProp[i].bcID = bcIndex;

        // bcIndex is 1 bias coming from attribute mapper
        //bcIndex = bcIndex -1;

        // Copy boundary condition name
        AIM_FREE(bcProps->surfaceProp[i].name);

        AIM_STRDUP(bcProps->surfaceProp[i].name, bcTuple[i].name, aimInfo, status);

        // Do we have a json string?
        if (strncmp(bcTuple[i].value, "{", 1) == 0) {

            //printf("\t\tJSON String - %s\n",bcTuple[i].value);

            /*! \page cfdBoundaryConditions
             * \section jsonStringCFDBoundary JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary
             * \if (FUN3D || SU2)
             *  (eg. "Value" = {"bcType": "Viscous", "wallTemperature": 1.1})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * <ul>
             * <li> <B>bcType = "Inviscid"</B> </li> <br>
             *      Boundary condition type. Options:
             *
             * \if (FUN3D)
             *  - Inviscid
             *  - Viscous
             *  - Farfield
             *  - Extrapolate
             *  - Freestream
             *  - BackPressure
             *  - Symmetry
             *  - SubsonicInflow
             *  - SubsonicOutflow
             *  - MassflowIn
             *  - MassflowOut
             *  - MachOutflow
             *
             * \elseif (SU2)
             *  - Inviscid
             *  - Viscous
             *  - Farfield
             *  - Freestream
             *  - BackPressure
             *  - Symmetry
             *  - SubsonicInflow
             *  - SubsonicOutflow

             * \else
             *  - Inviscid
             *  - Viscous
             *  - Farfield
             *  - Extrapolate
             *  - Freestream
             *  - BackPressure
             *  - Symmetry
             *  - SubsonicInflow
             *  - SubsonicOutflow
             *  - MassflowIn
             *  - MassflowOut
             *  - MachOutflow
             *  - FixedInflow
             *  - FixedOutflow
             *
             * \endif
             * </ul>
             */

            // Get property Type
            keyWord = "bcType";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                //{UnknownBoundary, Inviscid, Viscous, Farfield, Extrapolate, Freestream,
                // BackPressure, Symmetry, SubsonicInflow, SubsonicOutflow,
                // MassflowIn, MassflowOut, FixedInflow, FixedOutflow, MachOutflow}
                if      (strcasecmp(keyValue, "\"Inviscid\"")        == 0) bcProps->surfaceProp[i].surfaceType = Inviscid;
                else if (strcasecmp(keyValue, "\"Viscous\"")         == 0) bcProps->surfaceProp[i].surfaceType = Viscous;
                else if (strcasecmp(keyValue, "\"Farfield\"")        == 0) bcProps->surfaceProp[i].surfaceType = Farfield;
                else if (strcasecmp(keyValue, "\"Extrapolate\"")     == 0) bcProps->surfaceProp[i].surfaceType = Extrapolate;
                else if (strcasecmp(keyValue, "\"Freestream\"")      == 0) bcProps->surfaceProp[i].surfaceType = Freestream;
                else if (strcasecmp(keyValue, "\"BackPressure\"")    == 0) bcProps->surfaceProp[i].surfaceType = BackPressure;
                else if (strcasecmp(keyValue, "\"Symmetry\"")        == 0) bcProps->surfaceProp[i].surfaceType = Symmetry;
                else if (strcasecmp(keyValue, "\"SubsonicInflow\"")  == 0) bcProps->surfaceProp[i].surfaceType = SubsonicInflow;
                else if (strcasecmp(keyValue, "\"SubsonicOutflow\"") == 0) bcProps->surfaceProp[i].surfaceType = SubsonicOutflow;
                else if (strcasecmp(keyValue, "\"MassflowIn\"")      == 0) bcProps->surfaceProp[i].surfaceType = MassflowIn;
                else if (strcasecmp(keyValue, "\"MassflowOut\"")     == 0) bcProps->surfaceProp[i].surfaceType = MassflowOut;
                else if (strcasecmp(keyValue, "\"MachOutflow\"")     == 0) bcProps->surfaceProp[i].surfaceType = MachOutflow;
                else if (strcasecmp(keyValue, "\"FixedInflow\"")     == 0) bcProps->surfaceProp[i].surfaceType = FixedInflow;
                else if (strcasecmp(keyValue, "\"FixedOutflow\"")    == 0) bcProps->surfaceProp[i].surfaceType = FixedOutflow;
                else {

                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Boundary_Condition tuple %s, current options (not all options "
                            "are valid for this analysis tool - see AIM documentation) are "
                            "\" Inviscid, Viscous, Farfield, Extrapolate, Freestream, BackPressure, Symmetry, "
                            "SubsonicInflow, SubsonicOutflow, MassflowIn, MassflowOut, MachOutflow, "
                            "FixedInflow, FixedOutflow"
                            "\"",
                            keyWord, keyValue, bcTuple[i].name);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else {
                AIM_ERROR(aimInfo, "\tNo \"%s\" specified for tuple %s in json string", keyWord,
                        bcTuple[i].name);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            AIM_FREE(keyValue);

            // Wall specific properties

            /*! \page cfdBoundaryConditions
             *  \subsection cfdBoundaryWallProp Wall Properties
             *
             * \if FUN3D
             *  <ul>
             *	<li> <B>wallTemperature = 0.0</B> </li> <br>
             *  The ratio of wall temperature to reference temperature for inviscid and viscous surfaces. Adiabatic wall = -1
             *  </ul>
             * \elseif SU2
             *
             *  <ul>
             *  <li> <B>wallTemperature = 0.0</B> </li> <br>
             *  Dimensional wall temperature for inviscid and viscous surfaces
             *  </ul>
             *
             * \else
             * <ul>
             * <li>  <B>wallTemperature = 0.0</B> </li> <br>
             *  Temperature on inviscid and viscous surfaces.
             * </ul>
             * \endif
             */
            keyWord = "wallTemperature";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                bcProps->surfaceProp[i].wallTemperatureFlag = (int) true;

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].wallTemperature);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * \elseif SU2
             *
             * <ul>
             * <li>  <B>wallHeatFlux = 0.0</B> </li> <br>
             *  Heat flux on viscous surfaces.
             * </ul>
             *
             * \else
             * <ul>
             * <li>  <B>wallHeatFlux = 0.0</B> </li> <br>
             *  Heat flux on inviscid and viscous surfaces.
             * </ul>
             * \endif
             */

            keyWord = "wallHeatFlux";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                bcProps->surfaceProp[i].wallTemperatureFlag = (int) true;
                bcProps->surfaceProp[i].wallTemperature = -10;

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].wallHeatFlux);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            // Stagnation quantities

            /*! \page cfdBoundaryConditions
             * \subsection cfdBoundaryStagProp Stagnation Properties
             *
             *
             * \if FUN3D
             * <ul>
             *  <li>  <B>totalPressure = 0.0</B> </li> <br>
             *  Ratio of total pressure to reference pressure on a boundary surface.
             *
             * </ul>
             *
             * \elseif SU2
             * <ul>
             * <li>  <B>totalPressure = 0.0</B> </li> <br>
             *  Dimensional total pressure on a boundary surface.
             * </ul>
             *
             * \else
             * <ul>
             * <li>  <B>totalPressure = 0.0</B> </li> <br>
             *  Total pressure on a boundary surface.
             * </ul>
             * \endif
             *
             */
            keyWord = "totalPressure";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].totalPressure);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * <ul>
             * <li>  <B>totalTemperature = 0.0</B> </li> <br>
             *  Ratio of total temperature to reference temperature on a boundary surface.
             * </ul>
             *
             * \elseif SU2
             * <ul>
             * <li>  <B>totalTemperature = 0.0</B> </li> <br>
             *  Dimensional total temperature on a boundary surface.
             * </ul>
             *
             * \else
             * <ul>
             * <li>  <B>totalTemperature = 0.0</B> </li> <br>
             *  Total temperature on a boundary surface.
             * </ul>
             * \endif
             *
             */
            keyWord = "totalTemperature";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].totalTemperature);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             * <ul>
             *  <li> <B>totalDensity = 0.0</B> </li> <br>
             *  Total density of boundary.
             * </ul>
             * \endif
             *
             */
            keyWord = "totalDensity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].totalDensity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            // Static quantities

            /*! \page cfdBoundaryConditions
             * \subsection cfdBoundaryStaticProp Static Properties
             *
             * \if FUN3D
             *  <ul>
             *  <li> <B>staticPressure = 0.0</B> </li> <br>
             *  Ratio of static pressure to reference pressure on a boundary surface.
             *  </ul>
             * \elseif SU2
             *  <ul>
             *  <li> <B>staticPressure = 0.0</B> </li> <br>
             *  Dimensional static pressure on a boundary surface.
             *  </ul>
             * \else
             * <ul>
             * <li> <B>staticPressure = 0.0</B> </li> <br>
             *  Static pressure of boundary.
             * </ul>
             * \endif
             */
            keyWord = "staticPressure";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].staticPressure);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li> <B>staticTemperature = 0.0</B> </li> <br>
             *  Static temperature of boundary.
             *  </ul>
             * \endif
             */
            keyWord = "staticTemperature";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].staticTemperature);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li>  <B>staticDensity = 0.0</B> </li> <br>
             *  Static density of boundary.
             *  </ul>
             * \endif
             */
            keyWord = "staticDensity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].staticDensity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }


            // Velocity components
            /*! \page cfdBoundaryConditions
             * \subsection cfdBoundaryVelocity Velocity Components
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             * <ul>
             *  <li>  <B>uVelocity = 0.0</B> </li> <br>
             *  X-velocity component on boundary.
             * </ul>
             * \endif
             */
            keyWord = "uVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].uVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li>  <B>vVelocity = 0.0</B> </li> <br>
             *  Y-velocity component on boundary.
             *  </ul>
             * \endif
             *
             */
            keyWord = "vVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].vVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li>  <B>wVelocity = 0.0</B> </li> <br>
             *  Z-velocity component on boundary.
             * </ul>
             * \endif
             */
            keyWord = "wVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].wVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cfdBoundaryConditions
             *
             * \if FUN3D
             * <ul>
             *  <li> <B>machNumber = 0.0</B> </li> <br>
             *  Mach number on boundary.
             * </ul>
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li> <B>machNumber = 0.0</B> </li> <br>
             *  Mach number on boundary.
             * </ul>
             * \endif
             */
            keyWord = "machNumber";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].machNumber);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            // Massflow
            /*! \page cfdBoundaryConditions
             * \subsection cfdBoundaryMassflow Massflow Properties
             *
             * \if FUN3D
             *  <ul>
             *	<li> <B>massflow = 0.0</B> </li> <br>
             *  Massflow through the boundary in units of grid units squared.
             *  </ul>
             * \elseif SU2
             *
             * \else
             *  <ul>
             *  <li> <B>massflow = 0.0</B> </li> <br>
             *  Massflow through the boundary.
             *  </ul>
             * \endif
             */
            keyWord = "massflow";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[i].massflow);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

        } else {
            //printf("\t\tkeyString value - %s\n", bcTuple[i].value);

            /*! \page cfdBoundaryConditions
             * \section keyStringCFDBoundary Single Value String
             *
             * If "Value" is a single string the following options maybe used:
             *
             * \if (SU2)
             * - "Inviscid" (default)
             * - "Viscous"
             * - "Farfield"
             * - "Freestream"
             * - "SymmetryX"
             * - "SymmetryY"
             * - "SymmetryZ"
             * \else
             * - "Inviscid" (default)
             * - "Viscous"
             * - "Farfield"
             * - "Extrapolate"
             * - "Freestream"
             * - "SymmetryX"
             * - "SymmetryY"
             * - "SymmetryZ"
             * \endif
             * */
            //{UnknownBoundary, Inviscid, Viscous, Farfield, Extrapolate, Freestream,
            // BackPressure, Symmetry, SubsonicInflow, SubsonicOutflow,
            // MassflowIn, MassflowOut, FixedInflow, FixedOutflow, MachOutflow}

            keyValue = string_removeQuotation(bcTuple[i].value);
            if (keyValue == NULL) {
              status = EGADS_MALLOC;
              goto cleanup;
            }

            if      (strcasecmp(keyValue, "Inviscid" ) == 0) bcProps->surfaceProp[i].surfaceType = Inviscid;
            else if (strcasecmp(keyValue, "Viscous"  ) == 0) bcProps->surfaceProp[i].surfaceType = Viscous;
            else if (strcasecmp(keyValue, "Farfield" ) == 0) bcProps->surfaceProp[i].surfaceType = Farfield;
            else if (strcasecmp(keyValue, "Extrapolate" ) == 0) bcProps->surfaceProp[i].surfaceType = Extrapolate;
            else if (strcasecmp(keyValue, "Freestream" ) == 0)  bcProps->surfaceProp[i].surfaceType = Freestream;
            else if (strcasecmp(keyValue, "SymmetryX") == 0) {
                bcProps->surfaceProp[i].surfaceType = Symmetry;
                bcProps->surfaceProp[i].symmetryPlane = 1;
            }
            else if (strcasecmp(keyValue, "SymmetryY") == 0) {
                bcProps->surfaceProp[i].surfaceType = Symmetry;
                bcProps->surfaceProp[i].symmetryPlane = 2;
            }
            else if (strcasecmp(keyValue, "SymmetryZ") == 0) {
                bcProps->surfaceProp[i].surfaceType = Symmetry;
                bcProps->surfaceProp[i].symmetryPlane = 3;
            }
            else {
                AIM_ERROR(aimInfo, "\tUnrecognized bcType (%s) for BC %s!\n", keyValue,
                                                                              bcProps->surfaceProp[i].name);
                status = CAPS_BADVALUE;
                goto cleanup;
            }
            AIM_FREE(keyValue);
        }
    }

    for (i = 0; i < attrMap->numAttribute; i++) {
        for (j = 0; j < bcProps->numSurfaceProp; j++) {

            found = (int) false;
            if (attrMap->attributeIndex[i] == bcProps->surfaceProp[j].bcID) {
                found = (int) true;
                break;
            }
        }

        if (found == (int) false) {
            printf("\tWarning: No boundary condition specified for capsGroup %s!\n", attrMap->attributeName[i]);
        }
    }

    printf("\tDone getting CFD boundary conditions\n");

    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(keyValue);
    return status;
}


// Fill modalAeroelastic in a cfdModalAeroelasticStruct format with modal aeroelastic data from Modal Tuple
int cfd_getModalAeroelastic(int numTuple,
                            capsTuple modalTuple[],
                            cfdModalAeroelasticStruct *modalAeroelastic) {

    /*! \if FUN3D
     * \page cfdModalAeroelastic CFD Modal Aeroelastic
     * Structure for the modal aeroelastic tuple  = ("EigenVector_#", "Value").
     * The tuple name "EigenVector_#" defines the eigen-vector in which the supplied information corresponds to, where "#" should be
     * replaced by the corresponding mode number for the eigen-vector (eg. EigenVector_3 would correspond to the third mode, while
     * EigenVector_6 would be the sixth mode). This notation is the same as found in \ref dataTransferFUN3D. The "Value" must
     * be a JSON String dictionary (see Section \ref jsonStringcfdModalAeroelastic).
     * \endif
     */

    int status; //Function return

    int i; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    int stringLength = 0;

    // Destroy our cfdModalAeroelasticStruct structures coming in if it isn't 0 and NULL already
    status = destroy_cfdModalAeroelasticStruct(modalAeroelastic);
    if (status != CAPS_SUCCESS) return status;

    printf("Getting CFD modal aeroelastic information\n");

    modalAeroelastic->numEigenValue = numTuple;

    if (modalAeroelastic->numEigenValue > 0) {

        modalAeroelastic->eigenValue = (cfdEigenValueStruct *) EG_alloc(modalAeroelastic->numEigenValue * sizeof(cfdEigenValueStruct));
        if (modalAeroelastic->eigenValue == NULL) return EGADS_MALLOC;

    } else {

        printf("Warning: Number of modal aeroelastic tuples is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < modalAeroelastic->numEigenValue; i++) {
        status = initiate_cfdEigenValueStruct(&modalAeroelastic->eigenValue[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    //printf("Number of tuple pairs = %d\n", numTuple);

    for (i = 0; i < modalAeroelastic->numEigenValue; i++) {

        printf("\tEigen-Vector name - %s\n", modalTuple[i].name);

        stringLength = strlen(modalTuple[i].name);

        modalAeroelastic->eigenValue[i].name = (char *) EG_alloc((stringLength + 1)*sizeof(char));
        if (modalAeroelastic->eigenValue[i].name == NULL) return EGADS_MALLOC;

        memcpy(modalAeroelastic->eigenValue[i].name,
               modalTuple[i].name,
               stringLength*sizeof(char));

        modalAeroelastic->eigenValue[i].name[stringLength] = '\0';

        status = sscanf(modalTuple[i].name, "EigenVector_%d", &modalAeroelastic->eigenValue[i].modeNumber);
        if (status != 1) {
            printf("Unable to determine which EigenVector mode number - Defaulting the mode 1!!!\n");
            modalAeroelastic->eigenValue[i].modeNumber = 1;
        }

        // Do we have a json string?
        if (strncmp(modalTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n",bcTuple[i].value);

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *  \section jsonStringcfdModalAeroelastic JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary (eg. "Value" = {"generalMass": 1.0, "frequency": 10.7})
             *  the following keywords ( = default values) may be used:
             *  \endif
             */


            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>frequency = 0.0</B> </li> <br>
             *  This is the frequency of specified mode, in rad/sec.
             *  </ul>
             * \endif
             */
            keyWord = "frequency";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].frequency);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>damping = 0.0</B> </li> <br>
             *  The critical damping ratio of the mode.
             *  </ul>
             * \endif
             */
            keyWord = "damping";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].damping);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>generalMass = 0.0</B> </li> <br>
             *  The generalized mass of the mode.
             *  </ul>
             * \endif
             */
            keyWord = "generalMass";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].generalMass);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>generalDisplacement = 0.0</B> </li> <br>
             *  The generalized displacement used at the starting time step to perturb the mode and excite a dynamic response.
             *  </ul>
             * \endif
             */
            keyWord = "generalDisplacement";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].generalDisplacement);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>generalVelocity = 0.0</B> </li> <br>
             *  The generalized velocity used at the starting time step to perturb the mode and excite a dynamic response.
             *  </ul>
             * \endif
             */
            keyWord = "generalVelocity";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].generalVelocity);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \if (FUN3D)
             *  \page cfdModalAeroelastic
             *
             *  <ul>
             *  <li> <B>generalForce = 0.0</B> </li> <br>
             *  The generalized force used at the starting time step to perturb the mode and excite a dynamic response.
             *  </ul>
             * \endif
             */
            keyWord = "generalForce";
            status = search_jsonDictionary( modalTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &modalAeroelastic->eigenValue[i].generalForce);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

        } else {
            printf("\tA JSON string was NOT provided for tuple %s\n!!!", modalTuple[i].name);
            return CAPS_NOTFOUND;

        }

    }

    printf("Done getting CFD boundary conditions\n");

    return CAPS_SUCCESS;
}

static int _setDesignVariable(void *aimInfo,
                              const char *name,
                 /*@unused@*/ int numAnalysisVal, capsValue *analysisVal,
                              cfdDesignVariableStruct *variable) {

    int status;

    int found = (int) false;
    int i, j, index;

    int numGeomIn;
    capsValue *geomInVal;
    const char *geomInName;

    numGeomIn = aim_getIndex(aimInfo, NULL, GEOMETRYIN);
    if (numGeomIn <= 0) {
        AIM_ERROR(aimInfo, "No DESPMTR available!");
        return CAPS_NULLVALUE;
    }

    // Loop through geometry
    for (i = 0; i < numGeomIn; i++) {

        status = aim_getName(aimInfo, i+1, GEOMETRYIN, &geomInName);
        AIM_STATUS(aimInfo, status);

        if (strcasecmp(name, geomInName) != 0)  continue;

        if(aim_getGeomInType(aimInfo, i+1) == EGADS_OUTSIDE) {
            printf("GeometryIn value %s is a configuration parameter and not a valid design parameter - can't get sensitivity.\n", geomInName);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        status = aim_getValue(aimInfo, i+1, GEOMETRYIN, &geomInVal);
        AIM_STATUS(aimInfo, status);

        status = allocate_cfdDesignVariableStruct(name, geomInVal->length, variable);
        if (status != CAPS_SUCCESS) goto cleanup;

        variable->type = DesignVariableGeometry;

        for (j = 0; j< geomInVal->length; j++) {

            if (geomInVal->type == Double) {

                if (geomInVal->length == 1) {
                    variable->initialValue[j] = geomInVal->vals.real;
                } else {
                    variable->initialValue[j] = geomInVal->vals.reals[j];
                }

                variable->upperBound[j] = geomInVal->limits.dlims[0];
                variable->lowerBound[j] = geomInVal->limits.dlims[1];

            }

            if (geomInVal->type == Integer) {
                if (geomInVal->length == 1) {
                    variable->initialValue[j] = (double) geomInVal->vals.integer;
                } else {
                    variable->initialValue[j] = (double) geomInVal->vals.integers[j];
                }

                variable->upperBound[j] = (double) geomInVal->limits.ilims[0];
                variable->lowerBound[j] = (double) geomInVal->limits.ilims[1];
            }
        }

        found = (int) true;
        break;
    }

    // Analysis
    if (found == (int) false) {

        index = aim_getIndex(aimInfo, name, ANALYSISIN);
        if (index < CAPS_SUCCESS) {
            status = index;
            goto cleanup;
        }

        index -= 1;

        status = allocate_cfdDesignVariableStruct(name,  analysisVal[index].length, variable);
        if (status != CAPS_SUCCESS) goto cleanup;

        variable->type = DesignVariableAnalysis;

        if (analysisVal[index].nullVal != IsNull) { // Can not set upper/lower bounds from capsValue

            for (j = 0; j< analysisVal[index].length; j++) {

                if (analysisVal[index].type == Double) {

                    if (analysisVal[index].length == 1) {
                        variable->initialValue[j] = analysisVal[index].vals.real;
                    } else {
                        variable->initialValue[j] = analysisVal[index].vals.reals[j];
                    }

                    variable->upperBound[j] = analysisVal[index].limits.dlims[0];
                    variable->lowerBound[j] = analysisVal[index].limits.dlims[1];
                }

                if (analysisVal[index].type == Integer) {

                    if (analysisVal[index].length == 1) {
                        variable->initialValue[j] = (double) analysisVal[index].vals.integer;
                    } else {
                        variable->initialValue[j] = (double) analysisVal[index].vals.integers[j];
                    }

                    variable->upperBound[j] = (double) analysisVal[index].limits.ilims[0];
                    variable->lowerBound[j] = (double) analysisVal[index].limits.ilims[1];
                }
            }
        } else {
            printf("Warning: No initial value set for %s\n", variable->name);
        }

        found = (int) true;
    }

    if (found == (int) false) {
        // If we made it this far we haven't found the variable
        printf("Warning: Variable %s is neither a GeometryIn or an AnalysisIn variable.\n", name);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Premature exit in _setDesignVariable status = %d\n", status);

    return status;
}

// Get the design variables from a capsTuple
int cfd_getDesignVariable(void *aimInfo,
                          int numDesignVariableTuple,
                          capsTuple designVariableTuple[],
                          int numAnalysisVal, capsValue *analysisVal,
                          int *numDesignVariable,
                          cfdDesignVariableStruct *variable[]) {

    /*! \if FUN3D
     *
     * \page cfdDesignVariable CFD Design Variable
     * Structure for the design variable tuple  = ("DesignVariable Name", "Value").
     * "DesignVariable Name" defines the reference name for the design variable being specified.
     * The "Value" may be a JSON String dictionary (see Section \ref jsonStringDesignVariable) or just
     * a blank string (see Section \ref stringDesignVariable).
     *
     * \endif
     */
    int status; //Function return

    int i, j; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    int length = 0;
    double *tempArray = NULL;

    cfdDesignVariableStruct *var = NULL;

    // Destroy our design variables structures coming in if aren't 0 and NULL already
    if (*variable != NULL) {
        for (i = 0; i < *numDesignVariable; i++) {
            status = destroy_cfdDesignVariableStruct(&(*variable)[i]);
            AIM_STATUS(aimInfo, status);
        }
    }

    AIM_FREE(*variable);
    *variable = NULL;
    *numDesignVariable = 0;

    printf("\nGetting CFD design variables.......\n");

    *numDesignVariable = numDesignVariableTuple;

    printf("\tNumber of design variables - %d\n", numDesignVariableTuple);

    if (numDesignVariableTuple > 0) {
        AIM_ALLOC(*variable, numDesignVariableTuple, cfdDesignVariableStruct, aimInfo, status);
    } else {
        printf("\tNumber of design variable values in input tuple is 0\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }


    for (i = 0; i < numDesignVariableTuple; i++) {
        status = initiate_cfdDesignVariableStruct(&(*variable)[i]);
        if (status != CAPS_SUCCESS) {
            status = EGADS_MALLOC;
            goto cleanup;
        }
    }

    for (i = 0; i < numDesignVariableTuple; i++) {

        var = &(*variable)[i];

        printf("\tDesign Variable name - %s\n", designVariableTuple[i].name);

        status = _setDesignVariable(aimInfo,
                                    designVariableTuple[i].name,
                                    numAnalysisVal, analysisVal,
                                    var);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Do we have a json string?
        if (strncmp(designVariableTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", designVariableTuple[i].value);

            /*! \if FUN3D
             *  \page cfdDesignVariable
             * \section jsonStringDesignVariable JSON String Dictionary
             *
             * If "Value" is JSON string dictionary (eg. "Value" = {"upperBound": 10.0})
             * the following keywords ( = default values) may be used:
             *
             * \endif
             */

            //Fill up designVariable properties

            /*! \if FUN3D
             *
             *  \page cfdDesignVariable
             *
             * <ul>
             *  <li> <B>lowerBound = 0.0 or [0.0, 0.0,...]</B> </li> <br>
             *  Lower bound for the design variable.
             * </ul>
             * \endif
             */
            keyWord = "lowerBound";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                //status = string_toDouble(keyValue, &(*variable)[i].lowerBound);
                status =  string_toDoubleDynamicArray(keyValue, &length, &tempArray);
                string_free(&keyValue);
                if (status != CAPS_SUCCESS) goto cleanup;

                if (var->length != length) {
                    printf("Inconsistent lower bound lengths for %s!\n", var->name);
                    status = CAPS_MISMATCH;
                    goto cleanup;
                } else {
                    for (j = 0; j < length; j++) var->lowerBound[j] = tempArray[j];
                }

                if (tempArray != NULL) EG_free(tempArray);
                tempArray = NULL;
            }

            /*! \if FUN3D
             *
             * \page cfdDesignVariable
             *
             * <ul>
             *  <li> <B>upperBound = 0.0</B> </li> <br>
             *  Upper bound for the design variable.
             * </ul>
             * \endif
             */
            keyWord = "upperBound";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                // status = string_toDouble(keyValue, &(*variable)[i].upperBound);
                status =  string_toDoubleDynamicArray(keyValue, &length, &tempArray);
                     string_free(&keyValue);
                     if (status != CAPS_SUCCESS) goto cleanup;

                     if (var->length != length) {
                         printf("Inconsistent upper bound lengths for %s!\n", var->name);
                         status = CAPS_MISMATCH;
                         goto cleanup;
                     } else {
                         for (j = 0; j < length; j++) var->upperBound[j] = tempArray[j];
                     }

                     if (tempArray != NULL) EG_free(tempArray);
                     tempArray = NULL;
            }

        } else {

            /*! \if FUN3D
             *
             * \page cfdDesignVariable
             * \section keyStringDesignVariable Single Value String
             *
             * If "Value" is a string, the string value will be ignored.
             *
             * \endif
             */

//            printf("\tError: Design_Variable tuple value is expected to be a JSON string\n");
//            status = CAPS_BADVALUE;
//            goto cleanup;
        }
    }

    printf("\tDone getting CFD design variables\n");

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:

        if (status != CAPS_SUCCESS) printf("Error: Premature exit in cfd_getDesignVariable status = %d\n", status);

        if (tempArray != NULL) EG_free(tempArray);
        string_free(&keyValue);

        return status;
}

// Fill objective in a cfdDesignObjectiveStruct format with objective data from Objective Tuple
int cfd_getDesignObjective(int numObjectiveTuple,
                           capsTuple objectiveTuple[],
                           int *numObjective,
                           cfdDesignObjectiveStruct *objective[]) {

    /*! \if FUN3D
     *
     * \page cfdDesignObjective CFD Objective
     * Structure for the design objective tuple  = ("Objective Name", "Value").
     * "Objective Name" defines the objective function being specified.
     *  The "Value" must be a JSON String dictionary (see Section \ref jsonStringcfdDesignObjective).
     *
     * For FUN3D, the objective function in which the adjoint will be taken with respect to can be, "Cl",
     * "Cd". "Cmx", "Cmy", "Cmz", "ClCd", "Cx", "Cy", "Cz" (eg. "cl" would correspond to the lift coefficient,
     * while "clcd" would be the aerodynamic efficiency).
     *
     * FUN3D calculates the objective function using the following form:
     * f = sum( w_i *(C_i - C_i^*) ^ p_i)
     * Where:
     * f : Objective function
     * w_i : Weighting of objective function
     * C_i : Objective function type (cl, cd, etc.)
     * C_i^* : Objective function target
     * p_i : Exponential factor of objective function
     *
     * \endif
     */

    int status; //Function return

    int i; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    // Destroy our design objective structures coming in if aren't 0 and NULL already
     if (*objective != NULL) {
         for (i = 0; i < *numObjective; i++) {
             status = destroy_cfdDesignObjectiveStruct(&(*objective)[i]);
             if (status != CAPS_SUCCESS) goto cleanup;
         }
     }

     if (*objective != NULL) EG_free(*objective);
     *objective = NULL;
     *numObjective = 0;

    printf("Getting CFD objective.......\n");

    *numObjective = numObjectiveTuple;

    printf("\tNumber of design variables - %d\n", numObjectiveTuple);

    if (numObjectiveTuple > 0) {

        *objective = (cfdDesignObjectiveStruct *) EG_alloc(numObjectiveTuple * sizeof(cfdDesignObjectiveStruct));
        if (objective == NULL) {
            *numObjective = 0;
            status = EGADS_MALLOC;
            goto cleanup;
        }

    } else {
        printf("\tNumber of objective values in input tuple is 0\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    for (i = 0; i < numObjectiveTuple; i++) {
        status = initiate_cfdDesignObjectiveStruct(&(*objective)[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    for (i = 0; i < numObjectiveTuple; i++) {

        printf("\tObjective name - %s\n", objectiveTuple[i].name);

        (*objective)[i].name = (char *) EG_alloc(((strlen(objectiveTuple[i].name)) + 1)*sizeof(char));
        if ((*objective)[i].name == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        memcpy((*objective)[i].name, objectiveTuple[i].name, strlen(objectiveTuple[i].name)*sizeof(char));
        (*objective)[i].name[strlen(objectiveTuple[i].name)] = '\0';

        // Get adjoint function type
        if      ( strcasecmp(objectiveTuple[i].name, "cl"  ) == 0) (*objective)[i].objectiveType = ObjectiveCl;
        else if ( strcasecmp(objectiveTuple[i].name, "cd"  ) == 0) (*objective)[i].objectiveType = ObjectiveCd;
        else if ( strcasecmp(objectiveTuple[i].name, "cmx" ) == 0) (*objective)[i].objectiveType = ObjectiveCmx;
        else if ( strcasecmp(objectiveTuple[i].name, "cmy" ) == 0) (*objective)[i].objectiveType = ObjectiveCmy;
        else if ( strcasecmp(objectiveTuple[i].name, "cmz" ) == 0) (*objective)[i].objectiveType = ObjectiveCmz;
        else if ( strcasecmp(objectiveTuple[i].name, "clcd") == 0) (*objective)[i].objectiveType = ObjectiveClCd;
        else if ( strcasecmp(objectiveTuple[i].name, "cx"  ) == 0) (*objective)[i].objectiveType = ObjectiveCx;
        else if ( strcasecmp(objectiveTuple[i].name, "cy"  ) == 0) (*objective)[i].objectiveType = ObjectiveCy;
        else if ( strcasecmp(objectiveTuple[i].name, "cz"  ) == 0) (*objective)[i].objectiveType = ObjectiveCz;
        else {
            printf("\tFunction name not recognized: '%s'\n", objectiveTuple[i].name);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // Do we have a json string?
        if (strncmp(objectiveTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n",bcTuple[i].value);

            /*! \if FUN3D
             * \page cfdDesignObjective
             *
             *  \section jsonStringcfdDesignObjective JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary (eg. "Value" = {"weight": 1.0, "target": 10.7})
             *  the following keywords ( = default values) may be used:
             *
             * \endif
             */


            /*! \if FUN3D
             *
             *  \page cfdDesignObjective
             *
             *  <ul>
             *  <li> <B>weight = 1.0</B> </li> <br>
             *  This weighting of the objective function.
             *  </ul>
             * \endif
             */
            keyWord = "weight";
            status = search_jsonDictionary( objectiveTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*objective)[i].weight);
                string_free(&keyValue);

                if (status != CAPS_SUCCESS) goto cleanup;
            }

            /*! \if FUN3D
             *
             *  \page cfdDesignObjective
             *
             *  <ul>
             *  <li> <B>target = 0.0</B> </li> <br>
             *  This is the target value of the objective function.
             *  </ul>
             * \endif
             */
            keyWord = "target";
            status = search_jsonDictionary( objectiveTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*objective)[i].target);
                string_free(&keyValue);

                if (status != CAPS_SUCCESS) goto cleanup;
            }

            /*! \if FUN3D
             *
             * \page cfdDesignObjective
             *
             *  <ul>
             *  <li> <B>power = 1.0</B> </li> <br>
             *  This is the user defined power operator for the objective function.
             *  </ul>
             * \endif
             */
            keyWord = "power";
            status = search_jsonDictionary( objectiveTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*objective)[i].power);
                string_free(&keyValue);

                if (status != CAPS_SUCCESS) goto cleanup;
            }

        } else {
            printf("\tA JSON string was NOT provided for tuple %s\n!!!", objectiveTuple[i].name);
            printf("\tTuple value is expected to be a JSON string\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting CFD objective\n");

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:

        if (status != CAPS_SUCCESS) printf("Error: Premature exit in cfd_getDesignObjective status = %d\n", status);

        string_free(&keyValue);

        return status;
}


// Finish filling out the design variable information


// Initiate (0 out all values and NULL all pointers) of surfaceProps in the cfdSurfaceStruct structure format
int initiate_cfdSurfaceStruct(cfdSurfaceStruct *surfaceProps) {

    if (surfaceProps == NULL) return CAPS_NULLVALUE;

    surfaceProps->surfaceType = UnknownBoundary;

    surfaceProps->name = NULL;

    surfaceProps->bcID = 0; // ID of boundary

    // Wall specific properties
    surfaceProps->wallTemperatureFlag = (int) false; // Wall temperature flag
    surfaceProps->wallTemperature = -1.0;            // Wall temperature value -1 = adiabatic ; >0 = isothermal
    surfaceProps->wallHeatFlux = 0.0;                // Wall heat flux. to use Temperature flag should be true and wallTemperature < 0

    // Symmetry plane
    surfaceProps->symmetryPlane = 0; // Symmetry flag / plane

    // Stagnation quantities
    surfaceProps->totalPressure = 0;    // Total pressure
    surfaceProps->totalTemperature = 0; // Total temperature
    surfaceProps->totalDensity = 0;     // Total density

    // Static quantities
    surfaceProps->staticPressure = 0;   // Static pressure
    surfaceProps->staticTemperature = 0;// Static temperature
    surfaceProps->staticDensity = 0;    // Static temperature

    // Velocity components
    surfaceProps->uVelocity = 0;  // x-component of velocity
    surfaceProps->vVelocity = 0;  // y-component of velocity
    surfaceProps->wVelocity = 0;  // z-compoentn of velocity
    surfaceProps->machNumber = 0; // Mach number

    // Massflow
    surfaceProps->massflow = 0; // Mass flow through a boundary

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of surfaceProps in the cfdSurfaceStruct structure format
int destroy_cfdSurfaceStruct(cfdSurfaceStruct *surfaceProps) {

    if (surfaceProps == NULL) return CAPS_NULLVALUE;

    surfaceProps->surfaceType = UnknownBoundary;

    if (surfaceProps->name != NULL) EG_free(surfaceProps->name);
    surfaceProps->name = NULL;

    surfaceProps->bcID = 0;             // ID of boundary

    // Wall specific properties
    surfaceProps->wallTemperatureFlag = (int) false; // Wall temperature flag
    surfaceProps->wallTemperature = -1;              // Wall temperature value -1 = adiabatic ; >0 = isothermal

    // Symmetry plane
    surfaceProps->symmetryPlane = 0; // Symmetry flag / plane

    // Stagnation quantities
    surfaceProps->totalPressure = 0;    // Total pressure
    surfaceProps->totalTemperature = 0; // Total temperature
    surfaceProps->totalDensity = 0;     // Total density

    // Static quantities
    surfaceProps->staticPressure = 0;   // Static pressure
    surfaceProps->staticTemperature = 0;// Static temperature
    surfaceProps->staticDensity = 0;    // Static temperature

    // Velocity components
    surfaceProps->uVelocity = 0;  // x-component of velocity
    surfaceProps->vVelocity = 0;  // y-component of velocity
    surfaceProps->wVelocity = 0;  // z-compoentn of velocity
    surfaceProps->machNumber = 0; // Mach number

    // Massflow
    surfaceProps->massflow = 0; // Mass flow through a boundary

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of bcProps in the cfdBoundaryConditionStruct structure format
int initiate_cfdBoundaryConditionStruct(cfdBoundaryConditionStruct *bcProps) {

    if (bcProps == NULL) return CAPS_NULLVALUE;

    bcProps->name = NULL;

    bcProps->numSurfaceProp = 0;

    bcProps->surfaceProp = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of bcProps in the cfdBoundaryConditionStruct structure format
int destroy_cfdBoundaryConditionStruct(cfdBoundaryConditionStruct *bcProps) {

    int status; // Function return
    int i; // Indexing

    if (bcProps == NULL) return CAPS_NULLVALUE;

    AIM_FREE(bcProps->name);
 
    for  (i = 0; i < bcProps->numSurfaceProp; i++) {
        status = destroy_cfdSurfaceStruct(&bcProps->surfaceProp[i]);
        if (status != CAPS_SUCCESS) printf("Error in destroy_cfdBoundaryConditionStruct, status = %d\n", status);
    }

    bcProps->numSurfaceProp = 0;
    AIM_FREE(bcProps->surfaceProp);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of eigenValue in the cfdEigenValueStruct structure format
int initiate_cfdEigenValueStruct(cfdEigenValueStruct *eigenValue) {

    if (eigenValue == NULL) return CAPS_NULLVALUE;

    eigenValue->name = NULL;

    eigenValue->modeNumber = 0;

    eigenValue->frequency = 0.0;
    eigenValue->damping = 0.0;

    eigenValue->generalMass = 0.0;
    eigenValue->generalDisplacement = 0.0;
    eigenValue->generalVelocity = 0.0;
    eigenValue->generalForce = 0.0;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of eigenValue in the cfdEigenValueStruct structure format
int destroy_cfdEigenValueStruct(cfdEigenValueStruct *eigenValue) {

    if (eigenValue == NULL) return CAPS_NULLVALUE;

    if (eigenValue->name != NULL) EG_free(eigenValue->name);
    eigenValue->name = NULL;

    eigenValue->modeNumber = 0;

    eigenValue->frequency = 0.0;
    eigenValue->damping = 0.0;

    eigenValue->generalMass = 0.0;
    eigenValue->generalDisplacement = 0.0;
    eigenValue->generalVelocity = 0.0;
    eigenValue->generalForce = 0.0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of modalAeroelastic in the cfdModalAeroelasticStruct structure format
int initiate_cfdModalAeroelasticStruct(cfdModalAeroelasticStruct *modalAeroelastic) {

    if (modalAeroelastic == NULL) return CAPS_NULLVALUE;

    modalAeroelastic->surfaceID = 0;

    modalAeroelastic->numEigenValue = 0;
    modalAeroelastic->eigenValue = NULL;

    modalAeroelastic->freestreamVelocity = 0.0;
    modalAeroelastic->freestreamDynamicPressure = 0.0;
    modalAeroelastic->lengthScaling = 1.0;

    return CAPS_SUCCESS;
}


// Destroy (0 out all values and NULL all pointers) of modalAeroelastic in the cfdModalAeroelasticStruct structure format
int destroy_cfdModalAeroelasticStruct(cfdModalAeroelasticStruct *modalAeroelastic) {

    int status; // Function return status;
    int i; // Indexing

    if (modalAeroelastic == NULL) return CAPS_NULLVALUE;

    modalAeroelastic->surfaceID = 0;

    if ( modalAeroelastic->eigenValue != NULL) {

        for (i = 0; i < modalAeroelastic->numEigenValue; i++) {

            status = destroy_cfdEigenValueStruct(&modalAeroelastic->eigenValue[i]);

            if (status != CAPS_SUCCESS) printf("Error in destroy_cfdModalAeroelasticStruct, status = %d\n", status);
        }

        EG_free(modalAeroelastic->eigenValue);
    }

    modalAeroelastic->numEigenValue = 0;
    modalAeroelastic->eigenValue = NULL;

    modalAeroelastic->freestreamVelocity = 0.0;
    modalAeroelastic->freestreamDynamicPressure = 0.0;
    modalAeroelastic->lengthScaling = 1.0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of designVariable in the cfdDesignVariableStruct structure format
int initiate_cfdDesignVariableStruct(cfdDesignVariableStruct *designVariable) {

    if (designVariable == NULL ) return CAPS_NULLVALUE;

    designVariable->name = NULL;

    designVariable->type = DesignVariableUnknown;

    designVariable->length = 0;
    designVariable->initialValue = NULL;
    designVariable->lowerBound = NULL;
    designVariable->upperBound = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of designVariable in the cfdDesignVariableStruct structure format
int destroy_cfdDesignVariableStruct(cfdDesignVariableStruct *designVariable) {

    if (designVariable == NULL ) return CAPS_NULLVALUE;

    EG_free(designVariable->name);
    designVariable->name = NULL;

    designVariable->type = DesignVariableUnknown;

    designVariable->length = 0;

    EG_free(designVariable->initialValue);
    designVariable->initialValue = NULL;

    EG_free(designVariable->lowerBound);
    designVariable->lowerBound = NULL;

    EG_free(designVariable->upperBound);
    designVariable->upperBound = NULL;

    return CAPS_SUCCESS;
}

// Allocate cfdDesignVariableStruct structure
int allocate_cfdDesignVariableStruct(const char *name,  int length, cfdDesignVariableStruct *designVariable) {

    int status;

    int i;

    designVariable->name = EG_strdup(name);
    if (designVariable->name == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    designVariable->length = length;

    designVariable->initialValue = (double *) EG_alloc(length*sizeof(double));
    if (designVariable->initialValue == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    designVariable->lowerBound = (double *) EG_alloc(length*sizeof(double));
    if (designVariable->lowerBound == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    designVariable->upperBound = (double *) EG_alloc(length*sizeof(double));
    if (designVariable->upperBound == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < length; i++) {
        designVariable->initialValue[i] = 0.0;
        designVariable->lowerBound[i] = 0.0;
        designVariable->upperBound[i] = 0.0;
    }

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("Error: Premature exit in allocate_cfdDesignVariableStruct status = %d\n", status);

        return status;
}

// Initiate (0 out all values and NULL all pointers) of objective in the cfdDesignObjectiveStruct structure format
int initiate_cfdDesignObjectiveStruct(cfdDesignObjectiveStruct *objective) {

    if (objective == NULL) return CAPS_NULLVALUE;

    objective->name = NULL;

    objective->objectiveType = ObjectiveUnknown;

    objective->target = 0.0;
    objective->weight = 1.0;
    objective->power = 1.0;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of objective in the cfdDesignObjectiveStruct structure format
int destroy_cfdDesignObjectiveStruct(cfdDesignObjectiveStruct *objective) {

    if (objective == NULL) return CAPS_NULLVALUE;

    EG_free(objective->name);
    objective->name = NULL;

    objective->objectiveType = ObjectiveUnknown;

    objective->target = 0.0;
    objective->weight = 1.0;
    objective->power = 1.0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of objective in the cfdDesignStruct structure format
int initiate_cfdDesignStruct(cfdDesignStruct *design) {

    if (design == NULL ) return CAPS_NULLVALUE;

    design->numDesignObjective = 0;
    design->designObjective = NULL; // [numObjective]

    design->numDesignVariable = 0;
    design->designVariable = NULL; // [numDesignVariable]

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of objective in the cfdDesignStruct structure format
int destroy_cfdDesignStruct(cfdDesignStruct *design) {

    int i;

    if (design == NULL ) return CAPS_NULLVALUE;

    for (i = 0; i < design->numDesignObjective; i++)
      destroy_cfdDesignObjectiveStruct(design->designObjective + i);

    design->numDesignObjective = 0;
    EG_free(design->designObjective); // [numObjective]
    design->designObjective = NULL;

    for (i = 0; i < design->numDesignVariable; i++)
      destroy_cfdDesignVariableStruct(design->designVariable + i);

    design->numDesignVariable = 0;
    EG_free(design->designVariable); // [numDesignVariable]
    design->designVariable = NULL;

    return CAPS_SUCCESS;

}

// Initiate (0 out all values and NULL all pointers) of objective in the cfdUnitsStruct structure format
int initiate_cfdUnitsStruct(cfdUnitsStruct *units)
{
    if (units == NULL ) return CAPS_NULLVALUE;

    // base units
    units->length       = NULL;
    units->mass         = NULL;
    units->time         = NULL;
    units->temperature  = NULL;

    // derived units
    units->density      = NULL;
    units->pressure     = NULL;
    units->speed        = NULL;
    units->acceleration = NULL;
    units->force        = NULL;
    units->viscosity    = NULL;
    units->area         = NULL;

    // coefficient units
    units->Cpressure    = NULL;
    units->Cforce       = NULL;
    units->Cmoment      = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of objective in the cfdUnitsStruct structure format
int destroy_cfdUnitsStruct(cfdUnitsStruct *units)
{
    if (units == NULL ) return CAPS_NULLVALUE;

    // base units
    AIM_FREE(units->length     );
    AIM_FREE(units->mass       );
    AIM_FREE(units->time       );
    AIM_FREE(units->temperature);

    // derived units
    AIM_FREE(units->density     );
    AIM_FREE(units->pressure    );
    AIM_FREE(units->speed       );
    AIM_FREE(units->acceleration);
    AIM_FREE(units->force       );
    AIM_FREE(units->viscosity   );
    AIM_FREE(units->area        );

    // coefficient units
    AIM_FREE(units->Cpressure  );
    AIM_FREE(units->Cforce     );
    AIM_FREE(units->Cmoment    );

    return CAPS_SUCCESS;
}

// Compute derived units from bas units
int cfd_cfdDerivedUnits(void *aimInfo, cfdUnitsStruct *units)
{

    int status = CAPS_SUCCESS;
    char *volume=NULL, *tmpUnit=NULL;

    if (units == NULL) return CAPS_NULLVALUE;
    if (units->length      == NULL ) return CAPS_NULLVALUE;
    if (units->mass        == NULL ) return CAPS_NULLVALUE;
    if (units->time        == NULL ) return CAPS_NULLVALUE;
    if (units->temperature == NULL ) return CAPS_NULLVALUE;

    AIM_FREE(units->density     );
    AIM_FREE(units->pressure    );
    AIM_FREE(units->speed       );
    AIM_FREE(units->acceleration);
    AIM_FREE(units->force       );
    AIM_FREE(units->viscosity   );
    AIM_FREE(units->area        );

    // construct area unit
    status = aim_unitRaise(aimInfo, units->length, 2, &units->area); // length^2
    AIM_STATUS(aimInfo, status);

    // construct volume unit
    status = aim_unitRaise(aimInfo, units->length, 3, &volume); // length^3
    AIM_STATUS(aimInfo, status);
  
    // construct density unit
    status = aim_unitDivide(aimInfo, units->mass, volume, &units->density); // mass/length^3
    AIM_STATUS(aimInfo, status);

    // construct speed unit
    status = aim_unitDivide(aimInfo, units->length, units->time, &units->speed); // length/time
    AIM_STATUS(aimInfo, status);

    // construct acceleration unit
    status = aim_unitDivide(aimInfo, units->speed, units->time, &units->acceleration); // length/time^2
    AIM_STATUS(aimInfo, status);

    // construct force unit
    status = aim_unitMultiply(aimInfo, units->mass, units->acceleration, &units->force); // mass * length / time^2 (force)
    AIM_STATUS(aimInfo, status);

    // construct pressure unit
    status = aim_unitDivide(aimInfo, units->force, units->area, &units->pressure); // force/length^2
    AIM_STATUS(aimInfo, status);
 
    // construct viscosity unit
    status = aim_unitMultiply(aimInfo, units->time, units->pressure, &units->viscosity); // time * force/length^2
    AIM_STATUS(aimInfo, status);

cleanup:
    AIM_FREE(volume);
    AIM_FREE(tmpUnit);

    return status;
}

// Compute coefficient units from reference quantities
int cfd_cfdCoefficientUnits(void *aimInfo,
                            double length  , const char *lengthUnit,
                            double area    , const char *areaUnit,
                            double density , const char *densityUnit,
                            double speed   , const char *speedUnit,
                            double pressure, const char *pressureUnit,
                            cfdUnitsStruct *units)
{
    int status = CAPS_SUCCESS;

    char buffer[32], *refLengthUnit=NULL, *refAreaUnit=NULL, *refDynPressureUnit=NULL;
    char *refPressureUnit=NULL, *refOffsetPressure=NULL, *speedUnit2=NULL, *dyn_pressureUnit=NULL;
    char *tmp=NULL;
    double dynamic_pressure = 0.0, offset_pressure = 0.0;

    if (units == NULL) return CAPS_NULLVALUE;

    AIM_FREE(units->Cpressure);
    AIM_FREE(units->Cforce);
    AIM_FREE(units->Cmoment);

    // construct reference length/area units
    snprintf(buffer, 32, "%.16e", length);
    status = aim_unitMultiply(aimInfo, buffer, lengthUnit, &refLengthUnit);
    AIM_STATUS(aimInfo, status);

    snprintf(buffer, 32, "%.16e", area);
    status = aim_unitMultiply(aimInfo, buffer, areaUnit, &refAreaUnit);
    AIM_STATUS(aimInfo, status);


    // construct reference dynamic pressure unit
    status = aim_unitRaise(aimInfo, speedUnit, 2, &speedUnit2); // speed^2
    AIM_STATUS(aimInfo, status);

    status = aim_unitMultiply(aimInfo, densityUnit, speedUnit2, &dyn_pressureUnit);
    AIM_STATUS(aimInfo, status);

    dynamic_pressure = 0.5 * density * pow(speed, 2.0);
    snprintf(buffer, 32, "%.16e", dynamic_pressure);
    status = aim_unitMultiply(aimInfo, buffer, dyn_pressureUnit, &refDynPressureUnit);
    AIM_STATUS(aimInfo, status);

    // construct reference pressure units
    snprintf(buffer, 32, "%.16e", pressure);
    status = aim_unitMultiply(aimInfo, buffer, pressureUnit, &refPressureUnit);
    AIM_STATUS(aimInfo, status);

    // construct pressure coefficient offset
    status = aim_unitDivide(aimInfo, refPressureUnit, refDynPressureUnit, &refOffsetPressure); // Reference_Pressure/dynamic_pressure
    AIM_STATUS(aimInfo, status);
    offset_pressure = strtod(refOffsetPressure, &tmp);

    // construct pressure coefficient unit
    status = aim_unitOffset(aimInfo, refDynPressureUnit, offset_pressure, &units->Cpressure); // Pressure coefficient
    AIM_STATUS(aimInfo, status);

    // construct force coefficient unit
    status = aim_unitMultiply(aimInfo, refDynPressureUnit, refAreaUnit, &units->Cforce); // Force coefficient
    AIM_STATUS(aimInfo, status);

    // construct moment coefficient unit
    status = aim_unitMultiply(aimInfo, units->Cforce, refLengthUnit, &units->Cmoment); // Moment coefficient
    AIM_STATUS(aimInfo, status);

cleanup:

    AIM_FREE(refLengthUnit);
    AIM_FREE(refAreaUnit);
    AIM_FREE(refDynPressureUnit);
    AIM_FREE(refPressureUnit);
    AIM_FREE(refOffsetPressure);
    AIM_FREE(speedUnit2);
    AIM_FREE(dyn_pressureUnit);

    return status;
}
