#include <string.h>

#include "aimUtil.h"
#include "miscUtils.h"
#include "jsonUtils.h"
#include "zaeroUtils.h"
#include "zaeroCards.h"
#include "zaeroDiscipline.h"


/* Get Trim variables from input tuples */
int zaero_getTrimVariables(void *aimInfo,
                           int numtrimVarTuple,
                           capsTuple trimVarTuple[],
                           zaeroProblemStruct *zaeroProblem)
{
  /*! \page zaeroTRIMVAR ZAero Trim_Variable
   * Structure for the Trim_Variable tuple  = ("Trim Variable Name", "Value").
   * "Trim Variable Name" defines the reference name for the Trim being specified.
   * The "Value" must be a JSON String dictionary.
   */

  int i, status;

  zaeroTrimVarStruct *trimvar;

  // Ensure we are starting with no trim variables
  if (zaeroProblem->trimVariables != NULL) {

    for (i = 0; i < zaeroProblem->numTrimVariables; i++) {

      status = destroy_zaeroTrimVarStruct(&zaeroProblem->trimVariables[i]);
      if (status != CAPS_SUCCESS) {
        return status;
      }

    }
    AIM_FREE(zaeroProblem->trimVariables);
  }
  zaeroProblem->trimVariables = NULL;

  printf("\nGetting ZAERO Trim Variables.......\n");

  zaeroProblem->numTrimVariables = numtrimVarTuple;
  printf("\tNumber of Trim Variables - %d\n", zaeroProblem->numTrimVariables);

  if (zaeroProblem->numTrimVariables > 0) {

    AIM_ALLOC(zaeroProblem->trimVariables, zaeroProblem->numTrimVariables, zaeroTrimVarStruct, aimInfo, status);

    // initiate TrimVar structure
    for (i = 0; i < zaeroProblem->numTrimVariables; i++) {
      status = initiate_zaeroTrimVarStruct(&zaeroProblem->trimVariables[i]);
      AIM_STATUS(aimInfo, status);
    }
  } else {
    AIM_ERROR(aimInfo, "Number of analysis trimVariables in 'Trim_Variable' "
              "tuple is %d", zaeroProblem->numTrimVariables);
    return CAPS_NOTFOUND;
  }

  // for each analysis TrimVar tuple
  for (i = 0; i < zaeroProblem->numTrimVariables; i++) {

    trimvar = &zaeroProblem->trimVariables[i];

    // set name
    AIM_STRDUP(trimvar->name, trimVarTuple[i].name, aimInfo, status);

    // make sure TrimVar tuple value is json string
    if (!json_isDict(trimVarTuple[i].value)) {
      AIM_ERROR(aimInfo, "'Trim_Variable' tuple value must be a JSON dictionary");
      return CAPS_BADVALUE;
    }

    /*! \page zaeroTRIMVAR
     * \section jsonTRIMVAR Trim_Variable JSON String Dictionary
     *
     * For the JSON string "Value" dictionary
     *  (e.g. "Value" = {"label": "ALPHA", "value": "free"})
     * \endif
     *  the following keywords ( = default values) may be used:
     *
     * <ul>
     * <li> <B>label</B> </li> <br>
     *    The trim label
     * </ul>
     */

    status = json_getString(trimVarTuple[i].value, "label", &trimvar->label);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"label\" "
                         "in 'Trim_Variable' tuple value");
      return status;
    }

    /*! \page zaeroTRIMVAR
     * <ul>
     * <li> <B>value</B> </li> <br>
     *   Trime value.
     * </ul>
     */

    status = json_getString(trimVarTuple[i].value, "value", &trimvar->value);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"value\" "
                         "in 'Trim_Variable' tuple value");
      return status;
    }

    // #ifdef DEBUG
    printf("\n\tTrim Variable: %s\n", trimvar->name);
    printf("\t- Label : %s\n", trimvar->label);
    printf("\t- Value : %s\n", trimvar->value);
    // #endif
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

int zaero_getTrimDiscipline(void *aimInfo,
                            const char *jsonDict,
                            const cfdUnitsStruct *units,
                            zaeroSubcaseStruct *subcase)
{
  int status; // Function return

  char *keyValue = NULL;
  char *keyWord = NULL;

  double I[6];

  zaeroTrimStruct *trim;

  if (subcase->discipline != NULL) {
    AIM_ERROR(aimInfo, "subcase discipline is already set");
    return CAPS_BADVALUE;
  }

  AIM_ALLOC(subcase->discipline, 1, zaeroTrimStruct, aimInfo, status);

  trim = (zaeroTrimStruct *) subcase->discipline;
  initiate_zaeroTrimStruct(trim);

  /*! \page zaeroAnalysis
   * \section jsonTRIM StaticAeroelastic or Trim Discipline JSON String Dictionary
   *
   * For the JSON string "Value" dictionary
   *  (e.g. "Value" = {"dynamicPressure": rho*V^2/2, "gravityAcceleration": 9.81 m/s^2})
   * \endif
   *  the following keywords ( = default values) may be used:
   *
   * <ul>
   * <li> <B>dynamicPressure</B> </li> <br>
   *    Dynamic pressure
   * </ul>
   */

  keyWord = "dynamicPressure";
  status = search_jsonDictionary( jsonDict, keyWord, &keyValue);
  if (status == CAPS_SUCCESS) {
    if ( units->pressure != NULL ) {
        status = string_toDoubleUnits(aimInfo, keyValue, units->pressure, &trim->dynamicPressure);
    } else {
        status = string_toDouble(keyValue, &trim->dynamicPressure);
    }
    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
    AIM_FREE(keyValue);
  } else {
    AIM_ERROR(aimInfo, "missing \"%s\" entry required "
                       "by Trim discipline", keyWord);
    goto cleanup;
  }

  /*! \page zaeroAnalysis
   * <ul>
   * <li> <B>vectorToCG = [0,0,0]</B> </li> <br>
   *    Vector to offset the CG
   * </ul>
   */

  keyWord = "vectorToCG";
  status = search_jsonDictionary( jsonDict, keyWord, &keyValue);
  if (status == CAPS_SUCCESS) {
    if ( units->length != NULL ) {
        status = string_toDoubleArrayUnits(aimInfo, keyValue, units->length, 3, trim->vectorToCG);
    } else {
        status = string_toDoubleArray(keyValue, 3, trim->vectorToCG);
    }
    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
    AIM_FREE(keyValue);
  }

  /*! \page zaeroAnalysis
   *  <ul>
   * <li> <B>gravityAcceleration = NULL</B> </li> <br>
   *    Gravity acceleration
   * </ul>
   */

  keyWord = "gravityAcceleration";
  status = search_jsonDictionary( jsonDict, keyWord, &keyValue);
  if (status == CAPS_SUCCESS) {
    if ( units->acceleration != NULL ) {
        status = string_toDoubleUnits(aimInfo, keyValue, units->acceleration, &trim->gravityAcceleration);
    } else {
        status = string_toDouble(keyValue, &trim->gravityAcceleration);
    }
    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
    AIM_FREE(keyValue);
  } else {
    AIM_ERROR(aimInfo, "missing \"%s\" entry required "
                       "by Trim discipline", keyWord);
    goto cleanup;
  }

  /*! \page zaeroAnalysis
   * <ul>
   * <li> <B>weight = NULL</B> </li> <br>
   *    Weight of the overall aircraft if specified
   * </ul>
   */

  keyWord = "weight";
  status = search_jsonDictionary( jsonDict, keyWord, &keyValue);
  if (status == CAPS_SUCCESS) {
    if ( units->force != NULL ) {
        status = string_toDoubleUnits(aimInfo, keyValue, units->force, &trim->weight);
    } else {
        status = string_toDouble(keyValue, &trim->weight);
    }
    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
    AIM_FREE(keyValue);
  }

  /*! \page zaeroAnalysis
   * <ul>
   * <li> <B>weightMomentOfInertia = NULL</B> </li> <br>
   *    Moments of inertia of the full aircraft if specified: [IXX, IXY, IYY, IXZ, IYZ, IZZ]
   * </ul>
   */


  // IXX, IXY, IYY, IXZ, IYZ, IZZ
  keyWord = "weightMomentOfInertia";
  status = search_jsonDictionary( jsonDict, keyWord, &keyValue);
  if (status == CAPS_SUCCESS) {
    if ( units->momentOfInertia != NULL ) {
        status = string_toDoubleArrayUnits(aimInfo, keyValue, units->momentOfInertia, 6, I);
    } else {
        status = string_toDoubleArray(keyValue, 6, I);
    }
    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
    AIM_FREE(keyValue);

    trim->weightMomentOfInertia[I11] = I[0];
    trim->weightMomentOfInertia[I21] = I[1];
    trim->weightMomentOfInertia[I22] = I[2];
    trim->weightMomentOfInertia[I31] = I[3];
    trim->weightMomentOfInertia[I32] = I[4];
    trim->weightMomentOfInertia[I33] = I[5];
  }

  status = json_getString(jsonDict, "accUnit", &trim->accUnit);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "missing \"accUnit\" entry required "
                       "by Trim discipline");
    return status;
  }

  status = json_getStringDynamicArray(jsonDict, "variables", &trim->numVariables, &trim->variables);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "missing \"variables\" entry required "
                       "by Trim discipline");
    return status;
  }

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(keyValue);
  return status;
}

static int _findTrimVarByName(void *aimInfo, const char *matchName,
                              const zaeroProblemStruct *zaeroProblem,
                              zaeroTrimVarStruct **matchTrimVar)
{
  int i;

  for (i = 0; i < zaeroProblem->numTrimVariables; i++) {
    if (strcmp(zaeroProblem->trimVariables[i].name, matchName) == 0) {
      *matchTrimVar = &zaeroProblem->trimVariables[i];
      return CAPS_SUCCESS;
    }
  }

  AIM_ERROR(aimInfo, "No trim variable '%s' in Trim_Variable!", matchName);
  AIM_ADDLINE(aimInfo, "Trim_Variable names:");
  for (i = 0; i < zaeroProblem->numTrimVariables; i++)
    AIM_ADDLINE(aimInfo, "    %s", zaeroProblem->trimVariables[i].name);

  // printf("Error: No Trim variable found with name: %s\n", matchName);
  return CAPS_NOTFOUND;
}

// create (write cards for) trim analysis from a zaeroSubcase
int zaero_trimAnalysis(void *aimInfo,
                       FILE *fp,
                       const zaeroSubcaseStruct *subcase,
                       const zaeroProblemStruct *zaeroProblem,
                       const feaFileTypeEnum formatType)
{
  int i, status;

  int trimID, mkaerozID, numVars, *varIDs=NULL;

  /*@-observertrans@*/
  char *none = "none", *varLabel = NULL, **varValues=NULL;
  char *longAccel = none, *latAccel = none, *vertAccel = none;
  char *rollAccel = none, *pitchAccel = none, *yawAccel = none;
  char *dc[6] = {none, none, none, none, none, none};
  /*@+observertrans@*/

  double defaultLower = -1.0e30;
  double defaultUpper =  1.0e30;

  zaeroUAICStruct *uaic = NULL;
  zaeroTrimVarStruct *trimvar = NULL;

  zaeroTrimStruct *trim = (zaeroTrimStruct *) subcase->discipline;

  status = zaero_findUAICByName(aimInfo, subcase->uaicName, zaeroProblem, &uaic);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(uaic, aimInfo, status);

  mkaerozID = uaic->id;
  trimID = subcase->analysisID;

  AIM_ALLOC(varIDs   , trim->numVariables, int  , aimInfo, status );
  AIM_ALLOC(varValues, trim->numVariables, char*, aimInfo, status );

  // write trim variables
  numVars = 0;
  for (i = 0; i < trim->numVariables; i++) {

    status = _findTrimVarByName(aimInfo, trim->variables[i], zaeroProblem, &trimvar);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(trimvar, aimInfo, status);

    varLabel = trimvar->label;

    if (strcmp(varLabel, "NX") == 0) {
      longAccel = trimvar->value;
    } else if (strcmp(varLabel, "NY") == 0) {
      latAccel = trimvar->value;
    } else if (strcmp(varLabel, "NZ") == 0) {
      vertAccel = trimvar->value;
    } else if (strcmp(varLabel, "PDOT") == 0) {
      rollAccel = trimvar->value;
    } else if (strcmp(varLabel, "QDOT") == 0) {
      pitchAccel = trimvar->value;
    } else if (strcmp(varLabel, "RDOT") == 0) {
      yawAccel = trimvar->value;
    } else {

      varIDs[i] = trimID + numVars++;
      varValues[i] = trimvar->value;

      status = zaero_card_trimvar(aimInfo,
          fp,
          varIDs[i], // idvar
          varLabel, // label
          defaultLower, // lower
          defaultUpper, // upper
          0,  // trimlnk, TODO: hardcoded
          NULL, // dmi, TODO: hardcoded
          "sym", // sym, TODO: hardcoded
          (int) false, // guess `initial` ?
          0, // initial, ignored if `guess` is false, TODO: hardcoded
          dc, // dcd, dcy, dcl, dcr, dcm, dcn
          formatType);
      AIM_STATUS(aimInfo, status);
    }
  }

  status = zaero_card_trim(aimInfo,
       fp,
       trimID, // trimid
       mkaerozID, // idmk
       trim->dynamicPressure, // qinf
       0, // idobj, TODO: hardcoded
       0, // idcons, TODO: hardcoded
       trim->vectorToCG, // rhox, rhoy, rhoz
       1 / trim->gravityAcceleration, // wtmass
       trim->weight, // weight
       trim->weightMomentOfInertia,// ixx, ixy, iyy, ixz, iyz, izz
       trim->accUnit, // trnacc
       longAccel, latAccel, vertAccel, // nx, ny, nz
       rollAccel, pitchAccel, yawAccel, // pdot, qdot, rdot
       0, // loadset, TODO: hardcoded
       numVars,
       varIDs, // idvari
       varValues, // vali
       formatType);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(varIDs);
  AIM_FREE(varValues);

  return status;
}
