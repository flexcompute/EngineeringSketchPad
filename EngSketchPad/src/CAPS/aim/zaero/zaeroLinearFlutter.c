#include <string.h>

#include "aimUtil.h"
#include "miscUtils.h"
#include "jsonUtils.h"
#include "zaeroUtils.h"
#include "zaeroCards.h"
#include "zaeroDiscipline.h"

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#define PATH_MAX   _MAX_PATH
#endif

static int _getFlutterFunction(void *aimInfo,
                               char *jsonDict, char *jsonKey,
                               zaeroFlutterFuncEnum *function) {

    int status = CAPS_SUCCESS;
    char *functionStr = NULL;

    status = json_getString(jsonDict, jsonKey, &functionStr);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(functionStr, aimInfo, status);

    if (strcasecmp(functionStr, "FixMach") == 0) {

        *function = FixMach;
    }
    else if (
        (strcasecmp(functionStr, "FixMachDensity") == 0) ||
        (strcasecmp(functionStr, "fixmden") == 0)) {

        *function = FixMachDensity;
    }
    else if (
        (strcasecmp(functionStr, "FixMachAtmos") == 0) ||
        (strcasecmp(functionStr, "fixmatm") == 0)) {

        *function = FixMachAtmos;
    }
    else if (
        (strcasecmp(functionStr, "FixAltAtmos") == 0) ||
        (strcasecmp(functionStr, "fixhatm") == 0)) {

        *function = FixAltAtmos;
    }
    else {
        printf("Error: Unrecognized flutter function.\n");
        return CAPS_BADVALUE;
    }

    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(functionStr);
    return status;
}

static int _getFixMachDensity(char *jsonDict,
                              zaeroLinearFlutterStruct *flutter) {

    int status;

    status = json_getDouble(jsonDict, "density", &flutter->density);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getDouble(jsonDict, "refVelocity", &flutter->refVelocity);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getDoubleDynamicArray(
        jsonDict, "velocities",
        &flutter->numVelocities, &flutter->velocities);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getInteger(
        jsonDict, "printFlutter", &flutter->printFlag);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    return CAPS_SUCCESS;
}

static int _getModalDamping(char *jsonDict, zaeroLinearFlutterStruct *flutter) {

    int status;

    status = json_getDoubleDynamicArray(
        jsonDict, "dampingFreqs",
        &flutter->numDampingFreq, &flutter->dampingFreq);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getDoubleDynamicArray(
        jsonDict, "dampingValues",
        &flutter->numDampingValues, &flutter->dampingValues);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getString(
        jsonDict, "dampingType", &flutter->dampingUnits);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    } else {
        string_toUpperCase(flutter->dampingUnits);
    }

    return CAPS_SUCCESS;
}

int zaero_getLinearFlutterDiscipline(void *aimInfo,
                                     char *jsonDict,
                                     zaeroSubcaseStruct *subcase) {

    int status; // Function return

    zaeroLinearFlutterStruct *flutter;

    // Ensure we are starting with fresh discipline
    // if (subcase->discipline != NULL) {
    //   EG_free(subcase->discipline);
    // }
    if (subcase->discipline != NULL) {
        printf("Error: subcase discipline must be NULL\n");
        return CAPS_BADVALUE;
    }

    subcase->discipline = ((zaeroLinearFlutterStruct *)
                           EG_alloc(sizeof(zaeroLinearFlutterStruct)));

    if (subcase->discipline == NULL) {
        return EGADS_MALLOC;
    }

    flutter = (zaeroLinearFlutterStruct *) subcase->discipline;
    initiate_zaeroLinearFlutterStruct(flutter);

    status = json_getInteger(jsonDict, "numModes", &flutter->numModes);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    }

    status = json_getString(
        jsonDict, "boundaryCondition", &flutter->boundaryCondition);
    if (status != CAPS_SUCCESS) {
        // TODO: is it required?
    } else {
        string_toUpperCase(flutter->boundaryCondition);
    }

    status = _getFlutterFunction(aimInfo, jsonDict, "function", &flutter->function);
    if (status != CAPS_SUCCESS) {
        // required parameter
        AIM_ERROR(aimInfo, "Required LinearFlutter parameter \"function\" not specified.\n");
        return CAPS_NOTFOUND;
    }

    switch (flutter->function) {
      case FixMach:
        status = CAPS_NOTIMPLEMENT;
        break;
      case FixMachDensity:
        status = _getFixMachDensity(jsonDict, flutter);
        break;
      case FixMachAtmos:
        status = CAPS_NOTIMPLEMENT;
        break;
      case FixAltAtmos:
        status = CAPS_NOTIMPLEMENT;
        break;
    // should not be unknown here, including to avoid switch warning
      case UnknownFlutterFunction:
        status = CAPS_BADVALUE;
        break;
    }
    AIM_STATUS(aimInfo, status);

    status = _getModalDamping(jsonDict, flutter);
    if (status != CAPS_SUCCESS) return status;

    status = CAPS_SUCCESS;
cleanup:
    return status;
}

// create (write cards for) flutter analysis from a zaeroSubcase
int zaero_flutterAnalysis(void *aimInfo,
                          FILE *fp,
                          const zaeroSubcaseStruct *subcase,
                          const zaeroProblemStruct *zaeroProblem,
                          const cfdUnitsStruct *units,
                          const feaFileTypeEnum formatType)
{
  int status;

  int flutterID, mkaerozID, fixID, tabdmp1ID, mlistID, conmlstID;

  zaeroUAICStruct *uaic;

  zaeroLinearFlutterStruct *flutter = (
      (zaeroLinearFlutterStruct *) subcase->discipline);

  status = zaero_findUAICByName(aimInfo, subcase->uaicName, zaeroProblem, &uaic);
  AIM_STATUS(aimInfo, status);

  mkaerozID = uaic->id;
  flutterID = subcase->analysisID;
  fixID = subcase->analysisID;
  tabdmp1ID = subcase->analysisID;
  mlistID = -1; //blank for now
  conmlstID = -1; // blank for now

  status = zaero_card_flutter(
      aimInfo,
      fp,
      flutterID, // setid
      flutter->boundaryCondition, // sym
      fixID, // fix
      flutter->numModes, // nmode
      tabdmp1ID, // tabdmp
      mlistID, // mlist
      conmlstID, // conmlst
      -1, // nkstep, TODO: hardcoded
      formatType
  );
  AIM_STATUS(aimInfo, status);

  if (flutter->function == FixMach) {
    return CAPS_NOTIMPLEMENT;
  }
  else if (flutter->function == FixMachDensity) {
    status = zaero_card_fixmden(
        aimInfo,
        fp,
        fixID, // setid
        mkaerozID, // idmk
        flutter->density, // den
        units->mass != NULL ? units->mass : "kg", // ftmunit
        units->length != NULL ? units->length : "m", //ftlunit
        flutter->refVelocity, // vref
        0, // fluttf, TODO: hardcoded
        flutter->printFlag, // print
        flutter->numVelocities,
        flutter->velocities, // Vi
        formatType
    );
    AIM_STATUS(aimInfo, status);
  }
  else if (flutter->function == FixMachAtmos) {
    return CAPS_NOTIMPLEMENT;
  }
  else if (flutter->function == FixMachAtmos) {
    return CAPS_NOTIMPLEMENT;
  }
  else {
    AIM_ERROR(aimInfo, "Unknown flutter function.\n");
    return CAPS_BADVALUE;
  }

  status = zaero_card_tabdmp1(
      aimInfo,
      fp,
      tabdmp1ID, // tid
      flutter->dampingUnits, // type
      flutter->numDampingFreq,
      flutter->dampingFreq, // Fi
      flutter->dampingValues, // Gi
      formatType
  );
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;
cleanup:
  return status;
}
