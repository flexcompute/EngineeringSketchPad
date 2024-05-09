#include <string.h>

#include "aimUtil.h"
#include "miscUtils.h"
#include "jsonUtils.h"
#include "zaeroUtils.h"
#include "zaeroCards.h"
#include "zaeroGraphical.h"


int zaero_textFileGenerationAero(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                    /*@unused@*/ const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID)
{
  int status;
  int i;

  int numFilenames=0, numForms=0;
  char **filenames = NULL, **forms = NULL;

  // make sure jsonDict is json string
  if (!json_isDict(jsonDict)) {
    AIM_ERROR(aimInfo, "Output Aero tuple value is expected to be a JSON dictionary string");
    return CAPS_BADVALUE;
  }

  status = json_getStringDynamicArray(jsonDict, "filename", &numFilenames, &filenames);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "Output Aero tuple missing 'filename'");
    return CAPS_BADVALUE;
  }
  AIM_NOTNULL(filenames, aimInfo, status);

  status = json_getStringDynamicArray(jsonDict, "form", &numForms, &forms);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "Output Aero tuple missing 'form'");
    return CAPS_BADVALUE;
  }
  AIM_NOTNULL(forms, aimInfo, status);

  if (!zaero_allEqual(Integer, &numFilenames, &numForms, NULL)) {
    AIM_ERROR(aimInfo, "Error in Output Aero input, not all values have same length.");
    return CAPS_BADVALUE;
  }

  for (i = 0; i < numFilenames; i++) {

    (*setID)++;

    status = zaero_card_pltaero(
        aimInfo,
        fp,
        (*setID), // setid
        "YES", // femgrid, TODO: hardcoded
        0, // offset, TODO: hardcoded
        forms[i], // form
        filenames[i], // filenm
        "NO", // cell, TODO: hardcoded
        "NO", // vct, TODO: hardcoded
        formatType
    );
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;
cleanup:
  string_freeArray(numFilenames, &filenames);
  string_freeArray(numForms, &forms);

  return status;
}

int zaero_textFileGenerationFlutter(void *aimInfo,
                                    FILE *fp,
                                    const char *jsonDict,
                                    const zaeroProblemStruct *zaeroProblem,
                                    const feaFileTypeEnum formatType,
                                    int *setID)
{
  int status;
  int i;

  int flutterID;
  int numAnalysis, numModes, numFilenames, numAeronames, numForms;
  int *modes;
  char **analysis = NULL, **filenames = NULL, **aeronames, **forms = NULL;
  zaeroSubcaseStruct *subcase;

  // make sure jsonDict is json string
  if (!json_isDict(jsonDict)) {
    AIM_ERROR(aimInfo, "Output Flutter tuple value is expected to be a JSON dictionary string\n");
    return CAPS_BADVALUE;
  }

  status = json_getStringDynamicArray(jsonDict, "analysis", &numAnalysis, &analysis);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getIntegerDynamicArray(jsonDict, "mode", &numModes, &modes);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "filename", &numFilenames, &filenames);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "form", &numForms, &forms);
  if (status != CAPS_SUCCESS) {
    numForms = numAnalysis;
    zaero_setBlankStringArray(numForms, &forms);
  }

  status = json_getStringDynamicArray(jsonDict, "aeroname", &numAeronames, &aeronames);
  if (status != CAPS_SUCCESS) {
    numAeronames = numAnalysis;
    zaero_setBlankStringArray(numAeronames, &aeronames);
  }

  if (!zaero_allEqual(Integer, &numFilenames, &numAeronames, &numForms, &numAnalysis, NULL)) {
    AIM_ERROR(aimInfo, "Error in Output Flutter input, not all values have same length.\n");
    return CAPS_BADVALUE;
  }

  for (i = 0; i < numFilenames; i++) {
    AIM_NOTNULL(analysis, aimInfo, status);
    AIM_NOTNULL(forms, aimInfo, status);
    AIM_NOTNULL(filenames, aimInfo, status);

    (*setID)++;
    status = zaero_findSubcaseByName(aimInfo, analysis[i], zaeroProblem, &subcase);
    AIM_STATUS(aimInfo, status);
    flutterID = subcase->analysisID;
    // TODO: assert(subcase->analysisType == AerolasticFlutter);

    status = zaero_card_pltflut(
        aimInfo,
        fp,
        (*setID), // setid
        flutterID, // idflut
        modes[i], // mode
        1, // ntime, TODO: hardcoded
        1.0, // maxdisp, TODO: hardcoded
        forms[i], // form
        filenames[i], // filenm
        aeronames[i], // aeronm
        formatType
    );
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;
cleanup:
  string_freeArray(numAnalysis, &analysis);
  string_freeArray(numFilenames, &filenames);
  string_freeArray(numAeronames, &aeronames);
  string_freeArray(numForms, &forms);
  AIM_FREE(modes);

  return status;
}

int zaero_textFileGenerationMode(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                    /*@unused@*/ const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID)
{
  int status;
  int i;

  int numModes, numFilenames, numForms, numAeronames;
  int *modes = NULL;
  char **filenames = NULL, **forms = NULL, **aeronames;

  // make sure jsonDict is json string
  if (!json_isDict(jsonDict)) {
    AIM_ERROR(aimInfo, "Output Mode tuple value is expected to be a JSON dictionary string\n");
    return CAPS_BADVALUE;
  }

  status = json_getIntegerDynamicArray(jsonDict, "mode", &numModes, &modes);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "filename", &numFilenames, &filenames);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "form", &numForms, &forms);
  if (status != CAPS_SUCCESS) {
    numForms = numModes;
    zaero_setBlankStringArray(numForms, &forms);
  }

  status = json_getStringDynamicArray(jsonDict, "aeroname", &numAeronames, &aeronames);
  if (status != CAPS_SUCCESS) {
    numAeronames = numModes;
    zaero_setBlankStringArray(numAeronames, &aeronames);
  }

  if (!zaero_allEqual(Integer, &numFilenames, &numForms, &numModes, &numAeronames, NULL)) {
    AIM_ERROR(aimInfo, "Error in Output Mode input, not all values have same length.\n");
    return CAPS_BADVALUE;
  }

  for (i = 0; i < numFilenames; i++) {
    AIM_NOTNULL(modes, aimInfo, status);
    AIM_NOTNULL(forms, aimInfo, status);
    AIM_NOTNULL(filenames, aimInfo, status);

    (*setID)++;

    status = zaero_card_pltmode(
        aimInfo,
        fp,
        (*setID), // setid
        "SYM", // sym, TODO: hardcoded
        modes[i], // mode
        1.0, // maxdisp, TODO: hardcoded
        forms[i], // form
        filenames[i], // filenm
        aeronames[i], // aeronm
        formatType
    );
    AIM_STATUS(aimInfo, status);
  }

cleanup:
  string_freeArray(numFilenames, &filenames);
  string_freeArray(numAeronames, &aeronames);
  string_freeArray(numForms, &forms);
  AIM_FREE(modes);

  return status;
}

int zaero_textFileGenerationTrim(void *aimInfo,
                                 FILE *fp,
                                 const char *jsonDict,
                                 const zaeroProblemStruct *zaeroProblem,
                                 const feaFileTypeEnum formatType,
                                 int *setID)
{
  int status;
  int i;

  int trimID;
  int numAnalysis, numFlexs, numTypes, numFilenames, numForms, numScales, numAeronames;
  double *scales = NULL;
  char **analysis, **filenames = NULL, **flexs = NULL, **types = NULL,
       **forms = NULL, **aeronames = NULL;
  zaeroSubcaseStruct *subcase;

  // make sure jsonDict is json string
  if (!json_isDict(jsonDict)) {
    AIM_ERROR(aimInfo, "Output Trim tuple value is expected to be a JSON dictionary string\n");
    return CAPS_BADVALUE;
  }

  status = json_getStringDynamicArray(jsonDict, "analysis", &numAnalysis, &analysis);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "flex", &numFlexs, &flexs);
  if (status != CAPS_SUCCESS) {
    numFlexs = numAnalysis;
    zaero_setBlankStringArray(numFlexs, &flexs);
  }

  status = json_getStringDynamicArray(jsonDict, "type", &numTypes, &types);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "Output Trim \"type\" value is required\n");
    return status;
  }

  status = json_getStringDynamicArray(jsonDict, "form", &numForms, &forms);
  if (status != CAPS_SUCCESS) {
    numForms = numAnalysis;
    zaero_setBlankStringArray(numForms, &forms);
  }

  status = json_getStringDynamicArray(jsonDict, "filename", &numFilenames, &filenames);
  if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "Output Trim \"filename\" value is required\n");
    return status;
  }

  status = json_getDoubleDynamicArray(jsonDict, "scale", &numScales, &scales);
  if (status != CAPS_SUCCESS) {
    numScales = numAnalysis;
    zaero_setBlankRealArray(numScales, &scales);
  }

  status = json_getStringDynamicArray(jsonDict, "aeroname", &numAeronames, &aeronames);
  if (status != CAPS_SUCCESS) {
    numAeronames = numAnalysis;
    zaero_setBlankStringArray(numAeronames, &aeronames);
  }

  if (!zaero_allEqual(Integer, &numAnalysis, &numScales, &numFilenames,
                      &numForms, &numFlexs, &numTypes, &numAeronames, NULL)) {
    AIM_ERROR(aimInfo, "Error in Output Trim input, not all values have same length.\n");
    return CAPS_BADVALUE;
  }

  for (i = 0; i < numFilenames; i++) {
    AIM_NOTNULL(flexs, aimInfo, status);
    AIM_NOTNULL(types, aimInfo, status);
    AIM_NOTNULL(forms, aimInfo, status);
    AIM_NOTNULL(filenames, aimInfo, status);
    AIM_NOTNULL(scales, aimInfo, status);
    AIM_NOTNULL(aeronames, aimInfo, status);

    (*setID)++;
    status = zaero_findSubcaseByName(aimInfo, analysis[i], zaeroProblem, &subcase);
    AIM_STATUS(aimInfo, status);
    trimID = subcase->analysisID;
    // TODO: assert(subcase->analysisType == AerolasticFlutter);

    status = zaero_card_plttrim(
        aimInfo,
        fp,
        (*setID), // setid
        trimID, // idtrim
        flexs[i], // flex
        types[i], // type
        forms[i], // form
        filenames[i], // filenm
        scales[i], // scale
        aeronames[i], // aeronm
        formatType
    );
    AIM_STATUS(aimInfo, status);

  }

  status = CAPS_SUCCESS;
cleanup:
  string_freeArray(numAnalysis, &analysis);
  string_freeArray(numFlexs, &flexs);
  string_freeArray(numTypes, &types);
  string_freeArray(numForms, &forms);
  string_freeArray(numFilenames, &filenames);
  string_freeArray(numAeronames, &aeronames);
  AIM_FREE(scales);
  // if (aeronames != NULL) EG_free(aeronames); // blank strings are NULL pointers

  return status;
}

int zaero_textFileGenerationVG(void *aimInfo,
                               FILE *fp,
                               const char *jsonDict,
                               const zaeroProblemStruct *zaeroProblem,
                               const feaFileTypeEnum formatType,
                               int *setID)
{
  int status;
  int i;

  int flutterID;
  int numAnalysis, numFilenames, numForms, numXAxis;
  char **analysis = NULL, **filenames = NULL, **forms = NULL, **xAxis = NULL;
  zaeroSubcaseStruct *subcase;

  // make sure jsonDict is json string
  if (!json_isDict(jsonDict)) {
    AIM_ERROR(aimInfo, "Output VG tuple value is expected "
                       "to be a JSON dictionary string\n");
    return CAPS_BADVALUE;
  }

  status = json_getStringDynamicArray(jsonDict, "analysis", &numAnalysis, &analysis);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "filename", &numFilenames, &filenames);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  status = json_getStringDynamicArray(jsonDict, "form", &numForms, &forms);
  if (status != CAPS_SUCCESS) {
    numForms = numAnalysis;
    zaero_setBlankStringArray(numForms, &forms);
  }

  status = json_getStringDynamicArray(jsonDict, "xAxis", &numXAxis, &xAxis);
  if (status != CAPS_SUCCESS) {
    // TODO: is it required?
  }

  if (!zaero_allEqual(Integer, &numFilenames, &numForms,
                      &numAnalysis, &numXAxis, NULL)) {
    AIM_ERROR(aimInfo, "Error in Output VG input, not all values have same length.\n");
    return CAPS_BADVALUE;
  }

  for (i = 0; i < numFilenames; i++) {
    AIM_NOTNULL(analysis, aimInfo, status);
    AIM_NOTNULL(xAxis, aimInfo, status);
    AIM_NOTNULL(forms, aimInfo, status);
    AIM_NOTNULL(filenames, aimInfo, status);

    (*setID)++;
    status = zaero_findSubcaseByName(aimInfo, analysis[i], zaeroProblem, &subcase);
    AIM_STATUS(aimInfo, status);
    flutterID = subcase->analysisID;
    // TODO: assert(subcase->analysisType == AerolasticFlutter);

    status = zaero_card_pltvg(
        aimInfo,
        fp,
        (*setID), // setid
        flutterID, // idflut
        -1, // nmode (blank), TODO: hardcoded
        xAxis[i], // ntime
        forms[i], // form
        filenames[i], // filenm
        -1.0, // refrho (blank), TODO: hardcoded
        formatType
    );
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;
cleanup:
  string_freeArray(numAnalysis, &analysis);
  string_freeArray(numFilenames, &filenames);
  string_freeArray(numXAxis, &xAxis);
  string_freeArray(numForms, &forms);

  return status;
}
