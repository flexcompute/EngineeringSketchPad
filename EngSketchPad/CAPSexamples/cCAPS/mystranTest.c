/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             mystran AIM tester
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "caps.h"
#include <string.h>

#ifdef WIN32
#define unlink     _unlink
#define getcwd     _getcwd
#define PATH_MAX   _MAX_PATH
#define strcasecmp  stricmp
#define snprintf   _snprintf
#else
#include <limits.h>
#include <unistd.h>
#endif


static void
printErrors(int nErr, capsErrs *errors)
{
  int         i, j, stat, eType, nLines;
  char        **lines;
  capsObj     obj;
  static char *type[5] = {"Cont:   ", "Info:   ", "Warning:", "Error:  ",
                          "Status: "};

  if (errors == NULL) return;

  for (i = 1; i <= nErr; i++) {
    stat = caps_errorInfo(errors, i, &obj, &eType, &nLines, &lines);
    if (stat != CAPS_SUCCESS) {
      printf(" printErrors: %d/%d caps_errorInfo = %d\n", i, nErr, stat);
      continue;
    }
    for (j = 0; j < nLines; j++) {
      if (j == 0) {
        printf(" CAPS %s ", type[eType+1]);
      } else {
        printf("               ");
      }
      printf("%s\n", lines[j]);
    }
  }
  
  caps_freeError(errors);
}


int main(int argc, char *argv[])
{

    int status; // Function return status;
    int i; // Indexing
    int outLevel = 1;

    // CAPS objects
    capsObj  problemObj, mystranObj, tempObj;
    capsErrs *errors;

    // CAPS return values
    int   nErr;

    // Input values
    capsTuple        *material = NULL, *constraint = NULL, *property = NULL, *load = NULL;
    int               numMaterial = 1, numConstraint = 1, numProperty = 2, numLoad = 1;

    int               intVal;
    //double            doubleVal;
    enum capsBoolean  boolVal;
    char *stringVal = NULL;

    int major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char *analysisPath, *unitSystem, *intents, **fnames;
    char currentPath[PATH_MAX];

    printf("\n\nAttention: mystranTest is hard coded to look for ../csmData/aeroelasticDataTransferSimple.csm\n");
    printf("To not make system calls to the mystran executable the third command line "
           "option may be supplied - 0 = no analysis , >0 run analysis (default).\n\n");

    if (argc > 2) {
        printf(" usage: mystranTest outLevel!\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }

    status = caps_open("MyStran_Example", NULL, 0,
                       "../csmData/aeroelasticDataTransferSimple.csm", outLevel,
                       &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Load the AIMs
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "mystranAIM", NULL, NULL, NULL, &exec,
                               &mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Set Mystran inputs - Materials
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Material", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    material = (capsTuple *) EG_alloc(numMaterial*sizeof(capsTuple));
    material[0].name  = EG_strdup("Madeupium");
    material[0].value = EG_strdup("{\"youngModulus\": 2.2E6, \"density\": 7850}");

    status = caps_setValue(tempObj, Tuple, numMaterial, 1,  (void **) material,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                       - Properties
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Property", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    property = (capsTuple *) EG_alloc(numProperty*sizeof(capsTuple));
    property[0].name  = EG_strdup("Skin");
    property[0].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.1}");
    property[1].name  = EG_strdup("Rib_Root");
    property[1].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.2}");

    status = caps_setValue(tempObj, Tuple, numProperty, 1,  (void **) property,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                       - Constraints
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Constraint", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    constraint = (capsTuple *) EG_alloc(numConstraint*sizeof(capsTuple));
    constraint[0].name  = EG_strdup("edgeConstraint");
    constraint[0].value = EG_strdup("{\"groupName\": \"Rib_Root\", \"dofConstraint\": 123456}");

    status = caps_setValue(tempObj, Tuple, numConstraint, 1,
                           (void **) constraint, NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                       - Analysis
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Load", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    load = (capsTuple *) EG_alloc(numLoad*sizeof(capsTuple));
    load[0].name  = EG_strdup("appliedLoad");
    load[0].value = EG_strdup("{\"groupName\": \"Skin\", \"loadType\": \"Pressure\", \"pressureForce\": 2.0E6}");

    status = caps_setValue(tempObj, Tuple, numLoad, 1,  (void **) load,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Edge_Point_Max", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 3;
    status = caps_setValue(tempObj, Integer, 1, 1, (void *) &intVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Edge_Point_Min", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 2;
    status = caps_setValue(tempObj, Integer, 1, 1, (void *) &intVal, NULL, NULL,
                           &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Analysis_Type", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    stringVal = EG_strdup("Static");
    status = caps_setValue(tempObj, String, 1, 1, (void *) stringVal, NULL, NULL,
                           &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Quad_Mesh", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = (int) true;
    status = caps_setValue(tempObj, Boolean, 1, 1, (void *) &boolVal, NULL, NULL,
                           &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Run mystran pre-analysis
    status = caps_preAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* get analysis info */
    status = caps_analysisInfo(mystranObj, &analysisPath, &unitSystem, &major,
                               &minor, &intents, &nFields, &fnames, &ranks,
                               &fInOut, &exec, &dirty);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Execute Mystran via system call
    (void) getcwd(currentPath, PATH_MAX);

    if (chdir(analysisPath) != 0) {
        status = CAPS_DIRERR;
        printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
        goto cleanup;
    }

    printf("\n\nRunning mystran!\n\n");
    system("mystran.exe mystran_CAPS.dat > mystranOutput.txt");

    (void) chdir(currentPath);

    // Run Mystran post
    status = caps_postAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);

cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d\n",
                                       status);

    if (constraint != NULL) {
        for (i = 0; i < numConstraint; i++) {
            if (constraint[i].name  != NULL) EG_free(constraint[i].name);
            if (constraint[i].value != NULL) EG_free(constraint[i].value);
        }
        EG_free(constraint);
    }

    if (load != NULL) {
        for (i = 0; i < numLoad; i++) {
            if (load[i].name  != NULL) EG_free(load[i].name);
            if (load[i].value != NULL) EG_free(load[i].value);
        }
        EG_free(load);
    }

    if (property != NULL) {
        for (i = 0; i < numProperty; i++) {
            if (property[i].name  != NULL) EG_free(property[i].name);
            if (property[i].value != NULL) EG_free(property[i].value);
        }
        EG_free(property);
    }

    if (material != NULL) {
        for (i = 0; i < numMaterial; i++) {
            if (material[i].name  != NULL) EG_free(material[i].name);
            if (material[i].value != NULL) EG_free(material[i].value);
        }
        EG_free(material);
    }

    if (stringVal != NULL) EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
