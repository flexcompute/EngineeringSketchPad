/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Friction AIM tester
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
#define strcasecmp stricmp
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
    capsObj          problemObj, frictionObj, tempObj;
    capsErrs         *errors;

    // CAPS return values
    int   nErr, major, minor, nfields, exec, stat, *frank, *fInOut;
    char *usys, *intent, **fnames;
    enum capsvType   vtype;

    // Input values
    double            doubleVal[2];

    // Output values
    int nrow, ncol;
    const int  *partial;
    const char *valunits;
    const void *data;

    char *stringVal = NULL;

    char *analysisPath;
    char currentPath[PATH_MAX];

    printf("\n\nAttention: frictionTest is hard coded to look for ../csmData/frictionWingTailFuselage.csm\n");
    printf("To not make system calls to the friction executable the second command line option may be "
           "supplied - 0 = no analysis , >0 run analysis (default).\n\n");

    if (argc > 2) {
        printf(" usage: frictionTest outLevel!\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }


    status = caps_open("FRICTION_Example", NULL, 0,
                       "../csmData/frictionWingTailFuselage.csm", outLevel,
                       &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Now load the frictionAIM
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "frictionAIM", NULL, NULL, NULL,
                               &exec, &frictionObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;
  
    status = caps_analysisInfo(frictionObj, &analysisPath, &usys, &major, &minor,
                               &intent, &nfields, &fnames, &frank, &fInOut,
                               &exec, &stat);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Mach number
    status = caps_childByName(frictionObj, VALUE, ANALYSISIN, "Mach", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal[0] = 0.5;
    doubleVal[1] = 1.5;
    status = caps_setValue(tempObj, Double, 2, 1, (void *) &doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    // Find & set Altitude
    status = caps_childByName(frictionObj, VALUE, ANALYSISIN, "Altitude", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal[0] = 29.52756; // kft
    doubleVal[1] = 59.711286;
    status = caps_setValue(tempObj, Double, 2, 1, (void *) &doubleVal, NULL,
                           "kft", &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    // Do the analysis setup for FRICTION;
    status = caps_preAnalysis(frictionObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Execute friction via system call
    (void) getcwd(currentPath, PATH_MAX);

    if (chdir(analysisPath) != 0) {
        status = CAPS_DIRERR;
        printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
        goto cleanup;
    }

    printf("\n\nRunning FRICTION!\n\n");

#ifdef WIN32
    system("friction.exe frictionInput.txt frictionOutput.txt > Info.out");
#else
    system("friction frictionInput.txt frictionOutput.txt > Info.out");
#endif


    (void) chdir(currentPath);

    status = caps_postAnalysis(frictionObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Get  total Cl
    status = caps_childByName(frictionObj, VALUE, ANALYSISOUT, "CDfric", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_getValue(tempObj, &vtype, &nrow, &ncol, &data, &partial,
                           &valunits, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < nrow*ncol; i++) {
        printf("\nValue of CDfric = %lf\n", ((double *) data)[i]);
    }

cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d\n",
                                       status);

    if (stringVal != NULL) EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
