/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             mses AIM tester
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "caps.h"
#include <string.h>

#ifdef WIN32
#define unlink   _unlink
#define getcwd   _getcwd
#define PATH_MAX _MAX_PATH
#define strcasecmp  stricmp
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
    capsObj          problemObj, msesObj, tempObj;
    capsErrs         *errors=NULL;

    // CAPS return values
    int   nErr;
    enum capsvType   vtype;

    // Input values
    double    doubleVal[2];
    //enum capsBoolean  boolVal;

    // Output values
    int nrow, ncol;
    const char *valunits;
    const void *data;
    const int  *partial;

    char *stringVal = NULL;

    int exec;

    printf("\n\nAttention: msesTest is hard coded to look for ../csmData/airfoilSection.csm\n");

    if (argc > 2) {
        printf(" usage: msesTest outLevel\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }

    status = caps_open("mses_Example", NULL, 0, "../csmData/airfoilSection.csm",
                       outLevel, &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Now load the msesAIM (with auto execution)
    exec   = 1;
    status = caps_makeAnalysis(problemObj, "msesAIM", NULL, NULL, NULL, &exec,
                               &msesObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Mach number
    status = caps_childByName(msesObj, VALUE, ANALYSISIN, "Mach", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal[0] = 0.5;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set AoA
    status = caps_childByName(msesObj, VALUE, ANALYSISIN, "Alpha", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal[0] = 2.0;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // mses analysis executes automatically

    // Get CL
    status = caps_childByName(msesObj, VALUE, ANALYSISOUT, "CL", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_getValue(tempObj, &vtype, &nrow, &ncol, &data, &partial,
                           &valunits, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    printf("\nValue of CL = %f \n", ((double *) data)[0]);

cleanup:
    if (status != CAPS_SUCCESS) {
        printf("\n\nPremature exit - status = %d\n",
                                       status);
    }

    EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
