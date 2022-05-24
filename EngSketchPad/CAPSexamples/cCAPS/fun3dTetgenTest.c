/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             fun3d/tetgen AIM tester
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
    capsObj          problemObj, surfMeshObj, meshObj, fun3dObj, tempObj;
    capsErrs         *errors;

    // CAPS return values
    int   nErr;
    capsObj source, target;

    // Input values
    capsTuple        *tupleVal = NULL;
    int               tupleSize = 3;
    double            doubleVal;
    enum capsBoolean  boolVal;

    char *stringVal = NULL;

    int major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char *analysisPath, *unitSystem, *intents, **fnames;
    char currentPath[PATH_MAX];

    printf("\n\nAttention: fun3dTetgenTest is hard coded to look for ../csmData/cfdMultiBody.csm\n");

    if (argc > 2) {
        printf(" usage: fun3dTetgenTest outLevel!\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }

    status = caps_open("FUN3D_Tetgen_Example", NULL, 0,
                       "../csmData/cfdMultiBody.csm", outLevel, &problemObj,
                       &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Load the EGADS Tess AIM */
    exec   = 1;
    status = caps_makeAnalysis(problemObj, "egadsTessAIM", NULL, NULL, NULL,
                               &exec, &surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    // EGADS executes automatically

    // Load the Tetgen AIM */
    exec   = 1;
    status = caps_makeAnalysis(problemObj, "tetgenAIM", NULL, NULL, NULL,
                               &exec, &meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    // Link surface mesh from EGADS to TetGen
    status = caps_childByName(surfMeshObj, VALUE, ANALYSISOUT, "Surface_Mesh", &source,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("surfMeshObj childByName for Surface_Mesh = %d\n", status);
      goto cleanup;
    }
    status = caps_childByName(meshObj, VALUE, ANALYSISIN, "Surface_Mesh", &target,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("meshObj childByName for tessIn = %d\n", status);
      goto cleanup;
    }
    status = caps_linkValue(source, Copy, target, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf(" caps_linkValue = %d\n", status);
      goto cleanup;
    }

    // Tetgen executes automatically

    // Now load the fun3dAIM
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "fun3dAIM", NULL, NULL, NULL,
                               &exec, &fun3dObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Boundary_Conditions
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Boundary_Condition", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    tupleVal = (capsTuple *) EG_alloc(tupleSize*sizeof(capsTuple));
    tupleVal[0].name  = EG_strdup("Wing1");
    tupleVal[0].value = EG_strdup("{\"bcType\": \"Inviscid\", \"wallTemperature\": 1}");

    tupleVal[1].name  = EG_strdup("Wing2");
    tupleVal[1].value =  EG_strdup("{\"bcType\": \"Inviscid\", \"wallTemperature\": 1.2}");

    tupleVal[2].name  = EG_strdup("Farfield");
    tupleVal[2].value = EG_strdup("farfield");

    status = caps_setValue(tempObj, Tuple, tupleSize, 1,  (void **) tupleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Mach number input
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Mach", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = 0.4;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) &doubleVal, NULL,
                           NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Overwrite_NML */
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Overwrite_NML", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = true;
    status = caps_setValue(tempObj, Boolean, 1, 1, (void *) &boolVal, NULL,
                           NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Link the volume mesh from TetGen to Fun3D
    status = caps_childByName(meshObj, VALUE, ANALYSISOUT, "Volume_Mesh", &source,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("meshObj childByName for Volume_Mesh = %d\n", status);
      goto cleanup;
    }
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Mesh",  &target,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("fun3dObj childByName for Mesh = %d\n", status);
      goto cleanup;
    }
    status = caps_linkValue(source, Copy, target, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf(" caps_linkValue = %d\n", status);
      goto cleanup;
    }

    // Do the analysis setup for FUN3D;
    status = caps_preAnalysis(fun3dObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* get analysis info */
    status = caps_analysisInfo(fun3dObj, &analysisPath, &unitSystem, &major,
                               &minor, &intents, &nFields, &fnames, &ranks,
                               &fInOut, &exec, &dirty);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Execute FUN3D via system call
    (void) getcwd(currentPath, PATH_MAX);

    if (chdir(analysisPath) != 0) {
        status = CAPS_DIRERR;
        printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
        goto cleanup;
    }

    printf(" NOT Running fun3d!\n");
    //system("nodet_mpi > fun3dOutput.txt");

    (void) chdir(currentPath);

    status = caps_postAnalysis(fun3dObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);

cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d",
                                       status);

    if (tupleVal != NULL) {
        for (i = 0; i < tupleSize; i++) {
            if (tupleVal[i].name  != NULL) EG_free(tupleVal[i].name);
            if (tupleVal[i].value != NULL) EG_free(tupleVal[i].value);
        }
        EG_free(tupleVal);
    }
    if (stringVal != NULL) EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
