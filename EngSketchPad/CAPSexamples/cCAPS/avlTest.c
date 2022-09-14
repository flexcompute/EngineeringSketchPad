/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             avl AIM tester
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


static void
freeTuple(int tlen,          /* (in)  length of the tuple */
          capsTuple *tuple)  /* (in)  tuple freed */
{
    int i;

    if (tuple != NULL) {
        for (i = 0; i < tlen; i++) {

            EG_free(tuple[i].name);
            EG_free(tuple[i].value);
        }
    }
    EG_free(tuple);
}


int main(int argc, char *argv[])
{

    int status; // Function return status;
    int i; // Indexing
    int outLevel = 1;

    // CAPS objects
    capsObj          problemObj, avlObj, tempObj;
    capsErrs         *errors;

    // CAPS return values
    int   nErr;
    enum capsvType   vtype;

    // Input values
    capsTuple *surfaceTuple = NULL, *flapTuple = NULL;
    int       surfaceSize = 2, flapSize = 3;
    double    doubleVal;
    //enum capsBoolean  boolVal;

    // Output values
    int nrow, ncol;
    const char *valunits;
    const void *data;
    const int  *partial;

    char *stringVal = NULL;

    int major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char *analysisPath, *unitSystem, *intents, **fnames;
    char currentPath[PATH_MAX];

    printf("\n\nAttention: avlTest is hard coded to look for ../csmData/avlWingTail.csm\n");

    if (argc > 2) {
        printf(" usage: avlTest outLevel\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }

    status = caps_open("AVL_Example", NULL, 0, "../csmData/avlWingTail.csm",
                       outLevel, &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Now load the avlAIM (disabled auto execution)
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "avlAIM", NULL, NULL, NULL, &exec,
                               &avlObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set AVL_Surface
    status = caps_childByName(avlObj, VALUE, ANALYSISIN, "AVL_Surface", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    surfaceTuple = (capsTuple *) EG_alloc(surfaceSize*sizeof(capsTuple));
    surfaceTuple[0].name = EG_strdup("Wing");
    surfaceTuple[0].value = EG_strdup("{\"numChord\": 8, \"spaceChord\": 1, \"numSpanPerSection\": 12, \"spaceSpan\": 1}");

    surfaceTuple[1].name = EG_strdup("Vertical_Tail");
    surfaceTuple[1].value = EG_strdup("{\"numChord\": 5, \"spaceChord\": 1, \"numSpanTotal\": 10, \"spaceSpan\": 1}");

    status = caps_setValue(tempObj, Tuple, surfaceSize, 1,
                           (void **) surfaceTuple, NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    // Find & set AVL_Control
    status = caps_childByName(avlObj, VALUE, ANALYSISIN, "AVL_Control", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    flapTuple = (capsTuple *) EG_alloc(flapSize*sizeof(capsTuple));
    flapTuple[0].name  = EG_strdup("WingRightLE");
    flapTuple[0].value = EG_strdup("{\"controlGain\": 0.5, \"deflectionAngle\": 25}");

    flapTuple[1].name  = EG_strdup("WingRightTE");
    flapTuple[1].value = EG_strdup("{\"controlGain\": 1.0, \"deflectionAngle\": 15}");

    flapTuple[2].name  = EG_strdup("Tail");
    flapTuple[2].value = EG_strdup("{\"controlGain\": 1.0, \"deflectionAngle\": 15}");

    status = caps_setValue(tempObj, Tuple, flapSize, 1,  (void **) flapTuple,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Mach number
    status = caps_childByName(avlObj, VALUE, ANALYSISIN, "Mach", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = 0.5;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) &doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set AoA
    status = caps_childByName(avlObj, VALUE, ANALYSISIN, "Alpha", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = 1.0;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) &doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Do the analysis setup for AVL;
    status = caps_preAnalysis(avlObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* get analysis info */
    status = caps_analysisInfo(avlObj, &analysisPath, &unitSystem, &major,
                               &minor, &intents, &nFields, &fnames, &ranks,
                               &fInOut, &exec, &dirty);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Execute avl via system call
    (void) getcwd(currentPath, PATH_MAX);

    if (chdir(analysisPath) != 0) {
        status = CAPS_DIRERR;
        printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
        goto cleanup;
    }

    printf("\n\nRunning AVL!\n\n");

#ifdef WIN32
    system("avl.exe caps < avlInput.txt > avlOutput.txt");
#else
    system("avl caps < avlInput.txt > avlOutput.txt");
#endif

    (void) chdir(currentPath);

    status = caps_postAnalysis(avlObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Get  total Cl
    status = caps_childByName(avlObj, VALUE, ANALYSISOUT, "CLtot", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_getValue(tempObj, &vtype, &nrow, &ncol, &data, &partial,
                           &valunits, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    printf("\nValue of Cltot = %f\n", ((double *) data)[0]);

    // Get strip forces
    status = caps_childByName(avlObj, VALUE, ANALYSISOUT, "StripForces", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_getValue(tempObj, &vtype, &nrow, &ncol, &data, &partial,
                           &valunits, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    printf("\nStripForces\n\n");
    printf("%s = %s\n\n", ((capsTuple *) data)[0].name,
                          ((capsTuple *) data)[0].value);
    printf("%s = %s\n\n", ((capsTuple *) data)[1].name,
                          ((capsTuple *) data)[1].value);


cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d\n",
                                       status);

    freeTuple(surfaceSize, surfaceTuple);
    freeTuple(flapSize, flapTuple);

    EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
