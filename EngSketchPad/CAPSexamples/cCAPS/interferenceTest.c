/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             interference AIM tester
 *
 *      Copyright 2020-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "caps.h"
#include <string.h>


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
printResults(capsObj analysis)
{
  int            status, i, j, k, n, nErr, nrow, ncol;
  const int      *partial;
  const char     *names, *units;
  const double   *dists, *vols;
  capsObj        namesObj, distObj, volObj;
  capsErrs       *errors;
  enum capsvType vtype;
  
  status = caps_childByName(analysis, VALUE, ANALYSISOUT, "Names", &namesObj,
                            &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) {
    printf(" Error: Cannot get Names!\n");
    return;
  }
  status = caps_childByName(analysis, VALUE, ANALYSISOUT, "Distances", &distObj,
                            &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) {
    printf(" Error: Cannot get Distances!\n");
    return;
  }
  status = caps_childByName(analysis, VALUE, ANALYSISOUT, "Volumes", &volObj,
                            &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) {
    printf(" Error: Cannot get Volumes!\n");
    return;
  }
  
  status = caps_getValue(namesObj, &vtype, &nrow, &ncol, (const void **) &names,
                         &partial, &units, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) return;
  status = caps_getValue(distObj, &vtype, &nrow, &ncol, (const void **) &dists,
                         &partial, &units, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) return;
  n = nrow;
  status = caps_getValue(volObj, &vtype, &nrow, &ncol, (const void **) &vols,
                         &partial, &units, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) return;

  for (j = i = 0; i < n; i++) {
    printf(" %2d: %13.6le  %s\n", i+1, vols[i], &names[j]);
    j += strlen(&names[j]) + 1;
  }
  printf("\n");
  
  for (j = i = 0; i < n; i++) {
    for (k = 0; k < n; k++, j++)
      printf(" %13.6le ", dists[j]);
    printf("\n");
  }
  printf("\n");
}


int main(int argc, char *argv[])
{

  int      i, status, nErr, intVal, exec;
  capsObj  problemObj, interfereObj, tempObj;
  capsErrs *errors;

  printf("\n\nNote: interferenceTest uses ../csmData/interference.csm\n");

  status = caps_open("Interference_Example", NULL, 0,
                     "../csmData/interference.csm", 1, &problemObj,
                     &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;

  /* Load the AIM */
  exec   = 0;
  status = caps_makeAnalysis(problemObj, "interferenceAIM", NULL, NULL, NULL,
                             &exec, &interfereObj, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;
  
  /* first use all default inputs */
  status = caps_execute(interfereObj, &exec, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;
  printResults(interfereObj);
  
  /* now only do the inners */
  intVal = False;
  status = caps_childByName(interfereObj, VALUE, ANALYSISIN, "OML", &tempObj,
                            &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
   if (status != CAPS_SUCCESS) goto cleanup;
  status = caps_setValue(tempObj, Boolean, 1, 1, (void *) &intVal, NULL, NULL,
                         &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;
  status = caps_childByName(interfereObj, VALUE, ANALYSISIN, "Attr_Name", &tempObj,
                            &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
   if (status != CAPS_SUCCESS) goto cleanup;
  status = caps_setValue(tempObj, String, 1, 1, (void *) "Inner", NULL, NULL,
                         &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;
  
  status = caps_execute(interfereObj, &exec, &nErr, &errors);
  if (nErr != 0) printErrors(nErr, errors);
  if (status != CAPS_SUCCESS) goto cleanup;
  printResults(interfereObj);

cleanup:
  if (status != CAPS_SUCCESS)
    printf("\n\nPremature exit - status = %d\n", status);

  i = 0;
  if (status == CAPS_SUCCESS) i = 1;
  (void) caps_close(problemObj, i, NULL);

  return 0;
}
