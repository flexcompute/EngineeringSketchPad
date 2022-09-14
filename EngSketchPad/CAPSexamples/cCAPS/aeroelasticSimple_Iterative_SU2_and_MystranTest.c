/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             SU2, tetgen, mystran AIM tester
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


/***********************************************************************/
/*                                                                     */
/*   helper functions                                                  */
/*                                                                     */
/***********************************************************************/


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


static int
setValueByName_Double(capsObj probObj,              /* (in)  problem object */
                      int     type,                 /* (in)  variable subtype */
                      char    name[],               /* (in)  variable name */
                      int     nrow,                 /* (in)  number of rows */
                      int     ncol,                 /* (in)  number of columns */
                      double *value)                /* (in)  value to set */
{
    int  status = CAPS_SUCCESS;

    int       nErr;
    capsErrs  *errors;

    capsObj          valObj;
    const int        *partial=NULL;
    const char       *units=NULL;

    /* --------------------------------------------------------------- */

    status = caps_childByName(probObj, VALUE, type, name, &valObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
        printf("caps_childByName(%s) -> status=%d\n", name, status);
        return status;
    }

    status = caps_setValue(valObj, Double, nrow, ncol, (void *) value, partial,
                           units, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)
        printf("caps_setValue(%s) -> status=%d\n", name, status);

    return status;
}


static int
setValueByName_Integer(capsObj probObj,              /* (in)  problem object */
                       int     type,                 /* (in)  variable type */
                       char    name[],               /* (in)  variable name */
                       int     nrow,                 /* (in)  number of rows */
                       int     ncol,                 /* (in)  number of columns */
                       int    *value)                /* (in)  value to set */
{
    int  status = CAPS_SUCCESS;

    int       nErr;
    capsErrs  *errors;

    capsObj          valObj;
    const int        *partial=NULL;
    const char       *units=NULL;

    /* --------------------------------------------------------------- */

    status = caps_childByName(probObj, VALUE, type, name, &valObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
        printf("caps_childByName(%s) -> status=%d\n", name, status);
        return status;
    }

    status = caps_setValue(valObj, Integer, nrow, ncol, (void *) value, partial,
                           units, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)
        printf("caps_setValue(%s) -> status=%d\n", name, status);

    return status;
}


static int
setValueByName_Boolean(capsObj           probObj,              /* (in)  problem object */
                       int               type,                 /* (in)  variable type */
                       char              name[],               /* (in)  variable name */
                       int               nrow,                 /* (in)  number of rows */
                       int               ncol,                 /* (in)  number of columns */
                       enum capsBoolean *value)                /* (in)  value to set */
{
    int  status = CAPS_SUCCESS;

    int       nErr;
    capsErrs  *errors;

    capsObj          valObj;
    const int        *partial=NULL;
    const char       *units=NULL;

    /* --------------------------------------------------------------- */

    status = caps_childByName(probObj, VALUE, type, name, &valObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
        printf("caps_childByName(%s) -> status=%d\n", name, status);
        return status;
    }

    status = caps_setValue(valObj, Boolean, nrow, ncol, (void *) value, partial,
                           units, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)
        printf("caps_setValue(%s) -> status=%d\n", name, status);

    return status;
}


static int
setValueByName_Tuple(capsObj    probObj,              /* (in)  problem object */
                     int        type,                 /* (in)  variable type */
                     char       name[],               /* (in)  variable name */
                     int        nrow,                 /* (in)  number of rows */
                     int        ncol,                 /* (in)  number of columns */
                     capsTuple *value)                /* (in)  value to set */
{
    int  status = CAPS_SUCCESS;

    int       nErr;
    capsErrs  *errors;

    capsObj          valObj;
    const int        *partial=NULL;
    const char       *units=NULL;

    /* --------------------------------------------------------------- */

    status = caps_childByName(probObj, VALUE, type, name, &valObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
        printf("caps_childByName(%s) -> status=%d\n", name, status);
        return status;
    }

    status = caps_setValue(valObj, Tuple, nrow, ncol, (void *) value, partial,
                           units, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)
        printf("caps_setValue(%s) -> status=%d\n", name, status);

    return status;
}


static int
setValueByName_String(capsObj     probObj,              /* (in)  problem object */
                      int         type,                 /* (in)  variable type */
                      char        name[],               /* (in)  variable name */
                      int         nrow,                 /* (in)  number of rows */
                      int         ncol,                 /* (in)  number of columns */
                      const char *value)                /* (in)  value to set */
{
    int  status = CAPS_SUCCESS;

    int       nErr;
    capsErrs  *errors;

    capsObj          valObj;
    const int        *partial=NULL;
    const char       *units=NULL;

    /* --------------------------------------------------------------- */

    status = caps_childByName(probObj, VALUE, type, name, &valObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
        printf("caps_childByName(%s) -> status=%d\n", name, status);
        return status;
    }

    status = caps_setValue(valObj, String, nrow, ncol, (void *) value, partial,
                           units, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)
        printf("caps_setValue(%s) -> status=%d\n", name, status);

    return status;
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


/***********************************************************************/
/*                                                                     */
/*   main - main program                                               */
/*                                                                     */
/***********************************************************************/

int main(int argc, char *argv[])
{

    int status; // Function return status;
    int i, iter; // Indexing
    int outLevel = 1;

    // CAPS objects
    capsObj  problemObj, surfMeshObj, meshObj, su2Obj, mystranObj;
    capsErrs *errors;

    // Data transfer
    const char *transfers[3] = {"Skin_Top", "Skin_Bottom", "Skin_Tip"};
    capsObj  boundObj[3], vertexSU2Obj[3], vertexMystranObj[3], pressureSU2Obj[3], pressureMystranObj[3], displacementMystranObj[3], displacementSU2Obj[3];
#ifdef DEBUG
    char filename[256];
#endif

    // CAPS return values
    int   nErr;
    capsObj source, target;

    // Input values
    capsTuple         *su2BC = NULL, *material = NULL, *constraint = NULL, *property = NULL, *load = NULL;
    int               numSU2BC = 3, numMaterial = 1, numConstraint = 1, numProperty = 2, numLoad = 1;

    int               intVal;
    double            doubleVal, tessParams[3];
    enum capsBoolean  boolVal;
    double            displacement0[3] = {0,0,0};

    int major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char *analysisPath=NULL, *unitSystem, *intents, **fnames;

    const char projectName[] = "aeroelasticSimple_Iterative";
    char currentPath[PATH_MAX];

    double speedofSound = 340.0, refVelocity = 100.0, refDensity = 1.2;

    int numIteration = 2;

    printf("\n\nAttention: aeroelasticIterativeTest is hard coded to look for ../csmData/aeroelasticDataTransferSimple.csm\n");

    if (argc > 2) {
        printf(" usage: aeroelasticSimple_Iterative_SU2_and_Mystran outLevel!\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }


    status = caps_open("SU2_MyStran_Aeroelastic_Interative_Example", NULL, 0,
                       "../csmData/aeroelasticDataTransferSimple.csm", outLevel,
                       &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* --------------------------------------------------------------- */

    // Load the AIMs
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "egadsTessAIM", NULL, NULL, NULL,
                               &exec, &surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "tetgenAIM", NULL, NULL, NULL,
                               &exec, &meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "su2AIM", NULL, NULL, NULL,
                               &exec, &su2Obj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "mystranAIM", NULL, NULL, NULL,
                               &exec, &mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* --------------------------------------------------------------- */

    // Create data bounds
    for (i = 0; i < 3; i++) {
        status = caps_makeBound(problemObj, 2, transfers[i], &boundObj[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeVertexSet(boundObj[i], su2Obj, NULL, &vertexSU2Obj[i],
                                    &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeVertexSet(boundObj[i], mystranObj, NULL,
                                    &vertexMystranObj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeDataSet(vertexSU2Obj[i], "Pressure", FieldOut, 0,
                                  &pressureSU2Obj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeDataSet(vertexMystranObj[i], "Pressure", FieldIn, 0,
                                  &pressureMystranObj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeDataSet(vertexMystranObj[i], "Displacement", FieldOut, 0,
                                  &displacementMystranObj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_makeDataSet(vertexSU2Obj[i], "Displacement", FieldIn, 3,
                                  &displacementSU2Obj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_initDataSet(displacementSU2Obj[i], 3, displacement0,
                                  &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_linkDataSet(pressureSU2Obj[i], Conserve,
                                  pressureMystranObj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_linkDataSet(displacementMystranObj[i], Interpolate,
                                  displacementSU2Obj[i], &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = caps_closeBound(boundObj[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    /* --------------------------------------------------------------- */

    // Link surface mesh from EGADS to TetGen
    status = caps_childByName(surfMeshObj, VALUE, ANALYSISOUT, "Surface_Mesh",
                              &source, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("surfMeshObj childByName for Surface_Mesh = %d\n", status);
      goto cleanup;
    }
    status = caps_childByName(meshObj, VALUE, ANALYSISIN, "Surface_Mesh",
                              &target, &nErr, &errors);
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

    /* Link the volume mesh from TetGen to SU2 */
    status = caps_childByName(meshObj, VALUE, ANALYSISOUT, "Volume_Mesh",
                              &source, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) {
      printf("meshObj childByName for Volume_Mesh = %d\n", status);
      goto cleanup;
    }
    status = caps_childByName(su2Obj, VALUE, ANALYSISIN, "Mesh",  &target,
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

    /* --------------------------------------------------------------- */

    // Set parameters for SU2
    status = setValueByName_String(su2Obj, ANALYSISIN, "Proj_Name", 1, 1,
                                   projectName);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = setValueByName_String(su2Obj, ANALYSISIN, "SU2_Version", 1, 1,
                                   "Blackbird");
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal = refVelocity/speedofSound;
    status = setValueByName_Double(su2Obj, ANALYSISIN, "Mach", 1, 1, &doubleVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = setValueByName_String(su2Obj, ANALYSISIN, "Equation_Type", 1, 1,
                                   "compressible");
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 5;
    status = setValueByName_Integer(su2Obj, ANALYSISIN, "Num_Iter", 1, 1,
                                    &intVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = setValueByName_String(su2Obj, ANALYSISIN, "Output_Format", 1, 1,
                                   "Tecplot");
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = true;
    status = setValueByName_Boolean(su2Obj, ANALYSISIN, "Overwrite_CFG", 1, 1,
                                    &boolVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = 0.5*refDensity*refVelocity*refVelocity;
    status = setValueByName_Double(su2Obj, ANALYSISIN, "Pressure_Scale_Factor", 1, 1,
                                   &doubleVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Find & set Boundary_Conditions for SU2
    su2BC = (capsTuple *) EG_alloc(numSU2BC*sizeof(capsTuple));
    su2BC[0].name  = EG_strdup("Skin");
    su2BC[0].value = EG_strdup("{\"bcType\": \"Inviscid\"}");

    su2BC[1].name  = EG_strdup("SymmPlane");
    su2BC[1].value =  EG_strdup("SymmetryY");

    su2BC[2].name  = EG_strdup("Farfield");
    su2BC[2].value = EG_strdup("farfield");

    status = setValueByName_Tuple(su2Obj, ANALYSISIN, "Boundary_Condition",
                                  numSU2BC, 1, su2BC);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* --------------------------------------------------------------- */

    // Set Mystran inputs - Materials
    material = (capsTuple *) EG_alloc(numMaterial*sizeof(capsTuple));
    material[0].name  = EG_strdup("Madeupium");
    material[0].value = EG_strdup("{\"youngModulus\": 72.0E9, \"density\": 2.8E3}");

    status = setValueByName_Tuple(mystranObj, ANALYSISIN, "Material",
                                  numMaterial, 1, material);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                    - Properties
    property = (capsTuple *) EG_alloc(numProperty*sizeof(capsTuple));
    property[0].name  = EG_strdup("Skin");
    property[0].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.05}");
    property[1].name  = EG_strdup("Rib_Root");
    property[1].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.1}");

    status = setValueByName_Tuple(mystranObj, ANALYSISIN, "Property",
                                  numProperty, 1, property);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                   - Constraints
    constraint = (capsTuple *) EG_alloc(numConstraint*sizeof(capsTuple));
    constraint[0].name  = EG_strdup("edgeConstraint");
    constraint[0].value = EG_strdup("{\"groupName\": \"Rib_Root\", \"dofConstraint\": 123456}");

    status = setValueByName_Tuple(mystranObj, ANALYSISIN, "Constraint",
                                  numConstraint, 1, constraint);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = setValueByName_String(mystranObj, ANALYSISIN, "Proj_Name", 1, 1,
                                   projectName);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 3;
    status = setValueByName_Integer(mystranObj, ANALYSISIN, "Edge_Point_Min", 1, 1,
                                    &intVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 10;
    status = setValueByName_Integer(mystranObj, ANALYSISIN, "Edge_Point_Max", 1, 1,
                                    &intVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    tessParams[0] = 1.5;
    tessParams[1] = 0.1;
    tessParams[2] = 0.15;
    status = setValueByName_Double(mystranObj, ANALYSISIN, "Tess_Params", 3, 1,
                                   tessParams);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = setValueByName_String(mystranObj, ANALYSISIN, "Analysis_Type", 1, 1,
                                   "Static");
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = (int) false;
    status = setValueByName_Boolean(mystranObj, ANALYSISIN, "Quad_Mesh", 1, 1,
                                    &boolVal);
    if (status != CAPS_SUCCESS) goto cleanup;

    load = (capsTuple *) EG_alloc(numLoad*sizeof(capsTuple));
    load[0].name  = EG_strdup("pressureAero");
    load[0].value = EG_strdup("{\"loadType\": \"PressureExternal\", \"loadScaleFactor\": -1.0}");

    status = setValueByName_Tuple(mystranObj, ANALYSISIN, "Load", numLoad, 1,
                                  load);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* --------------------------------------------------------------- */

    // Do the analysis -- actually run EGADS
    status = caps_preAnalysis(surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Everything is done in preAnalysis, so we just do the post
    status = caps_postAnalysis(surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* --------------------------------------------------------------- */

    // Do the analysis -- actually run TetGen
    status = caps_preAnalysis(meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Everything is done in preAnalysis, so we just do the post
    status = caps_postAnalysis(meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Start iteration loop
    for (iter = 0; iter < numIteration; iter++) {

        /* --------------------------------------------------------------- */

#ifdef DEBUG
        /* Dump out vertex sets for debugging */
        for (i = 0; i < 3; i++) {
            sprintf(filename, "SU2_%d.vs", i);
            status = caps_outputVertexSet( vertexSU2Obj[i], filename );
            if (status != CAPS_SUCCESS) goto cleanup;
        }
#endif
        /* --------------------------------------------------------------- */

        // Do the pre-analysis  for SU2
        status = caps_preAnalysis(su2Obj, &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Execute SU2 via system call
        (void) getcwd(currentPath, PATH_MAX);

        /* get analysis info */
        status = caps_analysisInfo(su2Obj, &analysisPath, &unitSystem, &major,
                                   &minor, &intents, &nFields, &fnames, &ranks,
                                   &fInOut, &exec, &dirty);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (chdir(analysisPath) != 0) {
            printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
            status = CAPS_DIRERR;
            goto cleanup;
        }

        if (iter > 0) {
          printf("\n\nRunning SU2_DEF!\n\n");
          system("SU2_DEF aeroelasticSimple_Iterative.cfg > su2DEFOut.txt");
        }

        printf("\n\nRunning SU2_CFD!\n\n");
        system("SU2_CFD aeroelasticSimple_Iterative.cfg > su2CFDOut.txt");

        (void) chdir(currentPath);

        // Run su2 post - to kick off data transfer
        status = caps_postAnalysis(su2Obj, &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        /* --------------------------------------------------------------- */

#ifdef DEBUG
        /* Dump out vertex sets for debugging */
        for (i = 0; i < 3; i++) {
            sprintf(filename, "mystran_%d.vs", i);
            status = caps_outputVertexSet( vertexMystranObj[i], filename );
            if (status != CAPS_SUCCESS) goto cleanup;
        }
#endif

        /* --------------------------------------------------------------- */

        // Do the analysis  for Mystran
        status = caps_preAnalysis(mystranObj, &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Execute Mystran via system call
        (void) getcwd(currentPath, PATH_MAX);

        /* get analysis info */
        status = caps_analysisInfo(mystranObj, &analysisPath, &unitSystem,
                                   &major, &minor, &intents, &nFields, &fnames,
                                   &ranks, &fInOut, &exec, &dirty);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (chdir(analysisPath) != 0) {
            printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
            status = CAPS_DIRERR;
            goto cleanup;
        }

        printf("\n\nRunning mystran!\n\n");
        system("mystran.exe aeroelasticSimple_Iterative.dat > mystranOut.txt");

        (void) chdir(currentPath);

        // Run Mystran post
        status = caps_postAnalysis(mystranObj, &nErr, &errors);
        if (nErr != 0) printErrors(nErr, errors);
        if (status != CAPS_SUCCESS) goto cleanup;

        /* --------------------------------------------------------------- */
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d\n",
                                       status);

    freeTuple(numSU2BC, su2BC);
    freeTuple(numConstraint, constraint);
    freeTuple(numLoad, load);
    freeTuple(numProperty, property);
    freeTuple(numMaterial, material);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
