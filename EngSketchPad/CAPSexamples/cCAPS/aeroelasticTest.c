/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             fun3d, tetgen, mystran AIM tester
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


#define DEPENDENT


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
    capsObj  problemObj, surfMeshObj, meshObj, fun3dObj, fun3dObj2, mystranObj, tempObj;
    capsErrs *errors;

    // Data transfer
    capsObj topBoundObj, vertexSourceObj, vertexDestObj, dataSourceObj, dataDestObj;

    // CAPS return values
    int  nErr;
    capsObj source, target;

    int major, minor, nFields, *ranks, *fInOut, dirty, exec;
    char *analysisPath=NULL, *unitSystem, *intents, **fnames;

    // Input values
    capsTuple         *fun3dBC = NULL, *material = NULL, *constraint = NULL, *property = NULL, *load = NULL;
    int               numFUN3DBC = 3, numMaterial = 1, numConstraint = 1, numProperty = 2, numLoad = 1;

    int               intVal;
    double            doubleVal, tessParams[3];
    enum capsBoolean  boolVal;
    char *stringVal = NULL;

    const char projectName[] = "aeroelasticSimple";

    char currentPath[PATH_MAX];

    double speedofSound = 340.0, refVelocity = 100.0, refDensity = 1.2;

    int runAnalysis = (int) true;

    printf("\n\nAttention: aeroelasticTest is hard coded to look for ../csmData/aeroelasticDataTransferSimple.csm\n");

    if (argc > 2) {
        printf(" usage: aeroelasticTest outLevel\n");
        return 1;
    } else if (argc == 2) {
        outLevel = atoi(argv[1]);
    }

    status = caps_open("FUN3D_MyStran_Aeroelastic_Example", NULL, 0,
                       "../csmData/aeroelasticDataTransferSimple.csm", outLevel,
                       &problemObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Load the AIMs
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "egadsTessAIM", NULL,
                               NULL, NULL, &exec, &surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "tetgenAIM", NULL,
                               NULL, NULL, &exec, &meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS)  goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "fun3dAIM", NULL,
                               NULL, NULL, &exec, &fun3dObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    exec   = 0;
    status = caps_makeAnalysis(problemObj, "mystranAIM", NULL,
                               NULL, NULL, &exec, &mystranObj, &nErr, &errors);
   if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    /* unused analysis object */
    exec   = 0;
    status = caps_makeAnalysis(problemObj, "fun3dAIM", "DummyName",
                               NULL, NULL, &exec, &fun3dObj2, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Create data bounds
    status = caps_makeBound(problemObj, 2, "Skin_Top", &topBoundObj);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_makeVertexSet(topBoundObj, fun3dObj, NULL, &vertexSourceObj,
                                &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_makeVertexSet(topBoundObj, mystranObj, NULL, &vertexDestObj,
                                &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_makeDataSet(vertexSourceObj, "Pressure", FieldOut, 0,
                              &dataSourceObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_makeDataSet(vertexDestObj, "Pressure", FieldIn, 0,
                              &dataDestObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_closeBound(topBoundObj);
    if (status != CAPS_SUCCESS) goto cleanup;


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

    /* Link the volume mesh from TetGen to Fun3D */
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


    // Find & set Boundary_Conditions for FUN3D
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Boundary_Condition", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    fun3dBC = (capsTuple *) EG_alloc(numFUN3DBC*sizeof(capsTuple));
    fun3dBC[0].name = EG_strdup("Skin");
    fun3dBC[0].value = EG_strdup("{\"bcType\": \"Inviscid\", \"wallTemperature\": 1}");

    fun3dBC[1].name = EG_strdup("SymmPlane");
    fun3dBC[1].value =  EG_strdup("SymmetryY");

    fun3dBC[2].name = EG_strdup("Farfield");
    fun3dBC[2].value = EG_strdup("farfield");

    status = caps_setValue(tempObj, Tuple, numFUN3DBC, 1,  (void **) fun3dBC,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    // Set additional parameters for FUN3D
    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Mach", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = refVelocity/speedofSound;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) &doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Num_Iter", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal  = 10;
    status = caps_setValue(tempObj, Integer, 1, 1, (void *) &intVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Viscous", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    stringVal = EG_strdup("inviscid");
    status = caps_setValue(tempObj, String, 1, 1, (void *) stringVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    EG_free(stringVal); stringVal = NULL;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Restart_Read", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    stringVal = EG_strdup("off");
    status = caps_setValue(tempObj, String, 1, 1, (void *) stringVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;
    EG_free(stringVal); stringVal = NULL;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Overwrite_NML", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = true;
    status = caps_setValue(tempObj, Boolean, 1, 1, (void *) &boolVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Proj_Name", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_setValue(tempObj, String, 1, 1, (void *) &projectName,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(fun3dObj, VALUE, ANALYSISIN, "Pressure_Scale_Factor", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    doubleVal  = 0.5 * refDensity * refVelocity*refVelocity;
    status = caps_setValue(tempObj, Double, 1, 1, (void *) &doubleVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    // Set Mystran inputs - Materials
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Material", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    material = (capsTuple *) EG_alloc(numMaterial*sizeof(capsTuple));
    material[0].name = EG_strdup("Madeupium");
    material[0].value = EG_strdup("{\"youngModulus\": 72.0E9, \"density\": 2.8E3}");

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
    property[0].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.05}");
    property[1].name  = EG_strdup("Rib_Root");
    property[1].value = EG_strdup("{\"propertyType\": \"Shell\", \"membraneThickness\": 0.1}");

    status = caps_setValue(tempObj, Tuple, numProperty, 1,  (void **) property,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    //                      - Constraints
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Constraint", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    constraint = (capsTuple *) EG_alloc(numConstraint*sizeof(capsTuple));
    constraint[0].name = EG_strdup("edgeConstraint");
    constraint[0].value = EG_strdup("{\"groupName\": \"Rib_Root\", \"dofConstraint\": 123456}");

    status = caps_setValue(tempObj, Tuple, numConstraint, 1,
                           (void **) constraint, NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Proj_Name", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_setValue(tempObj, String, 1, 1, (void *) &projectName,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Edge_Point_Min", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 3;
    status = caps_setValue(tempObj, Integer, 1, 1, (void *) &intVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Edge_Point_Max", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    intVal = 10;
    status = caps_setValue(tempObj, Integer, 1, 1, (void *) &intVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Tess_Params", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    tessParams[0] = 0.5;
    tessParams[1] = 0.1;
    tessParams[2] = 0.15;
    status = caps_setValue(tempObj, Double, 3, 1, (void **) tessParams,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;


    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Analysis_Type", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    stringVal = EG_strdup("Static");
    status = caps_setValue(tempObj, String, 1, 1, (void *) stringVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Quad_Mesh", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    boolVal = (int) false;
    status = caps_setValue(tempObj, Boolean, 1, 1, (void *) &boolVal,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

#ifdef DEPENDENT
    // look at dependent Objects
    {
        int     nObjs;
        capsObj *objs;

        status = caps_dirtyAnalysis(problemObj, &nObjs, &objs);
        printf("  Problem dependencies = %d  (status = %d)\n", nObjs, status);
        if (objs != NULL) EG_free(objs);
        if (nObjs != 0) {
            status = caps_dirtyAnalysis(topBoundObj, &nObjs, &objs);
            printf("  Bound   dependencies = %d  (status = %d)\n", nObjs, status);
            if (objs != NULL) EG_free(objs);
        }
    }
#endif

    // =========================================================
    // Do the analysis -- actually run EGADS
    status = caps_preAnalysis(surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Everything is done in preAnalysis, so we just do the post
    status = caps_postAnalysis(surfMeshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;
    // =========================================================

    // =========================================================
    // Do the analysis -- actually run TetGen
    status = caps_preAnalysis(meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Everything is done in preAnalysis, so we just do the post
    status = caps_postAnalysis(meshObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;
    // =========================================================

    // Do the analysis  for FUN3D
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

    if (runAnalysis == (int) false) {
        printf("\n\nNot Running fun3d!\n\n");
    } else {
        printf("\n\nRunning fun3d!\n\n");
        system("nodet_mpi --write_aero_loads_to_file > fun3dOutput.txt");
    }

    (void) chdir(currentPath);

    // Run fun3d post - to kick off data transfer
    status = caps_postAnalysis(fun3dObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

#ifdef DEPENDENT
    // Look at dependent objects
    {
        int     nObjs;
        capsObj *objs;

        status = caps_dirtyAnalysis(problemObj, &nObjs, &objs);
        printf("  Problem dependencies = %d  (status = %d)\n", nObjs, status);
        if (objs != NULL) EG_free(objs);
        if (nObjs != 0) {
            status = caps_dirtyAnalysis(topBoundObj, &nObjs, &objs);
            printf("  Bound   dependencies = %d  (status = %d)\n", nObjs, status);
            if (objs != NULL) EG_free(objs);
        }
    }
#endif

    // Run mystran pre- and post- to get things enabled for data transfer
    status = caps_preAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = caps_postAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

#ifdef DEPENDENT
    // look at dependent Objects
    {
        int     nObjs;
        capsObj *objs;

        status = caps_dirtyAnalysis(problemObj, &nObjs, &objs);
        printf("  Problem dependencies = %d  (status = %d)\n", nObjs, status);
        if (objs != NULL) EG_free(objs);
        if (nObjs != 0) {
            status = caps_dirtyAnalysis(topBoundObj, &nObjs, &objs);
            printf("  Bound   dependencies = %d  (status = %d)\n", nObjs, status);
            if (objs != NULL) EG_free(objs);
        }
    }
#endif

    // Set load tuple for Mystran
    status = caps_childByName(mystranObj, VALUE, ANALYSISIN, "Load", &tempObj,
                              &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    load = (capsTuple *) EG_alloc(numLoad*sizeof(capsTuple));
    load[0].name  = EG_strdup("pressureAero");
    load[0].value = EG_strdup("{\"loadType\": \"PressureExternal\", \"loadScaleFactor\": -1.0}");

    status = caps_setValue(tempObj, Tuple, numLoad, 1,  (void **) load,
                           NULL, NULL, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Do the analysis  for Mystran
    status = caps_preAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Execute Mystran via system call
    (void) getcwd(currentPath, PATH_MAX);

    /* get analysis info */
    status = caps_analysisInfo(mystranObj, &analysisPath, &unitSystem, &major,
                               &minor, &intents, &nFields, &fnames, &ranks,
                               &fInOut, &exec, &dirty);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (chdir(analysisPath) != 0) {
        printf(" ERROR: Cannot change directory to -> %s\n", analysisPath);
        status = CAPS_DIRERR;
        goto cleanup;
    }

    if (runAnalysis == (int) false) {
        printf("\n\nNOT Running mystran!\n\n");
    } else {
        printf("\n\nRunning mystran!\n\n");
        system("mystran.exe aeroelasticSimple.dat > mystranOutput.txt");
    }

    if (chdir(currentPath) != 0) {
        printf(" ERROR: Cannot change directory to -> %s\n", currentPath);
        status = CAPS_DIRERR;
        goto cleanup;
    }

    // Run Mystran post
    status = caps_postAnalysis(mystranObj, &nErr, &errors);
    if (nErr != 0) printErrors(nErr, errors);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("\n\nPremature exit - status = %d\n",
                                       status);

    freeTuple(numFUN3DBC, fun3dBC);
    freeTuple(numConstraint, constraint);
    freeTuple(numLoad, load);
    freeTuple(numProperty, property);
    freeTuple(numMaterial, material);

    EG_free(stringVal);

    i = 0;
    if (status == CAPS_SUCCESS) i = 1;
    (void) caps_close(problemObj, i, NULL);

    return status;
}
