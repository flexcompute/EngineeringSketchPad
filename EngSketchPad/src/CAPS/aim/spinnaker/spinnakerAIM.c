/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Spinnaker AIM
 *
 *     Written by Dr. Marshall Galbraith MIT
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 */


/*!\mainpage Introduction
 *
 * \section overviewSPINNAKER Spinnaker AIM Overview
 * This module can be used to interface with the open-source Spinnaker code developed at Sandia National Laboratories
 *  with geometry in the CAPS system. For Spinnaker capabilities and related documentation, please refer to
 * https://github.com/spinnakerengine. Spinnaker expects a volume, surface, or area mesh and a corresponding
 * configuration file to perform the analysis.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsSPINNAKER and \ref aimOutputsSPINNAKER, respectively.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferSPINNAKER
 *
 * \subsection meshSPINNAKER Automatic generation of Spinnaker Exodus Mesh file
 * The mesh file from Spinnaker AIM is written in native Exodus
 * format ("filename.exo"). The description of the native Exodus mesh can be
 * found Exodus website (https://sandialabs.github.io/seacas-docs/html/index.html).
 * For the automatic generation of mesh file, Spinnaker AIM
 * depends on Mesh AIMs, for example, TetGen or AFLR4/3 AIM.
 *
 * \section examplesSPINNAKER Spinnaker Examples
 *  Here is an example that illustrated use of Spinnaker AIM \ref spinnakerExample. Note
 *  this AIM uses TetGen AIM for volume mesh generation.
 *
 */

#include <string.h>
#include <ctype.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"

#include "exodusWriter.h"
#include "libMeshb/sources/libmeshb7.h"

#include <exodusII.h>

/* these values requested to emulate Feflo.a behavior */
#define EXPORT_MESHB_VERTEX_ID (1)
#define EXPORT_MESHB_2D_ID (1)
#define EXPORT_MESHB_3D_ID (0)
#define EXPORT_MESHB_VERTEX_3 (10000000)
#define EXPORT_MESHB_VERTEX_4 (200000000)

#define MAX(A,B)         (((A) < (B)) ? (B) : (A))

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#else
#include <unistd.h>
#endif

#define MXCHAR  255

//static const char resultsFile[] = "results.exo";

//#define DEBUG


enum aimInputs
{
  Proj_Name = 1,        /* index is 1-based */
  inSolutionFile,
  inOutputScalarField,
  inRestartFile,
  inMesh_Morph,
  inMesh,
  NUMINPUT = inMesh       /* Total number of inputs */
};

enum aimOutputs
{
  outScalarFieldSolbFile = 1,         /* index is 1-based */
  outRestartSolbFile,
  NUMOUTPUT = outRestartSolbFile  /* Total number of inputs */
};

typedef struct {

  char* filename;

  int numField;
  char **field;

} outputFieldStruct;

int initiate_outputFieldStruct(outputFieldStruct *outputField)
{
  outputField->filename = NULL;
  outputField->numField = 0;
  outputField->field = NULL;
  return CAPS_SUCCESS;
}


int destroy_outputFieldStruct(outputFieldStruct *outputField)
{
  int i;
  AIM_FREE(outputField->filename);
  for (i = 0; i < outputField->numField; i++) {
    AIM_FREE(outputField->field[i]);
  }
  AIM_FREE(outputField->field);
  outputField->numField = 0;
  return CAPS_SUCCESS;
}

typedef struct {

  // Mesh reference obtained from meshing AIM
  aimMeshRef *meshRef, meshRefObj;

} aimStorage;


static int initiate_aimStorage(aimStorage *spinnakerInstance)
{
  // Mesh holders
  spinnakerInstance->meshRef = NULL;
  aim_initMeshRef(&spinnakerInstance->meshRefObj, aimUnknownMeshType);

  return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *spinnakerInstance)
{
  spinnakerInstance->meshRef = NULL;
  aim_freeMeshRef(&spinnakerInstance->meshRefObj);

  return CAPS_SUCCESS;
}

/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@null@*/ /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS; // Function return
    int i;
    //int *ints=NULL;
    //char **strs=NULL;
    //const char *keyWord;
    //char *keyValue = NULL;
    //double real = 1;

    aimStorage *spinnakerInstance=NULL;

    #ifdef DEBUG
        printf("\n spinnakerAIM/aimInitialize   inst = %d!\n", inst);
    #endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames = 0;

    /* specify the name of each field variable * /
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
    strs[0]  = EG_strdup("Pressure");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
    *fnames  = strs;


    / * specify the dimension of each field variable * /
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = 1;
    ints[1]  = 1;
    ints[2]  = 1;
    ints[3]  = 1;
    ints[4]  = 3;
    *franks   = ints;
    ints = NULL;

    / * specify if a field is an input field or output field * /
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldOut;
    ints[1]  = FieldOut;
    ints[2]  = FieldOut;
    ints[3]  = FieldOut;
    ints[4]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;
    */

    // Allocate spinnakerInstance
    AIM_ALLOC(spinnakerInstance, 1, aimStorage, aimInfo, status);
    *instStore = spinnakerInstance;

    // Set initial values for spinnakerInstance
    initiate_aimStorage(spinnakerInstance);

cleanup:
    if (status != CAPS_SUCCESS) {
        /* release all possibly allocated memory on error */
        if (*fnames != NULL)
          for (i = 0; i < *nFields; i++) AIM_FREE((*fnames)[i]);
        AIM_FREE(*franks);
        AIM_FREE(*fInOut);
        AIM_FREE(*fnames);
        AIM_FREE(*instStore);
        *nFields = 0;
    }

    return status;
}


int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
              char **ainame, capsValue *defval)
{
    int status = CAPS_SUCCESS;
    //aimStorage *spinnakerInstance;

#ifdef DEBUG
    printf(" spinnakerAIM/aimInputs  index = %d!\n", index);
#endif

    *ainame = NULL;

    // Spinnaker Inputs
    /*! \page aimInputsSPINNAKER AIM Inputs
     */

    //spinnakerInstance = (aimStorage *) instStore;
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        AIM_STRDUP(defval->vals.string, "spinnaker_CAPS", aimInfo, status);

        /*! \page aimInputsSPINNAKER
         * - <B> Proj_Name = "spinnaker_CAPS"</B> <br>
         * This corresponds to the project name used for file naming.
         */

    } else if (index == inSolutionFile) {
        *ainame              = EG_strdup("SolutionFile");
        defval->type         = String;
        defval->nullVal      = IsNull;

       /*! \page aimInputsSPINNAKER
        * - <B>SolutionFile = NULL</B> <br>
        * Spinnaker exodus solution file for generating ScalarFieldSolbFile
        */

    } else if (index == inOutputScalarField) {
        *ainame              = EG_strdup("OutputScalarField");
        defval->type         = String;
        defval->nullVal      = NotNull;
        AIM_STRDUP(defval->vals.string, "velocity_0", aimInfo, status);

       /*! \page aimInputsSPINNAKER
        * - <B>OutputScalarField = "velocity_0"</B> <br>
        * Scalar field quantity for the ScalarFieldSolbFile output.
        */

    } else if (index == inRestartFile) {
        *ainame              = EG_strdup("RestartFile");
        defval->type         = String;
        defval->nullVal      = IsNull;

       /*! \page aimInputsSPINNAKER
        * - <B>RestartFile = NULL</B> <br>
        * Restart file for spinnnaker. A libMeshb file will be converted to exodus.
        */

    } else if (index == inMesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        defval->type         = Boolean;
        defval->lfixed       = Fixed;
        defval->vals.integer = (int) false;
        defval->dim          = Scalar;
        defval->nullVal      = NotNull;

        /*! \page aimInputsSPINNAKER
         * - <B> Mesh_Morph = False</B> <br>
         * Project previous surface mesh onto new geometry and write out a 'Proj_Name'_body#.dat file.
         */

    } else if (index == inMesh) {
        *ainame             = AIM_NAME(Mesh);
        defval->type        = PointerMesh;
        defval->nrow        = 1;
        defval->lfixed      = Fixed;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;
        AIM_STRDUP(defval->meshWriter, MESHWRITER, aimInfo, status);

        /*! \page aimInputsSPINNAKER
         * - <B>Mesh = NULL</B> <br>
         * An Area_Mesh or Volume_Mesh link for 2D and 3D calculations respectively.
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(/*@unused@*/ void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status

    // AIM input bodies
    int  numBody;
    ego *bodies = NULL;
    const char *intents;

    aimStorage *spinnakerInstance;

    spinnakerInstance = (aimStorage *) instStore;
    AIM_NOTNULL(spinnakerInstance, aimInfo, status);
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Free our meshRef
    (void) aim_freeMeshRef(&spinnakerInstance->meshRefObj);

    if (aimInputs[inMesh-1].nullVal == IsNull &&
        aimInputs[inMesh_Morph-1].vals.integer == (int) false) {
        AIM_ANALYSISIN_ERROR(aimInfo, inMesh, "'Mesh' input must be linked to a 'Volume_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    // Get mesh
    spinnakerInstance->meshRef = (aimMeshRef *)aimInputs[inMesh-1].vals.AIMptr;

    if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
        spinnakerInstance->meshRef == NULL) { // If we are mighty morphing

        // Lets "load" the meshRef now since it's not linked
        status = aim_loadMeshRef(aimInfo, &spinnakerInstance->meshRefObj);
        AIM_STATUS(aimInfo, status);

        // Mightly Morph the mesh
        status = aim_morphMeshUpdate(aimInfo, &spinnakerInstance->meshRefObj, numBody, bodies);
        AIM_STATUS(aimInfo, status);
        /*@-immediatetrans@*/
        spinnakerInstance->meshRef = &spinnakerInstance->meshRefObj;
        /*@+immediatetrans@*/
    }
    AIM_NOTNULL(spinnakerInstance->meshRef, aimInfo, status);

    status = CAPS_SUCCESS;
cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(/*@unused@*/ const void *instStore, void *aimInfo, capsValue *aimInputs)
{
  // Function return flag
  int status;

  int i, j, imap;
  size_t slen = 0;

  // EGADS return values
  const char   *intents;
  ego body;
  int state, nGlobal, ptype, pindex;
  double xyz[3];
  FILE *fp = NULL;

  // Output filename
  char filename[PATH_MAX];
  char filepath[PATH_MAX];

  // AIM input bodies
  int  numBody;
  ego *bodies = NULL;

  // exodus file I/O
  int CPU_word_size = sizeof(double);
  int IO_word_size = sizeof(double);
  float version;
  int exoid = 0;
  ex_init_params par;
  double *x=NULL, *y=NULL, *z=NULL;
  int dim, nVertex;
  int num_vars=0, ivar;
  char **var_names=NULL;
  double **nodals=NULL;

  // libMeshb file I/O
  int meshbVersion;
  int64_t solbID=0;
  int numType, solSize, typTab[GmfMaxTyp];
  int *szfld = NULL;
  void *BegTab[42];
  void *EndTab[42];

  // Volume Mesh obtained from meshing AIM
  aimMeshRef *meshRef = NULL;

  const aimStorage *spinnakerInstance;

  spinnakerInstance = (const aimStorage *) instStore;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

#ifdef DEBUG
  printf(" spinnakerAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif
  if ((numBody <= 0) || (bodies == NULL)) {
    AIM_ERROR(aimInfo, "No Bodies!");
    return CAPS_SOURCEERR;
  }

  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get mesh
  meshRef = spinnakerInstance->meshRef;
  AIM_NOTNULL(meshRef, aimInfo, status);

  if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true) { // If we are mighty morphing
    if (aimInputs[inMesh-1].nullVal == NotNull) {
      // store the current mesh for future iterations
      status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[inMesh-1].vals.AIMptr, MESHEXTENSION);
      AIM_STATUS(aimInfo, status);

    } else {

      snprintf(filename, PATH_MAX, "%s%s", meshRef->fileName, MESHEXTENSION);

      exoid = ex_open(filename, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        return CAPS_IOERR;
      }

      status = ex_get_init_ext(exoid, &par);
      AIM_STATUS(aimInfo, status);

      dim     = par.num_dim;
      nVertex = par.num_nodes;

      AIM_ALLOC(x, nVertex, double, aimInfo, status);
      AIM_ALLOC(y, nVertex, double, aimInfo, status);
      if (dim == 3)
        AIM_ALLOC(z, nVertex, double, aimInfo, status);

      /* get all of the vertices */
/*@-nullpass@*/
      status = ex_get_coord(exoid, x, y, z);
      AIM_STATUS(aimInfo, status);
/*@+nullpass@*/

      ex_close(exoid);
      exoid = 0;

      for (imap = 0; imap < meshRef->nmap; imap++) {
        status = EG_statusTessBody(meshRef->maps[imap].tess, &body, &state, &nGlobal);
        AIM_STATUS(aimInfo, status);

        // Write the map
        for (i = 0; i < nGlobal; i++) {
          status = EG_getGlobal(meshRef->maps[imap].tess, i+1, &ptype, &pindex, xyz);
          AIM_STATUS(aimInfo, status);

          j = meshRef->maps[imap].map[i]-1;

          x[j] = xyz[0];
          y[j] = xyz[1];
          if (z != NULL)
            z[j] = xyz[2];
        }
      }

      exoid = ex_open(filename, EX_WRITE | EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        status = CAPS_IOERR;
        goto cleanup;
      }

      /* set all of the vertices */
/*@-nullpass@*/
      status = ex_put_coord(exoid, x, y, z);
      AIM_STATUS(aimInfo, status);
/*@+nullpass@*/

      AIM_FREE(x);
      AIM_FREE(y);
      AIM_FREE(z);

      ex_close(exoid);
      exoid = 0;
    }
  }

  /* create a symbolic link to the file name*/
  snprintf(filename, PATH_MAX, "%s%s", meshRef->fileName, MESHEXTENSION);
  snprintf(filepath, PATH_MAX, "%s%s", aimInputs[Proj_Name-1].vals.string, MESHEXTENSION);
  status = aim_symLink(aimInfo, filename, filepath);
  AIM_STATUS(aimInfo, status);

  if ( aimInputs[inRestartFile-1].nullVal != IsNull) {

    if (strlen(aimInputs[inRestartFile-1].vals.string) < 6) {
      AIM_ERROR(aimInfo, "Filename %s is too short", aimInputs[inRestartFile-1].vals.string);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    slen = strlen(aimInputs[inRestartFile-1].vals.string);

    if (strcmp(aimInputs[inRestartFile-1].vals.string+slen-4,".sol") == 0 ||
        strcmp(aimInputs[inRestartFile-1].vals.string+slen-5,".solb") == 0) {

      if (aimInputs[inSolutionFile-1].nullVal == IsNull) {
        AIM_ERROR(aimInfo, "SolutionFile must be specified!");
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      status = aim_file(aimInfo, aimInputs[inSolutionFile-1].vals.string, filename);
      AIM_STATUS(aimInfo, status);

      exoid = ex_open(filename, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        return CAPS_IOERR;
      }

      status = ex_get_init_ext(exoid, &par);
      AIM_STATUS(aimInfo, status);

      status = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
      AIM_STATUS(aimInfo, status);

      AIM_ALLOC(var_names, num_vars, char*, aimInfo, status);

      for (ivar = 0; ivar < num_vars; ivar++) var_names[ivar] = NULL;
      for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(var_names[ivar], EX_MAX_NAME, char, aimInfo, status);

      status = ex_get_variable_names(exoid, EX_NODAL, num_vars, (char**)var_names);
      AIM_STATUS(aimInfo, status);

      ex_close(exoid);
      exoid = 0;


      //ex_opts(EX_VERBOSE | EX_DEBUG | EX_NULLVERBOSE);

      snprintf(filepath, PATH_MAX, "%s%s", aimInputs[Proj_Name-1].vals.string, MESHEXTENSION);

      /* Open up the exodus restart file */
      status = aim_file(aimInfo, filepath, filename);
      AIM_STATUS(aimInfo, status);

      exoid = ex_open(filename, EX_WRITE | EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        return CAPS_IOERR;
      }

      status = ex_put_variable_param(exoid, EX_NODAL, num_vars);
      AIM_STATUS(aimInfo, status);

      status = ex_put_variable_names(exoid, EX_NODAL, num_vars, (char**)var_names);
      AIM_STATUS(aimInfo, status);


      /* Read meshb file */
      solbID = GmfOpenMesh(aimInputs[inRestartFile-1].vals.string, GmfRead, &meshbVersion, &dim);
      if (solbID == 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", aimInputs[inRestartFile-1].vals.string);
        status = CAPS_IOERR;
        goto cleanup;
      }

      status = GmfStatKwd(solbID, GmfSolAtVertices, &numType, &solSize, typTab);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      status = GmfGotoKwd(solbID, GmfSolAtVertices);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      if (solSize != num_vars) {
        AIM_ERROR(aimInfo, "Expecting %d vars in %s, found %d", num_vars, aimInputs[inRestartFile-1].vals.string, solSize);
        status = CAPS_IOERR;
        goto cleanup;
      }

      status = ex_get_init_ext(exoid, &par);
      AIM_STATUS(aimInfo, status);

      AIM_ALLOC(nodals, par.num_nodes, double*, aimInfo, status);
      for (ivar = 0; ivar < num_vars; ivar++) nodals[ivar] = NULL;
      for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(nodals[ivar], par.num_nodes, double, aimInfo, status);

      AIM_ALLOC(szfld, num_vars, int, aimInfo, status);

      for (ivar = 0; ivar < num_vars; ivar++) {
        szfld[ivar] = 1;

        BegTab[ivar] = &nodals[ivar][0];
        EndTab[ivar] = &nodals[ivar][par.num_nodes-1];
      }
      /*@-nullpass@*/
      status = GmfGetBlock(solbID, GmfSolAtVertices, 1, par.num_nodes, 0, NULL, NULL,
                           GmfArgTab, typTab, szfld, BegTab, EndTab);
      /*@+nullpass@*/
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      for (ivar = 0; ivar < num_vars; ivar++) {
        status = ex_put_var(exoid, 1, EX_NODAL, ivar+1, 1, par.num_nodes, nodals[ivar]);
        AIM_STATUS(aimInfo, status);
      }
    }
  }

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(x);
  AIM_FREE(y);
  AIM_FREE(z);

  if (var_names != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(var_names[ivar]);
  AIM_FREE(var_names);
  if (nodals != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(nodals[ivar]);
  AIM_FREE(nodals);
  AIM_FREE(szfld);

  if (fp != NULL) fclose(fp);
  if (exoid > 0) ex_close(exoid);

  return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *aimInputs)
{
  return CAPS_SUCCESS;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
  int status = CAPS_SUCCESS;

 /*! \page aimOutputsSPINNAKER AIM Outputs
   * Spinnaker outputs
   */

  #ifdef DEBUG
      printf(" spinnakerAIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
  #endif

  if (index == outScalarFieldSolbFile) {

    *aoname           = AIM_NAME(ScalarFieldSolbFile);
    form->type        = String;
    form->nullVal     = IsNull;

    /*! \page aimOutputsSPINNAKER
     * - <B> ScalarFieldSolbFile </B> <br>
     * String to file containing a scalar field in libMeshb solb format
     */

  } else if (index == outRestartSolbFile) {

    *aoname           = AIM_NAME(RestartSolbFile);
    form->type        = String;
    form->nullVal     = IsNull;

    /*! \page aimOutputsSPINNAKER
     * - <B> RestartSolbFile </B> <br>
     * String to file containing a scalar field
     */



  } else {
      status = CAPS_BADINDEX;
      AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);
  }

  AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
  return status;
}


static int
_getExodusNodalVar(void *aimInfo, int exoid, int num_vars, char *const* var_names, const char *var_name, int num_nodes, double *nodal)
{
  int status = CAPS_SUCCESS;
  int64_t time_step = ex_inquire_int(exoid, EX_INQ_TIME);

  int ivar = 0;

  while (ivar < num_vars && strcasecmp(var_names[ivar], var_name)) ivar++;
  if (ivar == num_vars) {
    AIM_ERROR(aimInfo, "Could not find '%s' in Spinnaker output file!", var_name);
    status = CAPS_IOERR;
    goto cleanup;
  }

  status = ex_get_var(exoid, time_step, EX_NODAL, ivar+1, 1, num_nodes, nodal);
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


// Calculate Spinnaker output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, capsValue *val)
{
  int status = CAPS_SUCCESS;

  int i, ivar;

  char spinFilename[PATH_MAX];
  char solFilename[PATH_MAX];
  char tmpFilename[PATH_MAX];

  //const char *analysisType = NULL;

  int CPU_word_size = sizeof(double);
  int IO_word_size = sizeof(double);
  float version;
  int exoid = 0;
  ex_init_params par;
  int64_t time_step = 0;


  int num_vars=0;
  char **var_names=NULL;
  //int whole_time_step = 1;
  double *nodal = NULL;
  double **nodals = NULL;
  int *node_id_map = NULL;
  double *nodal_ordered = NULL;
  double **nodals_ordered = NULL;
//  int num_elem_in_blk, num_nodes_per_elem, num_attr, iblk;
//  char elem_type[MAX_STR_LENGTH+1];
//  char name[2*MAX_STR_LENGTH+1];

  capsValue *SolutionFile = NULL;
  capsValue *OutputScalarField = NULL;
  FILE *fp = NULL;

  int meshbVersion;
  int64_t solbID=0;
  int *szfld = NULL;
  int *TypTab = NULL;
  void *BegTab[42];
  void *EndTab[42];

  // Mesh reference obtained from meshing AIM
  //aimMeshRef *meshRef = NULL;

  //aimStorage *spinnakerInstance = (aimStorage*)instStore;

#ifdef DEBUG
  printf(" spinnakerAIM/aimCalcOutput  index = %d!\n", index);
#endif

  // Get mesh
//  meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;
//  AIM_NOTNULL(meshRef, aimInfo, status);

  if (index == outScalarFieldSolbFile) {

    status = aim_getValue(aimInfo, inSolutionFile, ANALYSISIN, &SolutionFile);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(SolutionFile, aimInfo, status);

    status = aim_file(aimInfo, SolutionFile->vals.string, spinFilename);
    AIM_STATUS(aimInfo, status);

    exoid = ex_open(spinFilename, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                    &CPU_word_size, &IO_word_size, &version);
    if (exoid <= 0) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", spinFilename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    //ex_opts(EX_VERBOSE | EX_DEBUG | EX_NULLVERBOSE);

    status = ex_get_init_ext(exoid, &par);
    AIM_STATUS(aimInfo, status);

    status = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(var_names, num_vars, char*, aimInfo, status);
    AIM_ALLOC(nodal        , par.num_nodes, double, aimInfo, status);
    AIM_ALLOC(nodal_ordered, par.num_nodes, double, aimInfo, status);
    AIM_ALLOC(node_id_map  , par.num_nodes, int   , aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) var_names[ivar] = NULL;
    for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(var_names[ivar], EX_MAX_NAME, char, aimInfo, status);

    status = ex_get_variable_names(exoid, EX_NODAL, num_vars, (char**)var_names);
    AIM_STATUS(aimInfo, status);

    status = aim_getValue(aimInfo, inOutputScalarField, ANALYSISIN, &OutputScalarField);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(OutputScalarField, aimInfo, status);

    _getExodusNodalVar(aimInfo, exoid, num_vars, var_names, OutputScalarField->vals.string, par.num_nodes, nodal);

    /* get node mapping to undo partitioning reordering */
    status = ex_get_id_map(exoid, EX_NODE_MAP, node_id_map);
    AIM_STATUS(aimInfo, status);

    ex_close(exoid);
    exoid = 0;

    /* properly order the scalar value */
    for (i = 0; i < par.num_nodes; i++) {
      nodal_ordered[node_id_map[i]-1] = nodal[i];
    }

    /* Write binary meshb file */
    snprintf(tmpFilename, PATH_MAX, "%s.solb", OutputScalarField->vals.string);
    status = aim_file(aimInfo, tmpFilename, solFilename);
    AIM_STATUS(aimInfo, status);

    meshbVersion = 2;
    if (EXPORT_MESHB_VERTEX_3 < par.num_nodes) meshbVersion = 3;
    if (EXPORT_MESHB_VERTEX_4 < par.num_nodes) meshbVersion = 4;

    solbID = GmfOpenMesh(solFilename, GmfWrite, meshbVersion, par.num_dim);
    if (solbID == 0) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", solFilename);
      return CAPS_IOERR;
    }

    AIM_ALLOC(szfld, 1, int, aimInfo, status);
    szfld[0] = GmfSca;
    status = GmfSetKwd(solbID, GmfSolAtVertices, par.num_nodes, 1, szfld);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    /*@-nullpass@*/
    status = GmfSetBlock(solbID, GmfSolAtVertices, 1, par.num_nodes, 0, NULL, NULL,
                         GmfDouble, &nodal_ordered[0], &nodal_ordered[par.num_nodes-1]);
    /*@+nullpass@*/
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    AIM_STRDUP(val->vals.string, solFilename, aimInfo, status);

  } else if (index == outRestartSolbFile) {

    status = aim_getValue(aimInfo, inSolutionFile, ANALYSISIN, &SolutionFile);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(SolutionFile, aimInfo, status);

    status = aim_file(aimInfo, SolutionFile->vals.string, spinFilename);
    AIM_STATUS(aimInfo, status);

    exoid = ex_open(spinFilename, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                    &CPU_word_size, &IO_word_size, &version);
    if (exoid <= 0) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", spinFilename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    //ex_opts(EX_VERBOSE | EX_DEBUG | EX_NULLVERBOSE);

    status = ex_get_init_ext(exoid, &par);
    AIM_STATUS(aimInfo, status);

    status = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(var_names, num_vars, char*, aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) var_names[ivar] = NULL;
    for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(var_names[ivar], EX_MAX_NAME, char, aimInfo, status);

    status = ex_get_variable_names(exoid, EX_NODAL, num_vars, (char**)var_names);
    AIM_STATUS(aimInfo, status);


    AIM_ALLOC(nodals        , par.num_nodes, double*, aimInfo, status);
    AIM_ALLOC(nodals_ordered, par.num_nodes, double*, aimInfo, status);
    for (ivar = 0; ivar < num_vars; ivar++) {
      nodals[ivar] = NULL;
      nodals_ordered[ivar] = NULL;
    }
    for (ivar = 0; ivar < num_vars; ivar++) {
      AIM_ALLOC(nodals[ivar]        , par.num_nodes, double, aimInfo, status);
      AIM_ALLOC(nodals_ordered[ivar], par.num_nodes, double, aimInfo, status);
    }
    AIM_ALLOC(node_id_map  , par.num_nodes, int   , aimInfo, status);

    time_step = ex_inquire_int(exoid, EX_INQ_TIME);

    for (ivar = 0; ivar < num_vars; ivar++) {
      status = ex_get_var(exoid, time_step, EX_NODAL, ivar+1, 1, par.num_nodes, nodals[ivar]);
      AIM_STATUS(aimInfo, status);
    }


    /* get node mapping to undo partitioning reordering */
    status = ex_get_id_map(exoid, EX_NODE_MAP, node_id_map);
    AIM_STATUS(aimInfo, status);


    /* Write binary meshb file */
    status = aim_file(aimInfo, "restart.solb", solFilename);
    AIM_STATUS(aimInfo, status);

    meshbVersion = 2;
    if (EXPORT_MESHB_VERTEX_3 < par.num_nodes) meshbVersion = 3;
    if (EXPORT_MESHB_VERTEX_4 < par.num_nodes) meshbVersion = 4;

    solbID = GmfOpenMesh(solFilename, GmfWrite, meshbVersion, par.num_dim);
    if (solbID == 0) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", solFilename);
      return CAPS_IOERR;
    }

    AIM_ALLOC(szfld, num_vars, int, aimInfo, status);
    for (ivar = 0; ivar < num_vars; ivar++) szfld[ivar] = GmfSca;

    status = GmfSetKwd(solbID, GmfSolAtVertices, par.num_nodes, num_vars, szfld);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    AIM_ALLOC(TypTab, num_vars, int, aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) {
      TypTab[ivar] = GmfDouble;
      szfld[ivar] = 1;

      /* properly order the scalar value */
      for (i = 0; i < par.num_nodes; i++) {
        nodals_ordered[ivar][node_id_map[i]-1] = nodals[ivar][i];
      }

      BegTab[ivar] = &nodals_ordered[ivar][0];
      EndTab[ivar] = &nodals_ordered[ivar][par.num_nodes-1];
    }

    /*@-nullpass@*/
    status = GmfSetBlock(solbID, GmfSolAtVertices, 1, par.num_nodes, 0, NULL, NULL,
                         GmfArgTab, TypTab, szfld, BegTab, EndTab);
    /*@+nullpass@*/
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    ex_close(exoid);
    exoid = 0;

    AIM_STRDUP(val->vals.string, solFilename, aimInfo, status);
  }

  status = CAPS_SUCCESS;

cleanup:
  if (var_names != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(var_names[ivar]);
  AIM_FREE(var_names);
  AIM_FREE(nodal);
  AIM_FREE(nodal_ordered);
  if (nodals != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(nodals[ivar]);
  if (nodals_ordered != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(nodals_ordered[ivar]);
  AIM_FREE(nodals);
  AIM_FREE(nodals_ordered);
  AIM_FREE(node_id_map);
  AIM_FREE(szfld);
  AIM_FREE(TypTab);

  if (exoid != 0) ex_close(exoid);
  if (solbID != 0) GmfCloseMesh(solbID);
  if (fp != NULL) fclose(fp);

  return status;
}


void aimCleanup(void *instStore)
{
  aimStorage *spinnakerInstance;

#ifdef DEBUG
  printf(" spinnakerAIM/aimCleanup!\n");
#endif
  spinnakerInstance = (aimStorage *) instStore;

  destroy_aimStorage(spinnakerInstance);

  AIM_FREE(spinnakerInstance);
}


/************************************************************************/
// CAPS transferring functions
void aimFreeDiscrPtr(void *ptr)
{
  AIM_FREE(ptr);
}


/************************************************************************/
int aimDiscr(char *tname, capsDiscr *discr)
{
    int i; // Indexing

    int status; // Function return status

    int numBody;
    //aimStorage *spinnakerInstance;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;
    capsValue *meshVal;

    // Volume Mesh obtained from meshing AIM
    aimMeshRef *meshRef;

#ifdef DEBUG
    printf(" spinnakerAIM/aimDiscr: tname = %s!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;
    //spinnakerInstance = (aimStorage *) discr->instStore;

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    AIM_STATUS(discr->aInfo, status);

    if (bodies == NULL) {
        AIM_ERROR(discr->aInfo, " spinnakerAIM/aimDiscr: No Bodies!\n");
        return CAPS_NOBODIES;
    }

    // Get the mesh input Value
    status = aim_getValue(discr->aInfo, inMesh, ANALYSISIN, &meshVal);
    AIM_STATUS(discr->aInfo, status);

    if (meshVal->nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(discr->aInfo, inMesh, "'Mesh' input must be linked to an output 'Area_Mesh' or 'Volume_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    meshRef = (aimMeshRef *)meshVal->vals.AIMptr;
    AIM_NOTNULL(meshRef, discr->aInfo, status);

    if (meshRef->nmap == 0 || meshRef->maps == NULL) {
        AIM_ERROR(discr->aInfo, "No surface mesh map in volume mesh - data transfer isn't possible.\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Do we have an individual surface mesh for each body
    if (meshRef->nmap != numBody) {
        AIM_ERROR(  discr->aInfo, "Number of surface mesh in the linked volume mesh (%d) does not match the number");
        AIM_ADDLINE(discr->aInfo,"of bodies (%d) - data transfer is NOT possible.", meshRef->nmap,numBody);
        status = CAPS_MISMATCH;
        goto cleanup;
    }

    // To this point is doesn't appear that the volume mesh has done anything bad to our surface mesh(es)

    // Store away the tessellation now
    AIM_ALLOC(tess, numBody, ego, discr->aInfo, status);
    for (i = 0; i < numBody; i++) {
        tess[i] = meshRef->maps[i].tess;
    }

#ifdef DEBUG
    printf(" spinnakerAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
#endif

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(tess);

    return status;
}


int aimLocateElement(capsDiscr *discr, double *params, double *param,
                     int *bIndex, int *eIndex, double *bary)
{
    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


int aimTransfer(/*@unused@*/ capsDiscr *discr, /*@unused@*/ const char *dataName,
    /*@unused@*/int numPoint,
    /*@unused@*/ int dataRank, /*@unused@*/ double *dataVal, /*@unused@*/ char **units)
{

    int status; // Function return status

    status = CAPS_SUCCESS;

//cleanup:
    return status;
}


int aimInterpolation(capsDiscr *discr, /*@unused@*/ const char *name,
                     int bIndex, int eIndex, double *bary, int rank,
                     double *data, double *result)
{
#ifdef DEBUG
    printf(" spinnakerAIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


int aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                      int bIndex, int eIndex, double *bary, int rank,
                      double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" spinnakerAIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


int aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                   int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" spinnakerAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                    int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" spinnakerAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
