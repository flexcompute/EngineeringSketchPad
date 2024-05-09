/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Sierra Structural Dynamics (SD) AIM
 *
 *     Written by Dr. Marshall Galbraith MIT
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 */


/*!\mainpage Introduction
 *
 * \section overviewSIERRASD Sierra AIM Overview
 * This module can be used to interface with the Sandia National Laboratories Sierra Mechanics structural analysis
 *  with geometry in the CAPS system. Sierra expects a mesh file and a corresponding
 * configuration file to perform the analysis.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsSIERRASD and \ref aimOutputsSU2, respectively.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferSIERRASD
 *
 * \subsection meshSIERRASD Automatic generation of Sierra Exodus mesh file
 * The mesh file from Sierra AIM is written in native Exodus
 * format ("filename.exo"). The description of the native Exodus mesh can be
 * found Exodus website (https://sandialabs.github.io/seacas-docs/html/index.html).
 * For the automatic generation of mesh file, Sierra AIM
 * depends on Mesh AIMs, for example, TetGen or AFLR4/3 AIM.
 *
 * \subsection configSIERRASD Automatic generation of Sierra input file
 * The input file ("input.i") is automatically
 * created by using the boundary conditions that were set in
 * the driver program as a user input. For the rest of the input
 * variables, default set of values are provided for a general execution.
 * If desired, a user has freedom to manually (a) change these variables based
 * on  personal preference, or (b) override the configuration file with unique
 * configuration variables.
 *
 * \section examplesSIERRASD Sierra Examples
 *  Here is an example that illustrated use of Sierra AIM \ref sierraExample. Note
 *  this AIM uses TetGen AIM for volume mesh generation.
 *
 */

#include <string.h>
#include <ctype.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h"    // Meshing utilities
#include "miscUtils.h"    // Miscellaneous utilities
#include "feaUtils.h"     // FEA utilities
#include "nastranUtils.h" // Nastran utilities

#include "exodusWriter.h"

#include <exodusII.h>

#define MAX(A,B)         (((A) < (B)) ? (B) : (A))

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#else
#include <unistd.h>
#endif

#define MXCHAR  255

//#define DEBUG


enum aimInputs
{
  Property = 1,        /* index is 1-based */
  Material,
  //Constraint,
  Load,
  Analysis_Type,
  Outputs,
  Input_String,
  OutputFieldFiles,
  Mesh_Morph,
  Mesh,
  NUMINPUT = Mesh       /* Total number of inputs */
};

#define NUMOUTPUT  0

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

  feaUnitsStruct units; // units system

  feaProblemStruct feaProblem;

  // Attribute to capsGroup index map
  mapAttrToIndexStruct groupMap;

  // Attribute to constraint index map
  mapAttrToIndexStruct constraintMap;

  // Attribute to load index map
  mapAttrToIndexStruct loadMap;

  // Attribute to transfer map
  mapAttrToIndexStruct transferMap;

  // Attribute to connect map
  mapAttrToIndexStruct connectMap;

  // Attribute to response map
  mapAttrToIndexStruct responseMap;

  // Mesh holders
  int numMesh;
  meshStruct *feaMesh;

  // Mesh reference obtained from meshing AIM
  aimMeshRef *meshRef, meshRefObj;

  int numOutputFieldFile;
  outputFieldStruct *outputFieldFiles;

} aimStorage;


static int initiate_aimStorage(aimStorage *sierraInstance)
{

    int status;

    status = initiate_feaUnitsStruct(&sierraInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->groupMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->constraintMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->loadMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for transfer to index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->transferMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for connect to index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->connectMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for response to index map
    status = initiate_mapAttrToIndexStruct(&sierraInstance->responseMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&sierraInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    // Mesh holders
    sierraInstance->numMesh = 0;
    sierraInstance->feaMesh = NULL;

    sierraInstance->meshRef = NULL;
    aim_initMeshRef(&sierraInstance->meshRefObj, aimUnknownMeshType);

    sierraInstance->numOutputFieldFile = 0;
    sierraInstance->outputFieldFiles = NULL;

    return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *sierraInstance)
{
    int i, status;

    status = destroy_feaUnitsStruct(&sierraInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->constraintMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->loadMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Transfer to index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->transferMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Connect to index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->connectMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Response to index map
    status = destroy_mapAttrToIndexStruct(&sierraInstance->responseMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (sierraInstance->feaMesh != NULL) {

        for (i = 0; i < sierraInstance->numMesh; i++) {
            status = destroy_meshStruct(&sierraInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        EG_free(sierraInstance->feaMesh);
    }

    sierraInstance->feaMesh = NULL;
    sierraInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&sierraInstance->feaProblem);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    sierraInstance->meshRef = NULL;
    aim_freeMeshRef(&sierraInstance->meshRefObj);


    for (i = 0; i < sierraInstance->numOutputFieldFile; i++) {
      (void) destroy_outputFieldStruct(&sierraInstance->outputFieldFiles[i]);
    }
    AIM_FREE(sierraInstance->outputFieldFiles);
    sierraInstance->numOutputFieldFile = 0;

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *sierraInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  status = fea_createMesh(aimInfo,
                          NULL,
                          0,
                          0,
                          (int)false,
                          &sierraInstance->groupMap,
                          &sierraInstance->constraintMap,
                          &sierraInstance->loadMap,
                          &sierraInstance->transferMap,
                          &sierraInstance->connectMap,
                          &sierraInstance->responseMap,
                          NULL,
                          &sierraInstance->numMesh,
                          &sierraInstance->feaMesh,
                          &sierraInstance->feaProblem );
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


static int
writeMaterial(/*@unused@*/ void *aimInfo, FILE *fp, const feaMaterialStruct *feaMaterial)
{
  int status = CAPS_SUCCESS;

  fprintf(fp, "MATERIAL %s\n", feaMaterial->name);
  fprintf(fp, "  density %16.12e\n", feaMaterial->density);
  fprintf(fp, "  E       %16.12e\n", feaMaterial->youngModulus);
  fprintf(fp, "  Nu      %16.12e\n", feaMaterial->poissonRatio);
  fprintf(fp, "END\n");
  fprintf(fp, "\n");

//cleanup:
  return status;
}

static int
writeBlock(/*@unused@*/ void *aimInfo, FILE *fp, const feaPropertyStruct *feaProperty)
{
  int status = CAPS_SUCCESS;
  char name[PATH_MAX];

  strcpy(name, feaProperty->name);
  string_toLowerCase(name);

  fprintf(fp, "BLOCK %s\n", name);
  fprintf(fp, "  material  %s\n", feaProperty->materialName);
  fprintf(fp, "  thickness %16.12e\n", feaProperty->membraneThickness);
  fprintf(fp, "END\n");
  fprintf(fp, "\n");

//cleanup:
  return status;
}

#if 0
static int
writeLoad(void *aimInfo, FILE *fp, const feaPropertyStruct *feaProperty)
{
  int status = CAPS_SUCCESS;
  char name[PATH_MAX];

  strcpy(name, feaProperty->name);

  fprintf(fp, "BLOCK %s\n", string_toLowerCase(name));
  fprintf(fp, "  material  %s\n", feaProperty->materialName);
  fprintf(fp, "  thickness %16.12e\n", feaProperty->membraneThickness);
  fprintf(fp, "END\n");
  fprintf(fp, "\n");

//cleanup:
  return status;
}
#endif

static int
writeInput(void *aimInfo, const aimStorage *sierraInstance, capsValue *aimInputs)
{
  int status = CAPS_SUCCESS;
  int i, j, k;
  int found;
  size_t slen;
  const char inputFile[] = "input.inp";
  FILE *fp = NULL;

  // Output filename
  char meshRefFilename[PATH_MAX];
  char meshFilename[PATH_MAX];

  const char *analysisType = NULL;

  // Mesh obtained from meshing AIM
  aimMeshRef *meshRef = NULL;

  fp = aim_fopen(aimInfo, inputFile, "w");
  if (fp == NULL) {
      AIM_ERROR(aimInfo, "Unable to open file: %s", inputFile);
      status = CAPS_IOERR;
      goto cleanup;
  }

  // Get mesh
  meshRef = sierraInstance->meshRef;
  AIM_NOTNULL(meshRef, aimInfo, status);

  /* create a symbolic link to the file name*/
  snprintf(meshRefFilename, PATH_MAX, "%s%s", meshRef->fileName, MESHEXTENSION);
  snprintf(meshFilename, PATH_MAX, "sierraMesh%s", MESHEXTENSION);
  status = aim_symLink(aimInfo, meshRefFilename, meshFilename);
  AIM_STATUS(aimInfo, status);

  // Analysis type
  analysisType = aimInputs[Analysis_Type-1].vals.string;

  fprintf(fp, "SOLUTION\n");
  fprintf(fp, "\n");
  fprintf(fp, "  title 'CAPS Sierra Structural Dynamics'\n");
  fprintf(fp, "\n");
  fprintf(fp, "  case 'analysis'\n");
  // Analysis type
  if(strcasecmp(analysisType, "Static") == 0) fprintf(fp, "    statics\n");
  else {
      AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s", analysisType);
      status = CAPS_BADVALUE;
      goto cleanup;
  }
  fprintf(fp, "  solver gdsw\n");
  fprintf(fp, "END\n");
  fprintf(fp, "\n");
  fprintf(fp, "FILE\n");
  fprintf(fp, "  geometry_file '%s'\n", meshFilename);
  fprintf(fp, "END\n");
  fprintf(fp, "\n");
  fprintf(fp, "PARAMETERS\n");
  fprintf(fp, "  wtMass 1.0\n");
  fprintf(fp, "END\n");
  fprintf(fp, "\n");

  if (aimInputs[Outputs-1].nullVal != IsNull) {
      fprintf(fp, "OUTPUTS\n");
      // write all outputs requested
      for (slen = i = 0; i < aimInputs[Outputs-1].length; i++) {
          found = (int) false;
          for (j = 0; j < sierraInstance->numOutputFieldFile && found == (int) false; j++)
            for (k = 0; k < sierraInstance->outputFieldFiles[j].numField && found == (int) false; k++)
              if (strcasecmp(sierraInstance->outputFieldFiles[j].field[k], aimInputs[Outputs-1].vals.string + slen))
                found = (int) true;

          if (found == (int) false)
            fprintf(fp,"  %s\n", aimInputs[Outputs-1].vals.string + slen);
          slen += strlen(aimInputs[Outputs-1].vals.string + slen) + 1;
      }
      // Write all outputs requested in files
      for (j = 0; j < sierraInstance->numOutputFieldFile; j++)
        for (k = 0; k < sierraInstance->outputFieldFiles[j].numField; k++)
          fprintf(fp,"  %s\n", sierraInstance->outputFieldFiles[j].field[k]);
      fprintf(fp, "END\n");
  }

  // Materials
  for (i = 0; i < sierraInstance->feaProblem.numMaterial; i++) {
      if (i == 0) printf("\tWriting material\n");
      status = writeMaterial(aimInfo, fp, &sierraInstance->feaProblem.feaMaterial[i]);
      AIM_STATUS(aimInfo, status);
  }

  // Block
  for (i = 0; i < sierraInstance->feaProblem.numProperty; i++) {
      if (i == 0) printf("\tWriting property cards\n");
      status = writeBlock(aimInfo, fp, &sierraInstance->feaProblem.feaProperty[i]);
      AIM_STATUS(aimInfo, status);
  }

  if (aimInputs[Input_String-1].nullVal != IsNull) {
      fprintf(fp, "\n");
      fprintf(fp,"// CAPS Input_String\n\n");
      for (slen = i = 0; i < aimInputs[Input_String-1].length; i++) {
          fprintf(fp,"%s\n", aimInputs[Input_String-1].vals.string + slen);
          slen += strlen(aimInputs[Input_String-1].vals.string + slen) + 1;
      }
  }
  fprintf(fp, "\n");

cleanup:

  if (fp != NULL) fclose(fp);

  return status;
}


// Get the support properties from a capsTuple
int output_getFieldFiles(void *aimInfo,
                         int numFieldFile,
                         capsTuple fieldFiles[],
                         int *numOutput,
                         outputFieldStruct **outputFields) {

    /*! \page outputField Field output files
     * Structure for the output filed tuple  = {"File Name": "Value"}.
     * "File Name" defines the file name for where filed quantities will be written.
     *  The "Value" can either be a single string a list of single strings for the field quantities to be written to the file.
     */

    int status; //Function return

    int i;

    // Destroy our support structures coming in if aren't 0 and NULL already
    for (i = 0; i < *numOutput; i++) {
      status = destroy_outputFieldStruct(*outputFields + i);
      AIM_STATUS(aimInfo, status);
    }
    AIM_FREE(*outputFields);

    // Allocate
    AIM_ALLOC((*outputFields), numFieldFile, outputFieldStruct, aimInfo, status);
    *numOutput = numFieldFile;

    // Initiate supports to default values
    for (i = 0; i < *numOutput; i++) {
        status = initiate_outputFieldStruct(*outputFields + i);
        AIM_STATUS(aimInfo, status);
    }

    // Loop through tuples and fill out the support structures
    for (i = 0; i < numFieldFile; i++) {

        // Set file name to tuple attribute name
        AIM_STRDUP((*outputFields)[i].filename, fieldFiles[i].name, aimInfo, status);

        status = string_toStringDynamicArray(fieldFiles[i].value, &(*outputFields)[i].numField, &(*outputFields)[i].field);
        AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;
cleanup:
    return status;
}



/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@null@*/ /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS; // Function return
    int i;
    int *ints=NULL;
    char **strs=NULL;
    //const char *keyWord;
    //char *keyValue = NULL;
    //double real = 1;

    aimStorage *sierraInstance=NULL;

    #ifdef DEBUG
        printf("\n sierraAIM/aimInitialize   inst = %d!\n", inst);
    #endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 1;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
    strs[0]  = EG_strdup("Pressure");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
    *fnames  = strs;


    /* specify the dimension of each field variable */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = 1;
    *franks   = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate sierraInstance
    AIM_ALLOC(sierraInstance, 1, aimStorage, aimInfo, status);
    *instStore = sierraInstance;

    // Set initial values for sierraInstance

    initiate_aimStorage(sierraInstance);

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
    //aimStorage *sierraInstance;

#ifdef DEBUG
    printf(" sierraAIM/aimInputs  index = %d!\n", index);
#endif

    *ainame = NULL;

    // Sierra Inputs
    /*! \page aimInputsSIERRASD AIM Inputs
     */

    //sierraInstance = (aimStorage *) instStore;
    if (index == Property) {
        *ainame              = EG_strdup("Property");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsSIERRASD
         * - <B> Property = NULL</B> <br>
         * Property tuple used to input property information for the model, see \ref feaProperty for additional details.
         */
    } else if (index == Material) {
        *ainame              = EG_strdup("Material");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsSIERRASD
         * - <B> Material = NULL</B> <br>
         * Material tuple used to input material information for the model, see \ref feaMaterial for additional details.
         */
    } else if (index == Load) {
        *ainame              = EG_strdup("Load");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsSIERRASD
         * - <B> Load = NULL</B> <br>
         * Load tuple used to input load information for the model, see \ref feaLoad for additional details.
         */

    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Static");
        defval->lfixed       = Change;

        /*! \page aimInputsSIERRASD
         * - <B> Analysis_Type = "Static"</B> <br>
         * Type of analysis to generate files for, options include "Modal" or "Static".
         */

    } else if (index == Outputs) {
        *ainame               = EG_strdup("Outputs");
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->lfixed       = Change;
        defval->nrow         = 0;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

       /*! \page aimInputsSIERRASD
         * - <B>Outputs = NULL</B> <br>
         * Array of input strings of output quantities.
        */

    } else if (index == Input_String) {
        *ainame               = EG_strdup("Input_String");
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->lfixed       = Change;
        defval->nrow         = 0;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

       /*! \page aimInputsSIERRASD
         * - <B>Input_String = NULL</B> <br>
         * Array of input strings that will be written as is to the end of the Sierra input.i file.
        */

    } else if (index == OutputFieldFiles) {
        *ainame               = EG_strdup("OutputFieldFiles");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

       /*! \page aimInputsSIERRASD
         * - <B>OutputFieldFiles = NULL</B> <br>
         * Tuple of output file names and field quantities to be written to those files.
        */

    } else if (index == Mesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        defval->type         = Boolean;
        defval->lfixed       = Fixed;
        defval->vals.integer = (int) false;
        defval->dim          = Scalar;
        defval->nullVal      = NotNull;

        /*! \page aimInputsSIERRASD
         * - <B> Mesh_Morph = False</B> <br>
         * Project previous surface mesh onto new geometry and write out a 'Proj_Name'_body#.dat file.
         */

    } else if (index == Mesh) {
        *ainame             = AIM_NAME(Mesh);
        defval->type        = PointerMesh;
        defval->nrow        = 1;
        defval->lfixed      = Fixed;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;
        AIM_STRDUP(defval->meshWriter, MESHWRITER, aimInfo, status);

        /*! \page aimInputsSIERRASD
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

    aimStorage *sierraInstance;

    sierraInstance = (aimStorage *) instStore;
    AIM_NOTNULL(sierraInstance, aimInfo, status);
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Free our meshRef
    (void) aim_freeMeshRef(&sierraInstance->meshRefObj);

    if (aimInputs[Mesh-1].nullVal == IsNull &&
        aimInputs[Mesh_Morph-1].vals.integer == (int) false) {
        AIM_ANALYSISIN_ERROR(aimInfo, Mesh, "'Mesh' input must be linked to a 'Volume_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    // Get mesh
    sierraInstance->meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;

    if ( aimInputs[Mesh_Morph-1].vals.integer == (int) true &&
        sierraInstance->meshRef == NULL) { // If we are mighty morphing

        // Lets "load" the meshRef now since it's not linked
        status = aim_loadMeshRef(aimInfo, &sierraInstance->meshRefObj);
        AIM_STATUS(aimInfo, status);

        // Mightly Morph the mesh
        status = aim_morphMeshUpdate(aimInfo, &sierraInstance->meshRefObj, numBody, bodies);
        AIM_STATUS(aimInfo, status);
        /*@-immediatetrans@*/
        sierraInstance->meshRef = &sierraInstance->meshRefObj;
        /*@+immediatetrans@*/
    }
    AIM_NOTNULL(sierraInstance->meshRef, aimInfo, status);

    // Get FEA mesh if we don't already have one
    status = checkAndCreateMesh(aimInfo, sierraInstance);
    AIM_STATUS(aimInfo, status);

    // Note: Setting order is important here.
    // 1. Materials should be set before properties.
    // 2. Coordinate system should be set before mesh and loads
    // 3. Mesh should be set before loads, constraints, supports, and connections
    // 4. Constraints and loads should be set before analysis
    // 5. Optimization should be set after properties, but before analysis

    // Set material properties
    if (aimInputs[Material-1].nullVal == NotNull) {
        status = fea_getMaterial(aimInfo,
                                 aimInputs[Material-1].length,
                                 aimInputs[Material-1].vals.tuple,
                                 &sierraInstance->units,
                                 &sierraInstance->feaProblem.numMaterial,
                                 &sierraInstance->feaProblem.feaMaterial);
        AIM_STATUS(aimInfo, status);
    } else printf("Material tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &sierraInstance->groupMap,
                                 &sierraInstance->units,
                                 &sierraInstance->feaProblem);
        AIM_STATUS(aimInfo, status);


        // Assign element "subtypes" based on properties set
        status = fea_assignElementSubType(sierraInstance->feaProblem.numProperty,
                                          sierraInstance->feaProblem.feaProperty,
                                          &sierraInstance->feaProblem.feaMesh);
        AIM_STATUS(aimInfo, status);

    } else printf("Property tuple is NULL - No properties set\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInfo,
                             aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &sierraInstance->loadMap,
                             &sierraInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Load tuple is NULL - No loads applied\n");

    // Set output files
    if (aimInputs[OutputFieldFiles-1].nullVal == NotNull) {
        status = output_getFieldFiles(aimInfo,
                                      aimInputs[OutputFieldFiles-1].length,
                                      aimInputs[OutputFieldFiles-1].vals.tuple,
                                      &sierraInstance->numOutputFieldFile,
                                      &sierraInstance->outputFieldFiles);
        AIM_STATUS(aimInfo, status);
    } else printf("OutputFieldFiles tuple is NULL - No output field files set\n");


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

  // EGADS return values
  const char   *intents;
  ego body;
  int state, nGlobal, ptype, pindex;
  double xyz[3];
  FILE *fp = NULL;

  // Output filename
  char meshfilename[PATH_MAX];

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
  double *pressure=NULL;
  int dim, nVertex;

  char *const variable_names[] = {"pressure"};
  int whole_time_step = 1;

  feaLoadStruct *feaLoad=NULL;

  const aimStorage *sierraInstance;

  sierraInstance = (const aimStorage *) instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

#ifdef DEBUG
  printf(" sierraAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif
  if ((numBody <= 0) || (bodies == NULL)) {
    AIM_ERROR(aimInfo, "No Bodies!");
    return CAPS_SOURCEERR;
  }

  // write input.i
  status = writeInput(aimInfo, sierraInstance, aimInputs);
  AIM_STATUS(aimInfo, status);

  snprintf(meshfilename, PATH_MAX, "%s%s", sierraInstance->meshRef->fileName, MESHEXTENSION);

  if ( aimInputs[Mesh_Morph-1].vals.integer == (int) true) { // If we are mighty morphing
    if (aimInputs[Mesh-1].nullVal == NotNull) {
      // store the current mesh for future iterations
      status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[Mesh-1].vals.AIMptr, MESHEXTENSION);
      AIM_STATUS(aimInfo, status);

    } else {

      exoid = ex_open(meshfilename, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", meshfilename);
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

      for (imap = 0; imap < sierraInstance->meshRef->nmap; imap++) {
        status = EG_statusTessBody(sierraInstance->meshRef->maps[imap].tess, &body, &state, &nGlobal);
        AIM_STATUS(aimInfo, status);

        // Write the map
        for (i = 0; i < nGlobal; i++) {
          status = EG_getGlobal(sierraInstance->meshRef->maps[imap].tess, i+1, &ptype, &pindex, xyz);
          AIM_STATUS(aimInfo, status);

          j = sierraInstance->meshRef->maps[imap].map[i]-1;

          x[j] = xyz[0];
          y[j] = xyz[1];
          if (z != NULL)
            z[j] = xyz[2];
        }
      }

      exoid = ex_open(meshfilename, EX_WRITE | EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC,
                      &CPU_word_size, &IO_word_size, &version);
      if (exoid <= 0) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", meshfilename);
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

  if (sierraInstance->feaProblem.numLoad > 0) {
      AIM_ALLOC(feaLoad, sierraInstance->feaProblem.numLoad, feaLoadStruct, aimInfo, status);
      for (i = 0; i < sierraInstance->feaProblem.numLoad; i++) initiate_feaLoadStruct(&feaLoad[i]);
      for (i = 0; i < sierraInstance->feaProblem.numLoad; i++) {
          status = copy_feaLoadStruct(aimInfo, &sierraInstance->feaProblem.feaLoad[i], &feaLoad[i]);
          AIM_STATUS(aimInfo, status);

          if (feaLoad[i].loadType == PressureExternal) {

              // Transfer external pressures from the AIM discrObj
              status = fea_transferExternalPressureNode(aimInfo,
                                                        &feaLoad[i]);
              AIM_STATUS(aimInfo, status);

              if (pressure == NULL) {
                AIM_ALLOC(pressure, sierraInstance->feaProblem.feaMesh.numNode, double, aimInfo, status);
                for (j = 0; j < sierraInstance->feaProblem.feaMesh.numNode; j++) pressure[j] = 0;
              }

              for (j = 0; j < feaLoad[i].numGridID; j++) {
                pressure[feaLoad[i].gridIDSet[j]-1] = feaLoad->pressureMultiDistributeForce[j];
              }

          } else if (feaLoad[i].loadType == ThermalExternal) {

              // Transfer external temperature from the AIM discrObj
              status = fea_transferExternalTemperature(aimInfo,
                                                       &feaLoad[i]);
              AIM_STATUS(aimInfo, status);
          }
      }

      if (pressure != NULL) {

        exoid = ex_open(meshfilename, EX_WRITE | EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC,
                        &CPU_word_size, &IO_word_size, &version);
        if (exoid <= 0) {
          AIM_ERROR(aimInfo, "Cannot open file: %s\n", meshfilename);
          status = CAPS_IOERR;
          goto cleanup;
        }

        status = ex_put_variable_param(exoid, EX_NODAL, 1);
        AIM_STATUS(aimInfo, status);

        status = ex_put_variable_names(exoid, EX_NODAL, 1, variable_names);
        AIM_STATUS(aimInfo, status);

        status = ex_put_var(exoid, whole_time_step, EX_NODAL, 1, 1, sierraInstance->feaProblem.feaMesh.numNode, pressure);
        AIM_STATUS(aimInfo, status);
        AIM_FREE(pressure);

        ex_close(exoid);
        exoid = 0;
      }

  }


  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(x);
  AIM_FREE(y);
  AIM_FREE(z);
  AIM_FREE(pressure);

  if (feaLoad != NULL) {
      for (i = 0; i < sierraInstance->feaProblem.numLoad; i++) {
          destroy_feaLoadStruct(&feaLoad[i]);
      }
      AIM_FREE(feaLoad);
  }

  if (fp != NULL) fclose(fp);
  if (exoid > 0) ex_close(exoid);

  return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *aimInputs)
{
  int status = CAPS_SUCCESS;
  int i, j, ivar;

  char meshFilename[PATH_MAX];
  char aimFile[PATH_MAX];

  const char *analysisType = NULL;

  int CPU_word_size = sizeof(double);
  int IO_word_size = sizeof(double);
  float version;
  int exoid = 0;
  ex_init_params par;

  int num_vars=0;
  char **var_names=NULL;
  int whole_time_step = 1;
  double *elemental = NULL;
  int *ids=NULL, *nodeconn=NULL, *numElemAvg=NULL;
  int num_elem_in_blk, num_nodes_per_elem, num_attr, iblk;
  char elem_type[MAX_STR_LENGTH+1];
  char name[2*MAX_STR_LENGTH+1];


  capsValue val;
 // FILE *fp = NULL;

  // Mesh reference obtained from meshing AIM
  //aimMeshRef *meshRef = NULL;

  //aimStorage *sierraInstance = (aimStorage*)instStore;

  status = aim_initValue(&val);
  AIM_STATUS(aimInfo, status);

  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get mesh
//  meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;
//  AIM_NOTNULL(meshRef, aimInfo, status);

  if ((aimInputs[Outputs-1].nullVal == NotNull ||
       aimInputs[OutputFieldFiles-1].nullVal == NotNull) &&
       restart == 0) {

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Analysis type
    if(strcasecmp(analysisType, "Static") == 0)
      snprintf(meshFilename, PATH_MAX, "sierraMesh-analysis%s", MESHEXTENSION);
    else {
      AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s", analysisType);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    status = aim_file(aimInfo, meshFilename, aimFile);
    AIM_STATUS(aimInfo, status);

    exoid = ex_open(aimFile, EX_READ | EX_NETCDF4 | EX_NOCLASSIC,
                    &CPU_word_size, &IO_word_size, &version);
    if (exoid <= 0) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", meshFilename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    ex_opts(EX_VERBOSE | EX_DEBUG | EX_NULLVERBOSE);

    status = ex_get_init_ext(exoid, &par);
    AIM_STATUS(aimInfo, status);

    status = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(var_names, num_vars, char*, aimInfo, status);
    AIM_ALLOC(val.vals.reals, par.num_nodes, double, aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) var_names[ivar] = NULL;
    for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(var_names[ivar], EX_MAX_NAME, char, aimInfo, status);

    status = ex_get_variable_names(exoid, EX_NODAL, num_vars, (char**)var_names);
    AIM_STATUS(aimInfo, status);
#if 0
    for (i = 0; i < sierraInstance->numOutputFieldFile; i++)
    {
      fp = aim_fopen(aimInfo, sierraInstance->outputFieldFiles[i].filename, "w");
      if (fp == NULL) {
        AIM_ERROR(aimInfo, "Could not open '%s'!", sierraInstance->outputFieldFiles[i].filename);
        status = CAPS_IOERR;
        goto cleanup;
      }

      fprintf(fp, "MeshVersionFormatted 2\n");
      fprintf(fp, "Dimension 3\n");
      fprintf(fp, "SolAtVertices %lld 1 1\n", (long long int)par.num_nodes);

      for (j = 0; j < sierraInstance->outputFieldFiles[i].numField; j++)
      {
        ivar = 0;
        while (ivar < num_vars && strcasecmp(var_names[ivar], sierraInstance->outputFieldFiles[i].field[j])) ivar++;
        if (ivar == num_vars) {
          AIM_ERROR(aimInfo, "Could not find '%s' in Sierra SD output file!", sierraInstance->outputFieldFiles[i].field[j]);
          status = CAPS_IOERR;
          goto cleanup;
        }

        status = ex_get_var(exoid, whole_time_step, EX_NODAL, ivar+1, 1, par.num_nodes, val.vals.reals);
        AIM_STATUS(aimInfo, status);

        for (k = 0; k < par.num_nodes; k++)
          fprintf(fp, "%16.12e\n", val.vals.reals[k]);
      }

      fclose(fp);
      fp = NULL;
    }
#endif
    AIM_FREE(val.vals.reals);

    for (ivar = 0; ivar < num_vars; ivar++) {

      val.dim  = Vector;
      val.type = Double;
      val.nrow = par.num_nodes;
      val.ncol = 1;
      val.length = val.nrow*val.ncol;

      AIM_ALLOC(val.vals.reals, par.num_nodes, double, aimInfo, status);

      status = ex_get_var(exoid, whole_time_step, EX_NODAL, ivar+1, 1, par.num_nodes, val.vals.reals);
      AIM_STATUS(aimInfo, status);

      /* create the dynamic output */
      status = aim_makeDynamicOutput(aimInfo, var_names[ivar], &val);
      AIM_STATUS(aimInfo, status);
    }

    for (ivar = 0; ivar < num_vars; ivar++) AIM_FREE(var_names[ivar]);
    AIM_FREE(var_names);

    status = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_vars);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(var_names, num_vars, char*, aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) var_names[ivar] = NULL;
    for (ivar = 0; ivar < num_vars; ivar++) AIM_ALLOC(var_names[ivar], EX_MAX_NAME, char, aimInfo, status);

    status = ex_get_variable_names(exoid, EX_ELEM_BLOCK, num_vars, (char**)var_names);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(elemental, par.num_elem, double, aimInfo, status);
    AIM_ALLOC(numElemAvg, par.num_nodes, int, aimInfo, status);

    AIM_ALLOC(ids, par.num_elem_blk, int, aimInfo, status);

    status = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
    AIM_STATUS(aimInfo, status);

    for (ivar = 0; ivar < num_vars; ivar++) {

      AIM_ALLOC(val.vals.reals, par.num_nodes, double, aimInfo, status);

      for (i = 0; i < par.num_nodes; i++) {
        numElemAvg[i] = 0;
        val.vals.reals[i] = 0;
      }

      for (iblk = 0; iblk < par.num_elem_blk; iblk++)
      {
        /*@-nullpass@*/
        status = ex_get_block(exoid, EX_ELEM_BLOCK, ids[iblk], elem_type,
                              &num_elem_in_blk, &num_nodes_per_elem, NULL, NULL, &num_attr);
        AIM_STATUS(aimInfo, status);

        AIM_ALLOC(nodeconn, num_elem_in_blk*num_nodes_per_elem, int, aimInfo, status);

        status = ex_get_conn(exoid, EX_ELEM_BLOCK, ids[iblk], nodeconn, NULL, NULL);
        AIM_STATUS(aimInfo, status);
        /*@+nullpass@*/

        status = ex_get_var(exoid, whole_time_step, EX_ELEM_BLOCK, ivar+1, ids[iblk], num_elem_in_blk, elemental);
        AIM_STATUS(aimInfo, status);

        for (i = 0; i < num_elem_in_blk; i++)
        {
          for (j = 0; j < num_nodes_per_elem; j++) {
            val.vals.reals[nodeconn[i*num_nodes_per_elem + j]-1] += elemental[i];
            numElemAvg[nodeconn[i*num_nodes_per_elem + j]-1] += 1;
          }
        }

        AIM_FREE(nodeconn);
      }

      for (i = 0; i < par.num_nodes; i++) {
        if (numElemAvg[i] > 0)
          val.vals.reals[i] /= numElemAvg[i];
      }

      snprintf(name, 2*MAX_STR_LENGTH+1, "%s_Grid", var_names[ivar]);

      //printf("%s\n", name);

      val.dim  = Vector;
      val.type = Double;
      val.nrow = par.num_nodes;
      val.ncol = 1;
      val.length = val.nrow*val.ncol;

      /* create the dynamic output */
      status = aim_makeDynamicOutput(aimInfo, name, &val);
      AIM_STATUS(aimInfo, status);
    }

    ex_close(exoid);
    exoid = 0;
  }

  status = CAPS_SUCCESS;

cleanup:
  if (var_names != NULL)
    for (ivar = 0; ivar < num_vars; ivar++)
      AIM_FREE(var_names[ivar]);
  AIM_FREE(var_names);
  AIM_FREE(val.vals.reals);
  AIM_FREE(elemental);
  AIM_FREE(ids);
  AIM_FREE(nodeconn);
  AIM_FREE(numElemAvg);

  if (exoid != 0) ex_close(exoid);
  //if (fp != NULL) fclose(fp);

  return status;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
	// SU2 Outputs
    /*! \page aimOutputsSIERRASD AIM Outputs
     * Sierra outputs
     */

    #ifdef DEBUG
        printf(" sierraAIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
    #endif

    return CAPS_SUCCESS;
}


// Calculate Sierra output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, capsValue *val)
{
    int status = CAPS_SUCCESS;

    //aimStorage *sierraInstance;

#ifdef DEBUG
    printf(" sierraAIM/aimCalcOutput  index = %d!\n", index);
#endif
    //sierraInstance = (aimStorage *) instStore;

    val->vals.real = 0.0; // Set default value

    status = CAPS_SUCCESS;
//cleanup:

    return status;
}


void aimCleanup(void *instStore)
{
    aimStorage *sierraInstance;

#ifdef DEBUG
    printf(" sierraAIM/aimCleanup!\n");
#endif
    sierraInstance = (aimStorage *) instStore;

    destroy_aimStorage(sierraInstance);

    AIM_FREE(sierraInstance);
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
    //aimStorage *sierraInstance;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;
    capsValue *meshVal;

    // Volume Mesh obtained from meshing AIM
    aimMeshRef *meshRef;
    aimStorage *sierraInstance;


#ifdef DEBUG
    printf(" sierraAIM/aimDiscr: tname = %s!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;
    sierraInstance = (aimStorage *) discr->instStore;

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    AIM_STATUS(discr->aInfo, status);

    if (bodies == NULL) {
        AIM_ERROR(discr->aInfo, " sierraAIM/aimDiscr: No Bodies!\n");
        return CAPS_NOBODIES;
    }

    // Get the mesh input Value
    status = aim_getValue(discr->aInfo, Mesh, ANALYSISIN, &meshVal);
    AIM_STATUS(discr->aInfo, status);

    if (meshVal->nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(discr->aInfo, Mesh, "'Mesh' input must be linked to an output 'Area_Mesh' or 'Volume_Mesh'");
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

    status = mesh_fillDiscr(tname, &sierraInstance->groupMap, sierraInstance->numMesh, tess, discr);
    AIM_STATUS(discr->aInfo, status);

#ifdef DEBUG
    printf(" sierraAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
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
    printf(" sierraAIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


int aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                      int bIndex, int eIndex, double *bary, int rank,
                      double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" sierraAIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


int aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                   int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" sierraAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                    int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" sierraAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
