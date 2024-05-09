/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Plato AIM
 *
 *     Written by Dr. Marshall Galbraith MIT
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 */


/*!\mainpage Introduction
 *
 * \section overviewPLATO Plato AIM Overview
 * This module can be used to interface with the open-source Plato code developed at Sandia National Laboratories
 *  with geometry in the CAPS system. For Plato capabilities and related documentation, please refer to
 * https://github.com/platoengine. Plato expects a volume, surface, or area mesh and a corresponding
 * configuration file to perform the analysis.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsPLATO and \ref aimOutputsPLATO, respectively.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferPLATO
 *
 * \subsection meshPLATO Automatic generation of Plato Exodus Mesh file
 * The mesh file from Plato AIM is written in native Exodus
 * format ("filename.exo"). The description of the native Exodus mesh can be
 * found Exodus website (https://sandialabs.github.io/seacas-docs/html/index.html).
 * For the automatic generation of mesh file, Plato AIM
 * depends on Mesh AIMs, for example, TetGen or AFLR4/3 AIM.
 *
 * \section examplesPLATO Plato Examples
 *  Here is an example that illustrated use of Plato AIM \ref platoExample. Note
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
  Proj_Name = 1,        /* index is 1-based */
  Mesh_Morph,
  Mesh,
  NUMINPUT = Mesh       /* Total number of inputs */
};

#define NUMOUTPUT  0

typedef struct {

  // Mesh reference obtained from meshing AIM
  aimMeshRef *meshRef, meshRefObj;

} aimStorage;

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

    aimStorage *platoInstance=NULL;

    #ifdef DEBUG
        printf("\n platoAIM/aimInitialize   inst = %d!\n", inst);
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

    // Allocate platoInstance
    AIM_ALLOC(platoInstance, 1, aimStorage, aimInfo, status);
    *instStore = platoInstance;

    // Set initial values for platoInstance

    platoInstance->meshRef = NULL;
    aim_initMeshRef(&platoInstance->meshRefObj, aimUnknownMeshType);

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
    //aimStorage *platoInstance;

#ifdef DEBUG
    printf(" platoAIM/aimInputs  index = %d!\n", index);
#endif

    *ainame = NULL;

    // Plato Inputs
    /*! \page aimInputsPLATO AIM Inputs
     */

    //platoInstance = (aimStorage *) instStore;
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("plato_CAPS");

        /*! \page aimInputsPLATO
         * - <B> Proj_Name = "plato_CAPS"</B> <br>
         * This corresponds to the project name used for file naming.
         */

    } else if (index == Mesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        defval->type         = Boolean;
        defval->lfixed       = Fixed;
        defval->vals.integer = (int) false;
        defval->dim          = Scalar;
        defval->nullVal      = NotNull;

        /*! \page aimInputsPLATO
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

        /*! \page aimInputsPLATO
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

    aimStorage *platoInstance;

    platoInstance = (aimStorage *) instStore;
    AIM_NOTNULL(platoInstance, aimInfo, status);
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Free our meshRef
    (void) aim_freeMeshRef(&platoInstance->meshRefObj);

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
    platoInstance->meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;

    if ( aimInputs[Mesh_Morph-1].vals.integer == (int) true &&
        platoInstance->meshRef == NULL) { // If we are mighty morphing

        // Lets "load" the meshRef now since it's not linked
        status = aim_loadMeshRef(aimInfo, &platoInstance->meshRefObj);
        AIM_STATUS(aimInfo, status);

        // Mightly Morph the mesh
        status = aim_morphMeshUpdate(aimInfo, &platoInstance->meshRefObj, numBody, bodies);
        AIM_STATUS(aimInfo, status);
        /*@-immediatetrans@*/
        platoInstance->meshRef = &platoInstance->meshRefObj;
        /*@+immediatetrans@*/
    }
    AIM_NOTNULL(platoInstance->meshRef, aimInfo, status);

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
  int state, np, nGlobal, ptype, pindex;
  double xyz[3];
  FILE *fp = NULL;

  // Output filename
  char meshfilename[PATH_MAX];
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

  // Volume Mesh obtained from meshing AIM
  aimMeshRef *meshRef = NULL;

  const aimStorage *platoInstance;

  platoInstance = (const aimStorage *) instStore;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

#ifdef DEBUG
  printf(" platoAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif
  if ((numBody <= 0) || (bodies == NULL)) {
    AIM_ERROR(aimInfo, "No Bodies!");
    return CAPS_SOURCEERR;
  }

  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get mesh
  meshRef = platoInstance->meshRef;
  AIM_NOTNULL(meshRef, aimInfo, status);

  if ( aimInputs[Mesh_Morph-1].vals.integer == (int) true) { // If we are mighty morphing
    if (aimInputs[Mesh-1].nullVal == NotNull) {
      // store the current mesh for future iterations
      status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[Mesh-1].vals.AIMptr, MESHEXTENSION);
      AIM_STATUS(aimInfo, status);

    } else {

      snprintf(meshfilename, PATH_MAX, "%s%s", meshRef->fileName, MESHEXTENSION);

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

  /* create a symbolic link to the file name*/
  snprintf(meshfilename, PATH_MAX, "%s%s", meshRef->fileName, MESHEXTENSION);
  snprintf(filepath, PATH_MAX, "%s%s", aimInputs[Proj_Name-1].vals.string, MESHEXTENSION);
  status = aim_symLink(aimInfo, meshfilename, filepath);
  AIM_STATUS(aimInfo, status);

  fp = aim_fopen(aimInfo, "sensMap.txt", "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Failed to open sensMap.txt'");
    status = CAPS_IOERR;
    goto cleanup;
  }

  // Write the number of maps
  fprintf(fp, "%d\n", meshRef->nmap);

  for (imap = 0; imap < meshRef->nmap; imap++) {
    status = EG_statusTessBody(meshRef->maps[imap].tess, &body, &state, &np);
    AIM_STATUS(aimInfo, status);

    // Write number of points in the map
    fprintf(fp, "%d\n", np);

    // Write the map
    for (i = 0; i < np; i++)
      fprintf(fp, "%d\n", meshRef->maps[imap].map[i]);


    snprintf(meshfilename, PATH_MAX, "%s_%d.eto", aimInputs[Proj_Name-1].vals.string, imap+1);
    status = aim_file(aimInfo, meshfilename, filepath);
    AIM_STATUS(aimInfo, status);

    remove(filepath);
    status = EG_saveTess(meshRef->maps[imap].tess, filepath);
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(x);
  AIM_FREE(y);
  AIM_FREE(z);

  if (fp != NULL) fclose(fp);
  if (exoid > 0) ex_close(exoid);

  return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *aimInputs)
{
  int status = CAPS_SUCCESS;

  // Mesh reference obtained from meshing AIM
  //aimMeshRef *meshRef = NULL;

  //aimStorage *platoInstance;

  //platoInstance = (aimStorage*)instStore;

//  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get mesh
//  meshRef = (aimMeshRef *)aimInputs[Mesh-1].vals.AIMptr;
//  AIM_NOTNULL(meshRef, aimInfo, status);

//cleanup:
  return status;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
	// SU2 Outputs
    /*! \page aimOutputsPLATO AIM Outputs
     * Plato outputs
     */

    #ifdef DEBUG
        printf(" platoAIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
    #endif

    return CAPS_SUCCESS;
}


// Calculate Plato output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, capsValue *val)
{
    int status = CAPS_SUCCESS;

    //aimStorage *platoInstance;

#ifdef DEBUG
    printf(" platoAIM/aimCalcOutput  index = %d!\n", index);
#endif
    //platoInstance = (aimStorage *) instStore;

    val->vals.real = 0.0; // Set default value

    status = CAPS_SUCCESS;
//cleanup:

    return status;
}


void aimCleanup(void *instStore)
{
    aimStorage *platoInstance;

#ifdef DEBUG
    printf(" platoAIM/aimCleanup!\n");
#endif
    platoInstance = (aimStorage *) instStore;

    aim_freeMeshRef(&platoInstance->meshRefObj);
    platoInstance->meshRef = NULL;

    AIM_FREE(platoInstance);
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
    //aimStorage *platoInstance;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;
    capsValue *meshVal;

    // Volume Mesh obtained from meshing AIM
    aimMeshRef *meshRef;

#ifdef DEBUG
    printf(" platoAIM/aimDiscr: tname = %s!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;
    //platoInstance = (aimStorage *) discr->instStore;

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    AIM_STATUS(discr->aInfo, status);

    if (bodies == NULL) {
        AIM_ERROR(discr->aInfo, " platoAIM/aimDiscr: No Bodies!\n");
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

#ifdef DEBUG
    printf(" platoAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
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
    printf(" platoAIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


int aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                      int bIndex, int eIndex, double *bary, int rank,
                      double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" platoAIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


int aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                   int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" platoAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                    int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" platoAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
