/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             metris AIM
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 *     Written by Marshall Galbraith (MIT)
 */

/*!\mainpage Introduction
 *
 * \section overviewMETRIS metris AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with the
 * unstructured mesh adaptation software <a href="https://github.com/LucienRochery/Metris">metris</a>.
 *
 * Currently Metris only supports 2D (Area_Mesh) adaptation, but this will be expeanded in the future.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsMETRIS and \ref aimOutputsMETRIS, respectively.
 *
 * The metris AIM can automatically execute metris, with details provided in \ref aimExecuteMETRIS.
 * The specific executable can be changed with the 'metris' input string.
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"       // Collection of helper functions for miscellaneous analysis

#include "libmeshbWriter.h"
#include "libMeshb/sources/libmeshb7.h"

#define EXPORT_MESHB_VERTEX_ID (1)
#define EXPORT_MESHB_2D_ID (1)
#define EXPORT_MESHB_3D_ID (0)
#define EXPORT_MESHB_VERTEX_3 (10000000)
#define EXPORT_MESHB_VERTEX_4 (200000000)


#ifdef WIN32
/* needs Advapi32.lib & Ws2_32.lib */
#include <Windows.h>
#define strcasecmp stricmp
#define snprintf   _snprintf
#define access     _access
#define F_OK       0
#else
#include <unistd.h>
#endif

//#define DEBUG

#define CROSS(a,b,c)      a[0] = ((b)[1]*(c)[2]) - ((b)[2]*(c)[1]);\
                          a[1] = ((b)[2]*(c)[0]) - ((b)[0]*(c)[2]);\
                          a[2] = ((b)[0]*(c)[1]) - ((b)[1]*(c)[0])
#define DOT(a,b)         ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

enum aimInputs
{
  inMetris = 1,                    /* index is 1-based */
  inPasses,
  inMesh_Format,
  inMesh,
  inMetricFieldFile,
  NUMINPUT = inMetricFieldFile  /* Total number of inputs */
};

enum aimOutputs
{
  outMesh = 1,                    /* index is 1-based */
  outxyz,
  NUMOUT = outxyz                 /* Total number of outputs */
};

static char egadsFileName[] = "metris_in.egads";
static char metris_out_pre[] = "metris_out";
static char metrisInput[] = "metrisInput.txt";
static char metricFileName[] = "metric.solb";


typedef struct {
  int    npts;
  double *xyz;
  double *t;
  int    *ivp;   // volume node index
} edgeData;


typedef struct {
  int    npts;
  double *xyz;
  double *uv;
  int    ntri;
  int    *tris;
  int    *ivp;   // volume node index
} faceData;


typedef struct {
  double **rvec;
  ego    *surfaces;
  ego    body;
  ego    *faces;
  ego    *edges;
  ego    *nodes;
  int    nfaces;
  int    nedges;
  int    nnodes;

  edgeData *tedges;
  faceData *tfaces;
} bodyData;


typedef struct {
  int ivp;         // global index into volume vertexes
  int egadsType;   // egads type, NODE, EDGE, FACE
  int egadsID;     // type-index
  double param[2]; // parametric coordinates of the vertex
} srfVertex;


typedef struct {
  // Mesh reference obtained from meshing AIM
  aimMeshRef *meshRefIn;

  aimMeshRef meshRefOut;

  // Attribute to index map
  mapAttrToIndexStruct groupMap;

  int *edgeID;
  int *faceID;

} aimStorage;


static int
readlibMeshb(void *aimInfo, aimStorage *metrisInstance, aimMesh *mesh);

static int
writeLibMeshbSurface(void *aimInfo, aimMesh *mesh);


static int initiate_bodyData(int numBody, bodyData *bodydata)
{
  int i;

  for (i = 0; i < numBody; i++) {
    bodydata[i].rvec = NULL;
    bodydata[i].surfaces = NULL;
    bodydata[i].faces = NULL;
    bodydata[i].edges = NULL;
    bodydata[i].nodes = NULL;
    bodydata[i].nfaces = 0;
    bodydata[i].nedges = 0;
    bodydata[i].nnodes = 0;
    bodydata[i].tedges = NULL;
    bodydata[i].tfaces = NULL;
  }

  return CAPS_SUCCESS;
}

static int destroy_bodyData(int numBody, bodyData *bodydata)
{
  int i, j;

  if (bodydata == NULL) return CAPS_SUCCESS;

  for (i = 0; i < numBody; i++) {
    for (j = 0; j < bodydata[i].nfaces; j++) {
      if (bodydata[i].surfaces != NULL)
        if (bodydata[i].surfaces[j+bodydata[i].nfaces] != NULL)
          EG_deleteObject(bodydata[i].surfaces[j+bodydata[i].nfaces]);
      if (bodydata[i].rvec != NULL)
        EG_free(bodydata[i].rvec[j]);
    }
    EG_free(bodydata[i].nodes);
    EG_free(bodydata[i].edges);
    EG_free(bodydata[i].faces);
    EG_free(bodydata[i].surfaces);
    EG_free(bodydata[i].rvec);

    if (bodydata[i].tedges != NULL) {
      for (j = 0; j < bodydata[i].nedges; j++) {
        EG_free(bodydata[i].tedges[j].xyz);
        EG_free(bodydata[i].tedges[j].t);
        EG_free(bodydata[i].tedges[j].ivp);
      }
      EG_free(bodydata[i].tedges);
    }

    if (bodydata[i].tfaces != NULL) {
      for (j = 0; j < bodydata[i].nfaces; j++) {
        EG_free(bodydata[i].tfaces[j].xyz);
        EG_free(bodydata[i].tfaces[j].uv);
        EG_free(bodydata[i].tfaces[j].tris);
        EG_free(bodydata[i].tfaces[j].ivp);
      }
      EG_free(bodydata[i].tfaces);
    }
  }

  return CAPS_SUCCESS;
}


static int initiate_aimStorage(aimStorage *metrisInstance)
{

  int status = CAPS_SUCCESS;

  metrisInstance->meshRefIn = NULL;

  // Mesh reference passed to solver
  status = aim_initMeshRef(&metrisInstance->meshRefOut, metrisInstance->meshRefOut.type);
  if (status != CAPS_SUCCESS) return status;

  status = initiate_mapAttrToIndexStruct(&metrisInstance->groupMap);
  if (status != CAPS_SUCCESS) return status;

  metrisInstance->edgeID = NULL;
  metrisInstance->faceID = NULL;

  return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *metrisInstance)
{
  int status; // Function return status

  metrisInstance->meshRefIn = NULL;

  // Free the meshRef
  aim_freeMeshRef(&metrisInstance->meshRefOut);

  status = destroy_mapAttrToIndexStruct(&metrisInstance->groupMap);
  if (status != CAPS_SUCCESS)
    printf("Status = %d, metrisAIM attributeMap group cleanup!!!\n", status);

  AIM_FREE(metrisInstance->edgeID);
  AIM_FREE(metrisInstance->faceID);

  return status;
}


/****************** exposed AIM entry points -- Analysis **********************/

/* aimInitialize: Initialization Information for the AIM */
int
aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
              /*@unused@*/ void **instStore, /*@unused@*/ int *major,
              /*@unused@*/ int *minor, int *nIn, int *nOut,
              int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int        i, status = CAPS_SUCCESS;
  aimStorage *metrisInstance = NULL;

  /* specify the number of analysis inputs  defined in aimInputs
   *     and the number of analysis outputs defined in aimOutputs */
  *nIn    = NUMINPUT;
  *nOut   = NUMOUT;

  /* return if "query" only */
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 0;

  /* specify the name of each field variable */
  *fnames = NULL;

  /* specify the rank of each field variable */
  *franks = NULL;

  /* specify if a field is an input field or output field */
  *fInOut = NULL;

  /* setup our AIM specific state */
  AIM_ALLOC(metrisInstance, 1, aimStorage, aimInfo, status);
  *instStore = metrisInstance;

  status = initiate_aimStorage(metrisInstance);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;

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


// ********************** AIM Function Break *****************************
/* aimInputs: Input Information for the AIM */
int
aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
          char **ainame, capsValue *defval)
{
  int status = CAPS_SUCCESS;

  /*! \page aimInputsMETRIS AIM Inputs
   * The following list outlines the metris inputs along with their default value available
   * through the AIM interface.
   */

  /* fill in the required members based on the index */
  if (index == inMetris) {
    *ainame             = AIM_NAME(metris);
    defval->type        = String;
    defval->lfixed      = Fixed;
    AIM_STRDUP(defval->vals.string, "metris", aimInfo, status);

    /*! \page aimInputsMETRIS
     * - <B>metris = "metris"</B> <br>
     * metris executable
     */

  } else if (index == inPasses) {
    *ainame              = AIM_NAME(Passes);
    defval->type         = Integer;
    defval->dim          = Scalar;
    defval->vals.integer = 30;

    /*! \page aimInputsMETRIS
     * - <B> Passes = 30</B> <br>
     * Number of metris internal adaptation iterations
     */

  } else if (index == inMesh) {
    *ainame             = AIM_NAME(Mesh);
    defval->type        = PointerMesh;
    defval->nrow        = 1;
    defval->lfixed      = Fixed;
    defval->vals.AIMptr = NULL;
    defval->nullVal     = IsNull;
    AIM_STRDUP(defval->meshWriter, MESHWRITER, aimInfo, status);

    /*! \page aimInputsMETRIS
     * - <B>Mesh = NULL</B> <br>
     * An Area_Mesh link for mesh adaptation
     */

  } else if (index == inMesh_Format) {
    *ainame               = EG_strdup("Mesh_Format");
    defval->type          = String;
    defval->vals.string   = NULL;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;
    defval->lfixed        = Change;

    /*! \page aimInputsMETRIS
     * \include{doc} Mesh_Format.dox
     */

  } else if (index == inMetricFieldFile) {
    *ainame           = AIM_NAME(MetricFieldFile);
    defval->type      = String;
    defval->lfixed    = Fixed;
    defval->dim       = Scalar;
    defval->nullVal   = IsNull;

    /*! \page aimInputsMETRIS
     * - <B> MetricFieldFile = NULL</B> <br>
     * Metric field file in libMeshb format. <br>
     */

  } else {
      status = CAPS_BADINDEX;
      AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimUpdateState: The always the first call in the execution sequence */
int
aimUpdateState(void *instStore, void *aimInfo, capsValue *inputs)
{
  int        nBody, status = CAPS_SUCCESS;
  const char *intents;
  ego        *bodies;

  char aimFile[PATH_MAX], metris_out[PATH_MAX];

  const char *groupName = NULL;
  aimStorage *metrisInstance;

  int nEdges, nFaces, iedge, iface, cID;
  ego *edges=NULL, *faces=NULL;

  int oclass, mtype, faceSense, ncild, *senses;
  double limits[4];
  ego surface, geom, *childs;

  int *header=NULL;
  double *gdata=NULL;
  double norm[3];

  metrisInstance = (aimStorage *) instStore;
  AIM_NOTNULL(inputs, aimInfo, status);

  snprintf(metris_out, PATH_MAX, "%s%s", metris_out_pre, MESHEXTENSION);

  if (inputs[inMesh-1].nullVal == IsNull &&
      aim_isFile(aimInfo, metris_out) != CAPS_SUCCESS) {
    AIM_ANALYSISIN_ERROR(aimInfo, inMesh, "'Mesh' input must be linked to generate the initial mesh!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if ( inputs[inMesh-1].nullVal == IsNull &&
       inputs[inMetricFieldFile-1].nullVal == IsNull ) {
    AIM_ANALYSISIN_ERROR(aimInfo, inMetricFieldFile, "'MetricFieldFile' input must be specified!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  AIM_STATUS(aimInfo, status);
  if (nBody != 1) {
    AIM_ERROR(aimInfo, "metris only supports a single body: numBody = %d", nBody);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  AIM_NOTNULL(bodies, aimInfo, status);

  if (inputs[inMesh-1].nullVal == NotNull) {
    // Get mesh
    metrisInstance->meshRefIn = (aimMeshRef *) inputs[inMesh-1].vals.AIMptr;

    // Get attribute to index mapping
    status = create_MeshRefToIndexMap(aimInfo, metrisInstance->meshRefIn, &metrisInstance->groupMap);
    AIM_STATUS(aimInfo, status);
  } else {
    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    if (metrisInstance->groupMap.numAttribute == 0) {
      status = create_CAPSGroupAttrToIndexMap(nBody,
                                              bodies,
                                              2, // Body, Faces, and Edges
                                              &metrisInstance->groupMap);
      AIM_STATUS(aimInfo, status);
    }
    metrisInstance->meshRefIn = NULL;
  }

  // clear the previous mesh
  aim_freeMeshRef(&metrisInstance->meshRefOut);

  status = aim_file(aimInfo, metris_out_pre, aimFile);
  AIM_STATUS(aimInfo, status);
  AIM_STRDUP(metrisInstance->meshRefOut.fileName, aimFile, aimInfo, status);

  status = EG_getBodyTopos(bodies[0], NULL, EDGE, &nEdges, &edges);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(edges, aimInfo, status);

  AIM_FREE(metrisInstance->edgeID);
  AIM_ALLOC(metrisInstance->edgeID, nEdges, int, aimInfo, status);

  for (iedge = 0; iedge < nEdges; iedge++) {

    // Look for component/boundary ID for attribute mapper based on capsGroup
    status = retrieve_CAPSGroupAttr(edges[iedge], &groupName);
    if (status == EGADS_NOTFOUND) {
      metrisInstance->edgeID[iedge] = iedge;
      continue;
    }
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "No capsGroup attribute found on Edge %d, unable to assign a boundary index value",
                iedge+1);
      print_AllAttr( aimInfo, edges[iedge] );
      goto cleanup;
    }

    /*@-nullpass@*/
    status = get_mapAttrToIndexIndex(&metrisInstance->groupMap, groupName, &cID);
    AIM_STATUS(aimInfo, status, "Unable to retrieve index from capsGroup: %s",
               groupName);
    /*@+nullpass@*/

    metrisInstance->edgeID[iedge] = cID;
  }

  status = EG_getBodyTopos(bodies[0], NULL, FACE, &nFaces, &faces);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(faces, aimInfo, status);

  AIM_FREE(metrisInstance->faceID);
  AIM_ALLOC(metrisInstance->faceID, nFaces, int, aimInfo, status);

  if (nFaces == 1) {

    status = EG_getTopology(faces[0], &surface, &oclass, &faceSense,
                            limits, &ncild, &childs, &senses);
    AIM_STATUS(aimInfo, status);

    status = EG_getGeometry(surface, &oclass, &mtype, &geom, &header, &gdata);
    AIM_STATUS(aimInfo, status);

    if (mtype == PLANE) {
      CROSS(norm, gdata+3, gdata+6);
      norm[2] *= faceSense;
      if (fabs(norm[2] - 1) > 1e-7) {
        AIM_ERROR(aimInfo, "metris requires 2D meshes in the x-y plane with the normal in positive z!");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    metrisInstance->faceID[0] = 1;

  } else {

    for (iface = 0; iface < nFaces; iface++) {

      // Look for component/boundary ID for attribute mapper based on capsGroup
      status = retrieve_CAPSGroupAttr(faces[iface], &groupName);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup attribute found on Face %d, unable to assign a boundary index value",
                  iface+1);
        print_AllAttr( aimInfo, faces[iface] );
        goto cleanup;
      }

      /*@-nullpass@*/
      status = get_mapAttrToIndexIndex(&metrisInstance->groupMap, groupName, &cID);
      AIM_STATUS(aimInfo, status, "Unable to retrieve boundary index from capsGroup: %s",
                 groupName);
      /*@+nullpass@*/

      metrisInstance->faceID[iface] = cID;
    }
  }

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(edges);
  AIM_FREE(faces);
  AIM_FREE(header);
  AIM_FREE(gdata);
  return status;
}


// ********************** AIM Function Break *****************************
/* aimPreAnalysis: Parse Inputs, Generate Input File(s) */
int
aimPreAnalysis(/*@unused@*/ const void *instStore, void *aimInfo,
               capsValue *inputs)
{
  int status = CAPS_SUCCESS;
  int i;

  int nBody=0;
  const char *intents;
  ego *bodies=NULL;

  ego *bodyCopy=NULL, context=NULL, model=NULL;

  // Output filename
  char metris_in[PATH_MAX];
  char metris_out[PATH_MAX];
  char aimEgadsFile[PATH_MAX];
  char aimFile[PATH_MAX];
  char relPath[PATH_MAX];
  FILE *fp = NULL;

  const aimStorage *metrisInstance;

  AIM_NOTNULL(inputs, aimInfo, status);
  metrisInstance = (const aimStorage *) instStore;

  snprintf(metris_out, PATH_MAX, "%s%s", metris_out_pre, MESHEXTENSION);

  if (metrisInstance->meshRefIn != NULL) {

    /* create a symbolic link to the file name*/
    snprintf(metris_in, PATH_MAX, "%s%s", metrisInstance->meshRefIn->fileName, MESHEXTENSION);
    status = aim_symLink(aimInfo, metris_in, metris_out);
    AIM_STATUS(aimInfo, status);

  } else {

    if (aim_isFile(aimInfo, metris_out) != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "'%s' does not exist!", metris_out);
      status = CAPS_IOERR;
      goto cleanup;
    }

    if (inputs[inMetricFieldFile-1].nullVal == NotNull) {

      if (access(inputs[inMetricFieldFile-1].vals.string, F_OK) != 0) {
        AIM_ERROR(aimInfo, "'%s' does not exist!", inputs[inMetricFieldFile-1].vals.string);
        status = CAPS_IOERR;
        goto cleanup;
      }

      snprintf(metris_in, PATH_MAX, "metris_in%s", MESHEXTENSION);

      // copy over metris_out to metris_in
      status = aim_file(aimInfo, metris_out, aimFile);
      AIM_STATUS(aimInfo, status);
      status = aim_cpFile(aimInfo, aimFile, metris_in);
      AIM_STATUS(aimInfo, status);

      // remove metris_out
      status = aim_rmFile(aimInfo, metris_out);
      AIM_STATUS(aimInfo, status);

      // get the relative path
      status = aim_relPath(aimInfo, inputs[inMetricFieldFile-1].vals.string, metricFileName, relPath);
      AIM_STATUS(aimInfo, status);

      if (strncmp(relPath, metricFileName, PATH_MAX) != 0) {
        // Simply create a link to the file
        status = aim_symLink(aimInfo, inputs[inMetricFieldFile-1].vals.string, metricFileName);
        AIM_STATUS(aimInfo, status);
      }

    } else {
      AIM_ERROR(aimInfo, "Developer error!");
      status = CAPS_NOTIMPLEMENT;
      goto cleanup;
    }

    // remove previous meshes after renaming metris_out to metris_in
    status = aim_deleteMeshes(aimInfo, &metrisInstance->meshRefOut);
    AIM_STATUS(aimInfo, status);

    status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    AIM_ALLOC(bodyCopy, nBody, ego, aimInfo, status);
    for (i = 0; i < nBody; i++) bodyCopy[i] = NULL;

    // Get context
    status = EG_getContext(bodies[0], &context);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(context, aimInfo, status);

    // Make a copy of the bodies
    for (i = 0; i < nBody; i++) {
      status = EG_copyObject(bodies[i], NULL, &bodyCopy[i]);
      AIM_STATUS(aimInfo, status);
    }

    // Create a model from the copied bodies
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, nBody, bodyCopy,
                             NULL, &model);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(model, aimInfo, status);

    status = aim_file(aimInfo, egadsFileName, aimEgadsFile);
    AIM_STATUS(aimInfo, status);

    //printf("Writing egads file '%s'....\n", aimEgadsFile);
    remove(aimEgadsFile);
    status = EG_saveModel(model, aimEgadsFile);
    AIM_STATUS(aimInfo, status);


    status = aim_mkDir(aimInfo, "tmp");
    AIM_STATUS(aimInfo, status);

    fp = aim_fopen(aimInfo, metrisInput, "w");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open %s", metrisInput);
      status = CAPS_IOERR;
      goto cleanup;
    }

    fprintf(fp,
            " -in %s -met %s -cad %s -adapt %d --refine-conventions-inp --refine-conventions-out -adp-opt-niter -1 -verb 1 -prefix tmp/ -out %s",
            metris_in,
            metricFileName,
            egadsFileName,
            inputs[inPasses-1].vals.integer,
            metris_out);

    fclose(fp); fp = NULL;
  }

  status = CAPS_SUCCESS;

cleanup:

  if (model != NULL) {
    EG_deleteObject(model);
  } else {
    if (bodyCopy != NULL) {
      for (i = 0;  i < nBody; i++) {
        if (bodyCopy[i] != NULL)  {
          (void) EG_deleteObject(bodyCopy[i]);
        }
      }
    }
  }
  AIM_FREE(bodyCopy);
  if (fp != NULL) fclose(fp);

  return status;
}


// ********************** AIM Function Break *****************************
/* aimExecute: runs the Analysis & specifies the AIM does the execution */
int
aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimInfo,
           int *state)
{
  /*! \page aimExecuteMETRIS AIM Execution
   *
   * If auto execution is enabled when creating an metris AIM,
   * the AIM will execute metris just-in-time with the command line:
   *
   * \code{.sh}
   * metris $(cat metrisInput.txt) > metrisOutput.txt
   * \endcode
   *
   * where preAnalysis generated the file "metrisInput.txt" which contains commandline arguments for metris.
   *
   * The metris analysis directory is assumed to contain a metric.meshb file. This file will
   * be generated automatically with preAnalysis using ScalarFieldFile or HessianFieldFile inputs, or
   * can be generated manually via system calls to metris and setting MetricFieldFile.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an metris AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or metris can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * metris.preAnalysis()
   *
   * print ("\n\nRunning......")
   * metris.system("metris $(cat metrisInput.txt) > metrisOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * metris.postAnalysis()
   * \endcode
   */

  int status = CAPS_SUCCESS;
  char command[PATH_MAX];
  capsValue *metris;

  const aimStorage *metrisInstance;

  metrisInstance = (const aimStorage *) instStore;

  *state = 0;

  if (metrisInstance->meshRefIn != NULL)
    return CAPS_SUCCESS;

  status = aim_getValue(aimInfo, inMetris, ANALYSISIN, &metris);
  AIM_STATUS(aimInfo, status);

  // execute metris adapt
  snprintf(command, PATH_MAX,
           "%s $(cat %s) > metrisOutput.txt",
           metris->vals.string, metrisInput);

  status = aim_system(aimInfo, NULL, command);
  AIM_STATUS(aimInfo, status, "Failed to execute: %s", command);

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimPostAnalysis: Perform any processing after the Analysis is run */
int
aimPostAnalysis(void *instStore, void *aimInfo,
                /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;
  aimMesh    mesh;
  aimStorage *metrisInstance;

  metrisInstance = (aimStorage *) instStore;

  /*@-immediatetrans@*/
  mesh.meshData = NULL;
  mesh.meshRef = &metrisInstance->meshRefOut;
  /*@+immediatetrans@*/

  //Only read the surface tessellation
  status = readlibMeshb(aimInfo, metrisInstance, &mesh);
  AIM_STATUS(aimInfo, status);

  if (mesh.meshRef->type == aimSurfaceMesh) {
    status = writeLibMeshbSurface(aimInfo, &mesh);
    AIM_STATUS(aimInfo, status);
  }

  status = aim_freeMeshData(mesh.meshData);
  AIM_STATUS(aimInfo, status);
  AIM_FREE(mesh.meshData);


  /* Explicitly write out any requested meshes */
  status = aim_queryMeshes( aimInfo, inMesh_Format, ANALYSISIN, &metrisInstance->meshRefOut );
  if (status > 0) {
    /*@-immediatetrans@*/
    mesh.meshData = NULL;
    mesh.meshRef = &metrisInstance->meshRefOut;
    /*@+immediatetrans@*/

    status = readlibMeshb(aimInfo, metrisInstance, &mesh);
    AIM_STATUS(aimInfo, status);

    status = aim_writeMeshes(aimInfo, inMesh_Format, ANALYSISIN, &mesh);
    AIM_STATUS(aimInfo, status);

    status = aim_freeMeshData(mesh.meshData);
    AIM_STATUS(aimInfo, status);
    AIM_FREE(mesh.meshData);
  }
  else
    AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimOutputs: Output Information for the AIM */
int
aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
           /*@unused@*/ int index, char **aoname, capsValue *form)
{
  int status = CAPS_SUCCESS;
#ifdef DEBUG
  printf(" metrisAIM/aimOutputs instance = %d  index  = %d!\n",
         aim_getInstance(aimInfo), index);
#endif

  /*! \page aimOutputsMETRIS AIM Outputs
   * List of available outputs from the metris AIM
   */

  if (index == outMesh) {

      *aoname           = AIM_NAME(Mesh);
      form->type        = PointerMesh;
      form->dim         = Scalar;
      form->lfixed      = Fixed;
      form->sfixed      = Fixed;
      form->vals.AIMptr = NULL;
      form->nullVal     = IsNull;

      /*! \page aimOutputsMETRIS
       * - <B> Mesh </B> <br>
       * The output Area_Mesh or Volume_Mesh for a link
       */

  } else if (index == outxyz) {

      *aoname           = AIM_NAME(xyz);
      form->type        = Double;
      form->dim         = Array2D;

      /*! \page aimOutputsMETRIS
       * - <B> xyz </B> <br>
       * Grid coordinates. Useful for constructing scalar, hessian, or metric fields
       */

  } else {
      status = CAPS_BADINDEX;
      AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);
  }

  AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimCalcOutput: Calculate/Retrieve Output Information */
int
aimCalcOutput(void *instStore, void *aimInfo,
              int index, capsValue *val)
{
  int        i, j, status = CAPS_SUCCESS;
  aimStorage *metrisInstance;
  aimMesh    mesh;

  metrisInstance = (aimStorage *) instStore;

   if (outMesh == index) {

     status = aim_queryMeshes( aimInfo, outMesh, ANALYSISOUT, &metrisInstance->meshRefOut );
     if (status > 0) {
       /*@-immediatetrans@*/
       mesh.meshData = NULL;
       mesh.meshRef = &metrisInstance->meshRefOut;
       /*@+immediatetrans@*/

       status = readlibMeshb(aimInfo, metrisInstance, &mesh);
       AIM_STATUS(aimInfo, status);

       status = aim_writeMeshes(aimInfo, outMesh, ANALYSISOUT, &mesh);
       AIM_STATUS(aimInfo, status);

       status = aim_freeMeshData(mesh.meshData);
       AIM_STATUS(aimInfo, status);
       AIM_FREE(mesh.meshData);
     }
     else
       AIM_STATUS(aimInfo, status);

     /*@-immediatetrans@*/
     // Return the volume mesh references
     val->nrow        = 1;
     val->vals.AIMptr = &metrisInstance->meshRefOut;
     /*@+immediatetrans@*/

   } else if (index == outxyz) {

     /*@-immediatetrans@*/
     mesh.meshData = NULL;
     mesh.meshRef = &metrisInstance->meshRefOut;
     /*@+immediatetrans@*/

     status = readlibMeshb(aimInfo, metrisInstance, &mesh);
     AIM_STATUS(aimInfo, status);
     AIM_NOTNULL(mesh.meshData, aimInfo, status);

     AIM_ALLOC(val->vals.reals, mesh.meshData->dim*mesh.meshData->nVertex, double, aimInfo, status);
     val->nrow = mesh.meshData->nVertex;
     val->ncol = mesh.meshData->dim;

     for (i = 0; i < mesh.meshData->nVertex; i++) {
       for (j = 0; j < mesh.meshData->dim; j++) {
         val->vals.reals[mesh.meshData->dim*i+j] = mesh.meshData->verts[i][j];
       }
     }

     status = aim_freeMeshData(mesh.meshData);
     AIM_STATUS(aimInfo, status);
     AIM_FREE(mesh.meshData);

   } else {

     status = CAPS_BADINDEX;
     AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

   }

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimCleanup: Free up the AIMs storage */
void aimCleanup(/*@unused@*/ void *instStore)
{
  aimStorage *metrisInstance = NULL;

  /* clean up any allocated data */

  metrisInstance = (aimStorage *) instStore;
  if (metrisInstance == NULL) return;

  destroy_aimStorage(metrisInstance);
  AIM_FREE(metrisInstance);
}


// ********************** AIM Function Break *****************************
static void swapd(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

static void swapi(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}


// ********************** AIM Function Break *****************************
// A function to implement bubble sort
static void
bubbleSortEdge(edgeData *tedge)
{

  int i, j;
  for (i = 0; i < tedge->npts-1; i++)
    // Last i elements are already in place
    for (j = 0; j < tedge->npts-i-1; j++)
      if (tedge->t[j] > tedge->t[j+1]) {
        swapd(&tedge->t[j]      , &tedge->t[j+1]        );
        swapd(&tedge->xyz[3*j+0], &tedge->xyz[3*(j+1)+0]);
        swapd(&tedge->xyz[3*j+1], &tedge->xyz[3*(j+1)+1]);
        swapd(&tedge->xyz[3*j+2], &tedge->xyz[3*(j+1)+2]);
        swapi(&tedge->ivp[j]    , &tedge->ivp[j+1]      );
      }
}


static void
bubbleSortFace(faceData *tface)
{

  int i, j;
  for (i = 0; i < tface->npts-1; i++)
    // Last i elements are already in place
    for (j = 0; j < tface->npts-i-1; j++)
      if (tface->ivp[j] > tface->ivp[j+1]) {
        swapd(&tface->uv[2*j+0] , &tface->uv[2*(j+1)+0] );
        swapd(&tface->uv[2*j+1] , &tface->uv[2*(j+1)+1] );
        swapd(&tface->xyz[3*j+0], &tface->xyz[3*(j+1)+0]);
        swapd(&tface->xyz[3*j+1], &tface->xyz[3*(j+1)+1]);
        swapd(&tface->xyz[3*j+2], &tface->xyz[3*(j+1)+2]);
        swapi(&tface->ivp[j]    , &tface->ivp[j+1]      );
      }
}

// use bisection to find the face Index
static int
faceIndex(const int ivp, faceData *tface)
{
  int i0 = 0;
  int i1 = tface->npts/2;
  int i2 = tface->npts;

  while(tface->ivp[i1] != ivp) {
    if (ivp > tface->ivp[i1]) {
      i0 = i1;
      i1 = (i1 + i2)/2;
    } else {
      i2 = i1;
      i1 = (i0 + i1)/2;
    }
  }

  return i1+1;
}

// ********************** AIM Function Break *****************************
static int
readlibMeshb(void *aimInfo, aimStorage *metrisInstance, aimMesh *mesh)
{
  int    status = CAPS_SUCCESS;

  char   attrname[128];
  int    nBody=0, nLine=0, nTri=0, nTet=0;
  int    i, j, elementIndex, nPoint, igroup, iglobal, nglobal, localIndex, topoIndex;
  int    elem[4], ivp, id, npts, ntri, iedge, iface;
  int    meshVersion, nEdgeVerts, nFaceVerts, *faceGroups=NULL, *edgeGroups=NULL, tetGroup;
  int    oclass, mtype, *faceVertID=NULL, *face_tris;
  double reference, t, uv[2], double_gref;
  double result[18], params[3];
  double     *face_uv = NULL, *face_xyz = NULL;
  double     v1[3], v2[3], faceNormal[3], triNormal[3], ndot;
  const int    *tris = NULL, *tric = NULL, *ptype = NULL, *pindex = NULL;
  const double *pxyz = NULL, *puv = NULL;
  const char   *groupName = NULL;
  enum aimMeshElem elementTopo;
  char filename[PATH_MAX];
  int64_t fileID=0;
  aimMeshData *meshData = NULL;
  bodyData bodydata;
  const char *intents;
  ego        *bodies, body, tess, ref, prev, next;

  if (mesh           == NULL) return CAPS_NULLOBJ;
  if (mesh->meshRef  == NULL) return CAPS_NULLOBJ;
  if (mesh->meshRef->fileName  == NULL) return CAPS_NULLOBJ;

  status = initiate_bodyData(1, &bodydata);
  AIM_STATUS(aimInfo, status);

  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  AIM_STATUS(aimInfo, status);

  body = bodies[0];

  status = EG_getBodyTopos(body, NULL, NODE, &bodydata.nnodes, NULL);
  AIM_STATUS(aimInfo, status);

  status = EG_getBodyTopos(body, NULL, EDGE, &bodydata.nedges, &bodydata.edges);
  AIM_STATUS(aimInfo, status);

  status = EG_getBodyTopos(body, NULL, FACE, &bodydata.nfaces, &bodydata.faces);
  AIM_STATUS(aimInfo, status);

  status = aim_freeMeshData(mesh->meshData);
  AIM_STATUS(aimInfo, status);
  AIM_FREE(mesh->meshData);

  AIM_ALLOC(meshData, 1, aimMeshData, aimInfo, status);
  status = aim_initMeshData(meshData);
  AIM_STATUS(aimInfo, status);

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  fileID = GmfOpenMesh(filename, GmfRead, &meshVersion, &meshData->dim);

  if (fileID == 0) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  AIM_ALLOC(edgeGroups, bodydata.nedges, int, aimInfo, status);
  AIM_ALLOC(faceGroups, bodydata.nfaces, int, aimInfo, status);
  for (i = 0; i < bodydata.nedges; i++) edgeGroups[i] = -1;
  for (i = 0; i < bodydata.nfaces; i++) faceGroups[i] = -1;
  tetGroup = -1;

  meshData->nVertex = GmfStatKwd(fileID, GmfVertices);
  AIM_ALLOC(meshData->verts, meshData->nVertex, aimMeshCoords, aimInfo, status);

  status = GmfGotoKwd(fileID, GmfVertices);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nLine = GmfStatKwd(fileID, GmfEdges);
  nTri  = GmfStatKwd(fileID, GmfTriangles);
  nTet  = GmfStatKwd(fileID, GmfTetrahedra);

  // Real nodal coordinates
  if (meshData->dim == 2) {
    for (i = 0; i < meshData->nVertex; i++) {
      status = GmfGetLin(fileID, GmfVertices, &meshData->verts[i][0],
                                              &meshData->verts[i][1], &reference);
      meshData->verts[i][2] = 0;
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

    meshData->nTotalElems = nLine + nTri;

    mesh->meshRef->type = aimAreaMesh;

  } else {
    for (i = 0; i < meshData->nVertex; i++) {
      status = GmfGetLin(fileID, GmfVertices, &meshData->verts[i][0],
                                              &meshData->verts[i][1],
                                              &meshData->verts[i][2], &reference);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

    meshData->nTotalElems =  nTri + nTet;

    if (nTet > 0)
      mesh->meshRef->type = aimVolumeMesh;
    else
      mesh->meshRef->type = aimSurfaceMesh;
  }


  // allocate the element map that maps back to the original element numbering
  AIM_ALLOC(meshData->elemMap, meshData->nTotalElems, aimMeshIndices, aimInfo, status);

  /* Start of element index */
  elementIndex = 0;

  if (meshData->dim == 2 || mesh->meshRef->type == aimSurfaceMesh) {

    // Elements Line
    status = GmfGotoKwd(fileID, GmfEdges);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    nPoint = 2;
    elementTopo = aimLine;
    for (i = 0; i < nLine; i++) {

      /* read the element and group */
      status = GmfGetLin(fileID, GmfEdges, &elem[0],
                                           &elem[1], &igroup);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      if (igroup <= 0 || igroup > bodydata.nedges) {
        AIM_ERROR(aimInfo, "Edge ID %d is out of range [1, %d]!=", igroup, bodydata.nedges);
        status = CAPS_IOERR;
        goto cleanup;
      }
      igroup -= 1; // make zero based

      igroup = metrisInstance->edgeID[igroup]-1;

      /* add the group if necessary */
      if (edgeGroups[igroup] == -1) {
        status = get_mapAttrToIndexKeyword(&metrisInstance->groupMap, igroup+1, &groupName);
        AIM_NOTFOUND(aimInfo, status);

        status = aim_addMeshElemGroup(aimInfo, groupName, igroup+1, elementTopo, 1, nPoint, meshData);
        AIM_STATUS(aimInfo, status);
        edgeGroups[igroup] = meshData->nElemGroup-1;
      }
      igroup = edgeGroups[igroup];

      /* add the element to the group */
      status = aim_addMeshElem(aimInfo, 1, &meshData->elemGroups[igroup]);
      AIM_STATUS(aimInfo, status);

      /* set the element connectivity */
      for (j = 0; j < nPoint; j++)
        meshData->elemGroups[igroup].elements[nPoint*(meshData->elemGroups[igroup].nElems-1) + j] = elem[j];

      meshData->elemMap[elementIndex][0] = igroup;
      meshData->elemMap[elementIndex][1] = meshData->elemGroups[igroup].nElems-1;

      elementIndex += 1;
    }
  }

  AIM_ALLOC(bodydata.tfaces, bodydata.nfaces, faceData, aimInfo, status);
  for (j = 0; j < bodydata.nfaces; j++) {
    bodydata.tfaces[j].npts  = 0;
    bodydata.tfaces[j].xyz   = NULL;
    bodydata.tfaces[j].uv    = NULL;
    bodydata.tfaces[j].ntri  = 0;
    bodydata.tfaces[j].tris  = NULL;
    bodydata.tfaces[j].ivp   = NULL;
  }

  /* Elements triangles */
  status = GmfGotoKwd(fileID, GmfTriangles);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nPoint = 3;
  elementTopo = aimTri;
  for (i = 0; i < nTri; i++) {

    /* read the element and group */
    status = GmfGetLin(fileID, GmfTriangles, &elem[0],
                                             &elem[1],
                                             &elem[2], &igroup);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    if (igroup <= 0) {
      AIM_ERROR(aimInfo, "Group must be a positive number: %d!", igroup);
      status = CAPS_IOERR;
      goto cleanup;
    }
    igroup -= 1; // make zero based

    if (mesh->meshRef->type != aimAreaMesh) {
      ntri = bodydata.tfaces[igroup].ntri;
      AIM_REALL(bodydata.tfaces[igroup].tris, 3*(ntri+1), int, aimInfo, status);
      bodydata.tfaces[igroup].tris[3*ntri+0] = elem[0];
      bodydata.tfaces[igroup].tris[3*ntri+1] = elem[1];
      bodydata.tfaces[igroup].tris[3*ntri+2] = elem[2];
      bodydata.tfaces[igroup].ntri++;

      igroup = metrisInstance->faceID[igroup]-1;
    }

    /* add the group if necessary */
    if (faceGroups[igroup] == -1) {
      status = get_mapAttrToIndexKeyword(&metrisInstance->groupMap, igroup+1, &groupName);
      if (status != CAPS_SUCCESS && bodydata.nfaces == 1)
        groupName = NULL;
      else
        AIM_STATUS(aimInfo, status);

      status = aim_addMeshElemGroup(aimInfo, groupName, igroup+1, elementTopo, 1, nPoint, meshData);
      AIM_STATUS(aimInfo, status);
      faceGroups[igroup] = meshData->nElemGroup-1;
    }
    igroup = faceGroups[igroup];

    /* add the element to the group */
    status = aim_addMeshElem(aimInfo, 1, &meshData->elemGroups[igroup]);
    AIM_STATUS(aimInfo, status);

    /* set the element connectivity */
    for (j = 0; j < nPoint; j++)
      meshData->elemGroups[igroup].elements[nPoint*(meshData->elemGroups[igroup].nElems-1) + j] = elem[j];

    meshData->elemMap[elementIndex][0] = igroup;
    meshData->elemMap[elementIndex][1] = meshData->elemGroups[igroup].nElems-1;

    elementIndex += 1;
  }

  // If the surface mesh has been processed, read Tets
  if (mesh->meshRef->maps != NULL) {

    if (mesh->meshRef->type == aimVolumeMesh) {

      // Elements Tetrahedral
      status = GmfGotoKwd(fileID, GmfTetrahedra);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      igroup = 1;
      nPoint = 4;
      elementTopo = aimTet;

      /* add the group tetGroup */
      if (tetGroup == -1) {
        status = aim_addMeshElemGroup(aimInfo, NULL, igroup, elementTopo, 1, nPoint, meshData);
        AIM_STATUS(aimInfo, status);
        tetGroup = meshData->nElemGroup-1;
      }
      igroup = tetGroup;

      /* add the element to the group */
      status = aim_addMeshElem(aimInfo, nTet, &meshData->elemGroups[igroup]);
      AIM_STATUS(aimInfo, status);

      for (i = 0; i < nTet; i++) {

        /* read the element and group */
        status = GmfGetLin(fileID, GmfTetrahedra, &elem[0],
                                                  &elem[1],
                                                  &elem[2],
                                                  &elem[3], &igroup);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

        if (igroup != 0) {
          AIM_ERROR(aimInfo, "Tetrahedra group must be 0: %d!", igroup);
          status = CAPS_IOERR;
          goto cleanup;
        }

        igroup = tetGroup;

        /* set the element connectivity */
        for (j = 0; j < nPoint; j++)
          meshData->elemGroups[igroup].elements[nPoint*i + j] = elem[j];

        meshData->elemMap[elementIndex][0] = igroup;
        meshData->elemMap[elementIndex][1] = i;

        elementIndex += 1;
      }
    }

  } else {
    // generate tessellation

    // read parametric coordinates

    AIM_ALLOC(bodydata.tedges, bodydata.nedges, edgeData, aimInfo, status);
    for (j = 0; j < bodydata.nedges; j++) {
      bodydata.tedges[j].npts = 0;
      bodydata.tedges[j].xyz  = NULL;
      bodydata.tedges[j].t    = NULL;
      bodydata.tedges[j].ivp  = NULL;
    }

    // Read EDGEs
    nEdgeVerts = GmfStatKwd(fileID, GmfVerticesOnGeometricEdges);
    status = GmfGotoKwd(fileID, GmfVerticesOnGeometricEdges);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    // first count points on each edge
    for (j = 0; j < nEdgeVerts; j++) {
      status = GmfGetLin(fileID, GmfVerticesOnGeometricEdges,
                         &ivp,
                         &id,
                         &t,
                         &double_gref);  // refine abuse of dist
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      if (id <= 0 || id > bodydata.nedges) {
        AIM_ERROR(aimInfo, "Edge ID %d is out of range [1, %d]", id, bodydata.nedges);
        status = CAPS_IOERR;
        goto cleanup;
      }
      bodydata.tedges[id-1].npts++;
    }

    for (j = 0; j < bodydata.nedges; j++) {
      npts = bodydata.tedges[j].npts;
      AIM_ALLOC(bodydata.tedges[j].xyz, 3*npts, double, aimInfo, status);
      AIM_ALLOC(bodydata.tedges[j].t  ,   npts, double, aimInfo, status);
      AIM_ALLOC(bodydata.tedges[j].ivp,   npts, int   , aimInfo, status);
      bodydata.tedges[j].npts = 0;
    }

    status = GmfGotoKwd(fileID, GmfVerticesOnGeometricEdges);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    // read the data
    for (j = 0; j < nEdgeVerts; j++) {
      status = GmfGetLin(fileID, GmfVerticesOnGeometricEdges,
                         &ivp,
                         &id,
                         &t,
                         &double_gref);  // refine abuse of dist
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      id = id-1;

      npts = bodydata.tedges[id].npts;

      bodydata.tedges[id].t[npts] = t;

      bodydata.tedges[id].xyz[3*npts+0] = meshData->verts[ivp-1][0];
      bodydata.tedges[id].xyz[3*npts+1] = meshData->verts[ivp-1][1];
      bodydata.tedges[id].xyz[3*npts+2] = meshData->verts[ivp-1][2];

      bodydata.tedges[id].ivp[npts] = ivp;

      bodydata.tedges[id].npts++;
    }

    for (j = 0; j < bodydata.nedges; j++) {
      bubbleSortEdge(&bodydata.tedges[j]);
    }

    if (mesh->meshRef->type != aimAreaMesh) {
      // Count face points first
      nFaceVerts = GmfStatKwd(fileID, GmfVerticesOnGeometricTriangles);
      status = GmfGotoKwd(fileID, GmfVerticesOnGeometricTriangles);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      for (j = 0; j < nFaceVerts; j++) {

        status = GmfGetLin(fileID, GmfVerticesOnGeometricTriangles,
                           &ivp,
                           &id,
                           &uv[0], &uv[1],
                           &double_gref); // refine abuse of dist
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

        if (id <= 0 || id > bodydata.nfaces) {
          AIM_ERROR(aimInfo, "Face ID %d is out of range [1, %d]", id, bodydata.nfaces);
          status = CAPS_IOERR;
          goto cleanup;
        }
        id = id-1;
        bodydata.tfaces[id].npts++;
      }

      for (j = 0; j < bodydata.nfaces; j++) {
        npts = bodydata.tfaces[j].npts;
        AIM_ALLOC(bodydata.tfaces[j].xyz , 3*npts, double, aimInfo, status);
        AIM_ALLOC(bodydata.tfaces[j].uv  , 2*npts, double, aimInfo, status);
        AIM_ALLOC(bodydata.tfaces[j].ivp ,   npts, int   , aimInfo, status);
        bodydata.tfaces[j].npts = 0;
      }

      status = GmfGotoKwd(fileID, GmfVerticesOnGeometricTriangles);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      // read the data
      for (j = 0; j < nFaceVerts; j++) {
        status = GmfGetLin(fileID, GmfVerticesOnGeometricTriangles,
                           &ivp,
                           &id,
                           &uv[0], &uv[1],
                           &double_gref); // refine abuse of dist
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
        id = id-1;

        npts = bodydata.tfaces[id].npts;

        bodydata.tfaces[id].uv[2*npts+0] = uv[0];
        bodydata.tfaces[id].uv[2*npts+1] = uv[1];

        bodydata.tfaces[id].xyz[3*npts+0] = meshData->verts[ivp-1][0];
        bodydata.tfaces[id].xyz[3*npts+1] = meshData->verts[ivp-1][1];
        bodydata.tfaces[id].xyz[3*npts+2] = meshData->verts[ivp-1][2];

        bodydata.tfaces[id].ivp[npts] = ivp;

        bodydata.tfaces[id].npts++;
      }

      for (j = 0; j < bodydata.nfaces; j++) {
        bubbleSortFace(&bodydata.tfaces[j]);
        // get the face triangulation
        for (i = 0; i < bodydata.tfaces[j].ntri; i++) {
          bodydata.tfaces[j].tris[3*i+0] = faceIndex(bodydata.tfaces[j].tris[3*i+0], &bodydata.tfaces[j]);
          bodydata.tfaces[j].tris[3*i+1] = faceIndex(bodydata.tfaces[j].tris[3*i+1], &bodydata.tfaces[j]);
          bodydata.tfaces[j].tris[3*i+2] = faceIndex(bodydata.tfaces[j].tris[3*i+2], &bodydata.tfaces[j]);
        }
      }
    }

    // Allocate surfaceMesh from number of bodies
    AIM_ALLOC(mesh->meshRef->maps, 1, aimMeshTessMap, aimInfo, status);
    mesh->meshRef->nmap = 1;
    mesh->meshRef->maps[0].tess = NULL;
    mesh->meshRef->maps[0].map = NULL;

    // Build up the body tessellation object
    status = EG_initTessBody(body, &tess);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(tess, aimInfo, status);

    for ( iedge = 0; iedge < bodydata.nedges; iedge++ ) {

      // Check if the edge is degenerate
      if (bodydata.edges[iedge]->mtype == DEGENERATE) continue;

      status = EG_setTessEdge(tess, iedge+1, bodydata.tedges[iedge].npts,
                              bodydata.tedges[iedge].xyz,
                              bodydata.tedges[iedge].t);
      AIM_STATUS(aimInfo, status, "Failed to set tessellation on Edge %d!", iedge+1);

      // Add the unique indexing of the edge tessellation
      snprintf(attrname, 128, "edgeVertID_%d",iedge+1);
      status = EG_attributeAdd(tess, attrname, ATTRINT,
                               bodydata.tedges[iedge].npts,
                               bodydata.tedges[iedge].ivp, NULL, NULL);
      AIM_STATUS(aimInfo, status);
    }

    if (mesh->meshRef->type == aimAreaMesh) {

      params[0] = -1;
      params[1] = 1;
      params[2] = 20;

      status = EG_finishTess( tess, params );
      AIM_STATUS(aimInfo, status);

    } else {
      for (iface = 0; iface < bodydata.nfaces; iface++) {

        ntri  = bodydata.tfaces[iface].ntri;

        face_tris = bodydata.tfaces[iface].tris;
        face_uv   = bodydata.tfaces[iface].uv;
        face_xyz  = bodydata.tfaces[iface].xyz;

        AIM_NOTNULL(face_tris, aimInfo, status);
        AIM_NOTNULL(face_uv  , aimInfo, status);
        AIM_NOTNULL(face_xyz , aimInfo, status);

        // check the normals of the elements match the geometry normals
        // only need to check one element per face to decide for all
        elem[0] = face_tris[0]-1;
        elem[1] = face_tris[1]-1;
        elem[2] = face_tris[2]-1;

        // get the uv centroid
        uv[0] = (face_uv[2*elem[0]+0] + face_uv[2*elem[1]+0] + face_uv[2*elem[2]+0])/3.;
        uv[1] = (face_uv[2*elem[0]+1] + face_uv[2*elem[1]+1] + face_uv[2*elem[2]+1])/3.;

        // get the normal of the face
        status = EG_evaluate(bodydata.faces[iface], uv, result);
        AIM_STATUS(aimInfo, status);

        // use cross dX/du x dX/dv to get geometry normal
        v1[0] = result[3]; v1[1] = result[4]; v1[2] = result[5];
        v2[0] = result[6]; v2[1] = result[7]; v2[2] = result[8];
        CROSS(faceNormal, v1, v2);

        // get mtype=SFORWARD or mtype=SREVERSE for the face to get topology normal
        status = EG_getInfo(bodydata.faces[iface], &oclass, &mtype, &ref, &prev, &next);
        AIM_STATUS(aimInfo, status);
        faceNormal[0] *= mtype;
        faceNormal[1] *= mtype;
        faceNormal[2] *= mtype;

        // get the normal of the mesh triangle
        v1[0] = face_xyz[3*elem[1]+0] - face_xyz[3*elem[0]+0];
        v1[1] = face_xyz[3*elem[1]+1] - face_xyz[3*elem[0]+1];
        v1[2] = face_xyz[3*elem[1]+2] - face_xyz[3*elem[0]+2];

        v2[0] = face_xyz[3*elem[2]+0] - face_xyz[3*elem[0]+0];
        v2[1] = face_xyz[3*elem[2]+1] - face_xyz[3*elem[0]+1];
        v2[2] = face_xyz[3*elem[2]+2] - face_xyz[3*elem[0]+2];

        CROSS(triNormal, v1, v2);

        // get the dot product between the triangle and face normals
        ndot = DOT(faceNormal,triNormal);

        // if the normals are opposite, swap all triangles
        if (ndot < 0) {
          // swap two vertices to reverse the normal
          for (i = 0; i < ntri; i++) {
            swapi(&face_tris[3*i+0], &face_tris[3*i+2]);
          }
        }

        status = EG_setTessFace(tess,
                                iface+1,
                                bodydata.tfaces[iface].npts,
                                face_xyz,
                                face_uv,
                                ntri,
                                face_tris);
        AIM_STATUS(aimInfo, status);


        // The points get reindexed to be consistent with loops in EG_setTessFace
        // This uses the new triangulation to map that index change.
        status = EG_getTessFace(tess, iface+1, &npts, &pxyz, &puv, &ptype,
                                &pindex, &ntri, &tris, &tric);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(tris, aimInfo, status);

        AIM_ALLOC(faceVertID, npts, int, aimInfo, status);

        for (i = 0; i < ntri; i++) {
          for (j = 0; j < 3; j++) {
            faceVertID[tris[3*i+j]-1] = bodydata.tfaces[iface].ivp[face_tris[3*i+j]-1];
          }
        }

        // Add the unique indexing of the tessellation
        snprintf(attrname, 128, "faceVertID_%d",iface+1);
        status = EG_attributeAdd(tess, attrname, ATTRINT,
                                 bodydata.tfaces[iface].npts,
                                 faceVertID, NULL, NULL);
        AIM_STATUS(aimInfo, status);

        // replace the shuffled volume ID's
        AIM_FREE(bodydata.tfaces[iface].ivp);
        bodydata.tfaces[iface].ivp = faceVertID;
        faceVertID = NULL;
      }
    }

    // finalize the tessellation
    status = EG_statusTessBody(tess, &body, &i, &nglobal);
    AIM_STATUS(aimInfo, status, "Tessellation object was not built correctly!!!");

    // save the tessellation with caps
    status = aim_newTess(aimInfo, tess);
    AIM_STATUS(aimInfo, status);

    /*@-kepttrans@*/
    // reference the surface mesh object
    mesh->meshRef->maps[0].tess = tess;
    tess = NULL;
    /*@+kepttrans@*/

    // Create the map from the tessellation global vertex index to the volume mesh vertex index
    AIM_ALLOC(mesh->meshRef->maps[0].map, nglobal, int, aimInfo, status);

    if (mesh->meshRef->type == aimSurfaceMesh) {

      /* construct global vertex indices
       * EGADS re-ordres vertexes, so the meshb file must be written
       * back out to be consistent with the tess object. This implies
       * the mapping is simply an identity.
       */
      for (i = 0; i < nglobal; i++) {
        mesh->meshRef->maps[0].map[i] = i+1;
      }

    } else if (mesh->meshRef->type == aimAreaMesh) {

      // Find the boundary mesh in the global tessellation
      for (i = 0; i < nglobal; i++) {

        // Get the local indexes from the boundary mesh
        status = EG_getGlobal(mesh->meshRef->maps[0].tess, i+1, &localIndex, &topoIndex, NULL);
        AIM_STATUS(aimInfo, status);

        // Get the global index in the full 2D mesh
        if (localIndex == 0) {
          status = EG_localToGlobal(mesh->meshRef->maps[0].tess, 0, topoIndex, &iglobal);
          AIM_STATUS(aimInfo, status);
        } else if (topoIndex > 0) {
          status = EG_localToGlobal(mesh->meshRef->maps[0].tess, -topoIndex, localIndex, &iglobal);
          AIM_STATUS(aimInfo, status);
        } else {
          AIM_ERROR(aimInfo, "Developer exception! Should not find Face index!");
          status = CAPS_NOTIMPLEMENT;
          goto cleanup;
        }

        mesh->meshRef->maps[0].map[i] = iglobal;
      }

    } else {

      for (iface = 0; iface < bodydata.nfaces; iface++) {
        status = EG_getTessFace(mesh->meshRef->maps[0].tess, iface+1, &npts, &pxyz, &puv, &ptype,
                                &pindex, &ntri, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        /* construct global vertex indices */
        for (i = 0; i < npts; i++) {
          status = EG_localToGlobal(mesh->meshRef->maps[0].tess, iface+1, i+1, &iglobal);
          AIM_STATUS(aimInfo, status);
          mesh->meshRef->maps[0].map[iglobal-1] = bodydata.tfaces[iface].ivp[i];
        }
      }
    }

    AIM_ALLOC(mesh->meshRef->bnds, metrisInstance->groupMap.numAttribute, aimMeshBnd, aimInfo, status);
    mesh->meshRef->nbnd = metrisInstance->groupMap.numAttribute;
    for (i = 0; i < mesh->meshRef->nbnd; i++) {
      status = aim_initMeshBnd(mesh->meshRef->bnds + i);
      AIM_STATUS(aimInfo, status);
    }

    for (i = 0; i < mesh->meshRef->nbnd; i++) {
      AIM_STRDUP(mesh->meshRef->bnds[i].groupName, metrisInstance->groupMap.attributeName[i], aimInfo, status);
      mesh->meshRef->bnds[i].ID = metrisInstance->groupMap.attributeIndex[i];
    }
  }

  mesh->meshData = meshData;
  meshData = NULL;

  status = CAPS_SUCCESS;

cleanup:
  if (status != CAPS_SUCCESS) {
    aim_freeMeshData(meshData);
    AIM_FREE(meshData);
  }

  if (fileID != 0) GmfCloseMesh(fileID);

  destroy_bodyData(1, &bodydata);

  AIM_FREE(faceGroups);
  AIM_FREE(edgeGroups);

  return status;
}


static
int writeLibMeshbSurface(void *aimInfo, aimMesh *mesh)
{
  int status; // Function return status
  int  i, j;
  int state, nglobal, id;
  int nNode, nEdge, nFace;
  int iedge, iface;
  int nNodeOffset, nEdgeOffset, nFaceOffset;
  int nNodeVerts, nEdgeVerts, nFaceVerts;
  int nLine, nTri;
  int plen, tlen, iglobal;
  int elem[3];
  const double *points, *uv, *t;
  const int *ptypes, *pindexs, *tris, *tric;
  double xyz[3];
  ego body, *edges=NULL;
  char filename[PATH_MAX];
  int64_t fileID=0;
  int meshVersion;
  aimMeshRef *meshRef = NULL;
  aimMeshData *meshData = NULL;

  if (mesh == NULL) return CAPS_NULLVALUE;
  if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
  if (mesh->meshData == NULL) return CAPS_NULLVALUE;

  meshRef  = mesh->meshRef;
  meshData = mesh->meshData;

  if (meshData->dim != 2 && meshData->dim != 3) {
    AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", mesh->meshData->dim);
    return CAPS_BADVALUE;
  }

  snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

  meshVersion = 2;
  if (EXPORT_MESHB_VERTEX_3 < meshData->nVertex) meshVersion = 3;
  if (EXPORT_MESHB_VERTEX_4 < meshData->nVertex) meshVersion = 4;

  fileID = GmfOpenMesh(filename, GmfWrite, meshVersion, meshData->dim);

  if (fileID == 0) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    return CAPS_IOERR;
  }

  status = GmfSetKwd(fileID, GmfVertices, meshData->nVertex);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  // Write nodal coordinates
  for (i = 0; i < meshData->nVertex; i++) {
    int localIndex, topoIndex;
    status = EG_getGlobal(meshRef->maps[0].tess, i+1,
                          &localIndex, &topoIndex, xyz);
    AIM_STATUS(aimInfo, status);

    status = GmfSetLin(fileID, GmfVertices, xyz[0],
                                            xyz[1],
                                            xyz[2], EXPORT_MESHB_VERTEX_ID);
    if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  }

  // write out elements

  // count the number of EDGE/FACE elements
  nTri = nLine = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(edges, aimInfo, status)

    for (iedge = 0; iedge < nEdge; iedge++) {
      if (edges[iedge]->mtype == DEGENERATE) continue;
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      AIM_STATUS(aimInfo, status);
      nLine += plen-1;
    }
    AIM_FREE(edges);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);
      nTri += tlen;
    }
  }

  // Write out EDGE line elements
  status = GmfSetKwd(fileID, GmfEdges, nLine);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nEdgeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      if (status == EGADS_DEGEN) continue;
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < plen-1; j++) {

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 1, &elem[0]);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 2, &elem[1]);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        status = GmfSetLin(fileID, GmfEdges, elem[0],
                                             elem[1],
                                             nEdgeOffset + iedge + 1);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }
    nEdgeOffset += nEdge;
  }

  // Write FACE triangle elements
  status = GmfSetKwd(fileID, GmfTriangles, nTri);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nFaceOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < tlen; j++) {
        /* triangle orientation flipped, per metris convention */
        status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 0], &elem[1/*0*/]);
        AIM_STATUS(aimInfo, status);
        status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 1], &elem[0/*1*/]);
        AIM_STATUS(aimInfo, status);
        status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, tris[3*j + 2], &elem[2]);
        AIM_STATUS(aimInfo, status);

        status = GmfSetLin(fileID, GmfTriangles, elem[0],
                                                 elem[1],
                                                 elem[2],
                                                 nFaceOffset + iface + 1);
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    nFaceOffset += nFace;
  }

  nNodeVerts = nEdgeVerts = nFaceVerts = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, NODE, &nNode, NULL);
    AIM_STATUS(aimInfo, status);
    nNodeVerts += nNode;

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(edges, aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      if (edges[iedge]->mtype == DEGENERATE) continue;
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      AIM_STATUS(aimInfo, status);
      nEdgeVerts += plen;
    }
    AIM_FREE(edges);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);
      nFaceVerts += plen;
    }
  }

  // write out parametric coordinates

  // Write NODEs
  status = GmfSetKwd(fileID, GmfVerticesOnGeometricVertices, nNodeVerts);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nNodeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, NODE, &nNode, NULL);
    AIM_STATUS(aimInfo, status);

    for (j = 0; j < nNode; j++) {
      status = EG_localToGlobal(meshRef->maps[i].tess, 0, j + 1, &iglobal);
      AIM_STATUS(aimInfo, status);
      status = GmfSetLin(fileID, GmfVerticesOnGeometricVertices,
                         iglobal,
                         nNodeOffset + j + 1);
      if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    }

    nNodeOffset += nNode;
  }

  // Write EDGEs
  status = GmfSetKwd(fileID, GmfVerticesOnGeometricEdges, nEdgeVerts);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nEdgeOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, NULL);
    AIM_STATUS(aimInfo, status);

    for (iedge = 0; iedge < nEdge; iedge++) {
      status = EG_getTessEdge(meshRef->maps[i].tess, iedge + 1, &plen, &points, &t);
      if (status == EGADS_DEGEN) continue;
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < plen; j++) {

        status = EG_localToGlobal(meshRef->maps[i].tess, -(iedge + 1), j + 1, &iglobal);
        if (status == EGADS_DEGEN) continue;
        AIM_STATUS(aimInfo, status);

        id = nEdgeOffset + iedge + 1;
        status = GmfSetLin(fileID, GmfVerticesOnGeometricEdges,
                           iglobal,
                           id,
                           t[j],
                           (double)id);  // metris abuse of dist
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    nEdgeOffset += nEdge;
  }


  // Write FACEs
  status = GmfSetKwd(fileID, GmfVerticesOnGeometricTriangles, nFaceVerts);
  if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  nFaceOffset = 0;
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      for (j = 0; j < plen; j++) {

        status = EG_localToGlobal(meshRef->maps[i].tess, iface + 1, j + 1, &iglobal);
        AIM_STATUS(aimInfo, status);

        id = nFaceOffset + iface + 1;
        status = GmfSetLin(fileID, GmfVerticesOnGeometricTriangles,
                           iglobal,
                           id,
                           uv[2*j+0], uv[2*j+1],
                           (double)id); // metris abuse of dist
        if (status <= 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
      }
    }

    nFaceOffset += nFace;
  }

  //printf("Finished writing libMeshb file\n\n");

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(edges);
  if (fileID != 0) GmfCloseMesh(fileID);
  return status;
}
