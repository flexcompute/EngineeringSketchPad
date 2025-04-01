/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             FlightStream AIM
 *
 *     Written by John Dannenhoffer (Syracuse) and
 *                Marshall Galbraith (MIT)
 *
 *
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewFLIGHTSTREAM FlightStream AIM Overview

 * A module in the Computational Aircraft Prototype Syntheses (CAPS)
 * has been developed to interact (primarily through input files) with
 * Research in Flight's FlightStream \cite FlightStream. FlightStream predicts
 * the aerodynamic performance of a vehicle via a panel methid, which
 * makes it very fast.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref
 * aimInputsFLIGHTSTREAM and \ref aimOutputsFLIGHTSTREAM, respectively.
 *
 * Geometric attributes recognized by the AIM are provided in \ref
 * attributeFLIGHTSTREAM.
 *
 * The accepted and expected geometric representation are detailed in
 * \ref geomRepIntentFLIGHTSTREAM.
 *
 */
 /*! \page attributeFLIGHTSTREAM Attribution
 *
 * The following list of attributes drives the FlightStream geometric definition.
 *
 *  - <b> capsLength</b> This attribute defines the length units that
 *  the *.csm file is generated in and is not optional for FligtStream.
 *  The FlightStream input grid will be scaled to either the default length of METER
 *  or the user specified length unit (see \ref aimUnitsFLIGHTSTREAM).
 *
 *  - <b> capsReferenceChord</b> and <b> capsReferenceSpan</b>
 * [Optional] These attributes may exist on any <em> Body</em>. Their
 * value will be used as the reference moment lengths in FlightStream's
 * input file with their units assumed to be consistent with the
 * attribute "capsLength". These values may be alternatively set through an input
 * value, "ReferenceChord" (see \ref aimInputsFLIGHTSTREAM)
 *
 *  - <b> capsReferenceArea</b> [Optional] This attribute may exist on
 * any <em> Body</em>.  Its value will be used as the reference area
 * in FlightStream's input file with its units assumed to be consistent with
 * the attribute "capsLength". This value may be alternatively set
 * through an input value, "ReferenceArea" (see \ref aimInputsFLIGHTSTREAM)
 *
 *  - <b> capsReferenceX</b>, <b> capsReferenceY</b>,  and <b> capsReferenceZ</b>
 * [Optional] These attribute may exist on any <em> Body</em>. Their
 * value will be used as the reference moment lengths in FlightStream's
 * input file with their units assumed to be consistent with the
 * attribute "capsLength". These values may be alternatively set through an input
 * value, "ReferenceX" (see \ref aimInputsFLIGHTSTREAM)
 *

 *
 */

//#define DEBUG

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "miscUtils.h"
#include "meshUtils.h"
#include "cfdUtils.h"

#include "wavefrontWriter.h"

#ifdef WIN32
   #define snprintf    _snprintf
   #define strcasecmp  stricmp
   #define strncasecmp _strnicmp
   #define strtok_r    strtok_s
#else
   #include <unistd.h>
   #include <limits.h>
#endif

/* define the indicies of the inputs and outputs */
enum aimInputs
{
    inProj_Name = 1,             /* index is 1-biased */
    inFlightStream,
    inMach,
    inAlpha,
    inBeta,
    inFluid_Properties,
    inAltitude,
    inReferenceChord,
    inReferenceSpan,
    inReferenceArea,
    inReferenceX,
    inReferenceY,
    inReferenceZ,
   // inReferenceLength,
    inReferenceVelocity,
    inPressure_Scale_Factor,
    inConvergenceTol,
    inMaxIterations,
    inExport_Solver_Analysis,
    inFlightScript,
    inMesh_Morph,
    inSurface_Mesh,
    NUMINPUT = inSurface_Mesh    /* Total number of inputs */
};

enum aimOutputs
{
    outCx = 1,                   /* index is 1-based */
    outCy,
    outCz,
    outCL,
    outCDi,
    outCDo,
    outCMx,
    outCMy,
    outCMz,
    NUMOUTPUT = outCMz           /* Total number of outputs */
};

#define MXCHAR  255

#define ANALYSIS_VTK

#ifdef ANALYSIS_PLOAD_BDF
static const char* loadbdf = "surface_load.bdf";
#elif defined(ANALYSIS_VTK)
static const char* cpvtk = "surface_cp.vtk";
#elif defined(ANALYSIS_CSV)
static const char* dpcsv = "surface_dp.csv";
#else
#error "must define ANALYSIS_PLOAD_BDF, ANALYSIS_VTK, or ANALYSIS_CSV"
#endif

typedef struct {

  // Fluid properties in specific units
  double rhoref;        // kg/m^3
  double pref;          // Pa
  double sonicspeedref; // m/s
  double tref;          // K
  double muref;         // Pa-sec

  // variables that could be obtained via attributes on the Body
  double Cref;       // capsReferenceChord
  double Bref;       // capsReferenceSpan
  double Sref;       // capsReferenceArea
  double Xref;       // capsReferenceX   (cg)
  double Yref;       // capsReferenceY   (cg)
  double Zref;       // capsReferenceZ   (cg)

  // outputs from FlightStream
  double Cx;         // X-force coefficient
  double Cy;         // Y-force coefficient
  double Cz;         // Z-force coefficient
  double CL;         // lift coefficient
  double CDi;        // drag coefficient (inviscid)
  double CDo;        // drag coefficient (total)
  double CMx;        // X-moment coefficient
  double CMy;        // Y-moment coefficient
  double CMz;        // Z-moment coefficient

  // Attribute to index map
  mapAttrToIndexStruct groupMap;

  // Units structure
  cfdUnitsStruct units;

  // Mesh reference obtained from meshing AIM
  const aimMeshRef *meshRefIn;
  aimMeshRef meshRefMorph, meshRefDisp;

} aimStorage;


static int initialize_aimStorage(aimStorage *flightstreamInstance)
{

  // Set initial values for flightstreamInstance
  flightstreamInstance->rhoref = 0;
  flightstreamInstance->pref = 0;
  flightstreamInstance->sonicspeedref = 0;
  flightstreamInstance->tref = 0;
  flightstreamInstance->muref = 0;

  flightstreamInstance->Cref = 0;
  flightstreamInstance->Bref = 0;
  flightstreamInstance->Sref = 0;
  flightstreamInstance->Xref = 0;
  flightstreamInstance->Yref = 0;
  flightstreamInstance->Zref = 0;

  // Container for attribute to index map
  (void)initiate_mapAttrToIndexStruct(&flightstreamInstance->groupMap);

  initiate_cfdUnitsStruct(&flightstreamInstance->units);

  flightstreamInstance->meshRefIn = NULL;
  aim_initMeshRef(&flightstreamInstance->meshRefMorph, aimUnknownMeshType);
  aim_initMeshRef(&flightstreamInstance->meshRefDisp, aimUnknownMeshType);

  return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *flightstreamInstance)
{
  // Attribute to index map
  (void) destroy_mapAttrToIndexStruct(&flightstreamInstance->groupMap);

  destroy_cfdUnitsStruct(&flightstreamInstance->units);

  // Surface Mesh
  aim_freeMeshRef(&flightstreamInstance->meshRefMorph);
  aim_freeMeshRef(&flightstreamInstance->meshRefDisp);
  flightstreamInstance->meshRefIn = NULL;

  initialize_aimStorage(flightstreamInstance);

  return CAPS_SUCCESS;
}


// Get the loads from a capsTuple
static int
_getFluid_Properties(void *aimInfo,
                     int numPropTuple,
                     capsTuple propTuple[],
                     aimStorage *flightstreamInstance) {

  /*! \page flightstreamFLUID_PROPERTIES FlightStream Fluid Properties
   * Structure for the load tuple  = ("Property", "Value").
   */

  int status; //Function return

  int i; // Indexing

  char *keyValue = NULL; // Key values from tuple searches
  const char *keyWord = NULL; // Key words to find in the tuples

  printf("\nGetting Fluid Properties.......\n");

  if (numPropTuple != 5) {
    AIM_ERROR(aimInfo, "All 5 fluid properties must be specified. Only %d provided:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s", propTuple[i].name);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Loop through tuples and get each property
  for (i = 0; i < numPropTuple; i++) {

    /*! \page flightstreamFLUID_PROPERTIES
     *
     * <ul>
     *  <li> <B>density = "(no default)"</B> </li> <br>
     *  Reference density
     * </ul>
     *
     */
    keyWord = "density";
    if (strcasecmp(propTuple[i].name, keyWord) == 0) {
      status = string_toDoubleUnits(aimInfo, propTuple[i].value, "kg/m^3", &flightstreamInstance->rhoref);
      AIM_STATUS(aimInfo, status, "while parsing '%s'", propTuple[i].value);
    }

    /*! \page flightstreamFLUID_PROPERTIES
     *
     * <ul>
     *  <li> <B>pressure = "(no default)"</B> </li> <br>
     *  Reference pressure
     * </ul>
     *
     */
    keyWord = "pressure";
    if (strcasecmp(propTuple[i].name, keyWord) == 0) {
      status = string_toDoubleUnits(aimInfo, propTuple[i].value, "Pa", &flightstreamInstance->pref);
      AIM_STATUS(aimInfo, status, "while parsing '%s'", propTuple[i].value);
    }

    /*! \page flightstreamFLUID_PROPERTIES
     *
     * <ul>
     *  <li> <B>sonic_velocity = "(no default)"</B> </li> <br>
     *  Reference speed of sound
     * </ul>
     *
     */
    keyWord = "sonic_velocity";
    if (strcasecmp(propTuple[i].name, keyWord) == 0) {
      status = string_toDoubleUnits(aimInfo, propTuple[i].value, "m/s", &flightstreamInstance->sonicspeedref);
      AIM_STATUS(aimInfo, status, "while parsing '%s'", propTuple[i].value);
    }

    /*! \page flightstreamFLUID_PROPERTIES
     *
     * <ul>
     *  <li> <B>temperature = "(no default)"</B> </li> <br>
     *  Reference temperature
     * </ul>
     *
     */
    keyWord = "temperature";
    if (strcasecmp(propTuple[i].name, keyWord) == 0) {
      status = string_toDoubleUnits(aimInfo, propTuple[i].value, "K", &flightstreamInstance->tref);
      AIM_STATUS(aimInfo, status, "while parsing '%s'", propTuple[i].value);
    }

    /*! \page flightstreamFLUID_PROPERTIES
     *
     * <ul>
     *  <li> <B>viscosity = "(no default)"</B> </li> <br>
     *  Reference viscosity
     * </ul>
     *
     */
    keyWord = "viscosity";
    if (strcasecmp(propTuple[i].name, keyWord) == 0) {
      status = string_toDoubleUnits(aimInfo, propTuple[i].value, "Pa*s", &flightstreamInstance->muref);
      AIM_STATUS(aimInfo, status, "while parsing '%s'", propTuple[i].value);
    }
  }

  if (flightstreamInstance->rhoref <= 0) {
    AIM_ERROR(aimInfo, "Positive 'density' must be specified. Provided Fluid_Properties:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s -> %s", propTuple[i].name, propTuple[i].value);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (flightstreamInstance->pref <= 0) {
    AIM_ERROR(aimInfo, "Positive 'pressure' must be specified. Provided Fluid_Properties:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s -> %s", propTuple[i].name, propTuple[i].value);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (flightstreamInstance->sonicspeedref <= 0) {
    AIM_ERROR(aimInfo, "Positive 'sonic_velocity' must be specified. Provided Fluid_Properties:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s -> %s", propTuple[i].name, propTuple[i].value);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (flightstreamInstance->tref <= 0) {
    AIM_ERROR(aimInfo, "Positive 'temperature' must be specified. Provided Fluid_Properties:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s -> %s", propTuple[i].name, propTuple[i].value);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (flightstreamInstance->muref <= 0) {
    AIM_ERROR(aimInfo, "Positive 'viscosity' must be specified. Provided Fluid_Properties:");
    for (i = 0; i < numPropTuple; i++)
      AIM_ADDLINE(aimInfo, "%s -> %s", propTuple[i].name, propTuple[i].value);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  printf("\tDone getting Fluid Properties\n");
  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(keyValue);

  return status;
}


// Get FlightStream data transfer files
static int
fs_dataTransfer(void *aimInfo,
                capsValue  aimInputs[],
                const aimStorage *fsInstance)
{

    /*! \page dataTransferFLIGHTSTREAM FlightStream Data Transfer
     *
     * \section dataToCart3D Data transfer to FlightStream (FieldIn)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements (as from a structural solver)
     *   and updates FlightStream surface mesh.
     * </ul>
     */

    int status; // Function return status
    int i, ibound, iglobal; // Indexing

    // Discrete data transfer variables
    capsDiscr *discr;
    char **boundNames = NULL;
    int numBoundName = 0;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char *units;

    const aimMeshRef *meshRefIn;

    aimMesh    mesh;

    // Data transfer variable
    double *dxyz = NULL;

    int foundDisplacement = (int) false;

    status = aim_getBounds(aimInfo, &numBoundName, &boundNames);
    if (status == CAPS_NOTFOUND) return CAPS_SUCCESS;
    AIM_STATUS(aimInfo, status);

    // Get mesh
    meshRefIn = (aimMeshRef *) aimInputs[inSurface_Mesh-1].vals.AIMptr;

    if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
        meshRefIn == NULL) { // If we are mighty morphing
      meshRefIn = &fsInstance->meshRefMorph;
    }

    foundDisplacement = (int) false;
    for (ibound = 0; ibound < numBoundName; ibound++) {
      AIM_NOTNULL(boundNames, aimInfo, status);

      status = aim_getDiscr(aimInfo, boundNames[ibound], &discr);
      if (status != CAPS_SUCCESS) continue;

      status = aim_getDataSet(discr,
                              "Displacement",
                              &dataTransferMethod,
                              &numDataTransferPoint,
                              &dataTransferRank,
                              &dataTransferData,
                              &units);
      if (status != CAPS_SUCCESS) continue;

      foundDisplacement = (int) true;

      // Is the rank correct?
      if (dataTransferRank != 3) {
        AIM_ERROR(aimInfo, "Displacement transfer data found however rank is %d not 3!!!!", dataTransferRank);
        status = CAPS_BADRANK;
        goto cleanup;
      }

    } // Loop through bound names

    if (foundDisplacement != (int) true ) {
        status = CAPS_SUCCESS;
        goto cleanup;
    }

    status = aim_freeMeshRef(&((aimStorage *)fsInstance)->meshRefDisp);
    AIM_STATUS(aimInfo, status);

    /* Create a local mesh file name */
    status = aim_localMeshRef(aimInfo, meshRefIn, &((aimStorage *)fsInstance)->meshRefDisp);
    AIM_STATUS(aimInfo, status);

    /*@-immediatetrans@*/
    mesh.meshData = NULL;
    mesh.meshRef = &((aimStorage *)fsInstance)->meshRefDisp;
    /*@+immediatetrans@*/

    status = mesh_surfaceMeshData(aimInfo, &fsInstance->groupMap, &mesh);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(mesh.meshData, aimInfo, status);

    /*@-immediatetrans@*/
    ((aimStorage *)fsInstance)->meshRefIn = &((aimStorage *)fsInstance)->meshRefDisp;
    /*@+immediatetrans@*/

    AIM_ALLOC(dxyz, 3*mesh.meshData->nVertex, double, aimInfo, status);
    memset(dxyz, 0, 3*mesh.meshData->nVertex*sizeof(double));

    // gather displacements first to avoid double counting edges/nodes
    for (ibound = 0; ibound < numBoundName; ibound++) {
      AIM_NOTNULL(boundNames, aimInfo, status);

      status = aim_getDiscr(aimInfo, boundNames[ibound], &discr);
      if (status != CAPS_SUCCESS) continue;

      status = aim_getDataSet(discr,
                              "Displacement",
                              &dataTransferMethod,
                              &numDataTransferPoint,
                              &dataTransferRank,
                              &dataTransferData,
                              &units);
      if (status != CAPS_SUCCESS) continue;

      if (numDataTransferPoint != discr->nPoints &&
          numDataTransferPoint > 1) {
        AIM_ERROR(aimInfo, "Developer error!! %d != %d", numDataTransferPoint, discr->nPoints);
        status = CAPS_MISMATCH;
        goto cleanup;
      }

      for (i = 0; i < discr->nPoints; i++) {
        iglobal = discr->tessGlobal[2*i+1];

        if (numDataTransferPoint == 1) {
          // A single point means this is an initialization phase
          status = aim_convert(aimInfo, 3, units, dataTransferData, fsInstance->units.length, &dxyz[3*(iglobal-1)+0]);
          AIM_STATUS(aimInfo, status);
        } else {
          // Apply delta displacements
          status = aim_convert(aimInfo, 3, units, &dataTransferData[3*i+0], fsInstance->units.length, &dxyz[3*(iglobal-1)+0]);
          AIM_STATUS(aimInfo, status);
        }
      }
    } // Loop through bound names

    // Apply the displacements
    for (i = 0; i < mesh.meshData->nVertex; i++) {
      mesh.meshData->verts[i][0] += dxyz[3*i+0];
      mesh.meshData->verts[i][1] += dxyz[3*i+1];
      mesh.meshData->verts[i][2] += dxyz[3*i+2];
    }

    /* write the mesh */
    status = aim_writeMesh(aimInfo, MESHWRITER, fsInstance->units.length, &mesh);
    AIM_STATUS(aimInfo, status);

    /* cleanup */
    status = aim_freeMeshData(mesh.meshData);
    AIM_STATUS(aimInfo, status);
    AIM_FREE(mesh.meshData);

    status = CAPS_SUCCESS;

// Clean-up
cleanup:
    AIM_FREE(dxyz);
    AIM_FREE(boundNames);

    return status;
}

/**********************************************************************/
/* aimInitialize - initialize the AIM                                 */
/**********************************************************************/

int
aimInitialize(int        qFlag,         /* (in)  -1 indiactes query only */
 /*@unused@*/ const char unitSys[],     /* (in)  unit system */
              void       *aimInfo,      /* (in)  AIM context */
              void       **instStore,   /* (out) AIM instance storage */
              int        *major,        /* (out) major version number */
              int        *minor,        /* (out) minor version number */
              int        *nIn,          /* (out) number of inputs to AIM */
              int        *nOut,         /* (out) number of outputs from AIM */
              int        *nFields,      /* (out) number of DataSet fields */
              char       ***fnames,     /* (out) array  of DataSet names */
              int        **franks,      /* (out) array  of DataSet ranks */
              int        **fInOut)      /* (out) array  of field flags */
{

    int status = CAPS_SUCCESS;

    int  *ints=NULL, i;
    char **strs=NULL;
    const char *keyWord;
    char *keyValue = NULL;
    double real = 1;
    cfdUnitsStruct *units=NULL;

    aimStorage *flightstreamInstance = NULL;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimInitialize(qFlag=%d, unitSys=%s)\n", qFlag, unitSys);
#endif

    /* specify the version number */
    *major = 1;
    *minor = 0;

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;

    /* if this is simply a query, return now */
    if (qFlag == -1) {
        return CAPS_SUCCESS;
    }

    /*! \page geomRepIntentFLIGHTSTREAM Geometry Representation The geometric
     * representation for the FlightStream AIM requires that the body be
     * either a solid body (SOLIDBODY) or a manifold sheet body
     * (SHEETBODY).
     */

    /* specify the field variables this analysis can generate and consume */
    *nFields = 2;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
    strs[0]  = EG_strdup("Pressure");
    strs[1]  = EG_strdup("Displacement");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
      }
    *fnames  = strs;

    /* specify the dimension of each field variable */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = 1;
    ints[1]  = 3;
    *franks   = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldOut;
    ints[1]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate flightstreamInstance
    AIM_ALLOC(flightstreamInstance, 1, aimStorage, aimInfo, status);

    initialize_aimStorage(flightstreamInstance);

    units = &flightstreamInstance->units;

    *instStore = flightstreamInstance;


    /*! \page aimUnitsFLIGHTSTREAM AIM Units
     *  FlightStream expects units for all inputs, and by default the AIM uses SI units, i.e.
     *   -  mass : kg
     *   - length : meter
     *   - time : seconds
     *   - temperature : K
     *  A unit system may be optionally specified during AIM instance initiation to use a different set of base units.
     *  A unit system may be specified via a JSON string dictionary for example:
     *  unitSys = "{"mass": "lb", "length": "feet", "time":"seconds", "temperature": "R"}"
     */
    if (unitSys != NULL) {

      // Do we have a json string?
      if (strncmp( unitSys, "{", 1) != 0) {
        AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
        return CAPS_BADVALUE;
      }

      /*! \page aimUnitsFLIGHTSTREAM
       *  \section jsonStringFLIGHTSTREAM JSON String Dictionary
       *  The key arguments of the dictionary are described in the following:
       *
       *  <ul>
       *  <li> <B>mass = "None"</B> </li> <br>
       *  Mass units - e.g. "kilogram", "k", "slug", ...
       *  </ul>
       */
      keyWord = "mass";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->mass = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->mass, &real, "kg", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->mass, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsFLIGHTSTREAM
       *  <ul>
       *  <li> <B>length = "None"</B> </li> <br>
       *  Length units - FlightStream support: "inch", "millimeter", "feet", "mile", "meter", "kilometer", "mils", "micron", "centimeter", "microinch"
       *  </ul>
       */
      keyWord = "length";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->length = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->length, &real, "m", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->length, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsFLIGHTSTREAM
       *  <ul>
       *  <li> <B>time = "None"</B> </li> <br>
       *  Time units - e.g. "second", "s", "minute", ...
       *  </ul>
       */
      keyWord = "time";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->time = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->time, &real, "s", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->time, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsFLIGHTSTREAM
       *  <ul>
       *  <li> <B>temperature = "None"</B> </li> <br>
       *  Temperature units - e.g. "Kelvin", "K", "degC", ...
       *  </ul>
       */
      keyWord = "temperature";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->temperature = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->temperature, &real, "K", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->temperature, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    } else {

      // default SI units
      AIM_STRDUP(units->mass, "kg", aimInfo, status);
      AIM_STRDUP(units->length, "meter", aimInfo, status);
      AIM_STRDUP(units->time, "seconds", aimInfo, status);
      AIM_STRDUP(units->temperature, "K", aimInfo, status);

    }

    status = cfd_cfdDerivedUnits(aimInfo, units);
    AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


/**********************************************************************/
/* aimInputs - return information about index'th analysis input       */
/**********************************************************************/

int
aimInputs(
/*@unused@*/void    *instStore,         /* (in)  AIM instance storage */
          void      *aimInfo,           /* (in)  AIM context */
          int       index,              /* (in)  input index (1-nIn] */
          char      **ainame,           /* (out) name of analysis input */
          capsValue *inval)             /* (out) description of analysis input */
{
    int status = CAPS_SUCCESS;
    aimStorage *flightstreamInstance;

    cfdUnitsStruct *units=NULL;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimInputs(index=%d)\n", index);
#endif

    *ainame = NULL;

    /*! \page aimInputsFLIGHTSTREAM AIM Inputs
     * The following list outlines the FlightStream inputs along with their
     * default values available through the AIM interface.
     */

    // entities to be defined for each aimInput
    // default values are in [brackets]

    // *ainame             = EG_strdup("name")
    // inval->type         = Boolean, Integer, Double, String, Tuple, Pointer, DoubleDeriv, or PointerMesh
    // inval->nullVal      = [NotNull], IsNull, NotAllowed, IsPartial
    // inval->units        = [NULL] or EG_strdup("bar")
    // inval->lfixed       = [Fixed] or Change (length is fixed)
    // inval->sfixed       = [Fixed] or Change (shape  is fixed)
    // inval->dim          = [0], 1, or 2 (maximum dimensions)
    // inval->nrow         = [1] or n (only if dim>0)
    // inval->ncol         = [1] or n (only if dim>1)
    // inval->length       = inval->nrow * inval->ncol
    // inval->units        = NULL or EG_strdup("meter") or EG_strdup("degree") or ...
    // inval->vals.real    = number or AIM_ALLOC(inval->vals.reals, inval->length, double, aimInfo, status)
    // inval->vals.integer = number or AIM_ALLOC(inval->vals.integer, inval->length, int, aimInfo, status)
    // inval->vals.string  = NULL or EG_strdup("flightstream_CAPS")


    flightstreamInstance = (aimStorage *) instStore;

    if (flightstreamInstance != NULL) units = &flightstreamInstance->units;

    // FlightStream Inputs
    if (index == inProj_Name) {
        *ainame             = EG_strdup("Proj_Name");
        inval->type         = String;
        inval->nullVal      = NotNull;
        inval->lfixed       = Change;
        AIM_STRDUP(inval->vals.string, "flightstream_CAPS", aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> ProjName = "flightstream_CAPS"</B> <br>
         * Name for files generated by fightstream AIM.
         */
    } else if (index == inFlightStream) {
        *ainame             = EG_strdup("FlightStream");
        inval->type         = String;
        inval->nullVal      = NotNull;
        inval->lfixed       = Change;
        AIM_STRDUP(inval->vals.string, "FlightStream", aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> FlightStream = "FlightStream"</B> <br>
         * The name of the FlightStream executable. May include full path.
         */
    } else if (index == inMach) {
        *ainame             = EG_strdup("Mach"); // Mach number
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->units        = NULL;
        inval->lfixed       = Fixed;
        inval->sfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Mach = 0.0 (default) </B> <br>
         *  Mach number
         */

    } else if (index == inAlpha) {
        *ainame             = EG_strdup("Alpha");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->units        = EG_strdup("degree");
        inval->lfixed       = Fixed;
        inval->sfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Alpha = 0.0 (default) </B> <br>
         *  Angle of attack [degree]
         */

    } else if (index == inBeta) {
        *ainame             = EG_strdup("Beta");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->units        = EG_strdup("degree");
        inval->lfixed       = Fixed;
        inval->sfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Beta = 0.0 (default) </B> <br>
         *  Sideslip angle
         */

    } else if (index == inFluid_Properties) {
        *ainame             = EG_strdup("Fluid_Properties");
        inval->type         = Tuple;
        inval->nullVal      = IsNull;
        inval->lfixed       = Change;
        inval->sfixed       = Fixed;
        inval->dim          = Vector;
        inval->vals.tuple   = NULL;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Fluid_Properties = 0.0 (default) </B> <br>
         *  Reference fluid properties. Altitude must be NULL. See \ref flightstreamFLUID_PROPERTIES.
         */

    } else if (index == inAltitude) {
        *ainame             = EG_strdup("Altitude");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->units        = EG_strdup("ft");
        inval->lfixed       = Fixed;
        inval->sfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Altitude = NULL (default) </B> <br>
         *  Altitude used to compute Fluid Properties. The Fluid_Properties input must be NULL.
         */

    } else if (index == inReferenceChord) {
        *ainame             = EG_strdup("ReferenceChord");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceChord = NULL </B> <br>
         * This sets the reference chord for used in force and moment
         *  calculations.  Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceChord" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceSpan) {
        *ainame             = EG_strdup("ReferenceSpan");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceSpan = NULL </B> <br>
         * This sets the reference span for used in force and moment
         *  calculations.  Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceSpan" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceArea) {
        *ainame             = EG_strdup("ReferenceArea");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->area, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceArea = NULL </B> <br>
         * This sets the reference area for used in force and moment
         *  calculations.  Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceArea" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceX) {
        *ainame             = EG_strdup("ReferenceX");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceX = NULL </B> <br>
         * This sets the reference X for moment calculations.
         * Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceX" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceY) {
        *ainame             = EG_strdup("ReferenceY");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceX = NULL </B> <br>
         * This sets the reference Y for moment calculations.
         * Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceY" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceZ) {
        *ainame             = EG_strdup("ReferenceZ");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceX = NULL </B> <br>
         * This sets the reference Z for moment calculations.
         * Alternatively, the geometry (body) attribute
         *  (see \ref attributeFLIGHTSTREAM) "capsReferenceZ" maybe used to
         *  specify this variable (note: values set through the AIM input
         *  will supersede the attribution value).
         */

    } else if (index == inReferenceVelocity) {
        *ainame             = EG_strdup("ReferenceVelocity");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 0.0;
        if (units != NULL)
            AIM_STRDUP(inval->units, units->speed, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ReferenceVelocity = NULL </B> <br>
         * This sets the reference velocity
         */

    } else if (index == inPressure_Scale_Factor) {
        *ainame             = EG_strdup("Pressure_Scale_Factor");
        inval->type         = Double;
        inval->vals.real    = 1.0;

        /*! \page aimInputsCART3D
         * - <B>Pressure_Scale_Factor = 1.0</B> <br>
         * Value to scale Pressure data when transferring data. Data is scaled based on Pressure = Pressure_Scale_Factor*Pressure.
         */

    } else if (index == inConvergenceTol) {
        *ainame             = EG_strdup("ConvergenceTol");
        inval->type         = Double;
        inval->nullVal      = IsNull;
        inval->units        = NULL;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.real    = 1.0e-5;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>ConvergenceTol = 1.0e-5 (default) </B> <br>
         * Solver convergence tolerance
         */

    } else if (index == inMaxIterations) {
        *ainame             = EG_strdup("MaxIterations");
        inval->type         = Integer;
        inval->nullVal      = IsNull;
        inval->units        = NULL;
        inval->lfixed       = Fixed;
        inval->dim          = Scalar;
        inval->vals.integer = 500;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>MaxIterations = 500 </B> <br>
         * Maximum number of solver iterations
         */

    } else if (index == inExport_Solver_Analysis) {
        *ainame             = EG_strdup("Export_Solver_Analysis");
        inval->type         = String;
        inval->dim          = Vector;
        inval->lfixed       = Change;
        inval->sfixed       = Fixed;
        inval->vals.string  = NULL;
        inval->nullVal      = IsNull;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>Export_Solver_Analysis = NULL</B> <br>
         * List of file types to export. Available options:<br>
         *   - Tecplot
         *   - VTK
         *   - CSV
         *   - BDF
         *   - Force_Distributions
         */

    } else if (index == inFlightScript) {
        *ainame             = EG_strdup("FlightScript");
        inval->type         = String;
        inval->dim          = Vector;
        inval->lfixed       = Change;
        inval->sfixed       = Fixed;
        inval->vals.string  = NULL;
        inval->nullVal      = IsNull;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>FlightScript = NULL</B> <br>
         * List of flight script commands to append at the end of script.txt
         */

    } else if (index == inMesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        inval->type          = Boolean;
        inval->lfixed        = Fixed;
        inval->vals.integer  = (int) false;
        inval->dim           = Scalar;
        inval->nullVal       = NotNull;

        /*! \page aimInputsFLIGHTSTREAM
         * - <B> Mesh_Morph = False</B> <br>
         * Project previous surface mesh onto new geometry.
         */

    } else if (index == inSurface_Mesh) {
        *ainame             = EG_strdup("Surface_Mesh");
        inval->type         = PointerMesh;
        inval->dim          = Vector;
        inval->lfixed       = Change;
        inval->sfixed       = Change;
        inval->vals.AIMptr  = NULL;
        inval->nullVal      = IsNull;
        AIM_STRDUP(inval->meshWriter, MESHWRITER, aimInfo, status);
        if (units != NULL)
            AIM_STRDUP(inval->units, units->length, aimInfo, status);

        /*! \page aimInputsFLIGHTSTREAM
         * - <B>Surface_Mesh = NULL</B> <br>
         * A Surface_Mesh link.
         */

    } else {
        AIM_ERROR(aimInfo, "Unknown input index $%d", index);
        status = CAPS_RANGEERR;
        goto cleanup;
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) {
        AIM_FREE(*ainame);
    }

    return status;
}


/**********************************************************************/
/* aimOutputs - return information about index'th analysis outputs    */
/**********************************************************************/

int
aimOutputs(/*@unused@*/ void *instStore,    /* (in)  AIM instance storage */
           /*@unused@*/ void *aimInfo,      /* (in)  AIM context */
           int index,                       /* (in)  output index (1-nOut) */
           char **aoname,                   /* (out) name of analysis output */
           capsValue *outval)               /* (out) description of analysis output */
{
    int status = CAPS_SUCCESS;

    /*! \page aimOutputsFLIGHTSTREAM AIM Outputs
     * The following list outlines the FlightStream outputs available through
     * the AIM interface. All variables currently correspond to values
     * found in the *.plt file
     */

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimOutputs(index=%d)\n", index);
#endif

    outval->type       = Double;
    outval->lfixed     = Fixed;
    outval->sfixed     = Fixed;
    outval->units      = NULL;
    outval->dim        = 0;
    outval->length     = 1;
    outval->nrow       = 1;
    outval->ncol       = 1;
    outval->vals.real  = 0.0;
    outval->vals.reals = NULL;

    if        (index == outCx) {
        *aoname = EG_strdup("Cx");
    } else if (index == outCy) {
        *aoname = EG_strdup("Cy");
    } else if (index == outCz) {
        *aoname = EG_strdup("Cz");
    } else if (index == outCL) {
        *aoname = EG_strdup("CL");
    } else if (index == outCDi) {
        *aoname = EG_strdup("CDi");
    } else if (index == outCDo) {
        *aoname = EG_strdup("CDo");
    } else if (index == outCMx) {
        *aoname = EG_strdup("CMx");
    } else if (index == outCMy) {
        *aoname = EG_strdup("CMy");
    } else if (index == outCMz) {
        *aoname = EG_strdup("CMz");
    } else {
        printf("flightstreamAIM/aimOutputs index = %d NOT Found\n", index);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    /*! \page aimOutputsFLIGHTSTREAM
     * Aerodynamic coefficients:
     * - <B>Cx</B> = X-force coefficient
     * - <B>Cx</B> = X-force coefficient
     * - <B>Cx</B> = X-force coefficient
     * - <B>CL</B> = Lift coefficient
     * - <B>CDi</B> = (?) drag coefficient
     * - <B>CDo</B> = (?) drag coefficient
     * - <B>CMx</B> = X-moment coefficient
     * - <B>CMy</B> = Y-moment coefficient
     * - <B>CMz</B> = Z-moment coefficient
     */

    AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
    return status;
}


/**********************************************************************/
/* aimUpdateState - update the AIM's internal state                   */
/**********************************************************************/

int
aimUpdateState(void      *instStore,    /* (in)  AIM instance storage */
               void      *aimInfo,      /* (in)  AIM context */
               capsValue aimInputs[])   /* (in)  array of analysis inputs */
{
    // Function return flag
    int status = CAPS_SUCCESS;

    int i;
    int foundSref=(int)false;
    int foundCref=(int)false;
    int foundBref=(int)false;
    int foundXref=(int)false;
    int foundYref=(int)false;
    int foundZref=(int)false;

    // AIM input bodies
    const char *intent;
    int  numBody;
    ego *bodies = NULL;

    // EGADS return values
    int          atype, alen;
    const int    *ints;
    const char   *string;
    const double *reals;

    cfdUnitsStruct *units = NULL;

    const char *lengthUnits=NULL;
    double scaleFactor = 1.0;

    int nbound = 0;
    char **boundNames = NULL;

    aimStorage *flightstreamInstance = (aimStorage *)instStore;
    units = &flightstreamInstance->units;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimUpdateState()\n");
#endif

    AIM_NOTNULL(aimInputs, aimInfo, status);

    if (aimInputs[inSurface_Mesh-1].nullVal == IsNull &&
        aimInputs[inMesh_Morph-1].vals.integer == (int) false) {
      AIM_ANALYSISIN_ERROR(aimInfo, inSurface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intent, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf("flightstreamAIM/aimUpdateState -> numBody=%d\n", numBody);
#endif

    if ((numBody <= 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "No body\n");
        status = CAPS_SOURCEERR;
        goto cleanup;

    } else if (numBody > 1) {
        AIM_ERROR(aimInfo, "FlightStream can only accept a single body\n");
        status = CAPS_SOURCEERR;
        goto cleanup;
    }

    if (aimInputs[inFluid_Properties-1].nullVal == aimInputs[inAltitude-1].nullVal) {
      AIM_ERROR(aimInfo, "Only one of 'Fluid_Properties' or 'Altitude' must be specified");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    if (aimInputs[inFluid_Properties-1].nullVal == NotNull) {
      status = _getFluid_Properties(aimInfo,
                                    aimInputs[inFluid_Properties-1].length,
                                    aimInputs[inFluid_Properties-1].vals.tuple,
                                    flightstreamInstance);
      AIM_STATUS(aimInfo, status);
    } else {
      status = aim_getBounds(aimInfo, &nbound, &boundNames);
      AIM_STATUS(aimInfo, status);
      if (nbound > 0) {
        AIM_ERROR(aimInfo, "Data transfer requires 'Fluid_Properties'");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = aim_capsLength(aimInfo, &lengthUnits);
    AIM_NOTFOUND(aimInfo, status);

    if (status == CAPS_NOTFOUND) {
        AIM_ERROR(aimInfo, "capsLength attribute must be specified for FlightStream");
        goto cleanup;
    }
    AIM_NOTNULL(lengthUnits, aimInfo, status);

    status = aim_convert(aimInfo, 1, lengthUnits, &scaleFactor, units->length, &scaleFactor);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf("flightstreamAIM/aimUpdateState -> scaleFactor=%f\n", scaleFactor);
#endif

    // Loop over bodies and look for reference quantity attributes
    for (i=0; i < numBody; i++) {
        status = EG_attributeRet(bodies[i], "capsReferenceChord",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {
            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Cref = reals[0] * scaleFactor;
                foundCref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceChord should be a scalar\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        status = EG_attributeRet(bodies[i], "capsReferenceSpan",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {
            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Bref = reals[0] * scaleFactor;
                foundBref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceSpan should be a scalar\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        status = EG_attributeRet(bodies[i], "capsReferenceArea",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {
            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Sref = reals[0] * scaleFactor * scaleFactor;
                foundSref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceArea should be a scalar\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        status = EG_attributeRet(bodies[i], "capsReferenceX",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {

            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Xref = reals[0] * scaleFactor;
                foundXref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceX should be followed by a single real value!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        status = EG_attributeRet(bodies[i], "capsReferenceY",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {

            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Yref = reals[0] * scaleFactor;
                foundYref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceY should be followed by a single real value!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        status = EG_attributeRet(bodies[i], "capsReferenceZ",
                                 &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS){

            if (atype == ATTRREAL && alen == 1) {
                flightstreamInstance->Zref = reals[0] * scaleFactor;
                foundZref = (int)true;
            } else {
                AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }
    }

    if (aimInputs[inReferenceArea-1].nullVal == NotNull) {
        flightstreamInstance->Sref = aimInputs[inReferenceArea-1].vals.real;
        foundSref = (int)true;
    }
    if (foundSref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceArea is not set on any body and 'ReferenceArea' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceChord-1].nullVal == NotNull) {
        flightstreamInstance->Cref = aimInputs[inReferenceChord-1].vals.real;
        foundCref = (int)true;
    }
    if (foundCref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceChord is not set on any body and 'ReferenceChord' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceSpan-1].nullVal == NotNull) {
        flightstreamInstance->Bref = aimInputs[inReferenceSpan-1].vals.real;
        foundBref = (int)true;
    }
    if (foundBref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceSpan is not set on any body and 'ReferenceSpan' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceX-1].nullVal == NotNull) {
        flightstreamInstance->Xref = aimInputs[inReferenceX-1].vals.real;
        foundXref = (int)true;
    }
    if (foundXref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceX is not set on any body and 'ReferenceX' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceY-1].nullVal == NotNull) {
        flightstreamInstance->Yref = aimInputs[inReferenceY-1].vals.real;
        foundYref = (int)true;
    }
    if (foundYref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceY is not set on any body and 'ReferenceY' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceZ-1].nullVal == NotNull) {
        flightstreamInstance->Zref = aimInputs[inReferenceZ-1].vals.real;
        foundZref = (int)true;
    }
    if (foundZref == (int)false) {
        AIM_ERROR(aimInfo, "capsReferenceZ is not set on any body and 'ReferenceZ' input not set!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[inReferenceVelocity-1].nullVal == IsNull) {
      AIM_ERROR(aimInfo, "Input ReferenceVelocity not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get mesh
    flightstreamInstance->meshRefIn = (aimMeshRef *) aimInputs[inSurface_Mesh-1].vals.AIMptr;

    if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
        flightstreamInstance->meshRefIn == NULL) { // If we are mighty morphing

      // Lets "load" the meshRef now since it's not linked
      status = aim_loadMeshRef(aimInfo, &flightstreamInstance->meshRefMorph);
      AIM_STATUS(aimInfo, status);

      // Mighty Morph the mesh
      status = aim_morphMeshUpdate(aimInfo, &flightstreamInstance->meshRefMorph, numBody, bodies);
      AIM_STATUS(aimInfo, status);

      /*@-immediatetrans@*/
      flightstreamInstance->meshRefIn = &flightstreamInstance->meshRefMorph;
      /*@+immediatetrans@*/
    }
    AIM_NOTNULL(flightstreamInstance->meshRefIn, aimInfo, status);

    // Get attribute to index mapping
    status = create_MeshRefToIndexMap(aimInfo, flightstreamInstance->meshRefIn, &flightstreamInstance->groupMap);
    AIM_STATUS(aimInfo, status);

#if 0
    status = cfd_cfdCoefficientUnits(aimInfo,
                                     flightstreamInstance->Cref, units->length,
                                     flightstreamInstance->Sref, units->area,
                                     aimInputs[Freestream_Density-1 ].vals.real, aimInputs[Freestream_Density-1 ].units,
                                     aimInputs[inReferenceVelocity-1].vals.real, aimInputs[inReferenceVelocity-1].units,
                                     aimInputs[Freestream_Pressure-1].vals.real, aimInputs[Freestream_Pressure-1].units,
                                     units);
    AIM_STATUS(aimInfo, status);
#endif

#ifdef DEBUG
    printf("flightstreamInstance->projectName = %s\n", flightstreamInstance->projectName);
    printf("flightstreamInstance->Cref        = %f\n", flightstreamInstance->Cref       );
    printf("flightstreamInstance->Bref        = %f\n", flightstreamInstance->Bref       );
    printf("flightstreamInstance->Sref        = %f\n", flightstreamInstance->Sref       );
    printf("flightstreamInstance->Xref        = %f\n", flightstreamInstance->Xref       );
    printf("flightstreamInstance->Yref        = %f\n", flightstreamInstance->Yref       );
    printf("flightstreamInstance->Zref        = %f\n", flightstreamInstance->Zref       );
#endif

cleanup:
    AIM_FREE(boundNames);
    return status;
}


/**********************************************************************/
/* aimPreAnalysis - generate FlightStream input file                  */
/**********************************************************************/

int
aimPreAnalysis(const void *instStore,   /* (in)  AIM instance storage */
               void       *aimInfo,     /* (in)  AIM context */
               capsValue  aimInputs[])  /* (in)  array of analysis inputs */
{
    // Function return flag
    int status = CAPS_SUCCESS;

    //double alt, temp, vel;
    FILE   *fp = NULL;

    int i;
    int nbound = 0;
    char **boundNames = NULL;

    const char *export = NULL;
    const char *flightscript = NULL;

    const aimStorage *flightstreamInstance = (const aimStorage *)instStore;
    const cfdUnitsStruct *units = &flightstreamInstance->units;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimPreAnalysis()\n");
#endif

    AIM_NOTNULL(aimInputs, aimInfo, status);

    if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
         aimInputs[inSurface_Mesh-1].nullVal == NotNull) { // If we are mighty morphing
        // store the current mesh for future iterations
        status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[inSurface_Mesh-1].vals.AIMptr, NULL);
        AIM_STATUS(aimInfo, status);
    }

    status = fs_dataTransfer(aimInfo, aimInputs, flightstreamInstance);
    AIM_STATUS(aimInfo, status);

    /*
    // compute velocity from standard atmosphere
    alt = 0.3048 * aimInputs[inAltitude-1].vals.real;
    if (alt <= 11000) {
        temp = 288.16 - 0.0065 * alt;
    } else if (alt <= 25000) {
        temp = 216.66;
    } else if (alt < 47000) {
        temp = 216.66 + 0.0030 * (alt - 25000);
    } else if (alt < 53000) {
        temp = 282.66;
    } else if (alt < 79000) {
        temp = 282.66 - 0.0045 * (alt - 53000);
    } else if (alt < 90000) {
        temp = 165.66;
    } else {
        temp = 165.66 + 0.0040 * (alt - 90000);
    }
    vel = aimInputs[inMach-1].vals.real * sqrt(1.40 * 287 * temp);
    */

    // write flightstream script file
    fp = aim_fopen(aimInfo, "script.txt", "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open \"script.txt\"");
        status = CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp, "SOLVER_CLEAR\n\n");

    fprintf(fp, "IMPORT\n");
    fprintf(fp, "UNITS        %s\n", units->length);
    fprintf(fp, "FILE_TYPE    OBJ\n");
    fprintf(fp, "FILE         %s%s\n", flightstreamInstance->meshRefIn->fileName, MESHEXTENSION);
    fprintf(fp, "CLEAR\n\n");


    // create a new csys for the moment reference (X-streamwise, Y-left, Z-up)
    fprintf(fp, "CREATE_NEW_COORDINATE_SYSTEM\n");
    fprintf(fp, "EDIT_COORDINATE_SYSTEM\n");
    fprintf(fp, "FRAME 2\n");
    fprintf(fp, "NAME Moment_Reference\n");
    fprintf(fp, "ORIGIN_X %15.5f %s\n", flightstreamInstance->Xref, units->length);
    fprintf(fp, "ORIGIN_Y %15.5f %s\n", flightstreamInstance->Yref, units->length);
    fprintf(fp, "ORIGIN_Z %15.5f %s\n", flightstreamInstance->Zref, units->length);
    fprintf(fp, "VECTOR_X_X 1\n");
    fprintf(fp, "VECTOR_X_Y 0\n");
    fprintf(fp, "VECTOR_X_Z 0\n");
    fprintf(fp, "VECTOR_Y_X 0\n");
    fprintf(fp, "VECTOR_Y_Y 1\n");
    fprintf(fp, "VECTOR_Y_Z 0\n");
    fprintf(fp, "VECTOR_Z_X 0\n");
    fprintf(fp, "VECTOR_Z_Y 0\n");
    fprintf(fp, "VECTOR_Z_Z 1\n");
    fprintf(fp, "SET_SOLVER_ANALYSIS_LOADS_FRAME 2\n");

    fprintf(fp, "AUTO_DETECT_BASE_REGIONS\n\n");

    //fprintf(fp, "PHYSICS\n");
    //fprintf(fp, "AUTO_TRAIL_EDGES\n");
    //fprintf(fp, "AUTO_WAKE_NODES\n");
    //fprintf(fp, "END\n\n");

    if (aimInputs[inFluid_Properties-1].nullVal == NotNull) {
      fprintf(fp, "FLUID_PROPERTIES\n");
      fprintf(fp, "DENSITY                     %15.12e\n",   flightstreamInstance->rhoref);
      fprintf(fp, "PRESSURE                    %15.12e\n",   flightstreamInstance->pref);
      fprintf(fp, "SONIC_VELOCITY              %15.12e\n",   flightstreamInstance->sonicspeedref);
      fprintf(fp, "TEMPERATURE                 %15.12e\n",   flightstreamInstance->tref);
      fprintf(fp, "VISCOSITY                   %15.12e\n",   flightstreamInstance->muref);
    } else {
      fprintf(fp, "AIR_ALTITUDE                %15.12e\n",   aimInputs[inAltitude-1].vals.real);
    }
    fprintf(fp, "SOLVER_SET_AOA              %15.12e\n",   aimInputs[inAlpha   -1].vals.real);
    fprintf(fp, "SOLVER_SET_SIDESLIP         %15.12e\n",   aimInputs[inBeta    -1].vals.real);
    fprintf(fp, "SOLVER_SET_MACH_NUMBER      %15.12e\n",   aimInputs[inMach    -1].vals.real);
    fprintf(fp, "#SOLVER_SET_REF_MACH_NUMBER  %15.12e\n",   aimInputs[inMach    -1].vals.real);
    // SOLVER_SET_REF_MACH_NUMBER does not currently work in FlightStream
    fprintf(fp, "SOLVER_SET_REF_VELOCITY     %15.12e\n",   aimInputs[inReferenceVelocity-1].vals.real);
    fprintf(fp, "SOLVER_SET_REF_AREA         %15.12e\n",   flightstreamInstance->Sref);
    fprintf(fp, "SOLVER_SET_REF_LENGTH       %15.12e\n\n", flightstreamInstance->Cref);

    fprintf(fp, "SET_SOLVER_MODEL            HYPERSONIC\n\n");

    fprintf(fp, "INITIALIZE_SOLVER\n");
//    fprintf(fp, "SURFACES                    -1\n");
    fprintf(fp, "SURFACES                    %d\n", flightstreamInstance->groupMap.numAttribute);
    for (i = 0; i < flightstreamInstance->groupMap.numAttribute; i++)
      fprintf(fp, "%d,0,DISABLE\n", i+1);
    fprintf(fp, "WAKE_TERMINATION_X          DEFAULT\n");
    fprintf(fp, "SYMMETRY_TYPE               NONE\n");
    fprintf(fp, "SYMMETRY_PERIODICITY        1\n");
    fprintf(fp, "LOAD_FRAME                  1\n");
    fprintf(fp, "PROXIMITY_AVOIDANCE         DISABLE\n");
    fprintf(fp, "STABILIZATION               ENABLE\n");
    fprintf(fp, "STABILIZATION_STRENGTH      1.0\n");
    fprintf(fp, "FAST_MULTIPOLE              ENABLE\n\n");

    fprintf(fp, "START_SOLVER\n\n");

    fprintf(fp, "EXPORT_SOLVER_ANALYSIS_SPREADSHEET\n");
    fprintf(fp, "results.txt\n\n");

    fprintf(fp, "EXPORT_LOG\n");
    fprintf(fp, "log.txt\n\n");

    // add custom exports
    if (aimInputs[inExport_Solver_Analysis-1].nullVal == NotNull) {
      export = aimInputs[inExport_Solver_Analysis-1].vals.string;
      for (i = 0; i < aimInputs[inExport_Solver_Analysis-1].length; i++) {
        if (strcasecmp(export, "TECPLOT") == 0) {
          fprintf(fp, "EXPORT_SOLVER_ANALYSIS_TECPLOT\n");
          fprintf(fp, "%s.tec\n\n", aimInputs[inProj_Name-1].vals.string);
        } else if (strcasecmp(export, "VTK") == 0) {
          fprintf(fp, "EXPORT_SOLVER_ANALYSIS_VTK\n");
          fprintf(fp, "%s.vtk\n\n", aimInputs[inProj_Name-1].vals.string);
        } else if (strcasecmp(export, "BDF") == 0) {
          fprintf(fp, "EXPORT_SOLVER_ANALYSIS_PLOAD_BDF\n");
          fprintf(fp, "%s.bdf\n\n", aimInputs[inProj_Name-1].vals.string);
        } else {
          AIM_ERROR(aimInfo, "Unknown Export_Solver_Analysis = %s", export);
          status = CAPS_BADVALUE;
          goto cleanup;
        }
        export += strlen(export)+1;
      }
    }

    status = aim_getBounds(aimInfo, &nbound, &boundNames);
    AIM_STATUS(aimInfo, status);

    if (nbound > 0) {
#ifdef ANALYSIS_PLOAD_BDF
      status = aim_rmFile(aimInfo, loadbdf);
      AIM_STATUS(aimInfo, status);

      fprintf(fp, "EXPORT_SOLVER_ANALYSIS_PLOAD_BDF\n");
      fprintf(fp, "%s\n", loadbdf);
      fprintf(fp, "SURFACES -1\n\n");
#elif defined(ANALYSIS_VTK)
      status = aim_rmFile(aimInfo, cpvtk);
      AIM_STATUS(aimInfo, status);

      fprintf(fp, "EXPORT_SOLVER_ANALYSIS_VTK\n");
      fprintf(fp, "%s\n", cpvtk);
      fprintf(fp, "SURFACES -1\n\n");
      fprintf(fp, "SET_VTK_EXPORT_VARIABLES 1 DISABLE\n");
      fprintf(fp, "CP\n\n");
#elif defined(ANALYSIS_CSV)
      status = aim_rmFile(aimInfo, dpcsv);
      AIM_STATUS(aimInfo, status);

      fprintf(fp, "EXPORT_SOLVER_ANALYSIS_CSV\n");
      fprintf(fp, "%s\n", dpcsv);
      fprintf(fp, "FORMAT DIFFERENCE-PRESSURE\n");
      fprintf(fp, "UNITS PASCALS\n");
      fprintf(fp, "SURFACES -1\n\n");
#else
#error "must define ANALYSIS_PLOAD_BDF or ANALYSIS_CSV"
#endif
    }

    if (aimInputs[inFlightScript-1].nullVal == NotNull) {
      flightscript = aimInputs[inFlightScript-1].vals.string;
      for (i = 0; i < aimInputs[inFlightScript-1].length; i++) {
        fprintf(fp, "%s\n", flightscript);
        flightscript += strlen(flightscript)+1;
      }
    }

    fprintf(fp, "CLOSE_FLIGHTSTREAM\n");

    fclose(fp); fp = NULL;

cleanup:
    AIM_FREE(boundNames);
    if (fp != NULL) fclose(fp);
    return status;
}


/**********************************************************************/
/* aimExecute - execute FlightStream                                  */
/**********************************************************************/
int aimExecute(/*@unused@*/ const void *instStore, void *aimInfo,
               int *state)
{
  /*! \page aimExecuteFLIGHTSTREAM AIM Execution
   *
   * If auto execution is enabled when creating an FlightStream AIM,
   * the AIM will execute FlightStream just-in-time on Linux with the command line:
   *
   * \code{.sh}
   * FlightStream script.txt > flightstreamOut.txt
   * \endcode
   *
   * and on Windows with the command:
   *
   * \code{.sh}
   * FlightStream -hidden -script script.txt > flightstreamOut.txt
   * \endcode
   *
   * In both cases the FlightStream executable is assumed to in the PATH environment variable.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an FlightStream AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or FlightStream can be executed by calling preAnalysis, system, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * flightstream.preAnalysis()
   *
   * print ("\n\nRunning......")
   * flightstream.system("FlightStream.exe -hidden -script script.txt"); # Run via system call in inputs analysis directory
   *
   * print ("\n\postAnalysis......")
   * flightstream.postAnalysis()
   * \endcode
   */
  int status = CAPS_SUCCESS;
  capsValue *flightstream = NULL;

#define fslen PATH_MAX+42
  char fscmd[fslen];

  *state = 0;

  aim_getValue(aimInfo, inFlightStream, ANALYSISIN, &flightstream);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(flightstream, aimInfo, status);

#ifdef WIN32
  snprintf(fscmd, fslen, "%s -hidden -script script.txt > flightstreamOut.txt", flightstream->vals.string);
#else
  snprintf(fscmd, fslen, "%s script.txt", flightstream->vals.string);
#endif
  printf(" Executing: %s\n", fscmd);
  status = aim_system(aimInfo, "", fscmd);
  AIM_STATUS(aimInfo, status);

#undef fslen
cleanup:
  return status;
}


/**********************************************************************/
/* aimPostAnalysis - read FlightStream output file                    */
/**********************************************************************/

int
aimPostAnalysis(/*@unused@*/ void *instStore,    /* (in)  AIM instance storage */
                /*@unused@*/ void *aimInfo,      /* (in)  AIM context */
                /*@unused@*/ int restart,        /* (in)  0=normal, 1=restart */
                /*@unused@*/ capsValue inputs[]) /* (in)  array of analysis inputs */
{
    int    status = CAPS_SUCCESS;

    size_t  linecap=0;
    char    *line=NULL;

    FILE    *fp=NULL;
    aimStorage *flightstreamInstance = (aimStorage *)instStore;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimPostAnalysis()\n");
#endif

    fp = aim_fopen(aimInfo, "results.txt", "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open \"results.txt\"\n");
        status = CAPS_IOERR;
        goto cleanup;
    }

    /* read the file until we find a line that starts with "     Total" */
    while (1) {
        status = getline(&line, &linecap, fp);
        if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));

        AIM_NOTNULL(line, aimInfo, status);

        if (strncmp(line, "     Total", 10) == 0) {

            /* now that we have found the line, get the various coefficients */
            sscanf(&line[13], "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &(flightstreamInstance->Cx),
                   &(flightstreamInstance->Cy),
                   &(flightstreamInstance->Cz),
                   &(flightstreamInstance->CL),
                   &(flightstreamInstance->CDi),
                   &(flightstreamInstance->CDo),
                   &(flightstreamInstance->CMx),
                   &(flightstreamInstance->CMy),
                   &(flightstreamInstance->CMz));

            status = CAPS_SUCCESS;
            fclose(fp);
            fp = NULL;
            goto cleanup;
        }
    }

    /* getting here means we did not find results */
    AIM_ERROR(aimInfo, "Did not find \"Total\" results\n");
    status = CAPS_IOERR;

cleanup:
    if (fp   != NULL) fclose(fp);
    if (line != NULL) free(line);

    return status;
}


/**********************************************************************/
/* aimCalcOutput - retreive FlightStream output information           */
/**********************************************************************/

int
aimCalcOutput(/*@unused@*/ void *instStore,    /* (in)  AIM instance storage */
              /*@unused@*/ void *aimInfo,      /* (in)  AIM context */
              /*@unused@*/ int index,          /* (in)  analysis output */
              /*@unused@*/ capsValue *outval)  /* (in)  pointer to capsValue to fill */
{
    int status = CAPS_SUCCESS;

    const   aimStorage *flightstreamInstance = (const aimStorage *)instStore;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimCalcOutput(index=%d)\n", index);
#endif

    if        (index == outCx) {
        outval->vals.real = flightstreamInstance->Cx;
    } else if (index == outCy) {
        outval->vals.real = flightstreamInstance->Cy;
    } else if (index == outCz) {
        outval->vals.real = flightstreamInstance->Cz;
    } else if (index == outCL) {
        outval->vals.real = flightstreamInstance->CL;
    } else if (index == outCDi) {
        outval->vals.real = flightstreamInstance->CDi;
    } else if (index == outCDo) {
        outval->vals.real = flightstreamInstance->CDo;
    } else if (index == outCMx) {
        outval->vals.real = flightstreamInstance->CMx;
    } else if (index == outCMy) {
        outval->vals.real = flightstreamInstance->CMy;
    } else if (index == outCMz) {
        outval->vals.real = flightstreamInstance->CMz;
    } else {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

cleanup:
    return status;
}


/**********************************************************************/
/* aimCleanup - free up memory allocated in aimInitialize             */
/**********************************************************************/

void
aimCleanup(void *instStore)             /* (in)  AIM instance storage */
{

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("flightstreamAIM/aimCleanup()\n");
#endif

    aimStorage *flightstreamInstance = (aimStorage *)instStore;

    destroy_aimStorage(flightstreamInstance);

    AIM_FREE(flightstreamInstance);
}


// ********************** AIM Function Break *****************************
int
aimDiscr(char *tname, capsDiscr *discr)
{

  int i; // Indexing

  int status; // Function return status

  int numBody;

  // EGADS objects
  ego *bodies = NULL, *tess = NULL;

  const char   *intents;

  // Volume Mesh obtained from meshing AIM
  const aimMeshRef *meshRef;

  aimStorage *flightstreamInstance;

  flightstreamInstance = (aimStorage *) discr->instStore;

#ifdef DEBUG
  printf(" flightstreamAIM/aimDiscr: tname = %s!\n", tname);
#endif

  if (tname == NULL) return CAPS_NOTFOUND;

  // Currently this ONLY works if the capsTranfer lives on single body!
  status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
  if (status != CAPS_SUCCESS) {
    printf(" flightstreamAIM/aimDiscr: aim_getBodies = %d!\n", status);
    return status;
  }
  if (bodies == NULL) {
    AIM_ERROR(discr->aInfo, "NULL Bodies!\n");
    return CAPS_NULLOBJ;
  }

  // Get mesh
  meshRef = flightstreamInstance->meshRefIn;
  AIM_NOTNULL(meshRef, discr->aInfo, status);

  if (aim_newGeometry(discr->aInfo) == CAPS_SUCCESS) {
    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            1, // Only search down to the face level of the EGADS body
                                            &flightstreamInstance->groupMap);
    AIM_STATUS(discr->aInfo, status);
  }

  // Lets check the volume mesh

  // Do we have an individual surface mesh for each body
  if (meshRef->nmap != numBody) {
    AIM_ERROR(  discr->aInfo, "Number of surface mesh in the linked volume mesh (%d) does not match the number");
    AIM_ADDLINE(discr->aInfo,"of bodies (%d) - data transfer is NOT possible.", meshRef->nmap,numBody);
    status = CAPS_MISMATCH;
    goto cleanup;
  }

  // Lets store away our tessellation now
  AIM_ALLOC(tess, meshRef->nmap, ego, discr->aInfo, status);
  for (i = 0; i < meshRef->nmap; i++) {
    tess[i] = meshRef->maps[i].tess;
  }

  status = mesh_fillDiscr(tname, &flightstreamInstance->groupMap, meshRef->nmap, tess, discr);
  AIM_STATUS(discr->aInfo, status);

#ifdef DEBUG
  printf(" flightstreamAIM/aimDiscr: Finished!!\n");
#endif

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(tess);
  return status;
}

// ********************** AIM Function Break *****************************
void aimFreeDiscrPtr(void *ptr)
{
  // Extra information to store into the discr void pointer - just a int array
  EG_free(ptr);
}

// ********************** AIM Function Break *****************************
int
aimLocateElement(capsDiscr *discr, double *params, double *param,
                 int *bIndex, int *eIndex, double *bary)
{
    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


// ********************** AIM Function Break *****************************
int
aimTransfer(capsDiscr *discr, const char *name, int npts, int rank, double *data,
            char **units)
{
  /*! \page dataTransferFLIGHTSTREAM AIM Data Transfer
   *
   * The FlightStream AIM has the ability to transfer surface data (e.g. pressure distributions) to and from the AIM
   * using the conservative and interpolative data transfer schemes in CAPS.
   *
   * \section dataFromFlightsteram Data transfer from FlightStream (FieldOut)
   *
   * <ul>
   * <li> <B>"Pressure" </B> </li> <br>
   *  Loads the pressure distribution from FlightStream vtk file.
   *  This distribution may be scaled based on
   *  Pressure = Pressure_Scale_Factor*Pressure, where "Pressure_Scale_Factor"
   *  is an AIM input (\ref aimInputsFLIGHTSTREAM)
   * </ul>
   *
   */
  int    i, j, global, status, bIndex;
  double **rvec=NULL, scale = 1.0;
  capsValue *Pressure_Scale_Factor_Value=NULL;
  int state, nglobal;
  ego        body;
  FILE *fp = NULL;
#ifdef ANALYSIS_PLOAD_BDF
  int ID, inode, jnode;
  char str[10];
  size_t linecap = 0;
  char *line = NULL; // Temporary line holder
  const char *PLOAD4 = "$$  PLOAD4 Data";
#elif defined(ANALYSIS_VTK)
  size_t linecap = 0;
  char *line = NULL; // Temporary line holder
  const char *CP = "SCALARS Cp FLOAT";
  double vref = 0, qref=0;
  capsValue *ReferenceVelocity=NULL;
#elif defined(ANALYSIS_CSV)
  double x, y, z;
#else
#error "must define ANALYSIS_PLOAD_BDF, ANALYSIS_VTK, or, ANALYSIS_CSV"
#endif

  const aimStorage *flightstreamInstance = (const aimStorage *)discr->instStore;
  const aimMeshRef *meshRef = flightstreamInstance->meshRefIn;

#ifdef DEBUG
  printf(" flightstreamAIM/aimTransfer name = %s  npts = %d/%d!\n",
         name, npts, len_wrt);
#endif

  if (strcmp(name, "Pressure") == 0) {

    if (rank != 1) {
      AIM_ERROR(discr->aInfo, "Rank (%d) must be 1 for 'Pressure'!", rank);
      status = CAPS_NOTIMPLEMENT;
      goto cleanup;
    }

#ifdef ANALYSIS_PLOAD_BDF
    fp = aim_fopen(discr->aInfo, loadbdf, "r");
    if (fp == NULL) {
      AIM_ERROR(discr->aInfo, "Cannot open \"%s\"", loadbdf);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* try and read the bdf file */
    while (getline(&line, &linecap, fp) > 0) {
      AIM_NOTNULL(line, discr->aInfo, status);
      if (strncmp(line, PLOAD4, strlen(PLOAD4)) == 0) {
        break;
      }
    }

    AIM_NOTNULL(line, discr->aInfo, status);
    if (strncmp(line, PLOAD4, strlen(PLOAD4)) != 0) {
      AIM_ERROR(discr->aInfo, "Could not find 'PLOAD4' data in \"%s\"", loadbdf);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // Skip '$$'
    status = getline(&line, &linecap, fp);
    if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));

    AIM_ALLOC(rvec, meshRef->nmap, double*, discr->aInfo, status);
    for (i = 0; i < meshRef->nmap; i++) {
      rvec[i] = NULL;
    }

    jnode = 1;
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(discr->aInfo, status);
      AIM_ALLOC(rvec[i], nglobal, double, discr->aInfo, status);

      for (j = 0; j < nglobal; j++, jnode++) {
        status = getline(&line, &linecap, fp);
        if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));
        status = sscanf(line, "%s %d %d %lf", str, &ID, &inode, &rvec[i][j]);
        if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));
        if (inode != jnode) {
          AIM_ERROR(discr->aInfo, "While reading %s, %d != %d!", loadbdf, inode, jnode);
          status = CAPS_NOTIMPLEMENT;
          goto cleanup;
        }
        // subtract off the reference pressure
        // both BDF pressure and pref are in Pa
        rvec[i][j] -= flightstreamInstance->pref;
      }
    }

    // Convert from Pa to working pressure
    scale = 1;
    status = aim_convert(discr->aInfo, 1, "Pa", &scale, flightstreamInstance->units.pressure, &scale);
    AIM_STATUS(discr->aInfo, status);

#elif defined(ANALYSIS_VTK)
    fp = aim_fopen(discr->aInfo, cpvtk, "r");
    if (fp == NULL) {
      AIM_ERROR(discr->aInfo, "Cannot open \"%s\"", cpvtk);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* try and read the vtk file */
    while (getline(&line, &linecap, fp) > 0) {
      AIM_NOTNULL(line, discr->aInfo, status);
      if (strncmp(line, CP, strlen(CP)) == 0) {
        break;
      }
    }

    AIM_NOTNULL(line, discr->aInfo, status);
    if (strncmp(line, CP, strlen(CP)) != 0) {
      AIM_ERROR(discr->aInfo, "Could not find 'PLOAD4' data in \"%s\"", cpvtk);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // Skip 'LOOKUP_TABLE default'
    status = getline(&line, &linecap, fp);
    if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));

    AIM_ALLOC(rvec, meshRef->nmap, double*, discr->aInfo, status);
    for (i = 0; i < meshRef->nmap; i++) {
      rvec[i] = NULL;
    }

    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(discr->aInfo, status);
      AIM_ALLOC(rvec[i], nglobal, double, discr->aInfo, status);

      for (j = 0; j < nglobal; j++) {
        status = getline(&line, &linecap, fp);
        if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));
        status = sscanf(line, "%lf", &rvec[i][j]);
        if (status <= 0) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));
      }
    }

    status = aim_getValue(discr->aInfo, inReferenceVelocity, ANALYSISIN, &ReferenceVelocity);
    AIM_STATUS(discr->aInfo, status);
    AIM_NOTNULL(ReferenceVelocity, discr->aInfo, status);

    // Compute dynamic pressure
    status = aim_convert(discr->aInfo, 1,
                         ReferenceVelocity->units,  &ReferenceVelocity->vals.real,
                         "m/s", &vref);
    AIM_STATUS(discr->aInfo, status);

    qref = 0.5 * flightstreamInstance->rhoref * vref*vref;

    status = aim_convert(discr->aInfo, 1, "Pa", &qref, flightstreamInstance->units.pressure, &qref);
    AIM_STATUS(discr->aInfo, status);
    scale *= qref;

#elif defined(ANALYSIS_CSV)

    fp = aim_fopen(discr->aInfo, dpcsv, "r");
    if (fp == NULL) {
      AIM_ERROR(discr->aInfo, "Cannot open \"%s\"", dpcsv);
      status = CAPS_IOERR;
      goto cleanup;
    }

    AIM_ALLOC(rvec, meshRef->nmap, double*, discr->aInfo, status);
    for (i = 0; i < meshRef->nmap; i++) {
      rvec[i] = NULL;
    }

    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(discr->aInfo, status);
      AIM_ALLOC(rvec[i], nglobal, double, discr->aInfo, status);

      for (j = 0; j < nglobal; j++) {
        status = fscanf(fp, "%lf, %lf, %lf, %lf\n", &x, &y, &z, &rvec[i][j]);
        if (status != 4) AIM_STATUS(discr->aInfo, (status = CAPS_IOERR));
      }
    }
#else
#error "must define ANALYSIS_PLOAD_BDF, ANALYSIS_VTK, or ANALYSIS_CSV"
#endif


    // Custom additional scale factor
    status = aim_getValue(discr->aInfo, inPressure_Scale_Factor, ANALYSISIN, &Pressure_Scale_Factor_Value);
    AIM_STATUS(discr->aInfo, status);
    AIM_NOTNULL(Pressure_Scale_Factor_Value, discr->aInfo, status);
    scale *= Pressure_Scale_Factor_Value->vals.real;

    // set the units
    AIM_STRDUP(*units, flightstreamInstance->units.pressure, discr->aInfo, status);
  }

  /* move the appropriate parts of the tessellation to data */
  AIM_NOTNULL(rvec, discr->aInfo, status);
  for (i = 0; i < npts; i++) {
    /* points might span multiple bodies */
    bIndex = discr->tessGlobal[2*i  ];
    global = discr->tessGlobal[2*i+1];
    for (j = 0; j < rank; j++)
      data[rank*i+j] = rvec[bIndex-1][rank*(global-1)+j] * scale;
  }

  status = CAPS_SUCCESS;
cleanup:
  if (rvec != NULL) {
    for (i = 0; i < meshRef->nmap; ++i) {
      AIM_FREE(rvec[i]);
    }
  }
  AIM_FREE(rvec);
#if defined(ANALYSIS_PLOAD_BDF) || defined(ANALYSIS_VTK)
  if (line != NULL) free(line);
#endif

  if (fp != NULL) fclose(fp);

  return status;
}


// ********************** AIM Function Break *****************************
int
aimInterpolation(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                 int eIndex, double *bary, int rank, double *data,
                 double *result)
{
#ifdef DEBUG
    printf(" flightstreamAIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


// ********************** AIM Function Break *****************************
int
aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                  int eIndex, double *bary, int rank, double *r_bar,
                  double *d_bar)
{
#ifdef DEBUG
    printf(" flightstreamAIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


// ********************** AIM Function Break *****************************
int
aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
               int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" flightstreamAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank,
                           data, result);
}


// ********************** AIM Function Break *****************************
int
aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" flightstreamAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
