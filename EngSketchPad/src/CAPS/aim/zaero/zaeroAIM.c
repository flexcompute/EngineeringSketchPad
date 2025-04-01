/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             ZAERO AIM
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewZAERO ZAero AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with Zona's <a href="https://www.zonatech.com/zaero.html">ZAero</a>.
 * ZAero is designed to multiplie different aeroelastic analysis disciplines. Currently only a subset of ZAero's
 * input options have been exposed in the analysis interface module (AIM), but features can easily be included
 * as future needs arise.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsZAERO and \ref aimOutputsZAERO, respectively.
 *
 * Geometric attributes recognized by the AIM are provided in \ref attributeZAERO.
 *
 */

 /*! \page attributeZAERO Attribution
 *
 * The following list of attributes drives the ZAero geometric definition.
 *
 *  - <b> capsLength</b> This attribute defines the length units that the *.csm file is generated in.  ZAero grids
 *  MUST be in units of meter, as such the geometry is scaled accordingly based on this value.
 *
 *  - <b> capsReferenceArea</b>  [Optional] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the reference area in ZAero's input file with its units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  This value may be alternatively set
 * through an input value, "ReferenceArea" (see \ref aimInputsZAERO)
 *
 *  - <b> capsReferenceChord</b> and <b> capsReferenceSpan</b> [Optional] These attribute may exist on any <em> Body</em>. Their
 * value will be used as the reference moment lengths in ZAero's input file with their units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  These values may be alternatively set
 * through an input value, "Moment_Length" (see \ref aimInputsZAERO)
 *
 *  - <b> capsReferenceX</b>, <b> capsReferenceY</b>, and <b>capsReferenceZ</b> [Optional]
 * These attribute may exist on any <em> Body</em>. Their
 * value will be used as the center of gravity (CG) location in ZAero's input file with their units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  These values may be alternatively set
 * through an input value, "Moment_Center" (see \ref aimInputsZAERO)
 *
 */

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "aimUtil.h"
#include "miscUtils.h"
#include "jsonUtils.h"
#include "feaUtils.h"
#include "vlmUtils.h"
#include "cfdUtils.h"


#include "zaeroUtils.h"
#include "zaeroGeneral.h"
#include "zaeroDiscipline.h"
#include "zaeroGraphical.h"

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#define PATH_MAX   _MAX_PATH
#else
#include <unistd.h>
#include <limits.h>
#endif

//#define DEBUG

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


enum aimInputs
{
  inProj_Name = 1,                 /* index is 1-based */
  inAnalysis,
  inFile_Format,
  inFEM_1,
  inF06_1,
  inFEM_2,
  inF06_2,
  inCPU,
  inMemory,
  inSmart_Restart,
  inEcho,
  inOutput,
  inHFG,
  inUAIC,
  inSpline,
  inVLM_Surface,
  inVLM_Control,
  inTrim_Variable,
  inReferenceArea,
  inReferenceChord,
  inReferenceSpan,
  inMoment_Center,
  inMassPropLink,
  NUMINPUT = inMassPropLink       /* Total number of inputs */
};

enum aimOutputs
{
  //outBeta = 1,      /* index is 1-based */
  NUMOUTPUT = 0       /* Total number of outputs */
};


/* ZAERO storage */
typedef struct {

  // whether to use Smart Restart
  int smartRestart;

  // Attribute to index map
  mapAttrToIndexStruct attrMap;

  // Units structure
  cfdUnitsStruct units;

  // Output formating
  feaFileTypeEnum feaFormatType;

  // ZAERO problem data
  zaeroProblemStruct zaeroProblem;

  // ZAERO artifact filenames, relative to analysis dir.
  zaeroArtifactsStruct artifacts;

} aimStorage;

// AIM storage constructor
static int initiate_aimStorage(aimStorage *zaeroInstance) {

  // Set initial values for zaeroInstances

  zaeroInstance->smartRestart = 0;

  zaeroInstance->feaFormatType = SmallField;

  // Container for attribute to index map
  (void)initiate_mapAttrToIndexStruct(&zaeroInstance->attrMap);

  (void)initiate_cfdUnitsStruct(&zaeroInstance->units);

  (void)initiate_zaeroProblemStruct(&zaeroInstance->zaeroProblem);

  (void)initiate_zaeroArtifactsStruct(&zaeroInstance->artifacts);

  return CAPS_SUCCESS;
}

// AIM storage destructor
static int destroy_aimStorage(aimStorage *zaeroInstances) {

  int status;

  // Attribute to index map
  status = destroy_mapAttrToIndexStruct(&zaeroInstances->attrMap);
  if (status != CAPS_SUCCESS)
    printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

  // Destroy units
  status = destroy_cfdUnitsStruct(&zaeroInstances->units);
  if (status != CAPS_SUCCESS)
    printf("Error: Status %d during destroy_cfdUnitsStruct!\n", status);

  // Destroy ZAERO problem struct
  status = destroy_zaeroProblemStruct(&zaeroInstances->zaeroProblem);
  if (status != CAPS_SUCCESS)  {
      printf("Error: Status %d during destroy_zaeroProblemStruct!\n", status);
  }

  // Destroy ZAERO artifacts struct
  status = destroy_zaeroArtifactsStruct(&zaeroInstances->artifacts);
  if (status != CAPS_SUCCESS)  {
    printf("Error: Status %d during destroy_zaeroArtifactsStruct!\n", status);
  }

  initiate_aimStorage(zaeroInstances);

  return CAPS_SUCCESS;
}

static int
_createVLM(void *aimInfo,
           capsValue *aimInputs,
           aimStorage *zaeroInstance)
{
  int status, status2; // Function return status

  int i, j, k, surfaceIndex = 0, sectionIndex; // Indexing

  // Bodies
  const char *intents;
  int   numBody; // Number of Bodies
  ego  *bodies;

  int foundSref=(int)false, foundCref=(int)false, foundBref=(int)false;
  int foundXref=(int)false, foundYref=(int)false, foundZref=(int)false;

  const char *lengthUnits=NULL;
  double scaleFactor = 1.0;

  // EGADS return values
  int          atype, alen;
  const int    *ints;
  const char   *string;
  const double *reals;

  // Aeroelastic information
  int numVLMSurface = 0;
  vlmSurfaceStruct *vlmSurface = NULL;
  int numSpanWise;

  int numVLMControl = 0;
  vlmControlStruct *vlmControl = NULL;

  const char *Lunits=NULL;
  const cfdUnitsStruct *units = &zaeroInstance->units;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  if (status != CAPS_SUCCESS) return status;

#ifdef DEBUG
  printf(" nastranAIM/createVLMMesh  nbody = %d!\n", numBody);
#endif

  if ((numBody <= 0) || (bodies == NULL)) {
#ifdef DEBUG
    printf(" nastranAIM/createVLMMesh No Bodies!\n");
#endif
    return CAPS_SOURCEERR;
  }

  if (units->length != NULL)
    Lunits = units->length;
  else
    Lunits = "m";

  // Get aerodynamic reference quantities

  status = aim_capsLength(aimInfo, &lengthUnits);
  AIM_NOTFOUND(aimInfo, status);
  if (status == CAPS_NOTFOUND) {
    AIM_ERROR(aimInfo, "capsLength attribute must be specified for ZAero");
    goto cleanup;
  }
  AIM_NOTNULL(lengthUnits, aimInfo, status);

  status = aim_convert(aimInfo, 1, lengthUnits, &scaleFactor, Lunits, &scaleFactor);
  AIM_STATUS(aimInfo, status);

  // Loop over bodies and look for reference quantity attributes
  for (i=0; i < numBody; i++) {
    status = EG_attributeRet(bodies[i], "capsReferenceArea",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        zaeroInstance->zaeroProblem.hfg.refArea = reals[0] * scaleFactor * scaleFactor;
        foundSref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceArea should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceChord",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        zaeroInstance->zaeroProblem.hfg.refChord = reals[0] * scaleFactor;
        foundCref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceChord should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceSpan",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        zaeroInstance->zaeroProblem.hfg.refSpan = reals[0] * scaleFactor;
        foundBref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceSpan should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceX",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {

      if (atype == ATTRREAL && alen == 1) {
        zaeroInstance->zaeroProblem.hfg.refCenter[0] = reals[0] * scaleFactor;
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
        zaeroInstance->zaeroProblem.hfg.refCenter[1] = reals[0] * scaleFactor;
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
        zaeroInstance->zaeroProblem.hfg.refCenter[2] = reals[0] * scaleFactor;
        foundZref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }
  }

  if (aimInputs[inReferenceArea-1].nullVal == NotNull) {
    zaeroInstance->zaeroProblem.hfg.refArea = aimInputs[inReferenceArea-1].vals.real;
    foundSref = (int)true;
  }
  if (aimInputs[inReferenceChord-1].nullVal == NotNull) {
    zaeroInstance->zaeroProblem.hfg.refChord = aimInputs[inReferenceChord-1].vals.real;
    foundCref = (int)true;
  }
  if (aimInputs[inReferenceSpan-1].nullVal == NotNull) {
    zaeroInstance->zaeroProblem.hfg.refSpan = aimInputs[inReferenceSpan-1].vals.real;
    foundBref = (int)true;
  }

  if (foundSref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceArea is not set on any body and 'ReferenceArea' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (foundCref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceChord is not set on any body and 'ReferenceChord' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (foundBref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceSpan is not set on any body and 'ReferenceSpan' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Check for moment reference overwrites
  if (aimInputs[inMoment_Center-1].nullVal == NotNull) {

    zaeroInstance->zaeroProblem.hfg.refCenter[0] = aimInputs[inMoment_Center-1].vals.reals[0];
    zaeroInstance->zaeroProblem.hfg.refCenter[1] = aimInputs[inMoment_Center-1].vals.reals[1];
    zaeroInstance->zaeroProblem.hfg.refCenter[2] = aimInputs[inMoment_Center-1].vals.reals[2];
  } else {
    if (foundXref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceX is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    if (foundYref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceY is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    if (foundZref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceZ is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  }


  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, inVLM_Surface) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, inVLM_Control) == CAPS_SUCCESS) {

    // Cleanup Aero storage first
    if (zaeroInstance->zaeroProblem.feaAero != NULL) {
      for (i = 0; i < zaeroInstance->zaeroProblem.numAero; i++) {
        status = destroy_feaAeroStruct(&zaeroInstance->zaeroProblem.feaAero[i]);
        AIM_STATUS(aimInfo, status);
      }
      AIM_FREE(zaeroInstance->zaeroProblem.feaAero);
    }
    zaeroInstance->zaeroProblem.numAero = 0;


    // Get capsGroup name and index mapping
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            3, //>2 - search the body, faces, edges, and all the nodes
                                            &zaeroInstance->attrMap);
    AIM_STATUS(aimInfo, status);


    // Get VLM surface information
    if (aimInputs[inVLM_Surface-1].nullVal != IsNull) {

      status = get_vlmSurface(aimInfo,
                              aimInputs[inVLM_Surface-1].length,
                              aimInputs[inVLM_Surface-1].vals.tuple,
                              &zaeroInstance->attrMap,
                              0.0, // default Cspace
                              &numVLMSurface,
                              &vlmSurface);
      AIM_STATUS(aimInfo, status);

    } else {
      AIM_ERROR(aimInfo, "No VLM_Surface tuple specified\n");
      status = CAPS_NOTFOUND;
      goto cleanup;
    }

    // Get VLM control surface information
    if (aimInputs[inVLM_Control-1].nullVal == NotNull) {

      status = get_vlmControl(aimInfo,
                              aimInputs[inVLM_Control-1].length,
                              aimInputs[inVLM_Control-1].vals.tuple,
                              zaeroInstance->units.length == NULL ? NULL : "degree",
                              &numVLMControl,
                              &vlmControl);

      AIM_STATUS(aimInfo, status);
    }

    printf("\nGetting FEA vortex lattice mesh\n");

    status = vlm_getSections(aimInfo,
                             numBody,
                             bodies,
                             "Aerodynamic",
                             zaeroInstance->attrMap,
                             vlmGENERIC,
                             numVLMSurface,
                             &vlmSurface);
    AIM_STATUS(aimInfo, status);
    if (vlmSurface == NULL) {
      status = CAPS_NULLVALUE;
      goto cleanup;
    }

    for (i = 0; i < numVLMSurface; i++) {

      // Compute equal spacing
      if (vlmSurface[i].NspanTotal > 0)
        numSpanWise = vlmSurface[i].NspanTotal;
      else if (vlmSurface[i].NspanSection > 0)
        numSpanWise = (vlmSurface[i].numSection-1)*vlmSurface[i].NspanSection;
      else {
        AIM_ERROR(aimInfo  , "Only one of numSpanTotal and numSpanPerSection can be non-zero!");
        AIM_ADDLINE(aimInfo, "    numSpanTotal      = %d", vlmSurface[i].NspanTotal);
        AIM_ADDLINE(aimInfo, "    numSpanPerSection = %d", vlmSurface[i].NspanSection);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      status = vlm_equalSpaceSpanPanels(aimInfo, numSpanWise,
                                        vlmSurface[i].numSection,
                                        vlmSurface[i].vlmSection);
      AIM_STATUS(aimInfo, status);
    }

    // Split the surfaces that have more than 2 sections into a new surface
    for (i = 0; i < numVLMSurface; i++) {

      if (vlmSurface->numSection < 2) {
        AIM_ERROR(aimInfo, "Surface '%s' has less than two-sections!", vlmSurface[i].name);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      for (j = 0; j < vlmSurface[i].numSection-1; j++) {

        // Increment the number of Aero surfaces
        zaeroInstance->zaeroProblem.numAero += 1;

        surfaceIndex = zaeroInstance->zaeroProblem.numAero - 1;

        // Allocate
        AIM_REALL(zaeroInstance->zaeroProblem.feaAero, zaeroInstance->zaeroProblem.numAero, feaAeroStruct, aimInfo, status);

        // Initiate feaAeroStruct
        status = initiate_feaAeroStruct(&zaeroInstance->zaeroProblem.feaAero[surfaceIndex]);
        AIM_STATUS(aimInfo, status);

        // Get surface Name - copy from original surface
        AIM_STRDUP(zaeroInstance->zaeroProblem.feaAero[surfaceIndex].name, vlmSurface[i].name, aimInfo, status);

        // Get surface ID - Multiple by 1000 !!
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].surfaceID =
            1000*zaeroInstance->zaeroProblem.numAero;

        // ADD something for coordinate systems

        // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
        sectionIndex = vlmSurface[i].vlmSection[j].sectionIndex;

        // Populate vmlSurface structure
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.Cspace = vlmSurface[i].Cspace;
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.Sspace = vlmSurface[i].Sspace;

        // use the section span count for the sub-surface
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.NspanTotal = vlmSurface[i].vlmSection[sectionIndex].Nspan;
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.Nchord     = vlmSurface[i].Nchord;

        // Copy section information
        zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.numSection = 2;

        AIM_ALLOC(zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.vlmSection, 2, vlmSectionStruct, aimInfo, status);

        for (k = 0; k < zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.numSection; k++) {

          // Add k to section indexing variable j to get j and j+1 during iterations

          // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
          sectionIndex = vlmSurface[i].vlmSection[j+k].sectionIndex;

          status = initiate_vlmSectionStruct(&zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
          AIM_STATUS(aimInfo, status);

          // Copy the section data - This also copies the control data for the section
          status = copy_vlmSectionStruct(&vlmSurface[i].vlmSection[sectionIndex],
                                         &zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
          AIM_STATUS(aimInfo, status);

          // Reset the sectionIndex that is keeping track of the section order.
          zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k].sectionIndex = k;
        }


        if (numVLMControl > 0) {
          AIM_NOTNULL(vlmControl, aimInfo, status);
          // transfer control surface data to sections
          status = get_ControlSurface(bodies,
                                      numVLMControl,
                                      vlmControl,
                                      &zaeroInstance->zaeroProblem.feaAero[surfaceIndex].vlmSurface);
          AIM_STATUS(aimInfo, status);
        }
      }
    }

    if (zaeroInstance->zaeroProblem.feaAero == NULL) {
      status = CAPS_NULLVALUE;
      goto cleanup;
    }
  }

  status = CAPS_SUCCESS;

cleanup:

  if (vlmSurface != NULL) {
    for (i = 0; i < numVLMSurface; i++) {
      status2 = destroy_vlmSurfaceStruct(&vlmSurface[i]);
      if (status2 != CAPS_SUCCESS)
        printf("\tdestroy_vlmSurfaceStruct status = %d\n", status2);
    }
  }

  AIM_FREE(vlmSurface);
  numVLMSurface = 0;

  // Destroy Control
  if (vlmControl != NULL) {
    for (i = 0; i < numVLMControl; i++) {
      (void) destroy_vlmControlStruct(&vlmControl[i]);
    }
    AIM_FREE(vlmControl);
  }
  numVLMControl = 0;

  return status;
}


static int
_getFEMModule(void *aimInfo, capsValue *assign_fem, capsValue *f06Input, zaeroProblemStruct *zaeroProblem) {

  /*! \page zaeroECS ZAero Executive Control Section
   * The following lists the input for the Executive Control Section of the ZAero input file
   */

  int i, j, status;

  int numSuports = 0;
  int printFlag = 0, asetFlag = 0;
  char *FEM = NULL, *form = NULL, *boundary=NULL;
  int *suports=NULL;

  feaSolFileStruct *f06 = NULL;

  if (assign_fem == NULL) return CAPS_NULLVALUE;

  char relPath[PATH_MAX];

  /*! \page zaeroECS
   * \section zaeroASSIGNFEM FEM JSON String Dictionary
   *
   * FEM input must be a JSON string dictionary
   *  (e.g. {"boundary": "sym", "suport": [1,3], "print": 1}
   * where the following keywords ( = default values) may be used:
   */

  for (i = 0; i < assign_fem->length; i++) {

    /*! \page zaeroECS
     *
     *  <ul>
     *  <li> <B>fem = "" </B> </li> <br>
     *   Optional Absolute file path to f06 file from a structural analysis. Cannot be specified
     *   if F06 is linked.
     *  </ul>
     *
     */
    if (strcasecmp("fem", assign_fem->vals.tuple[i].name) == 0) {

      if (f06Input->nullVal == NotNull) {
        AIM_ERROR(aimInfo, "Both \"fem\" is in 'FEM_%d' input and 'F06_%d' cannot be specified", zaeroProblem->numFEMs+1, zaeroProblem->numFEMs+1);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      FEM = string_removeQuotation(assign_fem->vals.tuple[i].value);

      /*! \page zaeroECS
       *  <ul>
       *  <li> <B>form = "" </B> </li> <br>
       *   Optional string describing solver used to generate FEM file. If "fem" is specified then "form" must be one of: <br>
       * 'MSC'     generated by MSC.NASTRAN or NX.NASTRAN
       * 'NE'      generated by NE/NASTRAN
       * 'ASTROS'  generated by ASTROS
       * 'IDEAS'   generated by I-DEAS
       * 'ELFINI'  generated by ELFINI
       * 'GENESIS' generated by GENESIS
       * 'ABAQUS'  generated by ABAQUS
       * 'ALTAIR'  generated by ALTAIR's RADIOSS
       * 'FREE '   stored according to the input instruction described in Remark 9 of ZAero manual
       *  </ul>
       */

      for (j = 0; j < assign_fem->length; j++) {
        if (strcasecmp("form", assign_fem->vals.tuple[j].name) == 0) {
          form = string_removeQuotation(assign_fem->vals.tuple[i].value);
          break;
        }
      }

      if (form == NULL) {
        AIM_ERROR(aimInfo, "Both \"fem\" and \"form\" must be specified", zaeroProblem->numFEMs+1, zaeroProblem->numFEMs+1);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    /*! \page zaeroECS
     *
     *  <ul>
     *  <li> <B>boundary = "SYM" </B> </li> <br>
     *   "boundary" indicates the boundary condition of the structural finite element model.<br>
     *   'SYM' for symmetric boundary condition <br>
     *   'ANTI' for anti-symmetric boundary condition <br>
     *   'ASYM' for asymmetric boundary condition <br>
     *  </ul>
     *
     */
    else if (strcasecmp("boundary", assign_fem->vals.tuple[i].name) == 0) {
      boundary = string_removeQuotation(assign_fem->vals.tuple[i].value);
    }

    /*! \page zaeroECS
     *
     *  <ul>
     *  <li> <B>support = [0, 0] </B> </li> <br>
     *   "support" up to length 2 array for the m/L ZAero support input. See ZAero manual for details.
     *  </ul>
     *
     */
    else if (strcasecmp("suport", assign_fem->vals.tuple[i].name) == 0) {
      status = string_toIntegerDynamicArray(assign_fem->vals.tuple[i].value, &numSuports, &suports);
      AIM_STATUS(aimInfo, status);
      if (numSuports > 2) {
        AIM_ERROR(aimInfo, "Only 2 suport numbers may be specified: numSuport = %d", numSuports);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    /*! \page zaeroECS
     *
     *  <ul>
     *  <li> <B>print = 0 </B> </li> <br>
     *    Print options to the standard output file; where n is an integer. <br>
     *    n = 0 \t no printout of the imported structural free vibration solution. <br>
     *    | n | >= 1 \t print out the structural grid point locations in the aerodynamic coordinate system. <br>
     *    n >= 2 \t print out the modal data (mode shapes) at the structural grid points in the aerodynamic coordinate system. <br>
     *    n <= -2 \t print out the interpolated modal data at the control points of the aerodynamic boxes in the aerodynamic coordinate system. <br>
     *    n = 3 \t print all of the above
     *  </ul>
     *
     */
    else if (strcasecmp("suport", assign_fem->vals.tuple[i].name) == 0) {
      status = string_toInteger(assign_fem->vals.tuple[i].value, &printFlag);
      AIM_STATUS(aimInfo, status);
    }

    /*! \page zaeroECS
     *
     *  <ul>
     *  <li> <B>print = 0 </B> </li> <br>
     *    Print options to the standard output file; where n is an integer. <br>
     *    n = 0 \t no printout of the imported structural free vibration solution. <br>
     *    | n | >= 1 \t print out the structural grid point locations in the aerodynamic coordinate system. <br>
     *    n >= 2 \t print out the modal data (mode shapes) at the structural grid points in the aerodynamic coordinate system. <br>
     *    n <= -2 \t print out the interpolated modal data at the control points of the aerodynamic boxes in the aerodynamic coordinate system. <br>
     *    n = 3 \t print all of the above
     *  </ul>
     *
     */
    else if (strcasecmp("aset", assign_fem->vals.tuple[i].name) == 0) {
      status = string_toInteger(assign_fem->vals.tuple[i].value, &asetFlag);
      AIM_STATUS(aimInfo, status);
    }
  }


  if (f06Input->nullVal == NotNull && FEM == NULL) {

    f06 = (feaSolFileStruct*) f06Input->vals.AIMptr;

    status = aim_relPath(aimInfo, f06->filename, ".", relPath);
    AIM_STATUS(aimInfo, status);

    AIM_STRDUP(FEM, relPath, aimInfo, status);

    switch (f06->fileForm) {
    case MSC_NASTRAN:
      AIM_STRDUP(form, "MSC", aimInfo, status);
      break;
    case NE_NASTRAN:
      AIM_STRDUP(form, "NE", aimInfo, status);
      break;
    case ASTROS:
      AIM_STRDUP(form, "ASTROS", aimInfo, status);
      break;
    case ABAQUS:
      AIM_STRDUP(form, "ABAQUS", aimInfo, status);
      break;
    default:
      AIM_ERROR(aimInfo, "Developer Error: Unknown fea solution file format %d", f06->fileForm);
      status = CAPS_NOTIMPLEMENT;
      goto cleanup;
    }
  }

  if (FEM == NULL) {
    AIM_ERROR(aimInfo, "\"fem\" missing 'FEM_%d' input and 'F06_%d' is not linked", zaeroProblem->numFEMs+1, zaeroProblem->numFEMs+1);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (boundary == NULL) {
    AIM_STRDUP(boundary, "SYM", aimInfo, status);
  }

  i = zaeroProblem->numFEMs;
  zaeroProblem->FEMs[i].filename    = FEM;      FEM = NULL;
  zaeroProblem->FEMs[i].form        = form;     form = NULL;
  zaeroProblem->FEMs[i].boundary    = boundary; boundary = NULL;
  for (j = 0; j < numSuports; j++) {
    AIM_NOTNULL(suports, aimInfo, status);
    zaeroProblem->FEMs[i].suport[j] = suports[j];
  }
  zaeroProblem->FEMs[i].printFlag   = printFlag;
  zaeroProblem->FEMs[i].asetFlag    = asetFlag;
  zaeroProblem->numFEMs++;

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(FEM);
  AIM_FREE(form);
  AIM_FREE(boundary);
  AIM_FREE(suports);

  return status;
}


/* Get HFG module information from input */
static int
_getHFGModule(void *aimInfo,
              capsValue *hfgInput,
              zaeroProblemStruct *zaeroProblem)
{

  int i, status = CAPS_SUCCESS; // Function return

  if (hfgInput == NULL) return CAPS_NULLVALUE;

  for (i = 0; i < hfgInput->length; i++) {

    if (strcasecmp("XZSymmetric", hfgInput->vals.tuple[i].name) == 0) {
      zaeroProblem->hfg.XZSymmetry = string_removeQuotation(hfgInput->vals.tuple[i].value);
      AIM_NOTNULL(zaeroProblem->hfg.XZSymmetry, aimInfo, status);
      string_toUpperCase(zaeroProblem->hfg.XZSymmetry);
    }

    else if (strcasecmp("flip", hfgInput->vals.tuple[i].name) == 0) {
      zaeroProblem->hfg.flip = string_removeQuotation(hfgInput->vals.tuple[i].value);
      AIM_NOTNULL(zaeroProblem->hfg.flip, aimInfo, status);
      string_toUpperCase(zaeroProblem->hfg.flip);
    }
  }

  if (zaeroProblem->hfg.XZSymmetry == NULL)
    AIM_STRDUP(zaeroProblem->hfg.XZSymmetry, "YES", aimInfo, status);

  if (zaeroProblem->hfg.flip == NULL)
    AIM_STRDUP(zaeroProblem->hfg.flip, "NO", aimInfo, status);

  status = CAPS_SUCCESS;
cleanup:
  return status;
}


/* Get UAIC Module configurations from inputs */
static int
_getUAICModuleConfigs(void *aimInfo,
                      int numUAICTuple,
                      capsTuple UAICTuple[],
                      zaeroProblemStruct *zaeroProblem) {

  /*! \page zaeroUAIC ZAero Unified Aerodynamic Influence Coefficients
   * Structure for the UAIC tuple  = ("UAIC Name", "Value").
   * "UAIC Name" defines the reference name for the UAIC being specified.
   * The "Value" must be a JSON String dictionary.
   */

  int i, status;

  zaeroUAICStruct *uaic;

  // Ensure we are starting with no subcases
  if (zaeroProblem->UAICs != NULL) {
    for (i = 0; i < zaeroProblem->numUAICs; i++) {
      status = destroy_zaeroUAICStruct(&zaeroProblem->UAICs[i]);
      AIM_STATUS(aimInfo, status);
    }
    AIM_FREE(zaeroProblem->UAICs);
  }
  zaeroProblem->numUAICs = 0;

  printf("\nGetting ZAERO UAIC Configurations.......\n");

  zaeroProblem->numUAICs = numUAICTuple;
  printf("\tNumber of UAIC Configurations - %d\n", zaeroProblem->numUAICs);

  if (zaeroProblem->numUAICs > 0) {

    AIM_ALLOC(zaeroProblem->UAICs, zaeroProblem->numUAICs, zaeroUAICStruct, aimInfo, status);

    // initiate UAIC structure
    for (i = 0; i < zaeroProblem->numUAICs; i++) {
      status = initiate_zaeroUAICStruct(&zaeroProblem->UAICs[i]);
      AIM_STATUS(aimInfo,status);
    }

  } else {
    printf("\tNumber of analysis UAICs in Analysis tuple is %d\n",
           zaeroProblem->numUAICs);
    return CAPS_NOTFOUND;
  }

  // for each analysis UAIC tuple
  for (i = 0; i < zaeroProblem->numUAICs; i++) {

    uaic = &zaeroProblem->UAICs[i];

    // set name
    AIM_STRDUP(uaic->name, UAICTuple[i].name, aimInfo, status);

    // make sure UAIC tuple value is json string
    if (!json_isDict(UAICTuple[i].value)) {
      AIM_ERROR(aimInfo, "'UAIC' input must be a JSON dictionary");
      return CAPS_BADVALUE;
    }

    // set UAIC ID
    uaic->id = i+1;

    // create the filename for the UAIC
    snprintf(uaic->aicFilename, sizeof(uaic->aicFilename), "UAIC%d", uaic->id);

    /*! \page zaeroUAIC
     * \section jsonUAIC UAIC JSON String Dictionary
     *
     * For the JSON string "Value" dictionary
     *  (e.g. "Value" = {"machNumber": 0.5, "method": 120000.0, "poissonRatio": 0.5, "materialType": "isotropic"})
     *
     *  the following keywords ( = default values) may be used:
     *
     * <ul>
     * <li> <B>machNumber</B> </li> <br>
     *    Mach number.
     * </ul>
     */

    status = json_getDouble(UAICTuple[i].value, "machNumber", &uaic->machNumber);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"machNumber\" in 'UAIC' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    /*! \page zaeroUAIC
     * <ul>
     * <li> <B>method</B> </li> <br>
     *   Integer aerodynamic method <br>
     *   method = 0 for the ZONA6/ZONA7/ZSAP method <br>
     *   method = 1 for the ZTAIC method <br>
     *   method = +/- 2 for the ZONA7U method <br>
     *   method = 3 for the ZTRAN method <br>
     * </ul>
     */

    status = json_getInteger(UAICTuple[i].value, "method", &uaic->methodFlag);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"method\" in 'UAIC' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    /*! \page zaeroUAIC
     * <ul>
     * <li> <B>reducedFreq</B> </li> <br>
     *   List of real reduced frequencies
     * </ul>
     */
    status = json_getDoubleDynamicArray(UAICTuple[i].value, "reducedFreq",
                                        &uaic->numReducedFreq, &uaic->reducedFreq);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"reducedFreq\" in 'UAIC' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    /*! \page zaeroUAIC
     * <ul>
     * <li> <B>print = 0</B> </li> <br>
     *   Integer controlling print output
     * </ul>
     */
    status = json_getInteger(UAICTuple[i].value, "print", &uaic->printFlag);
    AIM_NOTFOUND(aimInfo,status);

    printf("\n\tUAIC Configuration: %s, id = %d\n", uaic->name, uaic->id);
    printf("\t- Mach Number : %f\n", uaic->machNumber);
    printf("\t- Reduced Freq: [");
    if (uaic->numReducedFreq > 0)
      printf("%f", uaic->reducedFreq[0]);
    if (uaic->numReducedFreq > 1)
      printf(", ... , %f", uaic->reducedFreq[uaic->numReducedFreq-1]);
    printf("]\n");

  }

  status = CAPS_SUCCESS;

cleanup:
  return status;
}

/* Get SPLINE module information from inputs */
static int
_getSplineModule(void *aimInfo, capsValue *splineInput, zaeroProblemStruct *zaeroProblem) {

  /*! \page zaeroSpline ZAero Spline Module
   * Structure for the spline tuple is a tuple of key value pairs.
   *
   * The following keywords ( = default values) may be used:
   */

  int i, status = CAPS_SUCCESS; // Function return

  if (splineInput == NULL) return CAPS_NULLVALUE;
  if (splineInput->nullVal == IsNull) return CAPS_SUCCESS;

  // Set default values
  zaeroProblem->spline.method = 1;
  zaeroProblem->spline.attachFlex = 0.0;
  zaeroProblem->spline.eps = 1.0e-6;

  for (i = 0; i < splineInput->length; i++) {

    /*! \page zaeroSpline
     * <ul>
     * <li> <B>method = 1</B> </li> <br>
     *   0 : Imposes zero-displacement condition on aerodynamic boxes.
     *   1 : Defines a surface spline method (Infinite Plate Spline method) for CAERO7.
     *   2 : Defines a beam spline method for CAERO7 / BODY7.
     *   3 : Defines a 3-D spline (Thin Plate Spline method) for CAERO7 / BODY7.
     * </ul>
     */

    if (strcasecmp("method", splineInput->vals.tuple[i].name)) {
      status = string_toInteger(splineInput->vals.tuple[i].value, &zaeroProblem->spline.method);
      AIM_STATUS(aimInfo, status);
    }

    /*! \page zaeroSpline
     * <ul>
     * <li> <B>attachFlex = 0 </B> </li> <br>
     *   Linear attachment flexibility (for spline method 1)
     * </ul>
     */

    else if (strcasecmp("attachFlex", splineInput->vals.tuple[i].name)) {
      status = string_toDouble(splineInput->vals.tuple[i].value, &zaeroProblem->spline.attachFlex);
      AIM_STATUS(aimInfo, status);
    }

    /*! \page zaeroSpline
     * <ul>
     * <li> <B>eps = 1e-6 </B> </li> <br>
     *  Multiplication factor to obtain a small tolerance to detect any duplicated location of
     *  structural grid points. The tolerance is computed by EPS×REFC, where REFC is the
     *  reference chord defined in the AEROZ bulk data card (for spline method 1 or 3)
     * </ul>
     */

    else if (strcasecmp("eps", splineInput->vals.tuple[i].name)) {
      status = string_toDouble(splineInput->vals.tuple[i].value, &zaeroProblem->spline.eps);
      AIM_STATUS(aimInfo, status);
    }
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

static int
_getDisciplineType(void *aimInfo,
                   const char *disciplineStr,
                   zaeroDisciplineTypeEnum *disciplineType) {

  if (
    (strcasecmp(disciplineStr, "LinearFlutter") == 0) ||
    (strcasecmp(disciplineStr, "Flutter") == 0)) {
    *disciplineType = LinearFlutter;
  } else if (
    (strcasecmp(disciplineStr, "ParamFlutter") == 0) ||
    (strcasecmp(disciplineStr, "fltpram") == 0)) {
    *disciplineType = ParamFlutter;
  } else if (
    (strcasecmp(disciplineStr, "Aeroservoelastic") == 0) ||
    (strcasecmp(disciplineStr, "ase") == 0)) {
    *disciplineType = Aeroservoelastic;
  } else if (
    (strcasecmp(disciplineStr, "StaticAeroelastic") == 0) ||
    (strcasecmp(disciplineStr, "Trim") == 0)) {
    *disciplineType = StaticAeroelastic;
  } else if (
    (strcasecmp(disciplineStr, "EjectionLoads") == 0) ||
    (strcasecmp(disciplineStr, "eloads") == 0)) {
    *disciplineType = EjectionLoads;
  } else if (
    (strcasecmp(disciplineStr, "ManeuverLoads") == 0) ||
    (strcasecmp(disciplineStr, "mloads") == 0)) {
    *disciplineType = ManeuverLoads;
  } else if (
    (strcasecmp(disciplineStr, "GustLoads") == 0) ||
    (strcasecmp(disciplineStr, "gloads") == 0)) {
    *disciplineType = GustLoads;
  } else if (
    (strcasecmp(disciplineStr, "MFTGustLoads") == 0) ||
    (strcasecmp(disciplineStr, "mftgust") == 0)) {
    *disciplineType = MFTGustLoads;
  } else if (
    (strcasecmp(disciplineStr, "NonLinearFlutter") == 0) ||
    (strcasecmp(disciplineStr, "nlfltr") == 0)) {
    *disciplineType = NonLinearFlutter;
  } else {
    *disciplineType = UnknownDiscipline;
    AIM_ERROR(aimInfo, "Unknown discipline: %s", disciplineStr);
    return CAPS_BADVALUE;
  }

  return CAPS_SUCCESS;
}

/* Get analysis subcase information from inputs */
static int
_getAnalysisSubcases(void *aimInfo,
                     const cfdUnitsStruct* units,
                     int numAnalysisTuple,
                     capsTuple analysisTuple[],
                     zaeroProblemStruct *zaeroProblem) {

  /*! \page zaeroAnalysis ZAero Analysis
   * Structure for the Analysis tuple  = ("Case Name", "Value").
   * "Case Name" defines the reference name for the subcase being specified.
   * The "Value" must be a JSON String dictionary.
   */

  int i, status;
  int maxchar = 72; // case control strings can be at most 72 chars
  char *disciplineStr = NULL;

  zaeroSubcaseStruct *subcase;

  // Ensure we are starting with no subcases
  if (zaeroProblem->subcases != NULL) {
    for (i = 0; i < zaeroProblem->numSubcases; i++) {
      status = destroy_zaeroSubcaseStruct(&zaeroProblem->subcases[i]);
      AIM_STATUS(aimInfo, status);
    }
  }
  AIM_FREE(zaeroProblem->subcases);
  zaeroProblem->numSubcases = 0;

  printf("\nGetting ZAERO Analysis Subcases.......\n");

  zaeroProblem->numSubcases = numAnalysisTuple;
  printf("\tNumber of Analysis Subcases - %d\n", zaeroProblem->numSubcases);

  if (zaeroProblem->numSubcases > 0) {
    AIM_ALLOC(zaeroProblem->subcases, zaeroProblem->numSubcases, zaeroSubcaseStruct, aimInfo, status);

    // initiate subcase structure
    for (i = 0; i < zaeroProblem->numSubcases; i++) {
      status = initiate_zaeroSubcaseStruct(&zaeroProblem->subcases[i]);
      AIM_STATUS(aimInfo, status);
    }
  } else {
    AIM_ERROR(aimInfo, "Number of analysis subcases in 'Analysis' tuple is %d", zaeroProblem->numSubcases);
    return CAPS_NOTFOUND;
  }

  // for each analysis subcase tuple
  for (i = 0; i < zaeroProblem->numSubcases; i++) {

    subcase = &zaeroProblem->subcases[i];

    // set name
    AIM_STRDUP(subcase->name, analysisTuple[i].name, aimInfo, status);

    // make sure analysis tuple value is json string
    if (!json_isDict(analysisTuple[i].value)) {
      AIM_ERROR(aimInfo, "'Analysis' tuple value must be a JSON dictionary");
      return CAPS_BADVALUE;
    }

    // set subcase ID
    subcase->subcaseID = i+1;

    // set subcase analysis ID (at the moment
    // TODO: is there a better way to determine analysisID?
    subcase->analysisID = subcase->subcaseID * 100;

    /*! \page zaeroAnalysis
     * \section jsonSubcase Analysis JSON String Dictionary
     *
     * For the JSON string "Value" dictionary
     *  (e.g. "Value" = {"discipline": "LinearFlutter", "uaic": "cruise"})
     * \endif
     *  the following keywords ( = default values) may be used:
     *
     * <ul>
     * <li> <B>discipline = ""</B> </li> <br>
     *    Analysis discipline JSON string. Options : <br>
     *    LinearFlutter or flutter - Linear flutter<br>
     *    ParamFlutter or fltpram - Parametric flutter.<br>
     *    Aeroservoelastic or ase - Asymmetric parametric flutter.<br>
     *    StaticAeroelastic or trim - Static areoelastic (see \ref jsonTRIM).<br>
     *    EjectionLoads or eloads - Transient ejection loads.<br>
     *    ManeuverLoads or mloads - Transient manouver loads.<br>
     *    GustLoads or gloads - Discrete gust load.<br>
     *    MFTGustLoads or mftgust - Continuous gust load.<br>
     *    NonLinearFlutter or nlfltr - Nonlinear flutter.<br>
     * </ul>
     */

    // setup discipline according to "discipline" value
    status = json_getString( analysisTuple[i].value, "discipline", &disciplineStr);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "input for %s analysis is missing required entry"
                         "\"discipline\"", subcase->name);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    AIM_NOTNULL(disciplineStr, aimInfo, status);

    status = _getDisciplineType(aimInfo, disciplineStr, &subcase->disciplineType);
    AIM_STATUS(aimInfo, status);

    switch (subcase->disciplineType) {
      case LinearFlutter:
        status = zaero_getLinearFlutterDiscipline(aimInfo, analysisTuple[i].value, subcase);
        break;
      case ParamFlutter:
        status = CAPS_NOTIMPLEMENT;
        break;
      case Aeroservoelastic:
        status = CAPS_NOTIMPLEMENT;
        break;
      case StaticAeroelastic:
        status = zaero_getTrimDiscipline(aimInfo, analysisTuple[i].value, units, subcase);
        break;
      case ManeuverLoads:
        status = CAPS_NOTIMPLEMENT;
        break;
      case EjectionLoads:
        status = CAPS_NOTIMPLEMENT;
        break;
      case GustLoads:
        status = CAPS_NOTIMPLEMENT;
        break;
      case MFTGustLoads:
        status = CAPS_NOTIMPLEMENT;
        break;
      case NonLinearFlutter:
        status = CAPS_NOTIMPLEMENT;
        break;
      // should not be unknown at this point, including to avoid switch warning
      case UnknownDiscipline:
        status = CAPS_BADVALUE;
        break;
    }
    AIM_STATUS(aimInfo, status);

    /*! \page zaeroAnalysis
     * <ul>
     * <li> <B>uaic = ""</B> </li> <br>
     *    Name of the UAIC module for the subcase
     * </ul>
     */

    status = json_getString(analysisTuple[i].value, "uaic", &subcase->uaicName);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "Subcase '%s' is missing required \"uaic\"", subcase->name);
      return status;
    }

    AIM_ALLOC(subcase->subtitle, maxchar + 1, char, aimInfo, status);
    snprintf(subcase->subtitle, maxchar, "%s Analysis, using %s",
             disciplineStr, subcase->uaicName);

    status = json_getString(analysisTuple[i].value, "label", &subcase->label);
    AIM_NOTFOUND(aimInfo, status);

    printf("\n\tAnalysis Subcase: %s, id = %d\n",
            subcase->name, subcase->subcaseID);
    printf("\t- Discipline : %s\n", disciplineStr);
    printf("\t- UAIC       : %s\n", subcase->uaicName);
  }

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(disciplineStr);

  return status;
}

static void
_writeBanner(FILE *fp, char *msg) {
  fprintf(fp,
"$==============================================================================\n"
"$ %s\n"
"$==============================================================================\n",
    msg);
}

static void
_writeCaseBanner(FILE *fp, char *msg) {
  fprintf(fp,
"$------------------------------------------------------------------------------\n"
"$ %s\n"
"$------------------------------------------------------------------------------\n",
    msg);
}

/* Populate zaero problem struct with input information*/
static int
_getZaeroProblemData(void *aimInfo,
                     capsValue *aimInputs,
                     aimStorage *zaeroInstance) {

  int status = CAPS_SUCCESS;

  capsValue *hfgInput, *uaicInput, *splineInput;
  capsValue *analysisInput, *trimVarInput;

  const cfdUnitsStruct *units = &zaeroInstance->units;
  zaeroProblemStruct *zaeroProblem = &zaeroInstance->zaeroProblem;

  hfgInput      = &aimInputs[inHFG-1];
  uaicInput     = &aimInputs[inUAIC-1];
  splineInput   = &aimInputs[inSpline-1];
  analysisInput = &aimInputs[inAnalysis-1];
  trimVarInput  = &aimInputs[inTrim_Variable-1];

  if (aimInputs[inFEM_1-1].nullVal == IsNull) {
    AIM_ERROR(aimInfo, "'FEM_1' input is required");
    return CAPS_BADVALUE;
  }

  status = _getFEMModule(aimInfo, &aimInputs[inFEM_1-1], &aimInputs[inF06_1-1], zaeroProblem);
  AIM_STATUS(aimInfo, status);

  if (aimInputs[inFEM_2-1].nullVal == NotNull) {
    status = _getFEMModule(aimInfo, &aimInputs[inFEM_2-1], &aimInputs[inF06_2-1], zaeroProblem);
    AIM_STATUS(aimInfo, status);
  }

  status = _getHFGModule(aimInfo, hfgInput, zaeroProblem);
  AIM_STATUS(aimInfo, status);

  if (uaicInput->nullVal == NotNull) {
    status = _getUAICModuleConfigs(aimInfo,
                                   uaicInput->length,
                                   uaicInput->vals.tuple,
                                   zaeroProblem);
    AIM_STATUS(aimInfo, status);
  } else {
    printf("UAIC tuple is NULL - No UAIC configurations set");
  }

  status = _getSplineModule(aimInfo, splineInput, zaeroProblem);
  AIM_STATUS(aimInfo, status);

  if (analysisInput->nullVal == NotNull) {
    status = _getAnalysisSubcases(aimInfo,
                                  units,
                                  analysisInput->length,
                                  analysisInput->vals.tuple,
                                  zaeroProblem);
    AIM_STATUS(aimInfo, status);
  } else {
    printf("Analysis tuple is NULL - No analysis subcase set");
  }

  if (trimVarInput->nullVal == NotNull) {
    status = zaero_getTrimVariables(aimInfo,
                                    trimVarInput->length,
                                    trimVarInput->vals.tuple,
                                    zaeroProblem);
    AIM_STATUS(aimInfo, status);
  } else {
    printf("Trim_Variable tuple is NULL - No trim variables set\n");
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

static int
_writeAssignFEM(void *aimInfo, FILE *fp, const zaeroFEMStruct *fem) {

  int status = CAPS_SUCCESS;
  char *tempString = NULL;

  // FEM filename
  fprintf(fp, "ASSIGN FEM = %s,\n", fem->filename);

  // form, default is MSC
  if (fem->form == NULL) {
    AIM_STRDUP(tempString , "MSC", aimInfo, status);
  }
  else {
    AIM_STRDUP(tempString, fem->form, aimInfo, status);
    string_toUpperCase(tempString);
  }
  fprintf(fp, "FORM = %s, ", tempString);

  AIM_FREE(tempString);

  // boundary condition, default is SYM
  if (fem->boundary == NULL) {
    AIM_STRDUP(tempString, "SYM", aimInfo, status);
  }
  else {
    AIM_STRDUP(tempString, fem->boundary, aimInfo, status);
    string_toUpperCase(tempString);
  }
  fprintf(fp, "BOUNDARY = %s, ", tempString);

  AIM_FREE(tempString);

  // print
  fprintf(fp, "PRINT = %d, ", fem->printFlag);

  // suport
  if (fem->suport[0] != 0 && fem->suport[1] != 0) {
    fprintf(fp, "SUPORT = %d/%d, ", fem->suport[0], fem->suport[1]);
  }
  else if (fem->suport[0] != 0) {
    fprintf(fp, "SUPORT = %d, ", fem->suport[0]);
  }

  if (fem->asetFlag > 0) {
    fprintf(fp, "ASET = YES");
  }
  else {
    fprintf(fp, "ASET = NO");
  }

  fprintf(fp, "\n");

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(tempString);
  return status;
}

/* Executive Control Section
 * From the manual:
 * """
 * The Executive Control Section must be located at the beginning of the input file.
 * Its major functions are:
 *   *  to define the filename that contains the free vibration output from the
 *      structural finite element methods"
 *   *  to allow direct matrix input
 *   *  to turn on diagnostic routines
 * """
 */
static int
_writeExecutiveControlSection(void *aimInfo,
                              FILE *fp,
                              capsValue *aimInputs,
                              const aimStorage *zaeroInstance) {

  char *tempString = NULL;

  int i, tempInteger;

  int status; // Function return

  _writeBanner(fp, "Executive Control Section");

  // ASSIGN FEM
  for (i = 0; i < zaeroInstance->zaeroProblem.numFEMs; i++) {
    status = _writeAssignFEM(aimInfo, fp, &zaeroInstance->zaeroProblem.FEMs[i]);
    AIM_STATUS(aimInfo, status);
  }

  // CPU
  tempInteger = aimInputs[inCPU-1].vals.integer;
  fprintf(fp, "CPU %d\n", tempInteger);

  // MEMORY
  tempString = aimInputs[inMemory-1].vals.string;
  fprintf(fp, "MEMORY %s\n", tempString);

  // DOUBLE
  fprintf(fp, "DOUBLE\n");

  // CEND
  fprintf(fp, "CEND\n");

  status = CAPS_SUCCESS;

cleanup:

  return status;
}

// write analysis subcase in Case Control Section
static int
_writeCase(void *aimInfo, FILE *fp, zaeroSubcaseStruct *subcase) {

  fprintf(fp, "SUBCASE = %d\n", subcase->subcaseID);

  fprintf(fp, "SUBTITLE = %s\n", subcase->subtitle);

  if (subcase->label != NULL)
    fprintf(fp, "LABEL = %s\n", subcase->label);

  if (subcase->disciplineType == LinearFlutter) {
    fprintf(fp, "FLUTTER = %d\n", subcase->analysisID);
  }
  else if (subcase->disciplineType == ParamFlutter) {
    AIM_ERROR(aimInfo, "ParamFlutter discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == Aeroservoelastic) {
    AIM_ERROR(aimInfo, "Aeroservoelastic discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == StaticAeroelastic) {
    fprintf(fp, "TRIM = %d\n", subcase->analysisID);
  }
  else if (subcase->disciplineType == ManeuverLoads) {
    AIM_ERROR(aimInfo, "ManeuverLoads discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == EjectionLoads) {
    AIM_ERROR(aimInfo, "EjectionLoads discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == GustLoads) {
    AIM_ERROR(aimInfo, "GustLoads discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == MFTGustLoads) {
    AIM_ERROR(aimInfo, "MFTGustLoads discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else if (subcase->disciplineType == NonLinearFlutter) {
    AIM_ERROR(aimInfo, "NonLinearFlutter discipline not implemented");
    return CAPS_NOTIMPLEMENT;
  }
  else {
    AIM_ERROR(aimInfo, "Unknown discipline");
    return CAPS_BADVALUE;
  }

  return CAPS_SUCCESS;
}

/* Case Control Section
 * From the manual:
 * """
 * The Case Control Section must be located after the Executive Control Section and
 * before the Bulk Data Section. Its major functions are:
 *   *  to input title cards that describe the ZAERO analysis
 *   *  to select the disciplines (FLUTTER, FLTPRAM, ASE, TRIM, MLOADS, ELOADS, GLOADS
 *      or NLFLTR) for the analysis
 * """
 */
static int
_writeCaseControlSection(void *aimInfo,
                         FILE *fp,
                         capsValue *aimInputs,
                         const aimStorage *zaeroInstance) {
  int status = CAPS_SUCCESS;

  int i;

  char *tempString = NULL;

  _writeBanner(fp, "Case Control Section");

  // TITLE
  tempString = aimInputs[inProj_Name-1].vals.string;
  fprintf(fp, "TITLE = %s\n", tempString);

  // ECHO
  tempString = aimInputs[inEcho-1].vals.string;
  fprintf(fp, "ECHO = %s\n", tempString);

  // write each analysis subcase
  for (i = 0; i < zaeroInstance->zaeroProblem.numSubcases; i++) {
    status = _writeCase(aimInfo, fp, &zaeroInstance->zaeroProblem.subcases[i]);
    AIM_STATUS(aimInfo, status);
  }

cleanup:
  return status;
}

/* Flight Condition Definition */
static int
_writeFlightConditionCards(void *aimInfo,
                           FILE *fp,
              /*@unused@*/ capsValue *aimInputs,
                           const aimStorage *zaeroInstance) {

  int i, status;

  _writeCaseBanner(fp, "Flight Condition Definition");
  // printf("FLIGHT CONDITION DEFINITION\n");

  status = zaero_modelPhysicalData( aimInfo,
                                    fp,
                                    &zaeroInstance->zaeroProblem.hfg,
                                    &zaeroInstance->units,
                                    zaeroInstance->feaFormatType);
  AIM_STATUS(aimInfo, status);

  // write cards for each UAIC configuration
  for (i = 0; i < zaeroInstance->zaeroProblem.numUAICs; i++) {

    status = zaero_unsteadyAerodynamicsGenerator(aimInfo,
                                                 fp,
                                                 &zaeroInstance->zaeroProblem.UAICs[i],
                                                 zaeroInstance->feaFormatType);
    AIM_STATUS(aimInfo, status);

  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

/* Aerodynamic Model Definition */
static int
_writeAerodynamicModelCards(void *aimInfo,
                            FILE *fp,
               /*@unused@*/ capsValue *aimInputs,
                            const aimStorage *zaeroInstance) {
  int status, i;

  _writeCaseBanner(fp, "Aerodynamic Model Definition");

  for (i = 0; i < zaeroInstance->zaeroProblem.numAero; i++) {
    status = zaero_aerodynamicWingComponent(aimInfo,
                                            fp,
                                            &zaeroInstance->zaeroProblem.feaAero[i],
                                            i,
                                            zaeroInstance->feaFormatType);
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

/* Spline */
static int
_writeSplineCards(void *aimInfo,
                  FILE *fp,
     /*@unused@*/ capsValue *aimInputs,
                  const aimStorage *zaeroInstance) {
  int status, i;

  _writeCaseBanner(fp, "Spline Definition");

  for (i = 0; i < zaeroInstance->zaeroProblem.numAero; i++) {
    status = zaero_splineMethod(aimInfo,
                                fp,
                                &zaeroInstance->zaeroProblem.feaAero[i],
                                &zaeroInstance->zaeroProblem.spline,
                                zaeroInstance->feaFormatType);
    AIM_STATUS(aimInfo, status);
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

/* Analysis Definition */
static int
_writeAnalysisCards(void *aimInfo,
                    FILE *fp,
       /*@unused@*/ capsValue *aimInputs,
                    const aimStorage *zaeroInstance) {

  int i, status;

  zaeroSubcaseStruct *subcase = NULL;

  // write analysis cards for each subcase
  for (i = 0; i < zaeroInstance->zaeroProblem.numSubcases; i++) {

    subcase = &zaeroInstance->zaeroProblem.subcases[i];

    if (subcase->disciplineType == LinearFlutter) {

      _writeCaseBanner(fp, "Linear Flutter Analysis");

      status = zaero_flutterAnalysis(aimInfo,
                                     fp,
                                     subcase,
                                     &zaeroInstance->zaeroProblem,
                                     &zaeroInstance->units,
                                     zaeroInstance->feaFormatType);
      AIM_STATUS(aimInfo, status);
    }
    else if (subcase->disciplineType == StaticAeroelastic) {

      _writeCaseBanner(fp, "Static Aeroelastic/Trim Analysis");

      status = zaero_trimAnalysis(aimInfo,
                                  fp,
                                  subcase,
                                  &zaeroInstance->zaeroProblem,
                                  zaeroInstance->feaFormatType);
      AIM_STATUS(aimInfo, status);
    }
    else {
      AIM_ERROR(aimInfo, "Unknown subcase disciplineType.");
      return CAPS_BADVALUE;
    }
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

/* Output */
static int
_writeOutputCards(void *aimInfo,
                  FILE *fp,
                  const int numOutputTuple,
                  const capsTuple outputTuple[],
                  const aimStorage *zaeroInstance)
{
  int i, status, setID=9000;

  for (i = 0; i < numOutputTuple; i++) {

    if (strcmp(outputTuple[i].name, "Aero") == 0) {

      status = zaero_textFileGenerationAero(aimInfo,
                                            fp,
                                            outputTuple[i].value,
                                            &zaeroInstance->zaeroProblem,
                                            zaeroInstance->feaFormatType,
                                            &setID);
      AIM_STATUS(aimInfo, status);

    } else if (strcmp(outputTuple[i].name, "Flutter") == 0) {

      status = zaero_textFileGenerationFlutter(aimInfo,
                                               fp,
                                               outputTuple[i].value,
                                               &zaeroInstance->zaeroProblem,
                                               zaeroInstance->feaFormatType,
                                               &setID);
      AIM_STATUS(aimInfo, status);

    } else if (strcmp(outputTuple[i].name, "Mode") == 0) {

      status = zaero_textFileGenerationMode(aimInfo,
                                            fp,
                                            outputTuple[i].value,
                                            &zaeroInstance->zaeroProblem,
                                            zaeroInstance->feaFormatType,
                                            &setID);
      AIM_STATUS(aimInfo, status);

    } else if (strcmp(outputTuple[i].name, "Trim") == 0) {

      status = zaero_textFileGenerationTrim(aimInfo,
                                            fp,
                                            outputTuple[i].value,
                                            &zaeroInstance->zaeroProblem,
                                            zaeroInstance->feaFormatType,
                                            &setID);
      AIM_STATUS(aimInfo, status);

    }  else if (strcmp(outputTuple[i].name, "VG") == 0) {

      status = zaero_textFileGenerationVG(aimInfo,
                                          fp,
                                          outputTuple[i].value,
                                          &zaeroInstance->zaeroProblem,
                                          zaeroInstance->feaFormatType,
                                          &setID);
      AIM_STATUS(aimInfo, status);
    }
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

static int
_writeOutputCardsIfRequested(void *aimInfo,
                             FILE *fp,
                             capsValue *aimInputs,
                             const aimStorage *zaeroInstance) {
  int status;

  capsValue *outputInput;

  outputInput = &aimInputs[inOutput-1];

  // if output generation requested
  if (outputInput->nullVal == NotNull) {

    _writeBanner(fp, "OUTPUT");

    status = _writeOutputCards(aimInfo,
                               fp,
                               outputInput->length,
                               outputInput->vals.tuple,
                               zaeroInstance);
    AIM_STATUS(aimInfo, status);

  } else {
    printf("No graphical outputs requested.\n");
  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

/* Bulk Data Section
 * From the manual:
 * """
 * The Bulk Data Section begins right after the BEGIN BULK Case Control Command and
 * ends at a bulk data card ENDDATA. The Bulk Data Section contains data cards that
 * specify:
 *   *  the geometry of the aerodynamic model
 *   *  spline for displacement and force transversal between the structural finite
 *      element grid points and aerodynamic boxes
 *   *  the Mach numbers, natural frequencies and aerodynamic methods for unsteady
 *      aerodynamic data generation
 *   *  disciplines (FLUTTER, ASE, FLTPRAM, static aerolastic/TRIM, MLOADS, ELOADS,
 *      GLOADS, or NLFLTR) to be analyzed
 *   *  other miscellaneous inputs
 * """
*/
static int
_writeBulkDataSection(void *aimInfo,
                      FILE *fp,
                      capsValue *aimInputs,
                      const aimStorage *zaeroInstance) {

  int status;

  _writeBanner(fp, "$ Begin Bulk Data Section");

  fprintf(fp, "BEGIN BULK\n");

  status = _writeFlightConditionCards(aimInfo, fp, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeAerodynamicModelCards(aimInfo, fp, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeSplineCards(aimInfo, fp, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeAnalysisCards(aimInfo, fp, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeOutputCardsIfRequested(aimInfo, fp, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  fprintf(fp, "ENDDATA\n");

  status = CAPS_SUCCESS;
cleanup:
  return status;
}


static int
_setupArtifacts(void *aimInfo, capsValue *aimInputs,
                aimStorage *zaeroInstance) {

  int i, status = CAPS_SUCCESS;

  char *inpExt = ".inp", *outExt = ".out";
  zaeroUAICStruct *uaic = NULL;
  zaeroArtifactsStruct *artifacts = &zaeroInstance->artifacts;

  // ensure starting with fresh zaeroArtifactsStruct
  if (artifacts != NULL) {
    status = destroy_zaeroArtifactsStruct(artifacts);
    if (status != CAPS_SUCCESS) {
      return status;
    }
  }
  AIM_NOTNULL(artifacts, aimInfo, status);

  // setup input filename
  artifacts->input = string_format("%s%s", aimInputs[inProj_Name-1].vals.string, inpExt, NULL);

  // setup output filename
  artifacts->output = string_format("%s%s", aimInputs[inProj_Name-1].vals.string, outExt, NULL);

  // TODO: setup modal FEM filename(s) ?

  // setup AIC filename for each UAIC configuration

  artifacts->numAIC = zaeroInstance->zaeroProblem.numUAICs;
  AIM_ALLOC(artifacts->aic, artifacts->numAIC, capsTuple, aimInfo, status);
  for (i = 0; i < artifacts->numAIC; i++) {
    artifacts->aic[i].name = NULL;
    artifacts->aic[i].value = NULL;
  }

  for (i = 0; i < artifacts->numAIC; i++) {
    uaic = &zaeroInstance->zaeroProblem.UAICs[i];
    AIM_STRDUP(artifacts->aic[i].name, uaic->name, aimInfo, status);
    AIM_STRDUP(artifacts->aic[i].value, uaic->aicFilename, aimInfo, status);
  }

cleanup:
  return status;
}


static int
_propagateExtractedFEMVariables(void *aimInfo,
                                zaeroProblemStruct *zaeroProblem,
                                double mass,
                                double centerGravity[3],
                                double inertia[6]) {

  int i, j;//, status;

  zaeroSubcaseStruct *subcase;
  zaeroTrimStruct *trim;

  for (i = 0; i < zaeroProblem->numSubcases; i ++) {

    subcase = &zaeroProblem->subcases[i];

    if (subcase->disciplineType == StaticAeroelastic) {

      if (subcase->discipline == NULL) {
        AIM_ERROR(aimInfo, "subcase discipline is NULL");
        return CAPS_NULLVALUE;
      }

      trim = (zaeroTrimStruct *) subcase->discipline;
      printf("\n\tUpdating StaticAeroelastic subcase: %s\n", subcase->name);

      // calculate vectorToCG from aero moment center and center of gravity
      for (j = 0; j < 3; j++) {
        // only set if user did not already set this value
        if (trim->vectorToCG[j] != -1.0) {
          trim->vectorToCG[j] = (
              centerGravity[j] - zaeroProblem->hfg.refCenter[j]);
        }
      }
      printf("\t- vectorToCG -> [%lf, %lf, %lf]\n",
             trim->vectorToCG[0], trim->vectorToCG[1], trim->vectorToCG[2]);

      // set weight, if xzsymmetric then double
      // if (trim->gravityAcceleration == 0.0) {
      //   AIM_ERROR(aimInfo, "gravityAcceleration is 0.0 in subcase: %s", subcase->name);
      // }
      trim->weight = mass;// * trim->gravityAcceleration;
      if (strcasecmp(zaeroProblem->hfg.XZSymmetry, "YES") == 0) {
        trim->weight *= 2;
      }
      printf("\t- weight -> %f\n", trim->weight);

      // set weight moment of inertia, copy extracted inertia matrix
      for (j = 0; j < 6; j++) {
        trim->weightMomentOfInertia[j] = inertia[j];
      }
      printf("\t- weightMomentOfInertia -> [%lf, %lf, %lf, %lf, %lf, %lf]\n",
             trim->weightMomentOfInertia[0], trim->weightMomentOfInertia[1],
             trim->weightMomentOfInertia[2], trim->weightMomentOfInertia[3],
             trim->weightMomentOfInertia[4], trim->weightMomentOfInertia[5]);
    }
  }

  return CAPS_SUCCESS;
}


static int
_getMassProp(void *aimInfo,
                     capsValue *aimInputs,
                     aimStorage *zaeroInstance)
{
  int i, status;

  double mass = 0.0;
  double CG[3];
  double inertia[6];

  double Lunit=1.0;
  const char *Lunits, *Munits;
  char *Iunits = NULL, *tmpUnits = NULL;

  feaMassPropStruct *feaMassProp=NULL;

  const char *bodyLunits = NULL;
  const cfdUnitsStruct *units = &zaeroInstance->units;

  // initialize matrices to zero
  for (i = 0; i < 3; i++)
    CG[i] = 0.0;

  for (i = 0; i < 6; i++)
    inertia[i] = 0.0;

  if (units->length != NULL)
    Lunits = units->length;
  else
    Lunits = "m";

  if (units->mass != NULL)
    Munits = units->mass;
  else
    Munits = "kg";

  if (units->length != NULL)
  {
    // Get length units
    status = aim_capsLength(aimInfo, &bodyLunits);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "No units assigned *** capsLength is not set in *.csm file!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // conversion of the csm model units into units of Lunits
    Lunit = 1.0;
    status = aim_convert(aimInfo, 1, bodyLunits, &Lunit, Lunits, &Lunit);
    AIM_STATUS(aimInfo, status);

    status = aim_unitRaise(aimInfo, Lunits, 2, &tmpUnits ); // length^2
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(tmpUnits, aimInfo, status);
    status = aim_unitMultiply(aimInfo, Munits, tmpUnits, &Iunits ); // mass*length^2, e.g moment of inertia
    AIM_STATUS(aimInfo, status);
    AIM_FREE(tmpUnits);
    AIM_NOTNULL(Iunits, aimInfo, status);
  }

  if (aimInputs[inMassPropLink-1].nullVal == NotNull) {
    feaMassProp = (feaMassPropStruct *) aimInputs[inMassPropLink-1].vals.AIMptr;

    // Mass
    mass = feaMassProp->mass;

    if (units->mass != NULL) {
      status = aim_convert(aimInfo, 1, feaMassProp->massUnit, &mass, Munits, &mass);
      AIM_STATUS(aimInfo, status);
    }

    // Center of gravity
    CG[0] = feaMassProp->CG[0];
    CG[1] = feaMassProp->CG[1];
    CG[2] = feaMassProp->CG[2];

    if (units->length != NULL) {
      status = aim_convert(aimInfo, 3, feaMassProp->lengthUnit, CG, Lunits, CG);
      AIM_STATUS(aimInfo, status);
    }

    // Inertia order = Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    inertia[I11] = feaMassProp->massInertia[I11];
    inertia[I22] = feaMassProp->massInertia[I22];
    inertia[I33] = feaMassProp->massInertia[I33];
    inertia[I21] = feaMassProp->massInertia[I21];
    inertia[I31] = feaMassProp->massInertia[I31];
    inertia[I32] = feaMassProp->massInertia[I32];

    if (units->length != NULL) {
      status = aim_convert(aimInfo, 6, feaMassProp->momentOfInertiaUnit, inertia, Iunits, inertia);
      AIM_STATUS(aimInfo, status);
    }
  }

  status = _propagateExtractedFEMVariables(aimInfo,
                                           &zaeroInstance->zaeroProblem,
                                           mass, CG, inertia);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(Iunits);

  return status;
}


/****************** exposed AIM entry points -- Analysis **********************/

/* aimInitialize: Initialization Information for the AIM */
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int status;
  aimStorage *zaeroInstance = NULL;

  const char *keyWord;
  char *keyValue = NULL;
  double real = 1;
  cfdUnitsStruct *units=NULL;

#ifdef DEBUG
  printf("\n zaeroAIM/aimInitialize   inst = %d!\n", inst);
#endif

  /* specify the number of analysis inputs  defined in aimInputs
   *     and the number of analysis outputs defined in aimOutputs */
  *nIn    = NUMINPUT;
  *nOut   = NUMOUTPUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate */
  *nFields = 0;
  *franks  = NULL;
  *fnames  = NULL;
  *fInOut  = NULL;

  // Allocate zaeroInstance
  AIM_ALLOC(zaeroInstance, 1, aimStorage, aimInfo, status);

  // Initialize instance storage
  initiate_aimStorage(zaeroInstance);

  /*! \page aimUnitsZAERO AIM Units
   *  A unit system may be optionally specified during AIM instance initiation. If
   *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
   *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
   *  unit consistent. A unit system may be specified via a JSON string dictionary for example:
   *  unitSys = "{"mass": "kg", "length": "m", "time":"seconds", "temperature": "K"}"
   */
  if (unitSys != NULL) {
    units = &zaeroInstance->units;

    // Do we have a json string?
    if (strncmp( unitSys, "{", 1) != 0) {
      AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
      return CAPS_BADVALUE;
    }

    /*! \page aimUnitsZAERO
     *  \section jsonStringZAERO JSON String Dictionary
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

    /*! \page aimUnitsZAERO
     *  <ul>
     *  <li> <B>length = "None"</B> </li> <br>
     *  Length units - e.g. "meter", "m", "inch", "in", "mile", ...
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

    /*! \page aimUnitsZAERO
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

    /*! \page aimUnitsZAERO
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

    status = cfd_cfdDerivedUnits(aimInfo, units);
    AIM_STATUS(aimInfo, status);
  }

  *instStore = zaeroInstance;

  status = CAPS_SUCCESS;

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
/* aimInputs: Input Information for the AIM */
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
  /*! \page aimInputsZAERO AIM Inputs
   * The following list outlines the Zaero inputs along with their default value available
   * through the AIM interface.
   */
  int status = CAPS_SUCCESS;
  aimStorage *zaeroInstance;
  cfdUnitsStruct *units=NULL;

  *ainame = NULL;

  zaeroInstance = (aimStorage *) instStore;
  if (zaeroInstance == NULL) AIM_STATUS(aimInfo, CAPS_NULLVALUE);

  if (zaeroInstance != NULL) units = &zaeroInstance->units;

// #ifdef DEBUG
//   printf(" zaeroAIM/aimInputs  zaeroInstances = %d  index  = %d!\n", iIndex, index);
// #endif

  // Zaero Inputs
  if (index == inProj_Name) {
    *ainame              = EG_strdup("Proj_Name");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("zaero_CAPS");
    defval->lfixed       = Change;

    /*! \page aimInputsZAERO
     * - <B> Proj_Name = "zaero_CAPS"</B> <br>
     * This corresponds to the project name used for file naming.
     */

  } else if (index == inAnalysis) {
    *ainame              = EG_strdup("Analysis");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> Analysis = NULL</B> <br>
     * Analysis tuple used to input analysis/case information for the model.
     */

  } else if (index == inFile_Format) {
    *ainame              = EG_strdup("File_Format");
    defval->type         = String;
    defval->vals.string  = EG_strdup("Small"); // Small, Free

    /*! \page aimInputsZAERO/
     * - <B> File_Format = "Small"</B> <br>
     * Formatting type for the bulk file. Options: "Small", "Large", "Free".
     */

  } else if (index == inFEM_1) {
    *ainame              = EG_strdup("FEM_1");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> FEM_1 = NULL </B> <br>
     * JSON dictionary for first ZAero ASSIGN FEM inputs in Executive Control Section (See \ref zaeroECS)
     */

  } else if (index == inF06_1) {
    *ainame              = EG_strdup("F06_1");
    defval->type         = Pointer;
    defval->nullVal      = IsNull;
    AIM_STRDUP(defval->units, "feaSolFileStruct", aimInfo, status);

    /*! \page aimInputsZAERO
     * - <B> F06_1 = NULL </B> <br>
     * Link for F06 file from from a structural analysis AIM for first ZAero ASSIGN FEM.
     * zaeroAIM will attempt to extract / determine as many analysis parameters from the F06 file as possible.
     */

  } else if (index == inFEM_2) {
    *ainame              = EG_strdup("FEM_2");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> FEM_2 = NULL </B> <br>
     * JSON dictionary for first ZAero ASSIGN FEM inputs in Executive Control Section (See \ref zaeroECS)
     */

  } else if (index == inF06_2) {
    *ainame              = EG_strdup("F06_2");
    defval->type         = Pointer;
    defval->nullVal      = IsNull;
    AIM_STRDUP(defval->units, "feaSolFileStruct", aimInfo, status);

    /*! \page aimInputsZAERO
     * - <B> F06_2 = NULL </B> <br>
     * Link for F06 file from from a structural analysis AIM for second ZAero ASSIGN FEM.
     * zaeroAIM will attempt to extract / determine as many analysis parameters from the F06 file as possible.
     */

  } else if (index == inCPU) {
    *ainame              = EG_strdup("CPU");
    defval->type         = Integer;
    defval->vals.integer = 1;

    /*! \page aimInputsZAERO
     * - <B> CPU = 1 </B> <br>
     * Defines the number of processors for parallel computation.
     */

  } else if (index == inMemory) {
    *ainame              = EG_strdup("Memory");
    defval->type         = String;
    defval->vals.string  = EG_strdup("1600MB");
    defval->lfixed       = Change;

    /*! \page aimInputsZAERO
     * - <B> Memory = "1600MB" </B> <br>
     * Maximum memory in terms of megabytes that is allocable by ZAERO from
     * the heap space.
     */

  } else if (index == inSmart_Restart) {
    *ainame              = EG_strdup("Smart_Restart");
    defval->type         = Boolean;
    defval->vals.integer = (int) true;

    /*! \page aimInputsZAERO
     * - <B> Smart_Restart = True </B> <br>
     * If True, zaeroAIM will try to detect whether the ZAERO restart capability
     * can be used and configure the ZAERO input to load existing AIC matrices.
     * If False, new AIC matrices are always generated.
     */

  } else if (index == inEcho) {
    *ainame              = EG_strdup("Echo");
    defval->type         = String;
    defval->vals.string  = EG_strdup("sort");
    defval->lfixed       = Change;

    /*! \page aimInputsZAERO
     * - <B> Echo = "sort" </B> <br>
     * Controls echo (printout) of the Bulk Data Section
     */

  } else if (index == inOutput) {
    *ainame              = EG_strdup("Output");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> Output </B> <br>
     * Output tuple used to define analysis/case outputs for plotting, etc.
     */

   } else if (index == inHFG) {
     *ainame              = EG_strdup("HFG");
     defval->type         = Tuple;
     defval->nullVal      = IsNull;
     defval->lfixed       = Change;
     defval->dim          = Vector;

     /*! \page aimInputsZAERO
      * - <B> HFG </B> <br>
      * JSON dictionary used to define HFG module data
      */

   } else if (index == inUAIC) {
    *ainame              = EG_strdup("UAIC");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> UAIC </B> <br>
     * UAIC tuple used to define UAIC configurations for unsteady aerodynamics data generation, see \ref zaeroUAIC for additional details.
     */

  } else if (index == inSpline) {
    *ainame              = EG_strdup("Spline");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> Spline </B> <br>
     * JSON dictionary used to define SPLINE module data, see \ref zaeroSpline for additional details.
     */

  } else if (index == inVLM_Surface) {
    *ainame              = EG_strdup("VLM_Surface");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    //defval->units        = NULL;
    defval->dim          = Vector;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;

    /*! \page aimInputsZAERO
     * - <B>VLM_Surface = NULL </B> <br>
     * Vortex lattice method tuple input,  see \ref vlmSurface for additional details.
     */

  }  else if (index == inVLM_Control) {
    *ainame              = EG_strdup("VLM_Control");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    //defval->units        = NULL;
    defval->dim          = Vector;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;

    /*! \page aimInputsZAERO
     * - <B>VLM_Control = NULL </B> <br>
     * Vortex lattice method control surface tuple input, see \ref vlmControl for additional details.
     */

  } else if (index == inTrim_Variable) {
    *ainame              = EG_strdup("Trim_Variable");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsZAERO
     * - <B> Trim_Variable </B> <br>
     * Trim_Variable tuple used to define Trim variables and/or constraints, see \ref zaeroTRIMVAR for additional details.
     */

  } else if (index == inReferenceArea) {
    *ainame              = EG_strdup("ReferenceArea");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;
    if (units != NULL && units->area != NULL) {
        AIM_STRDUP(defval->units, units->area, aimInfo, status);
    }

    /*! \page aimInputsZAERO
     * - <B>ReferenceArea = NULL </B> <br>
     * This sets the reference area for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeZAERO) "capsReferenceArea" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */

  } else if (index == inReferenceChord) {
    *ainame              = EG_strdup("ReferenceChord");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;
    if (units != NULL && units->length != NULL) {
        AIM_STRDUP(defval->units, units->length, aimInfo, status);
    }

    /*! \page aimInputsZAERO
     * - <B>ReferenceChord = NULL </B> <br>
     * This sets the reference chord for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeZAERO) "capsReferenceChord" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */

  } else if (index == inReferenceSpan) {
    *ainame              = EG_strdup("ReferenceSpan");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;
    if (units != NULL && units->length != NULL) {
        AIM_STRDUP(defval->units, units->length, aimInfo, status);
    }

    /*! \page aimInputsZAERO
     * - <B>ReferenceSpan = NULL </B> <br>
     * This sets the reference span for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeZAERO) "capsReferenceSpan" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */

  } else if (index == inMoment_Center) {
    *ainame              = EG_strdup("Moment_Center");
    defval->type          = Double;
    defval->dim           = 1;
    defval->length        = 3;
    defval->nrow          = 3;
    defval->ncol          = 1;
    AIM_ALLOC(defval->vals.reals, defval->length, double, aimInfo, status);
    defval->vals.reals[0] = 0.0;
    defval->vals.reals[1] = 0.0;
    defval->vals.reals[2] = 0.0;
    defval->nullVal       = IsNull;
    defval->lfixed        = Fixed;
    if (units != NULL && units->length != NULL) {
        AIM_STRDUP(defval->units, units->length, aimInfo, status);
    }

    /*! \page aimInputsZAERO
     * - <B>Moment_Center = [0.0, 0.0, 0.0] (NULL)</B> <br>
     * Array values correspond to the x, y, and z center of gravity (CG) locations [meter].
     * Alternatively, the geometry (body) attributes (see \ref attributeZAERO) "capsReferenceX", "capsReferenceY",
     * and "capsReferenceZ" may be used to specify the center of gravity, respectively
     * (note: values set through the AIM input will supersede the attribution values).
     */

  } else if (index == inMassPropLink) {
    *ainame              = EG_strdup("MassPropLink");
    defval->type         = Pointer;
    defval->nullVal      = IsNull;
    AIM_STRDUP(defval->units, "feaMassPropStruct", aimInfo, status);

    /*! \page aimInputsZAERO
     * - <B>MassPropLink = NULL</B> <br>
     * Mass properties linked from structural analysis for eigen value analysis
     * Must be in units of kg, m, and kg*m^2 if unitSystem (see \ref aimUnitsZAERO) is not specified.
     */

  } else {
    AIM_ERROR(aimInfo, "Unknown input index $%d", index);
    status = CAPS_RANGEERR;
    goto cleanup;
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
  return CAPS_SUCCESS;

}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
  // Function return flag
  int status = CAPS_SUCCESS;
  int i;
  zaeroUAICStruct *uaic=NULL;

  aimStorage *zaeroInstance = (aimStorage *)instStore;
  AIM_NOTNULL(zaeroInstance, aimInfo, status);
  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get project name
  zaeroInstance->smartRestart = aimInputs[inSmart_Restart-1].vals.integer;

  if (zaeroInstance->smartRestart) {
    printf("\n'Smart_Restart' is ON\n");
  }
  else {
    printf("\n'Smart_Restart' is OFF\n");
  }

  /* Per the manual:
   * "In the following conditions, the restart process becomes inapplicable and
   *  the new AICs must be computed [...]
   *  - Any changes in the CAERO7 and/or BODY7 and their associated
   *    bulk data cards
   *  - Any changes in the MKAEROZ bulk data card"
   *
   * Therefore aicFilename is removed if any of the following has changed:
   *  - the MKAEROZ data (UAIC module)
   *  - the CAERO7 data
   *  - the BODY7
   */
  if (zaeroInstance->smartRestart == (int)false ||
      aim_newAnalysisIn(aimInfo, inUAIC       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, inVLM_Surface) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, inVLM_Control) == CAPS_SUCCESS ||
      aim_newGeometry(aimInfo) == CAPS_SUCCESS) {
    for (i = 0; i < zaeroInstance->zaeroProblem.numUAICs; i++) {

      uaic = &zaeroInstance->zaeroProblem.UAICs[i];

      status = aim_rmFile(aimInfo, uaic->aicFilename);
      AIM_STATUS(aimInfo, status);
    }
  }

  status = _createVLM(aimInfo, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _getZaeroProblemData(aimInfo, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _setupArtifacts(aimInfo, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  // Mass properties are propagated and must be done last
  status = _getMassProp(aimInfo, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  if (zaeroInstance->smartRestart == (int)true) {
    for (i = 0; i < zaeroInstance->zaeroProblem.numUAICs; i++) {
      uaic = &zaeroInstance->zaeroProblem.UAICs[i];
      if (aim_isFile(aimInfo, uaic->aicFilename) == CAPS_SUCCESS) {
        uaic->saveFlag = 1; // ACQUIRE
      }
    }
  }

  // Set file format type
  if        (strcasecmp(aimInputs[inFile_Format-1].vals.string, "Small") == 0) {
    zaeroInstance->feaFormatType = SmallField;
  } else if (strcasecmp(aimInputs[inFile_Format-1].vals.string, "Large") == 0) {
    zaeroInstance->feaFormatType = LargeField;
  } else if (strcasecmp(aimInputs[inFile_Format-1].vals.string, "Free") == 0)  {
    zaeroInstance->feaFormatType = FreeField;
  } else {
    AIM_ERROR(aimInfo, "Unrecognized \"File_Format\", valid choices are [Small, Large, or Free].");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

cleanup:
  return status;
}

// ********************** AIM Function Break *****************************
/* aimPreAnalysis: Generate Input File(s) */
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
  int status; // Status return

  FILE *inputFile = NULL; // input file to be generated

  const aimStorage *zaeroInstance = NULL;

  // Get pointer to current instance
  zaeroInstance = (const aimStorage*)instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

// #ifdef DEBUG
//   printf("\n zaeroAIM/aimPreAnalysis zaeroInstances = %d!\n", iIndex);
// #endif

  // open input file to be generated
  inputFile = aim_fopen(aimInfo, zaeroInstance->artifacts.input, "w");
  if (inputFile == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s", aimInfo, zaeroInstance->artifacts.input);
    status = CAPS_IOERR;
    goto cleanup;
  }

  status = _writeExecutiveControlSection(aimInfo, inputFile, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeCaseControlSection(aimInfo, inputFile, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = _writeBulkDataSection(aimInfo, inputFile, aimInputs, zaeroInstance);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;

cleanup:

  if (inputFile != NULL) {
    fclose(inputFile);
  }

  return status;
}

// ********************** AIM Function Break *****************************
/* aimPostAnalysis: Perform any processing after the Analysis is run */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;

  return status;
}

// ********************** AIM Function Break *****************************
/* aimOutputs: Output Information for the AIM */
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
// #ifdef DEBUG
//   printf(" zaeroAIM/aimOutputs zaeroInstances = %d  index  = %d!\n", iIndex, index);
// #endif

//  *aoname = EG_strdup("zaeroAIMout");
//  if (*aoname == NULL) return EGADS_MALLOC;
//  form->type = Double;

  return CAPS_SUCCESS;
}

// ********************** AIM Function Break *****************************
/* aimCalcOutput: Calculate/Retrieve Output Information */
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, /*@unused@*/ int index,
                  /*@unused@*/ capsValue *val)
{
// #ifdef DEBUG
//   int        status;
//   const char *name;
//
//   status = aim_getName(aimInfo, index, ANALYSISOUT, &name);
//   printf(" zaeroAIM/aimCalcOutput zaeroInstances = %d  index = %d %s %d!\n",
//          iIndex, index, name, status);
// #endif

//  *errors = NULL;
//  if ((iIndex < 0) || (iIndex >= numInstances)) return CAPS_BADINDEX;
//  val->vals.real = 3.1415926;

  return CAPS_SUCCESS;
}

// ********************** AIM Function Break *****************************
/* aimCleanup: Free up the AIMs storage */
void aimCleanup(void *instStore)
{
  // #ifdef DEBUG
  //   printf(" zaeroAIM/aimCleanup!\n");
  // #endif

  aimStorage *zaeroInstance = NULL;

  // Get pointer to current instance
  zaeroInstance = (aimStorage*)instStore;

  // Clean up zaeroInstances data
  destroy_aimStorage(zaeroInstance);

  AIM_FREE(zaeroInstance);

}
