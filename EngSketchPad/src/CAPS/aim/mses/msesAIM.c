/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             MESS AIM
 *
 *      Written by Marshall Galbraith and Robert Haimes MIT
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewMSES MSES AIM Overview
 *
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (through input
 * files) with the airfoil analysis tool MSES. MSES is not open-source and not freelay available. However,
 * a 'lite' version of MSES is provided with EngSketchPad that supports analysis of a single airfoil.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsMSES and \ref aimOutputsMSES, respectively.
 *
 * The AIM preAnalysis generates the mesh by calling the 'mset' executable, and generats an msesInput.txt file
 * containing the instructions for executing mses.
 *
 * The MSES AIM can automatically execute MSES, with details provided in \ref aimExecuteMSES.
 *
 * \section assumptionsMSES Assumptions
 * MSES inherently assumes the airfoil cross-section is in the x-y plane, if it isn't an attempt is made
 * to automatically rotate the provided body.
 *
 * Within <b> OpenCSM</b>, there are a number of airfoil generation UDPs (User Defined Primitives). These include NACA 4
 * series, a more general NACA 4/5/6 series generator, Sobieczky's PARSEC parameterization and Kulfan's CST
 * parameterization. All of these UDPs generate <b> EGADS</b> <em> FaceBodies</em> where the <em>Face</em>'s underlying
 * <em>Surface</em>  is planar and the bounds of the <em>Face</em> is a closed set of <em>Edges</em> whose
 * underlying <em>Curves</em> contain the airfoil shape. In all cases there is a <em>Node</em> that represents
 * the <em>Leading Edge</em> point and one or two <em>Nodes</em> at the <em>Trailing Edge</em> -- one if the
 * representation is for a sharp TE and the other if the definition is open or blunt. If there are 2 <em>Nodes</em>
 * at the back, then there are 3 <em>Edges</em> all together and closed, even though the airfoil definition
 * was left open at the TE. All of this information will be used to automatically fill in the MSES geometry
 * description.
 *
 * It should be noted that general construction in either <b> OpenCSM</b> or even <b> EGADS</b> will be supported
 * as long as the topology described above is used. But care should be taken when constructing the airfoil shape
 * so that a discontinuity (i.e.,  simply <em>C<sup>0</sup></em>) is not generated at the <em>Node</em> representing
 * the <em>Leading Edge</em>. This can be done by splining the entire shape as one and then intersecting the single
 *  <em>Edge</em> to place the LE <em>Node</em>.
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "egads_dot.h"
#include "aimUtil.h"

#include "miscUtils.h" // Bring in miscellaneous utilities
#include "meshUtils.h"
#include "vlmUtils.h"
#include "cfdUtils.h"

#include "msesUtils.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#define strtok_r   strtok_s
#endif

#define NUMPOINT  201
#define NMODE  40 // number of shape modes for MSES (must be even, max 40)

#define PI        3.1415926535897931159979635

#define MXCHAR  255

#define MIN(A,B)         (((A) < (B)) ? (A) : (B))
#define MAX(A,B)         (((A) < (B)) ? (B) : (A))

enum aimInputs
{
  inMach = 1,                      /* index is 1-based */
  inRe,
  inAlpha,
  inCL,
  inAcrit,
  inISMOM,
  inGridAlpha,
  inCoarse_Iteration,
  inFine_Iteration,
  inxGridRange,
  inyGridRange,
  inDesign_Variable,
  inCheby_Modes,
  NUMINPUT = inCheby_Modes          /* Total number of inputs */
};

enum aimOutputs
{
  outAlpha = 1,                     /* index is 1-based */
  outCL,
  outCD,
  outCD_p,
  outCD_v,
  outCD_w,
  outCM,
  outCheby_Modes,
  NUMOUTPUT = outCheby_Modes        /* Total number of outputs */
};

#define GCON_ALPHA 5
#define GCON_CL    6

typedef struct {

  capsValue Alpha;
  capsValue CL;
  capsValue CD;
  capsValue CDp;
  capsValue CDv;
  capsValue CDw;
  capsValue CM;
  capsValue Cheby_Modes;

  // Design information
  cfdDesignStruct design;

  int numBody;
  double **xCoord;
  double **yCoord;

} aimStorage;


static int destroy_aimStorage(aimStorage *msesInstance)
{
  int i;
  aim_freeValue(&msesInstance->Alpha);
  aim_freeValue(&msesInstance->CL);
  aim_freeValue(&msesInstance->CD);
  aim_freeValue(&msesInstance->CDp);
  aim_freeValue(&msesInstance->CDv);
  aim_freeValue(&msesInstance->CDw);
  aim_freeValue(&msesInstance->CM);
  aim_freeValue(&msesInstance->Cheby_Modes);

  // Design information
  (void) destroy_cfdDesignStruct(&msesInstance->design);

  for (i = 0; i < msesInstance->numBody; i++) {
    AIM_FREE(msesInstance->xCoord[i]);
    AIM_FREE(msesInstance->yCoord[i]);
  }
  msesInstance->numBody = 0;
  AIM_FREE(msesInstance->xCoord);
  AIM_FREE(msesInstance->yCoord);

  return CAPS_SUCCESS;
}


/* ********************** Exposed AIM Functions ***************************** */
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int status = CAPS_SUCCESS;
  aimStorage *msesInstance=NULL;

#ifdef DEBUG
  printf("\n msesAIM/aimInitialize   inst = %d!\n", inst);
#endif

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUTPUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 0;
  *fnames  = NULL;
  *franks  = NULL;
  *fInOut  = NULL;

  // Allocate msesInstance
  AIM_ALLOC(msesInstance, 1, aimStorage, aimInfo, status);
  *instStore = msesInstance;

/*@-uniondef@*/
  aim_initValue(&msesInstance->Alpha);
  aim_initValue(&msesInstance->CL);
  aim_initValue(&msesInstance->CD);
  aim_initValue(&msesInstance->CDp);
  aim_initValue(&msesInstance->CDv);
  aim_initValue(&msesInstance->CDw);
  aim_initValue(&msesInstance->CM);
  aim_initValue(&msesInstance->Cheby_Modes);
/*@+uniondef@*/

  // Design information
  status = initiate_cfdDesignStruct(&msesInstance->design);
  AIM_STATUS(aimInfo, status);

  msesInstance->numBody = 0;
  msesInstance->xCoord = NULL;
  msesInstance->yCoord = NULL;

cleanup:
  return status;
}


int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
  int status = CAPS_SUCCESS;
  /*! \page aimInputsMSES AIM Inputs
   * The following list outlines the xFoil inputs along with their default values available
   * through the AIM interface.
   */

#ifdef DEBUG
  printf(" msesAIM/aimInputs  index = %d!\n", index);
#endif

  if (index == inMach) {
    *ainame           = EG_strdup("Mach");
    defval->type      = Double;
    defval->vals.real = 0.0;
    defval->nullVal   = IsNull;

    /*! \page aimInputsMSES
     * - <B> Mach = NULL </B> <br>
     *  Mach number.
     */

  } else if (index == inRe) {
    *ainame           = EG_strdup("Re");
    defval->type      = Double;
    defval->vals.real = 0.0;

    /*! \page aimInputsMSES
     * - <B> Re = 0.0 </B> <br>
     *  Reynolds number. Use 0.0 for an inviscid calculation.
     */
  } else if (index == inAlpha) {
    *ainame           = EG_strdup("Alpha");
    defval->type      = Double;
    defval->dim       = Scalar;
    defval->vals.real = 0.0;
    defval->nullVal   = IsNull;
    //defval->units     = EG_strdup("degree");

    /*! \page aimInputsMSES
     * - <B> Alpha = NULL </B> <br>
     *  Angle of attack [degree].
     */
  } else if (index == inCL) {
    *ainame           = EG_strdup("CL");
    defval->type      = Double;
    defval->dim       = Scalar;
    defval->vals.real = 0.0;
    defval->nullVal   = IsNull;

    /*! \page aimInputsMSES
     * - <B> CL = NULL </B> <br>
     *  Prescribed coefficient of lift.
     */
  } else if (index == inAcrit) {
    *ainame           = EG_strdup("Acrit");
    defval->type      = Double;
    defval->dim       = Scalar;
    defval->vals.real = 9.0;
    defval->nullVal   = NotAllowed;

    /*! \page aimInputsMSES
     * - <B> Acrit =  9.0</B> <br>
     *  Critical amplification factor "n" for the e^n envelope transition model. 9.0 is the standard model.
     */

  } else if (index == inISMOM) {
    *ainame              = EG_strdup("ISMOM");
    defval->type         = Integer;
    defval->dim          = Scalar;
    defval->vals.integer = 4;
    defval->nullVal      = NotAllowed;

    /*! \page aimInputsMSES
     * - <B> ISMOM = 4</B> <br>
     *  MSES ISMOM input to select the momentum equation. Valid inputs: 1, 2, 3, 4
     */

  } else if (index == inGridAlpha) {
    *ainame           = EG_strdup("GridAlpha");
    defval->type      = Double;
    defval->dim       = Scalar;
    defval->vals.real = 0.0;
    defval->nullVal   = NotAllowed;

    /*! \page aimInputsMSES
     * - <B> GridAlpha =  0.0</B> <br>
     *  Angle of attack used to generate the grid.
     */
  } else if (index == inCoarse_Iteration) {
    *ainame              = EG_strdup("Coarse_Iteration");
    defval->type         = Integer;
    defval->vals.integer = 20;

    /*! \page aimInputsMSES
     * - <B> Coarse_Iteration = 20</B> <br>
     *  Maximum number of coarse mesh iterations (can help convergence).
     */
  } else if (index == inFine_Iteration) {
    *ainame              = EG_strdup("Fine_Iteration");
    defval->type         = Integer;
    defval->vals.integer = 200;

    /*! \page aimInputsMSES
     * - <B> Fine_Iteration = 200</B> <br>
     *  Maximum number of fine mesh iterations.
     */

  } else if (index == inxGridRange) {
    *ainame              = EG_strdup("xGridRange");
    defval->type         = Double;
    defval->dim          = Vector;
    defval->nrow         = 2;
    defval->lfixed       = Fixed;
    defval->sfixed       = Fixed;
    AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
    defval->vals.reals[0] = -1.75;
    defval->vals.reals[1] =  2.75;

    /*! \page aimInputsMSES
     * - <B> xGridRange = [-1.75, 2.75]</B> <br>
     *  x-min and x-max values for the grid domain size.
     */

  } else if (index == inyGridRange) {
    *ainame              = EG_strdup("yGridRange");
    defval->type         = Double;
    defval->dim          = Vector;
    defval->nrow         = 2;
    defval->lfixed       = Fixed;
    defval->sfixed       = Fixed;
    AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
    defval->vals.reals[0] = -2.0;
    defval->vals.reals[1] =  2.5;

    /*! \page aimInputsMSES
     * - <B> xGridRange = [-2.0, 2.5]</B> <br>
     *  x-min and x-max values for the grid domain size.
     */

  } else if (index == inDesign_Variable) {
    *ainame              = EG_strdup("Design_Variable");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    //defval->units        = NULL;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsMSES
     * - <B> Design_Variable = NULL</B> <br>
     * The design variable tuple is used to input design variable information for the model optimization, see \ref cfdDesignVariable for additional details.
     * Must be NULL of Cheby_Modes is not NULL.
     */

  } else if (index == inCheby_Modes) {
    *ainame              = EG_strdup("Cheby_Modes");
    defval->type         = Double;
    defval->lfixed       = Fixed;
    defval->dim          = Vector;
    defval->nrow         = NMODE;
    AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
    memset(defval->vals.reals, 0, defval->nrow*sizeof(double));
    defval->nullVal      = IsNull;

    /*! \page aimInputsMSES
     * - <B> Cheby_Modes = NULL</B> <br>
     * Chebyshev shape mode values for shape optimization (fixed length 40).
     * Must be NULL if Design_Variable is not NULL
     */

  }

  status = CAPS_SUCCESS;
cleanup:
  return status;
}


int aimPreAnalysis(/*@unused@*/ void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
  int status; // Function return status

  int i, ibody; // Indexing
  int gcon=0;

  // Bodies
  const char *intents;
  int numBody;
  ego *bodies = NULL;

  vlmSectionStruct vlmSection;
  char command[PATH_MAX];

  aimStorage *msesInstance=NULL;

  // File I/O
  FILE *fp = NULL;
  char inputMSES[]     = "msesInput.txt";
  char inputMSET[]     = "msetInput.txt";
  char msesFilename[]  = "mses.airfoil";
  char bladeFilename[] = "blade.airfoil";
  char modesFilename[] = "modes.airfoil";
  char paramsFilename[] = "params.airfoil";

  if (aimInputs == NULL) {
#ifdef DEBUG
      printf("\tmsesAIM/aimPreAnalysis aimInputs == NULL!\n");
#endif
      return CAPS_NULLVALUE;
  }

  msesInstance = (aimStorage*)instStore;
  destroy_aimStorage(msesInstance);


  if (aimInputs[inAlpha-1].nullVal == aimInputs[inCL-1].nullVal)
  {
    AIM_ERROR(aimInfo, "One of 'Alpha' and 'CL' inputs must be specified.");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (aimInputs[inMach-1].nullVal == IsNull)
  {
    AIM_ERROR(aimInfo, "'Mach' input must be specified.");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (aimInputs[inISMOM-1].vals.integer < 1 ||
      aimInputs[inISMOM-1].vals.integer > 4)
  {
    AIM_ERROR(aimInfo, "'ISMOM' must be in [1-4]: ISMOM = %d", aimInputs[inISMOM-1].vals.integer);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (aimInputs[inDesign_Variable-1].nullVal == NotNull &&
      aimInputs[inCheby_Modes-1].nullVal == NotNull)
  {
    AIM_ERROR(aimInfo, "Only one of 'Design_Variable' and 'Cheby_Modes' can be set");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Get design variables
  if (aimInputs[inDesign_Variable-1].nullVal == NotNull &&
      aim_newAnalysisIn(aimInfo, inDesign_Variable) == CAPS_SUCCESS) {
/*@-nullpass@*/
      status = cfd_getDesignVariable(aimInfo,
                                     aimInputs[inDesign_Variable-1].length,
                                     aimInputs[inDesign_Variable-1].vals.tuple,
                                     &msesInstance->design.numDesignVariable,
                                     &msesInstance->design.designVariable);
/*@+nullpass@*/
      AIM_STATUS(aimInfo, status);
  }

  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

  if (numBody == 0 || bodies == NULL) {
    AIM_ERROR(aimInfo, "No Bodies!");
    return CAPS_SOURCEERR;
  }

  status = initiate_vlmSectionStruct(&vlmSection);
  AIM_STATUS(aimInfo, status);

  AIM_ALLOC(msesInstance->xCoord, numBody, double*, aimInfo, status);
  AIM_ALLOC(msesInstance->yCoord, numBody, double*, aimInfo, status);
  msesInstance->numBody = numBody;

  // Accumulate cross coordinates of airfoil and write out data file
  for (ibody = 0; ibody < numBody; ibody++) {

    // Open and write the input to control the MSES session
    fp = aim_fopen(aimInfo, bladeFilename, "w");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Unable to open file %s\n!", bladeFilename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    vlmSection.ebody = bodies[ibody];

    status = finalize_vlmSectionStruct(aimInfo, &vlmSection);
    AIM_STATUS(aimInfo, status);

    status = vlm_getSectionCoord(aimInfo,
                                 &vlmSection,
                                 (int) true, // Normalize by chord (true/false)
                                 NUMPOINT,
                                 &msesInstance->xCoord[ibody],  //[numPoint]
                                 &msesInstance->yCoord[ibody]); //[numPoint] for upper and lower surface
    AIM_STATUS(aimInfo, status);

    fprintf(fp,"capsBody_%d\n",ibody+1);
    fprintf(fp, "%lf %lf %lf %lf\n", aimInputs[inxGridRange-1].vals.reals[0],
                                     aimInputs[inxGridRange-1].vals.reals[1],
                                     aimInputs[inyGridRange-1].vals.reals[0],
                                     aimInputs[inyGridRange-1].vals.reals[1]);
    for( i = 0; i < NUMPOINT; i++) {
      fprintf(fp, "%lf %lf\n", msesInstance->xCoord[ibody][i],
                               msesInstance->yCoord[ibody][i]);
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
    fp = NULL;

  } // End body loop

  // This writes the mses.xxx file
  // Described in MSES/pdf/mses.pdf
  fp = aim_fopen(aimInfo, msesFilename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file %s\n!", msesFilename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  if (aimInputs[inAlpha-1].nullVal == NotNull)
    gcon = GCON_ALPHA; // specify Alpha
  else
    gcon = GCON_CL; // specify CL

  fprintf(fp, "3 4 5 7 10 15 20 \n"); // First row are the variables, are unlikely to change
                                // These come from Drela's examples, see Pages 5,17,18,19
                                // The 20 enables the sensitivities
  fprintf(fp, "3 4 %d 7 15 17 20\n", gcon);  // Second row are constraints, also unlikely to change
                                       // In general they match the variables, see Pages 6,7,8,17,18,19
  fprintf(fp, "%lf %lf %lf     | MACHIN  CLIFIN ALFAIN\n", aimInputs[inMach-1].vals.real,
                                                           aimInputs[inCL-1].vals.real,
                                                           aimInputs[inAlpha-1].vals.real );
  fprintf(fp, "%d 2            | ISMOM   IFFBC\n", aimInputs[inISMOM-1].vals.integer); // These have to do with flow properties (Page 7) and should probably be exposed to the user
  fprintf(fp, "%lf %lf         | REYNIN  ACRIT\n", aimInputs[inRe-1].vals.real,
                                                   aimInputs[inAcrit-1].vals.real);
  fprintf(fp, "1.0 1.0         | XTRS    XTRP\n");    // Forced transition location on top and bottom, x/c
  fprintf(fp, "0.98 1.0        | MCRIT   MUCON\n");   //these have to do with shocks and should also be exposed as settings
  fprintf(fp, "1 1             | ISMOVE  ISPRES\n"); // These are ignored, but must be here
  fprintf(fp, "%d 0            | NMOD    NPOS\n", NMODE);   // Number of mode variables and position variables, could be a user setting

  if (fp != NULL) fclose(fp);
  fp = NULL;

  // This writes the modes file that determines the geometry
  // See pages 41-12 MSES/pdf/mses.pdf
  // Format is: Variable   Mode-Shape   (a multi element thing that is ignored)  Mode-lower-bound  Mode-upper-bound  (always 1)
  fp = aim_fopen(aimInfo, modesFilename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file %s\n!", modesFilename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  for (i = 0; i < NMODE/2; i++)
    fprintf(fp, "%d   %d   1.0   0.0    1.0   1\n", i+1, NMODE/2+1+i); // Upper Surface Modes
  for (i = 0; i < NMODE/2; i++)
    fprintf(fp, "%d   %d   1.0   0.0   -1.0   1\n", NMODE/2+1+i, NMODE/2+1+i); // Lower surface modes

  if (fp != NULL) fclose(fp);
  fp = NULL;

  if (aimInputs[inCheby_Modes-1].nullVal == NotNull) {
    fp = aim_fopen(aimInfo, paramsFilename, "w");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Unable to open file %s\n!", paramsFilename);
      status = CAPS_IOERR;
      goto cleanup;
    }
    fprintf(fp, "%d   0\n", NMODE);
    for (i = 0; i < NMODE; i++)
      fprintf(fp, "%16.12e\n", aimInputs[inCheby_Modes-1].vals.reals[i]); // Surface Modes

    fprintf(fp, "%lf %lf %lf \n", aimInputs[inAlpha-1].vals.real,
                                  aimInputs[inCL-1].vals.real,
                                  aimInputs[inMach-1].vals.real );
    fprintf(fp, "%lf\n", aimInputs[inRe-1].vals.real);

    if (fp != NULL) fclose(fp);
    fp = NULL;
  } else {
    status = aim_rmFile(aimInfo, paramsFilename);
    AIM_STATUS(aimInfo, status);
  }

  // Run MSET to set up the grid
  // This is where the interactive window comes up, we really need to figure out how to turn it off to make this practical
  fp = aim_fopen(aimInfo, inputMSET, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file %s\n!", inputMSET);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // The normal sequence is 1-2-3-4-0
  fprintf(fp, "1\n"); // Generates streamlines
  fprintf(fp, "%lf\n", aimInputs[inGridAlpha-1].vals.real); // Sets the alpha of the grid (this is not the flow alpha)
  fprintf(fp, "2\n"); // Grid spacing
  fprintf(fp, "\n"); // No changes desired, honestly not sure what these settings do
  fprintf(fp, "3\n"); // Does a smoothing
  fprintf(fp, "4\n"); // write the mdat file that MSES uses
  fprintf(fp, "0\n"); // quits mset
  fprintf(fp, "\n");
  fflush(fp);
  if (fp != NULL) fclose(fp);
  fp = NULL;
  snprintf(command, PATH_MAX, "mset airfoil noplot < %s > msetOutput.txt", inputMSET);
  status = aim_system(aimInfo, "", command);
  AIM_STATUS(aimInfo, status, "Failed to execute: %s", command);


  // Open and write the input to control the MSES session
  fp = aim_fopen(aimInfo, inputMSES, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file %s\n!", inputMSES);
    status = CAPS_IOERR;
    goto cleanup;
  }

  if (aimInputs[inCoarse_Iteration-1].vals.integer != 0)
    fprintf(fp, "-%d\n", abs(aimInputs[inCoarse_Iteration-1].vals.integer));
  fprintf(fp, "+%d\n", abs(aimInputs[inFine_Iteration-1].vals.integer));
  fprintf(fp, "0\n"); // Terminates mses

  if (fp != NULL) fclose(fp);
  fp = NULL;

  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);

  return status;
}


int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteMSES AIM Execution
   *
   * If auto execution is enabled when creating an MSES AIM,
   * the AIM will execute MSES just-in-time with the command line:
   *
   * \code{.sh}
   * mses airfoil < msesInput.txt > msesOutput.txt
   * \endcode
   *
   * where preAnalysis generated the file "msesInput.txt" which contains the input information.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an MSES AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or MSES can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * mses.preAnalysis()
   *
   * print ("\n\nRunning......")
   * mses.system("mses airfoil < msesInput.txt > msesOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * mses.postAnalysis()
   * \endcode
   */

  *state = 0;
  return aim_system(aimInfo, NULL, "mses airfoil < msesInput.txt > msesOutput.txt");
}


int aimPostAnalysis(void *instStore, void *aimInfo,
                    /*@unused@*/ int restart, capsValue *aimInputs)
{
  int status = CAPS_SUCCESS;
  int i, j, k, idv, igv, is, nis, ib, im, jm, index, nderiv;
  int numFunctional, irow, icol, gcon;
  vlmSectionStruct *vlmSections=NULL;
  aimStorage *msesInstance=NULL;

  msesSensx *sensx;
  double *xyz=NULL, coord[3], data[18], tm, tp, ism_dot[9], isp_dot[9];
  int sizes[2] = {NUMPOINT,0};

  double *M=NULL, *rhs=NULL, *dmod_dvar=NULL, ds, cl_alfa=1;

  int ibody;
  double **dxyz=NULL, *dx_dvar=NULL, *dy_dvar=NULL, functional_dvar;
  ego context=NULL, *curves=NULL;
  capsValue **values=NULL, *geomInVal;

  // Bodies
  const char *intents, *name;
  int numBody=0;
  ego *bodies = NULL;

  AIM_NOTNULL(instStore, aimInfo, status);
  AIM_NOTNULL(aimInputs, aimInfo, status);

  msesInstance = (aimStorage*)instStore;

  // check the MSES output file
//  if (aim_isFile(aimInfo, "mdat.airfoil") != CAPS_SUCCESS) {
//    AIM_ERROR(aimInfo, "mses execution did not produce mdat.airfoil");
//    return CAPS_EXECERR;
//  }
//
//  status = msesMdatRead(aimInfo, "mdat.airfoil", &mdat);
//  AIM_STATUS(aimInfo, status);
//
  if (aimInputs[inAlpha-1].nullVal == NotNull)
    gcon = GCON_ALPHA; // specify Alpha
  else
    gcon = GCON_CL; // specify CL

  if (aim_isFile(aimInfo, "sensx.airfoil") != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "mses execution did not produce sensx.airfoil!");
    return CAPS_EXECERR;
  }

  // read in the sensx.airfoil file
  status = msesSensxRead(aimInfo, "sensx.airfoil", &sensx);
  AIM_STATUS(aimInfo, status);

  numFunctional = 7;
  AIM_ALLOC(values, numFunctional, capsValue*, aimInfo, status);
  values[0] = &msesInstance->Alpha;
  values[1] = &msesInstance->CL;
  values[2] = &msesInstance->CD;
  values[3] = &msesInstance->CDp;
  values[4] = &msesInstance->CDv;
  values[5] = &msesInstance->CDw;
  values[6] = &msesInstance->CM;

  nderiv = 3;

  if (aimInputs[inDesign_Variable-1].nullVal == NotNull)
    nderiv += aimInputs[inDesign_Variable-1].length;
  else if (aimInputs[inCheby_Modes-1].nullVal == NotNull)
    nderiv += 1;

  for (i = 0; i < numFunctional; i++) {
    values[i]->type   = DoubleDeriv;
    values[i]->dim    = Scalar;
    values[i]->nderiv = nderiv;

    AIM_ALLOC(values[i]->derivs, nderiv, capsDeriv, aimInfo, status);
    for (j = 0; j < nderiv; j++) {
      values[i]->derivs[j].name    = NULL;
      values[i]->derivs[j].deriv   = NULL;
      values[i]->derivs[j].len_wrt = 1;
    }
    for (j = 0; j < 3; j++)
      AIM_ALLOC(values[i]->derivs[j].deriv, 1, double, aimInfo, status);

    if (gcon == GCON_ALPHA)
    { AIM_STRDUP(values[i]->derivs[0].name, "Alpha", aimInfo, status); }
    else
    { AIM_STRDUP(values[i]->derivs[0].name, "CL", aimInfo, status); }
    AIM_STRDUP(values[i]->derivs[1].name, "Mach", aimInfo, status);
    AIM_STRDUP(values[i]->derivs[2].name, "Re", aimInfo, status);
  }

/*@-nullderef@*/
  // Alpha -----------------------------------
  msesInstance->Alpha.vals.real = sensx->alfa * 180./PI;

  msesInstance->Alpha.derivs[0].deriv[0] = sensx->al_alfa;
  msesInstance->Alpha.derivs[1].deriv[0] = sensx->al_mach;
  msesInstance->Alpha.derivs[2].deriv[0] = sensx->al_reyn;

  // CL -----------------------------------
  msesInstance->CL.vals.real = sensx->cl;

  msesInstance->CL.derivs[0].deriv[0] = sensx->cl_alfa / 180. * PI;
  msesInstance->CL.derivs[1].deriv[0] = sensx->cl_mach;
  msesInstance->CL.derivs[2].deriv[0] = sensx->cl_reyn;

  // CD -----------------------------------
  msesInstance->CD.vals.real = sensx->cdv + sensx->cdw;

  msesInstance->CD.derivs[0].deriv[0] = (sensx->cdv_alfa + sensx->cdw_alfa) / 180. * PI;
  msesInstance->CD.derivs[1].deriv[0] =  sensx->cdv_mach + sensx->cdw_mach;
  msesInstance->CD.derivs[2].deriv[0] =  sensx->cdv_reyn + sensx->cdw_reyn;

  // CDp -----------------------------------
  msesInstance->CDp.vals.real = sensx->cdv + sensx->cdw - sensx->cdf;

  msesInstance->CDp.derivs[0].deriv[0] = (sensx->cdv_alfa + sensx->cdw_alfa - sensx->cdf_alfa) / 180. * PI;  // dCDp/dalpha
  msesInstance->CDp.derivs[1].deriv[0] =  sensx->cdv_mach + sensx->cdw_mach - sensx->cdf_mach;
  msesInstance->CDp.derivs[2].deriv[0] =  sensx->cdv_reyn + sensx->cdw_reyn - sensx->cdf_reyn;

  // CDv -----------------------------------
  msesInstance->CDv.vals.real = sensx->cdv;

  msesInstance->CDv.derivs[0].deriv[0] = sensx->cdv_alfa / 180. * PI;
  msesInstance->CDv.derivs[1].deriv[0] = sensx->cdv_mach;
  msesInstance->CDv.derivs[2].deriv[0] = sensx->cdv_reyn;

  // CDw -----------------------------------
  msesInstance->CDw.vals.real = sensx->cdw;

  msesInstance->CDw.derivs[0].deriv[0] = sensx->cdw_alfa / 180. * PI;
  msesInstance->CDw.derivs[1].deriv[0] = sensx->cdw_mach;
  msesInstance->CDw.derivs[2].deriv[0] = sensx->cdw_reyn;

  // CM -----------------------------------
  msesInstance->CM.vals.real = sensx->cm;

  msesInstance->CM.derivs[0].deriv[0] = sensx->cm_alfa / 180. * PI;
  msesInstance->CM.derivs[1].deriv[0] = sensx->cm_mach;
  msesInstance->CM.derivs[2].deriv[0] = sensx->cm_reyn;
  // --------------------------------------

  // Change alpha derivatives to CL derivatives when CL is specified
  if (gcon == GCON_CL)
  {
    cl_alfa = sensx->cl_alfa / 180. * PI;

    // Negative sign from Implicit Derivative (https://en.wikipedia.org/wiki/Implicit_function)
    msesInstance->Alpha.derivs[0].deriv[0] =              1./cl_alfa; // dalpha/dCL
    msesInstance->Alpha.derivs[1].deriv[0] = -sensx->cl_mach/cl_alfa; // dalpha/dMach
    msesInstance->Alpha.derivs[2].deriv[0] = -sensx->cl_reyn/cl_alfa; // dalpha/dRe

    msesInstance->CL.derivs[0].deriv[0]   = 1.0;     // dCL/dCL
    msesInstance->CD.derivs[0].deriv[0]  /= cl_alfa; // dCD/dCL
    msesInstance->CDp.derivs[0].deriv[0] /= cl_alfa; // dCDp/dCL
    msesInstance->CDv.derivs[0].deriv[0] /= cl_alfa; // dCDv/dCL
    msesInstance->CDw.derivs[0].deriv[0] /= cl_alfa; // dCDw/dCL
    msesInstance->CM.derivs[0].deriv[0]  /= cl_alfa; // dCM/dCL
  }


  // Cheby_Modes --------------------------------------
  if (aimInputs[inDesign_Variable-1].nullVal == NotNull) {
    msesInstance->Cheby_Modes.type   = DoubleDeriv;
    msesInstance->Cheby_Modes.nderiv = aimInputs[inDesign_Variable-1].length;

    AIM_ALLOC(msesInstance->Cheby_Modes.derivs, msesInstance->Cheby_Modes.nderiv, capsDeriv, aimInfo, status);
    for (j = 0; j < msesInstance->Cheby_Modes.nderiv; j++) {
      msesInstance->Cheby_Modes.derivs[j].name    = NULL;
      msesInstance->Cheby_Modes.derivs[j].deriv   = NULL;
      msesInstance->Cheby_Modes.derivs[j].len_wrt = 0;
    }
  } else {
    msesInstance->Cheby_Modes.type = Double;
  }
  msesInstance->Cheby_Modes.dim    = Vector;
  AIM_ALLOC(msesInstance->Cheby_Modes.vals.reals, sensx->nmod, double, aimInfo, status);
  msesInstance->Cheby_Modes.nrow   = sensx->nmod;
  msesInstance->Cheby_Modes.length = sensx->nmod;
  for (i = 0; i < sensx->nmod; i++)
    msesInstance->Cheby_Modes.vals.reals[i] = sensx->modn[i];
  // --------------------------------------
/*@+nullderef@*/

  // Only compute geometric sensitivities if requested
  if (aimInputs[inDesign_Variable-1].nullVal == NotNull) {

    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody == 0 || bodies == NULL) {
      AIM_ERROR(aimInfo, "No Bodies!");
      return CAPS_SOURCEERR;
    }

    AIM_ALLOC(vlmSections, numBody, vlmSectionStruct, aimInfo, status);
    for (ibody = 0; ibody < numBody; ibody++) {
      status = initiate_vlmSectionStruct(&vlmSections[ibody]);
      AIM_STATUS(aimInfo, status);
    }

    status = EG_getContext(bodies[0], &context);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(context, aimInfo, status);

    AIM_ALLOC(dxyz, numBody, double*, aimInfo, status);
    for (ibody = 0; ibody < numBody; ibody++) dxyz[ibody] = NULL;

    AIM_ALLOC(curves, numBody, ego, aimInfo, status);
    for (ibody = 0; ibody < numBody; ibody++) curves[ibody] = NULL;

    AIM_ALLOC(xyz, 3*NUMPOINT, double, aimInfo, status);

    for (ibody = 0; ibody < numBody; ibody++) {
      // get the tessellation of the body for sensitvities
      vlmSections[ibody].ebody = bodies[ibody];

      status = finalize_vlmSectionStruct(aimInfo, &vlmSections[ibody]);
      AIM_STATUS(aimInfo, status);
    }

    // allocate matrix data
    AIM_ALLOC(M, sensx->nmod*sensx->nmod, double, aimInfo, status);
    for (i = 0; i < sensx->nmod*sensx->nmod; i++) M[i] = 0;

    AIM_ALLOC(rhs, sensx->nmod, double, aimInfo, status);
    for (i = 0; i < sensx->nmod; i++) rhs[i] = 0;
    AIM_ALLOC(dmod_dvar, sensx->nmod, double, aimInfo, status);
    for (i = 0; i < sensx->nmod; i++) dmod_dvar[i] = 0;

    // construct the mass matrix for the projection
    // use simple trapezoidal integration
    for (im = 0; im < sensx->nmod; im++) {
      for (jm = 0; jm < sensx->nmod; jm++) {
        for (ibody = 0; ibody < numBody; ibody++) {

          nis = sensx->iteb[ibody] - sensx->ileb[ibody]+1;
          for (k = 0; k < 2; k++) {

            ib = 2*ibody+k;
            for (is = 0; is < nis-1; is++) {
              ds = sqrt(pow(sensx->xbi[ib][is+1]-sensx->xbi[ib][is],2.) +
                        pow(sensx->ybi[ib][is+1]-sensx->ybi[ib][is],2.));
              M[im*sensx->nmod + jm] += 0.5*(
                  sensx->xbi_mod[im][ib][is  ]*sensx->xbi_mod[jm][ib][is  ] +
                  sensx->ybi_mod[im][ib][is  ]*sensx->ybi_mod[jm][ib][is  ] +
                  sensx->xbi_mod[im][ib][is+1]*sensx->xbi_mod[jm][ib][is+1] +
                  sensx->ybi_mod[im][ib][is+1]*sensx->ybi_mod[jm][ib][is+1])*ds;
            }
          }
        }
      }
    }

    //#define PRINT_MATRIX
#ifdef PRINT_MATRIX
    for (im = 0; im < sensx->nmod; im++) {
      for (jm = 0; jm < sensx->nmod; jm++) {
        printf("%lf ",  M[im*sensx->nmod + jm]);
      }
      printf("\n");
    }
#endif

    // factorize the matrix in-place
    status = factorLU(sensx->nmod, M);
    AIM_STATUS(aimInfo, status);

    //#define PRINT_COORDINATES
#ifdef PRINT_COORDINATES
    printf("airfoil %d\n", sensx->nbl);
    for (ibody = 0; ibody < numBody; ibody++) {
      nis = sensx->iteb[ibody] - sensx->ileb[ibody]+1;
      for (k = 0; k < 2; k++) {
        printf("side %d\n", k);
        is = 2*ibody+k;
        for (int m = 0; m < nis; m++) {
          printf("%f %f\n", sensx->xbi[is][m],sensx->ybi[is][m]);
        }
      }
    }
    printf("done\n");
#endif

    /* allocate derivatives */
    for (idv = 0; idv < msesInstance->design.numDesignVariable; idv++) {

      name = msesInstance->design.designVariable[idv].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      AIM_STATUS(aimInfo, status);

      for (i = 0; i < numFunctional; i++) {
        AIM_STRDUP(values[i]->derivs[3+idv].name, name, aimInfo, status);

        AIM_ALLOC(values[i]->derivs[3+idv].deriv, geomInVal->length, double, aimInfo, status);
        values[i]->derivs[3+idv].len_wrt  = geomInVal->length;
        for (j = 0; j < geomInVal->length; j++)
          values[i]->derivs[3+idv].deriv[j] = 0;
      }

      AIM_STRDUP(msesInstance->Cheby_Modes.derivs[idv].name, name, aimInfo, status);
      AIM_ALLOC(msesInstance->Cheby_Modes.derivs[idv].deriv, sensx->nmod*geomInVal->length, double, aimInfo, status);
      msesInstance->Cheby_Modes.derivs[idv].len_wrt  = geomInVal->length;
      for (j = 0; j < sensx->nmod*geomInVal->length; j++)
        msesInstance->Cheby_Modes.derivs[idv].deriv[j] = 0;
    }

    /* set derivatives */
    for (idv = 0; idv < msesInstance->design.numDesignVariable; idv++) {

      name = msesInstance->design.designVariable[idv].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      AIM_STATUS(aimInfo, status);

      for (irow = 0; irow < geomInVal->nrow; irow++) {
        for (icol = 0; icol < geomInVal->ncol; icol++) {

          igv = geomInVal->ncol*irow + icol;

          // set the sensitivities of the spline fit for each blade
          for (ibody = 0; ibody < numBody; ibody++) {
            status = vlm_getSectionTessSens(aimInfo,
                                            &vlmSections[ibody],
                                            (int)true,
                                            name,
                                            irow+1, icol+1, // row, col
                                            NUMPOINT,
                                            &dx_dvar, &dy_dvar);
            AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", name);
            AIM_NOTNULL(dx_dvar, aimInfo, status);
            AIM_NOTNULL(dy_dvar, aimInfo, status);

            AIM_ALLOC(dxyz[ibody], 3*NUMPOINT, double, aimInfo, status);
            for( i = 0; i < NUMPOINT; i++) {
              dxyz[ibody][3*i+0] = dx_dvar[i];
              dxyz[ibody][3*i+1] = dy_dvar[i];
              dxyz[ibody][3*i+2] = 0.0;
              //printf("%lf %lf\n", dxyz[ibody][3*i+0], dxyz[ibody][3*i+1]);
            }
            AIM_FREE(dx_dvar);
            AIM_FREE(dy_dvar);
          }

          for (ibody = 0; ibody < numBody; ibody++) {
            // Create a spline representation of the airfoil coordinates
            // this must be don after calling aim_tessSensitivity as OpenCSM
            // calls EG_deleteObject(context) which removes the curve
            for( i = 0; i < NUMPOINT; i++) {
              xyz[3*i+0] = msesInstance->xCoord[ibody][i];
              xyz[3*i+1] = msesInstance->yCoord[ibody][i];
              xyz[3*i+2] = 0.0;
            }
            status = EG_approximate(context, 0, 1e-8, sizes, xyz, &curves[ibody]);
            AIM_STATUS(aimInfo, status);

            status = EG_approximate_dot(curves[ibody], 0, 1e-8, sizes, xyz, dxyz[ibody]);
            AIM_STATUS(aimInfo, status);
            AIM_FREE(dxyz[ibody]);
          }

          for (i = 0; i < sensx->nmod; i++) rhs[i] = 0;

          for (ibody = 0; ibody < numBody; ibody++) {

            nis = sensx->iteb[ibody] - sensx->ileb[ibody]+1;
            for (k = 0; k < 2; k++) {

              ib = 2*ibody+k;

              // get the t-value at is = 0
              coord[0] = sensx->xbi[ib][0];
              coord[1] = sensx->ybi[ib][0];
              coord[2] = 0;

              status = EG_invEvaluate(curves[ibody], coord, &tm, data);
              AIM_STATUS(aimInfo, status);

              status = EG_evaluate_dot(curves[ibody], &tm, NULL, data, ism_dot);
              AIM_STATUS(aimInfo, status);

              for (is = 0; is < nis-1; is++) {

                // get the t-value at is+1
                coord[0] = sensx->xbi[ib][is+1];
                coord[1] = sensx->ybi[ib][is+1];
                coord[2] = 0;

                status = EG_invEvaluate(curves[ibody], coord, &tp, data);
                AIM_STATUS(aimInfo, status);

                // get the spline sensitivity at is+1
                status = EG_evaluate_dot(curves[ibody], &tp, NULL, data, isp_dot);
                AIM_STATUS(aimInfo, status);

                //printf("%lf %lf %lf %lf\n", tp, isp_dot[0], isp_dot[1], isp_dot[2]);
                //printf("%lf %lf %lf %lf\n", tp, coord[0], coord[1], coord[2]);

                ds = sqrt(pow(sensx->xbi[ib][is+1]-sensx->xbi[ib][is],2.) +
                          pow(sensx->ybi[ib][is+1]-sensx->ybi[ib][is],2.));

                // integrate with trapezoidal rule
                for (im = 0; im < sensx->nmod; im++) {
                  rhs[im] += 0.5*(
                      sensx->xbi_mod[im][ib][is  ]*ism_dot[0] +
                      sensx->ybi_mod[im][ib][is  ]*ism_dot[1] +
                      sensx->xbi_mod[im][ib][is+1]*isp_dot[0] +
                      sensx->ybi_mod[im][ib][is+1]*isp_dot[1])*ds;
                }

                // cycle the data
                tm = tp;
                for (i = 0; i < 9; i++) ism_dot[i] = isp_dot[i];
              }
            }

            // free the curve
            EG_deleteObject(curves[ibody]);
          }

          // get the mode sensitivities w.r.t. design variables
          status = backsolveLU(sensx->nmod, M, rhs, dmod_dvar);
          AIM_STATUS(aimInfo, status);


          /* save of the mode sensitivities */
          for (j = 0; j < sensx->nmod; j++)
            msesInstance->Cheby_Modes.derivs[idv].deriv[sensx->nmod*igv+j] = dmod_dvar[j];


//#define DMOD_DVAR
#ifdef DMOD_DVAR
          for (j = 0; j < sensx->nmod; j++) {
            printf("%lf ", dmod_dvar[j]);
          }
          printf("\n");
#endif

          for (i = 0; i < numFunctional; i++) {
            functional_dvar = values[i]->derivs[3+idv].deriv[igv];

            switch (i+1) {
            case outAlpha:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += sensx->al_mod[j]*dmod_dvar[j];
              }
              break;
            case outCL:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += sensx->cl_mod[j]*dmod_dvar[j];
              }
              break;
            case outCD:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += (sensx->cdv_mod[j] + sensx->cdw_mod[j])*dmod_dvar[j];
              }
              break;
            case outCD_p:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += (sensx->cdv_mod[j] + sensx->cdw_mod[j] - sensx->cdf_mod[i])*dmod_dvar[j];
              }
              break;
            case outCD_v:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += sensx->cdv_mod[j]*dmod_dvar[j];
              }
              break;
            case outCD_w:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += sensx->cdw_mod[j]*dmod_dvar[j];
              }
              break;
            case outCM:
              for (j = 0; j < sensx->nmod; j++) {
                functional_dvar += sensx->cm_mod[j]*dmod_dvar[j];
              }
              break;
            default:
              AIM_ERROR(aimInfo, "Unknown functional %d", i+1);
              status = CAPS_NOTIMPLEMENT;
              goto cleanup;
            }

            values[i]->derivs[3+idv].deriv[igv] = functional_dvar;
          }
        }
      }
    }

    // correct for specifying CL
    if (gcon == GCON_CL) {
      for (idv = 0; idv < msesInstance->design.numDesignVariable; idv++) {

        name = msesInstance->design.designVariable[idv].name;

        // Loop over the geometry in values and compute sensitivities for all bodies
        index = aim_getIndex(aimInfo, name, GEOMETRYIN);
        status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
        AIM_STATUS(aimInfo, status);

        for (irow = 0; irow < geomInVal->nrow; irow++) {
          for (icol = 0; icol < geomInVal->ncol; icol++) {
            /*@-nullderef@*/
            msesInstance->Alpha.derivs[3+idv].deriv[geomInVal->ncol*irow + icol] = -msesInstance->CL.derivs[3+idv].deriv[geomInVal->ncol*irow + icol]/cl_alfa;
            /*@+nullderef@*/
          }
        }
      }
    }

  } // aimInputs[inDesign_Variable-1].nullVal == NotNull
  else if (aimInputs[inCheby_Modes-1].nullVal == NotNull) {

    for (i = 0; i < numFunctional; i++) {
      AIM_STRDUP(values[i]->derivs[3].name, "Cheby_Modes", aimInfo, status);

      AIM_ALLOC(values[i]->derivs[3].deriv, sensx->nmod, double, aimInfo, status);
      values[i]->derivs[3].len_wrt = sensx->nmod;

      switch (i+1) {
      case outAlpha:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = sensx->al_mod[j];
        }
        break;
      case outCL:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = sensx->cl_mod[j];
        }
        break;
      case outCD:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = (sensx->cdv_mod[j] + sensx->cdw_mod[j]);
        }
        break;
      case outCD_p:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = (sensx->cdv_mod[j] + sensx->cdw_mod[j] - sensx->cdf_mod[i]);
        }
        break;
      case outCD_v:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = sensx->cdv_mod[j];
        }
        break;
      case outCD_w:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = sensx->cdw_mod[j];
        }
        break;
      case outCM:
        for (j = 0; j < sensx->nmod; j++) {
          values[i]->derivs[3].deriv[j] = sensx->cm_mod[j];
        }
        break;
      default:
        AIM_ERROR(aimInfo, "Unknown functional %d", i+1);
        status = CAPS_NOTIMPLEMENT;
        goto cleanup;
      }
    }

    // correct for specifying CL
    if (gcon == GCON_CL) {
      for (j = 0; j < sensx->nmod; j++) {
        /*@-nullderef@*/
        msesInstance->Alpha.derivs[3].deriv[j] = -msesInstance->CL.derivs[3].deriv[j]/cl_alfa;
        /*@+nullderef@*/
      }
    }
  }

  status = CAPS_SUCCESS;

cleanup:
  msesSensxFree(&sensx);
  AIM_FREE(values);
  AIM_FREE(xyz);
  AIM_FREE(dx_dvar);
  AIM_FREE(dy_dvar);
  AIM_FREE(M);
  AIM_FREE(rhs);
  AIM_FREE(dmod_dvar);

  if (dxyz != NULL) {
    for (ibody = 0; ibody < numBody; ibody++) {
      AIM_FREE(dxyz[ibody]);
    }
    AIM_FREE(dxyz);
  }

  if (curves != NULL) {
    for (ibody = 0; ibody < numBody; ibody++) {
      EG_deleteObject(curves[ibody]);
    }
    AIM_FREE(curves);
  }

  AIM_FREE(vlmSections);

  return status;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsMSES AIM Outputs
     * The following list outlines the xFoil outputs available through the AIM interface.
     */

    int status = CAPS_SUCCESS;
#ifdef DEBUG
    printf(" msesAIM/aimOutputs index = %d!\n", index);
#endif

    form->type    = Double;
    form->dim     = Scalar;
    form->nullVal = IsNull;
    form->vals.real = 0;

    if (index == outAlpha) {
        *aoname = EG_strdup("Alpha");

        /*! \page aimOutputsMSES
         * - <B> Alpha = </B> Angle of attack value(s).
         */

    } else if (index == outCL) {
        *aoname = EG_strdup("CL");

        /*! \page aimOutputsMSES
         * - <B> CL = </B> Coefficient of lift value(s).
         */

    } else if (index == outCD) {
        *aoname = EG_strdup("CD");

        /*! \page aimOutputsMSES
         * - <B> CD = </B> Coefficient of drag value(s).
         */

    }  else if (index == outCD_p) {
        *aoname = EG_strdup("CD_p");

        /*! \page aimOutputsMSES
         * - <B> CD_p = </B> Coefficient of drag value, pressure contribution.
         */

    }  else if (index == outCD_v) {
        *aoname = EG_strdup("CD_v");

        /*! \page aimOutputsMSES
         * - <B> CD_v = </B> Coefficient of drag value(, viscous contribution.
         */

    }  else if (index == outCD_w) {
        *aoname = EG_strdup("CD_w");

        /*! \page aimOutputsMSES
         * - <B> CD_w = </B> Coefficient of drag value, inviscid (wave) drag from shock entropy wake.
         */

    }  else if (index == outCM) {
        *aoname = EG_strdup("CM");

        /*! \page aimOutputsMSES
         * - <B> CM = </B> Moment coefficient value(s).
         */

    } else if (index == outCheby_Modes) {
        *aoname              = EG_strdup("Cheby_Modes");
        form->type           = Double;
        form->lfixed         = Fixed;
        form->dim            = Vector;
        form->nrow           = NMODE;
        AIM_ALLOC(form->vals.reals, form->nrow, double, aimInfo, status);
        memset(form->vals.reals, 0, form->nrow*sizeof(double));

        /*! \page aimOutputsMSES
         * - <B> Cheby_Modes = </B> <br>
         * Chebyshev shape mode values for shape optimization (fixed length 40).
         */

    }

    status = CAPS_SUCCESS;
cleanup:

    return status;
}


int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/void *aimInfo,
                  int index, capsValue *val)
{
    int status; // Function return status

    aimStorage *msesInstance;

    msesInstance = (aimStorage*)instStore;

    if (index == outAlpha) {
      *val = msesInstance->Alpha;
      aim_initValue(&msesInstance->Alpha);
    } else if (index == outCL) {
      *val = msesInstance->CL;
      aim_initValue(&msesInstance->CL);
    } else if (index == outCD) {
      *val = msesInstance->CD;
      aim_initValue(&msesInstance->CD);
    } else if (index == outCD_p) {
      *val = msesInstance->CDp;
      aim_initValue(&msesInstance->CDp);
    } else if (index == outCD_v) {
      *val = msesInstance->CDv;
      aim_initValue(&msesInstance->CDv);
    } else if (index == outCD_w) {
      *val = msesInstance->CDw;
      aim_initValue(&msesInstance->CDw);
    } else if (index == outCM) {
      *val = msesInstance->CM;
      aim_initValue(&msesInstance->CM);
    } else if (index == outCheby_Modes) {
      *val = msesInstance->Cheby_Modes;
      aim_initValue(&msesInstance->Cheby_Modes);
    }

    status = CAPS_SUCCESS;

//cleanup:

    return status;
}


void aimCleanup(/*@unused@*/ void *instStore)
{
#ifdef DEBUG
  printf(" msesAIM/aimCleanup!\n");
#endif
  aimStorage *msesInstance;

  msesInstance = (aimStorage*)instStore;
  destroy_aimStorage(msesInstance);
  AIM_FREE(instStore);
}
