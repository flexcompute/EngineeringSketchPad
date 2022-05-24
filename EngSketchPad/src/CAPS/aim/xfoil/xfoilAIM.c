/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             XFOIL AIM
 *
 *      Written by Ryan Durscher AFLR/RQVC
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewXFOIL xFoil AIM Overview
 *
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (through input
 * files) with the subsonic airfoil analysis tool xFoil \cite Drela1989. xFoil is an open-source tool and
 * may be freely downloaded from http://web.mit.edu/drela/Public/web/xfoil/ . At this time only a subsection
 * of xFoil's capabilities are exposed through the AIM. Furthermore, only version 6.99 of xFoil
 * have been tested against.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsXFOIL and \ref aimOutputsXFOIL, respectively.
 *
 * Upon running preAnalysis the AIM generates two files: 1. "xfoilInput.txt" which contains instructions for
 * xFoil to execute and 2. "caps.xfoil" which contains the geometry to be analyzed.
 *
 * The xFoil AIM can automatically execute xfoil, with details provided in \ref aimExecuteXFOIL.
 *
 * \section assumptionsXFOIL Assumptions
 * xFoil inherently assumes the airfoil cross-section is in the x-y plane, if it isn't an attempt is made
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
 * was left open at the TE. All of this information will be used to automatically fill in the xFoil geometry
 * description.
 *
 * It should be noted that general construction in either <b> OpenCSM</b> or even <b> EGADS</b> will be supported
 * as long as the topology described above is used. But care should be taken when constructing the airfoil shape
 * so that a discontinuity (i.e.,  simply <em>C<sup>0</sup></em>) is not generated at the <em>Node</em> representing
 * the <em>Leading Edge</em>. This can be done by splining the entire shape as one and then intersecting the single
 *  <em>Edge</em> to place the LE <em>Node</em>.
 *
 * \section exampleXFOIL Examples
 * An example problem using the xFoil AIM may be found at \ref xfoilExample .
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "aimUtil.h"

#include "miscUtils.h" // Bring in miscellaneous utilities
#include "meshUtils.h"
#include "vlmUtils.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#define strtok_r   strtok_s
#endif

#define NUMPOINT  200

#define MXCHAR  255

enum aimInputs
{
  inMach = 1,                      /* index is 1-based */
  inRe,
  inAlpha,
  inAlpha_Increment,
  inCL,
  inCL_Increment,
  inCL_Inviscid,
  inAppend_PolarFile,
  inViscous_Iteration,
  inWrite_Cp,
  NUMINPUT = inWrite_Cp            /* Total number of inputs */
};

enum aimOutputs
{
  outAlpha = 1,                    /* index is 1-based */
  outCL,
  outCD,
  outCD_p,
  outCM,
  outCp_Min,
  outTransition_Top,
  outTransition_Bottom,
  NUMOUTPUT = outTransition_Bottom /* Total number of outputs */
};


/* ********************** Exposed AIM Functions ***************************** */
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{

#ifdef DEBUG
    printf("\n xfoilAIM/aimInitialize   inst = %d!\n", inst);
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

    return CAPS_SUCCESS;
}


int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
    int status = CAPS_SUCCESS;
    /*! \page aimInputsXFOIL AIM Inputs
     * The following list outlines the xFoil inputs along with their default values available
     * through the AIM interface.
     */

#ifdef DEBUG
    printf(" xfoilAIM/aimInputs  index = %d!\n", index);
#endif

    if (index == inMach) {
        *ainame           = EG_strdup("Mach");
        defval->type      = Double;
        defval->vals.real = 0.0;

        /*! \page aimInputsXFOIL
         * - <B> Mach = 0.0 </B> <br>
         *  Mach number.
         */

    } else if (index == inRe) {
        *ainame              = EG_strdup("Re");
        defval->type      = Double;
        defval->vals.real = 0.0;

        /*! \page aimInputsXFOIL
         * - <B> Re = 0.0 </B> <br>
         *  Reynolds number.
         */
    } else if (index == inAlpha) {
        *ainame           = EG_strdup("Alpha");
        defval->type      = Double;
        defval->dim       = Vector;
        defval->nrow      = 1;
        defval->ncol      = 1;
        defval->units     = NULL;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        defval->lfixed    = Change;
        //defval->units     = EG_strdup("degree");

        /*! \page aimInputsXFOIL
         * - <B> Alpha = NULL </B> <br>
         *  Angle of attack [degree], either a single value or an array of values ( [0.0, 4.0, ...] ) may be provided.
         */
    } else if (index == inAlpha_Increment) {
        *ainame               = EG_strdup("Alpha_Increment");
        defval->type          = Double;
        defval->dim           = Vector;
        defval->nrow          = 3;
        AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
        defval->vals.reals[0] = 0.0;
        defval->vals.reals[1] = 0.0;
        defval->vals.reals[2] = 0.0;
        defval->nullVal       = IsNull;
        defval->lfixed        = Fixed;
        //defval->units      = EG_strdup("degree");

        /*! \page aimInputsXFOIL
         * - <B> Alpha_Increment = NULL</B> <br>
         *  Angle of attack [degree] sequence - [first value, last value, increment].
         */
    }else if (index == inCL) {
        *ainame           = EG_strdup("CL");
        defval->type      = Double;
        defval->dim       = Vector;
        defval->nrow      = 1;
        defval->ncol      = 1;
        defval->vals.real = 0.0;
        defval->lfixed    = Change;
        defval->nullVal   = IsNull;

        /*! \page aimInputsXFOIL
         * - <B> CL =  NULL</B> <br>
         *  Prescribed coefficient of lift, either a single value or an array of values ( [0.1, 0.5, ...] ) may be provided.
         */
    } else if (index == inCL_Increment) {
        *ainame               = EG_strdup("CL_Increment");
        defval->type          = Double;
        defval->dim           = Vector;
        defval->nrow          = 3;
        defval->ncol          = 1;
        defval->units         = NULL;
        AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
        defval->vals.reals[0] = 0.0;
        defval->vals.reals[1] = 0.0;
        defval->vals.reals[2] = 0.0;
        defval->nullVal       = IsNull;
        defval->lfixed        = Fixed;

        /*! \page aimInputsXFOIL
         * - <B> CL_Increment = NULL</B> <br>
         *  Prescribed coefficient of lift sequence - [first value, last value, increment].
         */
    } else if (index == inCL_Inviscid) {
        *ainame           = EG_strdup("CL_Inviscid");
        defval->type      = Double;
        defval->dim       = Vector;
        defval->nrow      = 1;
        defval->ncol      = 1;
        defval->vals.real = 0.0;
        defval->lfixed    = Change;
        defval->nullVal   = IsNull;

        /*! \page aimInputsXFOIL
         * - <B> CL_Inviscid =  NULL</B> <br>
         *  Prescribed inviscid coefficient of lift, either a single value or an array of values ( [0.1, 0.5, ...] ) may be provided.
         */
    } else if (index == inAppend_PolarFile) {
        *ainame              = EG_strdup("Append_PolarFile");
        defval->type         = Boolean;
        defval->vals.integer = (int) false;

        /*! \page aimInputsXFOIL
         * - <B> Append_PolarFile =  False</B> <br>
         *  Append the file (xfoilPolar.dat) that polar data is written to.
         */
    } else if (index == inViscous_Iteration) {
        *ainame              = EG_strdup("Viscous_Iteration");
        defval->type         = Integer;
        defval->vals.integer = 20;

        /*! \page aimInputsXFOIL
         * - <B> Viscous_Iteration = 20</B> <br>
         *  Viscous solution iteration limit. Only set if a Re isn't 0.0 .
         */
    } else if (index == inWrite_Cp) {
        *ainame              = EG_strdup("Write_Cp");
        defval->type         = Boolean;
        defval->vals.integer = (int) false;

        /*! \page aimInputsXFOIL
         * - <B> Write_Cp =  False</B> <br>
         *  Have xFoil write the coefficient of pressure (Cp) distribution to a file - xfoilCp.dat.
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

    int i, bodyIndex; // Indexing

    // Bodies
    const char *intents;
    int numBody;
    ego *bodies = NULL;

    vlmSectionStruct vlmSection;

    // File I/O
    FILE *fp = NULL;
    char inputFilename[] = "xfoilInput.txt";
    char xfoilFilename[] = "caps.xfoil";
    char polarFilename[] = "xfoilPolar.dat";
    char cpFilename[]    = "xfoilCp.dat";

    if (aimInputs == NULL) {
#ifdef DEBUG
        printf("\txfoilAIM/aimPreAnalysis aimInputs == NULL!\n");
#endif
        return CAPS_NULLVALUE;
    }

    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody == 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }

    if (numBody != 1) {
        AIM_ERROR(aimInfo, "Only one body should be provided to the xfoilAIM! numBody = %d", numBody);
        return CAPS_SOURCEERR;
    }

    status = initiate_vlmSectionStruct(&vlmSection);
    AIM_STATUS(aimInfo, status);

    // Accumulate cross coordinates of airfoil and write out data file
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        // Open and write the input to control the XFOIL session
        fp = aim_fopen(aimInfo, xfoilFilename, "w");
        if (fp == NULL) {
            AIM_ERROR(aimInfo, "Unable to open file %s\n!", xfoilFilename);
            status = CAPS_IOERR;
            goto cleanup;
        }

        fprintf(fp,"capsBody_%d\n",bodyIndex+1);

        status = EG_copyObject(bodies[bodyIndex], NULL, &vlmSection.ebody);
        AIM_STATUS(aimInfo, status);

        status = finalize_vlmSectionStruct(aimInfo, &vlmSection);
        AIM_STATUS(aimInfo, status);

        // Write out the airfoil cross-section given an ego body
        status = vlm_writeSection(aimInfo,
                                  fp,
                                  &vlmSection,
                                  (int) false, // Normalize by chord (true/false)
                                  (int) NUMPOINT);
        AIM_STATUS(aimInfo, status);

        status = destroy_vlmSectionStruct(&vlmSection);
        AIM_STATUS(aimInfo, status);

        // Close file
        if (fp != NULL) {
            fclose(fp);
            fp = NULL;
        }

    } // End body loop

    // Open and write the input to control the XFOIL session
    if (fp != NULL) fclose(fp);
    fp = aim_fopen(aimInfo, inputFilename, "w");
    if (fp == NULL) {
        printf("Unable to open file %s\n!", inputFilename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Print the session file for the XFOIL run
    fprintf(fp, "PLOP\n"); // Enter PLOP
    fprintf(fp, "G F\n"); // Turn off graphics
    fprintf(fp, "\n"); // Exit PLOP

    fprintf(fp, "LOAD\n"); // Load airfoil
    fprintf(fp, "%s\n", xfoilFilename);

    fprintf(fp, "PANE\n"); // USE PANE Option to clean up airfoil mesh

    fprintf(fp, "OPER\n"); // Enter OPER

    fprintf(fp, "Mach %f\n", aimInputs[inMach-1].vals.real);    // Set Mach number

    if (aimInputs[inRe-1].vals.real > 0) {
        fprintf(fp, "Viscr\n");
        fprintf(fp,	"%f\n", aimInputs[inRe-1].vals.real); // Set Reynolds number

        fprintf(fp, "ITER\n");
        fprintf(fp,	"%d\n", aimInputs[inViscous_Iteration-1].vals.integer); // Set iteration limit
    }

    fprintf(fp,"CINC\n"); // Turn of minimum Cp inclusion in polar data
    fprintf(fp,"PACC\n"); // Turn of polar data accmulation

    fprintf(fp,"%s\n", polarFilename);
    fprintf(fp,"\n");

    // Check to see if polar file exist
    if (aim_isFile(aimInfo, polarFilename) == CAPS_SUCCESS){
        if (aimInputs[inAppend_PolarFile-1].vals.integer == (int) false) {
            status = aim_rmFile(aimInfo, polarFilename);
            AIM_STATUS(aimInfo, status);
        } else {
            fprintf(fp,"n\n");
        }
    }

    //fprintf(fp,"xfoilPolar.dump\n");

    // Alpha
    if (aimInputs[inAlpha-1].nullVal ==  NotNull) {
        if (aimInputs[inAlpha-1].length == 1) {
            fprintf(fp, "Alfa %f\n", aimInputs[inAlpha-1].vals.real);
        } else {
            for (i = 0; i < aimInputs[inAlpha-1].length; i++) {
                fprintf(fp, "Alfa %f\n", aimInputs[inAlpha-1].vals.reals[i]);
            }
        }
    }

    // Alpha sequence
    if (aimInputs[inAlpha_Increment-1].nullVal ==  NotNull) {

        fprintf(fp, "ASeq %f %f %f\n",
                aimInputs[inAlpha_Increment-1].vals.reals[0],
                aimInputs[inAlpha_Increment-1].vals.reals[1],
                aimInputs[inAlpha_Increment-1].vals.reals[2]);
    }

    // CL
    if (aimInputs[inCL-1].nullVal ==  NotNull) {

        if (aimInputs[inCL-1].length == 1) {
            fprintf(fp, "CL %f\n", aimInputs[inCL-1].vals.real);
        } else {
            for (i = 0; i < aimInputs[inCL-1].length; i++) {
                fprintf(fp, "CL %f\n", aimInputs[inCL-1].vals.reals[i]);
            }
        }
    }

    // CL sequence
    if (aimInputs[inCL_Increment-1].nullVal ==  NotNull) {

        fprintf(fp, "CSeq %f %f %f\n", aimInputs[inCL_Increment-1].vals.reals[0],
                                       aimInputs[inCL_Increment-1].vals.reals[1],
                                       aimInputs[inCL_Increment-1].vals.reals[2]);
    }

    // CL Invicid
    if (aimInputs[inCL_Inviscid-1].nullVal ==  NotNull) {

        if (aimInputs[inCL_Inviscid-1].length == 1) {
            fprintf(fp, "CLI %f\n", aimInputs[inCL_Inviscid-1].vals.real);
        } else {
            for (i = 0; i < aimInputs[inCL_Inviscid-1].length; i++) {
                fprintf(fp, "CLI %f\n", aimInputs[inCL_Inviscid-1].vals.reals[i]);
            }
        }
    }

    // Write Cp distribution to a file
    if (aimInputs[inWrite_Cp-1].vals.integer == (int) true) {

        fprintf(fp, "CPWR\n");
        fprintf(fp, "%s\n", cpFilename);
    }

    fprintf(fp,"\n"); // Exit OPER

    fprintf(fp, "Quit\n"); // Quit XFOIL

    fclose(fp);
    fp = NULL;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("xfoil/preAnalysis status = %d",
                                           status);
        if (fp != NULL) fclose(fp);

        return status;
}


int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteXFOIL AIM Execution
   *
   * If auto execution is enabled when creating an xFoil AIM,
   * the AIM will execute xFoil just-in-time with the command line:
   *
   * \code{.sh}
   * xfoil < xfoilInput.txt > xfoilOutput.txt
   * \endcode
   *
   * where preAnalysis generated the file "xFoilInput.txt" which contains the input information.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an xFoil AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or xFoil can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * xfoil.preAnalysis()
   *
   * print ("\n\nRunning......")
   * xfoil.system("xfoil < xfoilInput.txt > xfoilOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * xfoil.postAnalysis()
   * \endcode
   */

  *state = 0;
  return aim_system(aimInfo, NULL, "xfoil < xfoilInput.txt > xfoilOutput.txt");
}


int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  // check the XFOIL output file
  if (aim_isFile(aimInfo, "xfoilPolar.dat") != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "xfoil execution did not produce xfoilPolar.dat");
    return CAPS_EXECERR;
  }

  return CAPS_SUCCESS;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{

    /*! \page aimOutputsXFOIL AIM Outputs
     * The following list outlines the xFoil outputs available through the AIM interface.
     */

#ifdef DEBUG
    printf(" xfoilAIM/aimOutputs index = %d!\n", index);
#endif


    if (index == 1) {
        *aoname = EG_strdup("Alpha");

        /*! \page aimOutputsXFOIL
         * - <B> Alpha = </B> Angle of attack value(s).
         */

    } else if (index == 2) {
        *aoname = EG_strdup("CL");

        /*! \page aimOutputsXFOIL
         * - <B> CL = </B> Coefficient of lift value(s).
         */

    } else if (index == 3) {
        *aoname = EG_strdup("CD");

        /*! \page aimOutputsXFOIL
         * - <B> CD = </B> Coefficient of drag value(s).
         */

    }  else if (index == 4) {
        *aoname = EG_strdup("CD_p");

        /*! \page aimOutputsXFOIL
         * - <B> CD_p = </B> Coefficient of drag value(s), pressure contribution.
         */

    }  else if (index == 5) {
        *aoname = EG_strdup("CM");

        /*! \page aimOutputsXFOIL
         * - <B> CM = </B> Moment coefficient value(s).
         */

    } else if (index == 6) {
        *aoname = EG_strdup("Cp_Min");

        /*! \page aimOutputsXFOIL
         * - <B> Cp_Min = </B> Minimum coefficient of pressure value(s).
         */
    } else if (index == 7) {
        *aoname = EG_strdup("Transition_Top");

        /*! \page aimOutputsXFOIL
         * - <B> Transition_Top = </B> Value(s) of x- transition location on the top of the airfoil.
         */
    } else if (index == 8) {
        *aoname = EG_strdup("Transition_Bottom");

        /*! \page aimOutputsXFOIL
         * - <B> Transition_Bottom = </B> Value(s) of x- transition location on the bottom of the airfoil.
         */
    }

    form->type    = Double;
    form->dim     = Vector;
    form->nrow    = 1;
    form->ncol    = 1;
    form->units   = NULL;
    form->lfixed  = Change;
    form->nullVal = IsNull;

    form->vals.reals = NULL;

    return CAPS_SUCCESS;
}


int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/void *aimInfo,
                  int index, capsValue *val)
{
    int status; // Function return status
    int i, j, k; // Indexing
    int skipLine = (int) false;
    int valCount = 0; // Number of values in return array

    int valIndex = 0;

    // excess data entries in case future versions of xfoil change
#define MAX_DATA_ENTRY 20
    int numDataEntry = 0;
    double dataLine[MAX_DATA_ENTRY];

    size_t     linecap = 0;
    char       *line = NULL, *rest = NULL, *token = NULL;
    char       headers[MAX_DATA_ENTRY][20];
    const char *valHeader;
    FILE       *fp;

    // Open the XFOIL output file
    fp = aim_fopen(aimInfo, "xfoilPolar.dat", "r");
    if (fp == NULL) return CAPS_IOERR;

    // Move to beginning of data
    for (i = 0; i < 11; i++) {
        status = getline(&line, &linecap, fp);
        if (status < 0){
            AIM_ERROR(aimInfo, "Could not parse xfoilPolar.dat");
            status = CAPS_IOERR;
            goto cleanup;
        }
    }

    // Parse the header information,
    // both headers can come from xfoil 6.99 (lax version control), i.e.
    //    alpha    CL        CD       CDp       CM      Cpmin   XCpmin   Top_Xtr  Bot_Xtr
    // or
    //    alpha    CL        CD       CDp       CM      Cpmin   Top_Xtr  Bot_Xtr  Top_Itr  Bot_Itr
    rest = line;

    numDataEntry = 0;
    while ((token = strtok_r(rest, " ", &rest))) {
        strcpy(headers[numDataEntry], token);
        //printf("%s\n", headers[numDataEntry]);
        numDataEntry++;
        if (numDataEntry == MAX_DATA_ENTRY) {
            AIM_ERROR(aimInfo, "More than %d columns in xfoilPolar.dat is not expected!",
                      numDataEntry);
            status = CAPS_IOERR;
            goto cleanup;
        }
    }

    // line skip the heading just above the data, i.e. ---- ---- ----...
    status = getline(&line, &linecap, fp);
    if (status < 0){
        AIM_ERROR(aimInfo, "Could not parse xfoilPolar.dat");
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Headers expected in xfoilPolar.dat that correspond to the AIM output names
  
         if (index == outAlpha)
        valHeader = "alpha";
    else if (index == outCL)
        valHeader = "CL";
    else if (index == outCD)
        valHeader = "CD";
    else if (index == outCD_p)
        valHeader = "CDp";
    else if (index == outCM)
        valHeader = "CM";
    else if (index == outCp_Min)
        valHeader = "Cpmin";
    else if (index == outTransition_Top)
        valHeader = "Top_Xtr";
    else if (index == outTransition_Bottom)
        valHeader = "Bot_Xtr";
    else {
        AIM_ERROR(aimInfo, "Developer error: Unkown variable index %d", index);
        status = CAPS_BADINDEX;
        goto cleanup;
    }

    // Find which column contains the requested data
    valIndex = 0;
    while( valIndex < numDataEntry && strncmp(headers[valIndex], valHeader,
                                              strlen(valHeader)) != 0 ) valIndex++;
    if (valIndex == numDataEntry) {
        AIM_ERROR(aimInfo, "Could not find '%s' header in xfoilPolar.dat", valHeader);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Count how many line of data we have
    valCount = 0;
    while (getline(&line, &linecap, fp) >= 0) valCount += 1;

    // Reset file
    rewind(fp);

    if (valCount == 0) {
        AIM_ERROR(aimInfo, "No data in xfoilPolar.dat");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    if (valCount > 1) {
        AIM_ALLOC(val->vals.reals, valCount, double, aimInfo, status);
    }

    // Move back to the beginning of data
    for (i = 0; i < 12; i++) {
        (void) getline(&line, &linecap, fp);
    }

    k = 0;
    for (j = 0; j < valCount; j++) {
        skipLine = (int) false;

        for (i = 0; i < numDataEntry; i++) {
            status = fscanf(fp, "%lf", &dataLine[i]);

            if (status == CAPS_SUCCESS) {
                printf("Invalid data line, %d, skipping data entry!\n", j+1);
                (void) getline(&line, &linecap, fp);

                skipLine = (int) true;
                break;
            }
        }

        if (skipLine == (int) true) continue;

        if (valCount == 1) {
            val->vals.real = dataLine[valIndex];
        } else {
            //printf("k = %d, Val (%d) - %f\n", k, valIndex,dataLine[valIndex]);
            val->vals.reals[k] = dataLine[valIndex];
        }
        k += 1;
    }

    if (k != valCount) { // Realloc our array - shorten it
        AIM_REALL(val->vals.reals, valCount, double, aimInfo, status);
    }

    valCount = k;
    val->dim     = Vector;
    val->nrow    = valCount;
    val->ncol    = 1;
    val->nullVal = NotNull;

    status = CAPS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);

    AIM_FREE(line);

    return status;
}


void aimCleanup(/*@unused@*/ void *instStore)
{

#ifdef DEBUG
    printf(" xfoilAIM/aimCleanup!\n");
#endif

}
