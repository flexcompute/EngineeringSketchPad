/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             TSFOIL AIM
 *
 *      Written by Ryan Durscher AFLR/RQVC
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewTSFOIL TSFOIL AIM Overview
 *
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (through input
 * files) with the transonic airfoil analysis tool TSFOIL. TSFOIL can be downloaded from
 * http://www.dept.aoe.vt.edu/~mason/Mason_f/MRsoft.html .
 *
 * Note: In the tsfoil2.f file is may be necessary to comment out line 38 - "USE DFPORT"
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsTSFOIL and \ref aimOutputsTSFOIL, respectively.
 *
 * Upon running preAnalysis the AIM generates two files: 1. "tsfoilInput.txt" which contains instructions for
 * TSFOIL to execute and 2. "caps.tsfoil" which contains the geometry to be analyzed.
 * The TSFOIL AIM can automatically execute TSFOIL, with details provided in \ref aimExecuteTSFOIL.
 *
 * \section assumptionsTSFOIL Assumptions
 * TSFOIL inherently assumes the airfoil cross-section is in the x-y plane, if it isn't an attempt is made
 * to automatically rotate the provided body.
 *
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
#endif

#define MXCHAR  255


enum aimInputs
{
  inMach = 1,                  /* index is 1-based */
  inRe,
  inAlpha,
  NUMINPUT = inAlpha           /* Total number of inputs */
};

enum aimOutputs
{
  outCL = 1,                   /* index is 1-based */
  outCD,
  outCD_Wave,
  outCM,
  outCp_Critical,
  NUMOUTPUT = outCp_Critical   /* Total number of outputs */
};


/* ********************** Exposed AIM Functions ***************************** */
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{

#ifdef DEBUG
    printf("\n tsfoilAIM/aimInitialize   inst = %d!\n", inst);
#endif

    // Specify the number of analysis input and out "parameters"
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


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
    /*! \page aimInputsTSFOIL AIM Inputs
     * The following list outlines the TSFOIL inputs along with their default values available
     * through the AIM interface.
     */

#ifdef DEBUG
    printf(" tsfoilAIM/aimInputs  index = %d!\n", index);
#endif

    if (index == inMach) {
        *ainame           = EG_strdup("Mach");
        defval->type      = Double;
        defval->vals.real = 0.75;
        defval->units     = NULL;
        defval->limits.dlims[0] = 0.5; // Limits on Mach for TSFOIL
        defval->limits.dlims[1] = 2.0;

        /*! \page aimInputsTSFOIL
         * - <B> Mach = 0.75 </B> <br>
         *  Mach number. Valid range for TSFOIL is 0.5 to 2.0 .
         */

    } else if (index == inRe) {
        *ainame              = EG_strdup("Re");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsTSFOIL
         * - <B> Re = 0.0 </B> <br>
         *  Reynolds number based on chord length.
         */
    } else if (index == inAlpha) {
        *ainame           = EG_strdup("Alpha");
        defval->type      = Double;
        defval->lfixed    = Fixed;
        //defval->units     = EG_strdup("degree");

        /*! \page aimInputsTSFOIL
         * - <B> Alpha = 0.0 </B> <br>
         *  Angle of attack [degree].
         */
    }


    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(/*@unused@*/ void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status

    int i, bodyIndex, loopIndex, edgeIndex; // Indexing

    // Bodies
    const char *intents;
    int numBody;
    ego *bodies = NULL;

    // Topology
    int numLoop, numEdge;
    ego *loops = NULL, *edges = NULL, geom;
    int oclass, mtype;
    //double uvbox[4];
    int *edgeSense = NULL;

    // Attribute to index map
    //mapAttrToIndexStruct attrMap;

    // Tessellation
    double size, params[3], paramsScaled[3];
    int numEdgePoint = 40;
    double rPos[40]; //numEdgePoint
    ego egadsTess;
    int numPoints;
    const double *points = NULL, *uv = NULL;
    double  box[6];

    // Proper orientation
    int xMeshConstant = (int) true, yMeshConstant = (int) true, zMeshConstant= (int) true;
    int swapZX = (int) false, swapZY = (int) false;

    // File I/O
    FILE *fp = NULL;
    char inputFilename[] = "tsfoilInput.txt";
    char outputFilename[] = "tsfoilOutput.txt";
    char tsfoilFilename[] = "caps.tsfoil";

    if (aimInputs == NULL) {
#ifdef DEBUG
        printf("\ttsfoilAIM/aimPreAnalysis aimInputs == NULL!\n");
#endif
        return CAPS_NULLVALUE;
    }

    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    if (status != CAPS_SUCCESS) {
#ifdef DEBUG
        printf("\ttsfoilAIM/aimPreAnalysis getBodies = %d!\n", status);
#endif
        return status;
    }

    if (numBody == 0 || bodies == NULL) {
        printf("\tError: tsfoilAIM/aimPreAnalysis No Bodies!\n");
        return CAPS_SOURCEERR;
    }

    if (numBody != 1) {
        printf("\tError: Only one body should be provided to the tsfoilAIM at this time!!");
        return CAPS_SOURCEERR;
    }

    // Setup even distribution for grid points along edge
    for (i = 0; i < numEdgePoint; i++) {
        rPos[i] = (double) (i+1) / (double)(numEdgePoint-1);
    }

    //  Loop through bodies and set spacing
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        // Get edges on the body - there should only be 2 - checked later
        status = EG_getBodyTopos(bodies[bodyIndex], NULL, EDGE, &numEdge, &edges);
        if (status != EGADS_SUCCESS) goto cleanup;

        if ((numEdge != 2) || (edges == NULL)) {
            if (edges != NULL) EG_free(edges);
            edges = NULL;
            printf("\tError: The airfoil cross-section (of body %d) should consist of two edges!\n",
                   bodyIndex+1);
            goto cleanup;
        }

        for (i = 0; i < numEdge; i++) {
            status = EG_attributeAdd(edges[i], ".rPos", ATTRREAL,
                                     numEdgePoint-2, NULL, rPos, NULL);
        }

        if (edges != NULL) EG_free(edges);
        edges = NULL;

        if (status != EGADS_SUCCESS) goto cleanup;
    }

    // Open and write the input to control the TSFOIL session
    fp = aim_fopen(aimInfo, tsfoilFilename, "w");
    if (fp == NULL) {
        printf("\tUnable to open file %s\n!", tsfoilFilename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp,"CAPS generated TSFOIL input\n");
    fprintf(fp,"$INP\n"); // Start input namelist

    fprintf(fp, "EMACH=%f\n", aimInputs[inMach-1].vals.real);

    fprintf(fp, "DELTA=0.115\n");

    fprintf(fp, "ALPHA=%f\n", aimInputs[inAlpha-1].vals.real);

    fprintf(fp, "GAM=1.4\n");

    fprintf(fp, "WE=1.8,1.9,1.95\n");
    fprintf(fp, "EPS=0.2\n");

    fprintf(fp, "AMESH=T\n");

    fprintf(fp, "RIGF=0.0\n");
    fprintf(fp, "CVERGE=0.00001\n");
    fprintf(fp, "BCFOIL=4\n");
    fprintf(fp, "MAXIT=800\n");

    if (aimInputs[inRe-1].vals.real > 0) {
        fprintf(fp, "REYNLD=%e\n", aimInputs[inRe-1].vals.real);
    }

    fprintf(fp, "$END\n"); // End input namelist

    params[0] = 0.025;
    params[1] = 0.001;
    params[2] = 15.00;

    // Accumulate cross coordinates of airfoil and write out data file
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        // Reset orientation information for each body
        xMeshConstant = (int) true;
        yMeshConstant = (int) true;
        zMeshConstant= (int) true;
        swapZX = (int) false;
        swapZY = (int) false;

        fprintf(fp,"capsBody_%d\n",bodyIndex+1);
        fprintf(fp,"%10.5f%10.5f%10.5f\n", 0.0, (float) numEdgePoint,
                (float) numEdgePoint);

        // Check for x-y plane data
        status = EG_getBoundingBox(bodies[bodyIndex], box);
        if (status != EGADS_SUCCESS) goto cleanup;

        // Constant x?
        if (box[3] == box[0]) xMeshConstant = (int) false; // This is not that robust - may cause issues in the future

        // Constant y?
        if (box[4] == box[1]) yMeshConstant = (int) false;

        // Constant z?
        if (box[5] == box[2]) zMeshConstant = (int) false;

        //for (i = 0; i < 6; i++) printf("Box[%d] = %f\n", i, box[i]);

        if (zMeshConstant != (int) true) {
            printf("TSFOIL expects airfoil cross sections to be in the x-y plane... attempting to rotate body %d!\n",
                   bodyIndex);

            if (xMeshConstant == (int) true && yMeshConstant == (int) false) {
                printf("\tSwapping z and x coordinates!\n");
                swapZX = (int) true;

            } else if(xMeshConstant == (int) false && yMeshConstant == (int) true) {
                printf("\tSwapping z and y coordinates!\n");
                swapZY = (int) true;

            } else {
                printf("\tUnable to rotate mesh!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        // Get tessellation size
        size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                    (box[1]-box[4])*(box[1]-box[4]) +
                    (box[2]-box[5])*(box[2]-box[5]));

        paramsScaled[0] = params[0]*size;
        paramsScaled[1] = params[1]*size;
        paramsScaled[2] = params[2];

        // Make tessellation
        status = EG_makeTessBody(bodies[bodyIndex], paramsScaled, &egadsTess);
        if (status != EGADS_SUCCESS) {
            printf("\tProblem during tessellation of body %d\n", bodyIndex+1);
            goto cleanup;
        }

        // Get loops on the body - there should only be 1
        status = EG_getBodyTopos(bodies[bodyIndex], NULL, LOOP, &numLoop, &loops);
        if (status != EGADS_SUCCESS) goto cleanup;

        if ((numLoop != 1) || (loops == NULL)) {
            printf("\tError: The number of loops on body %d is more than 1!\n",
                   bodyIndex+1);
            goto cleanup;
        }

        for (loopIndex = 0; loopIndex < numLoop; loopIndex++) {

            // Get edges in loop
            status = EG_getTopology(loops[loopIndex], &geom, &oclass, &mtype,
                                    NULL, &numEdge, &edges, &edgeSense);
            if (status != EGADS_SUCCESS) {
                printf("\tEG_getTopology status = %d\n",status);
                goto cleanup;
            }

            if ((numEdge != 2) || (edges == NULL)) {
                printf("\tError: The airfoil cross-section (of body %d) should consist of two edges!\n",
                       bodyIndex+1);
                goto cleanup;
            }

            // Loop through edges on the loop
            for (edgeIndex = 0; edgeIndex < numEdge; edgeIndex++) {

                fprintf(fp,"Edge_%d\n", edgeIndex+1);

                // The particular body relative edge index
                //printf("edgeIndex = %d\n", edgeIndex);
                i = EG_indexBodyTopo(bodies[bodyIndex], edges[edgeIndex]);
                if (i < EGADS_SUCCESS) {
                    status = CAPS_BADINDEX;
                    goto cleanup;
                }

                // Get the edge tessellation
                //printf("EdgeBodyIndex = %d %d\n", i, edgeSense[edgeIndex]);
                status = EG_getTessEdge(egadsTess, i, &numPoints, &points, &uv);
                if ((status != EGADS_SUCCESS) || (points == NULL)) {
                    printf("\tEG_getTessEdge status = %d\n",status);
                    goto cleanup;
                }

                if (edgeIndex == 1) { // if (edgeSense[edgeIndex] > 0) {
                    for (i = 0; i < numPoints; i++) {

                        if (swapZX == (int) true) {
                            // x = z
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 2],
                                    points[3*i + 1]);

                        } else if (swapZY == (int) true) {
                            // y = z
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 0],
                                    points[3*i + 2]);
                        } else {
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 0],
                                    points[3*i + 1]);
                        }
                    }

                } else {

                    for (i = numPoints-1; i >= 0; i--) {

                        if (swapZX == (int) true) {
                            // x = z
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 2],
                                    points[3*i + 1]);

                        } else if (swapZY == (int) true) {
                            // y = z
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 0],
                                    points[3*i + 2]);
                        } else {
                            fprintf(fp, "%10.5f%10.5f\n", points[3*i + 0],
                                    points[3*i + 1]);
                        }
                    }
                } // End sense loop
            } // End edge loop
        } // End loop loop

        // Close file
        if (fp != NULL) {
            fclose(fp);
            fp = NULL;
        }

        // Free loop
        if (loops != NULL) EG_free(loops);
        loops = NULL;

    } // End body loop

    // Open and write the input to control the TSFOIL session
    if (fp != NULL) fclose(fp);
    fp = aim_fopen(aimInfo, inputFilename, "w");
    if (fp == NULL) {
        printf("Unable to open file %s\n!", inputFilename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Print the session file for the TSFOIL run
    fprintf(fp, "default\n");            // Set directory
    fprintf(fp, "%s\n", outputFilename); // Set outputFilename
    fprintf(fp, "%s\n", tsfoilFilename); // Set inputfile
    fprintf(fp, "\n"); // Exit
    fprintf(fp, "\n"); // 3 times
    fprintf(fp, "\n"); // just to be safe

    fclose(fp);
    fp = NULL;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS)
            printf("tsfoil/preAnalysis status = %d", status);
        if (fp != NULL) fclose(fp);

        if (loops != NULL) EG_free(loops);

        return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteTSFOIL AIM Execution
   *
   * If auto execution is enabled when creating an TSFOIL AIM,
   * the AIM will execute TSFOIL just-in-time with the command line:
   *
   * \code{.sh}
   * tsfoil2 < tsfoilInput.txt > Info.out
   * \endcode
   *
   * where preAnalysis generated the file "frictionInput.txt" which contains the input information.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an TSFOIL AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or TSFOIL can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * tsfoil.preAnalysis()
   *
   * print ("\n\nRunning......")
   * tsfoil.system("tsfoil2 < tsfoilInput.txt > Info.out"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * tsfoil.postAnalysis()
   * \endcode
   */

  *state = 0;
  return aim_system(aimInfo, NULL, "tsfoil2 < tsfoilInput.txt > Info.out");
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  // check the friction output file
  if (aim_isFile(aimInfo, "tsfoilOutput.txt") != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "tsfoil2 execution did not produce tsfoilOutput.txt");
    return CAPS_EXECERR;
  }

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int index, char **aoname, capsValue *form)
{

    /*! \page aimOutputsTSFOIL AIM Outputs
     * The following list outlines the TSFOIL outputs available through the AIM interface.
     */

#ifdef DEBUG
    printf(" tsfoilAIM/aimOutputs  index = %d!\n", index);
#endif


    if (index == outCL) {
        *aoname = EG_strdup("CL");

        /*! \page aimOutputsTSFOIL
         * - <B> CL = </B> Coefficient of lift value.
         */

    } else if (index == outCD) {
        *aoname = EG_strdup("CD");

        /*! \page aimOutputsTSFOIL
         * - <B> CD = </B> Coefficient of drag value. (Calculated from momentum integral)
         */

    } else if (index == outCD_Wave) {
        *aoname = EG_strdup("CD_Wave");

        /*! \page aimOutputsTSFOIL
         * - <B> CD_Wave = </B> Wave drag coefficient value.
         */

    }  else if (index == outCM) {
        *aoname = EG_strdup("CM");

        /*! \page aimOutputsTSFOIL
         * - <B> CM = </B> Moment coefficient value.
         */

    } else if (index == outCp_Critical) {
        *aoname = EG_strdup("Cp_Critical");

        /*! \page aimOutputsTSFOIL
         * - <B> Cp_Critical = </B> Critical pressure coefficient (M = 1).
         */
    }

    form->type   = Double;
    form->dim    = Vector;
    form->nrow   = 1;
    form->ncol   = 1;
    form->units  = NULL;
    form->lfixed = Fixed;

    form->vals.reals = NULL;
    form->vals.real = 0;

    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                  int index, capsValue *val)
{
    int status; // Function return status

    size_t     linecap = 0;
    char       *line = NULL;
    FILE       *fp;

    //char beginData[] = "                     FINAL MESH";

    char *valstr = NULL;

    // open the TSFOIL output file
    fp = aim_fopen(aimInfo, "tsfoilOutput.txt", "r");
    if (fp == NULL) return CAPS_DIRERR;

    status = CAPS_NOTFOUND;
    while (getline(&line, &linecap, fp) >= 0) {

        valstr = NULL;

        if (index == outCL || index == outCM || index == outCp_Critical) {

            if (line == NULL) {
                status = CAPS_NOTFOUND;
                goto cleanup;
            }
            valstr = strstr(line, "FINAL MESH"); // Look for final mesh values
            if (valstr != NULL) {

                status = getline(&line, &linecap, fp);
                if (status < 0) {
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }


                if (index == outCL) { // Get CL variable

                    valstr = strstr(line, "CL =");

                    sscanf(&valstr[5], "%lf", &val->vals.real);
                    status = CAPS_SUCCESS;
                    break;
                }

                status = getline(&line, &linecap, fp); // Skip to next line
                if (status < 0) {
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

                if (index == outCM) { // Get CM variable

                    valstr = strstr(line, "CM =");

                    sscanf(&valstr[5], "%lf", &val->vals.real);
                    status = CAPS_SUCCESS;
                    break;
                }

                status = getline(&line, &linecap, fp); // Skip to next line
                if (status < 0) {
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

                if (index == outCp_Critical) { // Get Cp critical

                    valstr = strstr(line, "CP* =");
                    sscanf(&valstr[5], "%lf", &val->vals.real);
                    status = CAPS_SUCCESS;
                    break;
                }

            }
        }

        if (index == outCD_Wave) { // Get wave drag
            if (line == NULL) {
                status = CAPS_NOTFOUND;
                goto cleanup;
            }
            valstr = strstr(line, "TOTAL CDWAVE =");

            if (valstr != NULL) {
                sscanf(&valstr[15], "%lf", &val->vals.real);
                status = CAPS_SUCCESS;
                break;
            }
        }

        if (index == outCD) { //Get drag value
            if (line == NULL) {
                status = CAPS_NOTFOUND;
                goto cleanup;
            }
            valstr = strstr(line, "CD     =");

            if (valstr != NULL) {
                sscanf(&valstr[8], "%lf", &val->vals.real);
                status = CAPS_SUCCESS;
                break;
            }
        }

    }

    goto cleanup;

    cleanup:
        if (fp   != NULL) fclose(fp);
        if (line != NULL) EG_free(line);

        return status;
}


// ********************** AIM Function Break *****************************
void aimCleanup(/*@unused@*/ void *instStore)
{

#ifdef DEBUG
    printf(" tsfoilAIM/aimCleanup!\n");
#endif

}
