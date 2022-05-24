#include <string.h>
#include <stdio.h>

#include "egads.h"     // Bring in egads utilss
#include "capsTypes.h" // Bring in CAPS types
#include "aimUtil.h"   // Bring in AIM utils
#include "aimMesh.h"// Bring in AIM meshing utils

#include "miscUtils.h" // Bring in misc. utility functions
#include "meshUtils.h" // Bring in meshing utility functions
#include "cfdTypes.h"  // Bring in cfd specific types
#include "su2Utils.h"  // Bring in su2 utility header

#ifdef WIN32
#define strcasecmp stricmp
#define strtok_r   strtok_s
#endif

// Extract flow variables from SU2 surface csv file (connectivity is ignored) - dataMatrix = [numVariable][numDataPoint]
int su2_readAeroLoad(void *aimInfo, char *filename, int *numVariable, char **variableName[],
                     int *numDataPoint, double ***dataMatrix)
{

    int status; // Function return
    int i, j; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder
    char *tempStr = NULL; // Temporary strings

    char comma;

    FILE *fp = NULL; // File pointer

    // Open file
    fp = aim_fopen(aimInfo, filename, "r");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        return CAPS_IOERR;
    }

    printf("Reading SU2 AeroLoad File - %s\n", filename);

    *numVariable = 0;
    *numDataPoint = 0;

    // Determine how many variables there are - get the first line
    status = getline(&line, &linecap, fp);
    if (status < 0) {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }
    if (line == NULL) {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Create a temp string to braket our line with '[' and ']'
    AIM_ALLOC(tempStr, strlen(line) + 4, char, aimInfo, status);

    strcpy(tempStr,"[");
    strncat(tempStr, line, strlen(line)-1);
    strcat(tempStr, "]");

    // Sort string into an array of strings
    status =  string_toStringDynamicArray(tempStr, numVariable, variableName);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Count the number of lines to determine number of data points
    while (getline(&line, &linecap, fp) >= 0) {
        *numDataPoint += 1;
    }

    // Go to the start of the file
    rewind(fp);

    printf("\tNumber of variables = %d\n", *numVariable);
    printf("\tNumber of data points = %d\n", *numDataPoint);

    // Re-Grab the first line
    status = getline(&line, &linecap, fp);
    if (status < 0) {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // If nozero values were found
    if (*numVariable != 0 && *numDataPoint != 0) {

        // Allocate dataMatrix array
        //if (*dataMatrix != NULL) EG_free(*dataMatrix);

        *dataMatrix = (double **) EG_alloc(*numVariable *sizeof(double *));
        if (*dataMatrix == NULL) {
            status = EGADS_MALLOC; // If allocation failed ....
            goto cleanup;
        }

        for (i = 0; i < *numVariable; i++) {

            (*dataMatrix)[i] = (double *) EG_alloc(*numDataPoint*sizeof(double));

            if ((*dataMatrix)[i] == NULL) { // If allocation failed ....

                for (j = 0; j < i; j++) {

                    if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
                }

                if ((*dataMatrix) != NULL) EG_free((*dataMatrix));
                *dataMatrix = NULL;

                status = EGADS_MALLOC;
                goto cleanup;
            }
        }

        // Loop through the file and fill up the data matrix
        for (j = 0; j < *numDataPoint; j++) {

            for (i = 0; i < *numVariable; i++) {
                fscanf(fp, "%lf", &(*dataMatrix)[i][j]);

                if (i != *numVariable-1) fscanf(fp, "%c", &comma);
            }
        }

        // Output the first two rows of the dataMatrix
        //for (i = 0; i < *numVariable; i++) printf("Variable %d - %.6f\n", i, (*dataMatrix)[i][0]);
        //for (i = 0; i < *numVariable; i++) printf("Variable %d - %.6f\n", i, (*dataMatrix)[i][1]);

    } else {

        AIM_ERROR(aimInfo, "\tNo data values extracted from file - %s",filename);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    status = CAPS_SUCCESS;

    cleanup:

        if (status != CAPS_SUCCESS) printf("\tError in su2_readAeroLoad, status = %d\n", status);

        if (fp != NULL) fclose(fp);

        if (tempStr != NULL) EG_free(tempStr);
        if (line != NULL) EG_free(line);

        if (status != CAPS_SUCCESS) {

            // Free up the dataMatrix
            if ((*dataMatrix) != NULL) {

                for (i = 0; i < *numVariable; i++) {
                    if ((*dataMatrix)[i] == NULL) EG_free((*dataMatrix)[i]);
                }

                EG_free((*dataMatrix));
                *dataMatrix = NULL;
            }

            // Free of variable array
            string_freeArray(*numVariable, variableName);

            *numVariable = 0;
            *numDataPoint = 0;
        }

        return status;
}


// Write SU2 surface motion file (connectivity is optional) - dataMatrix = [numVariable][numDataPoint], connectMatrix (optional) = [4*numConnect]
//  the formating of the data may be specified through dataFormat = [numVariable] (use capsTypes Integer and Double)- If NULL default to double
int su2_writeSurfaceMotion(void *aimInfo,
                           char *filename,
                           int numVariable,
                           int numDataPoint,
                           double **dataMatrix,
                           int *dataFormat,
                           int numConnect,
                           int *connectMatrix)
{

    int i, j; // Indexing

    FILE *fp; // File pointer

    printf("Writing SU2 Motion File - %s\n", filename);

    // Open file
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        return CAPS_IOERR;
    }

    if (connectMatrix == NULL) numConnect = 0;

    for (i = 0; i < numDataPoint; i++) {

        for (j = 0; j < numVariable; j++) {

            if (dataFormat != NULL) {

                if (dataFormat[j] == (int) Integer) {

                    fprintf(fp, "%d ", (int) dataMatrix[j][i]);

                } else if (dataFormat[j] == (int) Double) {

                    fprintf(fp, "%.18e ", dataMatrix[j][i]);

                } else {

                    AIM_ERROR(aimInfo, "Unrecognized data format requested - %d", (int) dataFormat[j]);
                    fclose(fp);
                    return CAPS_BADVALUE;
                }

            } else {

                fprintf(fp, "%.18e ", dataMatrix[j][i]);
            }
        }

        fprintf(fp, "\n");
    }

    if (connectMatrix != NULL) {

        for (i = 0; i < numConnect; i++ ) {

            for (j = 0; j < 4; j++) {
                fprintf(fp, "%d ", connectMatrix[4*i+j]);
            }

            fprintf(fp, "\n");
        }
    }

    fclose(fp);

    return CAPS_SUCCESS;
}

// Write SU2 data transfer files
int su2_dataTransfer(void *aimInfo,
                     char *projectName,
                     aimMeshRef *meshRef)
{

    /*! \page dataTransferSU2 SU2 Data Transfer
     *
     * \section dataToSU2 Data transfer to SU2 (FieldIn)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements (as from a structural solver)
     *   and updates SU2's surface mesh; a new [project_name]_motion.dat file is written out which may
     *   be loaded into SU2 to update the surface mesh/move the volume mesh.
     * </ul>
     */

    int status; // Function return status
    int i, j, ibound, ibody, iglobal; // Indexing

    int stringLength = 0;

    char *filename = NULL;

    // Discrete data transfer variables
    capsDiscr *discr;
    char **boundNames = NULL;
    int numBoundName = 0;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char *units;

    int state, nGlobal, *globalOffset=NULL;
    ego body;

    // Variables used in global node mapping
    int ptype, pindex;
    double xyz[3];

    // Data transfer Out variables

    double **dataOutMatrix = NULL;
    int dataOutFormat[] = {Integer, Double, Double, Double};

    int numOutVariable = 4; // ID and x,y,x
    int numOutDataPoint = 0;

    const char fileExt[] = "_motion.dat";

    int foundDisplacement = (int) false;

    AIM_ALLOC(dataOutMatrix, numOutVariable, double*, aimInfo, status);
    for (i = 0; i < numOutVariable; i++) dataOutMatrix[i] = NULL;

    status = aim_getBounds(aimInfo, &numBoundName, &boundNames);
    AIM_STATUS(aimInfo, status);

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
        printf("No recognized data transfer names found!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    printf("Writing SU2 data transfer files\n");

    // first construct the complete boundary mesh
    // SU2 will initailize all active MARKER boundaries to the motion file values,
    // so the safest thing to do is write out all boundary points
    AIM_ALLOC(globalOffset, meshRef->nmap+1, int, aimInfo, status);
    numOutDataPoint = 0;
    globalOffset[0] = 0;
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nGlobal);
      AIM_STATUS(aimInfo, status);

      // re-allocate data arrays
      for (j = 0; j < numOutVariable; j++) {
        AIM_REALL(dataOutMatrix[j], numOutDataPoint + nGlobal, double, aimInfo, status);
      }

      for (iglobal = 0; iglobal < nGlobal; iglobal++) {
        status = EG_getGlobal(meshRef->maps[i].tess,
                              iglobal+1, &ptype, &pindex, xyz);
        AIM_STATUS(aimInfo, status);

        // Volume mesh node ID (0-based for SU2)
        dataOutMatrix[0][globalOffset[i]+iglobal] = meshRef->maps[i].map[iglobal]-1;

        // First just set the Coordinates
        dataOutMatrix[1][globalOffset[i]+iglobal] = xyz[0];
        dataOutMatrix[2][globalOffset[i]+iglobal] = xyz[1];
        dataOutMatrix[3][globalOffset[i]+iglobal] = xyz[2];
      }

      numOutDataPoint += nGlobal;
      globalOffset[i+1] = globalOffset[i] + nGlobal;
    }

    // now apply the displacements
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

        ibody   = discr->tessGlobal[2*i+0];
        iglobal = discr->tessGlobal[2*i+1];

        status = EG_getGlobal(discr->bodys[ibody-1].tess,
                              iglobal, &ptype, &pindex, xyz);
        AIM_STATUS(aimInfo, status);

        // Find the disc tessellation in the original list of tessellations
        for (j = 0; j < meshRef->nmap; j++) {
          if (discr->bodys[ibody-1].tess == meshRef->maps[j].tess) {
            break;
          }
        }
        if (j == meshRef->nmap) {
          AIM_ERROR(aimInfo, "Could not find matching tessellation!");
          status = CAPS_MISMATCH;
          goto cleanup;
        }

        if (numDataTransferPoint == 1) {
          // A single point means this is an initialization phase
          dataOutMatrix[1][globalOffset[j]+iglobal-1] += dataTransferData[0];
          dataOutMatrix[2][globalOffset[j]+iglobal-1] += dataTransferData[1];
          dataOutMatrix[3][globalOffset[j]+iglobal-1] += dataTransferData[2];
        } else {
          // Apply delta displacements
          dataOutMatrix[1][globalOffset[j]+iglobal-1] += dataTransferData[3*i+0];
          dataOutMatrix[2][globalOffset[j]+iglobal-1] += dataTransferData[3*i+1];
          dataOutMatrix[3][globalOffset[j]+iglobal-1] += dataTransferData[3*i+2];
        }
      }
    } // Loop through bound names

    stringLength = strlen(projectName) + strlen(fileExt) + 1;
    AIM_ALLOC(filename, stringLength, char, aimInfo, status);
    strcpy(filename, projectName);
    strcat(filename, fileExt);

    // Write out displacement in tecplot file
/*@-nullpass@*/
    status = su2_writeSurfaceMotion(aimInfo,
                                    filename,
                                    numOutVariable,
                                    numOutDataPoint,
                                    dataOutMatrix,
                                    dataOutFormat,
                                    0, // numConnectivity
                                    NULL); // connectivity matrix
/*@+nullpass@*/
    AIM_STATUS(aimInfo, status);

    status = CAPS_SUCCESS;

// Clean-up
cleanup:

    if (dataOutMatrix != NULL) {
        for (i = 0; i < numOutVariable; i++) {
          AIM_FREE(dataOutMatrix[i]);
        }
    }

    AIM_FREE(dataOutMatrix);
    AIM_FREE(filename);
    AIM_FREE(boundNames);
    AIM_FREE(globalOffset);

    return status;

}

// Extract the boundary condition names that should added to the marker
int su2_marker(void *aimInfo, const char* iname, capsValue *aimInputs, FILE *fp,
               cfdBoundaryConditionStruct bcProps)
{

    int status = CAPS_SUCCESS;
    int counter = 0;
    int nmarker = 0;
    capsValue *markerValue;

    int i, j; // indexing

    char *marker = NULL;

    markerValue = &aimInputs[aim_getIndex(aimInfo, iname, ANALYSISIN)-1];

    // might not be anything to write in the list
    if (markerValue->nullVal == IsNull) {
        fprintf(fp," NONE )\n");
        return CAPS_SUCCESS;
    }

    nmarker = markerValue->length;
    marker = markerValue->vals.string;
    counter = 0;
    for (j = 0; j < nmarker; j++) {
        for (i = 0; i < bcProps.numSurfaceProp ; i++) {
            if (strcmp(bcProps.surfaceProp[i].name, marker) == 0) {

                if (counter > 0) fprintf(fp, ",");
                fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

                counter += 1;
                break;
            }
        }
        marker += strlen(marker)+1;
    }
    fprintf(fp," )\n");

    if (counter != nmarker || counter == 0) {
        AIM_ERROR(aimInfo, "Could not find all '%s' names:\n", iname);
        marker = markerValue->vals.string;
        for (j = 0; j < nmarker; j++) {
            AIM_ADDLINE(aimInfo, "\t%s", marker);
            marker += strlen(marker)+1;
        }
        AIM_ADDLINE(aimInfo, "");

        AIM_ADDLINE(aimInfo, "in the list of boundary condition names:\n");
        for (i = 0; i < bcProps.numSurfaceProp ; i++)
          AIM_ADDLINE(aimInfo, "\t%s", bcProps.surfaceProp[i].name);
        AIM_ADDLINE(aimInfo, "");

        status = CAPS_NOTFOUND;
    }

//cleanup:
    return status;
}

int
su2_unitSystem(const char *unitSystem,
               const char **length,
               const char **mass,
               const char **temperature,
               const char **force,
               const char **pressure,
               const char **density,
               const char **speed,
               const char **viscosity,
               const char **area)
{
  /*@-observertrans@*/
  if (strcasecmp(unitSystem, "SI") == 0) {
    *length      = "meters";
    *mass        = "kilograms";
    *temperature = "Kelvin";
    *force       = "Newton";
    *pressure    = "Pascal";
    *density     = "kg/m^3";
    *speed       = "m/s";
    *viscosity   = "N*s/m^2";
    *area        = "m^2";
  } else if (strcasecmp(unitSystem, "US") == 0) {
    *length      = "inches";
    *mass        = "slug";
    *temperature = "Rankines";
    *force       = "lbf";        // slug*ft/s^2
    *pressure    = "lbf/ft^2";   // psf
    *density     = "slug/ft^3";
    *speed       = "ft/s";
    *viscosity   = "lbf*s/ft^2";
    *area        = "ft^2";
  } else {
    return CAPS_UNITERR;
  }
  /*@+observertrans@*/

  return CAPS_SUCCESS;
}

#ifdef DEFINED_BUT_NOT_USED
// Use in the future to just the information to carry out the deform.
static int su2_writeConfig_Deform(void *aimInfo, capsValue *aimInputs,
                                  cfdBoundaryConditionStruct bcProps)
{
    int status; // Function return status

    int i; // Indexing

    int stringLength;

    // For SU2 boundary tagging
    int counter;

    char *filename = NULL;
    FILE *fp = NULL;
    char fileExt[] = "_Deform.cfg";

    printf("Write SU2 configuration file for grid deformation\n");
    stringLength = 1
                   + strlen(aimInputs[Proj_Name-1].vals.string)
                   + strlen(fileExt);

    filename = (char *) EG_alloc((stringLength +1)*sizeof(char));
    if (filename == NULL) {
        status =  EGADS_MALLOC;
        goto cleanup;
    }

    strcpy(filename, aimInputs[Proj_Name-1].vals.string);
    strcat(filename, fileExt);

    fp = aim_fopen(aimInfo, filename,"w");
    if (fp == NULL) {
        status =  CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    fprintf(fp,"%%                                                                              %%\n");
    fprintf(fp,"%% SU2 configuration file - for Grid Deformation                                %%\n");
    fprintf(fp,"%% Created by SU2AIM for Project: \"%s\"\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");

    fprintf(fp,"%% ----------------------- DYNAMIC MESH DEFINITION -----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Type of dynamic mesh (NONE, RIGID_MOTION, DEFORMING, ROTATING_FRAME,\n");
    fprintf(fp,"%%                       MOVING_WALL, STEADY_TRANSLATION, FLUID_STRUCTURE,\n");
    fprintf(fp,"%%                       AEROELASTIC, ELASTICITY, EXTERNAL,\n");
    fprintf(fp,"%%                       AEROELASTIC_RIGID_MOTION, GUST)\n");
    fprintf(fp,"GRID_MOVEMENT_KIND= DEFORMING\n");

    fprintf(fp,"%% Moving wall boundary marker(s) (NONE = no marker, ignored for RIGID_MOTION)\n");

    fprintf(fp,"MARKER_MOVING= (" );
    counter = 0;
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
            bcProps.surfaceProp[i].surfaceType == Viscous) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%% ----------------------- DESIGN VARIABLE PARAMETERS --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Kind of deformation (NO_DEFORMATION, TRANSLATION, ROTATION, SCALE,\n");
    fprintf(fp,"%%                      FFD_SETTING, FFD_NACELLE\n");
    fprintf(fp,"%%                      FFD_CONTROL_POINT, FFD_CAMBER, FFD_THICKNESS, FFD_TWIST\n");
    fprintf(fp,"%%                      FFD_CONTROL_POINT_2D, FFD_CAMBER_2D, FFD_THICKNESS_2D, FFD_TWIST_2D,\n");
    fprintf(fp,"%%                      HICKS_HENNE, SURFACE_BUMP)\n");
    fprintf(fp,"DV_KIND= SURFACE_FILE \n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker of the surface in which we are going apply the shape deformation\n");

    fprintf(fp,"DV_MARKER= (");
    counter = 0;
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
            bcProps.surfaceProp[i].surfaceType == Viscous) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    //fprintf(fp,"DV_MARKER= ( airfoil )\n");
    fprintf(fp,"MOTION_FILENAME=%s_motion.dat\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"\n");

    fprintf(fp,"%% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh input file\n");
    fprintf(fp,"MESH_FILENAME= %s.su2\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh input file format (SU2, CGNS)\n");
    fprintf(fp,"MESH_FORMAT= SU2\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Mesh output file\n");
    fprintf(fp,"MESH_OUT_FILENAME= %s.su2\n", aimInputs[Proj_Name-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW,\n");
    fprintf(fp,"%%                     FIELDVIEW, FIELDVIEW_BINARY)\n");
    string_toUpperCase(aimInputs[Output_Format-1].vals.string);
    fprintf(fp,"OUTPUT_FORMAT= %s\n", aimInputs[Output_Format-1].vals.string);
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Verbosity of console output: NONE removes minor MPI overhead (NONE, HIGH)\n");
    fprintf(fp,"CONSOLE_OUTPUT_VERBOSITY= HIGH\n");


    fprintf(fp,"%% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Euler wall boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_EULER= (" );

    counter = 0; // Euler boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {

        if (bcProps.surfaceProp[i].surfaceType == Inviscid) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( marker name, constant heat flux (J/m^2), ... )\n");
    fprintf(fp,"MARKER_HEATFLUX= (");

    counter = 0; // Viscous boundary w/ heat flux
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Viscous &&
            bcProps.surfaceProp[i].wallTemperatureFlag == (int) true &&
            bcProps.surfaceProp[i].wallTemperature < 0) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID, bcProps.surfaceProp[i].wallHeatFlux);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Navier-Stokes (no-slip), isothermal wall marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( marker name, constant wall temperature (K), ... )\n");
    fprintf(fp,"MARKER_ISOTHERMAL= (");

    counter = 0; // Viscous boundary w/ isothermal wall
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Viscous &&
            bcProps.surfaceProp[i].wallTemperatureFlag == (int) true &&
            bcProps.surfaceProp[i].wallTemperature >= 0) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID, bcProps.surfaceProp[i].wallTemperature);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Far-field boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_FAR= (" );

    counter = 0; // Farfield boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Farfield) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Symmetry boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_SYM= (" );

    counter = 0; // Symmetry boundary
    for (i = 0; i < bcProps.numSurfaceProp; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Symmetry) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Near-Field boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_NEARFIELD= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Zone interface boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_INTERFACE= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Actuator disk boundary type (VARIABLES_JUMP, NET_THRUST, BC_THRUST,\n");
    fprintf(fp,"%%                              DRAG_MINUS_THRUST, MASSFLOW, POWER)\n");
    fprintf(fp,"ACTDISK_TYPE= VARIABLES_JUMP\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Actuator disk boundary marker(s) with the following formats (NONE = no marker)\n");
    fprintf(fp,"%% Variables Jump: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%                   Takeoff pressure jump (psf), Takeoff temperature jump (R), Takeoff rev/min,\n");
    fprintf(fp,"%%                   Cruise  pressure jump (psf), Cruise temperature jump (R), Cruise rev/min )\n");
    fprintf(fp,"%% Net Thrust: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%               Takeoff net thrust (lbs), 0.0, Takeoff rev/min,\n");
    fprintf(fp,"%%               Cruise net thrust (lbs), 0.0, Cruise rev/min )\n");
    fprintf(fp,"%%BC Thrust: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%             Takeoff BC thrust (lbs), 0.0, Takeoff rev/min,\n");
    fprintf(fp,"%%             Cruise BC thrust (lbs), 0.0, Cruise rev/min )\n");
    fprintf(fp,"%%Drag-Thrust: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%               Takeoff Drag-Thrust (lbs), 0.0, Takeoff rev/min,\n");
    fprintf(fp,"%%               Cruise Drag-Thrust (lbs), 0.0, Cruise rev/min )\n");
    fprintf(fp,"%%MasssFlow: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%               Takeoff massflow (lbs/s), 0.0, Takeoff rev/min,\n");
    fprintf(fp,"%%               Cruise massflowt (lbs/s), 0.0, Cruise rev/min )\n");
    fprintf(fp,"%%Power: ( inlet face marker, outlet face marker,\n");
    fprintf(fp,"%%          Takeoff power (HP), 0.0, Takeoff rev/min\n");
    fprintf(fp,"%%          Cruise power (HP), 0.0, Cruise rev/min )\n");
    fprintf(fp,"MARKER_ACTDISK= ( NONE )\n");
    fprintf(fp,"%%\n");

    fprintf(fp,"%% Inlet boundary type (TOTAL_CONDITIONS, MASS_FLOW)\n");
    fprintf(fp,"INLET_TYPE= TOTAL_CONDITIONS\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Inlet boundary marker(s) with the following formats (NONE = no marker) \n");
    fprintf(fp,"%% Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x, \n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"%% Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x, \n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"%% Incompressible: (inlet marker, NULL, velocity magnitude, flow_direction_x,\n");
    fprintf(fp,"%%           flow_direction_y, flow_direction_z, ... ) where flow_direction is\n");
    fprintf(fp,"%%           a unit vector.\n");
    fprintf(fp,"MARKER_INLET= ( ");

    counter = 0; // Subsonic Inflow
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == SubsonicInflow) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f, %f, %f, %f, %f", bcProps.surfaceProp[i].bcID,
                                                  bcProps.surfaceProp[i].totalTemperature,
                                                  bcProps.surfaceProp[i].totalPressure,
                                                  bcProps.surfaceProp[i].uVelocity,
                                                  bcProps.surfaceProp[i].vVelocity,
                                                  bcProps.surfaceProp[i].wVelocity);
            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Supersonic inlet boundary marker(s) (NONE = no marker) \n");
    fprintf(fp,"%% Format: (inlet marker, temperature, static pressure, velocity_x, \n");
    fprintf(fp,"%%           velocity_y, velocity_z, ... ), i.e. primitive variables specified.\n");
    fprintf(fp,"MARKER_SUPERSONIC_INLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Outlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( outlet marker, back pressure (static), ... )\n");
    fprintf(fp,"MARKER_OUTLET= ( ");

    counter = 0; // Outlet boundary
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == BackPressure ||
            bcProps.surfaceProp[i].surfaceType == SubsonicOutflow) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d, %f", bcProps.surfaceProp[i].bcID,
                                  bcProps.surfaceProp[i].staticPressure);

            counter += 1;
        }
    }

    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    fprintf(fp,"%%\n");
    fprintf(fp,"%% Supersonic outlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_SUPERSONIC_OUTLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Periodic boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y, \n");
    fprintf(fp,"%% rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis, \n");
    fprintf(fp,"%% rotation_angle_z-axis, translation_x, translation_y, translation_z, ... )\n");
    fprintf(fp,"MARKER_PERIODIC= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine inflow boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (engine inflow marker, fan face Mach, ... )\n");
    fprintf(fp,"MARKER_ENGINE_INFLOW= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Engine exhaust boundary marker(s) with the following formats (NONE = no marker) \n");
    fprintf(fp,"%% Format: (engine exhaust marker, total nozzle temp, total nozzle pressure, ... )\n");
    fprintf(fp,"MARKER_ENGINE_EXHAUST= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Displacement boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( displacement marker, displacement value normal to the surface, ... )\n");
    fprintf(fp,"MARKER_NORMAL_DISPL= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Load boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( load marker, force value normal to the surface, ... )\n");
    fprintf(fp,"MARKER_NORMAL_LOAD= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Pressure boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: ( pressure marker )\n");
    fprintf(fp,"MARKER_PRESSURE= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Neumann bounday marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_NEUMANN= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Dirichlet boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"MARKER_DIRICHLET= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Riemann boundary marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (marker, data kind flag, list of data)\n");
    fprintf(fp,"MARKER_RIEMANN= ( NONE )\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Non Reflecting boundary conditions marker(s) (NONE = no marker)\n");
    fprintf(fp,"%% Format: (marker, data kind flag, list of data)\n");
    fprintf(fp,"MARKER_NRBC= ( NONE )\n");


    fprintf(fp,"%% ------------------------ SURFACES IDENTIFICATION ----------------------------%%\n");
    fprintf(fp,"%%\n");
    fprintf(fp,"%% Marker(s) of the surface in the surface flow solution file\n");
    fprintf(fp,"MARKER_PLOTTING= (" );

    counter = 0; // Surface marker
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {
        if (bcProps.surfaceProp[i].surfaceType == Inviscid ||
            bcProps.surfaceProp[i].surfaceType == Viscous) {

            if (counter > 0) fprintf(fp, ",");
            fprintf(fp," BC_%d", bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }
    if(counter == 0) fprintf(fp," NONE");
    fprintf(fp," )\n");

    // write monitoring information
    status = su2_monitor(aimInfo, aimInputs, fp, bcProps);
    if (status != CAPS_SUCCESS) goto cleanup;

    fprintf(fp,"%%\n");
    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:

        if (filename != NULL) EG_free(filename);
        if (fp != NULL) fclose(fp);

        return status;
}
#endif
