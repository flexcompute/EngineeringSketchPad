// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Tecplot  related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include <string.h>
#include <stdio.h>

#include "capsTypes.h"  // Bring in CAPS types
#include "aimUtil.h"

#include "errno.h"

#ifdef WIN32
#define strcasecmp  stricmp

// For _mkdir
#include "direct.h"

#else

// For mkdir
#include "sys/stat.h"

#endif


// Write FEPOINT Tecplot data (compatible with FUN3D) (connectivity is optional) - dataMatrix = [numVariable][numDataPoint],
// connectMatrix (optional) = [4*numConnect], the formating of the data may be specified through
// dataFormat = [numVariable] (use capsTypes Integer and Double)- If NULL default is a double
int tecplot_writeFEPOINT(void *aimInfo,
                         const char *filename,
                         const char *message,
                         const char *zoneTitle,
                         int numVariable,
                         char *variableName[],
                         int numDataPoint,
                         double **dataMatrix,
                         const int *dataFormat,
                         int numConnect,
                         int *connectMatrix,
                         double *solutionTime) {

    int i, j; // Indexing
    double time = 0.0;

    const char *zone;

    FILE *fp = NULL; // File pointer

    if (filename == NULL) return CAPS_NULLVALUE;
    if (variableName == NULL) return CAPS_NULLVALUE;
    if (dataMatrix == NULL) return CAPS_NULLVALUE;

    if (message != NULL) {
        printf("Writing %s File - %s\n", message, filename);
    } else {
        printf("Writing File - %s\n", filename);
    }

    // Open file
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        return CAPS_IOERR;
    }

    fprintf(fp,"title=\"%s\"\n","CAPS");

    fprintf(fp, "variables=");

    for (i = 0; i < numVariable;i++) {
        fprintf(fp, "\"%s\"",variableName[i]);
        if (i != numVariable-1) fprintf(fp, ",");
    }
    fprintf(fp,"\n");

    if (solutionTime != NULL) time = *solutionTime;

    if (connectMatrix == NULL) numConnect = 0;

    if (zoneTitle == NULL) zone = "caps";
    else zone = zoneTitle;

    fprintf(fp,"zone t=\"%s\", i=%d, j=%d, f=fepoint, solutiontime=%f, strandid=0\n", zone, numDataPoint, numConnect, time);

    for (i = 0; i < numDataPoint; i++) {

        for (j = 0; j < numVariable; j++) {

            if (dataFormat != NULL) {

                if (dataFormat[j] == (int) Integer) {

                    fprintf(fp, "%d ", (int) dataMatrix[j][i]);

                } else if (dataFormat[j] == (int) Double) {

                    fprintf(fp, "%e ", dataMatrix[j][i]);

                } else {

                    printf("Unrecognized data format requested - %d", (int) dataFormat[j]);
                    fclose(fp);
                    return CAPS_BADVALUE;
                }

            } else {

                fprintf(fp, "%e ", dataMatrix[j][i]);
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
