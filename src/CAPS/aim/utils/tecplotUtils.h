// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

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
                         double *solutionTime);
