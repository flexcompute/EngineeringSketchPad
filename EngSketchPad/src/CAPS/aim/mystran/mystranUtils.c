// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
#include "capsTypes.h" // Bring in CAPS types

#include "aimUtil.h"
#include "miscUtils.h" //Bring in misc. utility functions
#include "mystranUtils.h" // Bring in mystran utility header

#ifdef WIN32
#define strcasecmp  stricmp
#endif


// Read data from Mystran OUTPUT4 file and load it into a capsvalue
int mystran_readOutput4Data(FILE *fp, const char *keyword, capsValue *val)
{
    int matrix, i; // Indexing

    int sint = sizeof(int), sdouble = sizeof(double); // Size of variables

    int maxNumMatrix = 5; // Max number of matrices in a OUTPUT4 file

    int found = (int) false; // Boolean test

    // valueMatrix data information
    char valueMatrixHeader[9];
    int numCol=0, numRow=0, precision =0;
    int tempInt;
    double tempDouble;

    double *valueMatrix = NULL; // Temporary storage for valueMatrix

    //printf("Reading (unformatted) OUTPUT4 file....\n");

    /*
        ! Write matrix header

        WRITE(UNT) NCOLS, NROWS, FORM, precision, MAT_OUT_NAME(1:4), MAT_OUT_NAME(5:8)

        ! Write matrix data

        DO J=1,NCOLS
          WRITE(UNT) J, ROW_BEG, 2*NROWS, (MAT(I,J),I=1,NROWS)
        ENDDO

        ! Write matrix trailer

        WRITE(UNT) NCOLS+1, IROW, precision, (ZERO, I=1,precision)
     */

    for (matrix = 0; matrix < maxNumMatrix; matrix++ ) {

        /////   HEADER ///////

        // Number for bytes in line
        fread(&tempInt, sint, 1, fp);
        //printf("NumBytes = %d\n", tempInt);

        // Number of total columns in matrix
        fread(&numCol,   sint ,1,fp);
        //printf("Col = %d\n", numCol);

        // Number of total rows in matrix
        fread(&numRow,   sint ,1,fp);
        //printf("Row = %d\n", numRow);

        // Form
        fread(&tempInt, sint, 1, fp);
        //printf("Form = %d\n", tempInt);

        // Precision
        fread(&precision, sint, 1, fp);
        //printf("Precision = %d\n", precision);

        // Header title
        fread(valueMatrixHeader,  8*sizeof(char), 1, fp);
        valueMatrixHeader[8]= '\0';
        //printf("valueMatrixHeader - %s\n", valueMatrixHeader);

        // End for bytes in line
        fread(&tempInt, sint, 1, fp);
        //printf("NumBytesEnd = %d\n", tempInt);

        // Check the keyword header
        if (strcasecmp(valueMatrixHeader,keyword) == 0) {
            found = (int) true;
        }

        // Allocate valueMatrix array
        valueMatrix = (double *) EG_alloc(numCol*numRow*sdouble);
        if (valueMatrix == NULL) return EGADS_MALLOC;

        //////   Loop through matrix //////
        for (i = 0; i < numCol; i++) {

            // Number for bytes in line
            fread(&tempInt, sint, 1, fp);
            //printf("NumBytes = %d\n", tempInt);

            // Column number
            fread(&tempInt, sint, 1, fp);
            //printf("Column = %d\n", tempInt);

            // Row_Beg
            fread(&tempInt, sint, 1, fp);
            //printf("Row_Beg = %d\n", tempInt);

            // 2*NROWS
            fread(&tempInt, sint ,1,fp);
            //printf("2*NROWS = %d\n", tempInt);

            // Read row
            fread(valueMatrix + i, sdouble, numRow, fp);

            /*
            for (j = 0; j < numRow; j++) {
                Zprintf("valueMatrix = %f\n", valueMatrix[j+i]);
            }
             */

            // End for bytes in line
            fread(&tempInt, sint, 1, fp);
            //printf("NumBytesEnd = %d\n", tempInt);
        }


        ///////  Trailer  ///////

        // Number for bytes in line
        fread(&tempInt, sint, 1, fp);
        //printf("NumBytes = %d\n", tempInt);

        // Ncols+1
        fread(&tempInt, sint, 1, fp);
        //printf("NCOLS+1 = %d\n",tempInt);

        // IROW
        fread(&tempInt, sint, 1, fp);
        //printf("IROW = %d\n",tempInt);

        // Precision
        fread(&precision, sint, 1, fp);
        //printf("Precision = %d\n", precision);

        // Zeros
        for (i = 0; i < precision; i++) fread(&tempDouble, sdouble, 1, fp);

        // End for bytes in line
        fread(&tempInt, sint, 1, fp);
        //printf("NumBytesEnd = %d\n", tempInt);

        if (found == (int) true) {

            break;

        } else {

            if (valueMatrix != NULL) EG_free(valueMatrix);
            valueMatrix = NULL;
        }
    }

    /*
    // Print out matrix
    for (j = 0; j < numCol; j++) {
        printf("Col %d\n", j);

        for (i = 0; i < numRow; i++) {
            printf("valueMatrix= %e\n", valueMatrix[j+i]);
        }
    }
    */

    if ((found == (int) true) && (valueMatrix != NULL)) {
        // Copy valueMatrixs in capsvalueMatrix structure
        val->nrow   = numRow;
        val->ncol   = numCol;
        val->length = val->nrow*val->ncol;
        if (val->length == 1) val->dim = Scalar;
        else if (val->nrow == 1 || val->ncol == 1) val->dim = Vector;
        else val->dim = Array2D;

        if (val->length == 1) {
            val->vals.real = valueMatrix[0];

        } else {

            val->vals.reals = (double *) EG_alloc(val->length*sizeof(double));
            if (val->vals.reals == NULL) {
                if (valueMatrix != NULL) EG_free(valueMatrix);
                return EGADS_MALLOC;
            }

            memcpy(val->vals.reals, valueMatrix, val->length*sizeof(double));
        }

        if (valueMatrix != NULL) EG_free(valueMatrix);


        return CAPS_SUCCESS;
    } else {

        if (valueMatrix != NULL) EG_free(valueMatrix);
        return CAPS_NOTFOUND;

    }
}

// Read data from a Mystran F06 file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int mystran_readF06EigenVector(FILE *fp, int *numEigenVector, int *numGridPoint,
                               double ***dataMatrix)
{
    int status; // Function return

    int i, j, eigenValue; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *numEigenLine = "                                NUMBER OF EIGENVALUES EXTRACTED  . . . . . .";
    char *outputEigenLine = " OUTPUT FOR EIGENVECTOR        ";
    char *beginEigenLine=NULL;
    char *endEigenLine = "                         ------------- ------------- ------------- ------------- ------------- -------------";

    int numVariable = 8; // Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
    int intLength;

    printf("Reading Mystran FO6 file - extracting Eigen-Vectors!\n");

    *numEigenVector = 0;
    *numGridPoint = 0;

    // Loop through file line by line until we have determined how many Eigen-Values and grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // See how many Eigen-Values we have
        if (strncmp(numEigenLine, line, strlen(numEigenLine)) == 0) {
            sscanf(&line[strlen(numEigenLine)], "%d", numEigenVector);

            // Build begin Eigen-Value string
            beginEigenLine = (char *) EG_alloc((strlen(outputEigenLine)+2)*
                                               sizeof(char));
            if (beginEigenLine == NULL) {
                if (line != NULL) EG_free(line);
                return EGADS_MALLOC;
            }

            sprintf(beginEigenLine, "%s%d", outputEigenLine, 1);
            beginEigenLine[strlen(outputEigenLine)+1] = '\0';
        }

        // Once we know how many Eigen-Values we have, we need to determine how many grid points exist
        if ((*numEigenVector > 0) && (beginEigenLine != NULL)) {

            // Look for start of Eigen-Vector 1
            if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

                // Fast forward 5 lines
                for (i = 0; i < 5; i++) {
                    status = getline(&line, &linecap, fp);
                    if (status < 0) break;
                }

                // Loop through lines counting the number of grid points
                while (getline(&line, &linecap, fp) >= 0) {
                    *numGridPoint +=1;
                    if (strncmp(endEigenLine, line,strlen(endEigenLine)) == 0) {
                        *numGridPoint -= 1; // Get rid of last line
                        break;
                    }
                }
            }
        }
    }

    printf("\tNumber of Eigen-Vectors = %d\n", *numEigenVector);
    printf("\tNumber of Grid Points = %d for each Eigen-Vector\n",
           *numGridPoint);

    // Free begin Eigen-Value string
    if (beginEigenLine != NULL) EG_free(beginEigenLine);
    beginEigenLine = NULL;

    if (*numGridPoint == 0 || *numEigenVector == 0) {
        printf("\tEither the number of data points  = 0 and/or the number of Eigen-Values = 0!!!\n");
        printf("\tWas a modal analysis run?\n");
        if (line != NULL) EG_free(line);
        return CAPS_NOTFOUND;
    }
    // Rewind the file
    rewind(fp);

    // Allocate dataMatrix array
    if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numEigenVector *sizeof(double *));
    if (*dataMatrix == NULL) {
        if (line != NULL) EG_free(line);
        return EGADS_MALLOC; // If allocation failed ....
    }

    for (i = 0; i < *numEigenVector; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(*numGridPoint*numVariable*
                                               sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));

            if (line != NULL) EG_free(line);
            return EGADS_MALLOC;
        }
    }

    if      (*numEigenVector >= 1000) intLength = 4;
    else if (*numEigenVector >= 100) intLength = 3;
    else if (*numEigenVector >= 10) intLength = 2;
    else intLength = 1;

    // Loop through the file again and pull out data
    eigenValue = 1;
    while (eigenValue <= *numEigenVector) {

        if (beginEigenLine == NULL) {
            // Build begin Eigen-Value string
            beginEigenLine = (char *) EG_alloc((strlen(outputEigenLine)+
                                                intLength+1)*sizeof(char));
            if (beginEigenLine == NULL) {
                if (line != NULL) EG_free(line);
                return EGADS_MALLOC;
            }

            sprintf(beginEigenLine,"%s%d",outputEigenLine, eigenValue);
            beginEigenLine[strlen(outputEigenLine)+intLength] = '\0';
        }

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // Look for start of Eigen-Vector
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            printf("\tLoading Eigen-Vector = %d\n", eigenValue);

            // Fast forward 5 lines
            for (i = 0; i < 5; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through the file and fill up the data matrix
            for (j = 0; j < (*numGridPoint)*numVariable; j++) {
                // eigenValue is 1 bias
                fscanf(fp, "%lf", &(*dataMatrix)[eigenValue-1][j]);
            }

            if (beginEigenLine != NULL) EG_free(beginEigenLine);
            beginEigenLine = NULL;

            eigenValue += 1;
        }
    }

    if (beginEigenLine != NULL) EG_free(beginEigenLine);

    if (line != NULL) EG_free(line);

    return CAPS_SUCCESS;
}

// Read data from a Mystran F06 file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int mystran_readF06Displacement(FILE *fp, int subcaseId, int *numGridPoint,
                                double ***dataMatrix)
{
    int status; // Function return

    int i, j; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *outputSubcaseLine = " OUTPUT FOR SUBCASE        ";
    char *beginSubcaseLine=NULL;
    char *endSubcaseLine = "                         ------------- ------------- ------------- ------------- ------------- -------------";

    int numVariable = 8; // Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
    int intLength;
    int numDataRead = 0;

    printf("Reading Mystran FO6 file - extracting Displacements!\n");

    *numGridPoint = 0;

    if      (subcaseId >= 1000) intLength = 4;
    else if (subcaseId >= 100) intLength = 3;
    else if (subcaseId >= 10) intLength = 2;
    else intLength = 1;

    beginSubcaseLine = (char *) EG_alloc((strlen(outputSubcaseLine)+intLength+1)*
                                         sizeof(char));
    if (beginSubcaseLine == NULL) return EGADS_MALLOC;

    sprintf(beginSubcaseLine,"%s%d",outputSubcaseLine, subcaseId);
    beginSubcaseLine[strlen(outputSubcaseLine)+intLength] = '\0';

    // Loop through file line by line until we have determined how many grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // Look for start of subcaseId
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            // Fast forward 5 lines
            for (i = 0; i < 5; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through lines counting the number of grid points
            while (getline(&line, &linecap, fp) >= 0) {
                if (strncmp(endSubcaseLine, line,strlen(endSubcaseLine)) == 0) {
                    break;
                }
                *numGridPoint +=1;
            }
        }
    }

    printf("\tNumber of Grid Points = %d\n", *numGridPoint);

    if (*numGridPoint == 0) {
        printf("\tEither the number of data points  = 0 and/or subcase wasn't found!!!\n");

        if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
        if (line != NULL) EG_free(line);
        return CAPS_NOTFOUND;
    }

    // Rewind the file
    rewind(fp);

    // Allocate dataMatrix array
    if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numGridPoint *sizeof(double *));
    if (*dataMatrix == NULL) {
        if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
        return EGADS_MALLOC; // If allocation failed ....
    }

    for (i = 0; i < *numGridPoint; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));
            if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
            return EGADS_MALLOC;
        }
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {
        if (line == NULL) continue;

        // Look for start of Eigen-Vector
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            printf("\tLoading displacements for Subcase = %d\n", subcaseId);

            // Fast forward 5 lines
            for (i = 0; i < 5; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through the file and fill up the data matrix
            for (i = 0; i < (*numGridPoint); i++) {
                for (j = 0; j < numVariable; j++) {

                    if (fscanf(fp, "%lf", &(*dataMatrix)[i][j]) == 1)
                        numDataRead++;
                }
            }

            break;
        }
    }

    if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);

    if (line != NULL) EG_free(line);

    if (numDataRead/numVariable != *numGridPoint) {
        printf("Failed to read %d grid points. Only found %d.\n",
               *numGridPoint, numDataRead/numVariable);
        return CAPS_IOERR;
    }

    return CAPS_SUCCESS;
}
