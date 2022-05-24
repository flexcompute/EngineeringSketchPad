// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Miscellaneous related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#include <limits.h>
#endif

#include "egads.h"
#include "aimUtil.h"
#include "capsTypes.h"  // Bring in CAPS types
#include "miscTypes.h"  // Bring in misc. structures
#include "miscUtils.h"  // Bring in misc. utility header

#ifdef WIN32
#define snprintf _snprintf
#define strcasecmp stricmp
#endif

#define NINT(A)         (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define MAX(A,B)        (((A) > (B)) ? (A) : (B))

// Convert a string in tuple form to an array of strings - tuple is assumed to be bounded by '(' and ')' and comma separated
//  for example ("3.0", 5, "foo", ("f", 1, 4), [1,2,3]) - note the strings of the outer tuple should NOT contain commas, Tuple elements
//  consisting of internal tuples and arrays are okay (only 1 level deep).  If the string coming in is not a tuple the string
// is simply copied. Also quotations are removed from values elements of the (outer) tuple.
int json_parseTuple(/*@null@*/ const char *stringToParse, int *arraySize,
                    char **stringArray[])
{
    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    int haveArray = (int) false;
    char *quoteString = NULL; // Temporary string to hold the found string
    char *noQuoteString = NULL;  // Temporary string to hold the found string with quotation marks removed

    int foundInternalTuple = (int) false;
    int foundInternalArray = (int) false;

    // Debug function
    //printf("Search string - %s length %d\n", stringToParse, strlen(stringToParse));

    if (stringToParse == NULL) return CAPS_BADVALUE;

    // Check to make sure first and last of the incoming string denote an array
    if (stringToParse[0] == '[' ) {

        // Debug function
        //printf("[ found\n");

        if(stringToParse[strlen(stringToParse)-1] == ']') { // Lets count how many commas we have

            // Debug function
            //printf("] found\n");

            haveArray = (int) true;
            *arraySize = 1;
            for (i = 1; i < strlen(stringToParse)-1; i++) {

                if (stringToParse[i] == '(') foundInternalTuple = (int) true;
                if (stringToParse[i] == ')' ) foundInternalTuple = (int) false;

                if (stringToParse[i] == '[') foundInternalArray = (int) true;
                if (stringToParse[i] == ']') foundInternalArray = (int) false;

                if (foundInternalTuple == (int) true || foundInternalArray == (int) true ) continue;

                if (stringToParse[i] == ',') *arraySize = *arraySize + 1;
            }

        } else {

            *arraySize = 1;
        }

    } else {

        *arraySize = 1;
    }

    // Debug function
    //printf("Number of strings (arraySize) = %d\n", *arraySize);

    // Allocate stringArray
    *stringArray = (char **) EG_alloc(*arraySize*sizeof(char **));
    if (*stringArray == NULL) {
        *arraySize = 0;
        return EGADS_MALLOC;
    }

    if (*arraySize == 1 && haveArray == (int) false) {

        // Debug function
        //printf("We don't have an array\n");

        noQuoteString = string_removeQuotation(stringToParse);

        (*stringArray)[0] = (char *) EG_alloc((strlen(noQuoteString)+1) * sizeof(char));

        // Check for malloc error
        if ((*stringArray)[0] == NULL) {

            // Free string array
            if (*stringArray != NULL) EG_free(*stringArray);
            *stringArray = NULL;

            *arraySize = 0;

            // Free no quote array
            if (noQuoteString != NULL) EG_free(noQuoteString);

            return EGADS_MALLOC;
        }

        strncpy((*stringArray)[0], noQuoteString, strlen(noQuoteString));
        (*stringArray)[0][strlen(noQuoteString)] = '\0';

        if (noQuoteString != NULL) {
            EG_free(noQuoteString);
            noQuoteString = NULL;
        }

    } else {

        foundInternalTuple = (int) false;
        foundInternalArray = (int) false;

        startIndex = 1;
        arrayIndex = 0;
        // Parse string based on defined pattern
        for (i = 1; i < strlen(stringToParse); i++) {

            if (stringToParse[i] == '(') foundInternalTuple = (int) true;
            if (stringToParse[i] == ')') foundInternalTuple = (int) false;

            if (stringToParse[i] == '[') foundInternalArray = (int) true;
            if (stringToParse[i] == ']') foundInternalArray = (int) false;

            if (foundInternalTuple == (int) true || foundInternalArray == (int) true ) continue;

            if(stringToParse[i] == ',' ||
               i == strlen(stringToParse)-1) {//stringToParse[i] == ')' ) {

                matchLength = i-startIndex;

                // Debug function
                //printf("MatchLength = %d\n", matchLength);

                // Make sure we aren't exceeding our arrays bounds
                if (arrayIndex >= *arraySize) {
                    return CAPS_MISMATCH;
                }

                // If not just a blank space
                if (matchLength > 0) {

                    // Allocate the quoted string array
                    quoteString = (char *) EG_alloc( (matchLength+1) * sizeof(char));
                    if (quoteString == NULL) return EGADS_MALLOC;

                    strncpy(quoteString, stringToParse + startIndex, matchLength);
                    quoteString[matchLength] = '\0';

                    // Remove quotations
                    noQuoteString = string_removeQuotation(quoteString);

                    // Free quote string array
                    if (quoteString != NULL) {
                        EG_free(quoteString);
                        quoteString = NULL;
                    }

                    // Allocate string array element based on no quote string
                    (*stringArray)[arrayIndex] = (char *) EG_alloc( (strlen(noQuoteString)+1) * sizeof(char));

                    // Check for malloc error
                    if ((*stringArray)[arrayIndex] == NULL) {

                        // Free string array
                        EG_free(quoteString);
                        for (j = 0; j <  arrayIndex; j++) EG_free((*stringArray)[j]);
                        EG_free(*stringArray);

                        *arraySize = 0;

                        // Free no quote string
                        if (noQuoteString != NULL) EG_free(noQuoteString);

                        return EGADS_MALLOC;
                    }

                    // Copy no quote string into array
                    strncpy((*stringArray)[arrayIndex], noQuoteString, strlen(noQuoteString));
                    (*stringArray)[arrayIndex][strlen(noQuoteString)] = '\0';

                    // Free no quote array
                    if (noQuoteString != NULL) {
                        EG_free(noQuoteString);
                        noQuoteString = NULL;
                    }

                    // Increment start indexes
                    arrayIndex += 1;
                    startIndex = i+1; // +1 to skip the commas
                }
            }
        }
    }

    // Debug function - print out number array
    //for (i = 0; i < *arraySize; i++) printf("Value = %s\n", (*stringArray)[i]);

    // Free quote array
    if (quoteString != NULL) EG_free(quoteString);

    // Free no quote array
    if (noQuoteString != NULL) EG_free(noQuoteString);


    return CAPS_SUCCESS;
}

// Simple json dictionary parser - currently doesn't support nested arrays for keyValue
int search_jsonDictionary(const char *stringToSearch, const char *keyWord, char **keyValue) {

    int  i; // Indexing for loop

    int  keyIndexStart; // Starting index of keyValue

    int  keyLength = 0; // Length of keyValue
    int  nesting = 0;

    char pattern[255]; // Substring - keyWord

    char *patternMatch = NULL; // Mathced substring - keyWord

    // Build pattern
    snprintf(pattern, 255, "\"%s\":", keyWord);

    // Debug function
    //printf("Pattern - [%s]\n", pattern);
    //printf("Search string - [%s]\n", stringToSearch);

    // Search string for pattern - return pointer to match
    patternMatch = strstr(stringToSearch, pattern);

    // If pattern is found - get keyValue
    if (patternMatch != NULL) {

        // Get the starting index of the keyValue (using pointer arithmetic)
        keyIndexStart = (patternMatch - stringToSearch) + strlen(pattern);

        if (stringToSearch[keyIndexStart] == ' ') keyIndexStart += 1; // Skip whitespace after :, some JSON writers have spaces others don't

        // See how long our keyValue is - at most length of incoming string
        for (i = 0; i < (int) strlen(stringToSearch); i++) {

            if (stringToSearch[keyIndexStart+keyLength] == '[') nesting++;

            // Array possibly nested
            if (stringToSearch[keyIndexStart] == '[') {
                keyLength += 1;

                if (stringToSearch[keyIndexStart+keyLength] == ']') {
                    keyLength += 1;
                    nesting--;
                    if (nesting == 0) break;
                }

            } else {

                keyLength += 1;
                if (stringToSearch[keyIndexStart+keyLength] == ',' ||
                    stringToSearch[keyIndexStart+keyLength] == '}') break;
            }
        }

        // Debug function
        //printf("Key start %d\n"    , keyIndexStart);
        //printf("Size of keyLength - %d\n", keyLength);

        // Free keyValue (if not already null) and allocate
        AIM_FREE(*keyValue);

        if (keyLength > 0) {
            *keyValue = (char *) EG_alloc((keyLength+1) * sizeof(char));
            if (*keyValue == NULL) return EGADS_MALLOC;

        } else {
            //printf("No match found for %s\n",keyWord);
            *keyValue = NULL; //"\0";
            return CAPS_NOTFOUND;
        }

        // Copy matched expression to keyValue
        strncpy((*keyValue), stringToSearch + keyIndexStart, keyLength);
        (*keyValue)[keyLength] = '\0';

    } else {

        //printf("No match found for %s\n",keyWord);
        *keyValue = NULL; //"\0";
        return CAPS_NOTFOUND;
    }

    // Debug function
    // printf("keyValue = [%s]\n", *keyValue);

    return CAPS_SUCCESS;
}

// Free an array strings
int string_freeArray(int numString, char **strings[]) {

    int i; // Indexing

    // Debugging
    // printf("NumStrings = %d\n", numString);

    if (numString == 0) return CAPS_SUCCESS;

    if ((*strings) != NULL) {

        // Free string one by one
        for (i = 0; i < numString; i++) {

            // Debugging
            // printf("Strings %d = %s\n", i, (*strings)[i]);

            if ((*strings)[i] != NULL) EG_free((*strings)[i]);
        }

        EG_free(*strings);
        *strings = NULL;
    }

    return CAPS_SUCCESS;
}

// Remove quotation marks around a string
char * string_removeQuotation(/*@null@*/ const char *string) {

    int startIndex = 0, endIndex = 0; // Indexing

    int length; // String length with quotes removed

    char *stringNoQuotation;

    if (string == NULL) return NULL;

    for (startIndex = 0; startIndex < strlen(string); startIndex++) {
        if (string[startIndex] != ' ') break;
    }

    for (endIndex = strlen(string)-1; endIndex == 0; endIndex--) {
        if (string[endIndex] != ' ') break;
    }

    // Debug
    //printf("Start Index = %d\n", startIndex);
    //printf("End Index = %d\n", endIndex);
    //printf("Incoming string - %s\n", string);

    // If string doesn't have quatation marks simply copy the string
    if ((string[startIndex] == '\"' || string[startIndex] == '\'') &&
        (string[endIndex]   == '\"' || string[endIndex]   == '\'')) {

        length = endIndex - startIndex -1;

        stringNoQuotation = (char *) EG_alloc((length + 1) * sizeof(char));
        if (stringNoQuotation == NULL) {
            printf("Memory allocation error in string_removeQuotation\n");
            return stringNoQuotation;
        }

        memcpy(stringNoQuotation, &string[startIndex+1], length);
        stringNoQuotation[length] = '\0';

    } else {
        stringNoQuotation = EG_strdup(string);
    }

    // Debug
    //printf("Outgoing string - %s\n", stringNoQuotation);

    return stringNoQuotation;
}

// Convert a string to double
int string_toDouble(const char *string, double *number) {

    char *next = NULL;

    if (string == NULL) return CAPS_NULLVALUE;

    *number = strtod(string, &next);

    if (string == next) return CAPS_BADVALUE;

    return CAPS_SUCCESS;
}

// Convert a string tuple of double with data-units to double in units
int string_toDoubleUnits(void *aimInfo, const char *string, const char *units, double *number) {

    int status = CAPS_SUCCESS;

    int numString = 0;
    char **stringArray = NULL; // Freeable
    char *dataunits=NULL;

    status = json_parseTuple(string, &numString, &stringArray);
    AIM_STATUS(aimInfo, status);

    if (numString != 2 || stringArray == NULL) {
        AIM_ERROR(aimInfo, "Expected tuple with number and units: %s", string);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    dataunits = stringArray[1];
    status = string_toDouble(stringArray[0], number);
    AIM_STATUS(aimInfo, status);

    status = aim_convert(aimInfo, 1, dataunits, number, units, number);
    AIM_STATUS(aimInfo, status);

cleanup:
    (void) string_freeArray(numString, &stringArray);

    return status;
}

// Convert a string to double array - array is assumed to be bounded by '[' and ']' and comma separated
//  for example [3.0, 41, -4.53E2]
int string_toDoubleArray(char *stringToSearch, int arraySize, double numberArray[]) {

    int status = CAPS_SUCCESS;
    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    char *numberString = NULL; // Temporary string to hold the found number

    // Debug function
    // printf("Search string - [%s]\n", stringToSearch);

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] != '[' &&
        stringToSearch[strlen(stringToSearch)] != ']') {

        printf("Error (string_toDoubleArray): incoming string should be bounded by '[' and ']' and comma separated\n");
        return CAPS_BADVALUE;
    }

    // Set default array values
    for (i = 0; i < arraySize; i++) numberArray[i] = 0;

    startIndex = 1;
    arrayIndex = 0;
    // Parse string based on defined pattern
    for (i = 1; i < strlen(stringToSearch); i++) {

        if(stringToSearch[i] == ',' ||
           stringToSearch[i] == ']') {

            matchLength = i-startIndex;

            // Debug function
            // printf("MatchLength = %d\n", matchLength);

            // Make sure we aren't exceeding our arrays bounds
            if (arrayIndex >= arraySize) {
                printf("Warning (string_toDoubleArray): Array size mismatch - to many values found (arraySize = %d, found = %d)!\n",
                        arraySize, arrayIndex);
                return CAPS_MISMATCH;
            }

            // If not just a blank space
            if (matchLength > 0) {

                if (numberString != NULL) EG_free(numberString);

                numberString = (char *) EG_alloc( (matchLength+1) * sizeof(char));

                // Check for malloc error
                if (numberString == NULL) {

                    // If no match was found fill in the rest of the array with 0
                    for (j = arrayIndex; j < arraySize; j++) numberArray[j] = 0;

                    printf("string_toDoubleArray: MALLOC error!\n");
                    return EGADS_MALLOC;
                }

                // Copy sub-string, that is the elements of the array
                strncpy(numberString, stringToSearch + startIndex, matchLength);
                numberString[matchLength] = '\0';

                //Debug function
                //printf("Number string = [%s]\n",numberString);

                // Convert array element to double
                status = string_toDouble(numberString, &numberArray[arrayIndex]);
                if (status != CAPS_SUCCESS) {
                    printf("Error: Cannot convert '%s' to double!\n", numberString);
                    AIM_FREE(numberString);
                    return CAPS_BADVALUE;
                }

                // Increment start indexes
                arrayIndex += 1;
                startIndex = i+1; // +1 to skip the commas
            }
        }
    }

    // Check to see if our array is under-sized
    if (arrayIndex < arraySize) {
        printf("Warning (string_toDoubleArray): Array size mismatch - remaining values will be 0\n");
    }

    // Debug function - print out number array
    //for (i = 0; i < arraySize; i++) printf("Value = %f\n", numberArray[i]);

    // Free allocated numberString
    AIM_FREE(numberString);

    return CAPS_SUCCESS;
}

// Convert a string tuple of array double with data-units to array double in units
int string_toDoubleArrayUnits(void *aimInfo, const char *string,
                              const char *units, int arraySize, double *numberArray) {

    int status = CAPS_SUCCESS;

    int numString = 0;
    char **stringArray = NULL; // Freeable
    char *dataunits=NULL;

    status = json_parseTuple(string, &numString, &stringArray);
    AIM_STATUS(aimInfo, status);

    if (numString != 2 || stringArray == NULL) {
        AIM_ERROR(aimInfo, "Expected tuple with number and units: %s", string);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    dataunits = stringArray[1];
    status = string_toDoubleArray(stringArray[0], arraySize, numberArray);
    AIM_STATUS(aimInfo, status);

    status = aim_convert(aimInfo, arraySize, dataunits, numberArray, units, numberArray);
    AIM_STATUS(aimInfo, status);

cleanup:
    (void) string_freeArray(numString, &stringArray);

    return status;
}

// Convert a string to an array of doubles - array is assumed to be bounded by '[' and ']' and comma separated
//  for example [3.0, 1.0, 9.0]. If the string coming in is not an array the an array of one is created.
int string_toDoubleDynamicArray(char *stringToSearch, int *arraySize, double *numberArray[]) {

    int status;       // Function return
    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    int haveArray = (int) false;
    char *numberString = NULL; // Temporary string to hold the found element of array

    // Debug function
    //printf("Search string - %s\n", stringToSearch);

    if (stringToSearch == NULL) return CAPS_BADVALUE;

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] != '[' &&
            stringToSearch[strlen(stringToSearch)] != ']') {
        *arraySize = 1;

    } else { // Lets count how many commas we have
        haveArray = (int) true;
        *arraySize = 1;
        for (i = 1; i < strlen(stringToSearch); i++) {
            if (stringToSearch[i] == ',') *arraySize = *arraySize + 1;
        }
    }

    // Debug function
    //printf("Number of strings = %d\n", *arraySize);

    // Free array if already set
    //if (*numberArray != NULL) EG_free(*numberArray);

    // Allocate numberArray
    *numberArray = (double *) EG_alloc(*arraySize*sizeof(double));
    if (*numberArray == NULL) {
        printf("string_toDoubleDynamicArray: MALLOC error!\n");
        return EGADS_MALLOC;
    }

    if (*arraySize == 1 && haveArray == (int) false) {

        // Debug function
        //printf("We don't have an array\n");

        status = string_toDouble(stringToSearch, &(*numberArray)[0]);
        if (status != CAPS_SUCCESS) return status;


    } else {

        startIndex = 1;
        arrayIndex = 0;
        // Parse string based on defined pattern
        for (i = 1; i < strlen(stringToSearch); i++) {

            if(stringToSearch[i] == ',' ||
                    stringToSearch[i] == ']') {

                matchLength = i-startIndex;

                // Debug function
                //printf("MatchLength = %d\n", matchLength);

                // Make sure we aren't exceeding our arrays bounds
                if (arrayIndex >= *arraySize) {
                    return CAPS_MISMATCH;
                }

                if (matchLength > 0) {

                    if (numberString != NULL) EG_free(numberString);

                    numberString = (char *) EG_alloc( (matchLength+1) * sizeof(char));

                    // Check for malloc error
                    if (numberString == NULL) {

                        // If we have a malloc error fill in the rest of the array with 0
                        for (j = arrayIndex; j < *arraySize; j++) (*numberArray)[j] = 0;

                        printf("string_toDoubleDynamicArray: MALLOC error!\n");
                        return EGADS_MALLOC;
                    }

                    // Copy sub-string, that is the elements of the array
                    strncpy(numberString, stringToSearch + startIndex, matchLength);
                    numberString[matchLength] = '\0';

                    //Debug function
                    //printf("Number string = [%s]\n",numberString);

                    // Convert array element to double
                    status = string_toDouble(numberString, &(*numberArray)[arrayIndex]);
                    if (status != CAPS_SUCCESS) {
                      AIM_FREE(numberString);
                      return status;
                    }

                    // Increment start indexes
                    arrayIndex += 1;
                    startIndex = i+1; // +1 to skip the commas
                }
            }
        }
    }

    // Debug function - print out number array
    //for (i = 0; i < arraySize; i++) printf("Value = %f\n", numberArray[i]);

    // Free allocated numberString
    if (numberString != NULL) EG_free(numberString);

    return CAPS_SUCCESS;
}


// Convert a string tuple of array double with data-units to array double in units
int string_toDoubleDynamicArrayUnits(void *aimInfo, const char *string,
                                     const char *units, int *arraySize, double **numberArray) {

    int status = CAPS_SUCCESS;

    int numString = 0;
    char **stringArray = NULL; // Freeable
    char *dataunits=NULL;

    status = json_parseTuple(string, &numString, &stringArray);
    AIM_STATUS(aimInfo, status);

    if (numString != 2 || stringArray == NULL) {
        AIM_ERROR(aimInfo, "Expected tuple with number and units: %s", string);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    dataunits = stringArray[1];
    status = string_toDoubleDynamicArray(stringArray[0], arraySize, numberArray);
    AIM_STATUS(aimInfo, status);

    status = aim_convert(aimInfo, *arraySize, dataunits, *numberArray, units, *numberArray);
    AIM_STATUS(aimInfo, status);

cleanup:
    (void) string_freeArray(numString, &stringArray);

    return status;
}


// Convert a string to an array of strings - array is assumed to be bounded by '[' and ']' and comma separated
//  for example ["3.0", "hey", "foo"] - note if strings contain quotation marks, they should NOT contain commas.
int string_toStringArray(char *stringToSearch, int arraySize, char *stringArray[]) {

    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    // int haveArray = (int) false;
    char *quoteString = NULL; // Temporary string to hold the found string
    char *noQuoteString = NULL;  // Temporary string to hold the found string with quotation marks removed
    int insideQuotes = (int) false;

    // Debug function
    //printf("Search string - %s length %d\n", stringToSearch, strlen(stringToSearch));

    if (stringToSearch == NULL) return CAPS_BADVALUE;

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] != '[' &&
        stringToSearch[strlen(stringToSearch)] != ']') {

        printf("Error (string_toStringArray): incoming string should be bounded by '[' and ']' and comma separated\n");
        return CAPS_BADVALUE;
    }

    // Set default array values
    for (i = 0; i < arraySize; i++) stringArray[i] = NULL;

    startIndex = 1;
    arrayIndex = 0;
    // Parse string based on defined pattern
    for (i = 1; i < strlen(stringToSearch); i++) {

        if (stringToSearch[i] == '"') {
            insideQuotes = !insideQuotes;
        }

        if((stringToSearch[i] == ',' && (!insideQuotes)) ||
                stringToSearch[i] == ']') {
            matchLength = i-startIndex;

            // Debug function
            //printf("MatchLength = %d\n", matchLength);

            // Make sure we aren't exceeding our arrays bounds
            if (arrayIndex >= arraySize) {
                return CAPS_MISMATCH;
            }

            // If not just a blank space
            if (matchLength > 0) {

                // Allocate the quoted string array
                quoteString = (char *) EG_alloc( (matchLength+1) * sizeof(char));
                if (quoteString == NULL) return EGADS_MALLOC;

                strncpy(quoteString, stringToSearch + startIndex, matchLength);
                quoteString[matchLength] = '\0';

                // Remove quotations
                noQuoteString = string_removeQuotation(quoteString);

                // Free quote string array
                if (quoteString != NULL) {
                    EG_free(quoteString);
                    quoteString = NULL;
                }

                // Allocate string array element based on no quote string
                stringArray[arrayIndex] = (char *) EG_alloc( (strlen(noQuoteString)+1) * sizeof(char));

                // Check for malloc error
                if (stringArray[arrayIndex] == NULL) {
                    EG_free(quoteString);
                    for (j = 0; j <  arrayIndex; j++) EG_free(stringArray[j]);
                    EG_free(noQuoteString);
                    return EGADS_MALLOC;
                }

                // Copy no quote string into array
                strncpy(stringArray[arrayIndex], noQuoteString, strlen(noQuoteString));
                stringArray[arrayIndex][strlen(noQuoteString)] = '\0';

                // Free no quote array
                if (noQuoteString != NULL) {
                    EG_free(noQuoteString);
                    noQuoteString = NULL;
                }

                // Increment start indexes
                arrayIndex += 1;
                startIndex = i+1; // +1 to skip the commas
            }
        }
    }

    // Debug function - print out number array
    //for (i = 0; i < *arraySize; i++) printf("Value = %s\n", (*stringArray)[i]);

    EG_free(quoteString); // Free quote array
    EG_free(noQuoteString); // Free no quote array

    return CAPS_SUCCESS;
}

// Convert a string to an array of strings - array is assumed to be bounded by '[' and ']' and comma separated
//  for example ["3.0", "hey", "foo"] - note if strings contain quotation marks, they should NOT contain commas. If the
//  string coming in is not an array the string is simply copied. Also quotations are removed
int string_toStringDynamicArray(char *stringToSearch, int *arraySize, char **stringArray[]) {

    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex, nesting=0, length;
    int haveArray = (int) false;
    char *quoteString = NULL; // Temporary string to hold the found string
    int insideQuotes = (int) false;

    // Debug function
    //printf("Search string - %s length %d\n", stringToSearch, strlen(stringToSearch));

    if (stringToSearch == NULL) return CAPS_BADVALUE;

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] == '[' ) {

        // Debug function
        //printf("[ found\n");
        length = strlen(stringToSearch);
        if(stringToSearch[length-1] == ']') { // Lets count how many commas we have

            // Debug function
            //printf("] found\n");

            haveArray = (int) true;
            *arraySize = 1;
            for (i = 1; i < length; i++) {

                // Array possibly nested
                if (stringToSearch[i] == '[' ||
                    stringToSearch[i] == '{') {
                    nesting++;

                } else if (stringToSearch[i] == ']' ||
                           stringToSearch[i] == '}') {
                    nesting--;

                } else if (stringToSearch[i] == '"') {
                    insideQuotes = !insideQuotes;
                }

                if (stringToSearch[i] == ',' && (!insideQuotes) && nesting <= 0) {
                    *arraySize = *arraySize + 1;
                }
            }

        } else {
            *arraySize = 1;
        }

    } else {
        *arraySize = 1;
    }

    // Debug function
    //printf("Number of strings (arraySize) = %d\n", *arraySize);

    // Allocate stringArray
    *stringArray = (char **) EG_alloc(*arraySize*sizeof(char **));
    if (*stringArray == NULL) return EGADS_MALLOC;
    for (i = 0; i < *arraySize; i++) (*stringArray)[i] = NULL;

    if (*arraySize == 1 && haveArray == (int) false) {

        // Debug function
        //printf("We don't have an array\n");

        (*stringArray)[0] = string_removeQuotation(stringToSearch);

        // Check for malloc error
        if ((*stringArray)[0] == NULL) {
            // Free string array
            AIM_FREE(*stringArray);
            return EGADS_MALLOC;
        }

    } else {

        insideQuotes = (int) false;

        startIndex = 1;
        arrayIndex = 0;
        nesting = 0;
        length = strlen(stringToSearch);
        // Parse string based on defined pattern
        for (i = 1; i < length; i++) {

            // Array possibly nested
            if (stringToSearch[i] == '[' ||
                stringToSearch[i] == '{') {
                nesting++;

            } else if (stringToSearch[i] == ']' ||
                       stringToSearch[i] == '}') {
                nesting--;

            } else if (stringToSearch[i] == '"') {
                insideQuotes = !insideQuotes;
            }

            if ((stringToSearch[i] == ',' && (!insideQuotes) && nesting <= 0) ||
                i == length-1) {

                while (stringToSearch[startIndex] == ' ') startIndex++;
                matchLength = i-startIndex;

                // Debug function
                //printf("MatchLength = %d\n", matchLength);

                // Make sure we aren't exceeding our arrays bounds
                if (arrayIndex >= *arraySize) {
                    return CAPS_MISMATCH;
                }

                // If not just a blank space
                if (matchLength > 0) {

                    // Allocate the quoted string array
                    quoteString = (char *) EG_alloc( (matchLength+1) * sizeof(char));
                    if (quoteString == NULL) return EGADS_MALLOC;

                    strncpy(quoteString, stringToSearch + startIndex, matchLength);
                    quoteString[matchLength] = '\0';

                    // Remove quotations
                    (*stringArray)[arrayIndex] = string_removeQuotation(quoteString);

                    // Free quote string array
                    AIM_FREE(quoteString);

                    // Check for malloc error
                    if ((*stringArray)[arrayIndex] == NULL) {
                        for (j = 0; j <  arrayIndex; j++) AIM_FREE((*stringArray)[j]);
                        AIM_FREE(*stringArray);
                        return EGADS_MALLOC;
                    }

                    // Increment start indexes
                    arrayIndex += 1;
                    startIndex = i+1; // +1 to skip the commas
                }
            }
        }
    }

    // Debug function - print out number array
    //for (i = 0; i < *arraySize; i++) printf("Value = %s\n", (*stringArray)[i]);

    // Free quote array
    AIM_FREE(quoteString);


    return CAPS_SUCCESS;
}

// Convert a string to integer
int string_toInteger(const char *string, int *number) {

    char *next = NULL;

    if (string == NULL) return CAPS_NULLVALUE;

    *number = strtol(string, &next, 10);

    if (string == next) return CAPS_BADVALUE;

    return CAPS_SUCCESS;
}

// Convert a string to integer array - array is assumed to be bounded by '[' and ']' and comma separated
//  for example [3, 41, -4]
int string_toIntegerArray(char *stringToSearch, int arraySize, int numberArray[]) {

    int status = CAPS_SUCCESS;
    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    char *numberString = NULL; // Temporary string to hold the found number

    // Debug function
    // printf("Search string - [%s]\n", stringToSearch);

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] != '[' &&
        stringToSearch[strlen(stringToSearch)] != ']') {

        printf("Error (string_toIntegerArray): incoming string should be bounded by '[' and ']' and comma separated\n");
        return CAPS_BADVALUE;
    }

    // Set default array values
    for (i = 0; i < arraySize; i++) numberArray[i] = 0;

    startIndex = 1;
    arrayIndex = 0;
    // Parse string based on defined pattern
    for (i = 1; i < strlen(stringToSearch); i++) {

        if(stringToSearch[i] == ',' ||
           stringToSearch[i] == ']') {

            matchLength = i-startIndex;

            // Debug function
            // printf("MatchLength = %d\n", matchLength);

            // Make sure we aren't exceeding our arrays bounds
            if (arrayIndex >= arraySize) {
                printf("Warning (string_toIntegerArray): Array size mismatch - to many values found!!\n");
                return CAPS_MISMATCH;
            }

            // If not just a blank space
            if (matchLength > 0) {

                if (numberString != NULL) EG_free(numberString);

                numberString = (char *) EG_alloc( (matchLength+1) * sizeof(char));

                // Check for malloc error
                if (numberString == NULL) {

                    // If no match was found fill in the rest of the array with 0
                    for (j = arrayIndex; j < arraySize; j++) numberArray[j] = 0;

                    printf("string_toIntegerArray: MALLOC error!\n");
                    return EGADS_MALLOC;
                }

                // Copy sub-string, that is the elements of the array
                strncpy(numberString, stringToSearch + startIndex, matchLength);
                numberString[matchLength] = '\0';

                //Debug function
                //printf("Number string = [%s]\n",numberString);

                // Convert array element to double
                status = string_toInteger(numberString, &numberArray[arrayIndex]);
                if (status != CAPS_SUCCESS) {
                    printf("Error: Cannot convert '%s' to integer!\n", numberString);
                    AIM_FREE(numberString);
                    return CAPS_BADVALUE;
                }

                // Increment start indexes
                arrayIndex += 1;
                startIndex = i+1; // +1 to skip the commas
            }
        }
    }

    // Check to see if our array is under-sized
    if (arrayIndex != arraySize) {
        printf("Warning (string_toIntegerArray): Array size mismatch - remaining values will be 0\n");
        AIM_FREE(numberString);
        return CAPS_MISMATCH;
    }

    // Debug function - print out number array
    //for (i = 0; i < arraySize; i++) printf("Value = %f\n", numberArray[i]);

    // Free allocated numberString
    if (numberString != NULL) EG_free(numberString);

    return CAPS_SUCCESS;
}

// Convert a string to an array of doubles - array is assumed to be bounded by '[' and ']' and comma separated
//  for example [3, 1, 9]. If the string coming in is not an array the an array of one is created.
int string_toIntegerDynamicArray(char *stringToSearch, int *arraySize, int *numberArray[]) {

    int status;       // Function return
    int i, j;         // Array indexing

    int matchLength;  // String length of matching pattern
    int startIndex; // Keep track of where we are in the string array

    int arrayIndex;
    int haveArray = (int) false;
    char *numberString = NULL; // Temporary string to hold the found element of array

    // Debug function
    //printf("Search string - %s\n", stringToSearch);

    if (stringToSearch == NULL) return CAPS_BADVALUE;

    // Check to make sure first and last of the incoming string denote an array
    if (stringToSearch[0] != '[' &&
            stringToSearch[strlen(stringToSearch)] != ']') {
        *arraySize = 1;

    } else { // Lets count how many commas we have
        haveArray = (int) true;
        *arraySize = 1;
        for (i = 1; i < strlen(stringToSearch); i++) {
            if (stringToSearch[i] == ',') *arraySize = *arraySize + 1;
        }
    }

    // Debug function
    //printf("Number of strings = %d\n", *arraySize);

    // Free array if already set
    //if (*numberArray != NULL) EG_free(*numberArray);

    // Allocate numberArray
    *numberArray = (int *) EG_alloc(*arraySize*sizeof(int));
    if (*numberArray == NULL) {
        printf("string_toIntegerDynamicArray: MALLOC error!\n");
        return EGADS_MALLOC;
    }

    if (*arraySize == 1 && haveArray == (int) false) {

        // Debug function
        //printf("We don't have an array\n");

        status = string_toInteger(stringToSearch, &(*numberArray)[0]);
        if (status != CAPS_SUCCESS) return status;


    } else {

        startIndex = 1;
        arrayIndex = 0;
        // Parse string based on defined pattern
        for (i = 1; i < strlen(stringToSearch); i++) {

            if(stringToSearch[i] == ',' ||
                    stringToSearch[i] == ']') {

                matchLength = i-startIndex;

                // Debug function
                //printf("MatchLength = %d\n", matchLength);

                // Make sure we aren't exceeding our arrays bounds
                if (arrayIndex >= *arraySize) {
                    return CAPS_MISMATCH;
                }

                if (matchLength > 0) {

                    if (numberString != NULL) EG_free(numberString);

                    numberString = (char *) EG_alloc( (matchLength+1) * sizeof(char));

                    // Check for malloc error
                    if (numberString == NULL) {

                        // If we have a malloc error fill in the rest of the array with 0
                        for (j = arrayIndex; j < *arraySize; j++) (*numberArray)[j] = 0;

                        printf("string_toIntegerDynamicArray: MALLOC error!\n");
                        return EGADS_MALLOC;
                    }

                    // Copy sub-string, that is the elements of the array
                    strncpy(numberString, stringToSearch + startIndex, matchLength);
                    numberString[matchLength] = '\0';

                    //Debug function
                    //printf("Number string = [%s]\n",numberString);

                    // Convert array element to double
                    status = string_toInteger(numberString, &(*numberArray)[arrayIndex]);
                    if (status != CAPS_SUCCESS) {
                      AIM_FREE(numberString);
                      return status;
                    }

                    // Increment start indexes
                    arrayIndex += 1;
                    startIndex = i+1; // +1 to skip the commas
                }
            }
        }
    }

    // Debug function - print out number array
    //for (i = 0; i < arraySize; i++) printf("Value = %f\n", numberArray[i]);

    // Free allocated numberString
    AIM_FREE(numberString);

    return CAPS_SUCCESS;
}

// Convert a string to boolean
int string_toBoolean(char *string, int *number) {

    if (strncmp(string, "T", 1) == 0 ||
        strncmp(string, "t", 1) == 0 ||
        strncmp(string, "1", 1) == 0) {

        *number = (int) true;

    } else if (strncmp(string, "F", 1) == 0 ||
               strncmp(string, "f", 1) == 0 ||
               strncmp(string, "0", 1) == 0) {

        *number = (int) false;

    } else {

        printf("Error: Unrecognized boolean string - %s\n", string);
        return CAPS_BADVALUE;
    }

    return CAPS_SUCCESS;
}

// Convert a single string to a prog_argv array
int string_toProgArgs(char *meshInputString, int *prog_argc, char ***prog_argv)
{
    // Input:
    //        meshInputstring - (single) string as would be input at the command line
    // Output:
    //        prog_argc - pointer to number of command line inputs
    //        prog_argv - pointer to array of strings
    // Note: For a command line input the first argument is the program call itself - space is allocated
    //             to account for this but is isn't used.

    int startIndex, endIndex, stringIndex; //Indexing;
    int stringLength;

    //printf("Mesh inputString = %s\n", meshInputString);

    if (meshInputString != NULL) {

        // Count out how many space for allocation purposes
        *prog_argc = 2; // Argue 1 = program, Last argument will not have a space after it
        // *prog_argc = 1; // Last argument will not have a space after it

        for (startIndex = 0; startIndex < strlen(meshInputString); startIndex++) {
            if (meshInputString[startIndex] == ' ') *prog_argc += 1;
        }

        // Size argv array based on number of found spaces
        *prog_argv = (char **) EG_alloc(*prog_argc*sizeof(char *));
        if (*prog_argv == NULL) return EGADS_MALLOC;

        for (stringIndex = 0; stringIndex < *prog_argc; stringIndex++) (*prog_argv)[stringIndex] = NULL;

        // Loop through input string again and break at spaces
        startIndex = 0; // Index of start of string

        stringIndex = 1; // First index [0] is nothing

        for (endIndex = 0; endIndex < strlen(meshInputString); endIndex++) {

            if (meshInputString[endIndex] == ' ') {
                stringLength = endIndex-startIndex;

                (*prog_argv)[stringIndex] = (char *) EG_alloc((stringLength+1)*sizeof(char)); // Allocate array
                if ((*prog_argv)[stringIndex] == NULL) return EGADS_MALLOC;

                memcpy((*prog_argv)[stringIndex],
                        meshInputString+startIndex,
                        stringLength*sizeof(char));

                (*prog_argv)[stringIndex][stringLength] = '\0';

                stringIndex += 1; // Increment prog_argv index
                startIndex = endIndex+1; //Skip space, ' '
            }
        }

        // Get last value
        stringLength = endIndex-startIndex;

        (*prog_argv)[stringIndex] = (char *) EG_alloc((stringLength+1)*sizeof(char)); // Allocate array
        if ((*prog_argv)[stringIndex] == NULL) return EGADS_MALLOC;

        memcpy((*prog_argv)[stringIndex],
                meshInputString+startIndex,
                stringLength*sizeof(char));

        (*prog_argv)[stringIndex][stringLength] = '\0';

    } else {
        *prog_argc = 0;
        *prog_argv = NULL;
    }

    //printf("Number of args = %d\n",*prog_argc);
    //for (i = 0; i <*prog_argc ; i++) printf("Arg %d = %s\n", i, (*prog_argv)[i]);

    return CAPS_SUCCESS;

}

// Force a string to upper case
void string_toUpperCase ( char *sPtr )
{
    while ( *sPtr != '\0' ) {
        *sPtr = toupper ( ( unsigned char ) *sPtr );
        ++sPtr;
    }
}

// Return newly allocated string formatted by `format` and any variadic string args
// NOTE: sentinel denoting end of variadic args is NULL, if NULL is not provided,
//       behavior may be undefined
char *string_format(char *format, ...) {

    int i, narg, nspec, length;
    char c, *formatted = NULL, *arg = NULL;
    va_list args;

    // get number of expected format string args (number of %s specifiers)
    for (i = 0, nspec = 0; i < strlen(format) - 1; i++) {
        c = format[i];
        if (c == '%') {
            if ((i > 0) && format[i-1] == '\\') {
                continue;
            }
            if (format[i+1] == 's') {
                nspec++;
            }
        }
    }

    // get total string length
    va_start(args, format);
    length = strlen(format);
    narg = 0;
    arg = va_arg(args, char *);
    //   loop until arg is NULL, or narg == nspec, this safeguards a little
    //   against undefined behavior if NULL sentinel is not used
    while (arg != NULL && narg < nspec) {
        length += strlen(arg);
        arg = va_arg(args, char *);
        narg++;
    }
    va_end(args);

    // allocate
    formatted = EG_alloc(sizeof(char) * (length + 1));
    if (formatted == NULL) {
        return NULL;
    }

    // create formatted string
    va_start(args, format);
    vsprintf(formatted, format, args);
    va_end(args);

    return formatted;
}

// Return whether string `find` is in `array`
int string_isInArray(char *find, int arraySize, char **array) {
    int i;
    for (i = 0; i < arraySize; i++) {
        if (strcmp(array[i], find) == 0) {
            return (int) true;
        }
    }
    return (int) false;
}

// The max x,y,z coordinates where P(3*i + 0) = x_i, P(3*i + 1) = y_i, and P(3*i + 2) = z_i
void maxCoords(int sizeP, double *P, double *x, double *y, double *z) {

    //printf("Getting max bounds\n");
    int i;
    *x = 0;
    *y = 0;
    *z = 0;
    for (i = 0; i < sizeP ; i++){
        if (P[3*i  ] > *x) *x = P[3*i  ];
        if (P[3*i+1] > *y) *y = P[3*i+1];
        if (P[3*i+2] > *z) *z = P[3*i+2];
    }
}

// The min x,y,z coordinates where P(3*i + 0) = x_i, P(3*i + 1) = y_i, and P(3*i + 2) = z_i
void minCoords(int sizeP, double *P, double *x, double *y, double *z) {

    //printf("Getting min bounds\n");
    int i;
    *x = 0;
    *y = 0;
    *z = 0;
    for (i = 0; i < sizeP ; i++){
        if (P[3*i  ] < *x) *x = P[3*i  ];
        if (P[3*i+1] < *y) *y = P[3*i+1];
        if (P[3*i+2] < *z) *z = P[3*i+2];
    }
}

// Return the endianness of the machine we are working on
int get_MachineENDIANNESS(void) {
    // Return = 0 - little
    // Return = 1 - big

    int num = 1;

    if (*(char *)&num == 1) {

        //printf("\nLittle-Endian\n");
        return 0;

    } else {

        //printf("Big-Endian\n");
        return 1;
    }
}

#ifdef WIN32
int getline(char **linep, size_t *linecapp, FILE *stream)
{
    int  i = 0, cnt = 0;

    char *buffer;

    buffer = (char *) EG_reall(*linep, 1026*sizeof(char));
    if (buffer == NULL) return EGADS_MALLOC;

    *linep    = buffer;
    *linecapp = 1024;
    while (cnt < 1024 && !ferror(stream)) {
        i = fread(&buffer[cnt], sizeof(char), 1, stream);
        if (i == 0) break;
        cnt++;
        if ((buffer[cnt-1] == 10) || (buffer[cnt-1] == 13) || // 10 = new line, 13 = carriage return, 0 = NULL
            (buffer[cnt-1] == 0)) break;
    }

    buffer[cnt] = 0;
    return cnt == 0 ? -1 : cnt;
}
#endif

// Return the max value of two doubles
double max_DoubleVal(double x1, double x2)
{
    if (x1 >= x2) return x1;
    else          return x2;
}
// Return the min value of two doubles
double min_DoubleVal(double x1, double x2)
{
    if (x1 <= x2) return x1;
    else          return x2;
}

// Cross product axb = c
void cross_DoubleVal(double a[], double b[], double c[]){

    c[0] = (a[1]*b[2]) - (a[2]*b[1]);
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);
    c[2] = (a[0]*b[1]) - (a[1]*b[0]);
}

// Dot product a*b = c
double dot_DoubleVal(double a[], double b[]){
    double c;

    c = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);

    return c;
}

// Distance between two points sqrt(dot(a-b,a-b))
double dist_DoubleVal(double a[], double b[]){
    double d[3];

    d[0] = a[0]-b[0];
    d[1] = a[1]-b[1];
    d[2] = a[2]-b[2];

    return sqrt(dot_DoubleVal(d,d));
}

// Convert an integer to a string of a given field width and justification
char * convert_integerToString(int integerVal, int fieldWidth, int leftOrRight)
{

    // Input:
    //  integerVal  - Integer to convert to string
    //  fieldWidth  - Desired fieldwidth of string must be at least 1
    //  leftofRight - Left justified = 0, Right justified = anything else

    // Ouput:
    //  stringVal - Returned integer in string format. Returned as a const char *

    int inputTest = 0; // Input check

    char tmp[42];
    char *stringVal = NULL; // Returned string array


    // First check input parameters
    if( fieldWidth > 15 ) {
        inputTest = 2;
    }

    if (fieldWidth <= 0) inputTest = 3;

    // If checks failed report error and return an empty string
    if (inputTest > 0) {

        if (inputTest == 3) {
            printf("Error in convert_integerToString: fieldWidth <= 0 \n");
        } else if (inputTest == 2) {
            printf("Error in convert_integerToString: fieldWidth > 15 not verified\n");
        } else {
            printf("Error in convert_integerToString: Input %i is too large for the requested fieldWidth of %d\n", integerVal,
                                                                                                                   fieldWidth);
        }

        printf("\tReturning a 'NaN' string.\n");
        stringVal = EG_strdup("NaN");
        return stringVal;
    }


    // Populate output string array with blank spaces and integer depending on justification
    if (leftOrRight == 0) {

        inputTest = sprintf(tmp, "%-*d", fieldWidth, integerVal);

    } else {

        inputTest = sprintf(tmp, "%*d", fieldWidth, integerVal);
    }

    if (inputTest < 0 || inputTest > fieldWidth) {
        printf("Error in convert_integerToString: Input %i fails with requested fieldWidth of %d\n", integerVal,
                                                                                                     fieldWidth);
        printf("\tReturning a 'NaN' string.\n");
        stringVal = EG_strdup("NaN");
        return stringVal;
    }

    // Copy over the output
    stringVal = EG_strdup(tmp);

    // Return string array
    return stringVal;
}

// Convert an double to a string (scientific notation is used depending on the fieldwidth and value) of a given field width
char * convert_doubleToString(double doubleVal, int fieldWidth, int leftOrRight)
{

    // Input:
    //  integerVal  - Double to convert to string
    //  fieldWidth  - Desired fieldwidth of string must be at least 5
    //  leftofRight - Left justified = 0, Right justified = anything else

    // Ouput:
    //  stringVal - Returned integer in string format. Returned as a char * (freeable)

    int i; // Index

    int inputTest = 0; // Input check
    int scival, offset;
    int minfieldWidth = 0; // Minimum field width for input

    //int powerExpUpper, powerExpLower; // Exponent powers for max width number comparison

    char *stringVal = NULL;
    const char *nan = "NaN";
    char numString[255];
    char sci[10], tmp[42];

    // If checks inputs and return an empty string
    if (fabs(doubleVal) < 1 && fabs(doubleVal) != 0 && doubleVal < 0  && fieldWidth < 3 ) {
        minfieldWidth = 2;
        inputTest = 1;

    } else if ( fabs(doubleVal) >= 1  && doubleVal < 0 && fieldWidth < 3 ) {
        minfieldWidth = 2;
        inputTest = 1;

    } else if ( fabs(doubleVal) < 1 && fabs(doubleVal) != 0 && doubleVal > 0  && fieldWidth < 3 ) {
        minfieldWidth = 2;
        inputTest = 1;

    } else if ( fabs(doubleVal) >= 1  && doubleVal > 0 && fieldWidth < 2 ) {
        minfieldWidth = 1;
        inputTest = 1;

    } else if ( fieldWidth < 2 ) {
        minfieldWidth = 1;
        inputTest = 1;
    }

    if (inputTest == 1) {
        printf("Error in convert_doubleToString: Input fieldWidth of %d must be greater than %d for the input value %E\n", fieldWidth, minfieldWidth,doubleVal);
        printf("\tReturning a 'NaN' string.\n");
        stringVal = EG_strdup(nan);
        return stringVal;
    }

    // allocate the string
    stringVal = (char *) EG_alloc((fieldWidth+1)*sizeof(char));

    // If zero, and yes the sign matters!?!
    if (doubleVal == 0.0 || doubleVal == +0.0 || doubleVal == -0.0) {
        sprintf(stringVal, "%#.*f", fieldWidth-2, 0.0);
        return stringVal;
    }

    offset = 2;                  // the period and the 'E'
    if (doubleVal < 0) offset++; // account for the negative sign

    // extract the exponent with rounding to check for simple float formatting
    sprintf(tmp, "%1.*E", fieldWidth-offset, doubleVal);
    i = fieldWidth;
    tmp[i++] = '\0';
    scival = atoi(tmp+i);

    // check if a floating point number is more appropriate than scientific
    // sized to make sure the '.' is included in the string (otherwise Fortran isn't happy)
    if (scival > -2 && scival < fieldWidth-(offset-2)-1) {

        sprintf(numString, "%#*.*f", fieldWidth, fieldWidth-MAX(scival,0)-offset, doubleVal);

    } else {

        // do loop as the exponent might change due to rounding
        sprintf(sci, "%+d", scival);
        do {
            // check if the scientific number fits in the fieldWidth
            if ((int)(fieldWidth-offset-strlen(sci)-1) < 0) {
                printf("Error in convert_doubleToString: Cannot write %E with field with %d!\n", doubleVal, fieldWidth);
                printf("\tReturning a 'NaN' string.\n");
                EG_free(stringVal);
                stringVal = EG_strdup(nan);
                return stringVal;
            }

            // construct the format statement based on the available digits
            sprintf(tmp, "%#1.*E", (int)(fieldWidth-offset-strlen(sci)-1), doubleVal);

            // truncate at 'E'
            i = fieldWidth-strlen(sci)-1;
            tmp[i++] = '\0';

            // get the exponent again and minimize it's size
            scival = atoi(tmp+i);
            sprintf(sci, "%+d", scival);

            // print the final string with the exponent
            sprintf(numString, "%sE%s", tmp, sci);

        } while ( strlen(numString) != fieldWidth );
    }

    // Populate output string array with blank spaces and number string depending on justification
    // Fairly confident this is mute as the above padds everything with 0's
    if (leftOrRight == 0) {

        strncpy(stringVal, numString, fieldWidth);
        for (i = 0; i < fieldWidth - (int) strlen(numString); i++) {
            stringVal[strlen(numString)+i] = ' ';
        }

    } else {

        if (fieldWidth - (int) strlen(numString) > 0 ) {
            for (i = 0; i < fieldWidth - (int) strlen(numString); i++) {
                stringVal[i] = ' ';
            }

            for (i = 0; i < (int) strlen(numString); i++) {
                stringVal[fieldWidth - (int) strlen(numString) + i] = numString[i];
            }

        } else {

            strncpy(stringVal, numString, fieldWidth);
        }
    }

    // Add termination character to the end of the string
    stringVal[fieldWidth] = '\0';

    return stringVal;
}

// Factorizes in place the square linear system A x = b using simple LU decomposition
// Returns CAPS_BADVALUE for a singular matrix
int factorLU(int n, double A[] )
{
    int i,j,k;
    double y;

    // Compute the LU decomposition in place
    for(k = 0; k < n-1; k++) {
        if (A[k*n+k] == 0) return CAPS_BADVALUE;
        for(j = k+1; j < n; j++) {
            y = A[j*n+k]/A[k*n+k];
            for(i = k; i < n; i++) {
                A[j*n+i] = A[j*n+i]-y*A[k*n+i];
            }
            A[j*n+k] = y;
        }
    }

    return CAPS_SUCCESS;
}

// Solves the factorized square linear system LU x = b
int backsolveLU(int n, double LU[], double b[], double x[] )
{
    int i,j;
    double y;

    // Forward solve
    for(i = 0; i < n; i++) {
        y=0.0;
        for(j = 0 ;j < i;j++) {
            y += LU[i*n+j]*x[j];
        }
        x[i]=(b[i]-y);
    }

    // Back substitution
    for(i = n-1; i >=0; i--) {
        y = 0.0;
        for(j = i+1; j < n; j++) {
            y += LU[i*n+j]*x[j];
        }
        x[i] = (x[i]-y)/LU[i*n+i];
    }

    return CAPS_SUCCESS;
}

// Solves the square linear system A x = b using simple LU decomposition
// Returns CAPS_BADVALUE for a singular matrix
int solveLU(int n, double A[], double b[], double x[] )
{
    int status = CAPS_SUCCESS;

    status = factorLU(n, A);
    if (status != CAPS_SUCCESS) return status;
    backsolveLU(n, A, b, x);

    return CAPS_SUCCESS;
}

// Prints all attributes on an ego
int print_AllAttr( void *aimInfo, ego obj )
{
#define NBUF 512
    int          status;
    int          i, j, nattr, atype, alen;
    char         buffer[NBUF], tmp[NBUF];
    const int    *pints;
    const char   *name, *pstr;
    const double *preals;

    AIM_ADDLINE(aimInfo, "Available attributes are:");

    nattr = 0;
    status  = EG_attributeNum(obj, &nattr);
    if ((status == EGADS_SUCCESS) && (nattr != 0)) {
        for (i = 1; i <= nattr; i++) {
            status = EG_attributeGet(obj, i, &name, &atype, &alen,
                                   &pints, &preals, &pstr);
            if (status != EGADS_SUCCESS) continue;
            snprintf(buffer, NBUF, "   %s:", name);
            if (atype == ATTRINT) {
                for (j = 0; j < alen; j++) {
                  snprintf(tmp, NBUF, "%s %d", buffer, pints[j]);
                  strcpy(buffer, tmp);
                }
            } else if (atype == ATTRREAL) {
                for (j = 0; j < alen; j++) {
                  snprintf(tmp, NBUF, "%s %lf", buffer, preals[j]);
                  strcpy(buffer, tmp);
                }
            } else if (atype == ATTRSTRING) {
                snprintf(tmp, NBUF, "%s %s", buffer, pstr);
                strcpy(buffer, tmp);
            } else if (atype == ATTRCSYS) {
                snprintf(tmp, NBUF, "%s csys", buffer);
                strcpy(buffer, tmp);
                for (j = 0; j < alen; j++) {
                  snprintf(tmp, NBUF, "%s %lf", buffer, preals[j]);
                  strcpy(buffer, tmp);
                }
            } else if (atype == ATTRPTR) {
                snprintf(tmp, NBUF, "%s pointer", buffer);
                strcpy(buffer, tmp);
            } else {
                snprintf(tmp, NBUF, "%s unknown attribute type!", buffer);
                strcpy(buffer, tmp);
            }
            AIM_ADDLINE(aimInfo, buffer);
        }
    }

    status = CAPS_SUCCESS;

    return status;
#undef NBUF
}

// Search a mapAttrToIndex structure for a given keyword and set/return the corresponding index
int get_mapAttrToIndexIndex(mapAttrToIndexStruct *attrMap, const char *keyWord, int *index) {

    // If the keyword is not found a CAPS_NOTFOUND is returned

    int i; // Indexing

    *index = CAPSMAGIC; // Default value of return

    if (attrMap == NULL) return CAPS_NULLVALUE;
    if (keyWord == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < attrMap->numAttribute; i++) {

        if (strcmp(attrMap->attributeName[i], keyWord ) == 0) {

            *index = attrMap->attributeIndex[i];
            return CAPS_SUCCESS;
        }
    }


    return CAPS_NOTFOUND;
}

// Search a mapAttrToIndex structure for a given index and return the corresponding keyword
int get_mapAttrToIndexKeyword(mapAttrToIndexStruct *attrMap, int index, const char **keyWord) {

    // If the keyword is not found a CAPS_NOTFOUND is returned

    int i; // Indexing

    *keyWord = NULL; // Default value of return

    if (attrMap == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < attrMap->numAttribute; i++) {

        if (index == attrMap->attributeIndex[i]) {

            *keyWord = attrMap->attributeName[i];
            return CAPS_SUCCESS;
        }
    }

    return CAPS_NOTFOUND;
}

// Set the index of a given keyword in a mapAttrToIndex structure
int set_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap, const char *keyWord, int index) {

    int i = 0; // Indexing

    if (attrMap == NULL) return CAPS_NULLVALUE;
    if (keyWord == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < attrMap->numAttribute; i++) {

        if (strcmp(attrMap->attributeName[i], keyWord ) == 0) {

            attrMap->attributeIndex[i] = index;
            return CAPS_SUCCESS;
        }
    }

    return CAPS_NOTFOUND;
}

// Increment a mapAttrToIndex structure with the given keyword and set the default index (= numAttribute)
int increment_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap, const char *keyWord) {

    int dummyIndex; // Dummy index
    int status = 0; // Returning status code
    int i; // Indexing

    if (attrMap == NULL) return CAPS_NULLVALUE;
    if (keyWord == NULL) return CAPS_NULLVALUE;

    status = get_mapAttrToIndexIndex(attrMap, keyWord, &dummyIndex);
    if (status == CAPS_NOTFOUND){
        attrMap->numAttribute += 1;
    } else {
        //printf("Attribute name already exists in attrMap structure !! \n");
        return EGADS_EXISTS;
    }

    if (attrMap->numAttribute == 1 ) {

        /*
        printf("Size of char** - %d\n", (int) sizeof(char **));
        printf("Size of const char* - %d\n", (int) sizeof(const char *));
        printf("Size of char* - %d\n", (int) sizeof(char *));
        printf("Size of char - %d\n", (int) sizeof(char));
        printf("Size of int* - %d\n", (int) sizeof(int *));
        printf("Size of int - %d\n", (int) sizeof(int));
         */

        attrMap->attributeName  = (char **) EG_alloc(attrMap->numAttribute*sizeof(char *));
        attrMap->attributeIndex = (int   *) EG_alloc(attrMap->numAttribute*sizeof(int   ));
    } else {

        attrMap->attributeName  = (char **) EG_reall(attrMap->attributeName , attrMap->numAttribute*sizeof(char *));
        attrMap->attributeIndex = (int   *) EG_reall(attrMap->attributeIndex, attrMap->numAttribute*sizeof(int   ));
    }

    if (attrMap->attributeName  == NULL) return EGADS_MALLOC;
    if (attrMap->attributeIndex == NULL) return EGADS_MALLOC;

    attrMap->attributeName[attrMap->numAttribute-1] = (char *) EG_alloc((strlen(keyWord) + 1)*sizeof(char));
    if (attrMap->attributeName[attrMap->numAttribute-1] ==  NULL ) {
        for (i = 0; i < attrMap->numAttribute; i++) {
            EG_free(attrMap->attributeName[i]);
        }
        attrMap->attributeName = NULL;
        return EGADS_MALLOC;
    }

    //printf("KEY WORD = %s\n",keyWord);
    sprintf(attrMap->attributeName[attrMap->numAttribute-1], "%s", keyWord);

    status = set_mapAttrToIndexStruct(attrMap, keyWord, attrMap->numAttribute);
    if (status != CAPS_SUCCESS) return status;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) an attribute map in the mapAttrToIndexStruct structure format
int initiate_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap) {

    if (attrMap == NULL) return CAPS_NULLVALUE;

    attrMap->mapName = NULL;

    attrMap->numAttribute = 0;
    attrMap->attributeName = NULL;
    attrMap->attributeIndex = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) an attribute map in the mapAttrToIndexStruct structure format
int destroy_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap) {

    int i; // Indexing

    if (attrMap == NULL) return CAPS_NULLVALUE;

    if (attrMap->mapName != NULL) EG_free(attrMap->mapName);
    attrMap->mapName = NULL;

    if (attrMap->attributeName != NULL) {

        for (i = 0; i < attrMap->numAttribute; i++) {

            if (attrMap->attributeName[i] != NULL) EG_free(attrMap->attributeName[i]);
        }

        EG_free(attrMap->attributeName);
    }

    if (attrMap->attributeIndex != NULL) EG_free(attrMap->attributeIndex);

    attrMap->attributeName  = NULL;
    attrMap->attributeIndex = NULL;

    attrMap->numAttribute = 0;

    return CAPS_SUCCESS;
}

// Make a copy of attribute map (attrMapIn)
int copy_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMapIn, mapAttrToIndexStruct *attrMapOut) {

    int status; // Function return status
    int i, j; // Indexing
    char *keyWord = NULL;

    if (attrMapIn  == NULL) return CAPS_NULLVALUE;
    if (attrMapOut == NULL) return CAPS_NULLVALUE;

    // Destroy attrMapOut in case it is already allocated - this implies that it must have at least been initiated
    status =  destroy_mapAttrToIndexStruct(attrMapOut);
    if (status != CAPS_SUCCESS) return status;

    attrMapOut->mapName = EG_strdup(attrMapIn->mapName);

    attrMapOut->numAttribute = attrMapIn->numAttribute;

    attrMapOut->attributeIndex = (int *) EG_alloc(attrMapOut->numAttribute*sizeof(int));
    if (attrMapOut->attributeIndex == NULL) return EGADS_MALLOC;

    memcpy(attrMapOut->attributeIndex,
            attrMapIn->attributeIndex,
            attrMapOut->numAttribute*sizeof(int));

    attrMapOut->attributeName  = (char **) EG_alloc(attrMapOut->numAttribute*sizeof(char *));
    if (attrMapOut->attributeName == NULL) return EGADS_MALLOC;

    for (i = 0; i < attrMapOut->numAttribute; i++){

        keyWord = attrMapIn->attributeName[i];

        attrMapOut->attributeName[i] = (char *) EG_alloc((strlen(keyWord) + 1)*sizeof(char));

        if (attrMapOut->attributeName[i] ==  NULL ) {

            for (j = 0; j < i; j++) {
                EG_free(attrMapOut->attributeName[j]);
            }

            EG_free(attrMapOut->attributeName);
            EG_free(attrMapOut->attributeIndex);

            attrMapOut->attributeName = NULL;
            attrMapOut->attributeIndex = NULL;
            return EGADS_MALLOC;
        }

        sprintf(attrMapOut->attributeName[i], "%s", keyWord);
    }


    return CAPS_SUCCESS;
}

// Merge two attribute maps preserving the order (and name) of the first input map.
int merge_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap1, mapAttrToIndexStruct *attrMap2, mapAttrToIndexStruct *attrMapOut) {

    int status; // Function return status
    int i; // Indexing

    if (attrMap1  == NULL) return CAPS_NULLVALUE;
    if (attrMap2  == NULL) return CAPS_NULLVALUE;
    if (attrMapOut == NULL) return CAPS_NULLVALUE;

    // Destroy attrMapOut in case it is already allocated - this implies that it must have at least been initiated
    status =  destroy_mapAttrToIndexStruct(attrMapOut);
    if (status != CAPS_SUCCESS) return status;

    status = copy_mapAttrToIndexStruct(attrMap1, attrMapOut);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < attrMap2->numAttribute; i++) {
        status = increment_mapAttrToIndexStruct(attrMapOut, attrMap2->attributeName[i]);
        if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
    }

    status = CAPS_SUCCESS;
    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in merge_mapAttrToIndexStruct, status = %d\n", status);

        return status;
}

// Retrieve the string following a generic tag (given by attributeKey)
int retrieve_stringAttr(ego geomEntity, const char *attributeKey, const char **string) {

    int status; // Status return

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;


    status = EG_attributeRet(geomEntity, attributeKey, &atype, &alen, &ints, &reals, string);
    if (status != EGADS_SUCCESS) return status;

    if (atype != ATTRSTRING) {
        printf ("Error: Attribute %s should be followed by a single string\n", attributeKey);
        return EGADS_ATTRERR;
    }

    return CAPS_SUCCESS;
}

// Retrieve an int following a generic tag (given by attributeKey)
int retrieve_intAttrOptional(ego geomEntity, const char *attributeKey, int *val) {

    int status; // Status return

    int atype = ATTRPTR, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char   *string;

    if (val == NULL) {
        return CAPS_NULLVALUE;
    }

    status = EG_attributeRet(geomEntity, attributeKey, &atype, &alen, &ints, &reals, &string);
    if (status != EGADS_SUCCESS) return status;

    if (!((atype == ATTRINT || atype == ATTRREAL) && alen == 1)) {
        printf ("Error: Attribute %s should be a single integer or real\n", attributeKey);
        return EGADS_ATTRERR;
    }

    if (atype == ATTRINT ) *val = ints[0];
    if (atype == ATTRREAL) *val = NINT(reals[0]);

    return CAPS_SUCCESS;
}


// Retrieve a double following a generic tag (given by attributeKey)
int retrieve_doubleAttrOptional(ego geomEntity, const char *attributeKey, double *val) {

    int status; // Status return

    int atype = ATTRPTR, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char   *string;

    if (val == NULL) {
        return CAPS_NULLVALUE;
    }

    status = EG_attributeRet(geomEntity, attributeKey, &atype, &alen, &ints, &reals, &string);
    if (status != EGADS_SUCCESS) return status;

    if (atype != ATTRREAL || alen != 1) {
        printf ("Error: Attribute %s should be a single real\n", attributeKey);
        return EGADS_ATTRERR;
    }

    if (atype == ATTRREAL) {
      *val = reals[0];
      return CAPS_SUCCESS;
    }

    return CAPS_NOTFOUND;
}

// Retrieve the string following a capsGroup tag
int retrieve_CAPSGroupAttr(ego geomEntity, const char **string) {

    int status; // Status return

    char *attributeKey = "capsGroup";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsConstraint tag
int retrieve_CAPSConstraintAttr(ego geomEntity, const char **string) {

    int status; // Status return

    char *attributeKey = "capsConstraint";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsLoad tag
int retrieve_CAPSLoadAttr(ego geomEntity, const char **string) {

    int status; // Status return

    char *attributeKey = "capsLoad";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsBound tag
int retrieve_CAPSBoundAttr(ego geomEntity, const char **string) {

    int status; // Status return

    char *attributeKey = "capsBound";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsIgnore tag
int retrieve_CAPSIgnoreAttr(ego geomEntity, const char **string) {

    int status;
    char *attributeKey = "capsIgnore";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsConnect tag
int retrieve_CAPSConnectAttr(ego geomEntity, const char **string) {

    int status;
    char *attributeKey = "capsConnect";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}


// Retrieve the string following a capsConnectLink tag
int retrieve_CAPSConnectLinkAttr(ego geomEntity, const char **string) {

    int status;
    char *attributeKey = "capsConnectLink";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the string following a capsResponse tag
int retrieve_CAPSResponseAttr(ego geomEntity, const char **string) {

    int status;
    char *attributeKey = "capsResponse";

    status = retrieve_stringAttr(geomEntity, attributeKey, string);
    return status;
}

// Retrieve the value following a capsDiscipline
int retrieve_CAPSDisciplineAttr(ego geomEntity, const char **string) {

    char *attributeKey = "capsDiscipline";

    return retrieve_stringAttr(geomEntity, attributeKey, string);
}

// Retrieve the value following a capsMesh
int retrieve_CAPSMeshAttr(ego geomEntity, const char **string) {

    char *attributeKey = "capsMesh";

    return retrieve_stringAttr(geomEntity, attributeKey, string);
}

/*
// Retrieve the string following a CoordSystem tag
int retrieve_CoordSystemAttr(ego geomEntity, const char **name) {

    int status; // Status return
    int i; // Indexing

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;

    int numAttr;

    status = EG_attributeNum(geomEntity, &numAttr);
    if (status != EGADS_SUCCESS) return EGADS_NOTFOUND;

    for (i = 0; i < numAttr; i++) {
        status = EG_attributeGet(geomEntity, i+1, name, &atype, &alen, &ints, &reals, &string);

        if (atype == ATTRCSYS) return CAPS_SUCCESS;
    }

    return EGADS_NOTFOUND;
}
*/

// Create a mapping between unique, generic (specified via mapName) attribute names and an index value
int create_genericAttrToIndexMap(int numBody, ego bodies[], int attrLevelIn, const char *mapName, mapAttrToIndexStruct *attrMap) {

    // In:
    //    numBody   = Number of incoming bodies
    //    bodies    = Array of ego bodies
    //    attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                              1 - search the body and all the faces
    //                              2 - search the body, faces, and all the edges
    //                             >2 - search the body, faces, edges, and all the nodes
    // Out:
    //         attrMap = A filled mapAttrToIndex structure

    int attr, body, face, edge, node; // Indexing variables

    int status; // Function return integer

    const char *groupName; // Attribute name

    int numFace = 0, numEdge = 0, numNode = 0; // Number of egos
    ego *faces = NULL, *edges = NULL, *nodes = NULL; // Geometry

    int attrLevel;

    // EGADS function returns
    ego  eref, *echilds;
    int  oclass, nchild, *senses, bodySubType = 0; // Body classification
    double data[4];

    // Destroy attrMap in case it is already allocated - this implies that is must have at least been initiated
    status =  destroy_mapAttrToIndexStruct(attrMap);
    if (status != CAPS_SUCCESS) return status;

    if (mapName == NULL) return CAPS_NULLVALUE;

    printf("Mapping %s attributes ................\n", mapName);

    attrMap->mapName = EG_strdup(mapName);
    if (attrMap->mapName == NULL) return EGADS_MALLOC;

    // Search through bodies
    for (body = 0; body < numBody; body++) {

        attrLevel = attrLevelIn;

        status = EG_getTopology(bodies[body], &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;
        if (oclass == NODE) attrLevel = 0; // If we have a node body - change attrLevel to just the body

        // Get groupName following mapName
        status = retrieve_stringAttr(bodies[body], mapName, &groupName);
        if (status != EGADS_SUCCESS && status != EGADS_NOTFOUND) goto cleanup;

        // Set attribute map
        if (status == CAPS_SUCCESS) {
            status = increment_mapAttrToIndexStruct(attrMap, groupName);
            if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
        }

        // Search through faces
        if (attrLevel > 0) {

            // Determine the number of faces
            status = EG_getBodyTopos(bodies[body], NULL, FACE, &numFace, &faces);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through faces
            for (face = 0; face < numFace; face++) {

                // Get groupName following mapName
                status = retrieve_stringAttr(faces[face], mapName, &groupName);
                if (status != EGADS_SUCCESS && status != EGADS_NOTFOUND) goto cleanup;

                // Set attribute map
                if (status == CAPS_SUCCESS) {
                    status = increment_mapAttrToIndexStruct(attrMap, groupName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End face loop

        // Search through edges
        if (attrLevel > 1) {
            status = EG_getBodyTopos(bodies[body], NULL, EDGE, &numEdge, &edges);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through edges
            for (edge = 0; edge < numEdge; edge++) {

                // Get groupName following mapName
                status = retrieve_stringAttr(edges[edge], mapName, &groupName);
                if (status != EGADS_SUCCESS && status != EGADS_NOTFOUND) goto cleanup;

                // Set attribute map
                if (status == CAPS_SUCCESS) {
                    status = increment_mapAttrToIndexStruct(attrMap, groupName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End edge loop

        // Search through nodes
        if (attrLevel > 2) {
            status = EG_getBodyTopos(bodies[body], NULL, NODE, &numNode, &nodes);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through nodes
            for (node = 0; node < numNode; node++) {

                // Get groupName following mapName
                status = retrieve_stringAttr(nodes[node], mapName, &groupName);
                if (status != EGADS_SUCCESS && status != EGADS_NOTFOUND) goto cleanup;

                // Set attribute map
                if (status == CAPS_SUCCESS) {
                    status = increment_mapAttrToIndexStruct(attrMap, groupName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End node loop

        AIM_FREE(faces);
        AIM_FREE(edges);
        AIM_FREE(nodes);

    } // End body loop

    printf("\tNumber of unique %s attributes = %d\n", mapName, attrMap->numAttribute);
    for (attr = 0; attr < attrMap->numAttribute; attr++) {
        printf("\tName = %s, index = %d\n", attrMap->attributeName[attr], attrMap->attributeIndex[attr]);
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(faces);
    AIM_FREE(edges);
    AIM_FREE(nodes);

    return status;
}

// Create a mapping between unique capsGroup attribute names and an index value
int create_CAPSGroupAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //         attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    const char *mapName = "capsGroup";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique CoordSystem attribute names and an index value
int create_CoordSystemAttrToIndexMap(int numBody, ego bodies[], int attrLevelIn, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int attr, body, face, edge, node; // Indexing variables

    int status; // Function return integer

    const char *csysName; // Attribute name

    int numFace = 0, numEdge = 0, numNode = 0; // Number of egos
    ego *faces = NULL, *edges = NULL, *nodes = NULL; // Geometry

    int attrLevel;

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;
    int numAttr;

    ego  eref, *echilds;
    int  oclass, nchild, *senses, bodySubType = 0; // Body classification
    double data[4];


    // Destroy attrMap in case it is already allocated - this implies that is must have at least been initiated
    status =  destroy_mapAttrToIndexStruct(attrMap);
    if (status != CAPS_SUCCESS) return status;

    printf("Mapping Csys attributes ................\n");

    // Search through bodies
    for (body = 0; body < numBody; body++) {

        attrLevel = attrLevelIn;

        status = EG_getTopology(bodies[body], &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;
        if (bodySubType == NOMTYPE) attrLevel = 0; // If we have a node body - change attrLevel to just the body

        // Get csysName following CoordSystem
        status = EG_attributeNum(bodies[body], &numAttr);
        if (status != EGADS_SUCCESS) return EGADS_NOTFOUND;

        for (attr = 0; attr < numAttr; attr++) {

            status = EG_attributeGet(bodies[body], attr+1, &csysName, &atype, &alen, &ints, &reals, &string);
            if (status != EGADS_SUCCESS) goto cleanup;

            if (atype != ATTRCSYS) continue;

            // Set attribute map
            status = increment_mapAttrToIndexStruct(attrMap, csysName);
            if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
        }

        // Search through faces
        if (attrLevel > 0) {

            // Determine the number of faces
            status = EG_getBodyTopos(bodies[body], NULL, FACE, &numFace, &faces);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through faces
            for (face = 0; face < numFace; face++) {

                status = EG_attributeNum(faces[face], &numAttr);
                if (status != EGADS_SUCCESS) return EGADS_NOTFOUND;

                for (attr = 0; attr < numAttr; attr++) {

                    status = EG_attributeGet(faces[face], attr+1, &csysName, &atype, &alen, &ints, &reals, &string);
                    if (status != EGADS_SUCCESS) goto cleanup;

                    if (atype != ATTRCSYS) continue;

                    // Set attribute map
                    status = increment_mapAttrToIndexStruct(attrMap, csysName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End face loop

        // Search through edges
        if (attrLevel > 1) {
            status = EG_getBodyTopos(bodies[body], NULL, EDGE, &numEdge, &edges);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through edges
            for (edge = 0; edge < numEdge; edge++) {

                status = EG_attributeNum(edges[edge], &numAttr);
                if (status != EGADS_SUCCESS) return EGADS_NOTFOUND;

                for (attr = 0; attr < numAttr; attr++) {

                    status = EG_attributeGet(edges[edge], attr+1, &csysName, &atype, &alen, &ints, &reals, &string);
                    if (status != EGADS_SUCCESS) goto cleanup;

                    if (atype != ATTRCSYS) continue;

                    // Set attribute map
                    status = increment_mapAttrToIndexStruct(attrMap, csysName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End edge loop

        // Search through nodes
        if (attrLevel > 2) {
            status = EG_getBodyTopos(bodies[body], NULL, NODE, &numNode, &nodes);
            if (status != EGADS_SUCCESS) return status;

            // Loop through nodes
            for (node = 0; node < numNode; node++) {

                status = EG_attributeNum(nodes[node], &numAttr);
                if (status != EGADS_SUCCESS) return EGADS_NOTFOUND;

                for (attr = 0; attr < numAttr; attr++) {

                    status = EG_attributeGet(nodes[node], attr+1, &csysName, &atype, &alen, &ints, &reals, &string);
                    if (status != EGADS_SUCCESS) goto cleanup;

                    if (atype != ATTRCSYS) continue;

                    // Set attribute map
                    status = increment_mapAttrToIndexStruct(attrMap, csysName);
                    if (status != CAPS_SUCCESS && status != EGADS_EXISTS) goto cleanup;
                }
            }
        } // End node loop


        EG_free(faces); faces = NULL;
        EG_free(edges); edges = NULL;
        EG_free(nodes); nodes = NULL;

    } // End body loop

    printf("\tNumber of unique Csys attributes = %d\n", attrMap->numAttribute);
    for (attr = 0; attr < attrMap->numAttribute; attr++) {
        printf("\tName = %s, index = %d\n", attrMap->attributeName[attr], attrMap->attributeIndex[attr]);
    }

    status = CAPS_SUCCESS;

    cleanup:
        EG_free(faces);
        EG_free(edges);
        EG_free(nodes);

        return status;
}

// Create a mapping between unique capsConstraint attribute names and an index value
int create_CAPSConstraintAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsConstraint";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique capsLoad attribute names and an index value
int create_CAPSLoadAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsLoad";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique capsBound attribute names and an index value
int create_CAPSBoundAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsBound";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique capsConnect attribute names and an index value
int create_CAPSConnectAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsConnect";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique capsResponse attribute names and an index value
int create_CAPSResponseAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsResponse";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Create a mapping between unique capsMesh attribute names and an index value
int create_CAPSMeshAttrToIndexMap(int numBody, ego bodies[], int attrLevel, mapAttrToIndexStruct *attrMap) {

    // In:
    //      numBody   = Number of incoming bodies
    //      bodies    = Array of ego bodies
    //      attrLevel = Level of depth to traverse the body:  0 - search just body attributes
    //                                                        1 - search the body and all the faces
    //                                                        2 - search the body, faces, and all the edges
    //                                                       >2 - search the body, faces, edges, and all the nodes
    // Out:
    //      attrMap = A filled mapAttrToIndex structure

    int status; // Function return integer

    char *mapName = "capsMesh";

    status = create_genericAttrToIndexMap(numBody, bodies, attrLevel, mapName, attrMap);

    return status;
}

// Check capsLength consistency in the bodies and return the value. No check is done to make sure ALL bodies have
// a capsLength, just that if present it is consistent.
int check_CAPSLength(int numBody, ego bodies[], const char **lengthString) {

    int i;
    int status;

    int found = (int) false;

    // EGADS return values
    const char   *string;

    *lengthString = NULL;

    for (i =0; i < numBody; i++) {

        status = retrieve_stringAttr(bodies[i], "capsLength", &string);
        if (status != CAPS_SUCCESS)    continue;

        // Save the string
        if (found == (int) false) {
            found = true;
            *lengthString = string;

        } else { // Compare the strings

            if (strcasecmp(*lengthString, string) != 0) {
                printf("Inconsistent length units on bodies, capsLength %s found on one body, "
                        "while %s found on another\n", *lengthString, string);
                return CAPS_MISMATCH;
            }
        }
    }

    if (found == (int) false) return CAPS_NOTFOUND;

    return CAPS_SUCCESS;
}


// Check capsDiscipline consistency in the bodies and return the value.
// All bodies must either have or not have a capsDiscipline.
int check_CAPSDiscipline(int numBody, ego bodies[], const char **discpline) {

    int i;
    int status;

    const char *currentdiscpline = 0;

    *discpline = NULL;

    // Check bodies
    for (i = 0; i < numBody; i++) {

        status = retrieve_CAPSDisciplineAttr(bodies[i], &currentdiscpline);
        if (status != EGADS_SUCCESS) continue;

        if (*discpline == NULL) {
            *discpline = currentdiscpline;
        } else {

            if (strcasecmp(*discpline, currentdiscpline) != 0) {
                printf("All bodies don't have the same capsDiscipline value - one body "
                        "found with = %s, while another has %s!\n", *discpline, currentdiscpline);

                *discpline = NULL;
                return CAPS_MISMATCH;
            }
        }
    }

    // If we found at least one discipline, re-check bodies because the first bodies checked might not have had an discipline
    if (*discpline != NULL) {
        for (i = 0; i < numBody; i++) {

            status = retrieve_CAPSDisciplineAttr(bodies[i], &currentdiscpline);

            if (status != CAPS_SUCCESS) {
                printf("A capsDiscipline value of %s was found, all bodies must have this value as well!\n", *discpline);
                *discpline = NULL;
                return status;
            }
        }
    }

    return CAPS_SUCCESS;
}


// Check capsMeshLength consistency in the bodies and return the value. No check is done to make sure ALL bodies have
// a capsMeshLength, just that if present it is consistent.
int check_CAPSMeshLength(int numBody, ego bodies[], double *capsMeshLength) {

    int i;
    int status;

    int found = (int) false;

    // EGADS return values
    double val = 0;

    for (i =0; i < numBody; i++) {

        status = retrieve_doubleAttrOptional(bodies[i], "capsMeshLength", &val);
        if (status == EGADS_ATTRERR) return status;
        if (status != CAPS_SUCCESS)  continue;

        // Save the string
        if (found == (int) false) {
            found = true;
            *capsMeshLength = val;

        } else { // Compare the strings

            if (*capsMeshLength != val) {
                printf("Inconsistent mesh length on bodies, capsMeshLength %lf found on one body, "
                        "while %lf found on another\n", *capsMeshLength, val);
                return CAPS_MISMATCH;
            }
        }
    }

    if (found == (int) false) return CAPS_NOTFOUND;

    return CAPS_SUCCESS;
}

// Copy an integer array, out[] is freeable
int copy_intArray(int length, int *in, int *out[]) {

    int i;
    int *temp;

    *out = NULL;

    if (in != NULL && length > 0) {
        temp = (int *) EG_alloc(length*sizeof(int));
        if (temp == NULL) return EGADS_MALLOC;
        for (i = 0; i < length; i++) temp[i] = in[i];

    } else {
         temp = in;
     }

    *out = temp;

    return CAPS_SUCCESS;
}
