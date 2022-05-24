// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Miscellaneous related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include "stdarg.h"
#include "meshTypes.h"  // Bring in mesh structures
#include "capsTypes.h"  // Bring in CAPS types
#include "cfdTypes.h"   // Bring in cfd structures
#include "miscTypes.h"  // Bring in miscellaneous structures

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MSC_VER // Microsoft compilers
#if (_MSC_VER < 1900)
#define __func__  __FUNCTION__
#endif
#endif

#ifndef S_SPLINT_S
#define PRINT_ERROR(reason, ...) \
    fprintf(stderr, \
        "\n" \
        "--------------------------------------------------------------------------------\n" \
        "[ ERROR ]\n" \
        "File: %s\n" \
        "Line: %d\n" \
        "Function: %s\n" \
        "Reason: " reason "\n" \
        "--------------------------------------------------------------------------------\n" \
        "\n", \
        __FILE__, __LINE__, __func__, ##__VA_ARGS__);

#define PRINT_WARNING(fmt, ...) \
    fprintf(stderr, "[ WARNING ] " fmt "\n", ##__VA_ARGS__);
#else
extern void PRINT_ERROR(const char *reason, ...);
extern void PRINT_WARNING(const char *fmt, ...);
#endif

// Convert a string in tuple form to an array of strings - tuple is assumed to be bounded by '(' and ')' and comma separated
//  for example ["3.0", "hey", "foo"] - note the strings should NOT contain commas. If the
//  string coming in is not a tuple the string is simply copied. Also quotations are removed
int json_parseTuple(/*@null@*/ const char *stringToParse, int *arraySize,
                    char **stringArray[]);

// Simple json dictionary parser - currently doesn't support nested arrays for keyValue
int search_jsonDictionary(const char *stringToSearch, const char *keyWord,
                          char **keyValue);

// Convert a single string to a prog_argv array
int string_toProgArgs(char *meshInputString, int *prog_argc, char ***prog_argv);

// Convert a string to double
int string_toDouble(const char *string, double *number);

// Convert a string tuple of double with data-units to double in units
int string_toDoubleUnits(void *aimInfo, const char *string, const char *units, double *number);

// Convert a string to double array
int string_toDoubleArray(char *string, int arraySize, double numberArray[]);

// Convert a string tuple of array double with data-units to array double in units
int string_toDoubleArrayUnits(void *aimInfo, const char *string,
                              const char *units, int arraySize, double *numberArray);

// Convert a string to an array of doubles
int string_toDoubleDynamicArray(char *stringToSearch, int *arraySize,
                                double *numberArray[]);

// Convert a string tuple of array double with data-units to array double in units
int string_toDoubleDynamicArrayUnits(void *aimInfo, const char *string,
                                     const char *units, int *arraySize, double **numberArray);

// Convert a string to an array of strings
int string_toStringArray(char *stringToSearch, int arraySize,
                         char *stringArray[]);

// Convert a string to an array of strings
int string_toStringDynamicArray(char *stringToSearch, int *arraySize,
                                char **stringArray[]);

// Convert a string to boolean
int string_toBoolean(char *string, int *number);

// Convert a string to integer
int string_toInteger(const char *string, int *number);

// Convert a string to integer array
int string_toIntegerArray(char *stringToSearch, int arraySize,
                          int numberArray[]);

// Convert a string to an array of integer
int string_toIntegerDynamicArray(char *stringToSearch, int *arraySize,
                                 int *numberArray[]);

// Free an array strings
int string_freeArray(int numString, char **strings[]);

// Remove quotation marks around a string - Returning char * should be free'd after use
/*@null@*/
char * string_removeQuotation(/*@null@*/ const char *string);

// Force a string to upper case
void string_toUpperCase ( char *sPtr );

// Create formatted string from format string and variadic string args
// NUmber of variadic args determined by NULL sentinel or number of "%s" specifiers in format string
// NOTE: Returning char * should be free'd after use
/*@null@*/
char * string_format(char *format, ...);

// Return whether string `find` is in `array`
int string_isInArray(char *find, int arraySize, char **array);

// The max x,y,z coordinates where P(3*i + 0) = x_i, P(3*i + 1) = y_i, and P(3*i + 2) = z_i
void maxCoords(int sizeP, double *P, double *x, double *y, double *z);

// The min x,y,z coordinates where P(3*i + 0) = x_i, P(3*i + 1) = y_i, and P(3*i + 2) = z_i
void minCoords(int sizeP, double *P, double *x, double *y, double *z);

// Return the endianness of the machine we are working on - 0 = little, 1 = big
int get_MachineENDIANNESS(void);

#ifdef WIN32
int getline(char **linep, size_t *linecapp, FILE *stream);
#endif

// Return the max value of two doubles
double max_DoubleVal(double x1, double x2);

// Return the max value of two doubles
double min_DoubleVal(double x1, double x2);

// Cross product axb = c
void cross_DoubleVal(double a[], double b[], double c[]);

// Dot product a*b = c
double dot_DoubleVal(double a[], double b[]);

// Distance between two points sqrt(dot(a-b,a-b))
double dist_DoubleVal(double a[], double b[]);

// Convert an integer to a string of a given field width and justifacation - Returning char * should be free'd after use
/*@null@*/
char * convert_integerToString(int integerVal, int fieldWidth, int leftOrRight);

// Convert an double to a string (scientific notation is used depending on the fieldwidth and value) of a given field width
//      - Returning char * should be free'd after use
/*@null@*/
char * convert_doubleToString(double doubleVal, int fieldWidth, int leftOrRight);

// Factorizes in place the square linear system A using simple LU decomposition
int factorLU(int n, double A[] );

// Solves the factorized square linear system LU x = b
int backsolveLU(int n, double LU[], double b[], double x[] );

// Solves the square linear system A x = b using simple LU decomposition
int solveLU(int n, double A[], double b[], double x[] );

// Prints all attributes on an ego
int print_AllAttr( void *aimInfo, ego obj );

// Initiate (0 out all values and NULL all pointers) an attribute map in the mapAttrToIndexStruct structure format
int initiate_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap);

// Destroy (0 out all values and NULL all pointers) an attribute map in the mapAttrToIndexStruct structure format
int destroy_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap);

// Make a copy of attribute map (attrMapIn)
int copy_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMapIn,
                              mapAttrToIndexStruct *attrMapOut);

// Merge two attribute maps preserving the order (and name) of the first input map.
int merge_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap1,
                               mapAttrToIndexStruct *attrMap2,
                               mapAttrToIndexStruct *attrMapOut);

// Search a mapAttrToIndex structure for a given keyword and set/return the corresponding index
int get_mapAttrToIndexIndex(mapAttrToIndexStruct *attrMap, const char *keyWord,
                            int *index);

// Search a mapAttrToIndex structure for a given index and return the corresponding keyword
int get_mapAttrToIndexKeyword(mapAttrToIndexStruct *attrMap, int index,
                              const char **keyWord);

// Set the index of a given keyword in a mapAttrToIndex structure
int set_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap, const char *keyWord,
                             int index);

// Increment a mapAttrToIndex structure with the given keyword and set the default index (= numAttribute)
int increment_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap,
                                   const char *keyWord);

// Retrieve the string following a generic tag (given by attributeKey)
int retrieve_stringAttr(ego geomEntity, const char *attributeKey,
                        const char **string);

// Retrieve an int following a generic tag (given by attributeKey)
int retrieve_intAttrOptional(ego geomEntity, const char *attributeKey,
                             int *val);

// Retrieve a double following a generic tag (given by attributeKey)
int retrieve_doubleAttrOptional(ego geomEntity, const char *attributeKey,
                                double *val);

// Retrieve the string following a capsGroup tag
int retrieve_CAPSGroupAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsConstraint tag
int retrieve_CAPSConstraintAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsLoad tag
int retrieve_CAPSLoadAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsBound tag
int retrieve_CAPSBoundAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsIgnore tag
int retrieve_CAPSIgnoreAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsConnect tag
int retrieve_CAPSConnectAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsConnectLink tag
int retrieve_CAPSConnectLinkAttr(ego geomEntity, const char **string);

// Retrieve the string following a capsResponse tag
int retrieve_CAPSResponseAttr(ego geomEntity, const char **string);

// Retrieve the value following a capsDiscipline
int retrieve_CAPSDisciplineAttr(ego geomEntity, const char **string);

// Retrieve the value following a capsMesh
int retrieve_CAPSMeshAttr(ego geomEntity, const char **string);

/*
// Retrieve the string following a CoordSystem tag
int retrieve_CoordSystemAttr(ego geomEntity, const char **string);
*/

// Create a mapping between unique, generic (specified via mapName) attribute names and an index value
int create_genericAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                 const char *mapName,
                                 mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsGroup attribute names and an index value
int create_CAPSGroupAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                   mapAttrToIndexStruct *attrMap);

// Create a mapping between unique CoordSystem attribute names and an index value
int create_CoordSystemAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                     mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsConstraint attribute names and an index value
int create_CAPSConstraintAttrToIndexMap(int numBody, ego bodies[],
                                        int attrLevel,
                                        mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsLoad attribute names and an index value
int create_CAPSLoadAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                  mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsBound attribute names and an index value
int create_CAPSBoundAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                   mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsConnect attribute names and an index value
int create_CAPSConnectAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                     mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsResponse attribute names and an index value
int create_CAPSResponseAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                      mapAttrToIndexStruct *attrMap);

// Create a mapping between unique capsMesh attribute names and an index value
int create_CAPSMeshAttrToIndexMap(int numBody, ego bodies[], int attrLevel,
                                  mapAttrToIndexStruct *attrMap);

// Check capsLength consistency in the bodies and return the value. No check is done to make sure ALL bodies have
// a capsLength, just that if present it is consistent.
int check_CAPSLength(int numBody, ego bodies[], const char **lengthString);

// Check capsDiscipline consistency in the bodies and return the value.
// All bodies must either have or not have a capsDiscipline.
int check_CAPSDiscipline(int numBody, ego bodies[], const char  **discpline);

// Check capsMeshLength consistency in the bodies and return the value. No check is done to make sure ALL bodies have
// a capsMeshLength, just that if present it is consistent.
int check_CAPSMeshLength(int numBody, ego bodies[], double *capsMeshLength);

// Copy an integer array, out[] is freeable
int copy_intArray(int length, int *in, int *out[]);

#ifdef __cplusplus
}
#endif
