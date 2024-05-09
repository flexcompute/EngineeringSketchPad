#ifndef WIN32
#define _GNU_SOURCE // this ensures that POSIX `getline` is defined in stdio.h if supported
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
// #include <assert.h>
#include <ctype.h>

#include "aimUtil.h"
#include "miscUtils.h"
#include "zaeroUtils.h"

#ifdef WIN32
#define strcasecmp  stricmp
#define strtok_r   strtok_s
#endif

/* general utility */

int zaero_setBlankIntegerArray(int size, int **array) {
    int i;
    if (array == NULL) return CAPS_NULLVALUE;
    *array = (int *) EG_alloc(sizeof(int) * size);
    if (*array == NULL) return EGADS_MALLOC;
    for (i = 0; i < size; i++) (*array)[i] = -1;
    return CAPS_SUCCESS;
}

int zaero_setBlankRealArray(int size, double **array) {
    int i;
    if (array == NULL) return CAPS_NULLVALUE;
    *array = (double *) EG_alloc(sizeof(double) * size);
    if (*array == NULL) return EGADS_MALLOC;
    for (i = 0; i < size; i++) (*array)[i] = -1.0;
    return CAPS_SUCCESS;
}

int zaero_setBlankStringArray(int size, char ***array) {
    int i;
    if (array == NULL) return CAPS_NULLVALUE;
    *array = (char **) EG_alloc(sizeof(char *) * size);
    if (*array == NULL) return EGADS_MALLOC;
    for (i = 0; i < size; i++) (*array)[i] = NULL;
    return CAPS_SUCCESS;
}


/* file utility */

void zaero_readFile(FILE *fp, char **buffer) {

    long fsize;

    fseek(fp, 0, SEEK_END);
    fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);  /* same as rewind(f); */

    if (*buffer == NULL) {
        *buffer = malloc(fsize + 1);
    } else {
        *buffer = realloc(*buffer, (fsize + 1));
    }
    if (*buffer == NULL) return;

    fread(*buffer, 1, fsize, fp);

    (*buffer)[fsize] = 0;
}

// int zaero_fileExists(char *filename) {
//     FILE *fp;
//     if ((fp = fopen(filename, "r")) != NULL) {
//         // file exists
//         fclose(fp);
//         return (int) true;
//     }
//     else {
//         return (int) false;
//     }
// }

static inline int _areEqualIntegers(void *a, void *b) {
    return *((int *) a) == *((int *) b);
}

static inline int _areEqualDoubles(void *a, void *b) {
    return *((double *) a) == *((double *) b);
}

static inline int _areEqualStrings(void *a, void *b) {
    return strcmp(*((char **) a), *((char **) b)) == 0;
}

int zaero_allEqual(int type, void *a, void *b, ...) {

    typedef int (*_areEqualType) (void *, void *);

    _areEqualType _areEqual;
    va_list args;

    switch(type) {
        case Integer:
            _areEqual = _areEqualIntegers;
            break;
        case Double:
            _areEqual = _areEqualDoubles;
            break;
        case String:
            _areEqual = _areEqualStrings;
            break;
        default:
            printf("Error: Unrecognized `type` in zaero_allEqual. "
                   "Must be one of: Integer, Double, or String.\n");
            return (int) false;
    }

    va_start(args, b);

    while (b != NULL) {
        if (!_areEqual(a, b)) {
            return (int) false;
        }
        a = b;
        b = va_arg(args, void *);
    }
    va_end(args);

    return (int) true;
}

// finders

// find zaeroSubcaseStruct with name member equal `matchName`
int zaero_findSubcaseByName(void *aimInfo,
                            const char *matchName,
                            const zaeroProblemStruct *zaeroProblem,
                            zaeroSubcaseStruct **matchSubcase)
{
  int i;
  for (i = 0; i < zaeroProblem->numSubcases; i++) {
    if (strcmp(zaeroProblem->subcases[i].name, matchName) == 0) {
      *matchSubcase = &zaeroProblem->subcases[i];
      return CAPS_SUCCESS;
    }
  }

  AIM_ERROR(aimInfo, "No analysis subcase found with name: %s\n", matchName);
  return CAPS_NOTFOUND;
}

// find UAICStruct with name member equal `matchName`
int zaero_findUAICByName(void *aimInfo,
                         const char *matchName,
                         const zaeroProblemStruct *zaeroProblem,
                         zaeroUAICStruct **matchUAIC)
{
  int i;
  for (i = 0; i < zaeroProblem->numUAICs; i++) {
    if (strcmp(zaeroProblem->UAICs[i].name, matchName) == 0) {
      *matchUAIC = &zaeroProblem->UAICs[i];
      return CAPS_SUCCESS;
    }
  }

  AIM_ERROR(aimInfo, "No UAIC configuration found with name: %s\n", matchName);
  return CAPS_NOTFOUND;
}

