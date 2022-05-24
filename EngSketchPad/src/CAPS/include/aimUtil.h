#ifndef AIMUTIL_H
#define AIMUTIL_H
/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Utility Function Prototypes
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#define PATH_MAX _MAX_PATH
#else
#include <limits.h>
#endif

/* splint fixes */
#ifdef S_SPLINT_S
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
/*@-incondefs@*/
extern int fclose (/*@only@*/ FILE *__stream);
/*@+incondefs@*/
#endif

#include "capsTypes.h"


#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

__ProtoExt__ void
  aim_capsRev( int *major, int *minor );

__ProtoExt__ int
  aim_getRootPath( void *aimInfo, const char **fullPath );

__ProtoExt__ int
  aim_file( void *aimInfo, const char *file, char *aimFile );

__ProtoExt__ /*@null@*/ /*@out@*/ /*@only@*/ FILE *
  aim_fopen( void *aimInfo, const char *file, const char *mode );

__ProtoExt__ int
  aim_isFile(void *aimStruc, const char *file);

__ProtoExt__ int
  aim_rmFile(void *aimStruc, const char *file);

__ProtoExt__ int
  aim_cpFile(void *aimStruc, const char *src, const char *dst);

__ProtoExt__ int
  aim_relPath(void *aimStruc, const char *src,
              /*@null@*/ const char *dst, char *relPath);

__ProtoExt__ int
  aim_symLink(void *aimStruc, const char *src, /*@null@*/ const char *dst);

__ProtoExt__ int
  aim_isDir(void *aimStruc, const char *path);

__ProtoExt__ int
  aim_mkDir(void *aimStruc, const char *path);

__ProtoExt__ int
  aim_rmDir(void *aimStruc, const char *path);

__ProtoExt__ int
  aim_system( void *aimInfo, /*@null@*/ const char *rpath,
              const char *command );

__ProtoExt__ int
  aim_getBodies( void *aimInfo, const char **intent, int *nBody, ego **bodies );
  
__ProtoExt__ int
  aim_newGeometry( void *aimInfo );

__ProtoExt__ int
  aim_newAnalysisIn(void *aimStruc, int index);

__ProtoExt__ int
  aim_numInstance( void *aimInfo );

__ProtoExt__ int
  aim_getInstance( void *aimInfo );
  
__ProtoExt__ int
  aim_getUnitSys( void *aimInfo, char **unitSys );

__ProtoExt__ int
  aim_convert( void *aimInfo, const int count,
               /*@null@*/ const char  *inUnits, double  *inValue,
               /*@null@*/ const char *outUnits, double *outValue );

__ProtoExt__ int
  aim_unitMultiply( void *aimInfo, const char  *inUnits1, const char *inUnits2,
                    char **outUnits );

__ProtoExt__ int
  aim_unitDivide( void *aimInfo, const char  *inUnits1, const char *inUnits2,
                  char **outUnits );

__ProtoExt__ int
  aim_unitInvert( void *aimInfo, const char  *inUnits,
                  char **outUnits );

__ProtoExt__ int
  aim_unitRaise( void *aimInfo, const char  *inUnits, const int power,
                 char **outUnits );

__ProtoExt__ int
  aim_unitOffset( void *aimStruc, const char *inUnit, const double offset,
                  char **outUnits );

__ProtoExt__ int
  aim_getIndex( void *aimInfo, /*@null@*/ const char *name,
                enum capssType subtype );
  
__ProtoExt__ int
  aim_getValue( void *aimInfo, int index, enum capssType subtype,
                capsValue **value );

__ProtoExt__ int
  aim_initValue( capsValue *value );

__ProtoExt__ void
  aim_freeValue(capsValue *value);

__ProtoExt__ int
  aim_makeDynamicOutput( void *aimInfo, const char *dynObjName,
                         capsValue *value );
  
__ProtoExt__ int
  aim_getName( void *aimInfo, int index, enum capssType subtype,
               const char **name );

__ProtoExt__ int
  aim_getGeomInType( void *aimInfo, int index );
  
__ProtoExt__ int
  aim_newTess( void *aimInfo, ego tess );
  
__ProtoExt__ int
  aim_getDiscr( void *aimInfo, const char *bname, capsDiscr **discr );
  
__ProtoExt__ int
  aim_getDiscrState( void *aimInfo, const char *bname );
  
__ProtoExt__ int
  aim_getDataSet( capsDiscr *discr, const char *dname, enum capsdMethod *method,
                  int *npts, int *rank, double **data, char **units );

__ProtoExt__ int
  aim_getBounds( void *aimInfo, int *nTname, char ***tnames );

__ProtoExt__ int
  aim_valueAttrs( void *aimInfo, int index, enum capssType stype, int *nValue,
                  char ***names, capsValue **values );

__ProtoExt__ int
  aim_analysisAttrs( void *aimInfo, int *nValue, char ***names,
                     capsValue **values );

__ProtoExt__ void
  aim_freeAttrs( int nValue, char **names, capsValue *values );
  
__ProtoExt__ int
  aim_setSensitivity( void *aimInfo, const char *GIname, int irow, int icol );
  
__ProtoExt__ int
  aim_getSensitivity( void *aimInfo, ego body, int ttype, int index, int *npts,
                      double **dxyz );

__ProtoExt__ int
  aim_tessSensitivity( void *aimInfo, const char *GIname, int irow, int icol,
                       ego tess, int *npts, double **dxyz );

__ProtoExt__ int
  aim_isNodeBody( ego body, double *xyz );


  /* utility functions for nodal and cell centered data transfers */
__ProtoExt__ void
  aim_initBodyDiscr(capsBodyDiscr *discBody);

__ProtoExt__ int
  aim_nodalTriangleType(capsEleType *eletype);

__ProtoExt__ int
  aim_nodalQuadType(capsEleType *eletype);

__ProtoExt__ int
  aim_cellTriangleType(capsEleType *eletype);

__ProtoExt__ int
  aim_cellQuadType(capsEleType *eletype);

__ProtoExt__ int
  aim_locateElement( capsDiscr *discr, double *params,
                     double *param, int *bIndex, int *eIndex,
                     double *bary);
__ProtoExt__ int
  aim_interpolation(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                    double *bary, int rank, double *data, double *result);

__ProtoExt__ int
  aim_interpolateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                     double *bary, int rank, double *r_bar, double *d_bar);

__ProtoExt__ int
  aim_integration(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                  int rank, double *data, double *result);

__ProtoExt__ int
  aim_integrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                   int rank, double *r_bar, double *d_bar);

__ProtoExt__ void
  aim_status(void *aimInfo, const int status, const char *file, const int line,
             const char *func, const int narg, ...);

__ProtoExt__ void
  aim_message(void *aimInfo, enum capseType etype, int index, const char *file,
              const int line, const char *func, const char *format, ...);

__ProtoExt__ void
  aim_addLine(void *aimInfo, const char *format, ...);

__ProtoExt__ void
  aim_setIndexError(void *aimInfo, int index);

__ProtoExt__ void
  aim_removeError(void *aimInfo);


#ifndef S_SPLINT_S
/* Macro for counting varadaic arguments
 *
 * https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments/2124385#2124385
 */

#ifdef _MSC_VER // Microsoft compilers

#if (_MSC_VER < 1900)
#define __func__  __FUNCTION__
#endif

#  define GET_ARG_COUNT(...) INTERNAL_EXPAND_ARGS_PRIVATE(INTERNAL_ARGS_AUGMENTER(__VA_ARGS__))

#  define INTERNAL_ARGS_AUGMENTER(...) unused, __VA_ARGS__
#  define INTERNAL_EXPAND(x) x
#  define INTERNAL_EXPAND_ARGS_PRIVATE(...) \
    INTERNAL_EXPAND(INTERNAL_GET_ARG_COUNT_PRIVATE(__VA_ARGS__, \
                                                   69, 68, 67, 66, 65, 64, 63, 62, \
                                                   61, 60, 59, 58, 57, 56, 55, 54, \
                                                   53, 52, 51, 50, 49, 48, 47, 46, \
                                                   45, 44, 43, 42, 41, 40, 39, 38, \
                                                   37, 36, 35, 34, 33, 32, 31, 30, \
                                                   29, 28, 27, 26, 25, 24, 23, 22, \
                                                   21, 20, 19, 18, 17, 16, 15, 14, \
                                                   13, 12, 11, 10,  9,  8,  7,  6, \
                                                    5,  4,  3,  2,  1,  0))
#  define INTERNAL_GET_ARG_COUNT_PRIVATE(_1_, _2_, _3_, _4_, _5_, \
                                         _6_, _7_, _8_, _9_, _10_, \
                                         _11_, _12_, _13_, _14_, _15_, \
                                         _16_, _17_, _18_, _19_, _20_, \
                                         _21_, _22_, _23_, _24_, _25_, \
                                         _26_, _27_, _28_, _29_, _30_, \
                                         _31_, _32_, _33_, _34_, _35_, \
                                         _36, _37, _38, _39, _40, _41, \
                                         _42, _43, _44, _45, _46, _47, \
                                         _48, _49, _50, _51, _52, _53, \
                                         _54, _55, _56, _57, _58, _59, \
                                         _60, _61, _62, _63, _64, _65, \
                                         _66, _67, _68, _69, _70, count, ...) count

#else // Non-Microsoft compilers

#if defined(__clang__)
#  if __clang_major__ > 3 || (__clang_major__ == 3  && __clang_minor__ >= 4)
#    if __has_warning("-Wgnu-zero-variadic-macro-arguments")
#      pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#    endif
#  endif
#endif

#  define GET_ARG_COUNT(...) \
    INTERNAL_GET_ARG_COUNT_PRIVATE(0, ## __VA_ARGS__, 70, 69, 68, 67, 66, 65, 64, \
                                  63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, \
                                  51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, \
                                  39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, \
                                  27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, \
                                  15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0)
#  define INTERNAL_GET_ARG_COUNT_PRIVATE(_0, _1_, _2_, _3_, _4_, _5_, \
                                         _6_, _7_, _8_, _9_, _10_, \
                                         _11_, _12_, _13_, _14_, _15_, \
                                         _16_, _17_, _18_, _19_, _20_, \
                                         _21_, _22_, _23_, _24_, _25_, \
                                         _26_, _27_, _28_, _29_, _30_, \
                                         _31_, _32_, _33_, _34_, _35_, \
                                         _36, _37, _38, _39, _40, _41, \
                                         _42, _43, _44, _45, _46, _47, \
                                         _48, _49, _50, _51, _52, _53, \
                                         _54, _55, _56, _57, _58, _59, \
                                         _60, _61, _62, _63, _64, _65, \
                                         _66, _67, _68, _69, _70, count, ...) count

#endif

#define AIM_STATUS(aimInfo, status, ...) \
 if (status != CAPS_SUCCESS) { \
   aim_status(aimInfo, status, __FILE__, __LINE__, __func__, GET_ARG_COUNT(__VA_ARGS__), ##__VA_ARGS__); \
   goto cleanup; \
 }

#define AIM_ANALYSISIN_ERROR(aimInfo, index, format, ...) \
 { aim_message(aimInfo, CERROR, index, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); }

#define AIM_ERROR(aimInfo, format, ...) \
 { aim_message(aimInfo, CERROR, 0    , __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); }

#define AIM_WARNING(aimInfo, format, ...) \
 { aim_message(aimInfo, CWARN , 0    , __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); }

#define AIM_INFO(aimInfo, format, ...) \
 { aim_message(aimInfo, CINFO , 0    , __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); }

#define AIM_ADDLINE(aimInfo, format, ...) \
 { aim_addLine(aimInfo, format, ##__VA_ARGS__); }

#else

extern void AIM_STATUS(void *aimInfo, int status, ...);

extern void AIM_ANALYSISIN_ERROR(void *aimInfo, int index, const char *format, ...);

extern void AIM_ERROR(void *aimInfo, const char *format, ...);

extern void AIM_WARNING(void *aimInfo, const char *format, ...);

extern void AIM_INFO(void *aimInfo, const char *format, ...);

extern void AIM_ADDLINE(void *aimInfo, const char *format, ...);

extern ssize_t getline(char ** restrict linep, size_t * restrict linecapp,
                       FILE * restrict stream);

#endif

/* macro for crating names from enums */
#define AIM_NAME(a)       EG_strdup(#a)

#define AIM_ALLOC(ptr, size, type, aimInfo, status) \
 { \
   if (ptr != NULL) { \
      status = EGADS_MALLOC; \
      aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 1, "AIM_ALLOC: %s != NULL", #ptr); \
      goto cleanup; \
   } \
   size_t memorysize = size; \
   ptr = (type *) EG_alloc(memorysize*sizeof(type)); \
   if (ptr == NULL) { \
     status = EGADS_MALLOC; \
     aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 3, "AIM_ALLOC: %s size %zu type %s", #ptr, memorysize, #type); \
     goto cleanup; \
   } \
 }

#define AIM_REALL(ptr, size, type, aimInfo, status) \
 { \
   size_t memorysize = size;\
   ptr = (type *) EG_reall(ptr, memorysize*sizeof(type)); \
   if (ptr == NULL) { \
     status = EGADS_MALLOC; \
     aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 3, "AIM_REALL: %s size %zu type %s", #ptr, memorysize, #type); \
     goto cleanup; \
   } \
 }

#define AIM_STRDUP(ptr, str, aimInfo, status) \
 { \
   if (ptr != NULL) { \
     status = EGADS_MALLOC; \
     aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 1, "AIM_STRDUP: %s != NULL!", #ptr); \
     goto cleanup; \
   } \
   ptr = EG_strdup(str); \
   if (ptr == NULL) { \
     status = EGADS_MALLOC; \
     aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 2, "AIM_STRDUP: %s %s", #ptr, str); \
     goto cleanup; \
   } \
 }

#define AIM_FREE(ptr) { EG_free(ptr); ptr = NULL; }

#define AIM_NOTNULL(ptr, aimInfo, status) \
 { \
   if (ptr == NULL) { \
     status = CAPS_NULLVALUE; \
     aim_status(aimInfo, status, __FILE__, __LINE__, __func__, 1, "%s == NULL!", #ptr); \
     goto cleanup; \
   } \
 }

/*************************************************************************/
/* Prototypes of aim entry points to catch incorrect function signatures */

int
aimInitialize( int inst, /*@null@*/ const char *unitSys, void *aimInfo,
               void **instStore, int *maj, int *min, int *nIn, int *nOut,
               int *nFields, char ***fnames, int **franks, int **fInOut );

int
aimInputs( /*@null@*/ void *instStore, void *aimInfo, int index, char **ainame,
           capsValue *defval );

int
aimOutputs( /*@null@*/ void *instStore, void *aimInfo, int index, char **aoname,
            capsValue *form );

int
aimPreAnalysis( void *instStore, void *aimInfo, /*@null@*/ capsValue *inputs );

int
aimExecute( void *instStore, void *aimInfo, int *state );

int
aimPostAnalysis( void *instStore, void *aimInfo, int restart,
                 /*@null@*/ capsValue *inputs );

void
aimCleanup( /*@only@*/ void *instStore );

int
aimCalcOutput( void *instStore, void *aimInfo, int index, capsValue *val );

int
aimDiscr( char *tname, capsDiscr *discr );

int
aimFreeDiscr( capsDiscr *discr );

int
aimLocateElement( capsDiscr *discr, double *params, double *param,
                  int *bIndex, int *eIndex, double *bary );

int
aimTransfer( capsDiscr *discr, const char *fname, int npts, int rank,
             double *data, char **units );

int
aimInterpolation( capsDiscr *discr, const char *name, int bIndex, int eIndex,
                  double *bary, int rank, double *data, double *result );

int
aimInterpolateBar( capsDiscr *discr, const char *name, int bIndex, int eIndex,
                   double *bary, int rank, double *r_bar, double *d_bar );

int
aimIntegration( capsDiscr *discr, const char *name, int bIndex, int eIndex,
               int rank, double *data, double *result );

int
aimIntegrateBar( capsDiscr *discr, const char *name, int bIndex, int eIndex,
                 int rank, double *r_bar, double *d_bar );

int
aimBackdoor( void *instStore, void *aimInfo, const char *JSONin,
             char **JSONout );

#ifdef __cplusplus
}
#endif

#endif /* AIMUTIL_H */
