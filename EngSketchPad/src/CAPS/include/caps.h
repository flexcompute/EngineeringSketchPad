#ifndef CAPS_H
#define CAPS_H
/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Function Prototypes
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
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


/* base-level object functions */

__ProtoExt__ void
  caps_revision( int *major, int *minor );

__ProtoExt__ int
  caps_info( capsObj object, char **name, enum capsoType *type,
             enum capssType *subtype, capsObj *link, capsObj *parent,
             capsOwn *last );

__ProtoExt__ int
  caps_size( capsObj object, enum capsoType type, enum capssType stype,
             int *size, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_childByIndex( capsObj object, enum capsoType type,
                     enum capssType stype, int index, capsObj *child );

__ProtoExt__ int
  caps_childByName( capsObj object, enum capsoType typ,
                    enum capssType styp, const char *name, capsObj *child,
                    int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_bodyByIndex( capsObj pobject, int index, ego *body, char **units );

__ProtoExt__ int
  caps_ownerInfo( const capsObj pobject, const capsOwn owner, char **phase,
                  char **pname, char **pID, char **userID, int *nLines,
                  char ***lines, short *datetime, CAPSLONG *sNum );

__ProtoExt__ int
  caps_getHistory( capsObj object, int *nHist, capsOwn **history );

__ProtoExt__ int
  caps_delete( capsObj object );

__ProtoExt__ int
  caps_errorInfo( capsErrs *errs, int eIndex, capsObj *errObj, int *eType,
                  int *nLines, char ***lines );

__ProtoExt__ int
  caps_freeError( /*@only@*/ capsErrs *errs );

__ProtoExt__ int
  caps_printErrors( /*@null@*/ FILE *fp, int nErr, capsErrs *errs );

__ProtoExt__ void
  caps_freeValue( capsValue *value );


/* I/O functions */

__ProtoExt__ int
  caps_writeParameters( const capsObj pobject, char *filename );

__ProtoExt__ int
  caps_readParameters( const capsObj pobject, char *filename );

__ProtoExt__ int
  caps_writeGeometry(capsObj object, int flag, const char *filename,
                     int *nErr, capsErrs **errors);


/* attribute functions */

__ProtoExt__ int
  caps_attrByName( capsObj cobj, char *name, capsObj *attr );

__ProtoExt__ int
  caps_attrByIndex( capsObj cobj, int index, capsObj *attr );

__ProtoExt__ int
  caps_setAttr( capsObj cobj, /*@null@*/ const char *name, capsObj aval );

__ProtoExt__ int
  caps_deleteAttr( capsObj cobj, /*@null@*/ char *name );


/* problem functions */

__ProtoExt__ int
  caps_phaseState( const char *prName, /*@null@*/ const char *phName,
                   int *bitFlag );

__ProtoExt__ int
  caps_journalState( const capsObj pobject );

__ProtoExt__ int
  caps_open( const char *prName, /*@null@*/ const char *phName, int flag,
             /*@null@*/ void *ptr, int outLevel, capsObj *pobject, int *nErr,
             capsErrs **errors );

__ProtoExt__ int
  caps_brokenLink( /*@null@*/ void (*callBack)(capsObj problem, capsObj obj,
                                               enum capstMethod tmeth,
                                               char *name, enum capssType st) );

__ProtoExt__ int
  caps_close( capsObj pobject, int complete, /*@null@*/ const char *phName );

__ProtoExt__ int
  caps_outLevel( capsObj pobject, int outLevel );

__ProtoExt__ int
  caps_getRootPath( capsObj pobject, const char **fullPath );

__ProtoExt__ int
  caps_intentPhrase( capsObj pobject, int nLines,
                     /*@null@*/ const char **lines );

__ProtoExt__ int
  caps_debug( capsObj pobject );


/* analysis functions */

__ProtoExt__ int
  caps_queryAnalysis( capsObj pobj, const char *aname, int *nIn, int *nOut,
                      int *execution );

__ProtoExt__ int
  caps_getBodies( capsObj aobject, int *nBody, ego **bodies,
                  int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_execute( capsObj object, int *state, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_getInput( capsObj pobj, const char *aname, int index, char **ainame,
                 capsValue *defaults );

__ProtoExt__ int
  caps_getOutput( capsObj pobj, const char *aname, int index, char **aoname,
                  capsValue *form );

__ProtoExt__ int
  caps_AIMbackdoor( capsObj aobject, const char *JSONin, char **JSONout );

__ProtoExt__ int
  caps_makeAnalysis( capsObj pobj, const char *anam,
                     /*@null@*/ const char *nam, /*@null@*/ const char *unitSys,
                     /*@null@*/ const char *intents, int *exec,
                     capsObj *aobject, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_dupAnalysis( capsObj from, const char *name,capsObj *aobj );

__ProtoExt__ int
  caps_resetAnalysis( capsObj aobject, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_dirtyAnalysis( capsObj pobj, int *nAobj, capsObj **aobjs );

__ProtoExt__ int
  caps_analysisInfo( capsObj aobject, char **apath, char **unitSys, int *major,
                     int *minor, char **intents, int *nField, char ***fnames,
                     int **ranks, int **fInOut, int *execute, int *status );

__ProtoExt__ int
  caps_preAnalysis( capsObj aobject, int *nErr, capsErrs **errors );


#ifdef ASYNCEXEC
__ProtoExt__ int
  caps_checkAnalysis( capsObj aobject, int *phase, int *nErr, capsErrs **errs );
#endif

__ProtoExt__ int
  caps_system( capsObj aobject, /*@null@*/ const char *rpath,
               const char *command );

__ProtoExt__ int
  caps_postAnalysis( capsObj aobject, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_getTessels( capsObj aobject, int *nTessel, ego **tessels,
                   int *nErr, capsErrs **errors );


/* bound, vertexset and dataset functions */

__ProtoExt__ int
  caps_makeBound( capsObj pobject, int dim, const char *bname, capsObj *bobj );

__ProtoExt__ int
  caps_boundInfo( capsObj object, enum capsState *state, int *dim,
                  double *plims );

__ProtoExt__ int
  caps_closeBound( capsObj bobject );

__ProtoExt__ int
  caps_makeVertexSet( capsObj bobject, /*@null@*/ capsObj aobject,
                      /*@null@*/ const char *vname, capsObj *vobj,
                      int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_vertexSetInfo( capsObj vobject, int *nGpts, int *nDpts, capsObj *bobj,
                      capsObj *aobj );

__ProtoExt__ int
  caps_outputVertexSet( capsObj vobject, const char *filename );

__ProtoExt__ int
  caps_fillUnVertexSet( capsObj vobject, int npts, const double *xyzs );

__ProtoExt__ int
  caps_makeDataSet( capsObj vobject, const char *dname, enum capsfType ftype,
                    int rank, capsObj *dobj, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_dataSetInfo( capsObj dobject, enum capsfType *ftype, capsObj *link,
                    enum capsdMethod *dmeth );

__ProtoExt__ int
  caps_linkDataSet( capsObj link, enum capsdMethod method,
                    capsObj target, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_initDataSet( capsObj dobject, int rank, const double *data,
                    int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_setData( capsObj dobject, int nverts, int rank, const double *data,
                const char *units, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_getData( capsObj dobject, int *npts, int *rank, double **data,
                char **units, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_getDataSets( capsObj bobject, const char *dname, int *nobj,
                    capsObj **dobjs );

__ProtoExt__ int
  caps_getTriangles( capsObj vobject, int *nGtris,            int **gtris,
                                      int *nGsegs, /*@null@*/ int **gsegs,
                                      int *nDtris,            int **dtris,
                                      int *nDsegs, /*@null@*/ int **dsegs );


/* value functions */

__ProtoExt__ int
  caps_getValue( capsObj object, enum capsvType *vtype, int *nrow, int *ncol,
                 const void **data, const int **partial, const char **units,
                 int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_makeValue( capsObj pobject, const char *vname, enum capssType stype,
                  enum capsvType vtype, int nrow, int ncol,
                  const void *data, /*@null@*/ int *partial,
                  /*@null@*/ const char *units, capsObj *vobj );

__ProtoExt__ int
  caps_setValue( capsObj object, enum capsvType vtype, int nrow, int ncol,
                 /*@null@*/ const void *data, /*@null@*/ const int *partial,
                 /*@null@*/ const char *units, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_getLimits( const capsObj object, enum capsvType *vtype,
                  const void **limits, const char **units );

__ProtoExt__ int
  caps_setLimits( capsObj object, enum capsvType vtype,
                  void *limits, const char *units, int *nErr,
                  capsErrs **errors );

__ProtoExt__ int
  caps_getValueProps( capsObj object, int *dim, int *pmtr,
                      enum capsFixed *lfix, enum capsFixed *sfix,
                      enum capsNull *ntype );

__ProtoExt__ int
  caps_setValueProps( capsObj object, int dim, enum capsFixed lfixed,
                      enum capsFixed sfixed, enum capsNull ntype, int *nErr,
                      capsErrs **errors );

__ProtoExt__ int
  caps_convertValue( capsObj value, double inVal, const char *inUnit,
                     double *outVal );

__ProtoExt__ int
  caps_transferValues( capsObj source, enum capstMethod method,
                       capsObj target, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_linkValue( /*@null@*/ capsObj link, enum capstMethod method,
                  capsObj target, int *nErr, capsErrs **errors );

__ProtoExt__ int
  caps_hasDeriv(capsObject *vobj, int *nderiv, char ***names,
                int *nErr, capsErrs **errors);

__ProtoExt__ int
  caps_getDeriv( capsObj value, const char *name, int *len, int *len_wrt,
                 double **deriv, int *nErr, capsErrs **errors );

/* units */

__ProtoExt__ int
  caps_convert( const int count, const char *inUnit,  double  *inVal,
                                 const char *outUnit, double *outVal );

__ProtoExt__ int
  caps_unitParse( const char *unit );

__ProtoExt__ int
  caps_unitConvertible( const char *unit1, const char *unit2 );

__ProtoExt__ int
  caps_unitCompare( const char *unit1, const char *unit2, int *compare );

__ProtoExt__ int
  caps_unitMultiply( const char *inUnits1, const char *inUnits2,
                     char **outUnits );

__ProtoExt__ int
  caps_unitDivide( const char *inUnits1, const char *inUnits2,
                   char **outUnits );

__ProtoExt__ int
  caps_unitInvert( const char *inUnit, char **outUnits );

__ProtoExt__ int
  caps_unitRaise( const char *inUnit, const int power,
                  char **outUnits );

__ProtoExt__ int
  caps_unitOffset( const char *inUnit, const double offset,
                   char **outUnits );

/* others */

__ProtoExt__ void
  caps_externSignal();

__ProtoExt__ void
  caps_rmLock();

__ProtoExt__ void
  caps_printObjects( capsObj pobject, capsObj object, int indent );

#ifdef __cplusplus
}
#endif

#endif
