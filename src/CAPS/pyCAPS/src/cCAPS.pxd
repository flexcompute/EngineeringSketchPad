#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

cimport cEGADS

cdef extern from "capsTypes.h":

    cdef int CAPSMAGIC "CAPSMAGIC"
    
    # Define - need to check long long
    ctypedef unsigned long CAPSLONG
    
    # Structures
    ctypedef struct capsObject:
        int     type #;   /* object type */
        void    *blind #; /* blind pointer to object data */
        char    *name #; /* object name */
    
    ctypedef capsObject * capsObj
    
    ctypedef struct capsOwn:
        pass
    
    ctypedef struct capsErrs:
        pass
    
    ctypedef struct capsProblem:
        cEGADS.ego        context #; /* the EGADS context object */
    
    ctypedef struct capsAnalysis:
        int         nBody  #; /* number of Bodies for this Analysis */
        cEGADS.ego *bodies #; /* the bodies */    

    ctypedef struct capsTuple:
        char *name  #; /* the name */
        char *value #; /* the value for the pair */
    
    cdef union capsValueUnion:
        int        integer   # /* single int -- length == 1 */
        int        *integers # /* multiple ints */
        double     real      # /* single double -- length == 1 */
        double     *reals    # /* mutiple doubles */
        char       *string   # /* character string (no single char) */
        capsTuple  *tuple    # /* tuple (no single tuple) */
        capsObject *object   # /* single object */
        capsObject **objects # /* multiple objects */
    
    ctypedef struct capsValue: # structure for CAPS object -- VALUE
        int type   # /* value type */
        int length # /* number of values */
        int dim    # /* the dimension */
        int nrow   # /* number of rows */
        int ncol   # /* the number of columns */
        int lfixed # /* length is fixed */
        int sfixed # /* shape is fixed */
        int nullVal# /* NULL handling */
        capsValueUnion vals
        
    ctypedef struct capsDiscr: # defines a discretized collection of Elements
        pass
    
    # Enums
    cdef enum capsoType:
        BODIES, ATTRIBUTES, UNUSED, PROBLEM, VALUE, ANALYSIS, BOUND, VERTEXSET, DATASET
    
    cdef enum capssType:
        NONE, STATIC, PARAMETRIC, GEOMETRYIN, GEOMETRYOUT, BRANCH, PARAMETER, USER, 
        ANALYSISIN, ANALYSISOUT, CONNECTED, UNCONNECTED
    
    cdef enum capsBoolean:
        false, true
    
    cdef enum capsvType:
        Boolean, Integer, Double, String, Tuple, Value
    
    cdef enum capsFixed:
        Change, Fixed
    
    cdef enum capsNull:
        NotAllowed, NotNull, IsNull
    
    cdef enum capstMethod:
        Copy
    
    cdef enum capsdMethod:
        Analysis, Interpolate, Conserve
    
    cdef enum capsState:
        MultipleError, Open, Empty, Single, Multiple

cdef extern from "caps.h":
    
    #/* base-level object functions */
    int caps_info( const capsObj object, char **name, capsoType *type, 
                   capssType *subtype, capsObj *link, capsObj *parent, 
                   capsOwn *last )
    
    int caps_size( const capsObj object,  capsoType type,  capssType stype, 
                   int *size )
    
    int caps_childByIndex( const capsObj object,  capsoType type, 
                           capssType stype, int index, capsObj *child )
    
    int caps_childByName( const capsObj object,  capsoType typ,
                          capssType styp, const char *name, capsObj *child )
    
    int caps_bodyByIndex( const capsObj pobject, int index, cEGADS.ego *body, char **units )
    
    int caps_ownerInfo( const capsOwn owner, char **pname, char **pID, char **userID, 
                        short *datetime, CAPSLONG *sNum )
    
    int caps_delete( capsObj object )
    
    int caps_errorInfo( capsErrs *errs, int eIndex, capsObj *errObj, int *nLines,
                        char ***lines )
    
    int caps_freeError( capsErrs *errs )
    
    void caps_freeValue( capsValue *value )

    #/* attribute functions */
    int caps_attrByName( capsObj cobj, char *name, capsObj *attr )
    
    int caps_attrByIndex( capsObj cobj, int index, capsObj *attr )
    
    int caps_setAttr( capsObj cobj, const char *name, capsObj aval )
    
    int caps_deleteAttr( capsObj cobj, char *name )
    
    #/* problem functions */
    int caps_open( const char *filename, const char *pname, capsObj *pobject )
  
    int caps_save( capsObj pobject, const char *filename )
  
    int caps_close( capsObj pobject )
    
    int caps_outLevel( capsObj pobject, int outLevel )
 
    #/* analysis functions */
    
    int caps_queryAnalysis( capsObj pobj, const char *aname, int *nIn, int *nOut,
                            int *execution )

    int caps_getBodies( const capsObject *aobject, int *nBody, cEGADS.ego **bodies )

    int caps_getInput( capsObj pobj, const char *aname, int index, char **ainame,
                       capsValue *defaults )
  
    int caps_getOutput( capsObj pobj, const char *aname, int index, char **aoname,
                        capsValue *form )
  
    int caps_AIMbackdoor( const capsObject *aobject, const char *JSONin, char **JSONout )
                    
    int caps_load( capsObj pobj, const char *anam, const char *apath, const char *unitSys, 
                   char *intent, int nparent, capsObj *parents, capsObj *aobj )
    
    int caps_dupAnalysis( capsObj fromObj, const char *apath, int nparent, 
                          capsObj *parents,capsObj *aobj )
    
    int caps_dirtyAnalysis( capsObj pobj, int *nAobj, capsObj **aobjs )
    
    int caps_analysisInfo( const capsObj aobject, char **apath, char **unitSys,
                           char **intent, int *nparent, capsObj **parents, 
                           int *nField, char ***fnames, int **ranks, 
                           int *execution, int *status)
    
    int caps_preAnalysis( capsObj aobject, int *nErr, capsErrs **errors )
    
    int caps_postAnalysis( capsObj aobject, capsOwn current, int *nErr,
                           capsErrs **errors )
    
    #/* bound, vertexset and dataset functions */
    int caps_makeBound( capsObj pobject, int dim, const char *bname, capsObj *bobj )
    
    int caps_boundInfo( const capsObj object,  capsState *state, int *dim, 
                        double *plims )
    
    int caps_closeBound( capsObj bobject )
    
    int caps_fillVertexSets( capsObj bobject, int *nErr, capsErrs **errors )

    int caps_makeVertexSet( capsObj bobject, capsObj aobject, 
                            const char *vname, capsObj *vobj )
    
    int caps_vertexSetInfo( const capsObj vobject, int *nGpts, int *nDpts,
                            capsObj *bobj, capsObj *aobj )
    
    int caps_fillUnVertexSet( capsObj vobject, int npts, const double *xyzs )
    
    int caps_makeDataSet( capsObj dobject, const char *dname,  capsdMethod meth,
                          int rank, capsObj *dobj )
    
    int caps_initDataSet( capsObj dobject, int rank, double *startup )

    int caps_setData( capsObj dobject, int nverts, int rank, const double *data,
                      const char *units )
    
    int caps_getData( capsObj dobject, int *npts, int *rank, double **data,
                      char **units )
    
    int caps_getHistory( const capsObj dobject, capsObj *vobject, int *nHist,
                         capsOwn **history )
    
    int caps_getDataSets( const capsObj bobject, const char *dname, int *nobj,
                          capsObj **dobjs )
    
    int caps_triangulate( const capsObj vobject, int *nGtris, int **gtris,
                          int *nDtris, int **dtris )
    
    #/* value functions */
    int caps_getValue( capsObj object,  capsvType *type, int *vlen,
                       const void **data, const char **units, int *nErr,
                       capsErrs **errors )
    
    int caps_makeValue( capsObj pobject, const char *vname,  capssType stype,
                        capsvType vtype, int nrow, int ncol, const void *data,
                        const char *units, capsObj *vobj )
    
    int caps_setValue( capsObj object, int nrow, int ncol, const void *data )
    
    int caps_getLimits( const capsObj object, const void **limits )
    
    int caps_setLimits( capsObj object, void *limits )
    
    int caps_getValueShape( const capsObj object, int *dim,  capsFixed *lfix, 
                            capsFixed *sfix,  capsNull *ntype,
                            int *nrow, int *ncol )
    
    int caps_setValueShape( capsObj object, int dim,  capsFixed lfixed,
                            capsFixed sfixed,  capsNull ntype )
    
    int caps_convert( const capsObj object, const char *inUnits, double inVal, 
                                            const char *outUnits, double *outVal )
    
    int caps_transferValues( capsObj source,  capstMethod method, 
                             capsObj target, int *nErr, capsErrs **errors )
    
    int caps_linkValue( capsObj link,  capstMethod method,
                          capsObj target )
    
cdef extern from "capsErrors.h":
    cdef int CAPS_SUCCESS "CAPS_SUCCESS"
    cdef int CAPS_BADRANK
    cdef int CAPS_BADDSETNAME
    cdef int CAPS_NOTFOUND
    cdef int CAPS_BADINDEX 
    cdef int CAPS_NOTCHANGED
    cdef int CAPS_BADTYPE  
    cdef int CAPS_NULLVALUE
    cdef int CAPS_NULLNAME
    cdef int CAPS_NULLOBJ  
    cdef int CAPS_BADOBJECT "CAPS_BADOBJECT" 
    cdef int CAPS_BADVALUE
    cdef int CAPS_PARAMBNDERR
    cdef int CAPS_NOTCONNECT
    cdef int CAPS_NOTPARMTRIC
    cdef int CAPS_READONLYERR
    cdef int CAPS_FIXEDLEN
    cdef int CAPS_BADNAME
    cdef int CAPS_BADMETHOD
    cdef int CAPS_CIRCULARLINK
    cdef int CAPS_UNITERR
    cdef int CAPS_NULLBLIND
    cdef int CAPS_SHAPEERR
    cdef int CAPS_LINKERR
    cdef int CAPS_MISMATCH
    cdef int CAPS_NOTPROBLEM
    cdef int CAPS_RANGEERR
    cdef int CAPS_DIRTY
    cdef int CAPS_HIERARCHERR
    cdef int CAPS_STATEERR
    cdef int CAPS_SOURCEERR
    cdef int CAPS_EXISTS
    cdef int CAPS_IOERR
    cdef int CAPS_DIRERR
    cdef int CAPS_NOTIMPLEMENT
    cdef int CAPS_EXECERR
    cdef int CAPS_CLEAN
    cdef int CAPS_BADINTENT "CAPS_BADINTENT"
    cdef int CAPS_BADINIT      
    cdef int CAPS_NOTNEEDED    
    cdef int CAPS_NOSENSITVTY
