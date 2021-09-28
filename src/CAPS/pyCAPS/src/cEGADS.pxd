#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.
cdef extern from "egadsTypes.h":
    ctypedef struct ego:
        pass
    
    cdef int FACE "FACE"
    cdef int EDGE "EDGE"
    cdef int BODY "BODY"
    
    cdef int MODEL "MODEL"
    
    cdef int NODE       "NODE"
    cdef int WIREBODY   "WIREBODY"
    cdef int FACEBODY   "FACEBODY"
    cdef int SHEETBODY  "SHEETBODY"
    cdef int SOLIDBODY  "SOLIDBODY"
    
    cdef int ATTRINT    "ATTRINT"
    cdef int ATTRREAL   "ATTRREAL"
    cdef int ATTRSTRING "ATTRSTRING"
    cdef int ATTRCSYS   "ATTRCSYS"
            
cdef extern from "egads.h":

    #/* memory functions */

    void *EG_alloc( int nbytes )
    
    void *EG_calloc( int nele, int size )
    
    void *EG_reall( void *ptr,
                    int nbytes )
    
    void EG_free( void *pointer )
    
    #/* base-level object functions */
    
    int  EG_open( ego *context );
    int  EG_loadModel( ego context, int bflg, const char *name, 
                       ego *model );
                                
    int  EG_saveModel( const ego model, const char *name )
    
    int  EG_deleteObject( ego object )
    
    int  EG_copyObject( const ego object, void *oform, ego *copy )
    
    int  EG_getContext( ego object, ego *context )

    int  EG_close( ego context );
    
    #/* attribute functions */
    
    int  EG_attributeNum( const ego obj, int *num )
    
    int  EG_attributeGet( const ego obj, int index, const char **name,
                          int *atype, int *len, 
                          const int    **ints,
                          const double **reals, 
                          const char   **str )
    
    int  EG_attributeRet( const ego obj, const char *name, int *atype, 
                          int *len, 
                          const int    **ints,
                          const double **reals,
                          const char   **str )
    
    #/* topology functions */

    int EG_getTopology( const ego topo, ego *geom, int *oclass, 
                        int *type, double *limits, 
                        int *nChildren, ego **children, int **sense )
                                  
    int EG_makeTopology( ego context, ego geom, int oclass, int mtype, 
                         double *limits, int nChildren, ego *children,
                         int *senses, ego *topo )
    
    int EG_getBodyTopos( const ego body, ego src,
                         int oclass, int *ntopo, ego **topos )
    
    int EG_indexBodyTopo( const ego body, const ego src )
                                   
    int EG_getBoundingBox( const ego topo, double *bbox )
    
    #/* tessellation functions */
    
    int  EG_setTessParam( ego context, int iparam, double value,
                          double *oldvalue )
    
    int  EG_makeTessGeom( ego obj, double *params, int *sizes, 
                          ego *tess )
    
    int  EG_getTessGeom( const ego tess, int *sizes, double **xyz )
    
    int  EG_makeTessBody( ego object, double *params, ego *tess )
    
    int  EG_remakeTess( ego tess, int nobj, ego *objs, 
                        double *params )
    
    int  EG_mapTessBody( ego tess, ego body, ego *mapTess )
    
    int  EG_locateTessBody( const ego tess, int npt, const int *ifaces,
                            const double *uv,  int *itri, double *results )
    
    int  EG_getTessEdge( const ego tess, int eIndex, int *len, 
                         const double **xyz, const double **t )
    
    int EG_getTessFace( const ego tess, int fIndex, int *len, 
                         const double **xyz, const double **uv, 
                         const int **ptype, const int **pindex, 
                         int *ntri, const int **tris, 
                         const int **tric )
    
    int  EG_getTessLoops( const ego tess, int fIndex, int *nloop,
                          const int **lIndices )
    
    int  EG_getTessQuads( const ego tess, int *nquad,
                          int **fIndices )
    
    int  EG_makeQuads( ego tess, double *params, int fIndex )
    
    int  EG_getQuads( const ego tess, int fIndex, int *len, 
                      const double **xyz, const double **uv, 
                      const int **ptype, const int **pindex, 
                      int *npatch )
    
    int  EG_getPatch( const ego tess, int fIndex, int patch, 
                      int *nu, int *nv, const int **ipts, 
                      const int **bounds )
    
    int  EG_quadTess( const ego tess, ego *quadTess )
                                   
    int  EG_insertEdgeVerts( ego tess, int eIndex, int vIndex, 
                             int npts, double *t )
    
    int  EG_deleteEdgeVert( ego tess, int eIndex, int vIndex, 
                            int dir )
    
    int  EG_moveEdgeVert( ego tess, int eIndex, int vIndex, 
                          double t )
     
    int  EG_openTessBody( ego tess )
    
    int  EG_initTessBody( ego object, ego *tess )
    
    int  EG_statusTessBody( ego tess, ego *body, int *state, int *np )
    
    int  EG_setTessEdge( ego tess, int eIndex, int len,
                         const double *xyz, const double *t )
    
    int  EG_setTessFace( ego tess, int fIndex, int len,
                         const double *xyz, const double *uv,
                         int ntri, const int *tris )
    
    int  EG_localToGlobal( const ego tess, int index, int local,
                           int *globalID )
    
    int  EG_getGlobal( const ego tess, int globalID, int *ptype,
                       int *pindex,  double *xyz )

cdef extern from "egadsErrors.h":
    cdef int EGADS_CNTXTHRD
    cdef int EGADS_READERR
    cdef int EGADS_TESSTATE
    cdef int EGADS_EXISTS
    cdef int EGADS_ATTRERR
    cdef int EGADS_TOPOCNT
    cdef int EGADS_BADSCALE
    cdef int EGADS_NOTORTHO
    cdef int EGADS_DEGEN
    cdef int EGADS_CONSTERR
    cdef int EGADS_TOPOERR
    cdef int EGADS_GEOMERR
    cdef int EGADS_NOTBODY
    cdef int EGADS_WRITERR
    cdef int EGADS_NOTMODEL
    cdef int EGADS_NOLOAD
    cdef int EGADS_RANGERR
    cdef int EGADS_NOTGEOM
    cdef int EGADS_NOTTESS
    cdef int EGADS_EMPTY
    cdef int EGADS_NOTTOPO
    cdef int EGADS_REFERCE
    cdef int EGADS_NOTXFORM
    cdef int EGADS_NOTCNTX
    cdef int EGADS_MIXCNTX
    cdef int EGADS_NODATA
    cdef int EGADS_NONAME
    cdef int EGADS_INDEXERR
    cdef int EGADS_MALLOC
    cdef int EGADS_NOTOBJ
    cdef int EGADS_NULLOBJ
    cdef int EGADS_NOTFOUND
    cdef int EGADS_SUCCESS "EGADS_SUCCESS"
    cdef int EGADS_OUTSIDE

cdef extern from "prm.h":
    cdef int PRM_NOTCONVERGED  
    cdef int PRM_BADNUMVERTICES
    cdef int PRM_ZEROPIVOT     
    cdef int PRM_NEGATIVEAREAS 
    cdef int PRM_TOLERANCEUNMET
    cdef int PRM_NOGLOBALUV    
    cdef int PRM_BADPARAM      
    cdef int PRM_BADDIVISION   
    cdef int PRM_CANNOTFORMLOOP
    cdef int PRM_WIGGLEDETECTED
    cdef int PRM_INTERNAL      
