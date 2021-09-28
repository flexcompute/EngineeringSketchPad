#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

#cdef extern from "/home/dursch/Desktop/CAPS/EngSketchPad/src/wvServer/libwebsockets.h":
#    ctypedef struct libwebsocket:
#        pass
     
#cdef extern from "webViewer.c":
#    #void browserMessage(libwebsocket *wsi, char *text, int len)
#    void browserMessage(void *wsi, char *text, int len)

ctypedef wvContext *wvObj
    
cdef extern from "wsss.h":

    cdef int BUFLEN "BUFLEN"

    #/* Graphic Primitive Types 
    cdef int WV_POINT "WV_POINT"
    cdef int WV_LINE "WV_LINE"
    cdef int WV_TRIANGLE "WV_TRIANGLE"

    #/* Plotting Attribute Bits 
    cdef int WV_ON "WV_ON"
    cdef int WV_TRANSPARENT "WV_TRANSPARENT"
    cdef int WV_SHADING "WV_SHADING"
    cdef int WV_ORIENTATION "WV_ORIENTATION"
    cdef int WV_POINTS "WV_POINTS"
    cdef int WV_LINES "WV_LINES"
    
    #/* VBO & single data Bits 
    cdef int WV_VERTICES "WV_VERTICES"
    cdef int WV_INDICES "WV_INDICES"
    cdef int WV_COLORS "WV_COLORS"
    cdef int WV_NORMALS "WV_NORMALS"
    
    cdef int WV_PINDICES "WV_PINDICES"
    cdef int WV_LINDICES "WV_LINDICES"
    cdef int WV_PCOLOR "WV_PCOLOR"
    cdef int WV_LCOLOR "WV_LCOLOR"
    cdef int WV_BCOLOR "WV_BCOLOR"
    
    cdef int WV_DELETE "WV_DELETE"
    cdef int WV_DONE "WV_DONE"
    
    #/* Data Types 
    cdef int WV_UINT8  "WV_UINT8"
    cdef int WV_UINT16 "WV_UINT16"
    cdef int WV_INT32  "WV_INT32"
    cdef int WV_REAL32 "WV_REAL32"
    cdef int WV_REAL64 "WV_REAL64"

    ctypedef struct wvData:
        int    dataType         # VBO type 
        int    dataLen          # length of data 
        void  *dataPtr          # pointer to data 
        float  data[3]          # size = 1 data (not in ptr) 
    
    ctypedef struct wvStripe:
        int            nsVerts   #      number of vertices 
        int            nsIndices #      the number of top-level indices 
        int            nlIndices #      the number of line indices 
        int            npIndices #      the number of point indices 
        int            *gIndices #      global indices for the stripe -- NULL if 1 
        float          *vertices #      the stripes set of vertices (3XnsVerts) 
        float          *normals  #      the stripes suite of normals (3XnsVerts) 
        unsigned char  *colors   #      the stripes suite of colors (3XnsVerts) 
        unsigned short *sIndice2 #      the top-level stripe indices 
        unsigned short *lIndice2 #      the line indices 
        unsigned short *pIndice2 #      the point indices 
    
    ctypedef struct wvGPrim:
        int            gtype     #  graphic primitive type 
        int            updateFlg #  update flag 
        int            nStripe   #  number of stripes 
        int            attrs     #  initial plotting attribute bits 
        int            nVerts    #  number of vertices 
        int            nIndex    #  the number of top-level indices 
        int            nlIndex   #  the number of line indices 
        int            npIndex   #  the number of point indices 
        int            nameLen   #  length of name (modulo 4) 
        float          pSize     #  point size 
        float          pColor[3] #  point color 
        float          lWidth    #  line width 
        float          lColor[3] #  line color 
        float          fColor[3] #  triangle foreground color 
        float          bColor[3] #  triangle background color 
        float          normal[3] #  constant normal -- planar surfaces only 
        char           *name     #  gPrim name (nameLen is length+) 
        float          *vertices #  the complete set of vertices (3XnVerts) 
        float          *normals  #  the complete suite of normals (3XnVerts) 
        unsigned char  *colors   #  the complete suite of colors (3XnVerts) 
        int            *indices  #  the complete suite of top-level indices 
        int            *lIndices #  the complete suite of line indices 
        int            *pIndices #  the complete suite of point indices 
        wvStripe       *stripes  #  stripes 
    
    ctypedef struct wvContext:
         #wvCB    callback #     optional call-back 
        int     handShake  # larger scale handshaking 
        int     ioAccess   # IO currently has access 
        int     dataAccess # data routines have access 
        int     bias       # vaule subtracted from all indices 
        int     keepItems  # 0 - Items are transient, 1 - user control 
        float   fov        # angle field of view for perspective 
        float   zNear      # Z near value for perspective 
        float   zFar       # Z far value for perspective 
        float   eye[3]     # eye coordinates for "camera" 
        float   center[3]  # lookat center coordinates 
        float   up[3]      # the up direction -- usually Y 
        int     nGPrim     # number of graphics primitives 
        int     mGPrim     # number of malloc'd primitives 
        int     cleanAll   # flag for complete gPrim delete 
        int     sent       # send flag 
        int     nColor     # # Colors # -n sent, 0 delete, +n to be sent 
        int     titleLen   # titleLen (bytes) # neg - send delete key 
        
        float   *colors    # NULL or colors to be sent 
        float   lims[2]    # key limits 
        
        char    *title     # title for the key 
        
        int     tnWidth    # thumbnail width 
        int     tnHeight   # thumbnail height 
        
        char    *thumbnail # the thumbnail image (rgba) 
        
        wvGPrim *gPrims    # the graphics primitives 
        
cdef extern from "wsserver.h":
    
    void wv_freeItem(wvData *item) # Added 
    
    int  wv_startServer( int port, char *interface, 
                                   char *cert_path, 
                                   char *key_path, int opts,
                                   wvContext *WVcontext )
                                  
    int  wv_statusServer( int index )
  
    int  wv_nClientServer( int index )

    wvContext *wv_createContext( int bias,   float fov,  float zNear, 
                                 float zFar, float *eye, float *center, 
                                 float *up )
    
    void wv_destroyContext( wvContext **context )

    void wv_printGPrim( wvContext *cntxt, int index )

    int  wv_indexGPrim( wvContext *cntxt, char *name )

    int  wv_addGPrim( wvContext *cntxt, char *name, int gtype, 
                      int attrs, int nItems, wvData *items )
    
    void wv_removeGPrim( wvContext *cntxt, int index )
                               
    void wv_removeAll( wvContext *cntxt )
                               
    int  wv_modGPrim( wvContext *cntxt, int index, int nItems, wvData *items )
    
    int  wv_handShake( wvContext *cntxt )
                               
    int  wv_addArrowHeads( wvContext *cntxt, int index, float size, 
                           int nHeads, int *heads )
  
    int  wv_setKey( wvContext *context, int nCol, float *colors,
                    float beg, float end, char *title )

    int  wv_setData( int type, int len, void *data, int VBOtype, 
                     wvData *dstruct )

    # #ifndef STANDALONE
    void wv_sendText( void *wsi, char *text )
    # #endif

    void wv_broadcastText( char *text )
                              
    void wv_adjustVerts( wvData *dstruct, float *focus )

    void wv_cleanupServers() #( void )
