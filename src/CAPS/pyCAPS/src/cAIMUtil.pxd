#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

cimport cEGADS
cimport cCAPS

cdef extern from "aimUtil.h":

    int aim_file( void *aimInfo, const char *file, char *aimFile )

    int aim_getBodies( void *aimInfo, int *nBody, cEGADS.ego **bodies )
  
    int aim_newGeometry( void *aimInfo )

    int aim_convert( void *aimInfo, const char  *inUnits, double   inValue,
                     const char *outUnits, double *outValue )
    
    int aim_getIndex( void *aimInfo, const char *name, cCAPS.capssType subtype )
  
    int aim_getName( void *aimInfo, int index, cCAPS.capssType subtype,
                     const char **name )
  
    int aim_getData( void *aimInfo, const char *name, cCAPS.capsvType *vtype,
                     int *rank, int *nrow, int *ncol, void **data, char **units )
  
    int aim_getDiscr( void *aimInfo, const char *tname, cCAPS.capsDiscr **discr )
  
    int aim_getDataSet( cCAPS.capsDiscr *discr, const char *dname, cCAPS.capsdMethod *method,
                        int *npts, int *rank, double **data )

    int aim_getBounds( void *aimInfo, int *nTname, char ***tnames )
  
    int aim_setSensitivity( void *aimInfo, const char *GIname, int irow, int icol )
  
    int aim_getSensitivity( void *aimInfo, cEGADS.ego body, int ttype, int index, int *npts,
                            double **dxyz )

    int aim_sensitivity( void *aimInfo, const char *GIname, int irow, int icol,
                         cEGADS.ego tess, int *npts, double **dxyz )
