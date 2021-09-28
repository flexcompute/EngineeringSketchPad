#include "xddm.h"


extern int writeTrix(const char *fname, /*@null@*/ p_tsXddm p_xddm, int ibody,
                     int nvert, double *xyzs, int nv, /*@null@*/ double **dvar,
                     int ntri, int *tris);

extern int readTrix(const char *fname, const char *tag, int nvert, int dim,
                    double *data);
