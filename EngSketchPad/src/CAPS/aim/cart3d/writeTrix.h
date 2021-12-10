#include "xddm.h"
#include "egads.h"

extern int writeTrix(const char *fname, int nbody, ego *tess,
                     /*@null@*/ p_tsXddm p_xddm, int nv, /*@null@*/ double ***dvar);

extern int readTrix(const char *fname, const char *tag, int *dim,
                    double ***data_out);
