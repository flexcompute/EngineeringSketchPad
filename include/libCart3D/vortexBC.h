/*
 *           ----------
 *           vortexBC.h
 *           ----------
 *   $Id: vortexBC.h,v 1.2 2012/02/21 23:39:42 maftosmi Exp $
 *                          2011.02  A.Ning & M.Aftosmis
 */
#ifndef __VORTEX_BC_H__
#define __VORTEX_BC_H__

#include "basicTypes.h"  /* ...need dpoint3 for interface */
#include "stateVector.h" 

/* ---------- MACROS ---------- */

/* ---externBC_setupVortex()---load parameters from a file given by 'fileName'.
                     Uinf     - freestream velocity magnitude
                     returns: - 0 if success, FILE_ERROR  -------------------- */
int externBC_setupVortex(const char fileName[], double Uinf);



/* ---- externBC_getVortexStateV_prim() ----wrapper for vortex query using
                                             standard Cart3D datatyes --- */
tsState externBC_getVortexStateV_prim(dpoint3 pt,
                                      double* UinfPrim, bool y_is_spanwise);

#endif  /*  __VORTEX_BC_H__ */
