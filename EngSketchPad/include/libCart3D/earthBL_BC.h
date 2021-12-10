/*
 *           ----------
 *           earthBL_BC.h
 *           ----------
 *   $Id: earthBL_BC.h,v 1.2 2012/02/20 18:50:54 maftosmi Exp $
 *                                 2012.02 M.Aftosmis
 */
#ifndef __EARTHBL_BC_H__
#define __EARTHBL_BC_H__

#include "basicTypes.h"  /* ...need dpoint3 for interface */
#include "stateVector.h" 

/* ---------- MACROS ---------- */

/* ---externBC_setupEarthBL()---load parameters from a file given by 'fileName'.
                     Uinf     - freestream velocity magnitude
                     returns: - 0 if success, FILE_ERROR  -------------------- */
int externBC_setupEarthBL(const char fileName[]);



/* ---- externBC_getEarthBL_StateV_prim() ----wrapper for vortex query using
                                              standard Cart3D datatyes --- */
tsState externBC_getEarthBL_StateV_prim(dpoint3 pt,  double* UinfPrim, bool y_is_spanwise);

#endif  /*  __EARTHBL_BC_H__ */
