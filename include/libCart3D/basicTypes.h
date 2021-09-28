/*
 * $Id: basicTypes.h,v 1.2 2005/06/15 15:50:20 smurman Exp $
 *
 */
/* ------------
   basictypes.h
   ------------
                   ...elementary Type  definitions for general programming
*/

#ifndef __BASICTYPES_H_
#define __BASICTYPES_H_

#include "c3d_global.h"

#if defined(FALSE)
#  undef FALSE
#endif
#if defined(TRUE)
#  undef TRUE
#endif

#ifdef __cplusplus
# define TRUE true
# define FALSE false
#endif

/* we need to wrap the bool def or the swig pre-processor will
 * break.  if we need to get at bool types through swig wrappers
 * we will have to write special adaptors for each case -SM- */
#ifndef SWIG
#ifndef __cplusplus
typedef enum { FALSE, TRUE }       bool; /* --------- define Boolean Type -- */
#endif
#endif /* SWIG */

typedef enum { UNMARKED, MARKED } mtype; /* --------- define Mark Type ----  */
typedef unsigned char byte;              /* --------- define byte -----------*/
typedef double   dpoint3[DIM];           /* -- define double 3D point -------*/
typedef dpoint3 *p_dpoint3;

/* ---- */

#endif /* __BASICTYPES_H_ */
