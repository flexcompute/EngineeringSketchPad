/*
 *
 */

#ifndef __POSTPROCINFO_H_
#define __POSTPROCINFO_H_

/*
 * $Id: PostProcInfo.h,v 1.1 2004/03/22 18:53:14 smurman Exp $
 */

#include "c3d_global.h"

#ifndef   MAXSLICES
#  define MAXSLICES 30
#endif

typedef struct PostprocessingStructure tsPPinfo;
typedef tsPPinfo  *p_tsPPinfo;

struct  PostprocessingStructure {
  int     nSlices[DIM];
  float   a_slicesX[MAXSLICES], a_slicesY[MAXSLICES], a_slicesZ[MAXSLICES];
};

#endif /* __POSTPROCINFO_H_ */
