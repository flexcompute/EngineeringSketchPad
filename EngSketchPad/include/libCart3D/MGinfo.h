/*
 *
 */

#ifndef __MGINFO_H_
#define __MGINFO_H_

/*
 * $Id: MGinfo.h,v 1.2 2004/03/31 23:49:49 berger Exp $
 */

#include "basicTypes.h"
#include "IOinfo.h"

typedef struct MultiGridInfoStructure tsMGinfo;
typedef tsMGinfo *p_tsMGinfo;

struct MultiGridInfoStructure {
  int        nMGlev;             /*                  Num of MultiGrid levels */
  int        finestSoFar;        /*     for FMG, finest level up to so far   */
  int        cycleType;          /*   1 = "V", 2 = "W": often called 'gamma' */
  int        nPre;               /*   # of  pre-smoothing steps              */
  int        nPost;              /*   # of post-smoothing steps(0 = sawtooth)*/
  int        maxCycles;          /*                    max MG cycles allowed */
  int        nCycles;            /*     # of cycles with current finest Grid */
  int        firstCycle;         /*   1 if new case; old maxCycles if restart*/
  bool       gs;                 /* levels for grid seq. (init.) only, no mg */
};

#endif /* __MGINFO_H_ */
