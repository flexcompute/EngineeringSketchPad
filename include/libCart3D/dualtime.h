/*
 *  $Id: dualtime.h,v 1.4 2014/04/18 14:01:10 maftosmi Exp $
 *
 *  data structures and macros for the time-dependent
 *  dual-time scheme
 */

#ifndef __DUALTIME_H_
#define __DUALTIME_H_

#include <stdlib.h>

/*
 * Global data types
 */
typedef struct TDInfoStructure tsTDinfo;
typedef tsTDinfo *p_tsTDinfo;

struct TDInfoStructure {
  double     time_p;                   /* the total physical time            */
  int        nTDSteps;                 /* the number of physical timesteps   */
  double     dt_p;                     /* the physical timestep (!=0 for T-D)*/
  int        curTDStep;                /* the current physical timestep      */
  float      a_TDMethod[3];            /* the 2-step schemes' coefficients   */ 
};


/*
 * Global macros
 */
#define GET_PHYSICAL_DT(P_S) (P_S)->p_tdInfo->dt_p

#define IS_TIME_DEPENDENT(P_S) (0.0 < (P_S)->p_tdInfo->dt_p)

#define IS_TWOSTEP_SCHEME(THETA, XI, PHI) ((0.0 != (XI)) || (0.0 != (PHI)))

/*
 * public methods
 */

/*
 * this retrieves the current scheme from the tdInfo structure.  this
 * should be used rather than querying the structure directly so that
 * special cases such as start-up, or failure can be handled.  if any
 * of the parameter are set NULL, then they will be ignored and left
 * empty.
 */
void getTDMethod(const p_tsTDinfo p_tdInfo, double *theta, double *xi, double *phi);

/*
 * append the physical time to a filename when performing I/O
 * for time-dependent problems.
 */
void appendTDStamp(const double time, char *file_name, const size_t length);

#endif /* __DUALTIME_H_ */
