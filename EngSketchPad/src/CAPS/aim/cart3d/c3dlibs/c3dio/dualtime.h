
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 *  $Id: dualtime.h,v 1.8 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

/**
 * data structures and macros for the time-dependent dual-time scheme
 */

#ifndef __DUALTIME_H_
#define __DUALTIME_H_

#include <stdlib.h>
#include "basicTypes.h"

#define MAX_NUM_DT_MODS 64

typedef struct dtControlStruct tsDTcontrol;
typedef  tsDTcontrol *p_tsDTcontrol;

struct   dtControlStruct {
  int timestep[MAX_NUM_DT_MODS];   /* arrays to store contents of
                                      timestep.cntl */
  int action[MAX_NUM_DT_MODS];     /* and trigger DT mods */

  bool isDTstored;                 /* check if user-supplied DT mods
                                      are stored. If false, running
                                      -autoDT will use calc_timesteps
                                      and stored maxWaveSpeed */
  int iNextDTmod;                  /* DT mod timesteps should be
                                      monotonically increasing, so
                                      instead of searching for
                                      timestep in array, just keep
                                      track of the index of next
                                      timestep when we change DT */
};

/* NOTE: tsTDinfo struct DOES get written into unsteady checkpoint
 * files.  This WILL break unsteady checkpoints (although technically
 * should still work with Un and Unm1) WILL "break" adapt since it
 * uses offset(tsTDinfo) and error metric written after TDinfo.  Just
 * need to recompile libCart3D and everything that uses this header.
 */
typedef struct TDInfoStructure tsTDinfo;
typedef tsTDinfo *p_tsTDinfo;

struct TDInfoStructure {
  double       time_p;                   /* the total physical time            */
  int          nTDSteps;                 /* the number of physical timesteps   */
  double       dt_p;                     /* the physical timestep (!=0 for T-D)*/
  int          curTDStep;                /* the current physical timestep      */
  float        a_TDMethod[3];            /* the 2-step schemes' coefficients   */
  double       waveSpeed;         /* characteristic wave speed at current time */
  int          action;            /* current state,   use to trigger dT change */
  tsDTcontrol  dtControl;         /* struct to store user-specified autoDT info*/
};

/*
 * Global macros
 */
#define GET_PHYSICAL_DT(P_S) (P_S)->p_tdInfo->dt_p

#define IS_TIME_DEPENDENT(P_S) (0.0 < (P_S)->p_tdInfo->dt_p)

#define IS_TWOSTEP_SCHEME(THETA, XI, PHI) ((0.0 != (XI)) || (0.0 != (PHI)))

#define CONSTANT_DT   0  /* ...define states  TDinfo.status */
#define INCREASE_DT  +1
#define DECREASE_DT  -1


/*
 * public methods
 */

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* __DUALTIME_H_ */
