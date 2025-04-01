
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
 * $Id: SolverInfo.h,v 1.4 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __SOLVERINFO_H_
#define __SOLVERINFO_H_

#include "c3d_global.h"
#include "basicTypes.h"
#include "IOinfo.h"
#include "limiters.h"

#ifndef   MAXNUMSTAGES            
#  define MAXNUMSTAGES 10
#endif

typedef enum { NONE=-1, SCALAR, JACOBI} pretype;  /*-- precondition type --  */
                                                  /*- flux function type --  */
typedef enum { VANLEER, VLHANEL, COLELLA, HLLC, HCUSP, VLMOD} fftype;

typedef struct SolverInfoStructure tsSinfo;
typedef tsSinfo *p_tsSinfo;

struct SolverInfoStructure {
  int        nStage;                   /*  # of stages in Runge-Kutta Scheme */
  double     a_stageCoef[MAXNUMSTAGES];/*        Array of stage coefficients */
  bool       a_gradEval[MAXNUMSTAGES]; /* yes/no eval gradient at each stage */
  double     cfl;                      /*                         CFL number */
  double     rampUp;                   /* factor used to ramp up CFL         */
  double     rampedCfl;                /* need to remember cfl from last step*/
  limtype    limiter;                  /*   slope limiter for reconstruction */
  int        freezeAfter;              /*  freeze limitersr for convergence  */
  fftype     fluxFunction;             /*     ...inviscid flux function flag */
  pretype    pc;                       /*   scalar or matrix preconditioner? */
  int        bboxBCs[2*DIM];           /* BCs on  domain BBox [LoHiLoHiLoHi] */
  bool       doSubcellSurf;            /*  keep subcell tris or agglomerate? */
  bool       first_order;
  p_tsIOinfo p_fileInfo;               
}; 

#endif /* __SOLVERINFO_H_ */
