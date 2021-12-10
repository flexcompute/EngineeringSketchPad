/*
 *
 */

#ifndef __SOLVERINFO_H_
#define __SOLVERINFO_H_

/*
 * $Id: SolverInfo.h,v 1.3 2006/10/10 16:30:47 aftosmis Exp $
 */

#include "c3d_global.h"
#include "basicTypes.h"
#include "IOinfo.h"
#include "limiters.h"

#ifndef   MAXNUMSTAGES            
#  define MAXNUMSTAGES 10
#endif

/*                        ...these should not be defined in infoStructures.h */
/*                                             (they are not infoStructures) */
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
