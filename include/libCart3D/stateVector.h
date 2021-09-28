/*
 * $Id: stateVector.h,v 1.2 2004/11/05 23:52:44 aftosmis Exp $
 *
 */
/* -------------
   statevector.h
   -------------
                                 ...state vector of conservated quantities and
				    perfect gas relations
*/

#ifndef __STATEVECTOR_H_
#define __STATEVECTOR_H_

#include "c3d_global.h"   /* to get DIM */

/* the state vector may alternatively store a flow state in conservative
 * or primitive variables
 */
#define NSTATES                 5
#define RHO                     0       /*..keep                       */
#define XMOM                    1       /*      this                   */
#define YMOM                    2       /*          order              */
#define ZMOM                    3       /*               unchanged     */
#define RHOE                    4       /*                !!!!!!       */
#define PRESS                   4       /*                             */
#define XVEL                 XMOM       /*                             */
#define YVEL                 YMOM       /*                             */
#define ZVEL                 ZMOM       /*                             */

#define VEL2(A) ((A)[XVEL]*(A)[XVEL]+(A)[YVEL]*(A)[YVEL] + (A)[ZVEL]*(A)[ZVEL])
#define MOM2(A) ((A)[XMOM]*(A)[XMOM]+(A)[YMOM]*(A)[YMOM] + (A)[ZMOM]*(A)[ZMOM])

struct StateVectorStructure {        /* ...The vector of dependent variables */
  double  v[NSTATES];                /* 0 = density                          */
};                                   /* 1 = density * Xvel                   */
                                     /* 2 = density * Yvel                   */
                                     /* 3 = density * Zvel                   */
                                     /* 4 = density * E    <--called "rhoE"  */
/*                                      where E = total energy per unit mass */
/* Perfect Gas Relations: rhoH = rhoE + p                                    */
/* ---------------------  e = c_v T, h = c_p T                               */
/*                        E = e + q^2/2                                      */
/*                        H = h + q^2/2                                      */
/*                     rhoH = rhoE + p                                       */
/*                        p = (gam-1)(rhoE - rho*q^2/2)                      */
/*                    or  p = ((gam-1)/gam) (rhoH - rho*q^2/2)               */
/*                                                                           */
/* ------------------------------------------------------------------------- */


typedef struct StateVectorStructure  tsState;
typedef tsState *p_tsState;

typedef tsState  tsState3[DIM];                          /* type of gradient */
typedef tsState3 *p_tsState3;
			    
/*
  ----------
*/

#endif /* __STATEVECTOR_H_ */
