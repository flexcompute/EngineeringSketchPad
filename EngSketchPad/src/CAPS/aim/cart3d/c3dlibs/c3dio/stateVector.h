
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
 * $Id: stateVector.h,v 1.6 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

/**
 * state vector of conservated quantities and perfect gas relations
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

/* for accessing fields in Power BC state vector                       */
#define BACKPRESS               0       /* back pressure               */
#define NORMVEL                 1       /* normal velocity             */
#define TOTTEMP                 0       /* total temperature           */
#define MFR                     1       /* mass flow rate   != TOTTEMP */
#define TOTPRESS                1       /* total pressure   != TOTTEMP */

#define VEL2(A) ((A)[XVEL]*(A)[XVEL]+(A)[YVEL]*(A)[YVEL] + (A)[ZVEL]*(A)[ZVEL])
#define MOM2(A) ((A)[XMOM]*(A)[XMOM]+(A)[YMOM]*(A)[YMOM] + (A)[ZMOM]*(A)[ZMOM])

                        /* !!! left out factor of 2 in the linearization !!! */
#define LIN_VEL2(A,B) ((A)[XVEL]*(B)[XVEL] + (A)[YVEL]*(B)[YVEL] + (A)[ZVEL]*(B)[ZVEL])

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

#endif /* __STATEVECTOR_H_ */
