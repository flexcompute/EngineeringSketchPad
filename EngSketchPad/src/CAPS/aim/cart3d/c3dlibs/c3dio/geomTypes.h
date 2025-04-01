
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
 * $Id: geomTypes.h,v 1.3 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __GEOMTYPES_H_
#define __GEOMTYPES_H_

#include "basicTypes.h"

#define X  0                   /*      ...dimension names */
#define Y  1
#define Z  2

#define V0 0                   /* ...vertices in tsVertex */
#define V1 1
#define V2 2

typedef int      iquad[4];              /* --  define a indexed quad --------*/
typedef iquad   *p_iquad;

typedef struct VertexStructure tsVertex;/*--  a vertex ----------------------*/
typedef tsVertex * p_tsVertex;

struct VertexStructure {
  float x[DIM];
};

typedef struct DoublePrecisionVertexStructure tsDPVertex;
typedef tsDPVertex * p_tsDPVertex;
/* Double precision vertices - added so we can store and work with double
 * precision geometry
 */
struct DoublePrecisionVertexStructure {
  double x[DIM];
};

typedef struct TriStructure tsTri;      /*--   an annotated triangle --------*/
typedef tsTri * p_tsTri;

struct TriStructure { 
 int      vtx[3];                             /* limits compnums to +-32767  */
 int      Comp:16;                            /* 2 because mtype is an enum  */
 mtype    mark:2, mark2:2;                    /* mark/unmark... etc          */
 unsigned :0;                                 /* pad to word boundary        */
};

#endif /* __GEOMTYPES_H_ */
