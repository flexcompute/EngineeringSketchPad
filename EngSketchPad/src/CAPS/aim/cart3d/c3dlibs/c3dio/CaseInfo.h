
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
 *  $Id: CaseInfo.h,v 1.6 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __CASEINFO_H_
#define __CASEINFO_H_

#include "c3d_global.h"
#include "basicTypes.h"
#include "IOinfo.h"
#include "stateVector.h"

typedef struct CaseInfoStructure tsCinfo;
typedef tsCinfo *p_tsCinfo;

struct CaseInfoStructure {
  double           Minf;     /*                       freestream Mach number */
  double          alpha;     /* angle of attack  - measured in the X-Y plane */
  double           beta;     /* sideslip angle   - measured in the Z-X plane */
  double           roll;     /* body roll angle  - measured in the Y-Z plane */
  double         rhoinf;     /*                           freestream density */
  double           pinf;     /*                          freestream pressure */
  double           Hinf;     /*                    freestream total enthalpy */
  double           qinf;     /*                  freestream dynamic pressure */
  tsState          Uinf;     /*    ...collected state vector for free stream */
  tsState          UinfPrim; /*    ...primitive var state vector-free stream */
  double          gamma;     /*                    ..ratio of specific heats */
  double           ainf;     /*                ...free stream speed of sound */
  double         Froude;     /*     ...Froude number -- if using body forces */
  double   gravity[DIM];     /* ...gravity unit vect -- if using body forces */
  bool          restart;     /*   ... Restart into any # of partitions (T/F) */
  bool        restartMP;     /*   ... Restart for  static partitioning (T/F) */
  p_tsIOinfo p_fileInfo;
  char  cmdLine[STRING_LEN]; /*  <--2014.11.17 changed DIM from FILENAME_LEN */
};

#endif /* __CASEINFO_H_ */
