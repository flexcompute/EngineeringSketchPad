/*
 *  $Id: CaseInfo.h,v 1.5 2014/12/22 16:15:44 maftosmi Exp $
 */

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
