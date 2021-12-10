/*
 * $Id: infoStructures.h,v 1.5 2010/07/30 16:05:31 maftosmi Exp $
 *
 */
/* -----------------
   infoStructures.h
   -----------------
                          ...all of cart3Ds "InfoStructures"
			     get defined here as do pointers
			     to all of them

     these include:
       Structure             short name   pointer
       ---------            -----------  -------
     1. IOinfoStructure        tsIOinfo  p_tsIOinfo
     2. MultiGridInfoStructure tsMGinfo  p_tsMGinfo
     3. SolverInfoStructure    tsSinfo   p_tsSinfo 
     4. GridInfoStructure      tsGinfo   p_tsGinfo 
     5. CaseInfoStructure      tsCinfo   p_tsCinfo 
 */

#ifndef __INFOSTRUCTURES_H_
#define __INFOSTRUCTURES_H_

/* the infoStructures file as a large repository is now deprecated.  
 * This file is kept for backwards compatibility, but should no longer
 * be used.  The individual files should be used to get specific data
 * structures.  Some of these files themselves will likely be moved out
 * of libCart3D eventually. */

#include "IOinfo.h"
#include "MGinfo.h"
#include "SolverInfo.h"
#include "GridInfo.h"
#include "CaseInfo.h"
#include "IOinfo.h"
#include "ConvInfo.h"
#include "PostProcInfo.h"

#endif /* __INFOSTRUCTURES_H_ */
