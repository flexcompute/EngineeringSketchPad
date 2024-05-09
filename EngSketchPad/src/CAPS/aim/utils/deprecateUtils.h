// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_DEPRECATEUTILS_H_
#define _AIM_UTILS_DEPRECATEUTILS_H_

#include "meshTypes.h"  // Bring in mesh structures
#include "capsTypes.h"  // Bring in CAPS types
#include "cfdTypes.h"   // Bring in cfd structures
#include "feaTypes.h"   // Bring in fea structures
#include "miscTypes.h"  // Bring in miscellanous structures

#ifdef __cplusplus
extern "C" {
#endif

// Change the from the use of capsGroup for mesh sizing to capsMesh
int deprecate_SizingAttr(void *aimInfo,
                         int numTuple,
                         capsTuple meshBCTuple[],
                         const mapAttrToIndexStruct *meshMap,
                         const mapAttrToIndexStruct *groupMap);

#ifdef __cplusplus
}
#endif

#endif //_AIM_UTILS_DEPRECATEUTILS_H_
