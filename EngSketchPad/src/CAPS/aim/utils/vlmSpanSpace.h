// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_VLMSPANSPACE_H_
#define _AIM_UTILS_VLMSPANSPACE_H_

#include "vlmTypes.h"  // Bring in Vortex Lattice Method structures

#ifdef __cplusplus
extern "C" {
#endif

// Compute auto spanwise panel spacing based on equal spacing on either side of a section
int vlm_autoSpaceSpanPanels(void *aimInfo, int NspanTotal, int numSection, vlmSectionStruct vlmSection[]);

#ifdef __cplusplus
}
#endif

#endif // _AIM_UTILS_VLMSPANSPACE_H_
