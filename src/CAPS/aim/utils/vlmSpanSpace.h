
#include "vlmTypes.h"  // Bring in Vortex Lattice Method structures

#ifdef __cplusplus
extern "C" {
#endif

// Compute auto spanwise panel spacing based on equal spacing on either side of a section
int vlm_autoSpaceSpanPanels(int NspanTotal, int numSection, vlmSectionStruct vlmSection[]);

#ifdef __cplusplus
}
#endif
