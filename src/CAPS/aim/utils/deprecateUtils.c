// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "egads.h"
#include "aimUtil.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef WIN32
#define strcasecmp  stricmp
#endif

#include "deprecateUtils.h"
#include "miscUtils.h"

// Change from the use of capsGroup for mesh sizing to capsMesh
int deprecate_SizingAttr(void *aimInfo,
                         int numTuple,
                         capsTuple meshBCTuple[],
                         mapAttrToIndexStruct *meshMap,
                         mapAttrToIndexStruct *groupMap){
    int status;
    int i, j, attrIndex;

    for (i = 0; i < numTuple; i++) {

        status = get_mapAttrToIndexIndex(meshMap, meshBCTuple[i].name, &attrIndex);
        if (status == CAPS_NOTFOUND) {
            status = get_mapAttrToIndexIndex(groupMap, meshBCTuple[i].name, &attrIndex);
            if (status == CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "DEPRACATED: The capsGroup attribute (capsGroup=%s) is no longer used to specify mesh sizing parameters. Please "
                                   "use capsMesh instead.", meshBCTuple[i].name);
            } else {
                AIM_ERROR(aimInfo, "No attribute capsMesh == '%s'.", meshBCTuple[i].name);
                AIM_ADDLINE(aimInfo, "------------------------------");
                AIM_ADDLINE(aimInfo, "Available capsMesh attributes:");
                for (j = 0; j < meshMap->numAttribute; j++)
                  AIM_ADDLINE(aimInfo, "%s", meshMap->attributeName[j]);
                AIM_ADDLINE(aimInfo, "------------------------------");
            }
          
            return CAPS_BADVALUE;
        }
    }

    return CAPS_SUCCESS;
}
