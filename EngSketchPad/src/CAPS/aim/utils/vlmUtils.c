// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Vortex lattice analysis related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include <string.h>
#include <math.h>
#include "aimUtil.h"    // Bring in AIM utilities
#include "vlmTypes.h"   // Bring in Vortex Lattice Method structures
#include "vlmUtils.h"
#include "miscUtils.h"  // Bring in miscellaneous utilities

#ifdef WIN32
#define strcasecmp  stricmp
#endif

extern int EG_isPlanar(const ego object);

#define PI        3.1415926535897931159979635
#define NINT(A)   (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define MAX(A,B)  (((A) < (B)) ? (B) : (A))

#define CROSS(a,b,c)      a[0] = ((b)[1]*(c)[2]) - ((b)[2]*(c)[1]);\
                          a[1] = ((b)[2]*(c)[0]) - ((b)[0]*(c)[2]);\
                          a[2] = ((b)[0]*(c)[1]) - ((b)[1]*(c)[0])

// Tolerance for checking if a dot product between airfoil section normals is zero
#define DOTTOL 1.e-7

// Fill vlmSurface in a vlmSurfaceStruct format with vortex lattice information
// from an incoming surfaceTuple
int get_vlmSurface(void *aimInfo,
                   int numTuple,
                   capsTuple surfaceTuple[],
                   mapAttrToIndexStruct *attrMap,
                   double Cspace,
                   int *numVLMSurface,
                   vlmSurfaceStruct *vlmSurface[]) {

    /*! \page vlmSurface Vortex Lattice Surface
     * Structure for the Vortex Lattice Surface tuple  = ("Name of Surface", "Value").
     * "Name of surface defines the name of the surface in which the data should be applied.
     *  The "Value" can either be a JSON String dictionary (see Section \ref jsonStringVLMSurface)
     *  or a single string keyword string (see Section \ref keyStringVLMSurface).
     */

    int status; //Function return

    int i, groupIndex, attrIndex; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;

    // Clear out vlmSurface is it has already been allocated
    if (*vlmSurface != NULL) {
        for (i = 0; i < (*numVLMSurface); i++) {
            status = destroy_vlmSurfaceStruct(&(*vlmSurface)[i]);
            AIM_STATUS(aimInfo, status);
        }
        AIM_FREE(*vlmSurface);
    }

    printf("Getting vortex lattice surface data\n");

    if (numTuple <= 0){
        printf("\tNumber of VLM Surface tuples is %d\n", numTuple);
        return CAPS_NOTFOUND;
    }

    *numVLMSurface = numTuple;
    *vlmSurface = (vlmSurfaceStruct *) EG_alloc((*numVLMSurface) * sizeof(vlmSurfaceStruct));
    if (*vlmSurface == NULL) {
        *numVLMSurface = 0;
        return EGADS_MALLOC;
    }

    // Initiate vlmSurfaces
    for (i = 0; i < (*numVLMSurface); i++) {
        status = initiate_vlmSurfaceStruct(&(*vlmSurface)[i]);
        AIM_STATUS(aimInfo, status);

        // set default Cspace
        (*vlmSurface)[i].Cspace = Cspace;
    }

    for (i = 0; i < (*numVLMSurface); i++) {

        printf("\tVLM surface name - %s\n", surfaceTuple[i].name);

        // Copy surface name
        (*vlmSurface)[i].name = EG_strdup(surfaceTuple[i].name);

        // Do we have a json string?
        if (strncmp(surfaceTuple[i].value, "{", 1) == 0) {

            //printf("JSON String - %s\n",surfaceTuple[i].value);

            /*! \page vlmSurface
             * \section jsonStringVLMSurface JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary (eg. "Value" = {"numChord": 5, "spaceChord": 1.0, "numSpan": 10, "spaceSpan": 0.5})
             * the following keywords ( = default values) may be used:
             *
             * \if (AVL)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <em>capsGroup</em> names used to define the surface (e.g. "Name1" or ["Name1","Name2",...].
             *  If no groupName variable is provided an attempted will be made to use the tuple name instead;
             * </ul>
             * \else
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <em>capsGroup</em> names used to define the surface (e.g. "Name1" or ["Name1","Name2",...].
             *  If no groupName variable is provided an attempted will be made to use the tuple name instead;
             * </ul>
             * \endif
             *
             */

            // Get surface variables
            keyWord = "groupName";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                AIM_STATUS(aimInfo, status);

                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;

                // Determine how many capsGroups go into making the surface
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndex);

                    if (status == CAPS_NOTFOUND) {
                        printf("groupName name %s not found in map of capsGroups attributes!", groupName[groupIndex]);
                        continue;
                    }
                    AIM_STATUS(aimInfo, status);

                    (*vlmSurface)[i].numAttr += 1;
                    if ((*vlmSurface)[i].numAttr == 1) {
                        (*vlmSurface)[i].attrIndex = (int *) EG_alloc(1*sizeof(int));

                    } else{
                        (*vlmSurface)[i].attrIndex = (int *) EG_reall((*vlmSurface)[i].attrIndex,
                                (*vlmSurface)[i].numAttr*sizeof(int));
                    }

                    if ((*vlmSurface)[i].attrIndex == NULL) {
                        status = EGADS_MALLOC;
                        goto cleanup;
                    }

                    (*vlmSurface)[i].attrIndex[(*vlmSurface)[i].numAttr-1] = attrIndex;
                }

                status = string_freeArray(numGroupName, &groupName);
                AIM_STATUS(aimInfo, status);
                groupName = NULL;
                numGroupName = 0;

            } else {
                printf("\tNo \"groupName\" variable provided or no matches found, going to use tuple name\n");
            }

            if ((*vlmSurface)[i].numAttr == 0) {

                status = get_mapAttrToIndexIndex(attrMap, (const char *) (*vlmSurface)[i].name, &attrIndex);
                if (status == CAPS_NOTFOUND) {
                    AIM_ERROR(aimInfo, "Tuple name %s not found in map of capsGroups attributes!", (*vlmSurface)[i].name);
                    goto cleanup;
                }

                (*vlmSurface)[i].numAttr += 1;
                if ((*vlmSurface)[i].numAttr == 1) {
                    (*vlmSurface)[i].attrIndex = (int *) EG_alloc(1*sizeof(int));
                } else{
                    (*vlmSurface)[i].attrIndex = (int *) EG_reall((*vlmSurface)[i].attrIndex,
                    (*vlmSurface)[i].numAttr*sizeof(int));
                }

                if ((*vlmSurface)[i].attrIndex == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                (*vlmSurface)[i].attrIndex[(*vlmSurface)[i].numAttr-1] = attrIndex;
            }

            AIM_FREE(keyValue);

            /*! \page vlmSurface
             * \if (AVL)
             * <ul>
             * <li> <B>noKeyword = "(no default)"</B> </li> <br>
             *      "No" type. Options: NOWAKE, NOALBE, NOLOAD.
             * </ul>
             * \endif
             */

            keyWord = "noKeyword";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                if      (strcasecmp(keyValue, "\"NOWAKE\"") == 0) (*vlmSurface)[i].nowake = (int) true;
                else if (strcasecmp(keyValue, "\"NOALBE\"") == 0) (*vlmSurface)[i].noalbe = (int) true;
                else if (strcasecmp(keyValue, "\"NOLOAD\"") == 0) (*vlmSurface)[i].noload = (int) true;
                else {

                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for VLM Section tuple %s, current options are "
                            "\" NOWAKE, NOALBE, or  NOLOAD\"", keyWord, keyValue, surfaceTuple[i].name);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * <ul>
             * <li> <B>numChord = 10</B> </li> <br>
             *  The number of chordwise horseshoe vortices placed on the surface.
             * </ul>
             *
             */
            keyWord = "numChord";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &(*vlmSurface)[i].Nchord);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * \if ( AVL )
             * <ul>
             * <li> <B>spaceChord = 1.0</B> </li> <br>
             *  The chordwise vortex spacing parameter.<br>
             *  Note: The local spanwise spacing may be overridden using the <em>vlmSspace</em> BODY attribute on a section.
             * </ul>
             * \else
             * <ul>
             * <li> <B>spaceChord = 0.0</B> </li> <br>
             *  The chordwise vortex spacing parameter.
             * </ul>
             * \endif
             *
             */
            keyWord = "spaceChord";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*vlmSurface)[i].Cspace);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /* Check for lingering numSpan in old scripts */
            keyWord = "numSpan";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_FREE(keyValue);
                AIM_ERROR(aimInfo, "numSpan is depricated.");
                AIM_ADDLINE(aimInfo, "Please use numSpanTotal or numSpanPerSection instead.");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            /*! \page vlmSurface
             * <ul>
             * <li> <B>numSpanTotal = 0</B> </li> <br>
             *  Total number of spanwise horseshoe vortices placed on the surface.
             *  The vorticies are 'evenly' distributed across sections to minimize jumps in spacings.
             *  numpSpanPerSection must be zero if this is set.<br>
             *  Note: The local spanwise count may be overridden using the <em>vlmNumSpan</em> BODY attribute on a section.
             * </ul>
             *
             */
            keyWord = "numSpanTotal";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &(*vlmSurface)[i].NspanTotal);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * <ul>
             * <li> <B>numSpanPerSection = 0</B> </li> <br>
             *  The number of spanwise horseshoe vortices placed on each section the surface.
             *  The total number of spanwise vorticies are (numSection-1)*numSpanPerSection.
             *  The vorticies are 'evenly' distributed across sections to minimize jumps in spacings.
             *  numSpanTotal must be zero if this is set.<br>
             *  Note: The local spanwise count may be overridden using the <em>vlmNumSpan</em> BODY attribute on a section.
             * </ul>
             *
             */
            keyWord = "numSpanPerSection";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &(*vlmSurface)[i].NspanSection);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            if ((*vlmSurface)[i].NspanTotal != 0 && (*vlmSurface)[i].NspanSection != 0) {
                AIM_ERROR(aimInfo, "Only one of numSpanTotal and numSpanPerSection can be non-zero!");
                AIM_ADDLINE(aimInfo, "       numSpanTotal      = %d", (*vlmSurface)[i].NspanTotal);
                AIM_ADDLINE(aimInfo, "       numSpanPerSection = %d", (*vlmSurface)[i].NspanSection);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            /*! \page vlmSurface
             * \if ( ASTROS || NASTRAN || ZAERO )
             * <ul>
             * <li> <B>spaceSpan = 0.0</B> </li> <br>
             *  The spanwise vortex spacing parameter.<br>
             *  The local spacing may be overridden using the vlmSspace BODY attribute on a section.
             * </ul>
             * \endif
             */
            keyWord = "spaceSpan";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*vlmSurface)[i].Sspace);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * \if ( AVL )
             * <ul>
             * <li> <B>yMirror = False</B> </li> <br>
             *  Mirror surface about the y-direction.
             * </ul>
             * \endif
             */
            keyWord = "yMirror";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toBoolean(keyValue, &(*vlmSurface)[i].iYdup);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * \if ( AVL )
             * <ul>
             * <li> <B>component = 0</B> </li> <br>
             *  A positive number allows multiple input SURFACEs to be grouped together into a composite
             *  virtual surface.  Application examples are:
             *  - A wing component made up of a wing SURFACE and a winglet SURFACE
             *  - A T-tail component made up of horizontal and vertical tail SURFACEs.
             * </ul>
             * \endif
             *
             */
            keyWord = "component";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &(*vlmSurface)[i].compon);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * \if ( AVL )
             * <ul>
             * <li> <B>CLaf = 0.0</B> </li> <br>
             * This scales the effective dcl/da of the section airfoil as follows: <br>
             *  <br>
             *  dcl/da  =  2 pi ( 1 + CLaf t/c ) <br>
             * <br>
             * where t/c is the airfoil's thickness/chord ratio.
             * </ul>
             * \endif
             */
            keyWord = "CLaf";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*vlmSurface)[i].CLaf);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmSurface
             * \if ( ASTROS )
             * <ul>
             * <li> <B>surfaceType = "Wing"</B> </li> <br>
             *  Type of aerodynamic surface being described: "Wing", "Canard", "Tail".
             * </ul>
             *  \endif
             */
            keyWord = "surfaceType";
            status = search_jsonDictionary( surfaceTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                (*vlmSurface)[i].surfaceType = string_removeQuotation(keyValue);

                AIM_FREE(keyValue);
            } else {
                (*vlmSurface)[i].surfaceType = EG_strdup("Wing");
            }

        } else {

            /*! \page vlmSurface
             * \section keyStringVLMSurface Single Value String
             *
             * If "Value" is a single string the following options maybe used:
             * - (NONE Currently)
             *
             */
            AIM_ERROR(aimInfo, "No current defaults for get_vlmSurface, tuple value must be a JSON string");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting vortex lattice surface data\n");

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(keyValue);

    if (numGroupName != 0 && groupName != NULL){
        (void) string_freeArray(numGroupName, &groupName);
    }
    return status;
}


// Fill vlmControl in a vlmControlStruct format with vortex lattice information
// from an incoming controlTuple
int get_vlmControl(void *aimInfo,
                   int numTuple,
                   capsTuple controlTuple[],
                   const char *angleUnit,
                   int *numVLMControl,
                   vlmControlStruct *vlmControl[]) {

    /*! \page vlmControl Vortex Lattice Control Surface
     * Structure for the Vortex Lattice Control Surface tuple  = ("Name of Control Surface", "Value").
     * "Name of control surface defines the name of the control surface in which the data should be applied.
     *  The "Value" must be a JSON String dictionary (see Section \ref jsonStringVLMSection).
     */

    int status; //Function return

    int i; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    // Clear out vlmControl if it has already been allocated
    if (*vlmControl != NULL) {
        for (i = 0; i < (*numVLMControl); i++) {
            status = destroy_vlmControlStruct(&(*vlmControl)[i]);
            if (status != CAPS_SUCCESS) printf("destroy_vlmControlStruct status = %d\n", status);
        }
        AIM_FREE(*vlmControl);
    }

    printf("Getting vortex lattice control surface data\n");

    if (numTuple <= 0){
        printf("\tNumber of VLM Surface tuples is %d\n", numTuple);
        return CAPS_NOTFOUND;
    }

    AIM_ALLOC(*vlmControl, numTuple, vlmControlStruct, aimInfo, status);
    *numVLMControl = numTuple;

    // Initiate vlmControl
    for (i = 0; i < (*numVLMControl); i++) {
        status = initiate_vlmControlStruct(&(*vlmControl)[i]);
        AIM_STATUS(aimInfo, status);
    }

    for (i = 0; i < (*numVLMControl); i++) {

        printf("\tVLM control surface name - %s\n", controlTuple[i].name);

        // Copy surface name
        AIM_STRDUP((*vlmControl)[i].name, controlTuple[i].name, aimInfo, status);

        // Do we have a json string?
        if (strncmp(controlTuple[i].value, "{", 1) == 0) {

            //printf("JSON String - %s\n",surfaceTuple[i].value);

            /*! \page vlmControl
             * \section jsonStringVLMSection JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary (e.g. "Value" = {"deflectionAngle": 10.0}) the following keywords ( = default values) may be used:
             *
             */

            /*! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>deflectionAngle = 0.0</B> </li> <br>
             *      Deflection angle of the control surface.
             * </ul>
             * \endif
             */
            keyWord = "deflectionAngle";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                if (angleUnit != NULL)
                  status = string_toDoubleUnits(aimInfo, keyValue, angleUnit, &(*vlmControl)[i].deflectionAngle);
                else
                  status = string_toDouble(keyValue, &(*vlmControl)[i].deflectionAngle);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);

                AIM_FREE(keyValue);
            }

            /*
            ! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>percentChord = 0.0</B> </li> <br>
             *      Percentage along the airfoil chord the control surface's hinge line resides.
             * </ul>
             * \endif

            keyWord = "percentChord";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*vlmControl)[i].percentChord);
                AIM_STATUS(aimInfo, status);

                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;
            }
             */

            /*! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>leOrTe = (no default) </B> </li> <br>
             *      Is the control surface at the leading or trailing edge?
             *        - For leOrTe = 0, the control surface is always between the leading edge and the hinge line.
             *        - For leOrTe = 1, the control surface is always between the trailing edge and the hinge line.
             *
             *      The default behavior when leOrTe is not specified: If the hinge location percentage along
             *      the airfoil chord is < 50% a leading edge slat is assumed, while >= 50% indicates a
             *      trailing edge flap.
             * </ul>
             * \endif
             */
            keyWord = "leOrTe";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                (*vlmControl)[i].leOrTeOverride = (int) true;
                status = string_toInteger(keyValue, &(*vlmControl)[i].leOrTe);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>controlGain = 1.0</B> </li> <br>
             *      Control deflection gain, units:  degrees deflection / control variable
             * </ul>
             * \endif
             */
            keyWord = "controlGain";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*vlmControl)[i].controlGain);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>hingeLine = [0.0 0.0 0.0]</B> </li> <br>
             *      Alternative vector giving hinge axis about which surface rotates
             * </ul>
             * \endif
             */
            keyWord = "hingeLine";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue, 3, (*vlmControl)[i].xyzHingeVec);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmControl
             * \if (AVL)
             * <ul>
             * <li> <B>deflectionDup = 0 </B> </li> <br>
             *      Sign of deflection for duplicated surface
             * </ul>
             * \endif
             */
            keyWord = "deflectionDup";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &(*vlmControl)[i].deflectionDup);
                AIM_STATUS(aimInfo, status);

                AIM_FREE(keyValue);
            }

            /*! \page vlmControl
             * \if (ASTROS)
             * <ul>
             * <li> <B>surfaceSymmetry = "SYM" </B> </li> <br>
             *      The surface type of symmetry.
             * </ul>
             * \endif
             */
            keyWord = "surfaceSymmetry";
            status = search_jsonDictionary( controlTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                (*vlmControl)[i].surfaceSymmetry = string_removeQuotation(keyValue);
                if ((*vlmControl)[i].surfaceSymmetry == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;
            }




        } else {

            /*! \page vlmControl
             * \section keyStringVLMControl Single Value String
             *
             * If "Value" is a single string, the following options maybe used:
             * - (NONE Currently)
             *
             */
            AIM_ERROR(aimInfo, "No current defaults for get_vlmControl, tuple value must be a JSON string");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting vortex lattice control surface data\n");

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Premature exit in get_vlmControl, status = %d\n", status);

    AIM_FREE(keyValue);

    return status;
}


// Initiate (0 out all values and NULL all pointers) of a control in the vlmcontrol structure format
int initiate_vlmControlStruct(vlmControlStruct *control) {

    control->name = NULL; // Control surface name

    control->deflectionAngle = 0.0; // Deflection angle of the control surface
    control->controlGain = 1.0; //Control deflection gain, units:  degrees deflection / control variable

    control->percentChord = 0.0; // Percentage along chord

    control->xyzHinge[0] = 0.0; // xyz location of the hinge
    control->xyzHinge[1] = 0.0;
    control->xyzHinge[2] = 0.0;

    control->xyzHingeVec[0] = 0.0; // Vector of hinge line at xyzHinge
    control->xyzHingeVec[1] = 0.0;
    control->xyzHingeVec[2] = 0.0;


    control->leOrTeOverride = (int) false; // Does the user want to override the geometry set value?
    control->leOrTe = 0; // Leading = 0 or trailing > 0 edge control surface

    control->deflectionDup = 0; // Sign of deflection for duplicated surface

    control->surfaceSymmetry = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of a control in the vlmcontrol structure format
int destroy_vlmControlStruct(vlmControlStruct *control) {

    if (control->name != NULL) EG_free(control->name);
    control->name = NULL; // Control surface name

    control->deflectionAngle = 0.0; // Deflection angle of the control surface

    control->controlGain = 1.0; //Control deflection gain, units:  degrees deflection / control variable

    control->percentChord = 0.0; // Percentage along chord

    control->xyzHinge[0] = 0.0; // xyz location of the hinge
    control->xyzHinge[1] = 0.0;
    control->xyzHinge[2] = 0.0;

    control->xyzHingeVec[0] = 0.0; // Vector of hinge line at xyzHinge
    control->xyzHingeVec[1] = 0.0;
    control->xyzHingeVec[2] = 0.0;

    control->leOrTeOverride = (int) false; // Does the user want to override the geometry set value?
    control->leOrTe = 0; // Leading = 0 or trailing > 0 edge control surface

    control->deflectionDup = 0; // Sign of deflection for duplicated surface

    if (control->surfaceSymmetry != NULL) EG_free(control->surfaceSymmetry);
    control->surfaceSymmetry = NULL;

    return CAPS_SUCCESS;
}


// Initiate (0 out all values and NULL all pointers) of a section in the vlmSection structure format
int initiate_vlmSectionStruct(vlmSectionStruct *section) {

    section->name = NULL;

    section->ebody = NULL;
    section->type  = vlmSecUnknown;
    section->sectionIndex = 0;

    section->xyzLE[0] = 0.0;
    section->xyzLE[1] = 0.0;
    section->xyzLE[2] = 0.0;
    section->nodeIndexLE = 0;

    section->xyzTE[0] = 0.0;
    section->xyzTE[1] = 0.0;
    section->xyzTE[2] = 0.0;
    section->teObj = NULL;
    section->teClass = 0;

    section->chord = 0.0;
    section->ainc = 0.0;

    section->normal[0] = 0.0;
    section->normal[1] = 0.0;
    section->normal[2] = 0.0;

    section->Nspan = 0;
    section->Sspace = 0.0;
    section->Sset = (int)false;

    section->numControl = 0;
    section->vlmControl = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of a section in the vlmSection structure format
int destroy_vlmSectionStruct(vlmSectionStruct *section) {

    int status; // Function return status
    int i; // Indexing

    EG_free(section->name);
    section->name = NULL;

    EG_deleteObject(section->ebody);
    section->ebody = NULL;
    section->type  = vlmSecUnknown;
    section->sectionIndex = 0;

    section->xyzLE[0] = 0.0;
    section->xyzLE[1] = 0.0;
    section->xyzLE[2] = 0.0;
    section->nodeIndexLE = 0;

    section->xyzTE[0] = 0.0;
    section->xyzTE[1] = 0.0;
    section->xyzTE[2] = 0.0;
    section->teObj = NULL;
    section->teClass = 0;

    section->chord = 0.0;
    section->ainc = 0.0;

    section->normal[0] = 0.0;
    section->normal[1] = 0.0;
    section->normal[2] = 0.0;

    section->Nspan = 0;
    section->Sspace = 0.0;
    section->Sset = (int)false;

    if (section->vlmControl != NULL) {

        for (i = 0; i < section->numControl; i++) {
            status = destroy_vlmControlStruct(&section->vlmControl[i]);
            if (status != CAPS_SUCCESS) printf("destroy_vlmControlStruct %d\n", status);
        }

        EG_free(section->vlmControl);
    }

    section->vlmControl = NULL;
    section->numControl = 0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of a surface in the vlmSurface structure format
int initiate_vlmSurfaceStruct(vlmSurfaceStruct  *surface) {

    // Surface information
    surface->name   = NULL;

    surface->numAttr = 0;    // Number of capsGroup/attributes used to define a given surface
    surface->attrIndex = NULL; // List of attribute map integers that correspond given capsGroups

    surface->Nchord = 10;
    surface->Cspace = 0.0;

    surface->NspanTotal = 0;
    surface->NspanSection = 0;
    surface->Sspace = 0.0;

    surface->CLaf = 0.0;

    surface->nowake = (int) false;
    surface->noalbe = (int) false;;
    surface->noload = (int) false;;
    surface->compon = 0;
    surface->iYdup  = (int) false;

    // Section storage
    surface->numSection = 0;
    surface->vlmSection = NULL;

    surface->surfaceType = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of a surface in the vlmSurface structure format
int destroy_vlmSurfaceStruct(vlmSurfaceStruct  *surface) {

    int status; // Function return status
    int i; // Indexing

    // Surface information
    EG_free(surface->name);
    surface->name = NULL;

    surface->numAttr = 0;    // Number of capsGroup/attributes used to define a given surface
    EG_free(surface->attrIndex);
    surface->attrIndex = NULL; // List of attribute map integers that correspond given capsGroups

    surface->Nchord = 0;
    surface->Cspace = 0.0;

    surface->NspanTotal = 0;
    surface->NspanSection = 0;
    surface->Sspace = 0.0;

    surface->CLaf = 0.0;

    surface->nowake = (int) false;
    surface->noalbe = (int) false;;
    surface->noload = (int) false;;
    surface->compon = 0;
    surface->iYdup  = (int) false;

    // Section storage
    if (surface->vlmSection != NULL) {
        for (i = 0; i < surface->numSection; i++ ) {
            status = destroy_vlmSectionStruct(&surface->vlmSection[i]);
            if (status != CAPS_SUCCESS) printf("destroy_vlmSectionStruct status = %d", status);
        }

        EG_free(surface->vlmSection);
    }

    surface->vlmSection = NULL;
    surface->numSection = 0;

    EG_free(surface->surfaceType);
    surface->surfaceType = NULL;

    return CAPS_SUCCESS;
}

// Populate vlmSurface-section control surfaces from geometry attributes, modify control properties based on
// incoming vlmControl structures
int get_ControlSurface(void *aimInfo,
                       int numControl,
                       vlmControlStruct vlmControl[],
                       vlmSurfaceStruct *vlmSurface) {

    int status = 0; // Function status return

    int section, control, attr, index; // Indexing

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;

    const char *attributeKey = "vlmControl", *aName = NULL;
    const char *attrName = NULL;

    int numAttr = 0; // Number of attributes

    // Control related variables
    double chordPercent, chordVec[3]; //chordLength

    /* Remove previous control structures */
    for (section = 0; section < vlmSurface->numSection; section++) {
      for (control = 0; control < vlmSurface->vlmSection[section].numControl; control++) {
        status = destroy_vlmControlStruct(&vlmSurface->vlmSection[section].vlmControl[control]);
        AIM_STATUS(aimInfo, status);
      }
    }


    for (section = 0; section < vlmSurface->numSection; section++) {

        vlmSurface->vlmSection[section].numControl = 0;

        status = EG_attributeNum(vlmSurface->vlmSection[section].ebody, &numAttr);
        AIM_STATUS(aimInfo, status);

        // Control attributes
        for (attr = 0; attr < numAttr; attr++) {

            status = EG_attributeGet(vlmSurface->vlmSection[section].ebody,
                                     attr+1,
                                     &aName,
                                     &atype,
                                     &alen,
                                     &ints,
                                     &reals, &string);

            if (status != EGADS_SUCCESS) continue;
            if (atype  != ATTRREAL)      continue;

            if (strncmp(aName, attributeKey, strlen(attributeKey)) != 0) continue;

            if (alen == 0) {
                AIM_ERROR(aimInfo, "%s should be followed by a single value corresponding to the flap location "
                                   "as a function of the chord. 0 - 1 (fraction - %% / 100), 1-100 (%%)", aName);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            if (reals[0] > 100) {
                AIM_ERROR(aimInfo, "%s value (%f) must be less than 100", aName, reals[0]);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            //printf("Attribute name = %s\n", aName);

            if (aName[strlen(attributeKey)] == '\0')

                attrName = "Flap";

            else if (aName[strlen(attributeKey)] == '_') {

                attrName = &aName[strlen(attributeKey)+1];

            } else {

                attrName = &aName[strlen(attributeKey)];
            }

            //printf("AttrName = %s\n", attrName);

            vlmSurface->vlmSection[section].numControl += 1;
            AIM_REALL(vlmSurface->vlmSection[section].vlmControl,
                      vlmSurface->vlmSection[section].numControl, vlmControlStruct, aimInfo, status);

            index = vlmSurface->vlmSection[section].numControl-1; // Make copy to shorten the following lines of code

            status = initiate_vlmControlStruct(&vlmSurface->vlmSection[section].vlmControl[index]);
            AIM_STATUS(aimInfo, status);

            // Get name of control surface
            AIM_FREE(vlmSurface->vlmSection[section].vlmControl[index].name);

            AIM_STRDUP(vlmSurface->vlmSection[section].vlmControl[index].name, attrName, aimInfo, status);

            // Loop through control surfaces from input Tuple and see if defaults have be augmented
            for (control = 0; control < numControl; control++) {

                if (strcasecmp(vlmControl[control].name, attrName) != 0) continue;

                status = copy_vlmControlStruct(&vlmControl[control],
                                               &vlmSurface->vlmSection[section].vlmControl[index]);
                AIM_STATUS(aimInfo, status);
                break;
            }

            if (control == numControl) {
                printf("Warning: Control %s not found in controls tuple! Only defaults will be used.\n", attrName);
            }

            // Get percent of chord from attribute
            if (reals[0] < 0.0) {
                printf("Warning: Percent chord must > 0, converting to a positive number.\n");
                vlmSurface->vlmSection[section].vlmControl[index].percentChord = -1.0* reals[0];

            } else {

                vlmSurface->vlmSection[section].vlmControl[index].percentChord = reals[0];
            }

            // Was value given as a percentage or fraction
            if (vlmSurface->vlmSection[section].vlmControl[index].percentChord  >= 1.0) {
                vlmSurface->vlmSection[section].vlmControl[index].percentChord  = vlmSurface->vlmSection[section].vlmControl[index].percentChord / 100;
            }

            if (vlmSurface->vlmSection[section].vlmControl[index].leOrTeOverride == (int) false) {

                if (vlmSurface->vlmSection[section].vlmControl[index].percentChord  < 0.5) {
                    vlmSurface->vlmSection[section].vlmControl[index].leOrTe = 0;
                } else {
                    vlmSurface->vlmSection[section].vlmControl[index].leOrTe = 1;
                }
            }

            // Get xyz of hinge location
            chordVec[0] = vlmSurface->vlmSection[section].xyzTE[0] - vlmSurface->vlmSection[section].xyzLE[0];
            chordVec[1] = vlmSurface->vlmSection[section].xyzTE[1] - vlmSurface->vlmSection[section].xyzLE[1];
            chordVec[2] = vlmSurface->vlmSection[section].xyzTE[2] - vlmSurface->vlmSection[section].xyzLE[2];

            chordPercent = vlmSurface->vlmSection[section].vlmControl[index].percentChord;

            vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[0] = chordPercent*chordVec[0] +
                                                                            vlmSurface->vlmSection[section].xyzLE[0];

            vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[1] = chordPercent*chordVec[1] +
                                                                            vlmSurface->vlmSection[section].xyzLE[1];

            vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[2] = chordPercent*chordVec[2] +
                                                                            vlmSurface->vlmSection[section].xyzLE[2];

            /*
            printf("\nCheck hinge line \n");

            printf("ChordLength = %f\n", chordLength);
            printf("ChordPercent = %f\n", chordPercent);
            printf("LeadingEdge = %f %f %f\n", vlmSurface->vlmSection[section].xyzLE[0],
                                               vlmSurface->vlmSection[section].xyzLE[1],
                                               vlmSurface->vlmSection[section].xyzLE[2]);

            printf("TrailingEdge = %f %f %f\n", vlmSurface->vlmSection[section].xyzTE[0],
                                                vlmSurface->vlmSection[section].xyzTE[1],
                                                vlmSurface->vlmSection[section].xyzTE[2]);

            printf("Hinge location = %f %f %f\n", vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[0],
                                                  vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[1],
                                                  vlmSurface->vlmSection[section].vlmControl[index].xyzHinge[2]);
             */

            /*// Set the hinge vector to the plane normal - TOBE ADDED
            vlmSurface->vlmSection[section].vlmControl[index].xyzHingeVec[0] =
            vlmSurface->vlmSection[section].vlmControl[index].xyzHingeVec[1] =
            vlmSurface->vlmSection[section].vlmControl[index].xyzHingeVec[2] =
             */
        }
    }

    status = CAPS_SUCCESS;

cleanup:

    return status;
}


// Make a copy of vlmControlStruct (it is assumed controlOut has already been initialized)
int copy_vlmControlStruct(vlmControlStruct *controlIn, vlmControlStruct *controlOut) {

    int status;

    if (controlIn == NULL) return CAPS_NULLVALUE;
    if (controlOut == NULL) return CAPS_NULLVALUE;

    status = destroy_vlmControlStruct(controlOut);
    if (status != CAPS_SUCCESS) return status;

    if (controlIn->name != NULL) {
        controlOut->name = EG_strdup(controlIn->name);
        if (controlOut->name == NULL) return EGADS_MALLOC;
    }

    controlOut->deflectionAngle = controlIn->deflectionAngle; // Deflection angle of the control surface
    controlOut->controlGain = controlIn->controlGain; //Control deflection gain, units:  degrees deflection / control variable


    controlOut->percentChord = controlIn->percentChord; // Percentage along chord

    controlOut->xyzHinge[0] = controlIn->xyzHinge[0]; // xyz location of the hinge
    controlOut->xyzHinge[1] = controlIn->xyzHinge[1];
    controlOut->xyzHinge[2] = controlIn->xyzHinge[2];

    controlOut->xyzHingeVec[0] = controlIn->xyzHingeVec[0]; // xyz location of the hinge
    controlOut->xyzHingeVec[1] = controlIn->xyzHingeVec[1];
    controlOut->xyzHingeVec[2] = controlIn->xyzHingeVec[2];


    controlOut->leOrTeOverride = controlIn->leOrTeOverride; // Does the user want to override the geometry set value?
    controlOut->leOrTe = controlIn->leOrTe; // Leading = 0 or trailing > 0 edge control surface

    controlOut->deflectionDup = controlIn->deflectionDup; // Sign of deflection for duplicated surface

    if (controlIn->surfaceSymmetry != NULL) {
        controlOut->surfaceSymmetry = EG_strdup(controlIn->surfaceSymmetry);
        if (controlOut->surfaceSymmetry == NULL) return EGADS_MALLOC;
    }

    return CAPS_SUCCESS;
}

// Make a copy of vlmSectionStruct (it is assumed sectionOut has already been initialized)
int copy_vlmSectionStruct(vlmSectionStruct *sectionIn, vlmSectionStruct *sectionOut) {

    int status; // Function return status

    int i; // Indexing

    int numObj;
    ego *objs;

    if (sectionIn == NULL) return CAPS_NULLVALUE;
    if (sectionOut == NULL) return CAPS_NULLVALUE;

    status = destroy_vlmSectionStruct(sectionOut);
    if (status != CAPS_SUCCESS) return status;

    if (sectionIn->name != NULL) {
        sectionOut->name = EG_strdup(sectionIn->name);
        if (sectionOut->name == NULL) return EGADS_MALLOC;
    }

    status = EG_copyObject(sectionIn->ebody, NULL, &sectionOut->ebody);
    if (status != EGADS_SUCCESS) return status;
    sectionOut->type  = sectionIn->type;

    sectionOut->sectionIndex = sectionIn->sectionIndex; // Section index - 0 bias

    sectionOut->xyzLE[0] = sectionIn->xyzLE[0]; // xyz coordinates for leading edge
    sectionOut->xyzLE[1] = sectionIn->xyzLE[1];
    sectionOut->xyzLE[2] = sectionIn->xyzLE[2];

    sectionOut->nodeIndexLE = sectionIn->nodeIndexLE; // Leading edge node index with reference to xyzLE

    sectionOut->xyzTE[0] = sectionIn->xyzTE[0]; // xyz location of the trailing edge
    sectionOut->xyzTE[1] = sectionIn->xyzTE[1];
    sectionOut->xyzTE[2] = sectionIn->xyzTE[2];

    status = EG_getBodyTopos(sectionOut->ebody, NULL, sectionIn->teClass, &numObj, &objs);
    if (status != EGADS_SUCCESS) return status;
    for (i = 0; i < numObj; i++) {
      if ( EG_isEquivalent(objs[i], sectionIn->teObj) == EGADS_SUCCESS) {
        sectionOut->teObj = objs[i]; // Trailing edge object
        break;
      }
    }
    AIM_FREE(objs);

    sectionOut->teClass = sectionIn->teClass; // Trailing edge object class (NODE or EDGE)

    sectionOut->chord = sectionIn->chord; // Chord
    sectionOut->ainc  = sectionIn->ainc;  // Incidence angle

    sectionOut->normal[0] = sectionIn->normal[0]; // planar normal
    sectionOut->normal[1] = sectionIn->normal[1];
    sectionOut->normal[2] = sectionIn->normal[2];

    sectionOut->Nspan = sectionIn->Nspan;   // number of spanwise vortices (elements)
    sectionOut->Sspace = sectionIn->Sspace; // spanwise point distribution
    sectionOut->Sset = sectionIn->Sset;

    sectionOut->numControl = sectionIn->numControl;

    if (sectionOut->numControl != 0) {

        sectionOut->vlmControl = (vlmControlStruct *) EG_alloc(sectionOut->numControl*sizeof(vlmControlStruct));
        if (sectionOut->vlmControl == NULL) return EGADS_MALLOC;

        for (i = 0; i < sectionOut->numControl; i++) {

            status = initiate_vlmControlStruct(&sectionOut->vlmControl[i]);
            if (status != CAPS_SUCCESS) return status;

            status = copy_vlmControlStruct(&sectionIn->vlmControl[i], &sectionOut->vlmControl[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    return CAPS_SUCCESS;
}

// Make a copy of vlmSurfaceStruct (it is assumed surfaceOut has already been initialized)
// Also the section in vlmSurface are reordered based on a vlm_orderSections() function call
int copy_vlmSurfaceStruct(void *aimInfo, vlmSurfaceStruct *surfaceIn, vlmSurfaceStruct *surfaceOut) {

    int status; // Function return status

    int i, sectionIndex; // Indexing

    if (surfaceIn == NULL) return CAPS_NULLVALUE;
    if (surfaceOut == NULL) return CAPS_NULLVALUE;

    status = destroy_vlmSurfaceStruct(surfaceOut);
    if (status != CAPS_SUCCESS) return status;

    if (surfaceIn->name != NULL) {
        surfaceOut->name = EG_strdup(surfaceIn->name);
        if (surfaceOut->name == NULL) return EGADS_MALLOC;
    }

    surfaceOut->numAttr = surfaceIn->numAttr;

    if (surfaceIn->attrIndex != NULL) {

        surfaceOut->attrIndex = (int *) EG_alloc(surfaceOut->numAttr *sizeof(int));
        if (surfaceOut->attrIndex == NULL) return EGADS_MALLOC;

        memcpy(surfaceOut->attrIndex,
                surfaceIn->attrIndex,
                surfaceOut->numAttr*sizeof(int));
    }

    surfaceOut->Nchord = surfaceIn->Nchord;
    surfaceOut->Cspace = surfaceIn->Cspace;

    surfaceOut->NspanTotal   = surfaceIn->NspanTotal;
    surfaceOut->NspanSection = surfaceIn->NspanSection;
    surfaceOut->Sspace       = surfaceIn->Sspace;

    surfaceOut->nowake = surfaceIn->nowake;
    surfaceOut->noalbe = surfaceIn->noalbe;
    surfaceOut->noload = surfaceIn->noload;

    surfaceOut->compon = surfaceIn->compon;
    surfaceOut->iYdup = surfaceIn->iYdup;

    surfaceOut->numSection = surfaceIn->numSection;


    if (surfaceIn->vlmSection != NULL) {

        AIM_ALLOC(surfaceOut->vlmSection, surfaceOut->numSection, vlmSectionStruct, aimInfo, status);

        status = vlm_orderSections(aimInfo, surfaceIn->numSection, surfaceIn->vlmSection);
        AIM_STATUS(aimInfo, status);

        for (i = 0; i < surfaceOut->numSection; i++) {

            status = initiate_vlmSectionStruct(&surfaceOut->vlmSection[i]);
            AIM_STATUS(aimInfo, status);

            // Sections aren't necessarily stored in order coming out of vlm_GetSection, however sectionIndex is (after a
            // call to vlm_orderSection()) !

            sectionIndex = surfaceIn->vlmSection[i].sectionIndex;

            status = copy_vlmSectionStruct(&surfaceIn->vlmSection[sectionIndex], &surfaceOut->vlmSection[i]);
            AIM_STATUS(aimInfo, status);

            // Reset the sectionIndex that is keeping track of the section order.
            surfaceOut->vlmSection[i].sectionIndex = i;
        }
    }

    if (surfaceIn->surfaceType != NULL) {
        surfaceOut->surfaceType = EG_strdup(surfaceIn->surfaceType);
    }

    status = CAPS_SUCCESS;
cleanup:
    return status;
}


static
int vlm_findLeadingEdge(void *aimInfo, int numNode, ego *nodes, int *nodeIndexLE, double *xyzLE )
{
    int status; // Function return
    int  i; // Indexing

    // Node variables
    double xmin=0, xyz[3]={0.0, 0.0, 0.0};

    *nodeIndexLE = 0;

    // Assume the LE position is the most forward Node in X
    for (i = 0; i < numNode; i++) {
        status = EG_evaluate(nodes[i], NULL, xyz);
        AIM_STATUS(aimInfo, status);

        if (*nodeIndexLE == 0) {
            *nodeIndexLE = i+1;
            xmin = xyz[0];
        } else {

            if (xyz[0] < xmin) {
                *nodeIndexLE = i+1;
                xmin = xyz[0];
            }
        }
    }

    if (*nodeIndexLE == 0) {
        AIM_ERROR(aimInfo, "Body has no Leading Edge Node!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    status = EG_evaluate(nodes[*nodeIndexLE-1], NULL, xyzLE);
    AIM_STATUS(aimInfo, status);

    status = CAPS_SUCCESS;
cleanup:

    return status;
}

// Find the EGO object pertaining the to trailing edge
static
int vlm_findTrailingEdge(void *aimInfo,
                         int numNode, ego *nodes,
                         int numEdge, ego *edges,
            /*@unused@*/ int nodeIndexLE,
                         double *secnorm,
                         ego *teObj,
                         int *teClass,
                         double *xyzTE) {

    int status = CAPS_SUCCESS; // Function return status

    int i; // Indexing
    int nodeIndexTE;

    double t, norm, xmax;
    double trange[2], normEdge[3], result[18], vec1[3], vec2[3], xyz[3];

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    ego ref, *children = NULL;

    *teObj = NULL;

    // Find the node with the most rear X
    nodeIndexTE = 0;
    for (i = 0; i < numNode; i++) {
        status = EG_evaluate(nodes[i], NULL, xyz);
        AIM_STATUS(aimInfo, status);

        if (nodeIndexTE == 0) {
            nodeIndexTE = i+1;
            xmax = xyz[0];
        } else {
            if (xyz[0] > xmax) {
                nodeIndexTE = i+1;
                xmax = xyz[0];
            }
        }
    }

    if (nodeIndexTE == 0) {
        printf("Error in vlm_findTrailingEdge: Body has no TE node!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    for (i = 0; i < numEdge; i++) {

        status = EG_getTopology(edges[i], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);
        if (mtype == DEGENERATE) continue;

        if (status != EGADS_SUCCESS) {
            printf("Error in vlm_findTrailingEdge: Edge %d getTopology = %d!\n", i, status);
            goto cleanup;
        }

        if (numChildren != 2) {
            printf("Error in vlm_findTrailingEdge: Edge %d has %d nodes!\n", i, numChildren);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // If the edge doesn't at least contain the TE node pass it by
        if (children[0] != nodes[nodeIndexTE-1] && children[1] != nodes[nodeIndexTE-1]) continue;

        // evaluate at the edge mid point
        t = 0.5*(trange[0]+trange[1]);
        status = EG_evaluate(edges[i], &t, result);
        AIM_STATUS(aimInfo, status);

        // get the tangent vector
        vec1[0] = result[3+0];
        vec1[1] = result[3+1];
        vec1[2] = result[3+2];

        // cross it to get the 'normal' to the edge (i.e. in the airfoil section PLANE)
        cross_DoubleVal(vec1,secnorm,normEdge);

        // get the tangent vectors at the end points and make sure the dot product is near 1

        // get the tangent vector at t0
        t = trange[0];
        status = EG_evaluate(edges[i], &t, result);
        AIM_STATUS(aimInfo, status);

        vec1[0] = result[3+0];
        vec1[1] = result[3+1];
        vec1[2] = result[3+2];

        norm = sqrt(dot_DoubleVal(vec1,vec1));

        vec1[0] /= norm;
        vec1[1] /= norm;
        vec1[2] /= norm;

        // get the tangent vector at t1
        t = trange[1];
        status = EG_evaluate(edges[i], &t, result);
        AIM_STATUS(aimInfo, status);

        vec2[0] = result[3+0];
        vec2[1] = result[3+1];
        vec2[2] = result[3+2];

        norm = sqrt(dot_DoubleVal(vec2,vec2));

        vec2[0] /= norm;
        vec2[1] /= norm;
        vec2[2] /= norm;

        // compute the dot between the two tangents
        norm = fabs( dot_DoubleVal(vec1,vec2) );

        // if the x-component of the normal is larger, assume the edge is pointing in the streamwise direction
        // the tangent at the end points must also be pointing in the same direction
        if (fabs(normEdge[0]) > sqrt(normEdge[1]*normEdge[1] + normEdge[2]*normEdge[2]) && (1 - norm) < 1e-3) {

          if (*teObj != NULL) {
              printf("\tError in vlm_findTrailingEdge: Found multiple trailing edges!!\n");
              status = CAPS_SOURCEERR;
              goto cleanup;
          }

          *teObj = edges[i];
        }
    }

    // Assume a sharp trailing edge and use the Node
    if (*teObj == NULL) {
        *teObj = nodes[nodeIndexTE-1];

        // Get the class and coordinates
        status = EG_getTopology(*teObj, &ref, teClass, &mtype, xyzTE, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);

    } else {

        // Get the class and t-range for mid point evaluation
        status = EG_getTopology(*teObj, &ref, teClass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);

        t = 0.5*(trange[0]+trange[1]);

        status = EG_evaluate(*teObj, &t, result);
        AIM_STATUS(aimInfo, status);

        xyzTE[0] = result[0];
        xyzTE[1] = result[1];
        xyzTE[2] = result[2];
    }


cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Premature exit in vlm_findTrailingEdge, status = %d\n", status);

    return status;
}


// Get the normal to the airfoil cross-section plane
static
int vlm_secNormal(void *aimInfo, ego body,
                  double *secnorm, vlmSectionEnum *secType)
{
    int status; // Function return status

    int i, j; // Indexing

    int numEdge, numLoop, numFace;
    int *ivec=NULL;
    double *rvec=NULL;

    double norm, t, alen, dir[2][9];
    double trange[4], nodesXYZ[2][9], dX1[3], dX2[3];

    //EGADS returns
    int oclass, mtype, *sens = NULL, *esens=NULL, numChildren;
    ego ref, *children=NULL, *edges=NULL, *loops=NULL, *faces=NULL, refGeom;

    secnorm[0] = secnorm[1] = secnorm[2] = 0.;

    if (aim_isNodeBody(body, trange) == CAPS_SUCCESS) {
      *secType = vlmSecNode;
      return CAPS_SUCCESS;
    }

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, LOOP, &numLoop, &loops);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &faces);
    AIM_STATUS(aimInfo, status);

    *secType = vlmSecLine;
    for (i = 0; i < numEdge; i++) {
      status = EG_getTopology(edges[i], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
      AIM_STATUS(aimInfo, status);
      if (mtype == DEGENERATE) continue;

      status = EG_getGeometry(ref, &oclass, &mtype, &refGeom, NULL, NULL);
      AIM_STATUS(aimInfo, status);

      if (mtype == LINE) {
        if (i == 0) {
          status = EG_evaluate(edges[i], &trange[0], dir[0]);
          AIM_STATUS(aimInfo, status);

          norm = sqrt(dot_DoubleVal(dir[0]+3, dir[0]+3));
          dir[0][3] /= norm;
          dir[0][4] /= norm;
          dir[0][5] /= norm;

          continue;
        }
      }
      if (mtype == BSPLINE) {
        status = EG_arcLength(ref, trange[0], trange[1], &alen);
        AIM_STATUS(aimInfo, status);

        status = EG_evaluate(children[0], NULL, nodesXYZ[0]);
        AIM_STATUS(aimInfo, status);

        status = EG_evaluate(children[1], NULL, nodesXYZ[1]);
        AIM_STATUS(aimInfo, status);

        dX1[0] = nodesXYZ[1][0] - nodesXYZ[0][0];
        dX1[1] = nodesXYZ[1][1] - nodesXYZ[0][1];
        dX1[2] = nodesXYZ[1][2] - nodesXYZ[0][2];

        norm = sqrt(dot_DoubleVal(dX1, dX1));

        // check if the integrated lenght is the same as the linear distance between nodes
        if (fabs(norm - alen)/alen < 1e-6) {
          if (i == 0) {
            status = EG_evaluate(edges[i], &trange[0], dir[0]);
            AIM_STATUS(aimInfo, status);

            norm = sqrt(dot_DoubleVal(dir[0]+3, dir[0]+3));
            dir[0][3] /= norm;
            dir[0][4] /= norm;
            dir[0][5] /= norm;
            continue;
          }
        }
      }

      // Check that all Edges are parallel
      if (i > 0) {
        status = EG_evaluate(edges[i], &trange[0], dir[1]);
        AIM_STATUS(aimInfo, status);

        norm = sqrt(dot_DoubleVal(dir[1]+3, dir[1]+3));
        dir[1][2] /= norm;
        dir[1][3] /= norm;
        dir[1][4] /= norm;

        norm = fabs(dot_DoubleVal(dir[0]+3, dir[1]+3));
        if (fabs(norm - 1) < 1e-6) {
          continue;
        }
      }

      // Not a single line geometry
      *secType = vlmSecAirfoil;
      break;
    }

    if (*secType == vlmSecLine) {
      goto cleanup;
    }

    status = EG_isPlanar(body);
    if (status != EGADS_SUCCESS) {
        AIM_ERROR(aimInfo, "body is not planar!");
        status = CAPS_SOURCEERR;
        goto cleanup;
    }

    // get the PLANE normal vector for the airfoil section
    if (numFace == 1) {
        status = EG_getTopology(faces[0], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);
        i = mtype;

        status = EG_getGeometry(ref, &oclass, &mtype, &refGeom, &ivec, &rvec);
        AIM_STATUS(aimInfo, status);

        if (oclass == SURFACE && mtype == PLANE) {
            cross_DoubleVal(rvec+3, rvec+6, secnorm);
            secnorm[0] *= i;
            secnorm[1] *= i;
            secnorm[2] *= i;
            goto cleanup;
        }
        AIM_FREE(ivec);
        AIM_FREE(rvec);
    }
    AIM_FREE(faces);

    // get the edge senses from the loop
    status = EG_getTopology(loops[0], &ref, &oclass, &mtype, trange, &numChildren, &children, &esens);
    AIM_STATUS(aimInfo, status);

    for (i = 0; i < numEdge && (secnorm[0] == 0 && secnorm[1] == 0 && secnorm[2] == 0); i++) {

        status = EG_getTopology(edges[i], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);
        if (mtype == DEGENERATE) continue;

        if (esens[i] == SFORWARD) {
            status = EG_evaluate(children[0], NULL, nodesXYZ[0]);
        } else {
            status = EG_evaluate(children[1], NULL, nodesXYZ[0]);
        }
        AIM_STATUS(aimInfo, status);

        t = (trange[0] + trange[1])/2.;
        status = EG_evaluate(edges[i], &t, nodesXYZ[1]);
        AIM_STATUS(aimInfo, status);

        dX1[0] = nodesXYZ[1][0] - nodesXYZ[0][0];
        dX1[1] = nodesXYZ[1][1] - nodesXYZ[0][1];
        dX1[2] = nodesXYZ[1][2] - nodesXYZ[0][2];

        for (j = 0; j < numEdge; j++) {
            if (i == j) continue;

            status = EG_getTopology(edges[j], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
            AIM_STATUS(aimInfo, status);
            if (mtype == DEGENERATE) continue;

            t = (trange[0] + trange[1])/2.;
            status = EG_evaluate(edges[j], &t, nodesXYZ[1]);
            AIM_STATUS(aimInfo, status);

            dX2[0] = nodesXYZ[1][0] - nodesXYZ[0][0];
            dX2[1] = nodesXYZ[1][1] - nodesXYZ[0][1];
            dX2[2] = nodesXYZ[1][2] - nodesXYZ[0][2];

            if (fabs(dot_DoubleVal(dX1,dX2)) < 1e-7) continue;

            cross_DoubleVal(dX1, dX2, secnorm);
            break;
        }
    }
    if (secnorm[0] == 0 && secnorm[1] == 0 && secnorm[2] == 0) {
        AIM_ERROR(aimInfo, "Failed to determine airfoil section PLANE normal!");
        status = EGADS_GEOMERR;
        goto cleanup;
    }

    // normalize the section normal vector
    norm = sqrt(dot_DoubleVal(secnorm,secnorm));

    secnorm[0] /= norm;
    secnorm[1] /= norm;
    secnorm[2] /= norm;

cleanup:

    AIM_FREE(ivec);
    AIM_FREE(rvec);
    AIM_FREE(edges);
    AIM_FREE(loops);
    AIM_FREE(faces);
    return status;
}

// Finalizes populating vlmSectionStruct member data after the ebody is set
int finalize_vlmSectionStruct(void *aimInfo, vlmSectionStruct *vlmSection)
{
    int status; // Function return status

    int i; // Index
    int numNode, numEdge, numLoop;
    int numEdgeMinusDegenrate;
    ego ebody;
    ego *nodes = NULL, *edges = NULL;

    double xdot[3], X[3] = {1,0,0}, Y[3];

    //EGADS returns
    int oclass, mtype;
    ego ref, prev, next;

    ebody = vlmSection->ebody;

    if (aim_isNodeBody(ebody, X) == CAPS_SUCCESS) {
      vlmSection->type = vlmSecNode;
      vlmSection->nodeIndexLE = 0;
      vlmSection->teObj = NULL;
      vlmSection->teClass = 0;
      vlmSection->xyzLE[0] = vlmSection->xyzTE[0] = X[0];
      vlmSection->xyzLE[1] = vlmSection->xyzTE[1] = X[1];
      vlmSection->xyzLE[2] = vlmSection->xyzTE[2] = X[2];
      vlmSection->chord = 0;
      vlmSection->ainc = 0;

      return CAPS_SUCCESS;
    }

    status = EG_getBodyTopos(ebody, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(ebody, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(ebody, NULL, LOOP, &numLoop, NULL);
    AIM_STATUS(aimInfo, status);

    numEdgeMinusDegenrate = 0;
    for (i = 0; i < numEdge; i++) {
        status = EG_getInfo(edges[i], &oclass, &mtype, &ref, &prev, &next);
        AIM_STATUS(aimInfo, status);
        if (mtype == DEGENERATE) continue;
        numEdgeMinusDegenrate += 1;
    }

    if (vlmSection->type == vlmSecAirfoil) {
      // There must be at least 2 nodes and 2 edges
      if ((numEdgeMinusDegenrate != numNode) || (numNode < 2) || (numLoop != 1)) {
        AIM_ERROR  (aimInfo, "Body has %d Nodes, %d Edges and %d Loops!", numNode, numEdge, numLoop);
        AIM_ADDLINE(aimInfo, "The body must have at least one leading and one trailing edge Node and only one Loop!");
        status = CAPS_SOURCEERR;
        goto cleanup;
      }
    }

    // Get the section normal from the body
    status = vlm_secNormal(aimInfo,
                           ebody,
                           vlmSection->normal, &vlmSection->type);
    AIM_STATUS(aimInfo, status);

    // Find the leadinge edge Node
    status = vlm_findLeadingEdge(aimInfo, numNode, nodes, &vlmSection->nodeIndexLE, vlmSection->xyzLE);
    AIM_STATUS(aimInfo, status);

    // Find the trailing edge Object (Node or EDGE)
    status = vlm_findTrailingEdge(aimInfo,
                                  numNode, nodes,
                                  numEdge, edges,
                                  vlmSection->nodeIndexLE,
                                  vlmSection->normal,
                                  &vlmSection->teObj,
                                  &vlmSection->teClass,
                                  vlmSection->xyzTE);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(vlmSection->teObj, aimInfo ,status);

    // There must be at least 2 nodes and 2 nonTE edges
    if ((numEdgeMinusDegenrate == 2) && (vlmSection->teObj->oclass == EDGE)) {
        AIM_ERROR  (aimInfo, "Body has %d Nodes and %d Edges with one trailing edge Edge!", numNode, numEdge);
        AIM_ADDLINE(aimInfo, "The body must have at least one leading edge Node!");
        status = CAPS_SOURCEERR;
        goto cleanup;
    }

    xdot[0] = vlmSection->xyzTE[0] - vlmSection->xyzLE[0];
    xdot[1] = vlmSection->xyzTE[1] - vlmSection->xyzLE[1];
    xdot[2] = vlmSection->xyzTE[2] - vlmSection->xyzLE[2];

    vlmSection->chord = sqrt(dot_DoubleVal(xdot,xdot));

    xdot[0] /= vlmSection->chord;
    xdot[1] /= vlmSection->chord;
    xdot[2] /= vlmSection->chord;

    // cross with section PLANE normal to get perpendicular vector in the PLANE
    cross_DoubleVal(vlmSection->normal, X, Y);

    vlmSection->ainc = -atan2(dot_DoubleVal(xdot,Y), xdot[0])*180./PI;

cleanup:

    AIM_FREE(nodes);
    AIM_FREE(edges);

    return status;
}

// flips a section
static
int vlm_flipSection(void *aimInfo, ego body, ego *flipped)
{
    int status; // Function return status

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    double data[4];
    ego context, ref, *children = NULL, eflip;

    status = EG_getContext(body, &context);
    AIM_STATUS(aimInfo, status);

    // get the child of the body as the body it self cannot be flipped
    status = EG_getTopology(body, &ref, &oclass, &mtype, data, &numChildren, &children, &sens);
    AIM_STATUS(aimInfo, status);

    if (numChildren != 1) {
      AIM_ERROR(aimInfo, "Body has %d children (may only have 1)!", numChildren);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Flip the airfoil so the normals are consistent
    status = EG_flipObject(children[0], &eflip);
    AIM_STATUS(aimInfo, status);

    // create the new body with the flipped airfoil
    status = EG_makeTopology(context, NULL, oclass, mtype, data, 1, &eflip, sens, flipped);
    AIM_STATUS(aimInfo, status);

    // copy over the body attributes
    status = EG_attributeDup(body, *flipped);
    AIM_STATUS(aimInfo, status);

cleanup:
    return status;
}


// returns a list of bodies with normal vectors pointing in the negative y- or z-directions
// requires that all sections be in y- or z-constant planes
static
int vlm_getSectionYZ(void *aimInfo, ego body, ego *copy)
{
    int status;
    vlmSectionEnum secType;

    double secnorm[3] = {0,0,0};

    status = vlm_secNormal(aimInfo, body, secnorm, &secType);
    AIM_STATUS(aimInfo, status);

    if ( fabs(fabs(secnorm[1]) - 1.) > DOTTOL &&
         fabs(fabs(secnorm[2]) - 1.) > DOTTOL ) {
        AIM_ERROR(aimInfo, "Section is neither purely in the y- or the z-plane.");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if ( fabs(fabs(secnorm[1]) - 1.) < DOTTOL ) {
        if (secnorm[1] > 0.) {
            // Flip the body
            status = vlm_flipSection(aimInfo, body, copy);
            AIM_STATUS(aimInfo, status);
        } else {
            // Store a copy of the body
            status = EG_copyObject(body, NULL, copy);
            AIM_STATUS(aimInfo, status);
        }
    } else {
        if (secnorm[2] > 0.) {
            // Flip the body
            status = vlm_flipSection(aimInfo, body, copy);
            AIM_STATUS(aimInfo, status);
        } else {
            // Store a copy of the body
            status = EG_copyObject(body, NULL, copy);
            AIM_STATUS(aimInfo, status);
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    return status;
}


// returns a list of bodies with normal vectors pointing in the negative radial direction
static
int vlm_getSectionRadial(void *aimInfo, ego body, ego *copy)
{
    int status;

    int numNode;
    ego *nodes = NULL;

    int nodeIndexLE;
    double xyzLE[3] = {0,0,0}, radLE[3], norm;

    double secnorm[3] = {0,0,0};
    vlmSectionEnum secType;

    if (aim_isNodeBody(body, xyzLE) == CAPS_SUCCESS) {
      // Store a copy of the body
      status = EG_copyObject(body, NULL, copy);
      AIM_STATUS(aimInfo, status);
      return CAPS_SUCCESS;
    }

    status = vlm_secNormal(aimInfo, body, secnorm, &secType);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = vlm_findLeadingEdge(aimInfo, numNode, nodes, &nodeIndexLE, xyzLE );
    AIM_STATUS(aimInfo, status);


    radLE[0] = 0;
    radLE[1] = xyzLE[1];
    radLE[2] = xyzLE[2];

    norm = sqrt(dot_DoubleVal(radLE,radLE));

    if (norm < DOTTOL) {
        AIM_ERROR(aimInfo, "Section LE cannot be on y = 0 and z = 0!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    radLE[1] /= norm;
    radLE[2] /= norm;

    if ( fabs(secnorm[1]*radLE[1] + secnorm[2]*radLE[2]) < DOTTOL) {
        AIM_ERROR(aimInfo, "Section normal is not radial!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // make sure the radial dot product is negative
    if ( secnorm[1]*radLE[1] + secnorm[2]*radLE[2] > 0. ) {
        // Flip the body
        status = vlm_flipSection(aimInfo, body, copy);
        AIM_STATUS(aimInfo, status);
    } else {
        // Store a copy of the body
        status = EG_copyObject(body, NULL, copy);
        AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(nodes);

    return status;
}


// Accumulate VLM section data from a set of bodies. If disciplineFilter is not NULL
// bodies not found with disciplineFilter (case insensitive) for a capsDiscipline attribute
// will be ignored.
int vlm_getSections(void *aimInfo,
                    int numBody,
                    ego bodies[],
                    /*@null@*/ const char *disciplineFilter,
                    mapAttrToIndexStruct attrMap,
                    vlmSystemEnum sys,
                    int numSurface,
                    vlmSurfaceStruct *vlmSurface[]) {

    int status; // Function return status

    int i, k, body, surf, section; // Indexing
    int attrIndex;
    int Nspan;
    double Sspace;
    double xyz[3], mag;
    vlmSectionEnum secType;

    int found = (int) false; // Boolean tester

    const char *groupName = NULL;
    const char *discipline = NULL;

    // Loop through bodies
    for (body = 0; body < numBody; body++) {

        if (disciplineFilter != NULL) {

            status = retrieve_CAPSDisciplineAttr(bodies[body], &discipline);
            if (status == CAPS_SUCCESS) {
                //printf("capsDiscipline = %s, Body, %i\n", discipline, body);
                if (strcasecmp(discipline, disciplineFilter) != 0) continue;
            }
        }

        // Loop through surfaces
        for (surf = 0; surf < numSurface; surf++) {

            status = retrieve_CAPSGroupAttr(bodies[body], &groupName);
            if (status != CAPS_SUCCESS) {
                printf("Warning (vlm_getSections): No capsGroup value found on body %d, body will not be used\n", body+1);
                continue;
            }

            status = get_mapAttrToIndexIndex(&attrMap, groupName, &attrIndex);
            if (status == CAPS_NOTFOUND) {
                AIM_ERROR(aimInfo, "VLM Surface name \"%s\" not found in attrMap\n", groupName);
                goto cleanup;
            }

            found = (int) false;

            // See if attrIndex is in the attrIndex array for the surface
            for (i = 0; i < (*vlmSurface)[surf].numAttr; i++) {

                if (attrIndex == (*vlmSurface)[surf].attrIndex[i]) {
                    found = (int) true;
                    break;
                }
            }

            // If attrIndex isn't in the array the; body doesn't belong in this surface
            if (found == (int) false) continue;

            // realloc vlmSection size
            AIM_REALL((*vlmSurface)[surf].vlmSection, (*vlmSurface)[surf].numSection+1, vlmSectionStruct, aimInfo, status);
            (*vlmSurface)[surf].numSection += 1;

            // Get section index
            section = (*vlmSurface)[surf].numSection-1;

            status = initiate_vlmSectionStruct(&(*vlmSurface)[surf].vlmSection[section]);
            AIM_STATUS(aimInfo, status);

            // Store the section index
            (*vlmSurface)[surf].vlmSection[section].sectionIndex = section;

            // get the specified number of span points from the body
            Nspan = 0;
            status = retrieve_intAttrOptional(bodies[body], "vlmNumSpan", &Nspan);
            if (status == EGADS_ATTRERR) goto cleanup;
            (*vlmSurface)[surf].vlmSection[section].Nspan = Nspan;

            // get the specified span points distribution from the body
            Sspace = (*vlmSurface)[surf].Sspace;
            status = retrieve_doubleAttrOptional(bodies[body], "vlmSspace", &Sspace);
            if (status == EGADS_ATTRERR) goto cleanup;
            (*vlmSurface)[surf].vlmSection[section].Sspace = Sspace;
            (*vlmSurface)[surf].vlmSection[section].Sset = (int)(status == CAPS_SUCCESS);

            // Get the section normal
            status = vlm_secNormal(aimInfo, bodies[body],
                                   (*vlmSurface)[surf].vlmSection[section].normal, &secType);
            AIM_STATUS(aimInfo, status);

            if (secType == vlmSecLine) {

              // Store a copy of the body
              status = EG_copyObject(bodies[body], NULL, &(*vlmSurface)[surf].vlmSection[section].ebody);
              AIM_STATUS(aimInfo, status);

              if (section >= 1) {
                xyz[0] = 0;
                xyz[1] = (*vlmSurface)[surf].vlmSection[section].normal[1] - (*vlmSurface)[surf].vlmSection[section-1].normal[1];
                xyz[2] = (*vlmSurface)[surf].vlmSection[section].normal[2] - (*vlmSurface)[surf].vlmSection[section-1].normal[2];
                
                mag = sqrt(dot_DoubleVal(xyz,xyz));
                
                (*vlmSurface)[surf].vlmSection[section].normal[0] = 0;
                (*vlmSurface)[surf].vlmSection[section].normal[1] = xyz[1]/mag;
                (*vlmSurface)[surf].vlmSection[section].normal[2] = xyz[2]/mag;
                
                if (section == 1 && (*vlmSurface)[surf].vlmSection[section-1].type == vlmSecLine) {
                  (*vlmSurface)[surf].vlmSection[0].normal[0] = (*vlmSurface)[surf].vlmSection[1].normal[0];
                  (*vlmSurface)[surf].vlmSection[0].normal[1] = (*vlmSurface)[surf].vlmSection[1].normal[1];
                  (*vlmSurface)[surf].vlmSection[0].normal[2] = (*vlmSurface)[surf].vlmSection[1].normal[2];
                }
              }

            } else {

              // modify bodies as needed for the given coordinate system
              if (sys == vlmGENERIC) {

                  // For a generic system the section normal vectors must be consistent
                  k = 0;
                  while (fabs((*vlmSurface)[surf].vlmSection[k].normal[1]*(*vlmSurface)[surf].vlmSection[section].normal[1] +
                              (*vlmSurface)[surf].vlmSection[k].normal[2]*(*vlmSurface)[surf].vlmSection[section].normal[2]) < DOTTOL) {
                      k++;
                      if (k == section) {
                          AIM_ERROR(aimInfo, "Body %d is orthogonal to all other airfoils!\n", body+1);
                          status = CAPS_NOTFOUND;
                          goto cleanup;
                      }
                  }

                  if (section == 0) {

                    if ((*vlmSurface)[surf].vlmSection[section].normal[1] > 0) {
                        // Flip the body
                        status = vlm_flipSection(aimInfo, bodies[body], &(*vlmSurface)[surf].vlmSection[section].ebody);
                        AIM_STATUS(aimInfo, status);
                    } else {
                        // Store a copy of the body
                        status = EG_copyObject(bodies[body], NULL, &(*vlmSurface)[surf].vlmSection[section].ebody);
                        AIM_STATUS(aimInfo, status);
                    }

                  } else { // section == 0

                    if (((*vlmSurface)[surf].vlmSection[k].normal[1]*(*vlmSurface)[surf].vlmSection[section].normal[1] +
                         (*vlmSurface)[surf].vlmSection[k].normal[2]*(*vlmSurface)[surf].vlmSection[section].normal[2]) < 0) {
                        // Flip the body
                        status = vlm_flipSection(aimInfo, bodies[body], &(*vlmSurface)[surf].vlmSection[section].ebody);
                        AIM_STATUS(aimInfo, status);
                    } else {
                        // Store a copy of the body
                        status = EG_copyObject(bodies[body], NULL, &(*vlmSurface)[surf].vlmSection[section].ebody);
                        AIM_STATUS(aimInfo, status);
                    }

                  }

              } else if (sys == vlmPLANEYZ) {
                  status = vlm_getSectionYZ(aimInfo, bodies[body], &(*vlmSurface)[surf].vlmSection[section].ebody);
                  AIM_STATUS(aimInfo, status);
              } else if (sys == vlmRADIAL) {
                  status = vlm_getSectionRadial(aimInfo, bodies[body], &(*vlmSurface)[surf].vlmSection[section].ebody);
                  AIM_STATUS(aimInfo, status);

                  if (aim_isNodeBody(bodies[body], xyz) == CAPS_SUCCESS) {
                    xyz[0] = 0;
                    mag = sqrt(dot_DoubleVal(xyz,xyz));

                    (*vlmSurface)[surf].vlmSection[section].normal[0] = -xyz[0]/mag;
                    (*vlmSurface)[surf].vlmSection[section].normal[1] = -xyz[1]/mag;
                    (*vlmSurface)[surf].vlmSection[section].normal[2] = -xyz[2]/mag;
                  }
              } else {
                  AIM_ERROR(aimInfo, "Developer Error: Unknown coordinate system");
                  status = CAPS_BADVALUE;
                  goto cleanup;
              }
            }

            // Populate remaining data after the body is set
            status = finalize_vlmSectionStruct(aimInfo, &(*vlmSurface)[surf].vlmSection[section]);
            AIM_STATUS(aimInfo, status);
        }
    }

    // order the sections in the surfaces
    for (surf = 0; surf < numSurface; surf++) {
        status = vlm_orderSections(aimInfo, (*vlmSurface)[surf].numSection, (*vlmSurface)[surf].vlmSection);
        AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:
    return status;
}


// Order VLM sections increasing order
int vlm_orderSections(void *aimInfo, int numSection, vlmSectionStruct vlmSections[])
{
    int i1, i2, j, k, hit;
    double dot, vec[3];
    vlmSectionEnum secType;

    // one section is always in order
    if (numSection == 1) return CAPS_SUCCESS;

    if (numSection <= 0) {
        AIM_ERROR(aimInfo, "invalid number of sections %d!\n", numSection);
        return CAPS_BADVALUE;
    }

    secType = vlmSecLine;
    for (k = 0; k < numSection; k++) {
      if (vlmSections[k].type != vlmSecLine) {
        secType = vlmSections[k].type;
        break;
      }
    }

    // Section lines muts be sorted by the user
    if (secType == vlmSecLine) {
      return CAPS_SUCCESS;
    }

    // the loop below will get stuck in an infinite loop if the normals are not consistent
    for (k = 1; k < numSection; k++) {
        j = 0;
        while(fabs(vlmSections[j].normal[1]*vlmSections[k].normal[1] +
                   vlmSections[j].normal[2]*vlmSections[k].normal[2]) < DOTTOL) {
          j++;
          if(j == numSection) {
            AIM_ERROR(aimInfo, "One airfoil is orthogonal to all other airfoils!");
            return CAPS_NOTFOUND;
          }
        }

        if ((vlmSections[j].normal[1]*vlmSections[k].normal[1] +
             vlmSections[j].normal[2]*vlmSections[k].normal[2]) < 0) {
            AIM_ERROR(aimInfo, "Section normals are not consistent: (%lf,%lf).(%lf,%lf) < 0!",
                      vlmSections[j].normal[1], vlmSections[j].normal[2],
                      vlmSections[k].normal[1], vlmSections[k].normal[2]);
            return CAPS_BADVALUE;
        }
    }

    // order the sections so the dot product between the the section normals and the
    // distance vector between sections is negative
    do {
        hit = 0;
        for (k = 0; k < numSection-1; k++) {
            i1 = vlmSections[k  ].sectionIndex;
            i2 = vlmSections[k+1].sectionIndex;

            vec[0] = 0;
            vec[1] = vlmSections[i2].xyzLE[1] - vlmSections[i1].xyzLE[1];
            vec[2] = vlmSections[i2].xyzLE[2] - vlmSections[i1].xyzLE[2];

            dot = vec[1]*(vlmSections[i1].normal[1] + vlmSections[i2].normal[1])/2. +
                  vec[2]*(vlmSections[i1].normal[2] + vlmSections[i2].normal[2])/2.;

            if (dot > 0) {
              vlmSections[k  ].sectionIndex = i2;
              vlmSections[k+1].sectionIndex = i1;
              hit++;
            }
        }
    } while (hit != 0);

    return CAPS_SUCCESS;
}


// Compute spanwise panel spacing with close to equal spacing on each panel
int vlm_equalSpaceSpanPanels(void *aimInfo, int NspanTotal, int numSection, vlmSectionStruct vlmSection[])
{
    int status = CAPS_SUCCESS;
    int    i, j, sectionIndex1, sectionIndex2;
    double distLE, distLETotal = 0, *b = NULL;
    int Nspan;
    int NspanMax, imax;
    int NspanMin, imin;

    int numSeg = numSection-1;

    if (numSeg == 0) {
      AIM_ERROR(aimInfo, "VLM must have at least 2 sections\n");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // special case for just one segment (2 sections)
    if (numSeg == 1) {
        sectionIndex1 = vlmSection[0].sectionIndex;

        // use any specified counts
        if (vlmSection[sectionIndex1].Nspan >= 2)
          return CAPS_SUCCESS;

        // just set the total
        vlmSection[sectionIndex1].Nspan = NspanTotal;
        return CAPS_SUCCESS;
    }

    // length of each span section
    b = (double*) EG_alloc(numSeg*sizeof(double));

    // go over all but the last section
    for (i = 0; i < numSection-1; i++) {

        // get the section indices
        sectionIndex1 = vlmSection[i  ].sectionIndex;
        sectionIndex2 = vlmSection[i+1].sectionIndex;

        // skip sections explicitly specified
        if (vlmSection[sectionIndex1].Nspan > 1) continue;

        // use the y-z distance between leading edge points to scale the number of spanwise points
        distLE = 0;
        for (j = 1; j < 3; j++) {
            distLE += pow(vlmSection[sectionIndex2].xyzLE[j] - vlmSection[sectionIndex1].xyzLE[j], 2);
        }
        distLE = sqrt(distLE);

        b[i] = distLE;
        distLETotal += distLE;
    }

    // set the number of spanwise points
    for (i = 0; i < numSection-1; i++) {

        // get the section indices
        sectionIndex1 = vlmSection[i].sectionIndex;

        // skip sections explicitly specified
        if (vlmSection[sectionIndex1].Nspan > 1) continue;

        b[i] /= distLETotal;
        Nspan = NINT(b[i]*abs(NspanTotal));

        vlmSection[sectionIndex1].Nspan = Nspan > 1 ? Nspan : 1;
    }

    // make sure the total adds up
    do {

        Nspan = 0;
        NspanMax = 0;
        NspanMin = NspanTotal;
        imax = 0;
        imin = 0;
        for (i = 0; i < numSection-1; i++) {
            // get the section indices
            sectionIndex1 = vlmSection[i].sectionIndex;

            if ( vlmSection[sectionIndex1].Nspan > NspanMax ) {
                NspanMax = vlmSection[sectionIndex1].Nspan;
                imax = sectionIndex1;
            }
            if ( vlmSection[sectionIndex1].Nspan < NspanMin ) {
                NspanMin = vlmSection[sectionIndex1].Nspan;
                imin = sectionIndex1;
            }

            Nspan += vlmSection[sectionIndex1].Nspan;
        }

        if (Nspan > NspanTotal) {
            vlmSection[imax].Nspan--;

            if (vlmSection[imax].Nspan == 0) {
                AIM_ERROR(aimInfo, "Insufficient spanwise sections! Increase numSpanTotal or numSpanPerSection!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }
        if (Nspan < NspanTotal) {
            vlmSection[imin].Nspan++;
        }

    } while (Nspan != NspanTotal);

    status = CAPS_SUCCESS;
cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in vlm_equalSpaceSpanPanels, status = %d\n", status);

    EG_free(b); b = NULL;

    return status;
}


static int
curvatureArcLenSeg(void *aimInfo, const ego geom, double t1, double t2, double *arc)
{
  int     i, status = CAPS_SUCCESS;
  double  t, s, k, d, ur, mid, *d1, *d2, dir[3];
  double  result[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
/*
  static int     ngauss   = 5;
  static double  wg[2*5]  = { 0.5688888888888889,  0.0000000000000000,
                              0.4786286704993665, -0.5384693101056831,
                              0.4786286704993665,  0.5384693101056831,
                              0.2369268850561891, -0.9061798459386640,
                              0.2369268850561891,  0.9061798459386640 };
 */
/*
  // degree 23 polynomial; 12 points
  static int     ngauss   = 12;
  static double  wg[2*20] = {0.0471753363865118271946160, -0.9815606342467192506905491,
                             0.1069393259953184309602547, -0.9041172563704748566784659,
                             0.1600783285433462263346525, -0.7699026741943046870368938,
                             0.2031674267230659217490645, -0.5873179542866174472967024,
                             0.2334925365383548087608499, -0.3678314989981801937526915,
                             0.2491470458134027850005624, -0.1252334085114689154724414,
                             0.2491470458134027850005624,  0.1252334085114689154724414,
                             0.2334925365383548087608499,  0.3678314989981801937526915,
                             0.2031674267230659217490645,  0.5873179542866174472967024,
                             0.1600783285433462263346525,  0.7699026741943046870368938,
                             0.1069393259953184309602547,  0.9041172563704748566784659,
                             0.0471753363865118271946160,  0.9815606342467192506905491 };
*/
/*
  static int     ngauss   = 15;
  static double  wg[2*15] = { 0.2025782419255613,  0.0000000000000000,
                              0.1984314853271116, -0.2011940939974345,
                              0.1984314853271116,  0.2011940939974345,
                              0.1861610000155622, -0.3941513470775634,
                              0.1861610000155622,  0.3941513470775634,
                              0.1662692058169939, -0.5709721726085388,
                              0.1662692058169939,  0.5709721726085388,
                              0.1395706779261543, -0.7244177313601701,
                              0.1395706779261543,  0.7244177313601701,
                              0.1071592204671719, -0.8482065834104272,
                              0.1071592204671719,  0.8482065834104272,
                              0.0703660474881081, -0.9372733924007060,
                              0.0703660474881081,  0.9372733924007060,
                              0.0307532419961173, -0.9879925180204854,
                              0.0307532419961173,  0.9879925180204854 };
*/
  /* degree 39 polynomial; 20 points */
  static int     ngauss   = 20;
  static double  wg[2*20] = {0.0176140071391521183118620, -0.9931285991850949247861224,
                             0.0406014298003869413310400, -0.9639719272779137912676661,
                             0.0626720483341090635695065, -0.9122344282513259058677524,
                             0.0832767415767047487247581, -0.8391169718222188233945291,
                             0.1019301198172404350367501, -0.7463319064601507926143051,
                             0.1181945319615184173123774, -0.6360536807265150254528367,
                             0.1316886384491766268984945, -0.5108670019508270980043641,
                             0.1420961093183820513292983, -0.3737060887154195606725482,
                             0.1491729864726037467878287, -0.2277858511416450780804962,
                             0.1527533871307258506980843, -0.0765265211334973337546404,
                             0.1527533871307258506980843,  0.0765265211334973337546404,
                             0.1491729864726037467878287,  0.2277858511416450780804962,
                             0.1420961093183820513292983,  0.3737060887154195606725482,
                             0.1316886384491766268984945,  0.5108670019508270980043641,
                             0.1181945319615184173123774,  0.6360536807265150254528367,
                             0.1019301198172404350367501,  0.7463319064601507926143051,
                             0.0832767415767047487247581,  0.8391169718222188233945291,
                             0.0626720483341090635695065,  0.9122344282513259058677524,
                             0.0406014298003869413310400,  0.9639719272779137912676661,
                             0.0176140071391521183118620,  0.9931285991850949247861224 };

  *arc   = 0.0;
  ur     =      t2 - t1;
  mid    = 0.5*(t2 + t1);
  for (i = 0; i < ngauss; i++) {
    t    = 0.5*wg[2*i+1]*ur + mid;
    status = EG_evaluate(geom, &t, result);
    AIM_STATUS(aimInfo, status);

    // tangent magnitude
    s = sqrt(result[3]*result[3] + result[4]*result[4] + result[5]*result[5]);

    // curvature k
    d1        = &result[3];
    d2        = &result[6];
    CROSS(dir, d1, d2);
    d         = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    k         = d/(s*s*s);
    if (k == 0) k = 1;

    // cbrt curvature weighted arc-length
    *arc += cbrt(k)*s*wg[2*i];
  }
  *arc   *= 0.5*ur;

cleanup:
  return status;
}


typedef struct {
  double u; // curvature weighted arc length space
  double t; // parameter space
} CurvatureSpace;


static int
curvatureArcLen(void *aimInfo, const ego geom, double t1, double t2, int *nseg_out, CurvatureSpace **segs_out)
{
  int    status = CAPS_SUCCESS;
  int    i, nseg=0, degree, end;
  int    *header=NULL, oclass, mtype;
  double t, *data=NULL, arc;
  CurvatureSpace *segs=NULL;
  ego ref;

  *nseg_out = 0;
  *segs_out = NULL;

  if (geom->mtype != BSPLINE) {

    AIM_ALLOC(segs, 2, CurvatureSpace, aimInfo, status);
    segs[0].u = 0;
    segs[0].t = t1;
    nseg = 1;

    status = curvatureArcLenSeg(aimInfo, geom, t1, t2, &arc);
    AIM_STATUS(aimInfo, status);
    segs[nseg].u = arc + segs[nseg-1].u;
    segs[nseg].t = t2;
    nseg++;

  } else {

    status = EG_getGeometry(geom, &oclass, &mtype, &ref, &header, &data);
    AIM_STATUS(aimInfo, status);

    /* get length of each set of knots */
    t      = t1;
    degree = header[1];
    end    = header[3]-degree-1;

    AIM_ALLOC(segs, end - degree + 1, CurvatureSpace, aimInfo, status);
    segs[0].u = 0;
    segs[0].t = t1;
    nseg = 1;

    for (i = degree; i < end; i++) {
      if (data[i] <= t)  continue;
      if (data[i] >= t2) break;
      status = curvatureArcLenSeg(aimInfo, geom, t, data[i], &arc);
      AIM_STATUS(aimInfo, status);
      t = data[i];
      segs[nseg].u = arc + segs[nseg-1].u;
      segs[nseg].t = t;
      nseg++;
    }

    status = curvatureArcLenSeg(aimInfo, geom, t, t2, &arc);
    AIM_STATUS(aimInfo, status);
    segs[nseg].u = arc + segs[nseg-1].u;
    segs[nseg].t = t2;
    nseg++;
  }

  status = CAPS_SUCCESS;

  *nseg_out = nseg;
  *segs_out = segs;

cleanup:
  if (status != CAPS_SUCCESS) {
    AIM_FREE(*segs_out);
  }

  AIM_FREE(header);
  AIM_FREE(data);

  return status;
}


// Get curvature weighted arc-length based point counts on each edge of a section
static
int vlm_secEdgePoints(void *aimInfo,
                      int numPoint,
                      int numEdge, ego *edges,
                      ego teObj,
                      int **numEdgePointsOut, // Number of interior points on an Edge (i.e. excluding the Nodes)
                      int **numEdgeSegsOut,
                      CurvatureSpace ***edgeSegsOut)
{
    int status; // Function return status

    int i, j; // Indexing
    double totLen, arcLen;

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    double trange[4];
    ego ref, *children = NULL;

    int numPointTot;
    int *numEdgePoint = NULL, *numEdgeSegs = NULL;;
    CurvatureSpace **edgeSegs = NULL;

    *numEdgePointsOut = NULL;

    // weight the number of points on each edge based on the curvature weighted arc length
    AIM_ALLOC(numEdgePoint, numEdge, int, aimInfo, status);
    for (i = 0; i < numEdge; i++) numEdgePoint[i] = 0;

    AIM_ALLOC(edgeSegs, numEdge, CurvatureSpace*, aimInfo, status);
    for (i = 0; i < numEdge; i++) edgeSegs[i] = NULL;
    AIM_ALLOC(numEdgeSegs, numEdge, int, aimInfo, status);
    for (i = 0; i < numEdge; i++) numEdgeSegs[i] = 0;

    totLen = 0.0;
    for (i = 0; i < numEdge; i++) {
        numEdgePoint[i] = 0;

        if ( teObj == edges[i] ) continue; // don't count the trailing edge
        status = EG_getTopology(edges[i], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);
        if (mtype == DEGENERATE) continue;

        status = curvatureArcLen( aimInfo, ref, trange[0], trange[1], &numEdgeSegs[i], &edgeSegs[i] );
        AIM_STATUS(aimInfo, status);

        totLen += edgeSegs[i][numEdgeSegs[i]-1].u;
    }

    numPointTot = 1; // One because the airfoil coordinates are an open loop
    for (i = 0; i < numEdge; i++) {

        if ( teObj == edges[i] ) continue; // don't count the trailing edge
        if (edges[i]->mtype == DEGENERATE) continue;

        AIM_NOTNULL(edgeSegs[i], aimInfo, status);
        arcLen = edgeSegs[i][numEdgeSegs[i]-1].u;

        numEdgePoint[i] = MAX(0,numPoint*arcLen/totLen - 1); // Number of segments, and remove 1 to get Edge interior points
        numPointTot += numEdgePoint[i] + 1; // Account for the Node in total point count with +1
    }

    // adjust any rounding so the total number of points matches maxNumPoint
    while (numPointTot != numPoint) {
        j = 0;
        if (numPointTot > numPoint) {
            // remove one point from the largest count
            for (i = 0; i < numEdge; i++) {
                if ( teObj == edges[i] ) continue; // don't count the trailing edge
                if (edges[i]->mtype == DEGENERATE) continue;
                if (numEdgePoint[i] > numEdgePoint[j]) j = i;
            }
            numEdgePoint[j]--;
            numPointTot--;
        } else {
            // add one point to the smallest edge count
            for (i = 0; i < numEdge; i++) {
                if ( teObj == edges[i] ) continue; // don't count the trailing edge
                if (edges[i]->mtype == DEGENERATE) continue;
                if (numEdgePoint[i] < numEdgePoint[j]) j = i;
            }
            numEdgePoint[j]++;
            numPointTot++;
        }
    }


    *numEdgePointsOut = numEdgePoint;
    *numEdgeSegsOut = numEdgeSegs;
    *edgeSegsOut = edgeSegs;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
        AIM_FREE(numEdgePoint);
        AIM_FREE(numEdgeSegs);
        if (edgeSegs != NULL) {
            for (i = 0; i < numEdge; i++) {
                AIM_FREE(edgeSegs[i]);
            }
            AIM_FREE(edgeSegs);
        }
    }

    return status;
}

// Retrieve edge ordering such that the loop starts at the trailing edge NODE
// with the teObj last if it is an EDGE
int vlm_secOrderEdges(
         /*@unused@*/ int numNode, ego *nodes,
                      int numEdge, ego *edges,
                      ego body, ego teObj,
                      int **edgeLoopOrderOut,
                      int **edgeLoopSenseOut,
                      ego *nodeTEOut)
{
    int status; // Function return status

    int i; // Indexing

    int numLoop;
    ego *loops=NULL;

    int edgeIndex, *edgeLoopOrder=NULL, *edgeLoopSense = NULL;

    int sense = 0, itemp = 0, nodeIndexTE2[2];
    int teClass;
    ego nodeTE;

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    double trange[4];
    ego ref, *children = NULL;

    *edgeLoopOrderOut = NULL;
    *edgeLoopSenseOut = NULL;
    *nodeTEOut = NULL;


    status = EG_getBodyTopos(body, NULL, LOOP, &numLoop, &loops);
    if (status != EGADS_SUCCESS) {
        printf("\tError in vlm_secOrderEdges, getBodyTopos Loops = %d\n", status);
        goto cleanup;
    }

    // Get the NODE(s) indexing for the trailing edge
    status = EG_getTopology(teObj, &ref, &teClass, &mtype, trange, &numChildren, &children, &sens);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (teClass == NODE) {
      nodeIndexTE2[0] = EG_indexBodyTopo(body, teObj);
      nodeIndexTE2[1] = nodeIndexTE2[0];
    } else {
      nodeIndexTE2[0] = EG_indexBodyTopo(body, children[0]);
      nodeIndexTE2[1] = EG_indexBodyTopo(body, children[1]);
    }

    // Determine edge order
    edgeLoopSense = (int *) EG_alloc(numEdge*sizeof(int));
    edgeLoopOrder = (int *) EG_alloc(numEdge*sizeof(int));
    if ( (edgeLoopSense == NULL) || (edgeLoopOrder == NULL) ) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    status = EG_getTopology(loops[0], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
    if (status != EGADS_SUCCESS) goto cleanup;

    // Get the edge ordering in the loop
    // The first edge may not start at the trailing edge
    for (i = 0; i < numChildren; i++) {
        edgeIndex = EG_indexBodyTopo(body, children[i]);

        if (edgeIndex < EGADS_SUCCESS) {
            status = CAPS_BADINDEX;
            goto cleanup;
        }

        edgeLoopOrder[i] = edgeIndex;
        edgeLoopSense[i] = sens[i];
    }

    // Reorder edge indexing such that a trailing edge node is the first node in the loop
    while ( (int) true ) {

        // the first edge cannot be the TE edge
        if (edges[edgeLoopOrder[0]-1] != teObj) {

            status = EG_getTopology(edges[edgeLoopOrder[0]-1], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
            if (status != EGADS_SUCCESS) goto cleanup;

            if (mtype == DEGENERATE) continue;

            // Get the sense of the edge from the loop
            sense = edgeLoopSense[0];

            // check if the starting child node is one of the TE nodes
            if (sense == 1) {
                if ( children[0] == nodes[nodeIndexTE2[0]-1] || children[0] == nodes[nodeIndexTE2[1]-1] ) {
                    nodeTE = children[0];
                    break;
                }
            } else {
                if ( children[1] == nodes[nodeIndexTE2[0]-1] || children[1] == nodes[nodeIndexTE2[1]-1] ) {
                    nodeTE = children[1];
                    break;
                }
            }
        }

        // rotate the avl order and the edge sense to the left by one
        itemp = edgeLoopOrder[0];
        for (i = 0; i < numEdge - 1; i++)
            edgeLoopOrder[i] = edgeLoopOrder[i + 1];
        edgeLoopOrder[i] = itemp;

        itemp = edgeLoopSense[0];
        for (i = 0; i < numEdge - 1; i++)
          edgeLoopSense[i] = edgeLoopSense[i + 1];
        edgeLoopSense[i] = itemp;
    }

    if (teClass == EDGE && teObj != edges[edgeLoopOrder[numEdge-1]-1]) {
        printf("Developer ERROR: Found trailing edge but it's not the last edge in the loop!!!!\n");
        status = CAPS_SOURCEERR;
        goto cleanup;
    }

    *edgeLoopOrderOut = edgeLoopOrder;
    *edgeLoopSenseOut = edgeLoopSense;
    *nodeTEOut = nodeTE;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
        AIM_FREE(edgeLoopOrder);
        AIM_FREE(edgeLoopSense);
    }

    AIM_FREE(loops);

    return status;
}

// Get the airfoil cross-section tessellation ego given a vlmSectionStruct
int vlm_getSectionTessSens(void *aimInfo,
                           vlmSectionStruct *vlmSection,
                           int normalize,      // Normalize by chord (true/false)
                           const char *geomInName,
                           const int irow, const int icol,
                           ego tess,
                           double **dx_dvar_out,
                           double **dy_dvar_out)
{
    int status; // Function return status

    int i, j, counter=0; // Indexing

    int n, edgeIndex, *edgeLoopOrder=NULL, *edgeLoopSense = NULL;

    int numPoint, numEdge, numNode, numChildren;
    ego *nodes = NULL, *edges = NULL, *children = NULL;

    int sense = 0, state;
    double *ts=NULL, *xyz=NULL, trange[4], result[18];
    double *dxyz=NULL;
    double *dx_dvar=NULL, *dy_dvar=NULL;
    double *xyzLE, *xyzTE;

    double chord;
    double xdot[3], ydot[3], *secnorm;

    //EGADS returns
    int oclass, mtype, *sens = NULL;
    ego ref, nodeTE = NULL;
    ego teObj = NULL, body;

    *dx_dvar_out = NULL;
    *dy_dvar_out = NULL;

    body = vlmSection->ebody;
    chord = vlmSection->chord;
    secnorm = vlmSection->normal;
    xyzLE = vlmSection->xyzLE;
    xyzTE = vlmSection->xyzTE;
    teObj = vlmSection->teObj;

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    // Get the loop edge ordering so it starts at the trailing edge NODE
    status = vlm_secOrderEdges(numNode, nodes,
                               numEdge, edges,
                               body, teObj,
                               &edgeLoopOrder, &edgeLoopSense, &nodeTE);
    AIM_STATUS(aimInfo, status);

    // get the total number of points
    status = EG_statusTessBody(tess, &body, &state, &numPoint);
    AIM_STATUS(aimInfo, status);
    numPoint += 1; // Because the airfoil is an open loop

    // vector from LE to TE normalized
    xdot[0]  =  xyzTE[0] - xyzLE[0];
    xdot[1]  =  xyzTE[1] - xyzLE[1];
    xdot[2]  =  xyzTE[2] - xyzLE[2];

    if (normalize == (int) false) chord = 1;

    xdot[0] /=  chord;
    xdot[1] /=  chord;
    xdot[2] /=  chord;

    // cross with section PLANE normal to get perpendicular vector in the PLANE
    cross_DoubleVal(secnorm, xdot, ydot);

    AIM_ALLOC(dx_dvar, numPoint, double, aimInfo, status);
    AIM_ALLOC(dy_dvar, numPoint, double, aimInfo, status);

    status = aim_setSensitivity(aimInfo, geomInName, irow, icol);
    AIM_STATUS(aimInfo, status);

    // get the coordinate of the starting trailing edge node
//    status = EG_evaluate(nodeTE, NULL, result);
//    AIM_STATUS(aimInfo, status);
//
//    result[0] -= xyzLE[0];
//    result[1] -= xyzLE[1];
//    result[2] -= xyzLE[2];

    // get the sensitivity of the starting trailing edge node
    status = aim_getSensitivity(aimInfo, tess, 0, EG_indexBodyTopo(body,nodeTE),
                                &n, &dxyz);
    AIM_STATUS(aimInfo, status);

    result[0] = dxyz[0];
    result[1] = dxyz[1];
    result[2] = dxyz[2];

    dx_dvar[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
    dy_dvar[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
    counter++;
    AIM_FREE(dxyz);

    // Loop through edges based on order
    for (i = 0; i < numEdge; i++) {
        //printf("Edge order %d\n", edgeOrder[i]);

        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        if (edges[edgeIndex] == teObj) continue;

        // Get children for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        if (status != EGADS_SUCCESS) goto cleanup;

        // Get the sense of the edge from the loop
        sense = edgeLoopSense[edgeIndex];

        status = aim_getSensitivity(aimInfo, tess, 1, edgeIndex+1, &n, &dxyz);
        AIM_STATUS(aimInfo, status);

        for (j = 1; j < n-1; j++) {
          if (sense == SFORWARD) {
            result[0] = dxyz[3*j+0];
            result[1] = dxyz[3*j+1];
            result[2] = dxyz[3*j+2];
          } else {
            result[0] = dxyz[3*(n-1-j)+0];
            result[1] = dxyz[3*(n-1-j)+1];
            result[2] = dxyz[3*(n-1-j)+2];
          }

          dx_dvar[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
          dy_dvar[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
          counter++;
        }
        AIM_FREE(dxyz);

        // get the last Node on the Edge
        if (sense == SFORWARD) {
            status = aim_getSensitivity(aimInfo, tess, 0, EG_indexBodyTopo(body,children[1]),
                                        &n, &dxyz);
        } else {
            status = aim_getSensitivity(aimInfo, tess, 0, EG_indexBodyTopo(body,children[0]),
                                        &n, &dxyz);
        }
        AIM_STATUS(aimInfo, status);

        result[0] = dxyz[0];
        result[1] = dxyz[1];
        result[2] = dxyz[2];

        dx_dvar[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
        dy_dvar[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        counter++;
        AIM_FREE(dxyz);
    }

    *dx_dvar_out = dx_dvar;
    *dy_dvar_out = dy_dvar;

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(ts);
    AIM_FREE(xyz);
    AIM_FREE(dxyz);
    AIM_FREE(nodes);
    AIM_FREE(edges);

    AIM_FREE(edgeLoopSense);
    AIM_FREE(edgeLoopOrder);

    return status;
}

// Get the airfoil cross-section coordinates given a vlmSectionStruct
int vlm_getSectionCoord(void *aimInfo,
                        vlmSectionStruct *vlmSection,
                        int normalize,      // Normalize by chord (true/false)
                        int numPoint,       // number of points in airfoil
                        double **xCoordOut, // [numPoint]
                        double **yCoordOut, // [numPoint] for upper and lower surface
                        ego *tessOut)       // Tess object that created points
{
    int status; // Function return status

    int i, j, k; // Indexing

    int counter=0;

    int edgeIndex, *edgeLoopOrder=NULL, *edgeLoopSense = NULL;

    int *numEdgePoint = NULL, *numEdgeSegs = NULL;;
    CurvatureSpace **edgeSegs = NULL;

    double chord;
    double xdot[3], ydot[3], *secnorm;
    double params[3];

    int numEdge, numNode;
    ego *nodes = NULL, *edges = NULL, nodeTE = NULL;

    int sense = 0;
    double s, u, du, t, tt[2], trange[4], result[18];

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    ego ref, *children = NULL;

    ego teObj = NULL, body, tess=NULL;

    double xyzLE[3], xyzTE[3];
    int n;
    const double *xyz, *ts;
    double *xyzS=NULL, *tS=NULL;

    double *xCoord=NULL, *yCoord=NULL;

    *xCoordOut = NULL;
    *yCoordOut = NULL;
    *tessOut = NULL;

    body = vlmSection->ebody;
    chord = vlmSection->chord;
    secnorm = vlmSection->normal;
    xyzLE[0] = vlmSection->xyzLE[0];
    xyzLE[1] = vlmSection->xyzLE[1];
    xyzLE[2] = vlmSection->xyzLE[2];
    xyzTE[0] = vlmSection->xyzTE[0];
    xyzTE[1] = vlmSection->xyzTE[1];
    xyzTE[2] = vlmSection->xyzTE[2];
    teObj = vlmSection->teObj;

//#define DUMP_EGADS_SECTIONS
#ifdef DUMP_EGADS_SECTIONS
    static int ID = 0;
    int atype, alen;
    const int *ints=NULL;
    const double *reals=NULL;
    const char *string=NULL;
    char filename[256];

    status = EG_attributeRet(body, "_name", &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      snprintf(filename, 256, "section_%d_%s.egads", ID, string);
    } else {
      snprintf(filename, 256, "section%d.egads", ID);
    }
    /* make a model and write it out */
    remove(filename);
    printf(" EG_saveModel(%s) = %d\n", filename, EG_saveModel(body, filename));
    EG_deleteObject(newModel);
    ID++;
#endif

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    // Get the number of points on each edge
    status = vlm_secEdgePoints(aimInfo,
                               numPoint,
                               numEdge, edges,
                               teObj, &numEdgePoint,
                               &numEdgeSegs, &edgeSegs);
    AIM_STATUS(aimInfo, status);

    // Get the loop edge ordering so it starts at the trailing edge NODE
    status = vlm_secOrderEdges(numNode, nodes,
                               numEdge, edges,
                               body, teObj,
                               &edgeLoopOrder, &edgeLoopSense, &nodeTE);
    AIM_STATUS(aimInfo, status);

    // initialize the tessellation
    status = EG_initTessBody(body, &tess);
    AIM_STATUS(aimInfo, status);

      // Loop through edges
    for (edgeIndex = 0; edgeIndex < numEdge; edgeIndex++) {

        // The whole body needs to be tessellated, so just set 0 points on the TE
        if (edges[edgeIndex] == teObj) {
          status = EG_attributeAdd(edges[edgeIndex], ".nPos", ATTRINT, 1, &numEdgePoint[edgeIndex], NULL, NULL);
          AIM_STATUS(aimInfo, status);
          continue;
        }

        // Get t-range and nodes for the edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        AIM_STATUS(aimInfo, status);
        if (mtype == DEGENERATE) continue;

        // Adjust the interior edge points to total edge points by adding the end nodes
        numEdgePoint[edgeIndex] += 2;

        AIM_REALL(tS  ,   numEdgePoint[edgeIndex], double, aimInfo, status);
        AIM_REALL(xyzS, 3*numEdgePoint[edgeIndex], double, aimInfo, status);

        // Uniform spacing in curvature weighted arc length
        k = numEdgeSegs[edgeIndex]-1;
        du = edgeSegs[edgeIndex][k].u/(numEdgePoint[edgeIndex]-1);

        // Create in points along edge
        k = 0;
        for (j = 0; j < numEdgePoint[edgeIndex]; j++) {

            if (j == 0) {
                status = EG_evaluate(nodes[0], NULL, result);
                AIM_STATUS(aimInfo, status);
                t = trange[0];
            } else if (j == numEdgePoint[edgeIndex]-1) {
                status = EG_evaluate(nodes[1], NULL, result);
                AIM_STATUS(aimInfo, status);
                t = trange[1];
            } else {

                u = j*du;

                while (k < numEdgeSegs[edgeIndex]) {
                  if (u < edgeSegs[edgeIndex][k+1].u && u >= edgeSegs[edgeIndex][k].u)
                    break;
                  k++;
                }

                // interpolate t based on u-space
                tt[0] = edgeSegs[edgeIndex][k  ].t;
                tt[1] = edgeSegs[edgeIndex][k+1].t;
                s = ( u - edgeSegs[edgeIndex][k].u ) / ( edgeSegs[edgeIndex][k+1].u - edgeSegs[edgeIndex][k].u );

                t = tt[0] + s*(tt[1]-tt[0]);

                status = EG_evaluate(edges[edgeIndex], &t, result);
                AIM_STATUS(aimInfo, status);
           }

            tS[j] = t;
            xyzS[3*j+0] = result[0];
            xyzS[3*j+1] = result[1];
            xyzS[3*j+2] = result[2];
        }

        status = EG_setTessEdge(tess, EG_indexBodyTopo(body, edges[edgeIndex]), numEdgePoint[edgeIndex], xyzS, tS);
        AIM_STATUS(aimInfo, status);
    }

    // close the tessellation. Use 0 length to prevent face points
    params[0] = 0;
    params[1] = chord;
    params[2] = 20;
    status = EG_finishTess(tess, params);
    AIM_STATUS(aimInfo, status);


    if (normalize == (int) true) {

      // vector from LE to TE normalized
      xdot[0]  =  xyzTE[0] - xyzLE[0];
      xdot[1]  =  xyzTE[1] - xyzLE[1];
      xdot[2]  =  xyzTE[2] - xyzLE[2];

      xdot[0] /= chord;
      xdot[1] /= chord;
      xdot[2] /= chord;

    } else {

      // Remove offset, scaling, and set xdot to x-direction
      chord = 1;

      xyzLE[0] = 0.0;
      xyzLE[1] = 0.0;
      xyzLE[2] = 0.0;

      xdot[0] = 1.0;
      xdot[1] = 0.0;
      xdot[2] = 0.0;
    }

    // cross with section PLANE normal to get perpendicular vector in the PLANE
    cross_DoubleVal(secnorm, xdot, ydot);

//#define DUMP_TESS_SECTIONS
#ifdef DUMP_TESS_SECTIONS
    static int ID = 0;
    int atype, alen;
    const int *ints=NULL;
    const double *reals=NULL;
    const char *string=NULL;
    char filename[256];
    ego newModel, context, bodies[2];

    status = EG_attributeRet(body, "_name", &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      snprintf(filename, 256, "section_%d_%s.egads", ID, string);
    } else {
      snprintf(filename, 256, "section%d.egads", ID);
    }

    EG_getContext(body, &context);

    EG_copyObject(body, NULL, &bodies[0]);
    EG_copyObject(tess, bodies[0], &bodies[1]);
    status = EG_makeTopology(context, NULL, MODEL, 2, NULL, 1, bodies, NULL, &newModel);

    /* make a model and write it out */
    remove(filename);
    printf(" EG_saveModel(%s) = %d\n", filename, EG_saveModel(newModel, filename));
    EG_deleteObject(newModel);
    ID++;

    status = EG_openTessBody(tess);
    AIM_STATUS(aimInfo, status);

    params[0] = 10*chord;
    params[1] = chord;
    params[2] = 20;
    status = EG_finishTess(tess, params);
    AIM_STATUS(aimInfo, status);
#endif

    // set output points

    AIM_ALLOC(xCoord, numPoint, double, aimInfo, status);
    AIM_ALLOC(yCoord, numPoint, double, aimInfo, status);

    for (i = 0; i < numPoint; i++) {
        xCoord[i] = 0.0;
        yCoord[i] = 0.0;
    }

    // get the coordinate of the starting trailing edge node
    status = EG_evaluate(nodeTE, NULL, result);
    AIM_STATUS(aimInfo, status);

    result[0] -= xyzLE[0];
    result[1] -= xyzLE[1];
    result[2] -= xyzLE[2];

    xCoord[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
    yCoord[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
    counter += 1;

    // Loop through edges based on order
    for (i = 0; i < numEdge; i++) {
        //printf("Edge order %d\n", edgeLoopOrder[i]);

        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        if (edges[edgeIndex] == teObj) continue;

        // Get t- range for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        AIM_STATUS(aimInfo, status);

        // Get the sense of the edge from the loop
        sense = edgeLoopSense[i];

        // get the loop edge tessellation
        status = EG_getTessEdge(tess, edgeIndex+1, &n, &xyz, &ts);
        AIM_STATUS(aimInfo, status);

        // Write out in points along each edge
        for (j = 1; j < n-1; j++) {
            if (sense == SFORWARD) {
                result[0] = xyz[3*j+0];
                result[1] = xyz[3*j+1];
                result[2] = xyz[3*j+2];
            } else {
                result[0] = xyz[3*(n-1-j)+0];
                result[1] = xyz[3*(n-1-j)+1];
                result[2] = xyz[3*(n-1-j)+2];
            }

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            xCoord[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
            yCoord[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
            counter += 1;
        }

        // Write the last node of the edge
        if (sense == SFORWARD) {
            status = EG_evaluate(children[1], NULL, result);
        } else {
            status = EG_evaluate(children[0], NULL, result);
        }
        AIM_STATUS(aimInfo, status);

        result[0] -= xyzLE[0];
        result[1] -= xyzLE[1];
        result[2] -= xyzLE[2];

        xCoord[counter] = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
        yCoord[counter] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        counter += 1;
    }
    if (counter != numPoint) {
      AIM_ERROR(aimInfo, "Development error counter > *numPoint!");
      status = CAPS_NOTIMPLEMENT;
      goto cleanup;
    }


    for (i = 0; i < counter; i++) {
        if (fabs(xCoord[i]) < 1e-12) {
            xCoord[i] = 0.0;
        }
        if (fabs(yCoord[i]) < 1e-12) {
            yCoord[i] = 0.0;
        }
    }

    // counter may be lower than numPoint if there are points on the teObj
    *xCoordOut = xCoord;
    *yCoordOut = yCoord;
    *tessOut = tess;

    xCoord = NULL;
    yCoord = NULL;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
        AIM_FREE(xCoord);
        AIM_FREE(yCoord);
        EG_deleteObject(tess);
    }

    if (edgeSegs != NULL) {
      for (i = 0; i < numEdge; i++)
        AIM_FREE(edgeSegs[i]);
      AIM_FREE(edgeSegs);
    }
    AIM_FREE(numEdgeSegs);
    AIM_FREE(numEdgePoint);

    AIM_FREE(xyzS);
    AIM_FREE(tS);

    AIM_FREE(nodes);
    AIM_FREE(edges);

    AIM_FREE(edgeLoopSense);
    AIM_FREE(edgeLoopOrder);

    return status;
}


// Write out the airfoil cross-section given an ego body
int vlm_writeSection(void *aimInfo,
                     FILE *fp,
                     vlmSectionStruct *vlmSection,
                     int normalize, // Normalize by chord (true/false)
                     int numPoint)  // max number of points in airfoil
{
    int status; // Function return status

    int i; // Indexing

    ego tess=NULL;
    double *xCoord=NULL, *yCoord=NULL;

    status = vlm_getSectionCoord(aimInfo,
                                 vlmSection,
                                 normalize,
                                 numPoint,
                                 &xCoord,
                                 &yCoord,
                                 &tess);
    AIM_STATUS(aimInfo, status);

    for( i = 0; i < numPoint; i++) {
        fprintf(fp, "%16.12e %16.12e\n", xCoord[i], yCoord[i]);
    }

    fprintf(fp, "\n");

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(xCoord);
    AIM_FREE(yCoord);

    EG_deleteObject(tess);

    return status;
}

static
void spacer( int N, double pspace, double scale, double x[]) {
/*  modified from avl source sgutil.f
 *
 *...PURPOSE     TO CALCULATE A NORMALIZED (0<=X<=1) SPACING ARRAY.
 *
 *...INPUT       N      =  NUMBER OF DESIRED POINTS IN ARRAY.
 *               PSPACE =  SPACING PARAMETER (-3<=PSPACE<=3).
 *                         DEFINES POINT DISTRIBUTION
 *                         TO BE USED AS FOLLOWS:
 *                 PSPACE = 0  : EQUAL SPACING
 *                 PSPACE = 1  : COSINE SPACING.
 *                 PSPACE = 2  : SINE SPACING
 *                               (CONCENTRATING POINTS NEAR 0).
 *                 PSPACE = 3  : EQUAL SPACING.
 *
 *                 NEGATIVE VALUES OF PSPACE PRODUCE SPACING
 *                 WHICH IS REVERSED (AFFECTS ONLY SINE SPACING).
 *                 INTERMEDIATE VALUES OF PSPACE WILL PRODUCE
 *                 A SPACING WHICH IS A LINEAR COMBINATION
 *                 OF THE CORRESPONDING INTEGER VALUES.
 *
 *...OUTPUT      X      =  NORMALIZED SPACING ARRAY (0 <= X <= 1)
 *                         THE FIRST ELEMENT WILL ALWAYS BE  X[0  ] = 0.
 *                         THE LAST ELEMENT WILL ALWAYS BE   X[N-1] = scale
 *
 */

  double pabs, pequ, pcos, psin;
  double frac, theta;
  int nabs, i;

  pabs = fabs(pspace);
  if (pabs > 3.0) pabs = 3;
  nabs = (int) pabs + 1;

  if (nabs == 1) {
    pequ   = 1.-pabs;
    pcos   = pabs;
    psin   = 0.;
  } else if (nabs == 2) {
    pequ   = 0.;
    pcos   = 2.-pabs;
    psin   = pabs-1.;
  } else {
    pequ   = pabs-2.;
    pcos   = 0.;
    psin   = 3.-pabs;
  }

  for (i = 0; i < N; i++) {
    frac = (double)i/(double)(N-1.);
    theta =  frac * PI;
    if (pspace >= 0. )
      x[i] = (pequ * frac
           +  pcos * ( 1. - cos ( theta )      ) / 2.
           +  psin * ( 1. - cos ( theta / 2. ) ) )*scale;
    if (pspace < 0. )
       x[i] = (pequ * frac
            +  pcos * ( 1. - cos ( theta )     ) / 2.
            +  psin * sin ( theta / 2. )         )*scale;
  }

  x[0]   = 0.;
  x[N-1] = scale;
}
// Use Newton solve to refine the t-value
static int _refineT(void *aimInfo, double x1, double xCoord, double scale, double chord,
                    ego edge, double xdot[], double result[], double xyzLE[],
                    double *t) {

    int status;

    int tries;
    double tempT, x_t, residual, deltaT;

    while( fabs(x1 - xCoord) > 1e-7*scale ) {

        x_t = (xdot[0]*result[3]+xdot[1]*result[4]+xdot[2]*result[5])/chord;

        residual = x1 - xCoord;
        deltaT = residual / x_t;

        tries = 0;
        do {

            tempT = *t - deltaT;

            status = EG_evaluate(edge, &tempT, result);
            if (status != EGADS_SUCCESS) return status;

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            x1 = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;

            deltaT /= 2;

            if (tries++ > 20) {
                AIM_ERROR(aimInfo, "Newton solve did not converge.\n"
                            "There is likely something wrong with the geometry of the airfoil.");
                return CAPS_BADVALUE;
            }

        } while ( fabs(x1 - xCoord) > fabs(residual) );

        *t = tempT;
    }

    return CAPS_SUCCESS;
}

// Get the airfoil cross-section given a vlmSectionStruct
// where y-upper and y-lower correspond to the x-value
// Only works for sharp trailing edges
int vlm_getSectionCoordX(void *aimInfo,
                         vlmSectionStruct *vlmSection,
                         double Cspace,      // Chordwise spacing (see spacer)
                         int normalize,      // Normalize by chord (true/false)
                         int rotated,        // Leave airfoil rotated (true/false)
                         int numPoint,       // Number of points in airfoil
                         double **xCoordOut, // [numPoint] increasing x values
                         double **yUpperOut, // [numPoint] for upper surface
                         double **yLowerOut) // [numPoint] for lower surface
{
    int status; // Function return status

    int i, j, ipnt, jbeg, jend; // Indexing

    int nodeIndexLE;
    int edgeIndex, *edgeLoopOrder=NULL, *edgeLoopSense = NULL;

    int *numEdgePoint = NULL, *numEdgeSegs = NULL;;
    CurvatureSpace **edgeSegs = NULL;

    double chord, x1, x2, ux, uy, uz, c, s;
    double xdot[3], ydot[3], tmp[3], *secnorm;

    int numEdge, numNode;
    ego *nodes = NULL, *edges = NULL, nodeTE = NULL, nodeLE;

    int sense = 0;
    double t, trange[4], result[18], params[3], scale;

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    ego ref, *children = NULL;

    ego teObj = NULL, body, tess=NULL;

    int nlen;
    const double *pxyz, *pt;
    double *xyzLE, *xyzTE;

    double *xCoord=NULL, *yUpper=NULL, *yLower=NULL;

    *xCoordOut = NULL;
    *yUpperOut = NULL;
    *yLowerOut = NULL;

    if (vlmSection->teClass != NODE) {
      AIM_ERROR(aimInfo, "Trailing edge must be sharp!");
      status = CAPS_SHAPEERR;
      goto cleanup;
    }

    body = vlmSection->ebody;
    scale = chord = vlmSection->chord;
    secnorm = vlmSection->normal;
    nodeIndexLE = vlmSection->nodeIndexLE;
    xyzLE = vlmSection->xyzLE;
    xyzTE = vlmSection->xyzTE;
    teObj = vlmSection->teObj;

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    // get the leading edge node
    nodeLE = nodes[nodeIndexLE-1];

    // Get the number of points on each edge
    status = vlm_secEdgePoints(aimInfo, numPoint,
                               numEdge, edges,
                               teObj, &numEdgePoint,
                               &numEdgeSegs, &edgeSegs);
    AIM_STATUS(aimInfo, status);

    // Get the loop edge ordering so it starts at the trailing edge NODE
    status = vlm_secOrderEdges(numNode, nodes,
                               numEdge, edges,
                               body, teObj,
                               &edgeLoopOrder, &edgeLoopSense, &nodeTE);
    AIM_STATUS(aimInfo, status);

    // vector from LE to TE normalized
    xdot[0]  =  xyzTE[0] - xyzLE[0];
    xdot[1]  =  xyzTE[1] - xyzLE[1];
    xdot[2]  =  xyzTE[2] - xyzLE[2];

    if (normalize == (int) false) chord = 1;

    scale /= chord;
    xdot[0] /=  chord;
    xdot[1] /=  chord;
    xdot[2] /=  chord;

    if (rotated == (int)true) {
      c = cos(-vlmSection->ainc*PI/180.);
      s = sin(-vlmSection->ainc*PI/180.);

      ux = secnorm[0];
      uy = secnorm[1];
      uz = secnorm[2];

      // rotation matrix about an arbitrary axis
      tmp[0] = (    c + (1-c)*ux*ux)*xdot[0] + (-uz*s + (1-c)*uy*ux)*xdot[1] + ( uy*s + (1-c)*uz*ux)*xdot[2];
      tmp[1] = ( uz*s + (1-c)*ux*uy)*xdot[0] + (    c + (1-c)*uy*uy)*xdot[1] + (-ux*s + (1-c)*uz*uy)*xdot[2];
      tmp[2] = (-uy*s + (1-c)*ux*uz)*xdot[0] + ( ux*s + (1-c)*uy*uz)*xdot[1] + (    c + (1-c)*uz*uz)*xdot[2];

      xdot[0] = tmp[0];
      xdot[1] = tmp[1];
      xdot[2] = tmp[2];
    }

    // cross with section PLANE normal to get perpendicular vector in the PLANE
    cross_DoubleVal(secnorm, xdot, ydot);

    // create a tessellation object on the edges

    // loop over all edges and set the desired point count
    for (i = 0; i < numEdge; i++) {
      status = EG_attributeAdd(edges[i], ".nPos", ATTRINT, 1, &numEdgePoint[i], NULL, NULL);
      AIM_STATUS(aimInfo, status);
    }

    // Negating the first parameter triggers EGADS to only put vertexes on edges
    params[0] = -chord;
    params[1] =  chord;
    params[2] =  20;

    status = EG_makeTessBody(body, params, &tess);
    AIM_STATUS(aimInfo, status);



    // set output points

    xCoord = (double *) EG_alloc(numPoint*sizeof(double));
    yUpper = (double *) EG_alloc(numPoint*sizeof(double));
    yLower = (double *) EG_alloc(numPoint*sizeof(double));
    if (xCoord == NULL || yUpper == NULL || yLower == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // generate the x-coordinates using AVL algorithm
    spacer( numPoint, Cspace, scale, xCoord);

    // Loop over the upper surface
    ipnt = numPoint-2; // avoid last point due to rounding errors
    for (i = 0; i < numEdge; i++) {
        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        // Get t- range for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        AIM_STATUS(aimInfo, status);

        // Check if the first node is the LE edge, which means the current edge is on the lower surface
        sense = edgeLoopSense[i];
        if (sense == SFORWARD) {
            if (children[0] == nodeLE) break;
        } else {
            if (children[1] == nodeLE) break;
        }

        status = EG_getTessEdge(tess, edgeIndex+1, &nlen, &pxyz, &pt);
        AIM_STATUS(aimInfo, status);

        jbeg = sense == SFORWARD ? 0 : nlen-1;
        jend = sense == SFORWARD ? nlen-1 : 0;
        j = jbeg;
        for (; ipnt > 0; ipnt--) { // avoid first point due to rounding errors

            while (j != jend) {
              x2 = (xdot[0]*(pxyz[3*(j      )+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j      )+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j      )+2]-xyzLE[2]))/chord;
              x1 = (xdot[0]*(pxyz[3*(j+sense)+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j+sense)+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j+sense)+2]-xyzLE[2]))/chord;
              if (x1 <= xCoord[ipnt] && xCoord[ipnt] <= x2) {
                t = (pt[j] + pt[j+sense])/2.;
                break;
              }
              j += sense;
            }
            if (j == jend) { break; }

            status = EG_evaluate(edges[edgeIndex], &t, result);
            AIM_STATUS(aimInfo, status);

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            x1 = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
            // use Newton solve to refine the t-value
            status = _refineT(aimInfo, x1, xCoord[ipnt], scale, chord,
                              edges[edgeIndex], xdot, result, xyzLE, &t);
            AIM_STATUS(aimInfo, status);

            yUpper[ipnt] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        }
    }

    // Loop over the lower surface
    ipnt = 1;  // avoid first point due to rounding errors
    for (; i < numEdge; i++) {
        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        // Get t- range for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        AIM_STATUS(aimInfo, status);

        sense = edgeLoopSense[i];

        status = EG_getTessEdge(tess, edgeIndex+1, &nlen, &pxyz, &pt);
        AIM_STATUS(aimInfo, status);

        jbeg = sense == SFORWARD ? 0 : nlen-1;
        jend = sense == SFORWARD ? nlen-1 : 0;
        j = jbeg;
        for (; ipnt < numPoint; ipnt++) { // avoid first point due to rounding errors

            while (j != jend) {
              x1 = (xdot[0]*(pxyz[3*(j      )+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j      )+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j      )+2]-xyzLE[2]))/chord;
              x2 = (xdot[0]*(pxyz[3*(j+sense)+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j+sense)+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j+sense)+2]-xyzLE[2]))/chord;
              if (x1 <= xCoord[ipnt] && xCoord[ipnt] <= x2) {
                t = (pt[j] + pt[j+sense])/2.;
                break;
              }
              j += sense;
            }
            if (j == jend) { break; }

            status = EG_evaluate(edges[edgeIndex], &t, result);
            AIM_STATUS(aimInfo, status);

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            x1 = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
            // use Newton solve to refine the t-value
            status = _refineT(aimInfo, x1, xCoord[ipnt], scale, chord,
                              edges[edgeIndex], xdot, result, xyzLE, &t);
            AIM_STATUS(aimInfo, status);

            yLower[ipnt] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        }
    }

    // Enforce leading and trailing edge nodes and fill in first/last points
    status = EG_evaluate(nodeLE, NULL, result);
    AIM_STATUS(aimInfo, status);

    result[0] -= xyzLE[0];
    result[1] -= xyzLE[1];
    result[2] -= xyzLE[2];

    yUpper[0] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
    yLower[0] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;

    status = EG_evaluate(nodeTE, NULL, result);
    AIM_STATUS(aimInfo, status);

    result[0] -= xyzLE[0];
    result[1] -= xyzLE[1];
    result[2] -= xyzLE[2];

    yUpper[numPoint-1] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
    yLower[numPoint-1] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;

    for (i = 0; i < numPoint; i++) {
        if (fabs(xCoord[i]) < 1e-12) {
            xCoord[i] = 0.0;
        }
        if (fabs(yUpper[i]) < 1e-12) {
            yUpper[i] = 0.0;
        }
        if (fabs(yLower[i]) < 1e-12) {
            yLower[i] = 0.0;
        }
    }

    *xCoordOut = xCoord;
    *yUpperOut = yUpper;
    *yLowerOut = yLower;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
        AIM_FREE(xCoord);
        AIM_FREE(yUpper);
        AIM_FREE(yLower);
    }

    if (edgeSegs != NULL) {
      for (i = 0; i < numEdge; i++)
        AIM_FREE(edgeSegs[i]);
      AIM_FREE(edgeSegs);
    }
    AIM_FREE(numEdgeSegs);
    AIM_FREE(numEdgePoint);

    AIM_FREE(nodes);
    AIM_FREE(edges);

    AIM_FREE(edgeLoopSense);
    AIM_FREE(edgeLoopOrder);

    EG_deleteObject(tess);

    return status;
}

// Get the camber line for a set of x coordinates
int vlm_getSectionCamberLine(void *aimInfo,
                             vlmSectionStruct *vlmSection,
                             double Cspace,       // Chordwise spacing (see spacer)
                             int normalize,       // Normalize by chord (true/false)
                             int numPoint,        // Number of points in airfoil
                             double **xCoordOut,  // [numPoint] increasing x values
                             double **yCamberOut) // [numPoint] camber line y values
{

    int i, status;
    double *xCoord = NULL, *yUpper = NULL, *yLower = NULL, *yCamber = NULL;

    status = vlm_getSectionCoordX(aimInfo,
                                  vlmSection,
                                  Cspace, // Cosine distribution
                                  normalize, (int)false, numPoint,
                                  &xCoord, &yUpper, &yLower);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(yCamber, numPoint, double, aimInfo, status);

    for (i = 0; i < numPoint; i++) {
        yCamber[i] = (yUpper[i] + yLower[i]) / 2;
    }

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS) {
        AIM_FREE(xCoord);
        AIM_FREE(yCamber);
    }
    else {
        *xCoordOut = xCoord;
        *yCamberOut = yCamber;
    }

    AIM_FREE(yUpper);
    AIM_FREE(yLower);

    return status;
}


// Get the airfoil thickness over chord (ToC) given a vlmSectionStruct
int vlm_getSectionToC(void *aimInfo,
                      vlmSectionStruct *vlmSection,
                      double *ToC)    // Resulting Thickness over Chord
{
    int status; // Function return status

    int i, j, ipnt, jbeg, jend; // Indexing

    int nodeIndexLE;
    int edgeIndex, *edgeLoopOrder=NULL, *edgeLoopSense = NULL;

    int *numEdgePoint = NULL, *numEdgeSegs = NULL;;
    CurvatureSpace **edgeSegs = NULL;

    double chord, x1, x2;
    double xdot[3], ydot[3], *secnorm;

    int numEdge, numNode;
    ego *nodes = NULL, *edges = NULL, nodeTE = NULL, nodeLE;

    int sense = 0;
    double t, trange[4], result[18], params[3], scale;

    //EGADS returns
    int oclass, mtype, *sens = NULL, numChildren;
    ego ref, *children = NULL;

    ego teObj = NULL, body, tess=NULL;

    int nlen;
    const double *pxyz, *pt;
    double *xyzLE, *xyzTE;

    double *xCoord=NULL, *yUpper=NULL, *yLower=NULL;

    double thick = 0;
    const int numPoint = 201;
    const double Cspace = 0; // Equal spacing

    body = vlmSection->ebody;
    scale = chord = vlmSection->chord;
    secnorm = vlmSection->normal;
    nodeIndexLE = vlmSection->nodeIndexLE;
    xyzLE = vlmSection->xyzLE;
    xyzTE = vlmSection->xyzTE;
    teObj = vlmSection->teObj;

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);

    // get the leading edge node
    nodeLE = nodes[nodeIndexLE-1];

    // Get the number of points on each edge
    status = vlm_secEdgePoints(aimInfo, numPoint,
                               numEdge, edges,
                               teObj, &numEdgePoint,
                               &numEdgeSegs, &edgeSegs);
    AIM_STATUS(aimInfo, status);

    // Get the loop edge ordering so it starts at the trailing edge NODE
    status = vlm_secOrderEdges(numNode, nodes,
                               numEdge, edges,
                               body, teObj,
                               &edgeLoopOrder, &edgeLoopSense, &nodeTE);
    AIM_STATUS(aimInfo, status);

    // vector from LE to TE normalized
    xdot[0]  =  xyzTE[0] - xyzLE[0];
    xdot[1]  =  xyzTE[1] - xyzLE[1];
    xdot[2]  =  xyzTE[2] - xyzLE[2];

    scale /= chord;
    xdot[0] /=  chord;
    xdot[1] /=  chord;
    xdot[2] /=  chord;

    // cross with section PLANE normal to get perpendicular vector in the PLANE
    cross_DoubleVal(secnorm, xdot, ydot);

    // create a tessellation object on the edges

    // loop over all edges and set the desired point count
    for (i = 0; i < numEdge; i++) {
      status = EG_attributeAdd(edges[i], ".nPos", ATTRINT, 1, &numEdgePoint[i], NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
    }

    // Negating the first parameter triggers EGADS to only put vertexes on edges
    params[0] = -chord;
    params[1] =  chord;
    params[2] =  20;

    status = EG_makeTessBody(body, params, &tess);
    AIM_STATUS(aimInfo, status);

    // set output points

    xCoord = (double *) EG_alloc(numPoint*sizeof(double));
    yUpper = (double *) EG_alloc(numPoint*sizeof(double));
    yLower = (double *) EG_alloc(numPoint*sizeof(double));
    if (xCoord == NULL || yUpper == NULL || yLower == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // generate the x-coordinates using AVL algorithm
    spacer( numPoint, Cspace, scale, xCoord);

    // Loop over the upper surface
    ipnt = numPoint-2; // avoid last point due to rounding errors
    for (i = 0; i < numEdge; i++) {
        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        if (edges[edgeIndex] == teObj) continue;

        // Get t- range for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        AIM_STATUS(aimInfo, status);

        // Check if the first node is the LE edge, which means the current edge is on the lower surface
        sense = edgeLoopSense[i];
        if (sense == SFORWARD) {
            if (children[0] == nodeLE) break;
        } else {
            if (children[1] == nodeLE) break;
        }

        status = EG_getTessEdge(tess, edgeIndex+1, &nlen, &pxyz, &pt);
        AIM_STATUS(aimInfo, status);

        jbeg = sense == SFORWARD ? 0 : nlen-1;
        jend = sense == SFORWARD ? nlen-1 : 0;
        j = jbeg;
        for (; ipnt > 0; ipnt--) { // avoid first point due to rounding errors

            while (j != jend) {
              x2 = (xdot[0]*(pxyz[3*(j      )+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j      )+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j      )+2]-xyzLE[2]))/chord;
              x1 = (xdot[0]*(pxyz[3*(j+sense)+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j+sense)+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j+sense)+2]-xyzLE[2]))/chord;
              if (x1 <= xCoord[ipnt] && xCoord[ipnt] <= x2) {
                t = (pt[j] + pt[j+sense])/2.;
                break;
              }
              j += sense;
            }
            if (j == jend) { break; }

            status = EG_evaluate(edges[edgeIndex], &t, result);
            AIM_STATUS(aimInfo, status);

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            x1 = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
            // use Newton solve to refine the t-value
            status = _refineT(aimInfo, x1, xCoord[ipnt], scale, chord,
                              edges[edgeIndex], xdot, result, xyzLE, &t);
            AIM_STATUS(aimInfo, status);

            yUpper[ipnt] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        }
    }

    // Loop over the lower surface
    ipnt = 1;  // avoid first point due to rounding errors
    for (; i < numEdge; i++) {
        edgeIndex = edgeLoopOrder[i] - 1; // -1 indexing

        if (edges[edgeIndex] == teObj) continue;

        // Get t- range for edge
        status = EG_getTopology(edges[edgeIndex], &ref, &oclass, &mtype, trange, &numChildren, &children, &sens);
        if (mtype == DEGENERATE) continue;
        AIM_STATUS(aimInfo, status);

        sense = edgeLoopSense[i];

        status = EG_getTessEdge(tess, edgeIndex+1, &nlen, &pxyz, &pt);
        AIM_STATUS(aimInfo, status);

        jbeg = sense == SFORWARD ? 0 : nlen-1;
        jend = sense == SFORWARD ? nlen-1 : 0;
        j = jbeg;
        for (; ipnt < numPoint; ipnt++) { // avoid first point due to rounding errors

            while (j != jend) {
              x1 = (xdot[0]*(pxyz[3*(j      )+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j      )+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j      )+2]-xyzLE[2]))/chord;
              x2 = (xdot[0]*(pxyz[3*(j+sense)+0]-xyzLE[0])+
                    xdot[1]*(pxyz[3*(j+sense)+1]-xyzLE[1])+
                    xdot[2]*(pxyz[3*(j+sense)+2]-xyzLE[2]))/chord;
              if (x1 <= xCoord[ipnt] && xCoord[ipnt] <= x2) {
                t = (pt[j] + pt[j+sense])/2.;
                break;
              }
              j += sense;
            }
            if (j == jend) { break; }

            status = EG_evaluate(edges[edgeIndex], &t, result);
            AIM_STATUS(aimInfo, status);

            result[0] -= xyzLE[0];
            result[1] -= xyzLE[1];
            result[2] -= xyzLE[2];

            x1 = (xdot[0]*result[0]+xdot[1]*result[1]+xdot[2]*result[2])/chord;
            // use Newton solve to refine the t-value
            status = _refineT(aimInfo, x1, xCoord[ipnt], scale, chord,
                              edges[edgeIndex], xdot, result, xyzLE, &t);
            AIM_STATUS(aimInfo, status);

            yLower[ipnt] = (ydot[0]*result[0]+ydot[1]*result[1]+ydot[2]*result[2])/chord;
        }
    }

    // Find the approximate max ToC, skip LE and TE
    (*ToC) = 0;
    for( i = 1; i < numPoint-1; i++) {
        thick = yUpper[i] - yLower[i];
        if (thick > (*ToC)) {
            (*ToC) = thick;
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(xCoord);
    AIM_FREE(yUpper);
    AIM_FREE(yLower);

    if (edgeSegs != NULL) {
      for (i = 0; i < numEdge; i++)
        AIM_FREE(edgeSegs[i]);
      AIM_FREE(edgeSegs);
    }
    AIM_FREE(numEdgeSegs);
    AIM_FREE(numEdgePoint);

    AIM_FREE(nodes);
    AIM_FREE(edges);

    AIM_FREE(edgeLoopSense);
    AIM_FREE(edgeLoopOrder);

    EG_deleteObject(tess);

    return status;
}
