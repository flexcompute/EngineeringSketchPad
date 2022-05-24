// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include "capsTypes.h" // Bring in CAPS types
#include "vlmTypes.h"  // Bring in Vortex Lattice Method structures
#include "miscTypes.h"  // Bring in miscellaneous types

#ifdef __cplusplus
extern "C" {
#endif

// Fill vlmSurface in a vlmSurfaceStruct format with vortex lattice information
// from an incoming surfaceTuple
int get_vlmSurface(int numTuple,
                   capsTuple surfaceTuple[],
                   mapAttrToIndexStruct *attrMap,
                   double Cspace,
                   int *numVLMSurface,
                   vlmSurfaceStruct *vlmSurface[]);

// Fill vlmControl in a vlmControlStruct format with vortex lattice information
// from an incoming controlTuple
int get_vlmControl(void *aimInfo,
                   int numTuple,
                   capsTuple controlTuple[],
                   int *numVLMControl,
                   vlmControlStruct *vlmControl[]);


// Initiate (0 out all values and NULL all pointers) of a control in the vlmcontrol structure format
int initiate_vlmControlStruct(vlmControlStruct *control);

// Destroy (0 out all values and NULL all pointers) of a control in the vlmcontrol structure format
int destroy_vlmControlStruct(vlmControlStruct *control);

// Initiate (0 out all values and NULL all pointers) of a section in the vlmSection structure format
int initiate_vlmSectionStruct(vlmSectionStruct *section);

// Destroy (0 out all values and NULL all pointers) of a section in the vlmSection structure format
int destroy_vlmSectionStruct(vlmSectionStruct *section);

// Initiate (0 out all values and NULL all pointers) of a surface in the vlmSurface structure format
int initiate_vlmSurfaceStruct(vlmSurfaceStruct *surface);

// Destroy (0 out all values and NULL all pointers) of a surface in the vlmSurface structure format
int destroy_vlmSurfaceStruct(vlmSurfaceStruct  *surface);

// Retrieve the string following a vlmControl tag
//int retrieve_vlmControlAttr(ego geomEntity, const char **attrName, const char **string);
//int retrieve_vlmControlAttr(ego geomEntity, const char **string);

// Populate vlmSurface-section control surfaces from geometry attributes, modify control properties based on
// incoming vlmControl structures
int get_ControlSurface(void *aimInfo,
                       int numControl,
                       vlmControlStruct vlmControl[],
                       vlmSurfaceStruct *vlmSurface);

// Make a copy of vlmControlStruct (it is assumed controlOut has already been initialized)
int copy_vlmControlStruct(vlmControlStruct *controlIn, vlmControlStruct *controlOut);

// Make a copy of vlmSectionStruct (it is assumed sectionOut has already been initialized)
int copy_vlmSectionStruct(vlmSectionStruct *sectionIn, vlmSectionStruct *sectionOut);

// Make a copy of vlmSurfaceStruct (it is assumed surfaceOut has already been initialized)
// Also the section in vlmSurface are reordered based on a vlm_orderSections() function call
int copy_vlmSurfaceStruct(vlmSurfaceStruct *surfaceIn, vlmSurfaceStruct *surfaceOut);

// Finalizes populating vlmSectionStruct member data after the ebody is set
int finalize_vlmSectionStruct(void *aimInfo, vlmSectionStruct *vlmSection);

// Accumulate VLM section data from a set of bodies.If disciplineFilter is not NULL
// bodies not found with disciplineFilter (case insensitive) for a capsDiscipline attribute
// will be ignored.
int vlm_getSections(void *aimInfo,
                    int numBody,
                    ego bodies[],
                    /*@null@*/ const char *disciplineFilter,
                    mapAttrToIndexStruct attrMap,
                    vlmSystemEnum sys,
                    int numSurface,
                    vlmSurfaceStruct *vlmSurface[]);

// Order VLM sections increasing order
int vlm_orderSections(int numSection, vlmSectionStruct section[]);

// Compute spanwise panel spacing with close to equal spacing on each pane
int vlm_equalSpaceSpanPanels(void *aimInfo, int NspanTotal, int numSection, vlmSectionStruct vlmSection[]);

// Get the airfoil cross-section tessellation ego given a vlmSectionStruct
int vlm_getSectionTessSens(void *aimInfo,
                           vlmSectionStruct *vlmSection,
                           int normalize,      // Normalize by chord (true/false)
                           const char *geomInName,
                           const int irow, const int icol,
                           int numPoint,    // Number of points in airfoil
                           double **dx_dvar_out,
                           double **dy_dvar_out);

// Get the airfoil cross-section given a vlmSectionStruct
int vlm_getSectionCoord(void *aimInfo,
                        vlmSectionStruct *vlmSection,
                        int normalize,       // Normalize by chord (true/false)
                        int numPoint,        // Number of points in airfoil
                        double **xCoordOut,  //[maxNumPoint]
                        double **yCoordOut); //[maxNumPoint] for upper and lower surface

// Write out the airfoil cross-section given a vlmSectionStruct
int vlm_writeSection(void *aimInfo,
                     FILE *fp,
                     vlmSectionStruct *vlmSection,
                     int normalize, // Normalize by chord (true/false)
                     int numPoint); // Number of points in airfoil

// Get the airfoil cross-section given a vlmSectionStruct
// where y-upper and y-lower correspond to the x-value
// Only works for sharp trailing edges
int vlm_getSectionCoordX(vlmSectionStruct *vlmSection,
                         double Cspace,       // Chordwise spacing (see spacer)
                         int normalize,       // Normalize by chord (true/false)
                         int rotated,         // Leave airfoil rotated (true/false)
                         int numPoint,        // Number of points in airfoil
                         double **xCoordOut,  // [numPoint] increasing x values
                         double **yUpperOut,  // [numPoint] for upper surface
                         double **yLowerOut); // [numPoint] for lower surface

// Get the camber line for a set of x coordinates
int vlm_getSectionCamberLine(vlmSectionStruct *vlmSection,
                            double Cspace,       // Chordwise spacing (see spacer)
                            int normalize,       // Normalize by chord (true/false)
                            int numPoint,        // Number of points in airfoil
                            double **xCoordOut,  // [numPoint] increasing x values
                            double **yCamberOut);// [numPoint] camber line y values


#ifdef __cplusplus
}
#endif
