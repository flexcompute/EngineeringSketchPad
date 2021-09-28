// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Structures for vortex lattice analysis - Written by Dr. Ryan Durscher AFRL/RQVC

#ifndef VLMTYPES_H
#define VLMTYPES_H


typedef enum  {vlmGENERIC, vlmPLANEYZ, vlmRADIAL} vlmSystemEnum;

typedef struct {

    char *name; // Control surface name

    double deflectionAngle; // Deflection angle of the control surface
    double controlGain; //Control deflection gain, units:  degrees deflection / control variable

    double percentChord; // Percentage along chord
    double xyzHinge[3]; // xyz location of the hinge

    double xyzHingeVec[3]; // Vector of hinge line at xyzHinge

    int    leOrTeOverride; // Does the user want to override the geometry set value?
    int    leOrTe; // Leading = 0 or trailing > 0 edge control surface

    int  deflectionDup; // Sign of deflection for duplicated surface
} vlmControlStruct;

typedef struct {
    char *name; // Section Name

    ego ebody;       // Body of the section (might be flipped relative to original)

    int sectionIndex; // Section index - 0 bias

    double xyzLE[3]; // xyz coordinates for leading edge
    int nodeIndexLE; // Leading edge node (in geometry) index with reference to xyzLE - 1 bias

    double xyzTE[3]; // xyz coordinates for trailing edge (Node or Edge mid point)
    ego teObj;       // Trailing edge object in the body
    int teClass;     // Trailing edge object class (NODE or EDGE)

    double chord;    // section chord length
    double ainc;     // section incidence angle

    double normal[3]; // planar normal for the section

    int    Nspan;   // number of spanwise vortices (elements)
    double Sspace;  // spanwise point distribution
    int    Sset;    // spanwise point distribution set by user

    int numControl;
    vlmControlStruct *vlmControl;

} vlmSectionStruct;

typedef struct {
    char   *name; // Name of surface

    int numAttr;    // Number of capsGroup/attributes used to define a given surface
    int *attrIndex; // List of attribute map integers that correspond given capsGroups

    double Cspace;
    double Sspace;

    int    Nchord;
    int    NspanTotal;    // Total number of spanwise vortices (elements) on the surface
    int    NspanSection;  // Number of spanwise vortices on each section of the surface

    int    nowake;
    int    noalbe;
    int    noload;
    int    compon;
    int    iYdup;

    int numSection; // Number of sections that make up the surface
    vlmSectionStruct *vlmSection; // Section information

    char *surfaceType;

} vlmSurfaceStruct;

#endif
