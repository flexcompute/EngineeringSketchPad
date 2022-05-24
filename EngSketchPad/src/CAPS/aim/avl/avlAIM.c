/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AVL AIM
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 *      Modified by Ryan Durscher AFLR/RQVC May 10, 2016
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "caps.h"

#include "miscUtils.h" // Bring in miscellaneous utilities
#include "vlmUtils.h"
#include "vlmSpanSpace.h"
#include "cfdUtils.h"

#ifdef WIN32
#define strcasecmp  stricmp
#define strncasecmp _strnicmp
#define strtok_r    strtok_s
#endif

#define NUMPOINT  200
#define PI        3.1415926535897931159979635
#define NINT(A)         (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))

//#define DEBUG

enum aimInputs
{
  inMach = 1,                    /* index is 1-based */
  inAlpha,
  inBeta,
  inRollRate,
  inPitchRate,
  inYawRate,
  inCDp,
  inAVL_Surface,
  inAVL_Control,
  inCL,
  inMoment_Center,
  inMassProp,
  inGravity,
  inDensity,
  inVelocity,
  NUMINPUT = inVelocity        /* Total number of inputs */
};

enum aimOutputs
{
  Alpha = 1,                   /* index is 1-based */
  Beta,
  Mach,
  pbd2V,
  qcd2V,
  rbd2V,
  pPbd2V,
  rPbd2V,
  CXtot,
  CYtot,
  CZtot,
  Cltot,
  Cmtot,
  Cntot,
  ClPtot,
  CnPtot,
  CLtot,
  CDtot,
  CDvis,
  CLff,
  CYff,
  CDind,
  CDff,
  e,
  CLa,
  CYa,
  ClPa,
  Cma,
  CnPa,
  CLb,
  CYb,
  ClPb,
  Cmb,
  CnPb,
  CLpP,
  CYpP,
  ClPpP,
  CmpP,
  CnPpP,
  CLqP,
  CYqP,
  ClPqP,
  CmqP,
  CnPqP,
  CLrP,
  CYrP,
  ClPrP,
  CmrP,
  CnPrP,
  CXu,
  CYu,
  CZu,
  Clu,
  Cmu,
  Cnu,
  CXv,
  CYv,
  CZv,
  Clv,
  Cmv,
  Cnv,
  CXw,
  CYw,
  CZw,
  Clw,
  Cmw,
  Cnw,
  CXp,
  CYp,
  CZp,
  Clp,
  Cmp,
  Cnp,
  CXq,
  CYq,
  CZq,
  Clq,
  Cmq,
  Cnq,
  CXr,
  CYr,
  CZr,
  Clr,
  Cmr,
  Cnr,
  Xnp,
  Xcg,
  Ycg,
  Zcg,
  ControlStability,
  ControlBody,
  HingeMoment,
  StripForces,
  EigenValues,
  NUMOUT = EigenValues         /* Total number of outputs */
};

/*!\mainpage Introduction
 * \tableofcontents
 *
 * \section overviewAVL AVL AIM Overview
 * The use of lower-dimensional design tools is clearly desirable in a multidisciplinary/multi-fidelity aero
 * design optimization setting. This is the crux of the Computational Aircraft Prototype Syntheses (CAPS) program.
 * In many ways describing geometry appropriate for AVL (the Athena Vortex Lattice) code is more cumbersome than
 * higher fidelity codes that require an Outer Mold Line. The goal is to make a CAPS AIM (Analysis Input Module)
 * that directly feeds input to AVL and extracts the output quantities of interest from AVL's execution. This
 * needs to be consistent with a build description that is hierarchical and multi-fidelity. That is, the build
 * description that generates the geometric data at this level can be further enhanced to produce the complete
 * OML of the aircraft design under consideration. As for the geometric description, AVL requires airfoil section
 * data specified at the appropriate locations that describe the <em> skeleton </em> of the aircraft. These
 * sections when <em> lofted </em> as groups and finally <em> unioned </em> together builds the OML. Clearly,
 * intercepting the state of the geometry before these higher-level operations are applied provides the data
 * appropriate for AVL. This naturally constructs a hierarchical geometric view where a design can progress
 * into higher fidelities and feedback can be achieved where we can go back to this level of description when need be.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsAVL and \ref aimOutputsAVL, respectively. An
 * alternative to the AIM's outputs for retrieving sensitivity information is provided in \ref aimBackDoorAVL.
 *
 * Details on the use of units are outlined in \ref aimUnitsAVL.
 *
 * Geometric attribution that the AIM makes use is provided in \ref attributeAVL.
 *
 * The AVL AIM can automatically execute avl, with details provided in \ref aimExecuteAVL.
 *
 * To populate output data the AIM expects files, "capsTotalForce.txt", "capsStripForce.txt", "capsStatbilityDeriv.txt", "capsBodyAxisDeriv.txt", and
 * "capsHingeMoment.txt" to exist after running AVL
 * (see \ref aimOutputsAVL for additional information). An example execution for AVL looks like:
 *
 *
 * \section assumptionsAVL Assumptions
 * The AVL coordinate system assumption (X -- downstream, Y -- out the right wing, Z -- up) needs to be followed.
 *
 * Within <b> OpenCSM</b> there are a number of airfoil generation UDPs (User Defined Primitives). These include NACA 4
 * series, a more general NACA 4/5/6 series generator, Sobieczky's PARSEC parameterization and Kulfan's CST
 * parameterization. All of these UDPs generate <b> EGADS</b> <em> FaceBodies</em> where the <em>Face</em>'s underlying
 * <em>Surface</em>  is planar and the bounds of the <em>Face</em> is a closed set of <em>Edges</em> whose
 * underlying <em>Curves</em> contain the airfoil shape. In all cases, there is a <em>Node</em> that represents
 * the <em>Leading Edge</em> point and one or two <em>Nodes</em> at the <em>Trailing Edge</em> -- one if the
 * representation is for a sharp TE and the other if the definition is open or blunt. If there are 2 <em>Nodes</em>
 * at the back, then there are 3 <em>Edges</em> all together and closed, even though the airfoil definition
 * was left open at the TE. All of this information will be used to automatically fill in the AVL geometry description.
 *
 * The AVL Sections are automatically generated, one from each <em> FaceBody</em> and the details extracted from the geometry.
 * The <em> FaceBody</em> must contain at least two edges and two nodes, but may contain any number of <em> Edges</em> otherwise.
 * If the <em> FaceBody</em> contains more nodes, the node with the smallest <b> x</b> value is used to
 * define the leading edge, the node with the largest <b> x</b> defines the trailing edge. The airfoil may have a
 * single <em> Edge</em> that defines a straight blunt trailing edge. <b> Xle</b>, <b> Yle</b>, and <b> Zle</b>, are taken
 * from the <em> Node</em> associated with the <em> Leading Edge</em>. The <b> Chord</b> is computed by getting
 * the distance between the LE and TE (if there is a blunt trainling <em> Edge</em> in the <em> FaceBody</em>
 * the TE point is considered the mid-position on that <em> Edge</em>). <b> Ainc</b> is computed by
 * registering the chordal direction of the <em> FaceBody</em> against the X-Z plane. The airfoil shapes are generated by
 * sampling the <em> Curves</em> and put directly in the input file via the <b> AIRFOIL</b> keyword after being normalized.
 *
 * It should be noted that general construction in either <b> OpenCSM</b> or even <b> EGADS</b> will be supported
 * as long as the topology described above is used. But care should be taken when constructing the airfoil shape
 * so that a discontinuity (i.e.,  simply <em>C<sup>0</sup></em>) is not generated at the <em>Node</em> representing
 * the <em>Leading Edge</em>. This can be done by splining the entire shape as one and then intersecting the single
 *  <em>Edge</em> to place the LE <em>Node</em>.
 *
 *  The rest of the information and options required to fill out the AVL geometry input file (<b> xxx.avl</b>)
 *  will be found in the attributes attached to the <em>FaceBody</em> itself. The conventions used will be
 *  described in the next section.
 *
 * Also note that this first implementation is not intended to provide complete control over AVL. In particular,
 * there is no mention above of the <b> BODY</b>, <b> DESIGN</b>, <b> CLAF</b>, or <b> CDCL</b> AVL keywords.
 *
 * \section exampleAVL Examples
 * An example problem using the AVL AIM may be found at \ref avlExample, which
 * contains example *.csm input files and pyCAPS scripts designed to make use of the AVL AIM.  These example
 * scripts make extensive use of the \ref attributeAVL, \ref aimInputsAVL, and \ref aimOutputsAVL.
 *
 */

/*! \page avlExample AVL AIM Examples
 * This example contains a set of *.csm and pyCAPS (*.py) inputs that uses the AVL AIM.  A user should have knowledge
 * on the generation of parametric geometry in Engineering Sketch Pad (ESP) before attempting to integrate with any AIM.
 * Specifically, this example makes use of Design Parameters, Set Parameters, User Defined Primitive (UDP) and attributes
 * in ESP.
 *
 * The follow code details the process in a *.csm file that generates three airfoil sections to create a wing. Note
 * to execute in serveCSM a dictionary file must be included
 * "serveCSM $ESP_ROOT/CAPSexamples/csmData/avlWing.csm"
 *
 * The CSM script generates Bodies which are designed to be used by specific AIMs.
 * The AIMs that the Body is designed for is communicated to the CAPS framework via
 * the capsAIM string attribute. This is a semicolon-separated string with the list of
 * AIM names. Thus, the CSM author can give a clear indication to which AIMs should
 * use the Body. In this example, the list contains only the avlAIM:
 * \snippet avlWing.csm capsAIM
 *
 * Next we will define the design parameters to define the wing cross section and planform.
 * \snippet avlWing.csm designParameters
 *
 * The design parameters will then be used to set parameters for use internally to create geometry.
 * \snippet avlWing.csm setParameters
 *
 * Finally, the airfoils are created using the User Defined Primitive (UDP) naca. The inputs used for this
 * example to the UDP are Thickness and Camber. Cross sections are in the X-Y plane and are rotated to the X-Z plane.
 * Reference quantities must exist on any body, otherwise AVL defaults to 1.0 for Area, Span,
 * Chord and 0.0 for X,Y,Z moment References
 *
 * \snippet avlWing.csm createAirfoils
 *
 * An example pyCAPS script that uses the above csm file to run AVL is as follows.
 *
 * First the pyCAPS and os module needs to be imported.
 * \snippet avl_PyTest.py import
 *
 * Next the *.csm file is loaded and design parameter is changed - area in the geometry. Any despmtr from the
 * avlWing.csm file is available inside the pyCAPS script. They are: thick, camber, area, aspect, taper,
 * sweep, washout, dihedral...
 *
 * \snippet avl_PyTest.py geometry
 *
 * The AVL AIM is then loaded with:
 * \snippet avl_PyTest.py loadAIM
 *
 * After the AIM is loaded the Mach number and angle of attack are set, though all \ref aimInputsAVL are available.
 * \snippet avl_PyTest.py setInputs
 *
 * Once all the inputs have been set, outputs can be directly requested. The avl
 * analysis will be automatically executed just-in-time (\ref aimExecuteAvl).
 *
 * Any of the AIM's output variables (\ref aimOutputsAVL) are readily available; for example,
 * \snippet avl_PyTest.py output
 *
 * results in
 * \code
 * CXtot   0.00061
 * CYtot   -0.0
 * CZtot   -0.30129
 * Cltot   -0.0
 * Cmtot   -0.19449
 * Cntot   -0.0
 * Cl'tot  -0.0
 * Cn'tot  -0.0
 * CLtot   0.30126
 * CDtot   0.00465
 * CDvis   0.0
 * CLff    0.30096
 * CYff    -0.0
 * CDind   0.0046467
 * CDff    0.0049692
 * e       0.967
 * \endcode
 *
 * Additionally, besides making a call to the AIM outputs, sensitivity values may be obtained in the following manner,
 * \snippet avl_PyTest.py sensitivity
 *
 * The avlAIM supports the control surface modeling functionality inside AVL.  Trailing edge control surfaces can be
 * added to the above example by making use of the vlmControlName attribute (see \ref attributeAVL regarding the attribution
 * specifics).  To add a <b>RightFlap</b> and <b>LeftFlap</b> to the
 * previous example *.csm file the naca UDP entries are augmented with the following attributes.
 *
 * \code {.py}
 * # left tip
 * udprim    naca      Thickness thick     Camber    camber
 * attribute vlmControl_LeftFlap 80 # Hinge line is at 80% of the chord
 * ...
 *
 * # root
 * udprim    naca      Thickness thick     Camber    camber
 * attribute vlmControl_LeftFlap 80 # Hinge line is at 80% of the chord
 * attribute vlmControl_RightFlap 80 # Hinge line is at 80% of the chord
 * ...
 *
 * # right tip
 * udprim    naca      Thickness thick     Camber    camber
 * attribute vlmControl_RightFlap 80 # Hinge line is at 80% of the chord
 * ...
 * \endcode
 *
 * Note how the root airfoil contains two attributes for both the left and right flaps.
 *
 * In the pyCAPS script the \ref aimInputsAVL, <b>AVL_Control</b>, must be defined.
 * \code {.py}
 * flap = {"controlGain"     : 0.5,
 *         "deflectionAngle" : 10.0}
 *
 *  myAnalysis.input.AVL_Control = {"LeftFlap": flap, "RightFlap": flap}
 \endcode
 *
 * Notice how the information defined in the <b>flap</b> variable is assigned to the vlmControl<b>Name</b> portion of the
 * attributes added to the *.csm file.
 *
 */

/*! \page attributeAVL AIM Attributes
 * The following list of attributes drives the AVL geometric definition. Each <em>FaceBody</em> which relates to AVL
 * <b> Sections</b> will be marked up in an appropriate manner to drive the input file construction. Many attributes
 * are required and those that are optional are marked so in the description:
 *
 * - <b> capsReferenceArea</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the SREF entry in the AVL input.
 *
 * - <b> capsReferenceChord</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the CREF entry in the AVL input.
 *
 * - <b> capsReferenceSpan</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the BREF entry in the AVL input.
 *
 * - <b> capsReferenceX</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Xref entry in the AVL input.
 *
 * - <b> capsReferenceY</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Yref entry in the AVL input.
 *
 * - <b> capsReferenceZ</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 *  value will be used as the Zref entry in the AVL input.
 *
 * - <b> capsGroup</b> This string attribute labels the <em> FaceBody</em> as to which AVL Surface the section
 *  is assigned. This should be something like: <em> Main\_Wing</em>, <em> Horizontal\_Tail</em>, etc. This informs
 *  the AVL AIM to collect all <em> FaceBodies</em> that match this attribute into a single AVL Surface.
 *
 * - <b> vlmControl"Name"</b> This string attribute attaches a control surface to the <em> FaceBody</em>. The hinge
 * location is defined as the double value between 0 or 1.0.  The range as percentage from 0 to 100
 * will also work. The name of the
 * control surface is the string information after vlmControl (or vlmControl\_).  For Example, to define a control
 * surface named Aileron the following are identical (<em>attribute vlmControlAileron 0.8 </em> or
 * <em>attribute vlmControl\_Aileron 80</em>) . Multiple <em> vlmControl</em> attributes, with different
 * names, can be defined on a single <em> FaceBody</em>.
 * <br> <br>By default control surfaces with percentages less than 0.5 (< 50%) are considered leading edge flaps, while values
 * greater than or equal to 0.5 (>= 50%) are considered trailing edge flaps. This behavior may be overwritten when setting up
 * the control surface in "AVL_Control" (see \ref aimInputsAVL) with the keyword "leOrTe" (see \ref vlmControl for additional details).
 *
 * - <b> vlmNumSpan </b> This attribute may be set on any given airfoil cross-section to overwrite the number of
 * spanwise horseshoe vortices placed on the surface (globally set - see keyword "numSpanPerSection" and "numSpanTotal" in \ref vlmSurface )
 * between two sections. Note, that the AIM internally sorts the sections in ascending y (or z) order, so care
 * should be taken to select the correct section for the desired intent.
 *
 * Note: The attribute <b> avlNumSpan </b> has been deprecated in favor of <b> vlmNumSpan </b>
 *
 * - <b> vlmSspace </b> This attribute may be set on any given airfoil cross-section in the range [-1 .. 1] specify the spanwise distribution function.
 *
 *
 */


// To be changed
/*  - <b> avlComponent</em>  [optional]. This integer (or real) attribute is examined for the AVL Component index.
 *  This may be on any or all <em> FaceBodies</em> associated with the AVL Surface.

- <b> avlKeyword</em>  [optional]. This string attribute may contain the word (or words): {\tt NOWAKE</em>, {\tt NOALBE</em> and/or {\tt NOLOAD</em>. Multiple keywords need
to be delimited by a colon (`:'). This is associated with the AVL Surface.


\item {\bf avlSecValues.} This real valued attribute must be of at least 5 {\it doubles} in length and contain:
   \begin{enumerate}
   \item {\tt Nchord} -- The number of chordwise horseshoe vortices placed on the surface.
   \item {\tt Cspace} -- The chordwise vortex spacing parameter.
   \item {\tt Nspan} -- The number of spanwise horseshoe vortices placed on the surface.
   \item {\tt Sspace} -- The spanwise vortex spacing parameter.
   \item {\tt iYdup} -- The Y symmetry for this Surface (0 -- do nothing, 1 -- duplicate about the Y=0.0 plane).
   \end{enumerate}
   It should be noted that the first 2 (and last) values refer to the entire AVL Surface/capsGroup and therefore the first occurrence of the \textbf{\textit {capsGroup}} attribute in the list of {\it bodies} sets these values in the AVL input file (though {\tt iYdup} can be set in any Section).

\item {\bf avlCntrl}\textbf{\textit {Name}}  [optional]. This optional numerical attribute indicates that this {\it FaceBody} Section is a member of the control surface named as part of the attribute name (for example:  {\tt avlCntrlAileron} specifies that this Section is part of the {\it Aileron} control surface). This real valued attribute must be of at least 6 {\it doubles} in length and contain:
   \begin{enumerate}
   \item {\tt gain} -- The control deflection gain.
   \item {\tt Xhinge} -- The x/c location of hinge.
   \item {\tt Xhvec} -- The first component of the vector giving hinge axis about which surface rotates.
   \item {\tt Yhvec} -- The second component.
   \item {\tt Zhvec} -- The third component.
   \item {\tt SgnDup } -- The signed magnitude (or simply the scale factor) for the deflection of the duplicated surface.
   \end{enumerate}
\end{itemize}
 */


typedef struct {

    mapAttrToIndexStruct controlMap; // vlmControl* attribute to index map

    mapAttrToIndexStruct groupMap; // capsGroup attribute to index map

    // AVL surface information
    int numAVLSurface;
    vlmSurfaceStruct *avlSurface;

    // AVL control surface information
    int numAVLControl;
    vlmControlStruct *avlControl;

    // Units structure
    cfdUnitsStruct units;

} aimStorage;


static int initiate_aimStorage(void *aimInfo, aimStorage *avlInstance)
{

    int status;

    // Set initial values for avlInstance
    status = initiate_mapAttrToIndexStruct(&avlInstance->controlMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&avlInstance->groupMap);
    AIM_STATUS(aimInfo, status);

    // AVL surface information
    avlInstance->numAVLSurface = 0;
    avlInstance->avlSurface = NULL;

    avlInstance->numAVLControl = 0;
    avlInstance->avlControl = NULL;

    status = initiate_cfdUnitsStruct(&avlInstance->units);
    AIM_STATUS(aimInfo, status);

cleanup:
    return status;
}


static int destroy_aimStorage(aimStorage *avlInstance)
{
    int status;
    int i;

    status = destroy_mapAttrToIndexStruct(&avlInstance->controlMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    status = destroy_mapAttrToIndexStruct(&avlInstance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    if (avlInstance->avlSurface != NULL) {
        for (i = 0; i < avlInstance->numAVLSurface; i++) {
            status = destroy_vlmSurfaceStruct(&avlInstance->avlSurface[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_vlmSurfaceStruct!\n", status);
        }
        AIM_FREE(avlInstance->avlSurface);
        avlInstance->numAVLSurface = 0;
    }

    if (avlInstance->avlControl != NULL) {
        for (i = 0; i < avlInstance->numAVLControl; i++) {
            status = destroy_vlmControlStruct(&avlInstance->avlControl[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_vlmControlStruct!\n", status);
        }
        AIM_FREE(avlInstance->avlControl);
        avlInstance->numAVLControl = 0;
    }

    status = destroy_cfdUnitsStruct(&avlInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_cfdUnitsStruct!\n", status);

    return CAPS_SUCCESS;
}


/* ********************** AVL AIM Helper Functions ************************** */
static int writeSection(void *aimInfo, FILE *fp, vlmSectionStruct *vlmSection)
{
    int status; // Function return status

    ego body;
    double *xyzLE;
    int Nspan;
    double chord, ainc, Sspace;

    // Attribute variables
    int          atype, alen;
    const int    *ints;
    const char   *string;
    const double *reals;

    // Get data from the section
    body  = vlmSection->ebody;
    xyzLE = vlmSection->xyzLE;
    chord = vlmSection->chord;
    ainc  = vlmSection->ainc;
    Nspan = vlmSection->Nspan;
    Sspace = vlmSection->Sspace;

    status = EG_attributeRet(body, "avlNumSpan", &atype, &alen, &ints, &reals,
                             &string);
    if (status == EGADS_SUCCESS) {
        printf("*************************************************************\n");
        printf("Warning: avlNumSpan is DEPRICATED in favor of vlmNumSpan!!!\n");
        printf("         Please update the attribution.\n");
        printf("*************************************************************\n");

        if (atype != ATTRINT && atype != ATTRREAL && alen != 1) {
            printf ("Error: Attribute %s should be followed by a single integer\n",
                    "avlNumSpan");
        }

        if (atype == ATTRINT) {
            Nspan = ints[0];
        }

        if (atype == ATTRREAL) {
            Nspan = (int) reals[0];
        }
    }

    fprintf(fp, "#Xle     Yls       Zle       Chord    Ainc  Nspan  Sspace\n");
    fprintf(fp, "SECTION\n%lf %lf %lf  %lf %lf  %d %lf\n\n",
            xyzLE[0], xyzLE[1], xyzLE[2], chord, ainc, Nspan, Sspace);
    fprintf(fp, "AIRFOIL\n");

    status = vlm_writeSection(aimInfo,
                              fp,
                              vlmSection,
                              (int) true, // Normalize by chord (true/false)
                              (int) NUMPOINT);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in writeSection, status = %d\n", status);

    return status;
}


static int
writeMassFile(void *aimInfo, capsValue *aimInputs, cfdUnitsStruct *units,
              const char *bodyLunits, char massFilename[])
{
  int status; // Function return status

  const char *keyWord = NULL;
        char *keyValue = NULL; // Freeable

  int numString = 0;
  char **stringArray = NULL; // Freeable

  // Eigen value quantities
  double Lunit, Munit=1.0, Tunit=1.0;
  const char *Lunits, *Munits, *Tunits;
  double gravity, density;
  char *Iunits = NULL, *tmpUnits = NULL;
  double MLL;

  double mass, xyz[3], inertia[6]; // Inertia order = Ixx, Iyy, Izz, Ixy, Ixz, Iyz
  double *I= NULL; // Freeable

  char *value = NULL, *errName = NULL, *errValue = NULL;
  const char *errMsg = NULL;
  capsTuple *massProp;
  int massPropLen, inertiaLen;

  int i, j; // Indexing

  FILE *fp = NULL;

  /* open the file and write the mass data */
  fp = aim_fopen(aimInfo, massFilename, "w");
  if (fp == NULL) {
      printf("Unable to open file %s\n!", massFilename);

      status =  CAPS_IOERR;
      goto cleanup;
  }

  printf("Writing mass properties file: %s\n", massFilename);

  fprintf(fp, "#-------------------------------------------------\n");
  fprintf(fp, "#  Dimensional unit and parameter data.\n");
  fprintf(fp, "#  Mass & Inertia breakdown.\n");
  fprintf(fp, "#-------------------------------------------------\n");
  fprintf(fp, "\n");
  fprintf(fp, "#  Names and scalings for units to be used for trim and eigenmode calculations.\n");
  fprintf(fp, "#  The Lunit and Munit values scale the mass, xyz, and inertia table data below.\n");
  fprintf(fp, "#  Lunit value will also scale all lengths and areas in the AVL input file.\n");

  if (units->length != NULL)
    Lunits = units->length;
  else
    Lunits = "m";

  if (units->mass != NULL)
    Munits = units->mass;
  else
    Munits = "kg";

  if (units->time != NULL)
    Tunits = units->time;
  else
    Tunits = "s";

  // conversion of the csm model units into units of Lunits
  Lunit = 1.0;
  status = aim_convert(aimInfo, 1, bodyLunits, &Lunit, Lunits, &Lunit);
  AIM_STATUS(aimInfo, status);

  fprintf(fp, "Lunit = %lf %s\n", Lunit, Lunits);
  fprintf(fp, "Munit = %lf %s\n", Munit, Munits);
  fprintf(fp, "Tunit = %lf %s\n", Tunit, Tunits);
  fprintf(fp, "\n");
  fprintf(fp, "#-------------------------\n");
  fprintf(fp, "#  Gravity and density to be used as default values in trim setup.\n");
  fprintf(fp, "#  Must be in the units given above.\n");
/*@-nullpass@*/

  gravity = aimInputs[inGravity-1].vals.real;
  density = aimInputs[inDensity-1].vals.real;

  fprintf(fp, "g   = %lf\n", gravity/(Lunit/(Tunit*Tunit)));
  fprintf(fp, "rho = %lf\n", density/(Munit/pow(Lunit,3)));
  fprintf(fp, "\n");
  fprintf(fp, "#-------------------------\n");
  fprintf(fp, "#  Mass & Inertia breakdown.\n");
  fprintf(fp, "#  x y z  is location of item's own CG.\n");
  fprintf(fp, "#  Ixx... are item's inertias about item's own CG.\n");
  fprintf(fp, "#\n");
  fprintf(fp, "#  x,y,z system here must be exactly the same one used in the AVL input file\n");
  fprintf(fp, "#     (same orientation, same origin location, same length units)\n");
  fprintf(fp, "#\n");
  fprintf(fp, "#  mass     x     y     z       Ixx    Iyy    Izz   [ Ixy  Ixz  Iyz ]\n");
  fprintf(fp, "#\n");

  massProp    = aimInputs[inMassProp-1].vals.tuple;
  massPropLen = aimInputs[inMassProp-1].length;

  status = aim_unitRaise(aimInfo, Lunits, 2, &tmpUnits ); // length^2
  AIM_STATUS(aimInfo, status);
  status = aim_unitMultiply(aimInfo, Munits, tmpUnits, &Iunits ); // mass*length^2, e.g moment of inertia
  AIM_STATUS(aimInfo, status);
  AIM_FREE(tmpUnits);
  AIM_NOTNULL(Iunits, aimInfo, status);
/*@+nullpass@*/
  MLL = Munit * Lunit * Lunit;

  //status = writeMassProp(aimInfo, aimInputs, fp);
  //AIM_STATUS(aimInfo, status);

  printf("Parsing MassProp\n");
  for (i = 0; i< massPropLen; i++ ) {

      // Set error message strings
      errName  = massProp[i].name;
      errValue = massProp[i].value;

      // Do we have a json string?
       if (strncmp( massProp[i].value, "{", 1) != 0) {
           errMsg = "  MassProp tuple value is expected to be a JSON string dictionary";
           status = CAPS_BADVALUE;
           goto cleanup;
       }

       keyWord = "mass";
       status  = search_jsonDictionary(massProp[i].value, keyWord, &keyValue);
       AIM_STATUS(aimInfo, status);
       AIM_NOTNULL(keyValue, aimInfo, status);

       if (units->mass != NULL)
         status = string_toDoubleUnits(aimInfo, keyValue, Munits, &mass);
       else
         status = string_toDouble(keyValue, &mass);
       AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);

       AIM_FREE(keyValue);
       (void) string_freeArray(numString, &stringArray);

       keyWord = "CG";
       status  = search_jsonDictionary(massProp[i].value, keyWord, &keyValue);
       AIM_STATUS(aimInfo, status);
       AIM_NOTNULL(keyValue, aimInfo, status);

       if (units->length != NULL)
         status = string_toDoubleArrayUnits(aimInfo, keyValue, Lunits, 3, xyz);
       else
         status = string_toDoubleArray(keyValue, 3, xyz);
       AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);

       AIM_FREE(keyValue);
       (void) string_freeArray(numString, &stringArray);

       for (j = 0; j< 6; j++) inertia[j] = 0;

       keyWord = "massInertia";
       status  = search_jsonDictionary(massProp[i].value, keyWord, &keyValue);
       AIM_STATUS(aimInfo, status);
       AIM_NOTNULL(keyValue, aimInfo, status);

       if (units->length != NULL)
         status = string_toDoubleDynamicArrayUnits(aimInfo, keyValue, Iunits, &inertiaLen, &I);
       else
         status = string_toDoubleDynamicArray(keyValue, &inertiaLen, &I);
       AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);

       AIM_NOTNULL(I, aimInfo, status);

       // Inertia order = Ixx, Iyy, Izz, Ixy, Ixz, Iyz
       for (j = 0; j < inertiaLen; j++) inertia[j] = I[j];

       AIM_FREE(I);
       AIM_FREE(keyValue);
       (void) string_freeArray(numString, &stringArray);

      // Normalize away the Lunit, Munit, Tunit values as AVL will multiply everything by those values
      fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ! %s\n",
                  mass/Munit,
                  xyz[0]/Lunit, xyz[1]/Lunit, xyz[2]/Lunit,
                  inertia[0]/MLL, inertia[1]/MLL, inertia[2]/MLL, inertia[3]/MLL,
                  inertia[4]/MLL, inertia[5]/MLL,
                  massProp[i].name);
  }

  status = CAPS_SUCCESS;

cleanup:

  if (status != CAPS_SUCCESS && errName != NULL && errValue != NULL) {
      AIM_ERROR(  aimInfo, "Cannot parse mass properties for:\n");
      AIM_ADDLINE(aimInfo, "  (\"%s\", %s)\n", errName, errValue);
      if (errMsg != NULL && keyWord != NULL) {
          AIM_ADDLINE(aimInfo, "\n");
          AIM_ADDLINE(aimInfo, "%s for %s\n", errMsg, keyWord);
      } else if (errMsg != NULL) {
          AIM_ADDLINE(aimInfo, "\n");
          AIM_ADDLINE(aimInfo, "%s\n", errMsg);
      }
      AIM_ADDLINE(aimInfo, "\n");
      AIM_ADDLINE(aimInfo, "  The 'value' JSON string should be of the form:\n");
      AIM_ADDLINE(aimInfo, "\t{\"mass\":[mass,\"kg\"], \"CG\":[[x,y,z],\"m\"], \"massInertia\":[[Ixx, Iyy, Izz, Ixy, Ixz, Iyz], \"kg*m2\"]}\n");
      AIM_ADDLINE(aimInfo, "or\n");
      AIM_ADDLINE(aimInfo, "\t{\"mass\":mass, \"CG\":[x,y,z], \"massInertia\":[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]}\n");
      AIM_ADDLINE(aimInfo, "if no unitSystem is specified\n");
  }

  // close the file
  if (fp != NULL) fclose(fp);
  fp = NULL;

  AIM_FREE(I);
  AIM_FREE(value);
  AIM_FREE(Iunits);
  AIM_FREE(tmpUnits);

  AIM_FREE(keyValue);
  string_freeArray(numString, &stringArray);

  AIM_FREE(I);

  return status;
}


static int read_Data(void *aimInfo, const char *file, const char *key, double *data)
{

    int i; //Indexing

    size_t     linecap = 0;
    char       *valstr = NULL, *line = NULL;
    FILE       *fp = NULL;

    *data = 0.0;

    if (file == NULL) return CAPS_NULLVALUE;
    if (key  == NULL) return CAPS_NULLVALUE;

    // Open the AVL output file
    fp = aim_fopen(aimInfo, file, "r");
    if (fp == NULL) {
        return CAPS_DIRERR;
    }

    // Scan the file for the string
    valstr = NULL;

    while (getline(&line, &linecap, fp) >= 0) {

        if (line == NULL) continue;

        // printf("Key = %s\n", key);

        valstr = strstr(line, key);

        if (valstr != NULL) {

            // Check string for ***
            for (i = 0; i < strlen(valstr); i++) {
                if (valstr[i] == '*'){
                    *data = 0.0;

                    printf("Unreal value for variable %s - setting it to 0.0\n", key);
                    break;
                }
            }
            // Found it -- get the value
            sscanf(&valstr[strlen(key)+1], "%lf", data);
            break;
        }
    }

    // Restore the path we came in with
    fclose(fp);

    AIM_FREE(line);

    return CAPS_SUCCESS;
}


static int
read_StripForces(void *aimInfo, int *length, capsTuple **surfaces_out)
{

    int status = CAPS_SUCCESS;
    int i; //Indexing

    size_t     linecap = 0, vallen = 0, alen = 0;
    char       *str = NULL, *line = NULL, *rest = NULL, *token = NULL, *value = NULL;
    capsTuple  *tuples = NULL, *surfaces = NULL;
    int        numSurfaces = 0, numDataEntry = 0;
    FILE       *fp = NULL;

    // Open the AVL output file
    fp = aim_fopen(aimInfo, "capsStripForce.txt", "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Failed to open 'capsStripForce.txt'");
        return CAPS_IOERR;
    }

    while (true) {

        // Scan for the next surface
        while ((status = getline(&line, &linecap, fp)) >= 0) {
            if (line == NULL) continue;
            str = strstr(line, "Surface #");
            if (str != NULL) break;
        }
        if (status < 0 || str == NULL) {
          // reached end of file without finding another surface
          break;
        }

        // remove line feed charachters
        while ((str = strstr(line, "\n")) != NULL) str[0] = ' ';
        while ((str = strstr(line, "\r")) != NULL) str[0] = ' ';

        rest = strstr(line, "#");
       (void)   strtok_r(rest, " ", &rest); // skip the #
        token = strtok_r(rest, " ", &rest); // skip the surface number
        while(rest[0]              == ' ' && rest[0] != '\0') rest++;                      // remove leading spaces
        while(rest[strlen(rest)-1] == ' '                   ) rest[strlen(rest)-1] = '\0'; // remove trailing spaces
        if (rest[0] == '\0') {
            printf("ERROR: Could not find a strip force surface name for surface # %s!\n", token);
            status = CAPS_IOERR;
            goto cleanup;
        }

        // create a new surface and save off the name of the surface, and start the value JSON dictionary
        numSurfaces++;
        surfaces = (capsTuple *)EG_reall(surfaces, sizeof(capsTuple)*(numSurfaces));
        if (surfaces == NULL) { status = EGADS_MALLOC; goto cleanup; }

        surfaces[numSurfaces-1].name  = (char *) EG_alloc((strlen(rest)+1)*sizeof(char));
        surfaces[numSurfaces-1].value = (char *) EG_alloc((strlen("{")+1)*sizeof(char));
        if (surfaces[numSurfaces-1].name  == NULL) { status = EGADS_MALLOC; goto cleanup; }
        if (surfaces[numSurfaces-1].value == NULL) { status = EGADS_MALLOC; goto cleanup; }

        strcpy(surfaces[numSurfaces-1].name, rest);
        strcpy(surfaces[numSurfaces-1].value, "{");

        // Scan the file for the line just before the header string
        while ((status = getline(&line, &linecap, fp)) >= 0) {
            if (line == NULL) continue;
            str = strstr(line, "Strip Forces referred");
            if (str != NULL) break;
        }
        if (status < 0){
            status = CAPS_IOERR;
            goto cleanup;
        }

        // read the header, e.g.
        //   j      Yle    Chord     Area     c cl      ai      cl_norm  cl       cd       cdv    cm_c/4    cm_LE  C.P.x/c
        status = getline(&line, &linecap, fp);
        if (status < 0){
            status = CAPS_IOERR;
            goto cleanup;
        }

        // remove line feed charachters
        while ((str = strstr(line, "\n")) != NULL) str[0] = ' ';
        while ((str = strstr(line, "\r")) != NULL) str[0] = ' ';

        // change "c cl" to "c_cl" so it does not get treated as two headers
        str = strstr(line, "c cl");
        if (str == NULL){
            status = CAPS_IOERR;
            goto cleanup;
        }
        str[1] = '_';

        rest = line;
        // skip the "j" column
        (void) strtok_r(rest, " ", &rest);

        numDataEntry = 0;
        while ((token = strtok_r(rest, " ", &rest))) {
            // resize the number of tuples
            AIM_REALL(tuples, numDataEntry+1, capsTuple, aimInfo, status);
            tuples[numDataEntry].name = tuples[numDataEntry].value = NULL;

            AIM_ALLOC(tuples[numDataEntry].name , strlen(token)+1, char, aimInfo, status);
            AIM_ALLOC(tuples[numDataEntry].value, strlen("["  )+1, char, aimInfo, status);

            // restore the name with a space in it
            str = strstr(token, "c_cl");
            if (str != NULL){
              str[1] = ' ';
            }

            // copy over the header and start the list
            strcpy(tuples[numDataEntry].name, token);
            strcpy(tuples[numDataEntry].value, "[");
            numDataEntry++;
        }
        if ((numDataEntry == 0) || (tuples == NULL) || (surfaces == NULL)) {
            AIM_ERROR(aimInfo, "No Tuples/Surfaces!");
            status = CAPS_IOERR;
            goto cleanup;
        }

        // read in the data columns
        while ((status = getline(&line, &linecap, fp)) >= 0) {

          str = strstr(line, "--------------");
          if (str != NULL) break; // end of the file

          // remove line feed characters
          while ((str = strstr(line, "\n")) != NULL) str[0] = ' ';
          while ((str = strstr(line, "\r")) != NULL) str[0] = ' ';

          rest = line;
          // skip the "j" column
          token = strtok_r(rest, " ", &rest);

          if (token == NULL) break; // end of data for this surface

          i = 0;
          while ((token = strtok_r(rest, " ", &rest))) {
              vallen = strlen(tuples[i].value);
              AIM_REALL(tuples[i].value, vallen+strlen(token)+2, char, aimInfo, status);

              // append the values to the list
              sprintf(tuples[i].value + vallen,"%s,", token);
              i++;
          }
        }
        if (status < 0){
            status = CAPS_IOERR;
            goto cleanup;
        }

        for (i = 0; i < numDataEntry; i++) {

            // close the lists by replacing the ',' with ']'
            tuples[i].value[strlen(tuples[i].value)-1] = ']';

            // collapse down the tuples for the surface
            vallen = strlen(surfaces[numSurfaces-1].value);
            // \" + name                     + \": + value                 + ,\0
            alen = 1 + strlen(tuples[i].name) + 2 + strlen(tuples[i].value) + 2;
            AIM_REALL(surfaces[numSurfaces-1].value, vallen+alen, char, aimInfo, status);

            value = surfaces[numSurfaces-1].value + vallen;
            sprintf(value,"\"%s\":%s,", tuples[i].name, tuples[i].value);

            // release the memory now that it's been consumed
            AIM_FREE(tuples[i].name );
            AIM_FREE(tuples[i].value);
        }
        AIM_FREE(tuples);
        numDataEntry = 0;

        // close the JSON dict by replacing the ',' with '}'
        surfaces[numSurfaces-1].value[strlen(surfaces[numSurfaces-1].value)-1] = '}';
    }

    *surfaces_out = surfaces;
    surfaces      = NULL;
    *length       = numSurfaces;
    status        = CAPS_SUCCESS;

cleanup:

    if (fp != NULL) fclose(fp);

    if ((status != CAPS_SUCCESS) && (surfaces != NULL)) {
        for (i = 0; i < numSurfaces; i++) {
            AIM_FREE(surfaces[i].name);
            AIM_FREE(surfaces[i].value);
        }
        AIM_FREE(surfaces);
    }

    if (tuples != NULL) {
        for (i = 0; i < numDataEntry; i++) {
            AIM_FREE(tuples[i].name);
            AIM_FREE(tuples[i].value);
        }
        AIM_FREE(tuples);
    }

    AIM_FREE(line);

    return status;
}


static int
read_EigenValues(void *aimInfo, int *length, capsTuple **eigen_out)
{

    int status = CAPS_SUCCESS;
    int i; //Indexing

    size_t     linecap = 0, vallen = 0;
    char       *str = NULL, *line = NULL, *rest = NULL, *token = NULL;
    char       caseName[30];
    capsTuple  *eigen = NULL;
    int        numCase = 0, icase = 0;
    FILE       *fp = NULL;
    char eigenValueFile[] = "capsEigenValues.txt";

    // ignore if file does not exist
    if (aim_isFile(aimInfo, eigenValueFile) == CAPS_NOTFOUND) {
        return CAPS_SUCCESS;
    }

    // Open the eigen value output file
    fp = aim_fopen(aimInfo, eigenValueFile, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Failed to open '%s'", eigenValueFile);
        return CAPS_IOERR;
    }

    // Scan for header
    while ((status = getline(&line, &linecap, fp)) >= 0) {
        if (line == NULL) continue;
        str = strstr(line, "#   Run case");
        if (str != NULL) break;
    }
    if (status < 0 || str == NULL) {
      // reached end of file without finding any eigenvalues
      status = CAPS_SUCCESS;
      goto cleanup;
    }

    while (getline(&line, &linecap, fp) >= 0) {

        // remove line feed charachters
        while ((str = strstr(line, "\n")) != NULL) str[0] = ' ';
        while ((str = strstr(line, "\r")) != NULL) str[0] = ' ';

        rest = line;
        token = strtok_r(rest, " ", &rest); // get the case number

        icase = atoi(token);
        if (icase > numCase) {
            // create a new case to save off the eigen informaiton
            AIM_REALL(eigen, icase, capsTuple, aimInfo, status);
            for (i = numCase; i < icase; i++)
              eigen[icase-1].name = eigen[icase-1].value = NULL;

            sprintf(caseName,"case %d", icase );
            AIM_STRDUP(eigen[icase-1].name , caseName, aimInfo, status);
            AIM_STRDUP(eigen[icase-1].value, "["     , aimInfo, status);
            numCase = icase;
        }

        if (eigen != NULL)
            for (i = 0; i < 2; i++) {
                token = strtok_r(rest, " ", &rest); // get the real/imaginary part of the eigen value

                vallen = strlen(eigen[icase-1].value);
                AIM_REALL(eigen[icase-1].value, vallen+strlen(token)+3, char, aimInfo, status);

                // append the values to the list
                if (i == 0)
                    sprintf(eigen[icase-1].value + vallen,"[%s,", token );
                else
                    sprintf(eigen[icase-1].value + vallen,"%s],", token );
            }
    }

    // close the array by replacing the ',' with ']'
    if (eigen != NULL)
        for (icase = 0; icase < numCase; icase++)
            eigen[icase].value[strlen(eigen[icase].value)-1] = ']';

    *eigen_out = eigen;
    eigen      = NULL;
    *length    = numCase;
    status     = CAPS_SUCCESS;

cleanup:

    if (fp != NULL) fclose(fp);

    if ((status != CAPS_SUCCESS) && (eigen != NULL)) {
        for (i = 0; i < numCase; i++) {
            AIM_FREE(eigen[i].name);
            AIM_FREE(eigen[i].value);
        }
        AIM_FREE(eigen);
    }

    AIM_FREE(line);

    return status;
}

static int get_controlDeriv(void *aimInfo, int controlIndex, int outputIndex, double *data)
{

    int status; // Function return status
    char *fileToOpen;

    char *coeff = NULL;
    char key[42];

    // Stability axis
    if        (outputIndex == CLtot) {
        fileToOpen = "capsStatbilityDeriv.txt";
        coeff = "CL";

    } else if (outputIndex == CYtot) {
        fileToOpen = "capsStatbilityDeriv.txt";
        coeff = "CY";

    } else if (outputIndex == ClPtot) {
        fileToOpen = "capsStatbilityDeriv.txt";
        coeff = "Cl";

    } else if (outputIndex == Cmtot) {
        fileToOpen = "capsStatbilityDeriv.txt";
        coeff = "Cm";

    } else if (outputIndex == CnPtot) {
        fileToOpen = "capsStatbilityDeriv.txt";
        coeff = "Cn";

    // Body axis
    } else if (outputIndex == CXtot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "CX";

    } else if (outputIndex == CYtot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "CY";

    } else if (outputIndex == CZtot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "CZ";

    } else if (outputIndex == Cltot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "Cl";

    } else if (outputIndex == Cmtot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "Cm";

    } else if (outputIndex == Cntot) {
        fileToOpen = "capsBodyAxisDeriv.txt";
        coeff = "Cn";

    } else {
        AIM_ERROR(aimInfo, "Unrecognized output variable for control derivatives!");
        status = CAPS_MISMATCH;
        goto cleanup;
    }

    sprintf(key, "%sd%d =", coeff, controlIndex);

    status = read_Data(aimInfo, fileToOpen, key, data);
    AIM_STATUS (aimInfo, status);

    status = CAPS_SUCCESS;

cleanup:
    return status;
}


static int parse_controlName(aimStorage *avlInstance, char string[],
                             int *controlNumber)
{

    int status; // Function return status

    char control[] = "AVL_Control";

    char *controlName = NULL;

    *controlNumber = CAPSMAGIC;

    //printf("String = %s\n", string);

    if (strncasecmp(string, control, strlen(control)) != 0 ) {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    controlName = strstr(string, ":");
    if (controlName == NULL) {
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Increment pointer to remove ':'
    ++controlName;

    //printf("Temp = %s\n", temp);

    // Loop through and determine which control integer this name corresponds to?
    status = get_mapAttrToIndexIndex(&avlInstance->controlMap,
            (const char *) controlName, controlNumber);
    if (status != CAPS_SUCCESS) goto cleanup;

    /*
    printf("ControlName = %s\n", controlName);
    printf("ControlIndex = %d\n", *controlNumber);
     */

    status = CAPS_SUCCESS;
    goto cleanup;

cleanup:
    return status;
}


static int
stabilityAngleDerivatives(void *aimInfo, capsValue *val)
{
  int i, status = CAPS_SUCCESS;

  int nderiv = 2;

  AIM_REALL(val->derivs, val->nderiv+nderiv, capsDeriv, aimInfo, status);
  for (i = val->nderiv; i < val->nderiv+nderiv; i++) {
    val->derivs[i].name    = NULL;
    val->derivs[i].deriv   = NULL;
    val->derivs[i].len_wrt = 1;
    AIM_ALLOC(val->derivs[i].deriv, 1, double, aimInfo, status);
  }

  i = val->nderiv;
  AIM_STRDUP(val->derivs[i+0].name, "Alpha", aimInfo, status);
  AIM_STRDUP(val->derivs[i+1].name, "Beta", aimInfo, status);
  val->nderiv += nderiv;

cleanup:
  return status;
}


static int
bodyRateDerivatives(void *aimInfo, capsValue *val)
{
  int i, status = CAPS_SUCCESS;

  int nderiv = 3;

  AIM_REALL(val->derivs, val->nderiv+nderiv, capsDeriv, aimInfo, status);
  for (i = val->nderiv; i < val->nderiv+nderiv; i++) {
    val->derivs[i].name    = NULL;
    val->derivs[i].deriv   = NULL;
    val->derivs[i].len_wrt = 1;
    AIM_ALLOC(val->derivs[i].deriv, 1, double, aimInfo, status);
  }

  i = val->nderiv;
  AIM_STRDUP(val->derivs[i+0].name, "RollRate", aimInfo, status);
  AIM_STRDUP(val->derivs[i+1].name, "PitchRate", aimInfo, status);
  AIM_STRDUP(val->derivs[i+2].name, "YawRate", aimInfo, status);
  val->nderiv += nderiv;

cleanup:
  return status;
}


static int
controlDerivatives(void *aimInfo, int outIndex, aimStorage *avlInstance, capsValue *val)
{
  int i, status = CAPS_SUCCESS;

  int nderiv = avlInstance->controlMap.numAttribute;

  if (nderiv == 0) return CAPS_SUCCESS;

  AIM_REALL(val->derivs, val->nderiv+nderiv, capsDeriv, aimInfo, status);
  for (i = val->nderiv; i < val->nderiv+nderiv; i++) {
    val->derivs[i].name    = NULL;
    val->derivs[i].deriv   = NULL;
    val->derivs[i].len_wrt = 1;
    AIM_ALLOC(val->derivs[i].deriv, 1, double, aimInfo, status);
  }

  // Loop through control surfaces
  for (i = 0; i < avlInstance->controlMap.numAttribute; i++) {
    AIM_STRDUP(val->derivs[val->nderiv+i].name, avlInstance->controlMap.attributeName[i], aimInfo, status);

    status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                              outIndex, val->derivs[val->nderiv+i].deriv);
    AIM_STATUS(aimInfo, status);
  }
  val->nderiv += nderiv;

cleanup:
  return status;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS;
    const char *keyWord;
    char *keyValue = NULL;
    double real = 1;
    cfdUnitsStruct *units=NULL;

    aimStorage *avlInstance=NULL;

#ifdef DEBUG
    printf("\n avlAIM/aimInitialize   instance = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks  = NULL;
    *fInOut  = NULL;

    // Allocate avlInstance
    AIM_ALLOC(avlInstance, 1, aimStorage, aimInfo, status);
    *instStore = avlInstance;

    // Initiate storage
    status = initiate_aimStorage(aimInfo, avlInstance);
    AIM_STATUS(aimInfo, status);

    /*! \page aimUnitsAVL AIM Units
     *  A unit system may be optionally specified during AIM instance initiation. If
     *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
     *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
     *  unit consistent. A unit system may be specified via a JSON string dictionary for example:
     *  unitSys = "{"mass": "kg", "length": "m", "time":"seconds", "temperature": "K"}"
     */
    if (unitSys != NULL) {
      units = &avlInstance->units;

      // Do we have a json string?
      if (strncmp( unitSys, "{", 1) != 0) {
        AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
        return CAPS_BADVALUE;
      }

      /*! \page aimUnitsAVL
       *  \section jsonStringAVL JSON String Dictionary
       *  The key arguments of the dictionary are described in the following:
       *
       *  <ul>
       *  <li> <B>mass = "None"</B> </li> <br>
       *  Mass units - e.g. "kilogram", "k", "slug", ...
       *  </ul>
       */
      keyWord = "mass";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->mass = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->mass, &real, "kg", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->mass, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsAVL
       *  <ul>
       *  <li> <B>length = "None"</B> </li> <br>
       *  Length units - e.g. "meter", "m", "inch", "in", "mile", ...
       *  </ul>
       */
      keyWord = "length";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->length = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->length, &real, "m", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->length, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsAVL
       *  <ul>
       *  <li> <B>time = "None"</B> </li> <br>
       *  Time units - e.g. "second", "s", "minute", ...
       *  </ul>
       */
      keyWord = "time";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->time = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->time, &real, "s", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->time, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsAVL
       *  <ul>
       *  <li> <B>temperature = "None"</B> </li> <br>
       *  Temperature units - e.g. "Kelvin", "K", "degC", ...
       *  </ul>
       */
      keyWord = "temperature";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->temperature = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->temperature, &real, "K", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->temperature, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      status = cfd_cfdDerivedUnits(aimInfo, units);
      AIM_STATUS(aimInfo, status);
    }

cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
    /*! \page aimInputsAVL AIM Inputs
     * The following list outlines the AVL inputs along with their default value available
     * through the AIM interface.
     */

    int status = CAPS_SUCCESS;
    aimStorage *avlInstance;
    cfdUnitsStruct *units=NULL;

#ifdef DEBUG
    printf(" avlAIM/aimInputs  index = %d!\n", index);
#endif

    avlInstance = (aimStorage *) instStore;

    if (avlInstance != NULL) units = &avlInstance->units;

    if (index == inMach) {
        *ainame           = EG_strdup("Mach");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsAVL
         * - <B> Mach = 0.0 </B> <br>
         *  Mach number.
         */

    } else if (index == inAlpha) {
        *ainame           = EG_strdup("Alpha");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->nrow      = 1;
        defval->ncol      = 1;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        defval->lfixed    = Change;
        defval->sfixed    = Change;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, "degree", aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B> Alpha = NULL </B> <br>
         *  Angle of attack [degree]. Either CL or Alpha must be defined but not both.
         */

    } else if (index == inBeta) {
        *ainame           = EG_strdup("Beta");
        defval->type      = Double;
        defval->vals.real = 0.0;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, "degree", aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B> Beta = 0.0 </B> <br>
         *  Sideslip angle [degree].
         */

    } else if (index == inRollRate) {
        *ainame           = EG_strdup("RollRate");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsAVL
         * - <B> RollRate = 0.0 </B> <br>
         *  Non-dimensional roll rate.
         */

    } else if (index == inPitchRate) {
        *ainame           = EG_strdup("PitchRate");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsAVL
         * - <B> PitchRate = 0.0 </B> <br>
         *  Non-dimensional pitch rate.
         */

    } else if (index == inYawRate) {
        *ainame           = EG_strdup("YawRate");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsAVL
         * - <B> YawRate = 0.0 </B> <br>
         *  Non-dimensional yaw rate.
         */

    } else if (index == inCDp) {
        *ainame           = EG_strdup("CDp");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->units     = NULL;

        /*! \page aimInputsAVL
         * - <B> CDp = 0.0 </B> <br>
         *  A fixed value of profile drag to be added to all simulations.
         */

    } else if (index == inAVL_Surface) {
        *ainame              = EG_strdup("AVL_Surface");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAVL
         * - <B>AVL_Surface = NULL </B> <br>
         * See \ref vlmSurface for additional details.
         */

    } else if (index == inAVL_Control) {
        *ainame              = EG_strdup("AVL_Control");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAVL
         * - <B>AVL_Control = NULL </B> <br>
         * See \ref vlmControl for additional details.
         */

    } else if (index == inCL) {
        *ainame           = EG_strdup("CL");
        defval->type      = Double;
        defval->dim       = Vector;
        defval->nrow      = 1;
        defval->ncol      = 1;
        defval->units     = NULL;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        defval->lfixed    = Change;

        /*! \page aimInputsAVL
         * - <B> CL = NULL </B> <br>
         *  Coefficient of Lift.  AVL will solve for Angle of Attack.  Either CL or Alpha must be defined but not both.
         */

    } else if (index == inMoment_Center) {
        *ainame              = EG_strdup("Moment_Center");
        defval->type          = Double;
        defval->dim           = Vector;
        defval->nrow          = 3;
        defval->ncol          = 1;
        AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
        defval->vals.reals[0] = 0.0;
        defval->vals.reals[1] = 0.0;
        defval->vals.reals[2] = 0.0;
        defval->nullVal       = IsNull;
        defval->lfixed        = Fixed;
        defval->sfixed        = Fixed;
        if (units != NULL && units->length != NULL) {
            AIM_STRDUP(defval->units, units->length, aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B>Moment_Center = NULL, [0.0, 0.0, 0.0]</B> <br>
         * Array values correspond to the Xref, Yref, and Zref variables.
         * Alternatively, the geometry (body) attributes "capsReferenceX", "capsReferenceY",
         * and "capsReferenceZ" may be used to specify the X-, Y-, and Z- reference centers, respectively
         * (note: values set through the AIM input will supersede the attribution values).
         */

    } else if (index == inMassProp) {
        *ainame              = EG_strdup("MassProp");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAVL
         * - <B>MassProp = NULL</B> <br>
         * Mass properties used for eigen value analysis
         * Structure for the mass property tuple = ("Name", "Value").
         * The "Name" of the mass component used for documenting the xxx.mass file.
         * The value is a JSON dictionary with values with unit pairs for mass, CG, and moments of inertia information
         * (e.g. "Value" = {"mass" : [mass,"kg"], "CG" : [[x,y,z],"m"], "massInertia" : [[Ixx, Iyy, Izz, Ixy, Ixz, Iyz], "kg*m2"]})
         * The components Ixy, Ixz, and Iyz are optional may be omitted.
         * Must be in units of kg, m, and kg*m^2 if unitSystem (see \ref aimUnitsAVL) is not specified and no units should be specified in the JSON dictionary.
         */

    } else if (index == inGravity) {
        *ainame           = EG_strdup("Gravity");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        if (units != NULL && units->acceleration != NULL) {
            AIM_STRDUP(defval->units, units->acceleration, aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B> Gravity = NULL </B> <br>
         *  Magnitude of the gravitational force used for Eigen value analysis.
         *  Must be in units of m/s^2 if unitSystem (see \ref aimUnitsAVL) is not specified.
         */

    } else if (index == inDensity) {
        *ainame           = EG_strdup("Density");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        if (units != NULL && units->density != NULL) {
            AIM_STRDUP(defval->units, units->density, aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B> Density = NULL </B> <br>
         *  Air density used for Eigen value analysis.
         *  Must be in units of kg/m^3 if unitSystem (see \ref aimUnitsAVL) is not specified.
         */

    } else if (index == inVelocity) {
        *ainame           = EG_strdup("Velocity");
        defval->type      = Double;
        defval->vals.real = 0.0;
        defval->nullVal   = IsNull;
        if (units != NULL && units->speed != NULL) {
            AIM_STRDUP(defval->units, units->speed, aimInfo, status);
        }

        /*! \page aimInputsAVL
         * - <B> Velocity = NULL </B> <br>
         *  Velocity used for Eigen value analysis.
         *  Must be in units of m/s if unitSystem (see \ref aimUnitsAVL) is not specified.
         */
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Function return status

    int i, j, k, surf, section, control; // Indexing

    int numBody;

    int found = (int) false; // Boolean tester
    int eigenValues = (int) false;

    int atype, alen;

    aimStorage *avlInstance;

    double      Sref, Cref, Bref, Xref, Yref, Zref;
    int         foundSref=(int)false, foundCref=(int)false, foundBref=(int)false, foundXref=(int)false;
    const int    *ints;
    const char   *string, *intents;
    const double *reals;
    ego          *bodies;

    const char *bodyLunits = NULL;

    // File I/O
    FILE *fp = NULL;
    char inputFilename[] = "avlInput.txt";
    char avlFilename[] = "caps.avl";
    char massFilename[] = "caps.mass";

    char totalForceFile[] = "capsTotalForce.txt";
    char stripForceFile[] = "capsStripForce.txt";
    char stabilityFile[] = "capsStatbilityDeriv.txt";
    char bodyAxisFile[] = "capsBodyAxisDeriv.txt";
    char hingeMomentFile[] = "capsHingeMoment.txt";
    char eigenValueFile[] = "capsEigenValues.txt";

    int numControlName = 0; // Unique control names
    char **controlName = NULL;

    int numSpanWise = 0;

    // Control related variables
    //double xyzHingeMag, xyzHingeVec[3];


#ifdef DEBUG
    printf(" avlAIM/aimPreAnalysis\n");
#endif

    avlInstance = (aimStorage *) instStore;

    // Initialize reference values
    Sref = 1.0;
    Cref = 1.0;
    Bref = 1.0;

    Xref = 0.0;
    Yref = 0.0;
    Zref = 0.0;

    if (aimInputs == NULL) {
#ifdef DEBUG
        printf(" avlAIM/aimPreAnalysis aimInputs == NULL!\n");
#endif
        return CAPS_NULLVALUE;
    }

    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody == 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        status = CAPS_SOURCEERR;
        goto cleanup;
    }

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS) {
      // Destroy previous groupMap (in case it already exists)
      status = destroy_mapAttrToIndexStruct(&avlInstance->groupMap);
      AIM_STATUS(aimInfo, status);

      // Get capsGroup name and index mapping to make sure all bodies have a capsGroup value
      status = create_CAPSGroupAttrToIndexMap(numBody,
                                              bodies,
                                              0, // Only search down to the body level of the EGADS body
                                              &avlInstance->groupMap);
      AIM_STATUS(aimInfo, status);
    }


    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Surface) == CAPS_SUCCESS) {
        // Get AVL surface information
        if (aimInputs[inAVL_Surface-1].nullVal == NotNull) {

            if (avlInstance->avlSurface != NULL) {
                for (i = 0; i < avlInstance->numAVLSurface; i++) {
                    status = destroy_vlmSurfaceStruct(&avlInstance->avlSurface[i]);
                    AIM_STATUS(aimInfo, status);
                }
                AIM_FREE(avlInstance->avlSurface);
                avlInstance->numAVLSurface = 0;
            }

            status = get_vlmSurface(aimInputs[inAVL_Surface-1].length,
                                    aimInputs[inAVL_Surface-1].vals.tuple,
                                    &avlInstance->groupMap,
                                    1.0, // default Cspace
                                    &avlInstance->numAVLSurface,
                                    &avlInstance->avlSurface);
            AIM_STATUS(aimInfo, status);

        } else {
            AIM_ERROR(aimInfo, "No AVL_Surface specified!");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        // Accumulate section data
        status = vlm_getSections(aimInfo, numBody, bodies, NULL, avlInstance->groupMap, vlmGENERIC,
                                 avlInstance->numAVLSurface, &avlInstance->avlSurface);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(avlInstance->avlSurface, aimInfo, status);
    }

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Control) == CAPS_SUCCESS) {

        // Get AVL control surface information
        if (aimInputs[inAVL_Control-1].nullVal == NotNull) {

            status = get_vlmControl(aimInfo,
                                    aimInputs[inAVL_Control-1].length,
                                    aimInputs[inAVL_Control-1].vals.tuple,
                                    &avlInstance->numAVLControl,
                                    &avlInstance->avlControl);
            AIM_STATUS(aimInfo, status);
        }

        // Loop through surfaces and transfer control surface data to sections
        for (surf = 0; surf < avlInstance->numAVLSurface; surf++) {
    /*@-nullpass@*/
            status = get_ControlSurface(aimInfo,
                                        avlInstance->numAVLControl,
                                        avlInstance->avlControl,
                                        &avlInstance->avlSurface[surf]);
    /*@+nullpass@*/
            AIM_STATUS (aimInfo, status);
        }
    }

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Surface) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Control) == CAPS_SUCCESS) {

        // Compute auto spacing
        for (surf = 0; surf < avlInstance->numAVLSurface; surf++) {

            if      (avlInstance->avlSurface[surf].NspanTotal > 0  && avlInstance->avlSurface[surf].NspanSection == 0)
                numSpanWise = avlInstance->avlSurface[surf].NspanTotal;
            else if (avlInstance->avlSurface[surf].NspanTotal == 0 && avlInstance->avlSurface[surf].NspanSection > 0 )
                numSpanWise = (avlInstance->avlSurface[surf].numSection-1)*avlInstance->avlSurface[surf].NspanSection;
            else {
                AIM_ERROR  (aimInfo,"Only one of numSpanTotal and numSpanPerSection must be non-zero!");
                AIM_ADDLINE(aimInfo,"       numSpanTotal      = %d", avlInstance->avlSurface[surf].NspanTotal);
                AIM_ADDLINE(aimInfo,"       numSpanPerSection = %d", avlInstance->avlSurface[surf].NspanSection);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            status = vlm_autoSpaceSpanPanels(aimInfo,
                                             numSpanWise, avlInstance->avlSurface[surf].numSection,
                                                          avlInstance->avlSurface[surf].vlmSection);
            AIM_STATUS (aimInfo, status);
        }
    }

    /*
    // Determine the hinge line for any flaps
    for (surf = 0; surf < numAVLSurface; surf++) {

        for (section = 0; section < avlSurface[surf].numSection; section++) {

            for (control = 0; control < avlSurface[surf].vlmSection[section].numControl; control++) {

                for (i = 0; i < avlSurface[surf].numSection; i++) {

                    if (i == section) continue;

                    for (j = 0; j < avlSurface[surf].vlmSection[i].numControl; j++) {

                        if( strcmp(avlSurface[surf].vlmSection[section].vlmControl[control].name,
                                   avlSurface[surf].vlmSection[i].vlmControl[j].name) == 0) {

                            xyzHingeVec[0] = avlSurface[surf].vlmSection[section].vlmControl[control].xyzHinge[0] -
                                             avlSurface[surf].vlmSection[i].vlmControl[j].xyzHinge[0];

                            xyzHingeVec[1] = avlSurface[surf].vlmSection[section].vlmControl[control].xyzHinge[1] -
                                             avlSurface[surf].vlmSection[i].vlmControl[j].xyzHinge[1];

                            xyzHingeVec[2] = avlSurface[surf].vlmSection[section].vlmControl[control].xyzHinge[2] -
                                             avlSurface[surf].vlmSection[i].vlmControl[j].xyzHinge[2];

                            xyzHingeMag = sqrt(pow(xyzHingeVec[0], 2) +
                                               pow(xyzHingeVec[1], 2) +
                                               pow(xyzHingeVec[2], 2));

                            avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[0] = xyzHingeVec[0]/xyzHingeMag;
                            avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[1] = xyzHingeVec[1]/xyzHingeMag;
                            avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[2] = xyzHingeVec[2]/xyzHingeMag;

                        }
                    }
                }
            }
        }
    }
     */

    // look for eigen value analysis
    if ( (aimInputs[inMassProp-1].nullVal == NotNull ||
          aimInputs[inGravity-1].nullVal  == NotNull ||
          aimInputs[inDensity-1].nullVal  == NotNull ||
          aimInputs[inVelocity-1].nullVal == NotNull) ) {

        // Get length units
        status = check_CAPSLength(numBody, bodies, &bodyLunits);
        if (status != CAPS_SUCCESS) {
          AIM_ERROR(aimInfo, "No units assigned *** capsLength is not set in *.csm file!");
          status = CAPS_BADVALUE;
          goto cleanup;
        }

        if ( !((aimInputs[inMassProp-1].nullVal == NotNull) &&
               (aimInputs[inGravity-1].nullVal  == NotNull) &&
               (aimInputs[inDensity-1].nullVal  == NotNull) &&
               (aimInputs[inVelocity-1].nullVal == NotNull)) ) {
            AIM_ERROR(  aimInfo, " All inputs 'MassProp', 'Gravity', 'Density', and 'Velocity'\n");
            AIM_ADDLINE(aimInfo, " must be set for AVL eigen value analysis.\n");
            AIM_ADDLINE(aimInfo, " Missing values for:\n");
            if (aimInputs[inMassProp-1].nullVal == IsNull) AIM_ADDLINE(aimInfo, "    MassProp\n");
            if (aimInputs[inGravity-1].nullVal  == IsNull) AIM_ADDLINE(aimInfo, "    Gravity\n");
            if (aimInputs[inDensity-1].nullVal  == IsNull) AIM_ADDLINE(aimInfo, "    Density\n");
            if (aimInputs[inVelocity-1].nullVal == IsNull) AIM_ADDLINE(aimInfo, "    Velocity\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        eigenValues = (int) true; // compute eigen values
    }

    // Open and write the input to control the AVL session
    fp = aim_fopen(aimInfo, inputFilename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file %s\n!", inputFilename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Print the session file for the AVL run

    fprintf(fp, "PLOP\n"); // Start with disabling graphics
    fprintf(fp, "G\n");    // toggle graphics
    fprintf(fp, "\n");     // main menu

    if (eigenValues == (int)true) {
        // apply the mass properties to move the CG
        // WARNING: mset will move the CG to (0,0,0) if no xxx.mass file is present
        fprintf(fp, "mset 0\n");
    }

    // Set operation parameters
    fprintf(fp, "OPER\n");

    if (aimInputs[inAlpha-1].nullVal ==  NotNull) {
        fprintf(fp, "A A ");//Alpha
        fprintf(fp, "%lf\n", aimInputs[1].vals.real);
    }

    if (aimInputs[inCL-1].nullVal ==  NotNull) {
        fprintf(fp, "A C ");//CL
        fprintf(fp, "%lf\n", aimInputs[9].vals.real);
    }

    fprintf(fp, "B B ");//Beta
    fprintf(fp, "%lf\n", aimInputs[inBeta-1].vals.real);

    fprintf(fp, "R R ");//Roll Rate pb/2V
    fprintf(fp, "%lf\n", aimInputs[inRollRate-1].vals.real);

    fprintf(fp, "P P ");//Pitch Rate qc/2v
    fprintf(fp, "%lf\n", aimInputs[inPitchRate-1].vals.real);

    fprintf(fp, "Y Y ");//Yaw Rate rb/2v
    fprintf(fp, "%lf\n", aimInputs[inYawRate-1].vals.real);

    // Destroy previous controlMap (in case it already exists)
    status = destroy_mapAttrToIndexStruct(&avlInstance->controlMap);
    AIM_STATUS(aimInfo, status);

    // Check for control surface information
    j = 1;
    numControlName = 0;
    for (surf = 0; surf < avlInstance->numAVLSurface; surf++) {

        for (i = 0; i < avlInstance->avlSurface[surf].numSection; i++) {

            section = avlInstance->avlSurface[surf].vlmSection[i].sectionIndex;

            for (control = 0; control < avlInstance->avlSurface[surf].vlmSection[section].numControl; control++) {

                // Check to see if control surface hasn't already been written
                found = (int) false;

                for (k = 0; k < numControlName; k++) {
                    /*@-nullderef@*/
                    if (strcmp(controlName[k],
                               avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].name) == 0) {
                        found = (int) true;
                        break;
                    }
                    /*@+nullderef@*/
                }
                if (found == (int) true) continue;

                AIM_REALL(controlName, numControlName+1, char*, aimInfo, status);
                controlName[numControlName] = NULL;
                numControlName += 1;

                AIM_STRDUP(controlName[numControlName-1],
                           avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].name, aimInfo, status);

                // Store control map for later use
                status = increment_mapAttrToIndexStruct(&avlInstance->controlMap,
                                                        controlName[numControlName-1]);
                AIM_STATUS(aimInfo, status);

                fprintf(fp, "D%d D%d %f\n",j, j, avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].deflectionAngle);
                j += 1;
            }
        }
    }

    /*if(aimInputs[7].nullVal != IsNull) {
        // write control surface settings
        if (aimInputs[7].length == 1){
            fprintf(fp, "D1 D1 ");//Pitch Rate qc/2v
            fprintf(fp, "%lf\n", aimInputs[7].vals.real);
        } else {
            for (i=0; i < aimInputs[7].length; i++) {
                fprintf(fp,"D%d D%d %lf\n",i+1,i+1,aimInputs[7].vals.reals[i]);
            }
        }
    }*/

    fprintf(fp, "M\n"); // Modify parameters
    fprintf(fp, "MN\n");//Mach
    fprintf(fp, "%lf\n", aimInputs[inMach-1].vals.real);

    if (aimInputs[inVelocity-1].nullVal == NotNull) {
        fprintf(fp, "V\n");                            // Velocity
        fprintf(fp, "%lf\n", aimInputs[inVelocity-1].vals.real);  // set the value
    }
    fprintf(fp, "\n");                    // exit modify parameters

    fprintf(fp, "X\n"); // Execute the calculation

    fprintf(fp, "S\n\n"); // save caps.run file
    if (aim_isFile(aimInfo, "caps.run") == CAPS_SUCCESS) {
        fprintf(fp, "y\n");
    }

    // Get Total forces
    fprintf(fp,"FT\n");
    fprintf(fp,"%s\n", totalForceFile);
    if (aim_isFile(aimInfo, totalForceFile) == CAPS_SUCCESS) {
        fprintf(fp, "O\n");
    }

    // Get strip forces
    fprintf(fp,"FS\n");
    fprintf(fp,"%s\n", stripForceFile);
    if (aim_isFile(aimInfo, stripForceFile) == CAPS_SUCCESS) {
        fprintf(fp, "O\n");
    }

    // Get stability derivatives
    fprintf(fp,"ST\n");
    fprintf(fp,"%s\n", stabilityFile);
    if (aim_isFile(aimInfo, stabilityFile) == CAPS_SUCCESS) {
        fprintf(fp, "O\n");
    }

    // Get stability (body axis) derivatives
    fprintf(fp,"SB\n");
    fprintf(fp,"%s\n", bodyAxisFile);
    if (aim_isFile(aimInfo, bodyAxisFile) == CAPS_SUCCESS) {
        fprintf(fp, "O\n");
    }

    // Get hinge moments
    fprintf(fp,"HM\n");
    fprintf(fp,"%s\n", hingeMomentFile);
    if (aim_isFile(aimInfo, hingeMomentFile) == CAPS_SUCCESS) {
        fprintf(fp, "O\n");
    }

    fprintf(fp, "\n"); // back to main menu

    if (eigenValues == (int)true) {
        fprintf(fp, "mode\n");                // enter eigen value analysis

        fprintf(fp, "n\n");                   // compute eigen values
        fprintf(fp, "w\n");                   // write eigen values to file
        fprintf(fp, "%s\n", eigenValueFile);
        if (aim_isFile(aimInfo, eigenValueFile) == CAPS_SUCCESS) {
            fprintf(fp, "Y\n");
        }
        fprintf(fp, "\n"); // back to main menu
    }

    fprintf(fp, "Quit\n"); // Quit AVL
    fclose(fp);
    fp = NULL;

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Surface) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inAVL_Control) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inMoment_Center) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, inCDp) == CAPS_SUCCESS) {

        /* open the file and write the avl data */
        fp = aim_fopen(aimInfo, avlFilename, "w");
        if (fp == NULL) {
            AIM_ERROR(aimInfo, "Unable to open file %s\n!", avlFilename);
            status =  CAPS_IOERR;
            goto cleanup;
        }

        // Loop over bodies and look for reference quantity attributes
        for (i=0; i < numBody; i++) {
            status = EG_attributeRet(bodies[i], "capsReferenceArea",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL && alen == 1) {
                    Sref = (double) reals[0];
                    foundSref = (int)true;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceArea should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[i], "capsReferenceChord",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL && alen == 1) {
                    Cref = (double) reals[0];
                    foundCref = (int)true;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceChord should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[i], "capsReferenceSpan",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL && alen == 1) {
                    Bref = (double) reals[0];
                    foundBref = (int)true;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceSpan should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[i], "capsReferenceX",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL && alen == 1) {
                    Xref = (double) reals[0];
                    foundXref = (int)true;
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceX should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[i], "capsReferenceY",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS) {

                if (atype == ATTRREAL && alen == 1) {
                    Yref = (double) reals[0];
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceY should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            status = EG_attributeRet(bodies[i], "capsReferenceZ",
                                     &atype, &alen, &ints, &reals, &string);
            if (status == EGADS_SUCCESS){

                if (atype == ATTRREAL && alen == 1) {
                    Zref = (double) reals[0];
                } else {
                    AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }
        }

        if (foundSref == (int)false) {
            AIM_ERROR(aimInfo, "capsReferenceArea is not set on any body!");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
        if (foundCref == (int)false) {
            AIM_ERROR(aimInfo, "capsReferenceChord is not set on any body!");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
        if (foundBref == (int)false) {
            AIM_ERROR(aimInfo, "capsReferenceSpan is not set on any body!");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
        if (foundXref == (int)false) {
            AIM_ERROR(aimInfo, "capsReferenceX is not set on any body!");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // Check for moment reference overwrites
        if (aimInputs[inMoment_Center-1].nullVal ==  NotNull) {

            Xref = aimInputs[inMoment_Center-1].vals.reals[0];
            Yref = aimInputs[inMoment_Center-1].vals.reals[1];
            Zref = aimInputs[inMoment_Center-1].vals.reals[2];
        }

        fprintf(fp, "CAPS generated Configuration\n");
        fprintf(fp, "0.0         # Mach\n");                                   /* Mach */
        fprintf(fp, "0 0 0       # IYsym   IZsym   Zsym\n");                   /* IYsym   IZsym   Zsym */
        fprintf(fp, "%lf %lf %lf # Sref    Cref    Bref\n", Sref, Cref, Bref); /* Sref    Cref    Bref */
        fprintf(fp, "%lf %lf %lf # Xref    Yref    Zref\n", Xref, Yref, Zref); /* Xref    Yref    Zref */
        fprintf(fp, "%lf         # CDp\n", aimInputs[inCDp-1].vals.real);      /* CDp */

        // Write out the Surfaces, one at a time
        for (surf = 0; surf < avlInstance->numAVLSurface; surf++) {

            printf("Writing surface - %s (ID = %d)\n", avlInstance->avlSurface[surf].name, surf);

            if (avlInstance->avlSurface[surf].numSection < 2) {
                printf("Surface %s only has %d Sections - it will be skipped!\n",
                       avlInstance->avlSurface[surf].name, avlInstance->avlSurface[surf].numSection);
                continue;
            }

            fprintf(fp, "#\nSURFACE\n%s\n%d %lf\n\n", avlInstance->avlSurface[surf].name,
                                                      avlInstance->avlSurface[surf].Nchord,
                                                      avlInstance->avlSurface[surf].Cspace);

            if  (avlInstance->avlSurface[surf].compon != 0)
                fprintf(fp, "COMPONENT\n%d\n\n", avlInstance->avlSurface[surf].compon);

            if  (avlInstance->avlSurface[surf].iYdup  != 0)
                fprintf(fp, "YDUPLICATE\n0.0\n\n");

            if  (avlInstance->avlSurface[surf].nowake == (int) true) fprintf(fp, "NOWAKE\n");
            if  (avlInstance->avlSurface[surf].noalbe == (int) true) fprintf(fp, "NOALBE\n");
            if  (avlInstance->avlSurface[surf].noload == (int) true) fprintf(fp, "NOLOAD\n");

            if ( (avlInstance->avlSurface[surf].nowake == (int) true) ||
                 (avlInstance->avlSurface[surf].noalbe == (int) true) ||
                 (avlInstance->avlSurface[surf].noload == (int) true)) {

                fprintf(fp,"\n");
            }

            // Write the sections for each surface
            for (i = 0; i < avlInstance->avlSurface[surf].numSection; i++) {

                section = avlInstance->avlSurface[surf].vlmSection[i].sectionIndex;

                printf("\tSection %d of %d (ID = %d)\n",
                       i+1, avlInstance->avlSurface[surf].numSection, section);

                // Write section data
                status = writeSection(aimInfo, fp, &avlInstance->avlSurface[surf].vlmSection[section]);
                AIM_STATUS(aimInfo, status);

                // Write control information for each section
                for (control = 0; control < avlInstance->avlSurface[surf].vlmSection[section].numControl; control++) {

                    printf("\t  Control surface %d of %d \n", control + 1, avlInstance->avlSurface[surf].vlmSection[section].numControl);

                    fprintf(fp, "CONTROL\n");
                    fprintf(fp, "%s ", avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].name);

                    fprintf(fp, "%f ", avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].controlGain);

                    if (avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].leOrTe == 0) { // Leading edge (-)

                        fprintf(fp, "%f ", -avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].percentChord);

                    } else { // Trailing edge (+)

                        fprintf(fp, "%f ", avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].percentChord);
                    }

                    fprintf(fp, "%f %f %f ",
                            avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[0],
                            avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[1],
                            avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].xyzHingeVec[2]);

                    fprintf(fp, "%f\n", (double) avlInstance->avlSurface[surf].vlmSection[section].vlmControl[control].deflectionDup);
                }
                fprintf(fp, "\n");
            }
        }
    }

    // write mass data file for needed for eigen value analysis
    if (eigenValues == (int) true) {
/*@-nullpass@*/
      status = writeMassFile(aimInfo, aimInputs, &avlInstance->units, bodyLunits, massFilename);
/*@+nullpass@*/
      AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:

    if (fp != NULL) fclose(fp);

    // Free array of control name strings
    if (numControlName != 0 && controlName != NULL) {
        (void) string_freeArray(numControlName, &controlName);
#ifdef S_SPLINT_S
        AIM_FREE(controlName);
#endif
    }

    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteAVL AIM Execution
   *
   * If auto execution is enabled when creating an AVL AIM,
   * the AIM will execute avl just-in-time with the command line:
   *
   * \code{.sh}
   * avl caps < avlInput.txt > avlOutput.txt
   * \endcode
   *
   * where preAnalysis generated the two files: 1) "avlInput.txt" which contains the input information and
   * control sequence for AVL to execute and 2) "caps.avl" which contains the geometry to be analyzed.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an AVL AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or avl can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * avl.preAnalysis()
   *
   * print ("\n\nRunning......")
   * avl.system("avl caps < avlInput.txt > avlOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * avl.postAnalysis()
   * \endcode
   */

  *state = 0;
  return aim_system(aimInfo, NULL, "avl caps < avlInput.txt > avlOutput.txt");
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  // check an AVL output file
  if (aim_isFile(aimInfo, "capsTotalForce.txt") != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "avl execution did not produce capsTotalForce.txt");
    return CAPS_EXECERR;
  }

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
    int status = CAPS_SUCCESS;

    cfdUnitsStruct *units=NULL;
    aimStorage     *avlInstance;

    avlInstance = (aimStorage *) instStore;
    AIM_NOTNULL(avlInstance, aimInfo, status);

    units = &avlInstance->units;

    /*! \page aimOutputsAVL AIM Outputs
     * Optional outputs that echo the inputs.  These are parsed from the resulting output and can be used as a sanity check.
     */

    if (index >= 90) {
        form->type       = Tuple;
        form->dim        = Vector;
        form->vals.tuple = NULL;
        form->nullVal    = IsNull;
        form->lfixed     = form->sfixed = Change;
    } else {
        form->type = Double;
        form->vals.real = 0;
    }

#ifdef DEBUG
    printf(" avlAIM/aimOutputs instance = %d  index = %d!\n", inst, index);
#endif

    // Echo AVL flow conditions
    if (index == Alpha) {
        *aoname = EG_strdup("Alpha");

        /*! \page aimOutputsAVL
         * - <B> Alpha </B> = Angle of attack.
         */
    } else if (index == Beta) {
        *aoname = EG_strdup("Beta");

        /*! \page aimOutputsAVL
         * - <B> Beta </B> = Sideslip angle.
         */
    } else if (index == Mach) {
        *aoname = EG_strdup("Mach");

        /*! \page aimOutputsAVL
         * - <B> Mach </B> = Mach number.
         */

    } else if (index == pbd2V) {
        *aoname = EG_strdup("pb/2V");

        /*! \page aimOutputsAVL
         * - <B> pb/2V </B> = Non-dimensional roll rate.
         */
    } else if (index == qcd2V) {
        *aoname = EG_strdup("qc/2V");

        /*! \page aimOutputsAVL
         * - <B> qc/2V </B> = Non-dimensional pitch rate.
         */
    } else if (index == rbd2V) {
        *aoname = EG_strdup("rb/2V");

        /*! \page aimOutputsAVL
         * - <B> rb/2V </B> = Non-dimensional yaw rate.
         */

    } else if (index == pPbd2V) {
        *aoname = EG_strdup("p'b/2V");

        /*! \page aimOutputsAVL
         * - <B> p'b/2V </B> = Non-dimensional roll acceleration.
         */

    } else if (index == rPbd2V) {
        *aoname = EG_strdup("r'b/2V");

        /*! \page aimOutputsAVL
         * - <B> r'b/2V </B> = Non-dimensional yaw acceleration.
         */

        // Body control derivatives
    } else if (index == CXtot) {
        *aoname = EG_strdup("CXtot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * Forces and moments:
         * - <B> CXtot </B> = X-component of total force in body axis
         */

    } else if (index == CYtot) {
        *aoname = EG_strdup("CYtot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> CYtot </B> = Y-component of total force in body axis
         */
    } else if (index == CZtot) {
        *aoname = EG_strdup("CZtot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> CZtot </B> = Z-component of total force in body axis
         */
    } else if (index == Cltot) {
        *aoname = EG_strdup("Cltot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> Cltot </B> = X-component of moment in body axis
         */
    } else if (index == Cmtot) {
        *aoname = EG_strdup("Cmtot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> Cmtot </B> = Y-component of moment in body axis
         */
    } else if (index == Cntot) {
        *aoname = EG_strdup("Cntot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> Cntot </B> = Z-component of moment in body axis
         */
    } else if (index == ClPtot) {
        *aoname = EG_strdup("Cl'tot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> Cl'tot </B> = x-component of moment in stability axis
         */
    } else if (index == CnPtot) {
        *aoname = EG_strdup("Cn'tot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> Cn'tot </B> = z-component of moment in stability axis
         */
    } else if (index == CLtot) {
        *aoname = EG_strdup("CLtot");
        form->type = DoubleDeriv;

        /*! \page aimOutputsAVL
         * - <B> CLtot </B> = total lift in stability axis
         */
    } else if (index == CDtot) {
        *aoname = EG_strdup("CDtot");

        /*! \page aimOutputsAVL
         * - <B> CDtot </B> = total drag in stability axis
         */
    } else if (index == CDvis) {
        *aoname = EG_strdup("CDvis");

        /*! \page aimOutputsAVL
         * - <B> CDvis </B> = viscous drag component
         */
    } else if (index == CLff) {
        *aoname = EG_strdup("CLff");

        /*! \page aimOutputsAVL
         * - <B> CLff </B> = trefftz plane lift force
         */
    } else if (index == CYff) {
        *aoname = EG_strdup("CYff");

        /*! \page aimOutputsAVL
         * - <B> CYff </B> = trefftz plane side force
         */
    } else if (index == CDind) {
        *aoname = EG_strdup("CDind");

        /*! \page aimOutputsAVL
         * - <B> CDind </B> = induced drag force
         */
    } else if (index == CDff) {
        *aoname = EG_strdup("CDff");

        /*! \page aimOutputsAVL
         * - <B> CDff </B> = trefftz plane drag force
         */
    } else if (index == e) {
        *aoname = EG_strdup("e");

        /*! \page aimOutputsAVL
         * - <B> e = </B> Oswald Efficiency
         */
    } else if (index == CLa) { // Alpha stability derivatives
        *aoname = EG_strdup("CLa");

    } else if (index == CYa) {
        *aoname = EG_strdup("CYa");

    } else if (index == ClPa) {
        *aoname = EG_strdup("Cl'a");

    } else if (index == Cma) {
        *aoname = EG_strdup("Cma");

    } else if (index == CnPa) {
        *aoname = EG_strdup("Cn'a");

        /*! \page aimOutputsAVL
         * Stability-axis derivatives - Alpha:
         * - <B>CLa</B> = z' force, CL, with respect to alpha.
         * - <B>CYa</B> = y force, CY, with respect to alpha.
         * - <B>Cl'a</B> = x' moment, Cl', with respect to alpha.
         * - <B>Cma</B> = y moment, Cm, with respect to alpha.
         * - <B>Cn'a</B> = z' moment, Cn', with respect to alpha.
         */

    } else if (index == CLb) { // Beta stability derivatives
        *aoname = EG_strdup("CLb");

    } else if (index == CYb) {
        *aoname = EG_strdup("CYb");

    } else if (index == ClPb) {
        *aoname = EG_strdup("Cl'b");

    } else if (index == Cmb) {
        *aoname = EG_strdup("Cmb");

    } else if (index == CnPb) {
        *aoname = EG_strdup("Cn'b");

        /*! \page aimOutputsAVL
         * Stability-axis derivatives - Beta:
         * - <B>CLb</B> = z' force, CL, with respect to beta.
         * - <B>CYb</B> = y force, CY, with respect to beta.
         * - <B>Cl'b</B> = x' moment, Cl', with respect to beta.
         * - <B>Cmb</B> = y moment, Cm, with respect to beta.
         * - <B>Cn'b</B> = z' moment, Cn', with respect to beta.
         */
    } else if (index == CLpP) { // Roll rate stability derivatives
        *aoname = EG_strdup("CLp'");

    } else if (index == CYpP) {
        *aoname = EG_strdup("CYp'");

    } else if (index == ClPpP) {
        *aoname = EG_strdup("Cl'p'");

    } else if (index == CmpP) {
        *aoname = EG_strdup("Cmp'");

    } else if (index == CnPpP) {
        *aoname = EG_strdup("Cn'p'");

        /*! \page aimOutputsAVL
         * Stability-axis derivatives - Roll rate, p':
         * - <B>CLp'</B> = z' force, CL, with respect to roll rate, p'.
         * - <B>CYp'</B> = y force, CY, with respect to roll rate, p'.
         * - <B>Cl'p'</B> = x' moment, Cl', with respect to roll rate, p'.
         * - <B>Cmp'</B> = y moment, Cm, with respect to roll rate, p'.
         * - <B>Cn'p'</B> = z' moment, Cn', with respect to roll rate, p'.
         */
    } else if (index == CLqP) { // Pitch rate stability derivatives
        *aoname = EG_strdup("CLq'");

    } else if (index == CYqP) {
        *aoname = EG_strdup("CYq'");

    } else if (index == ClPqP) {
        *aoname = EG_strdup("Cl'q'");

    } else if (index == CmqP) {
        *aoname = EG_strdup("Cmq'");

    } else if (index == CnPqP) {
        *aoname = EG_strdup("Cn'q'");

        /*! \page aimOutputsAVL
         * Stability-axis derivatives - Pitch rate, q':
         * - <B>CLq'</B> = z' force, CL, with respect to pitch rate, q'.
         * - <B>CYq'</B> = y force, CY, with respect to pitch rate, q'.
         * - <B>Cl'q'</B> = x' moment, Cl', with respect to pitch rate, q'.
         * - <B>Cmq'</B> = y moment, Cm, with respect to pitch rate, q'.
         * - <B>Cn'q'</B> = z' moment, Cn', with respect to pitch rate, q'.
         */
    } else if (index == CLrP) { // Yaw rate stability derivatives
        *aoname = EG_strdup("CLr'");

    } else if (index == CYrP) {
        *aoname = EG_strdup("CYr'");

    } else if (index == ClPrP) {
        *aoname = EG_strdup("Cl'r'");

    } else if (index == CmrP) {
        *aoname = EG_strdup("Cmr'");

    } else if (index == CnPrP) {
        *aoname = EG_strdup("Cn'r'");

        /*! \page aimOutputsAVL
         * Stability-axis derivatives - Yaw rate, r':
         * - <B>CLr'</B> = z' force, CL, with respect to yaw rate, r'.
         * - <B>CYr'</B> = y force, CY, with respect to yaw rate, r'.
         * - <B>Cl'r'</B> = x' moment, Cl', with respect to yaw rate, r'.
         * - <B>Cmr'</B> = y moment, Cm, with respect to yaw rate, r'.
         * - <B>Cn'r'</B> = z' moment, Cn', with respect to yaw rate, r'.
         */
    } else if (index == CXu) { // Axial velocity stability derivatives
        *aoname = EG_strdup("CXu");

    } else if (index == CYu) {
        *aoname = EG_strdup("CYu");

    } else if (index == CZu) {
        *aoname = EG_strdup("CZu");

    } else if (index == Clu) {
        *aoname = EG_strdup("Clu");

    } else if (index == Cmu) {
        *aoname = EG_strdup("Cmu");

    } else if (index == Cnu) {
        *aoname = EG_strdup("Cnu");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Axial velocity, u:
         * - <B>CXu</B> = x force, CX, with respect to axial velocity, u.
         * - <B>CYu</B> = y force, CY, with respect to axial velocity, u.
         * - <B>CZu</B> = z force, CZ, with respect to axial velocity, u.
         * - <B>Clu</B> = x moment, Cl, with respect to axial velocity, u.
         * - <B>Cmu</B> = y moment, Cm, with respect to axial velocity, u.
         * - <B>Cnu</B> = z moment, Cn, with respect to axial velocity, u.
         */
    } else if (index == CXv) { // Sideslip velocity stability derivatives
        *aoname = EG_strdup("CXv");

    } else if (index == CYv) {
        *aoname = EG_strdup("CYv");

    } else if (index == CZv) {
        *aoname = EG_strdup("CZv");

    } else if (index == Clv) {
        *aoname = EG_strdup("Clv");

    } else if (index == Cmv) {
        *aoname = EG_strdup("Cmv");

    } else if (index == Cnv) {
        *aoname = EG_strdup("Cnv");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Sideslip velocity, v:
         * - <B>CXv</B> = x force, CX, with respect to sideslip velocity, v.
         * - <B>CYv</B> = y force, CY, with respect to sideslip velocity, v.
         * - <B>CZv</B> = z force, CZ, with respect to sideslip velocity, v.
         * - <B>Clv</B> = x moment, Cl, with respect to sideslip velocity, v.
         * - <B>Cmv</B> = y moment, Cm, with respect to sideslip velocity, v.
         * - <B>Cnv</B> = z moment, Cn, with respect to sideslip velocity, v.
         */
    } else if (index == CXw) { // Normal velocity stability derivatives
        *aoname = EG_strdup("CXw");
    } else if (index == CYw) {
        *aoname = EG_strdup("CYw");
    } else if (index == CZw) {
        *aoname = EG_strdup("CZw");
    } else if (index == Clw) {
        *aoname = EG_strdup("Clw");
    } else if (index == Cmw) {
        *aoname = EG_strdup("Cmw");
    } else if (index == Cnw) {
        *aoname = EG_strdup("Cnw");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Normal velocity, w:
         * - <B>CXw</B> = x force, CX, with respect to normal velocity, w.
         * - <B>CYw</B> = y force, CY, with respect to normal velocity, w.
         * - <B>CZw</B> = z force, CZ, with respect to normal velocity, w.
         * - <B>Clw</B> = x moment, Cl, with respect to normal velocity, w.
         * - <B>Cmw</B> = y moment, Cm, with respect to normal velocity, w.
         * - <B>Cnw</B> = z moment, Cn, with respect to normal velocity, w.
         */
    } else if (index == CXp) { // Roll rate stability derivatives
        *aoname = EG_strdup("CXp");
    } else if (index == CYp) {
        *aoname = EG_strdup("CYp");
    } else if (index == CZp) {
        *aoname = EG_strdup("CZp");
    } else if (index == Clp) {
        *aoname = EG_strdup("Clp");
    } else if (index == Cmp) {
        *aoname = EG_strdup("Cmp");
    } else if (index == Cnp) {
        *aoname = EG_strdup("Cnp");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Roll rate, p:
         * - <B>CXp</B> = x force, CX, with respect to roll rate, p.
         * - <B>CYp</B> = y force, CY, with respect to roll rate, p.
         * - <B>CZp</B> = z force, CZ, with respect to roll rate, p.
         * - <B>Clp</B> = x moment, Cl, with respect to roll rate, p.
         * - <B>Cmp</B> = y moment, Cm, with respect to roll rate, p.
         * - <B>Cnp</B> = z moment, Cn, with respect to roll rate, p.
         */
    } else if (index == CXq) { // Pitch rate stability derivatives
        *aoname = EG_strdup("CXq");
    } else if (index == CYq) {
        *aoname = EG_strdup("CYq");
    } else if (index == CZq) {
        *aoname = EG_strdup("CZq");
    } else if (index == Clq) {
        *aoname = EG_strdup("Clq");
    } else if (index == Cmq) {
        *aoname = EG_strdup("Cmq");
    } else if (index == Cnq) {
        *aoname = EG_strdup("Cnq");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Pitch rate, q:
         * - <B>CXq</B> = x force, CX, with respect to pitch rate, q.
         * - <B>CYq</B> = y force, CY, with respect to pitch rate, q.
         * - <B>CZq</B> = z force, CZ, with respect to pitch rate, q.
         * - <B>Clq</B> = x moment, Cl, with respect to pitch rate, q.
         * - <B>Cmq</B> = y moment, Cm, with respect to pitch rate, q.
         * - <B>Cnq</B> = z moment, Cn, with respect to pitch rate, q.
         */
    } else if (index == CXr) { // Yaw rate stability derivatives
        *aoname = EG_strdup("CXr");
    } else if (index == CYr) {
        *aoname = EG_strdup("CYr");
    } else if (index == CZr) {
        *aoname = EG_strdup("CZr");
    } else if (index == Clr) {
        *aoname = EG_strdup("Clr");
    } else if (index == Cmr) {
        *aoname = EG_strdup("Cmr");
    } else if (index == Cnr) {
        *aoname = EG_strdup("Cnr");

        /*! \page aimOutputsAVL
         * Body-axis derivatives - Yaw rate, r:
         * - <B>CXr</B> = x force, CX, with respect to yaw rate, r.
         * - <B>CYr</B> = y force, CY, with respect to yaw rate, r.
         * - <B>CZr</B> = z force, CZ, with respect to yaw rate, r.
         * - <B>Clr</B> = x moment, Cl, with respect to yaw rate, r.
         * - <B>Cmr</B> = y moment, Cm, with respect to yaw rate, r.
         * - <B>Cnr</B> = z moment, Cn, with respect to yaw rate, r.
         */

    } else if (index == Xnp) {
        *aoname = EG_strdup("Xnp");
        if (units->length != NULL) {
          AIM_STRDUP(form->units, units->length, aimInfo, status);
        }
    } else if (index == Xcg) {
        *aoname = EG_strdup("Xcg");
        if (units->length != NULL) {
          AIM_STRDUP(form->units, units->length, aimInfo, status);
        }
    } else if (index == Ycg) {
        *aoname = EG_strdup("Ycg");
        if (units->length != NULL) {
          AIM_STRDUP(form->units, units->length, aimInfo, status);
        }
    } else if (index == Zcg) {
        *aoname = EG_strdup("Zcg");
        if (units->length != NULL) {
          AIM_STRDUP(form->units, units->length, aimInfo, status);
        }

        /*! \page aimOutputsAVL
         * Geometric output:
         * - <B>Xnp</B> = Neutral Point
         * - <B>Xcg</B> = x CG location
         * - <B>Ycg</B> = y CG location
         * - <B>Zcg</B> = z CG location
         */

    } else if (index == ControlStability) {
        *aoname = EG_strdup("ControlStability");

        /*! \page aimOutputsAVL
         * Controls:
         * - <B>ControlStability</B> = a (or an array of) tuple(s) with a
         * structure of ("Control Surface Name", "JSON Dictionary") for all control surfaces in the stability axis frame.
         * The JSON dictionary has the form = {"CLtot":value,"CYtot":value,"Cl'tot":value,"Cmtot":value,"Cn'tot":value}
         */
    } else if (index == ControlBody) {
        *aoname = EG_strdup("ControlBody");

        /*! \page aimOutputsAVL
         * - <B>ControlBody</B> = a (or an array of) tuple(s) with a
         * structure of ("Control Surface Name", "JSON Dictionary") for all control surfaces in the body axis frame.
         * The JSON dictionary has the form = {"CXtot":value,"CYtot":value,"CZtot":value,"Cltot":value,"Cmtot":value,"Cntot":value}
         */
    } else if (index == HingeMoment) {
        *aoname = EG_strdup("HingeMoment");

        /*! \page aimOutputsAVL
         * - <B>HingeMoment</B> = a (or an array of) tuple(s) with a
         * structure of ("Control Surface Name", "HingeMoment")
         */
    } else if (index == StripForces) {
        *aoname = EG_strdup("StripForces");

        /*! \page aimOutputsAVL
         * - <B>StripForces</B> = a (or an array of) tuple(s) with a
         * structure of ("Surface Name", "JSON Dictionary") for all surfaces.
         * The JSON dictionary has the form = {"cl":[value0,value1,value2],"cd":[value0,value1,value2]...}
         */
    } else if (index == EigenValues) {
        *aoname = EG_strdup("EigenValues");

        /*! \page aimOutputsAVL
         * - <B>EigenValues</B> = a (or an array of) tuple(s) with a
         * structure of ("case #", "Array of eigen values").
         * The array of eigen values is of the form = [[real0,imaginary0],[real0,imaginary0],...]
         */

    } else {

        return CAPS_NOTFOUND;
    }

cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{

    int status; // Function return status

    int i; // Indexing

    const char *key = NULL;
    const char *fileToOpen = "capsTotalForce.txt";
    const char *derivFile = NULL;
    const char *deriv_keys[5] = {NULL,NULL,NULL,NULL,NULL};
    int nderiv = 0;
    aimStorage *avlInstance;

    double tempVal[6];

    char jsonOut[200];

    /*
    fileToOpen= "capsTotalForce.txt";
    fileToOpen = "capsStatbilityDeriv.txt";
    fileToOpen= "capsBodyAxisDeriv.txt";
    fileToOpen = "capsHingeMoment.txt";
     */

#ifdef DEBUG
    const char *name;

    status = aim_getName(aimInfo, index, ANALYSISOUT, &name);
    printf(" avlAIM/aimCalcOutput index =  %s %d!\n", index, name, status);
#endif
    avlInstance = (aimStorage *) instStore;

    val->vals.real = 0.0;
/*@-observertrans@*/
    switch (index) {
    case Alpha:
        fileToOpen= "capsTotalForce.txt";
        key = "Alpha =";
        break;
    case Beta:
        fileToOpen= "capsTotalForce.txt";
        key = "Beta  =";
        break;
    case Mach:
        fileToOpen= "capsTotalForce.txt";
        key = "Mach  =";
        break;
    case pbd2V:
        fileToOpen= "capsTotalForce.txt";
        key = "pb/2V =";
        break;
    case qcd2V:
        fileToOpen= "capsTotalForce.txt";
        key = "qc/2V =";
        break;
    case rbd2V:
        fileToOpen= "capsTotalForce.txt";
        key = "rb/2V =";
        break;
    case pPbd2V:
        fileToOpen= "capsTotalForce.txt";
        key = "p'b/2V =";
        break;
    case rPbd2V:
        fileToOpen= "capsTotalForce.txt";
        key = "r'b/2V =";
        break;
    case CXtot:
        fileToOpen= "capsTotalForce.txt";
        key = "CXtot =";

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile= "capsBodyAxisDeriv.txt";
        nderiv = 3;
        deriv_keys[0] = "CXp =";
        deriv_keys[1] = "CXq =";
        deriv_keys[2] = "CXr =";

        status = controlDerivatives(aimInfo, CXtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case CYtot:
        fileToOpen= "capsTotalForce.txt";
        key = "CYtot =";

        status = stabilityAngleDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        deriv_keys[0] = "CYa =";
        deriv_keys[1] = "CYb =";
        for (i = 0; i < 2; i++) {
          status = read_Data(aimInfo, "capsStatbilityDeriv.txt", deriv_keys[i], val->derivs[i].deriv);
          AIM_STATUS(aimInfo, status);
        }

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        deriv_keys[0] = "CYp =";
        deriv_keys[1] = "CYq =";
        deriv_keys[2] = "CYr =";
        for (i = 0; i < 3; i++) {
          status = read_Data(aimInfo, "capsBodyAxisDeriv.txt", deriv_keys[i], val->derivs[i+2].deriv);
          AIM_STATUS(aimInfo, status);
        }

        status = controlDerivatives(aimInfo, CYtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case CZtot:
        fileToOpen= "capsTotalForce.txt";
        key = "CZtot =";

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile= "capsBodyAxisDeriv.txt";
        nderiv = 3;
        deriv_keys[0] = "CZp =";
        deriv_keys[1] = "CZq =";
        deriv_keys[2] = "CZr =";

        status = controlDerivatives(aimInfo, CZtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case Cltot:
        fileToOpen= "capsTotalForce.txt";
        key = "Cltot =";

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile= "capsBodyAxisDeriv.txt";
        nderiv = 3;
        deriv_keys[0] = "Clp =";
        deriv_keys[1] = "Clq =";
        deriv_keys[2] = "Clr =";

        status = controlDerivatives(aimInfo, Cltot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case Cmtot:
        fileToOpen= "capsTotalForce.txt";
        key = "Cmtot =";

        status = stabilityAngleDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        deriv_keys[0] = "Cma =";
        deriv_keys[1] = "Cmb =";
        for (i = 0; i < 2; i++) {
          status = read_Data(aimInfo, "capsStatbilityDeriv.txt", deriv_keys[i], val->derivs[i].deriv);
          AIM_STATUS(aimInfo, status);
        }

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        deriv_keys[0] = "Cmp =";
        deriv_keys[1] = "Cmq =";
        deriv_keys[2] = "Cmr =";
        for (i = 0; i < 3; i++) {
          status = read_Data(aimInfo, "capsBodyAxisDeriv.txt", deriv_keys[i], val->derivs[i+2].deriv);
          AIM_STATUS(aimInfo, status);
        }

        status = controlDerivatives(aimInfo, Cmtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case Cntot:
        fileToOpen= "capsTotalForce.txt";
        key = "Cntot =";

        status = bodyRateDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile= "capsBodyAxisDeriv.txt";
        nderiv = 3;
        deriv_keys[0] = "Cnp =";
        deriv_keys[1] = "Cnq =";
        deriv_keys[2] = "Cnr =";

        status = controlDerivatives(aimInfo, Cntot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case ClPtot:
        fileToOpen= "capsTotalForce.txt";
        key = "Cl'tot =";

        status = stabilityAngleDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile = "capsStatbilityDeriv.txt";
        nderiv = 2;
        deriv_keys[0] = "Cla =";
        deriv_keys[1] = "Clb =";

        status = controlDerivatives(aimInfo, ClPtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case CnPtot:
        fileToOpen= "capsTotalForce.txt";
        key = "Cn'tot =";

        status = stabilityAngleDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile = "capsStatbilityDeriv.txt";
        nderiv = 2;
        deriv_keys[0] = "Cna =";
        deriv_keys[1] = "Cnb =";

        status = controlDerivatives(aimInfo, CnPtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case CLtot:
        fileToOpen= "capsTotalForce.txt";
        key = "CLtot =";

        status = stabilityAngleDerivatives(aimInfo, val);
        AIM_STATUS(aimInfo, status);

        derivFile = "capsStatbilityDeriv.txt";
        nderiv = 2;
        deriv_keys[0] = "CLa =";
        deriv_keys[1] = "CLb =";

        status = controlDerivatives(aimInfo, CLtot, avlInstance, val);
        AIM_STATUS(aimInfo, status);

        break;
    case CDtot:
        fileToOpen= "capsTotalForce.txt";
        key = "CDtot =";
        break;
    case CDvis:
        fileToOpen= "capsTotalForce.txt";
        key = "CDvis =";
        break;
    case CLff:
        fileToOpen= "capsTotalForce.txt";
        key = "CLff  =";
        break;
    case CYff:
        fileToOpen= "capsTotalForce.txt";
        key = "CYff  =";
        break;
    case CDind:
        fileToOpen= "capsTotalForce.txt";
        key = "CDind =";
        break;
    case CDff:
        fileToOpen= "capsTotalForce.txt";
        key = "CDff  =";
        break;
    case e:
        fileToOpen= "capsTotalForce.txt";
        key = "e =";
        break;

        // Alpha stability derivatives
    case CLa:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CLa =";
        break;

    case CYa:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CYa =";
        break;

    case ClPa:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cla =";
        break;

    case Cma:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cma =";
        break;

    case CnPa:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cna =";
        break;

        // Beta stability derivatives
    case CLb:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CLb =";
        break;

    case CYb:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CYb =";
        break;

    case ClPb:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Clb =";
        break;

    case Cmb:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cmb =";
        break;

    case CnPb:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cnb =";
        break;

        // Roll rate p' stability derivatives
    case CLpP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CLp =";
        break;

    case CYpP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CYp =";
        break;

    case ClPpP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Clp =";
        break;

    case CmpP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cmp =";
        break;

    case CnPpP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cnp =";
        break;

        // Pitch rate q' stability derivatives
    case CLqP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CLq =";
        break;

    case CYqP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CYq =";
        break;

    case ClPqP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Clq =";
        break;

    case CmqP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cmq =";
        break;

    case CnPqP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cnq =";
        break;

        // Yaw rate r' stability derivatives
    case CLrP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CLr =";
        break;

    case CYrP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "CYr =";
        break;

    case ClPrP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Clr =";
        break;

    case CmrP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cmr =";
        break;

    case CnPrP:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Cnr =";
        break;


        // Axial vel (body axis)  stability derivatives
    case CXu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXu =";
        break;

    case CYu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYu =";
        break;

    case CZu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZu =";
        break;

    case Clu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clu =";
        break;

    case Cmu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmu =";
        break;

    case Cnu:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnu =";
        break;


        // Slidslip vel (body axis) stability derivatives
    case CXv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXv =";
        break;

    case CYv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYv =";
        break;

    case CZv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZv =";
        break;

    case Clv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clv =";
        break;

    case Cmv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmv =";
        break;

    case Cnv:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnv =";
        break;

        // Normal vel (body axis)  stability derivatives
    case CXw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXw =";
        break;

    case CYw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYw =";
        break;

    case CZw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZw =";
        break;

    case Clw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clw =";
        break;

    case Cmw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmw =";
        break;

    case Cnw:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnw =";
        break;

        // Roll rate (body axis)  stability derivatives
    case CXp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXp =";
        break;

    case CYp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYp =";
        break;

    case CZp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZp =";
        break;

    case Clp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clp =";
        break;

    case Cmp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmp =";
        break;

    case Cnp:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnp =";
        break;

        // Pitch rate (body axis)  stability derivatives
    case CXq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXq =";
        break;

    case CYq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYq =";
        break;

    case CZq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZq =";
        break;

    case Clq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clq =";
        break;

    case Cmq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmq =";
        break;

    case Cnq:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnq =";
        break;

        // Yaw rate (body axis)  stability derivatives
    case CXr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CXr =";
        break;

    case CYr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CYr =";
        break;

    case CZr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "CZr =";
        break;

    case Clr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Clr =";
        break;

    case Cmr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cmr =";
        break;

    case Cnr:
        fileToOpen= "capsBodyAxisDeriv.txt";
        key = "Cnr =";
        break;

    case Xnp:
        fileToOpen = "capsStatbilityDeriv.txt";
        key = "Xnp =";
        break;

    case Xcg:
        fileToOpen= "caps.run";
        key = "X_cg      =";
        break;

    case Ycg:
        fileToOpen= "caps.run";
        key = "Y_cg      =";
        break;

    case Zcg:
        fileToOpen= "caps.run";
        key = "Z_cg      =";
        break;

    case ControlStability:
    case ControlBody:
    case HingeMoment:
        fileToOpen = "capsHingeMoment.txt";
        break;

    case StripForces:
        fileToOpen = "capsStripForce.txt";
        break;

    case EigenValues:
        fileToOpen = "capsEigenValues.txt";
        break;
    }
/*@+observertrans@*/

    if (index < ControlStability) {

      if (key == NULL) {
          AIM_ERROR(aimInfo, "Developer error: No string key found!");
          status = CAPS_NOTFOUND;
          goto cleanup;
      }

      status = read_Data(aimInfo, fileToOpen, key, &val->vals.real);
      AIM_STATUS(aimInfo, status);

      /* read derivatives */
      for (i = 0; i < nderiv; i++) {
        /*@-nullpass@*/
        status = read_Data(aimInfo, derivFile, deriv_keys[i], val->derivs[i].deriv);
        AIM_STATUS(aimInfo, status);
        /*@+nullpass@*/
      }

    } else if (index >= ControlStability && index <= HingeMoment) { // Need to build something for control output

        // Initiate tuple base on number of control surfaces
        val->nrow = avlInstance->controlMap.numAttribute;

        // nothing to do if there is no control attributes
        if (val->nrow == 0) {
            val->nullVal = IsNull;
            return CAPS_SUCCESS;
        }

        AIM_ALLOC(val->vals.tuple, val->nrow, capsTuple, aimInfo, status);

        for (i = 0; i < val->nrow; i++) val->vals.tuple[i].name = val->vals.tuple[i].value = NULL;

        // Loop through control surfaces
        for (i = 0; i < avlInstance->controlMap.numAttribute; i++) {

            AIM_STRDUP(val->vals.tuple[i].name, avlInstance->controlMap.attributeName[i], aimInfo, status);

            // Stability axis
            if (index == ControlStability) {
                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CLtot, &tempVal[0]);
                AIM_STATUS(aimInfo, status);

                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CYtot, &tempVal[1]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          ClPtot, &tempVal[2]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          Cmtot, &tempVal[3]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CnPtot, &tempVal[4]);
                AIM_STATUS(aimInfo, status);

                sprintf(jsonOut,"{\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f}",
                        "CLtot",  tempVal[0],
                        "CYtot",  tempVal[1],
                        "Cl'tot", tempVal[2],
                        "Cmtot",  tempVal[3],
                        "Cn'tot", tempVal[4]);
                AIM_STRDUP(val->vals.tuple[i].value, jsonOut, aimInfo, status);

            } else if (index == ControlBody) {

                // Body axis
                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CXtot, &tempVal[0]);
                AIM_STATUS(aimInfo, status);

                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CYtot, &tempVal[1]);
                AIM_STATUS(aimInfo, status);

                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          CZtot, &tempVal[2]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          Cltot, &tempVal[3]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          Cmtot, &tempVal[4]);
                AIM_STATUS(aimInfo, status);


                status = get_controlDeriv(aimInfo, avlInstance->controlMap.attributeIndex[i],
                                          Cntot, &tempVal[5]);
                AIM_STATUS(aimInfo, status);;

                sprintf(jsonOut,"{\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f,\"%s\":%7.6f}",
                        "CXtot", tempVal[0], "CYtot", tempVal[1], "CZtot", tempVal[2],
                        "Cltot", tempVal[3], "Cmtot", tempVal[4], "Cntot", tempVal[5]);
                val->vals.tuple[i].value = EG_strdup(jsonOut);

            } else if (index == HingeMoment) {

                status = read_Data(aimInfo, fileToOpen, val->vals.tuple[i].name, &tempVal[0]);
                AIM_STATUS(aimInfo, status);

                sprintf(jsonOut, "%5.4e", tempVal[0]);

                AIM_STRDUP(val->vals.tuple[i].value, jsonOut, aimInfo, status);

            } else {
                status = CAPS_NOTFOUND;
                goto cleanup;
            }

        }

    } else if (index == StripForces) {

      // Read in the strip forces
      status = read_StripForces(aimInfo, &val->nrow, &val->vals.tuple);
      AIM_STATUS(aimInfo, status);

    } else if (index == EigenValues) {

      // Read in the strip forces
      status = read_EigenValues(aimInfo, &val->nrow, &val->vals.tuple);
      AIM_STATUS(aimInfo, status);

    } else {

        printf("DEVELOPER Error! Unknown index %d\n", index);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    status = CAPS_SUCCESS;

    cleanup:
        return status;
}


// ********************** AIM Function Break *****************************
void aimCleanup(void *instStore)
{

    aimStorage *avlInstance;

#ifdef DEBUG
    printf(" avlAIM/aimCleanup!\n");
#endif
    avlInstance = (aimStorage *) instStore;

    // cleanup storage
    destroy_aimStorage(avlInstance);
    AIM_FREE(avlInstance);
}

// ********************** AIM Function Break *****************************
// Back function used to calculate sensitivities
int aimBackdoor(void *instStore, void *aimInfo, const char *JSONin,
                char **JSONout)
{

    /*! \page aimBackDoorAVL AIM Back Door
     * The back door function of this AIM may be used as an alternative to retrieve sensitivity information produced by the
     * AIM. The JSONin string should be of the following form '{"mode": "sensitivity", "inputVar": "Name of Input Variable",
     * "outputVar": "Name of Output Variable"}', while the JSONout string will look like '{"sensitivity": value}'. Important: the
     * JSONout string is freeable! Invalid combinations of input and output variables returns a CAPS_MISMATCH error code.
     *
     * Acceptable values for the "Name of Input Variable" are as follows (definitions are consistent, where appropriate,
     *  with \ref aimInputsAVL):
     * - "Alpha"
     * - "Beta"
     * - "RollRate"
     * - "PitchRate"
     * - "YawRate"
     * - "AxialVelocity"
     * - "SideslipVelocity"
     * - "NormalVelocity"
     * - "AVL_Control:Name_of_Control_Surface", where Name_of_Control_Surface should be replaced with name of the desired
     * control surface
     *
     * Acceptable values for the "Name of Output Variable" are as follows (definitions are consistent with \ref aimOutputsAVL):
     * - "CLtot"
     * - "CYtot"
     * - "Cl'tot"
     * - "Cmtot"
     * - "Cn'tot"
     * - "CXtot"
     * - "CZtot"
     * - "Cltot"
     * - "Cntot"
     */

    int status; // Function return status

    int index = -99;
    double data;
    capsValue val;
    aimStorage *avlInstance;

    char *keyValue = NULL;
    char *keyWord  = NULL;

    char *tempString = NULL; // Freeable

    int inputIndex = 0, outputIndex = 0, controlIndex = 0;

    char *outputJSON = NULL;

    *JSONout = NULL;

    keyWord  = "mode";
    status   = search_jsonDictionary(JSONin, keyWord, &keyValue);
    if (status != CAPS_SUCCESS) goto cleanup;
    if (keyValue == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    // Type of acceptable modes for this AIM backdoor
    if (strcasecmp(keyValue, "\"Sensitivity\"") != 0) {

        printf("Error: A valid mode wasn't found for AIMBackDoor!\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    avlInstance = (aimStorage *) instStore;

    // We are using this to get sensitivities
    if (strcasecmp(keyValue, "\"Sensitivity\"") == 0) {

        AIM_FREE(keyValue);
        keyValue = NULL;

        // Input variable
        keyWord = "inputVar";
        status = search_jsonDictionary(JSONin, keyWord, &keyValue);
        if (status != CAPS_SUCCESS) goto cleanup;
        if (keyValue == NULL) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

        tempString = string_removeQuotation(keyValue);
        AIM_NOTNULL(tempString, aimInfo, status);

        inputIndex = aim_getIndex(aimInfo, tempString, ANALYSISIN);
        if (inputIndex <= 0) {
            status = inputIndex;

            if (strcasecmp(tempString, "AxialVelocity") == 0) {
                inputIndex = NUMINPUT+1;

            } else if (strcasecmp(tempString, "SideslipVelocity") == 0) {
                inputIndex = NUMINPUT+2;

            } else if (strcasecmp(tempString, "NormalVelocity") == 0) {
                inputIndex = NUMINPUT+3;

                // Try to see if it is a control parameter
            } else if (parse_controlName(avlInstance, tempString,
                                         &controlIndex) == CAPS_SUCCESS) {
                inputIndex = CAPSMAGIC;

            } else {

                printf("Error: Unable to get index for inputVar = %s\n", tempString);
                goto cleanup;
            }
        }

        AIM_FREE(tempString);

        // Outout variable
        keyWord = "outputVar";
        status  = search_jsonDictionary(JSONin, keyWord, &keyValue);
        if (status != CAPS_SUCCESS) goto cleanup;

        tempString = string_removeQuotation(keyValue);
        AIM_NOTNULL(tempString, aimInfo, status);

        outputIndex = aim_getIndex(aimInfo, tempString, ANALYSISOUT);
        if (outputIndex <= 0) {
            status = outputIndex;
            printf("Error: Unable to get index for outputVar = %s\n", tempString);
            goto cleanup;
        }

        AIM_FREE(tempString);

        //printf("InputIndex = %d\n", inputIndex);
        //printf("OutputIndex = %d\n", outputIndex);

        // Stability axis - Alpha
        if (outputIndex == CLtot && inputIndex == inAlpha) {

            index = CLa;
        }

        else if (outputIndex == CYtot && inputIndex == inAlpha) {

            index = CYa;

        }

        else if (outputIndex == ClPtot && inputIndex == inAlpha) {

            index = ClPa;
        }

        else if (outputIndex == Cmtot && inputIndex == inAlpha) {

            index = Cma;
        }

        else if (outputIndex == CnPtot && inputIndex == inAlpha) {

            index = CnPa;
        }

        // Stability axis - Beta
        else if (outputIndex == CLtot && inputIndex == inBeta) {

            index = CLb;
        }

        else if (outputIndex == CYtot && inputIndex == inBeta) {

            index = CYb;

        }

        else if (outputIndex == ClPtot && inputIndex == inBeta) {

            index = ClPb;
        }

        else if (outputIndex == Cmtot && inputIndex == inBeta) {

            index = Cmb;
        }

        else if (outputIndex == CnPtot && inputIndex == inBeta) {

            index = CnPb;
        }

        // Stability axis - RollRate
        else if (outputIndex == CLtot && inputIndex == inRollRate) {

            index = CLpP;
        }

        else if (outputIndex == CYtot && inputIndex == inRollRate) {

            index = CYpP;

        }

        else if (outputIndex == ClPtot && inputIndex == inRollRate) {

            index = ClPpP;
        }

        else if (outputIndex == Cmtot && inputIndex == inRollRate) {

            index = CmpP;
        }

        else if (outputIndex == CnPtot && inputIndex == inRollRate) {

            index = CnPpP;
        }

        // Stability axis - PitchRate
        else if (outputIndex == CLtot && inputIndex == inPitchRate) {

            index = CLqP;
        }

        else if (outputIndex == CYtot && inputIndex == inPitchRate) {

            index = CYqP;

        }

        else if (outputIndex == ClPtot && inputIndex == inPitchRate) {

            index = ClPqP;
        }

        else if (outputIndex == Cmtot && inputIndex == inPitchRate) {

            index = CmqP;
        }

        else if (outputIndex == CnPtot && inputIndex == inPitchRate) {

            index = CnPqP;
        }


        // Stability axis - YawRate
        else if (outputIndex == CLtot && inputIndex == inYawRate) {

            index = CLrP;
        }

        else if (outputIndex == CYtot && inputIndex == inYawRate) {

            index = CYrP;

        }

        else if (outputIndex == ClPtot && inputIndex == inYawRate) {

            index = ClPrP;
        }

        else if (outputIndex == Cmtot && inputIndex == inYawRate) {

            index = CmrP;
        }

        else if (outputIndex == CnPtot && inputIndex == inYawRate) {

            index = CnPrP;
        }


        //////////////////////////////////////////////////////////

        // Body axis - AxialVelocity
        else if (outputIndex == CXtot && inputIndex == NUMINPUT+1) {

            index = CXu;
        }

        else if (outputIndex == CYtot && inputIndex == NUMINPUT+1) {

            index = CYu;

        }

        else if (outputIndex == CZtot && inputIndex == NUMINPUT+1) {

            index = CZu;

        }

        else if (outputIndex == Cltot && inputIndex == NUMINPUT+1) {

            index = Clu;
        }

        else if (outputIndex == Cmtot && inputIndex == NUMINPUT+1) {

            index = Cmu;
        }

        else if (outputIndex == Cntot && inputIndex == NUMINPUT+1) {

            index = Cnu;
        }

        // Body axis - Sideslip
        else if (outputIndex == CXtot && inputIndex == NUMINPUT+2) {

            index = CXv;
        }

        else if (outputIndex == CYtot && inputIndex == NUMINPUT+2) {

            index = CYv;

        }

        else if (outputIndex == CZtot && inputIndex == NUMINPUT+2) {

            index = CZv;

        }

        else if (outputIndex == Cltot && inputIndex == NUMINPUT+2) {

            index = Clv;
        }

        else if (outputIndex == Cmtot && inputIndex == NUMINPUT+2) {

            index = Cmv;
        }

        else if (outputIndex == Cntot && inputIndex == NUMINPUT+2) {

            index = Cnv;
        }

        // Body axis - Normal
        else if (outputIndex == CXtot && inputIndex == NUMINPUT+3) {

            index = CXw;
        }

        else if (outputIndex == CYtot && inputIndex == NUMINPUT+3) {

            index = CYw;

        }

        else if (outputIndex == CZtot && inputIndex == NUMINPUT+3) {

            index = CZw;

        }

        else if (outputIndex == Cltot && inputIndex == NUMINPUT+3) {

            index = Clw;
        }

        else if (outputIndex == Cmtot && inputIndex == NUMINPUT+3) {

            index = Cmw;
        }

        else if (outputIndex == Cntot && inputIndex == NUMINPUT+3) {

            index = Cnw;
        }

        // Body axis - RollRate
        else if (outputIndex == CXtot && inputIndex == inRollRate) {

            index = CXp;
        }

        else if (outputIndex == CYtot && inputIndex == inRollRate) {

            index = CYp;

        }

        else if (outputIndex == CZtot && inputIndex == inRollRate) {

            index = CZp;

        }

        else if (outputIndex == Cltot && inputIndex == inRollRate) {

            index = Clp;
        }

        else if (outputIndex == Cmtot && inputIndex == inRollRate) {

            index = Cmp;
        }

        else if (outputIndex == Cntot && inputIndex == inRollRate) {

            index = Cnp;
        }

        // Body axis - PitchRate
        else if (outputIndex == CXtot && inputIndex == inPitchRate) {

            index = CXq;
        }

        else if (outputIndex == CYtot && inputIndex == inPitchRate) {

            index = CYq;

        }

        else if (outputIndex == CZtot && inputIndex == inPitchRate) {

            index = CZq;

        }

        else if (outputIndex == Cltot && inputIndex == inPitchRate) {

            index = Clq;
        }

        else if (outputIndex == Cmtot && inputIndex == inPitchRate) {

            index = Cmq;
        }

        else if (outputIndex == Cntot && inputIndex == inPitchRate) {

            index = Cnq;
        }

        // Body axis - YawRate
        else if (outputIndex == CXtot && inputIndex == inYawRate) {

            index = CXr;
        }

        else if (outputIndex == CYtot && inputIndex == inYawRate) {

            index = CYr;

        }

        else if (outputIndex == CZtot && inputIndex == inYawRate) {

            index = CZr;

        }

        else if (outputIndex == Cltot && inputIndex == inYawRate) {

            index = Clr;
        }

        else if (outputIndex == Cmtot && inputIndex == inYawRate) {

            index = Cmr;
        }

        else if (outputIndex == Cntot && inputIndex == inYawRate) {

            index = Cnr;
        }
        else if (inputIndex == CAPSMAGIC) { // Control variable found
            index = CAPSMAGIC;
        }

        else {
            printf("Invalid combination of input and output variables.\n");
            status = CAPS_MISMATCH;
            goto cleanup;
        }

        // If we have a control variable
        if (index == CAPSMAGIC) {

            status = get_controlDeriv(aimInfo, controlIndex, outputIndex, &data);
            if (status != CAPS_SUCCESS) goto cleanup;

            val.vals.real = data;


        } else { // Make the standard call

            status = aimCalcOutput(instStore, aimInfo, index, &val);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        outputJSON = (char *) EG_alloc(50*sizeof(char));
        if (outputJSON == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        sprintf(outputJSON, "{\"sensitivity\": %7.6f}", val.vals.real);
    }

    *JSONout = outputJSON;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in aimBackdoor, status = %d\n", status);

    AIM_FREE(keyValue);
    AIM_FREE(tempString);

    return status;
}
