/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Awave AIM
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>

#include "aimUtil.h"
#include "miscUtils.h"


#define PI        3.1415926535897931159979635
#define NINT(A)         (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))

//#define DEBUG

enum aimInputs
{
  inMach = 1,                  /* index is 1-based */
  inAlpha,
  NUMINPUT = inAlpha           /* Total number of inputs */
};

enum aimOutputs
{
  outCDwave = 1,               /* index is 1-based */
  outMach,
  outAlpha,
  NUMOUTPUT = outAlpha         /* Total number of outputs */
};


/*!\mainpage Introduction
 *
 * \tableofcontents
 *
 * \section overviewAwave Awave AIM Overview
 *
 * Awave provides an estimation for wave drag at supersonic Mach numbers at various angles of attack.
 * Taken from the Awave manual \cite AWAVE :
 *
 * "Awave is a streamlined, modified version of the Harris far-field wave drag program
 * described in the reference.  It has all of the capabilities and accuracy of the original program plus
 * the ability to include the approximate effects of angle of attack.
 * It is an order of magnitude faster, and improvements to the integration schemes
 * have reduced numerical integration errors by an order of magnitude.
 * A formatted input echo has been added so that those not intimately familiar
 * with the code can tell what has been input.
 *
 * Reference: Harris, Roy V., Jr.  An Analysis and Correlation of Aircraft Wave Drag.  NASA TMX-947.  March 1964.
 * "
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsAwave and
 * \ref aimOutputsAwave and \ref attributeAwave, respectively.
 *
 *
 * Upon running preAnalysis the AIM generates a single file, "awaveInput.txt" which contains the input
 * information and control sequence for Awave to execute (see \ref AwaveModification).
 *
 * The Awave AIM can automatically execute Awave, with details provided in \ref aimExecuteAwave.
 *
 * \section AwaveModification Awave Modifications
 * The AIM assumes that a modified version of Awave is being used. The modified version allows for longer input
 * and output file name lengths, as well as other I/O modifications. This modified version of
 * Awave, awavemod.f, is not currently supplied with the AIM due to licensing issues, please contact the CAPS creators
 * for additional details. Once this source code is obtained, it is automatically built with the AIM.
 * During compilation, the source code is compiled into an executable with the
 * name \a awave (Linux and OSX) or \a awave.exe (Windows)
 *
 * \section AwaveExample Examples
 * An example problem using the Awave AIM may be found at
 * \ref awaveExample.
 *
 */

 /*
 * \section assumptionsAwave AIM Assumptions
 *   All configurations are assumed to be symmetric with
 * respect to the X-Z plane so that only the "positive y" half of the
 * vehicle is defined.  The Awave AIM automatically generates the halfspan
 * inputs from full span Engineering Sketch Pad models that are oriented
 * with the coordinate system assumption (X -- downstream, Y -- out the
 * right wing, Z -- up).
 *
 *    Within <b> OpenCSM</b> there are a number of airfoil generation UDPs (User Defined Primitives). These include NACA 4
 * series, a more general NACA 4/5/6 series generator, Sobieczky's PARSEC parameterization and Kulfan's CST
 * parameterization. All of these UDPs generate <b> EGADS</b> <em> FaceBodies</em> where the <em>Face</em>'s underlying
 * <em>Surface</em>  is planar and the bounds of the <em>Face</em> is a closed set of <em>Edges</em> whose
 * underlying <em>Curves</em> contain the airfoil shape. In all cases there is a <em>Node</em> that represents
 * the <em>Leading Edge</em> point and one or two <em>Nodes</em> at the <em>Trailing Edge</em> -- one if the
 * representation is for a sharp TE and the other if the definition is open or blunt. If there are 2 <em>Nodes</em>
 * at the back, then there are 3 <em>Edges</em> all together and closed, even though the airfoil definition
 * was left open at the TE. All of this information will be used to automatically fill in the Awave geometry description.
 *
 * It should be noted that general construction in either <b> OpenCSM</b> or even <b> EGADS</b> will be supported
 * as long as the topology described above is used. But care should be taken when constructing the airfoil shape
 * so that a discontinuity (i.e.,  simply <em>C<sup>0</sup></em>) is not generated at the <em>Node</em> representing
 * the <em>Leading Edge</em>. This can be done by splining the entire shape as one and then intersecting the single
 *  <em>Edge</em> to place the LE <em>Node</em>.
 *
 */

/*! \page attributeAwave AIM Attributes
 * The following list of attributes drives the Awave geometric definition. Aircraft components are defined as cross sections
 * in the low fidelity geometry definition. To be able to logically group the cross sections into wings, tails, fuselage, etc.
 * they must be given a grouping attribute. This attribute defines a logical group along with identifying a set of cross sections
 * as a lifting surface or a body of revolution. The format is as follows.
 *
 *  - <b> capsType</b> This string attribute labels the <em> FaceBody</em> as to which type the section
 *  is assigned. This information is also used to logically group sections together to create wings, tails, stores, etc. Because
 *  Awave is relatively rigid <b> capsType </b> attributes must be on of the following items:
 *
 *    <em> Lifting Surfaces: </em> Wing, Tail, HTail, VTail, Cannard, Fin
 *
 *    <em> Body of Revolution: </em> Fuselage, Fuse, Store
 *
 *  - <b> capsGroup </b> This string attribute is used to group like components together.  This is a user defined
 *  unique string that can be used to tie sections to one another.  Examples are tail1, tail2, etc.
 *
 *  - <b> capsReferenceArea</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the SREF entry in the Awave input.
 *
 */

typedef struct {
    const char *name;
    const char *attribute;
    int type; // 0 lifting surface, 1 body of revolution (fuselage | store | pod | etc ...)
    double xyz[3]; // leading edge (type 0) or centroid (type 1) location
    // values for type 1 inputs
    double area; // section area
    double radius; // PI * radius * radius = area
    // values for type 0 inputs
    double chordLength;
    int ndiv; // number of chord divisions
    double *x; // chord locations length(*x) = ndiv, in (x/c)*100
    double *camber; // delta z chamber at each x location
    double *halfThick; // half thickness at x location (100 * (t/c) / 2
} awaveSec;


typedef struct {
  int          nVert;
  const double *pxyz, *pt;
} tessStorage;

typedef struct {
    // Inidicate if input angle should have units
    int withUnits;
} aimStorage;


/* ****************** Awave AIM Helper Functions ************************** */

static void initiate_Section(awaveSec *section) {

    section->name = NULL;
    section->attribute = NULL;
    section->type = 0; // 0 lifting surface, 1 body of revolution (fuselage | store | pod | etc ...)
    section->xyz[0] = 0.0; // leading edge (type 0) or centroid (type 1) location
    section->xyz[1] = 0.0;
    section->xyz[2] = 0.0;

    // values for type 1 inputs
    section->area = 0.0; // section area
    section->radius = 0.0; // PI * radius * radius = area

    // values for type 0 inputs
    section->chordLength = 0.0;
    section->ndiv = 0; // number of chord divisions
    section->x = NULL; // chord locations length(*x) = ndiv, in (x/c)*100
    section->camber = NULL; // delta z chamber at each x location
    section->halfThick = NULL; // half thickness at x location (100 * (t/c) / 2
}

static void destroy_Section(awaveSec *section) {

    section->name = NULL;
    section->attribute = NULL;
    section->type = 0; // 0 lifting surface, 1 body of revolution (fuselage | store | pod | etc ...)
    section->xyz[0] = 0.0; // leading edge (type 0) or centroid (type 1) location
    section->xyz[1] = 0.0;
    section->xyz[2] = 0.0;

    // values for type 1 inputs
    section->area = 0.0; // section area
    section->radius = 0.0; // PI * radius * radius = area

    // values for type 0 inputs
    section->chordLength = 0.0;
    section->ndiv = 0; // number of chord divisions
    if (section->x != NULL) EG_free(section->x); // chord locations length(*x) = ndiv, in (x/c)*100
    if (section->camber != NULL) EG_free(section->camber); // delta z chamber at each x location
    if (section->halfThick != NULL) EG_free(section->halfThick); // half thickness at x location (100 * (t/c) / 2
}

// length - size of array
// array - array of doubles
// str - character array that is written in the comment section of the first line (column 73+)
// fp - file pointer to write to, NULL writes to standard out
static void printAwaveArray(int length, double array[], char *str,
                            /*@null@*/ FILE *fp)
{
    int i;
    char * tmp;

    for (i = 0; i < length; i++) {

        //d2s7fld(array[i], tmp);
        tmp = convert_doubleToString(array[i],7,0);

        if (fp == NULL) {
            printf("%s",tmp);
            if ((double) (i+1) / 10.0 == (double) ((i+1)/10)) {
                if (i == 9 & str != NULL) {
                    printf("  %s\n",str);
                } else {
                    printf("\n");
                }
            }
        } else {
            fprintf(fp,"%s",tmp);
            if ((double) (i+1) / 10.0 == (double) ((i+1)/10)) {
                if (i == 9 & str != NULL) {
                    fprintf(fp,"  %s\n",str);
                } else {
                    fprintf(fp,"\n");
                }
            }
        }

        if(tmp != NULL) EG_free(tmp);

    }

    if ((double) (length) / 10.0 != (double) ((length)/10)) {
        if (length > 10 | str == NULL) {
            if (fp == NULL) {
                printf("\n");
            } else {
                fprintf(fp,"\n");
            }
        } else {
            for (i = 0; i < 10-length; i++) {
                if (fp == NULL) {
                    printf("       ");
                } else {
                    fprintf(fp,"       ");
                }
            }
            if (fp == NULL) {
                printf("  %s\n",str);
            } else {
                fprintf(fp,"  %s\n",str);
            }
        }
    }

}


// ********************** AIM Function Break *****************************
static int
linInterp2D(int length, double xIn[], double yIn[], double xOut, double* yOut)
{
    // find yOut given vectors xIn, yIn and desired xOut
    int i, i1, i2;
    i1 = i2 = -1;

    // find location in xIn / yIn for interpolation
    if (xIn[0] < xIn[length-1]) {
        if (xOut > xIn[length-1] | xOut < xIn[0]) {
            printf("linInterp2D does not support extrapolation\n");
            //status = printAwaveArray(length,xIn,"XIN",NULL);
            //status = printAwaveArray(length,yIn,"YIN",NULL);
            //printf("X DESIRED %4.3f\n",xOut);
            *yOut = 0.0;
            return CAPS_IOERR;
        }
        for (i = 0; i < length-1; i++) {
            if (xOut <= xIn[i+1] & xOut >= xIn[i]) {
                i1 = i;
                i2 = i+1;
                break;
            }
        }
    } else {
        if (xOut > xIn[0] | xOut < xIn[length-1]) {
            printf("linInterp2D does not support extrapolation\n");
            //status = printAwaveArray(length,xIn,"XIN",NULL);
            //status = printAwaveArray(length,yIn,"YIN",NULL);
            //printf("X DESIRED %4.3f\n",xOut);
            *yOut = 0.0;
            return CAPS_IOERR;
        }
        for (i = 0; i < length-1; i++) {
            if (xOut >= xIn[i+1] & xOut <= xIn[i]) {
                i1 = i;
                i2 = i+1;
                break;
            }
        }
    }

    // If location is not found error out - no extrapolation
    if (i1 != -1) {
        *yOut = yIn[i1] + (xOut - xIn[i1]) * (yIn[i2] - yIn[i1]) / (xIn[i2] - xIn[i1]);
    }

    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
static int
tessPointReturn(int length, double xOut, const ego edge, const double *xyz_tess,
                const double *t_tess, double *yOut, double *zOut)
{
    double dx, dX, t0, t1, t, s, result[18];
    int    i, status;

    *yOut = 0.0;
    *zOut = 0.0;

    for (i = 0; i < length-1; i++) {
        if ((xyz_tess[3*(i+1)] >= xOut && xyz_tess[3*(i  )] <= xOut) ||
            (xyz_tess[3*(i  )] >= xOut && xyz_tess[3*(i+1)] <= xOut)) {
          t0 = t_tess[i];
          t1 = t_tess[i+1];
          dX = xyz_tess[3*(i+1)] - xyz_tess[3*i];
          dx = xOut - xyz_tess[3*i];
          s = dx/dX;

          t = t0*(1.0-s) + s*t1;

          status = EG_evaluate(edge, &t, result);
          if (status != EGADS_SUCCESS) return status;

          *yOut = result[1];
          *zOut = result[2];

          return CAPS_SUCCESS;
        }
    }

    return CAPS_NOTFOUND;
}


// ********************** AIM Function Break *****************************
static int
defineAwavePod(int nbody, awaveSec sec[], int nPts, int *podLoc, double *xyz,
               double *xPod, double *radPod)
{
    // INPUTS
    // nbody - number of bodies / size of sec
    // sec - nbody size array of awaveSec structures
    // nPts - NPODOR number
    // podLoc - nbody array of integers most of them 0 with 1 indicating that
    //          the sec is part of the pod of interest

    // OUTPUTS
    // xyz - location of the front of the pod
    // xPod - x location starting at 0.0 in dimensional units for the pod
    // radPod - radius at each x location

    // it is likely that the number or sections in the input geometry (podLoc)
    //          and the nPts values are NOT the same
    // this function interpolates to determine xPod and radPod

    int    i, nsec, *iloc;
    double dx;
    double *xIn, *x, *rad, *y, *z, *X, *R, *Xtmp, *Rtmp, *Y, *Z, *Ytmp, *Ztmp;

    /* allocate our temp space */
    xIn = (double *) EG_alloc(nPts*sizeof(double));
    if (xIn == NULL) return EGADS_MALLOC;
    iloc = (int *) EG_alloc(nbody*sizeof(int));
    if (iloc == NULL) {
        EG_free(xIn);
        return EGADS_MALLOC;
    }
    x = (double *) EG_alloc(12*nbody*sizeof(double));
    if (x == NULL) {
        EG_free(iloc);
        EG_free(xIn);
        return EGADS_MALLOC;
    }
    rad  = &x[   nbody];
    y    = &x[ 2*nbody];
    z    = &x[ 3*nbody];
    X    = &x[ 4*nbody];
    R    = &x[ 5*nbody];
    Xtmp = &x[ 6*nbody];
    Rtmp = &x[ 7*nbody];
    Y    = &x[ 8*nbody];
    Z    = &x[ 9*nbody];
    Ytmp = &x[10*nbody];
    Ztmp = &x[11*nbody];

    // determine the number of cross sections making the pod input
    // keep the xyz information and the location in the list of all sections
    nsec = 0;
    for (i = 0; i < nbody; i++) {
        if (podLoc[i]) {
            x[i] = sec[i].xyz[0];
            y[i] = sec[i].xyz[1];
            z[i] = sec[i].xyz[2];
            rad[i] = sec[i].radius;
            iloc[nsec] = i;
            nsec++;
        }
    }

    // populate section variables with xyz and radius information
    for (i = 0; i < nsec; i++) {
        X[i] = Xtmp[i] = x[iloc[i]];
        Y[i] = Ytmp[i] = y[iloc[i]];
        Z[i] = Ztmp[i] = z[iloc[i]];
        R[i] = Rtmp[i] = rad[iloc[i]];
    }
    // make sure that X is increasing else fip all the data
    if (X[0] > X[nsec-1]) {
        for (i = 0; i < nsec; i++) {
            X[i] = Xtmp[nsec-1-i];
            Y[i] = Ytmp[nsec-1-i];
            Z[i] = Ztmp[nsec-1-i];
            R[i] = Rtmp[nsec-1-i];
        }
    }

    // return the front xyz location of the pod
    xyz[0] = X[0];
    xyz[1] = Y[0];
    xyz[2] = Z[0];

    // detrmine the even delta x for an awave pod input containing nPts
    dx = (X[nsec-1] - X[0]) / (double) (nPts -1);

    // create x arrays for the pod
    // awave requires on starting at 0.0 xPod
    // linear interpolation will need the one starting at the xyz front point location (xIn)
    xPod[0] = 0.0;
    xIn[0] = X[0];
    for (i = 1; i < nPts; i++) {
        xPod[i] = xPod[i-1] + dx;
        xIn[i] = xIn[i-1] + dx;
    }

    // interpolate along xIn to determine the radius of the pod at each xIn & xPod location
    for (i = 0; i < nPts; i++) {
        linInterp2D(nsec, X, R, xIn[i], &radPod[i]);
    }

    EG_free(x);
    EG_free(iloc);
    EG_free(xIn);
    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
static int
findSectionData(ego body, int *nle, double *xle, int *nte, double *xte,
                double *chordLength, double *arcLength, double *thickOverChord)
{
    //Input to this function is an egads body object
    //Find leading edge location x,y,z
    //Find trailing edge location x,y,z
    //Find chord length of section (distance between two points),
    //      both should have same y value if an egads body representing
    //      a cross section || to the yplane is entered
    //Find arc length around the body

    int    i, n, status, nNode, nEdge, oclass, mtype, *sens,index;
    double xmin=1e300, xyz[3], xmax=-1e300, box[6], data[4];
    double massData[14];
    double thickness;
    ego    ref, *nodes, *objs, *edges;

    *nle = 0;
    *nte = 0;
    *arcLength = 0.0;

    // check body type, looking for a NODE
    status = aim_isNodeBody(body, data);
    if (status < EGADS_SUCCESS) {
        printf(" Awave AIM Warning: EG_getTopology Failure in findSectionData Code :: %d\n", status);
        return status;
    }

    if (status == EGADS_SUCCESS) {
        *nle = 0;
        *nte = 0;
        *chordLength = 0.0;
        *arcLength = 0.0;
        *thickOverChord = 0.0;
        xle[0] = xte[0] = data[0];
        xle[1] = xte[1] = data[1];
        xle[2] = xte[2] = data[2];
        return CAPS_SUCCESS;
    }

    // assume the LE position is the most forward Node in X
    // assume the TE position is the most rearward Node in X

    // Get Nodes from EGADS body, input to this function
    status = EG_getBodyTopos(body, NULL, NODE, &nNode, &nodes);
    if (status != EGADS_SUCCESS) {
        printf(" Awave AIM Warning: LE getBodyTopos Nodes = %d\n", status);
        return status;
    }

    status = EG_getBoundingBox(body, box);
    if (status != EGADS_SUCCESS) {
        printf(" Awave AIM Warning: LE getBoundingBox Nodes = %d\n", status);
        EG_free(nodes);
        return status;
    }
    // Estimate thickness value
    thickness = sqrt((box[1]-box[4])*(box[1]-box[4]) + (box[2]-box[5])*(box[2]-box[5]));
    // Sort throught the list of nodes to determine the LE and TE nodes
    //printf("Nodes %d\n",nNode);
    index = 0;

    for (i = 0; i < nNode; i++) {
        status = EG_getTopology(nodes[i], &ref, &oclass, &mtype, xyz,
                              &n, &objs, &sens);
        if (status != EGADS_SUCCESS) continue;

        // LE search
        if (*nle == 0) {
            *nle = i+1;
            xmin = xyz[index];
        } else {
            if (xyz[index] < xmin) {
                *nle = i+1;
                xmin = xyz[index];
            }
        }
        // TE search
        if (*nte == 0) {
            *nte = i+1;
            xmax = xyz[index];
        } else {
            if (xyz[index] > xmax) {
                *nte = i+1;
                xmax = xyz[index];
            }
        }
    }

    if (*nle == 0) {
        if (nodes != NULL)EG_free(nodes);
        return  CAPS_BADVALUE;
    }
    if (*nte == 0) {
        if (nodes != NULL) EG_free(nodes);
        return CAPS_BADVALUE;
    }

    // assign xyz location variables xle and xte
    EG_getTopology(nodes[*nle-1], &ref, &oclass, &mtype, xle, &n, &objs, &sens);
    EG_getTopology(nodes[*nte-1], &ref, &oclass, &mtype, xte, &n, &objs, &sens);
    // determine distance between LT and TE points
    *chordLength = sqrt((xle[0]-xte[0])*(xle[0]-xte[0]) +
                        (xle[1]-xte[1])*(xle[1]-xte[1]) +
                        (xle[2]-xte[2])*(xle[2]-xte[2]));

    if ( fabs(*chordLength) < 1.0e-8 ) {
        *chordLength = box[3] - box[0];
    }

    *thickOverChord = thickness / *chordLength;
    if (nodes != NULL) EG_free(nodes);

    // determine the arc length around the body
    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    if (status != EGADS_SUCCESS) {
        printf(" Awave AIM Warning: LE getBodyTopos Edges = %d\n", status);
        return status;
    }

    for (i = 0; i < nEdge; i++) {

        status = EG_getMassProperties(edges[i], massData);

        if (status != EGADS_SUCCESS) continue;
        // massData array population
        // volume, surface area (length), cg(3), iniria(9)
        *arcLength = *arcLength + massData[1];
        // printf("EDGE %d :: %8.6f Length\n",i,massData[1]);

    }
    if (edges != NULL) EG_free(edges);

    return CAPS_SUCCESS;

}


// ********************** AIM Function Break *****************************
static int
defineAwaveAirfoil(ego body, int nPts, double *xOut, double *cOut, double *hOut)
{
    // INPUTS
    // body - egads body object
    // nPts - number of points along the x axis to define the airfoil
    //        note: awave assumes x-axis flow direction, thus this function
    //              is not fully general but awave specific
    // OUTPUTS
    // xOut - nPts long array of 100*(x/c) where c is the chord length
    // cOut - camber values at each xOut location
    // hOut - half thickness values at each xOut location 100*(t/c) / 2.0

    int          i, j, n, status, nle, nte, nEdge, *loc = NULL;
    int          nVert;
    const double *pxyz, *pt;
    double       dx, xle[3], xte[3], *x = NULL, dtmp, clen;
    double       *foil_xyz = NULL, chord; // sign; Never used
    double       yOut, zOut, yAvg, zAvg;
    double       parms[3];
    ego          tess = NULL, *edges = NULL;

    status = findSectionData(body, &nle, xle, &nte, xte, &clen, &dtmp, &dtmp);
    if (status != CAPS_SUCCESS) {
        printf(" Awave AIM Warning: findSectionData = %d!\n", status);
        return status;
    }

    x = (double *) EG_alloc(nPts*sizeof(double));
    if (x == NULL) {
        printf(" Awave AIM Warning: EG_alloc of x (nPts = %d)!\n", nPts);
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Output X location 100*(x/c)
    dx = 100.0 / (double) (nPts-1);
    xOut[0] = 0.0;

    for (i = 1; i < nPts; i++) {
        xOut[i] = xOut[i-1] + dx;
    }

    // Physical Space X location
    chord = xte[0] - xle[0];
    dx = (chord) / (double)(nPts-1);
    x[0] = xle[0];
    for (i = 1; i < nPts-1; i++) {
        x[i] = xle[0] + dx*i;
    }
    x[nPts-1] = xte[0];

    // a NODE has zero camber or thickness
    if (nle == 0 & nte == 0) {
        for (i = 0; i < nPts; i++) {
            cOut[i] = hOut[i] = 0.0;
        }
        status = CAPS_SUCCESS;
        goto cleanup;
    }

    status = EG_getBodyTopos(body, NULL, EDGE, &nEdge, &edges);
    if ((status != EGADS_SUCCESS) || (edges == NULL)) {
        printf(" Awave AIM Warning: LE getBodyTopos EDGE = %d\n", status);
        goto cleanup;
    }

    // Negating the first parameter triggers EGADS to only put vertexes on edges
    parms[0] = -clen / 100.0; // edge maximum length in physical space
    parms[1] =  clen / 10.0;  // curvature - based value, refines areas where curvature is high
    parms[2] =  15.0;         // max angle constraint, not required for this

    status = EG_makeTessBody(body, parms, &tess);
    if ((status != EGADS_SUCCESS) || (tess == NULL)) {
        if (status == EGADS_SUCCESS) status = EGADS_NOTTESS;
        printf(" Awave AIM Warning: EG_makeTessBody = %d\n", status);
        goto cleanup;
    }

    loc = (int *) EG_alloc(2*nPts*sizeof(int));
    if (loc == NULL) {
        printf(" Awave AIM Warning: EG_alloc of loc (nPts = %d)!\n", nPts);
        status = EGADS_MALLOC;
        goto cleanup;
    }
    foil_xyz = (double *) EG_alloc(3*2*nPts*sizeof(double));
    if (foil_xyz == NULL) {
        printf(" Awave AIM Warning: EG_alloc of foil_xyz (nPts = %d)!\n", nPts);
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < nPts; i++) {
        loc[2*i] = loc[2*i+1] = -1;
    }

    j = 0;
    for (n = 0; n < nEdge; n++) {
        status = EG_getTessEdge(tess, n+1, &nVert, &pxyz, &pt);
        if (status != EGADS_SUCCESS) {
            printf(" Awave AIM Warning: EG_getTessEdge = %d, Edge Number %d\n", status, n+1);
            goto cleanup;
        }

        for (i = 0; i < nPts; i++) {
            status = tessPointReturn(nVert, x[i], edges[n], pxyz, pt, &yOut, &zOut);
            // look for CAPS_SUCCESS otherwise extrapolation is present and no data is kept
            if (status == CAPS_SUCCESS) {
                foil_xyz[3*j+0] = x[i];
                foil_xyz[3*j+1] = yOut;
                foil_xyz[3*j+2] = zOut;
                if (loc[2*i] == -1) {
                    loc[2*i] = j;
                } else {
                    loc[2*i+1] = j;
                }
                j++;
            }
        }
        if (j == 2*nPts) break;
    }

    for (i = 0; i < nPts; i++) {
        // chamber is only input for wings not fins
        // define + chamber as +z
        // reference is the leading edge location xle[]

        if (loc[2*i] == -1 & loc[2*i+1] > -1) {
            j = loc[2*i+1];

            cOut[i] = sqrt((xle[1] - foil_xyz[3*j+1]) * (xle[1] - foil_xyz[3*j+1]) +
                           (xle[2] - foil_xyz[3*j+2]) * (xle[2] - foil_xyz[3*j+2]));
            hOut[i] = 0.0;
        } else if (loc[2*i] > -1 & loc[2*i+1] == -1) {
            j = loc[2*i];

            cOut[i] = sqrt((xle[1] - foil_xyz[3*j+1]) * (xle[1] - foil_xyz[3*j+1]) +
                           (xle[2] - foil_xyz[3*j+2]) * (xle[2] - foil_xyz[3*j+2]));
            hOut[i] = 0.0;
        } else if (loc[2*i] == -1 & loc[2*i+1] == -1) {
            cOut[i] = 0.0;
            hOut[i] = 0.0;
        } else {
            j = loc[2*i];
            n = loc[2*i+1];
            yAvg = (foil_xyz[3*n+1] + foil_xyz[3*j+1]) / 2.0;
            zAvg = (foil_xyz[3*n+2] + foil_xyz[3*j+2]) / 2.0;

            cOut[i] = sqrt((xle[1] - yAvg) * (xle[1] - yAvg) + (xle[2] - zAvg) * (xle[2] - zAvg));
            hOut[i] = (100.0 / chord) * sqrt((foil_xyz[3*n+1] - foil_xyz[3*j+1]) *
                                             (foil_xyz[3*n+1] - foil_xyz[3*j+1]) +
                                             (foil_xyz[3*n+2] - foil_xyz[3*j+2]) *
                                             (foil_xyz[3*n+2] - foil_xyz[3*j+2])) / 2.0;
        }
    }


    status = CAPS_SUCCESS;

cleanup:
    if (loc != NULL) EG_free(loc);
    if (foil_xyz != NULL) EG_free(foil_xyz);
    if (x != NULL) EG_free(x);

    if (edges != NULL) EG_free(edges);

    if (tess != NULL) (void) EG_deleteObject(tess);

    return status;
}


/* ********************** Exposed AIM Functions ***************************** */

// ********************** AIM Function Break *****************************
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS;
    const char *keyWord;
    char *keyValue = NULL, *angle = NULL;
    double real = 1;

    aimStorage *awaveInstance=NULL;

    #ifdef DEBUG
        printf("\n awaveAIM/aimInitialize   instance = %d!\n", inst);
    #endif

    // specify the number of analysis input and out "parameters"
    *nIn     = NUMINPUT;    // Mach, Altitude
    *nOut    = NUMOUTPUT;   // CDwave, MachOut, AoAOut
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks  = NULL;
    *fInOut  = NULL;

    // Allocate su2Instance
    AIM_ALLOC(awaveInstance, 1, aimStorage, aimInfo, status);
    *instStore = awaveInstance;
    awaveInstance->withUnits = (int)false;

    /*! \page aimUnitsAwave AIM Units
     *  A unit system may be optionally specified during AIM instance initiation. If
     *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
     *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
     *  unit consistent. A unit system may be specified via a JSON string dictionary for example:
     *  unitSys = "{"angle": "deg"}"
     */
    if (unitSys != NULL) {

      // Do we have a json string?
      if (strncmp( unitSys, "{", 1) != 0) {
        AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
        return CAPS_BADVALUE;
      }

      /*! \page aimUnitsAwave
       *  \section jsonStringAwave JSON String Dictionary
       *  The key arguments of the dictionary are described in the following:
       *
       *  <ul>
       *  <li> <B>angle = "None"</B> </li> <br>
       *  Angle units - e.g. "degree", "deg", "rad", ...
       *  </ul>
       */
      keyWord = "angle";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        angle = string_removeQuotation(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, angle, &real, "degree", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, angle, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      awaveInstance->withUnits = (int)true;
    }

cleanup:
    AIM_FREE(keyValue);
    AIM_FREE(angle);

    return status;
}

// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
  aimStorage *awaveInstance;

  *ainame = NULL;

  awaveInstance = (aimStorage *) instStore;

  /*! \page aimInputsAwave AIM Inputs
   * The following list outlines the Awave inputs along with their default value available
   * through the AIM interface.  All inputs to the Awave AIM are variable length arrays.
   * <B> All inputs must be the same length</B>.
   */

#ifdef DEBUG
  printf(" awaveAIM/aimInputs  index = %d!\n",  index);
#endif

    if (index == inMach) {
        *ainame               = EG_strdup("Mach");
        defval->limits.dlims[0] = 1.0; // Limit of accepted values
        defval->limits.dlims[1] = 100.0;


        /*! \page aimInputsAwave
         * - <B> Mach = double </B> <br> OR
         * - <B> Mach = [double, ... , double] </B> <br>
         *  Mach number.
         */

    } else if (index == inAlpha) {
        *ainame               = EG_strdup("Alpha");
        if (awaveInstance != NULL && awaveInstance->withUnits == (int)true)
          defval->units       = EG_strdup("degree");

        /*! \page aimInputsAwave
         * - <B> Alpha = double </B> <br> OR
         * - <B> Alpha = [double, ... , double] </B> <br>
         *  Angle of attack [degree].
         */

    }

    defval->type          = Double;
    defval->lfixed        = Change;
    defval->sfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;
    defval->vals.real = 0.0;
    defval->vals.reals = NULL;

    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(/*@unused@*/ void *instStore, void *aimInfo,
                   /*@null@*/ capsValue *inputs)
{
    int status; // Function status return

    int         nbody, nFace, i, j, n, nPod, nFin, nCan;

    int         atype, alen, nle, nte;
    int         nAwaveType[5];
    int         J0, J1, J2, J3, J4, J5, J6;
    int         NWAF, NWAFOR, NFUS, NRADX[4], NFORX[4];
    int         NP, NPODOR, NF, NFINOR, NCAN, NCANOR;

    double Sref, dtmp, data[4];

    // Freeable storage values
    int *iloc = NULL;
    int *awaveType = NULL;
    int *locPod = NULL, *locFin = NULL, *locCan = NULL, *locFuse = NULL, *locWing = NULL;

    double *xFuse = NULL, *zFuse = NULL, *aFuse = NULL;
    double *xPod = NULL, *rPod = NULL;

    awaveSec    *surfaces = NULL; // awaveSec structure defined global

    // EGADS types
    ego *bodies, *faces = NULL, *nodes = NULL;

    // capsIntent
    //int intent;

    // File
    FILE *fp = NULL;

    // create data arrays for fuselage input
    // These are hard coded to size 30 inside Awave
    double      xFuse1[30], xFuse2[30], xFuse3[30], xFuse4[30];
    double      zFuse1[30], zFuse2[30], zFuse3[30], zFuse4[30];
    double      aFuse1[30], aFuse2[30], aFuse3[30], aFuse4[30];
    double      xte[3], massData[14];
    double      xyzPod[3];
    double      Mach, AoA;
    char        tmpName[50];
    char        *tmpChar;
    const int   *ints;
    const char  *string, *intents;
    const double *reals;

    #ifdef DEBUG
        printf(" awaveAIM/aimPreAnalysis\n");
    #endif

    Sref = 1.0; // initilize Sref


    // Get EGADS bodies
    status = aim_getBodies(aimInfo, &intents, &nbody, &bodies);
    if (status != CAPS_SUCCESS) {
        printf(" awaveAIM/aimPreAnalysis getBodies = %d!\n", status);
        return status;
    }

    if ((nbody <= 0) || (bodies == NULL)) {

        #ifdef DEBUG
            printf(" awaveAIM/aimPreAnalysis No Bodies!\n");
        #endif

       return CAPS_SOURCEERR;
    }

    // Check inputs
    if (inputs == NULL) {
        #ifdef DEBUG
            printf(" awaveAIM/aimPreAnalysis inputs == NULL!\n");
        #endif

        return CAPS_NULLVALUE;
    }

    if (inputs[inMach-1].nullVal  == IsNull ||
        inputs[inAlpha-1].nullVal == IsNull) {

        printf("Either input Mach or Alpha has not been set!\n");
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    if (inputs[inMach-1].length != inputs[inAlpha-1].length) {

        printf("Inputs Mach and Alpha must be the same length\n");
        status = CAPS_MISMATCH;
        goto cleanup;
    }

    // Allocate memory to surfaces structure and awaveType array of ints
    awaveType = (int *) EG_alloc(nbody*sizeof(int));
    if (awaveType == NULL) {
        printf(" awaveAIM/aimPreAnalysis Cannot allocate %d types!\n", nbody);
        status = EGADS_MALLOC;
        goto cleanup;
    }
    surfaces = (awaveSec *) EG_alloc(nbody*sizeof(awaveSec));
    if (surfaces == NULL) {
        printf(" awaveAIM/aimPreAnalysis Cannot allocate %d surfaces!\n", nbody);
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Define default inputs for each awaveSec structure
    for (i = 0; i < nbody; i++) {
        (void) initiate_Section(&surfaces[i]);
    }

    //  Memory
    locPod  = (int *)    EG_alloc(nbody*sizeof(int));
    locCan  = (int *)    EG_alloc(nbody*sizeof(int));
    locFin  = (int *)    EG_alloc(nbody*sizeof(int));
    locFuse = (int *)    EG_alloc(nbody*sizeof(int));
    locWing = (int *)    EG_alloc(nbody*sizeof(int));
    iloc    = (int *)    EG_alloc(nbody*sizeof(int));
    xFuse   = (double *) EG_alloc(nbody*sizeof(double));
    zFuse   = (double *) EG_alloc(nbody*sizeof(double));
    aFuse   = (double *) EG_alloc(nbody*sizeof(double));

    if ((locPod  == NULL) || (locCan  == NULL) || (locFin == NULL) ||
        (locFuse == NULL) || (locWing == NULL) || (iloc   == NULL) ||
        (xFuse   == NULL) || (zFuse   == NULL) || (aFuse  == NULL)) {

        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Initialize loc--- to zeros
    for (i = 0; i < nbody; i++) {
        locPod[i] = locCan[i] = locFuse[i] = locFin[i] = locWing[i] = iloc[i] = 0;
        xFuse[i] = zFuse[i] = aFuse[i] = 0.0;
    }

    // ****************************************************************** //
    // PRE ANALYSIS CODE //

    printf("Writing Awave input file\n");

    fp = aim_fopen(aimInfo, "awaveInput.txt","w");

    if (fp == NULL) {
      AIM_ERROR(aimInfo, " awaveAIM/aimPreAnalysis Cannot Open File: awaveInput.txt!\n");
      status = CAPS_IOERR;
      goto cleanup;
    }

    fprintf(fp,"CAPS Awave AIM GENERATED INPUT\n");

    // ****************************************************************** //

    //Prelimary inputs:
    J0=1;//          1- 3     = 0, No reference area
    //                        = 1, Reference area will be input and used
    //                        = 2, Reference area from previous
    //                             configuration will be used

    J1=0;          //4- 6     = 0, No wing
    //                        = 1, Cambered wing data will be input
    //                        =-1, Uncambered wing data will be input
    //                        = 2, Wing data from previous configuration
    //                             will be used

    J2=0;//          7- 9     = 0, No fuselage
    //                        = 1, Arbitrarily shaped fuselage
    //                        =-1, Circular fuselage
    //                        = 2, Fuselage data from previous
    //                             configuration will be used

    J3=0;//         10-12     = 0, No pod
    //                        = 1, Pod data will be input
    //                        = 2, Pod data from previous configuration
    //                             will be used

    J4=0;//         13-15     = 0, No fin (Vertical tail)
    //                        = 1, Fin data will be input
    //                        = 2, Fin data from previous configuration
    //                             will be used

    J5=0;//         16-18     = 0, No canard (Horizontal tail)
    //                        = 1, Canard data will be input
    //                        = 2, Canard data from previous
    //                             configuration will be used

    J6=0;//         19-21     = 1, Complete configuration is symmetrical
    //                             with respect to the X-Y plane, implies
    //                             an uncambered circular fuselage if a
    //                             fuselage exists.
    //                        =-1, The circular fuselage is symmetrical
    //                             (uncambered)
    //                        = 0, No X-Y plane symmetry

    //NWAF          22-24     Number of airfoil sections used to describe
    //                        the wing (minimum = 2, maximum = 20)

    NWAFOR=21;//     25-27     Number of stations at which ordinates are
    //                        input for each wing airfoil section
    //                        (minimum = 3, maximum = 30)
    //                        If given a negative sign, the program
    //                        will expect to read both upper and lower
    //                        ordinates.  If positive, the airfoil is
    //                        assumed to be symmetrical.

    //NFUS          28-30     Number of fuselage segments
    //                        (minimum = 1, maximum = 4)

    NRADX[0]=21;//   31-33     Number of points per station used to
    //                        represent a half-section of the first
    //                        fuselage segment.  For circular fuselages,
    //                        NRADX y and z coordinates are computed
    //                        (minimum = 3, maximum = 30)

    NFORX[0]=0;//   34-36     Number of stations for first fuselage
    //                        segment (minimum = 4, maximum = 30)

    NRADX[1]=21;//   37-39     Same as NRADX(1) but for second segment

    NFORX[1]=0;//   40-42     Same as NFORX(1) but for second segment

    NRADX[2]=21;//   43-45     Same as NRADX(1) but for third segment

    NFORX[2]=0;//   46-48     Same as NFORX(1) but for third segment

    NRADX[3]=21;//   49-51     Same as NRADX(1) but for fourth segment

    NFORX[3]=0;//   52-54     Same as NFORX(1) but for fourth segment

    //NP            55-57     Number of pods described (maximum = 9)

    NPODOR=10;//     58-60     Number of stations at which pod radii are
    //           input (minimum = 4, maximum = 30)

    //NF            61-63     Number of fins (vertical tails) described
    //                        (maximum = 6) If given a negative sign, the
    //                        program will expect to read ordinates for
    //                        both the outboard and inboard sides of the
    //                        fins.  If positive, the airfoils are
    //                        assumed to be symmetrical.

    NFINOR=10;//     64-66     Number of stations at which ordinates are
    //                        input for each fin airfoil
    //                        (minimum = 3, maximum = 10)

    //NCAN          67-69     Number of canards (horizontal tails)
    //                        defined (maximum = 2).  If given a negative
    //                        sign, the program will expect to read
    //                        ordinates for both the root and tip.  If
    //                        positive, the root and tip airfoils are
    //                        assumed to be the same.

    NCANOR=10;//     70-72     Number of stations at which ordinates are
    //                        input for each canard airfoil
    //                        (minimum = 3, maximum = 10)
    //                        If given a negative sign, the program
    //                        will expect to read both upper and lower
    //                        ordinates.  If positive, the airfoil is
    //                        assumed to be symmetrical.

    // loop over all the bodies and populate the surfaces structure
    // at this stage even surfaces located in -Y are retained
    // populate awaveType only if surface will be printed (Y >= 0.0)
    // 1 (wing), 2 (fuselage), 3 (pod), 4 (vertical fin), 5 (cannard)
    nAwaveType[0] = nAwaveType[1] = nAwaveType[2] = nAwaveType[3] = nAwaveType[4] = 0;

    for (i = 0; i < nbody; i++) {

        awaveType[i] = 0; // default is 0, nothing will be written to awave input, section is -Y

        // search for the "capsReferenceArea" attribute
        status = EG_attributeRet(bodies[i], "capsReferenceArea", &atype, &alen,
                                 &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {

            if (atype != ATTRREAL) {
                printf("capsReferenceArea should be followed by a single real value!\n");
                status = EGADS_ATTRERR;
                goto cleanup;
            }

            Sref = (double) reals[0];
        }

        // Determine type of body
        // Look for "capsType" attribute
        // Wing, Tail, VTail, HTail, Cannard, Fin are all lifting surfaces
        // Fuse, Fuselage, Store are all bodies of revolution
        status = EG_attributeRet(bodies[i], "capsType", &atype, &alen, &ints,
                                 &reals, &surfaces[i].attribute);
        if (status != EGADS_SUCCESS) {

            printf(" *** WARNING AwaveAIM: capsType not found on body %d - defaulting to 'Wing'!\n", i+1);
/*@-observertrans@*/
            surfaces[i].attribute = "Wing";
/*@+observertrans@*/

        } else {

            if (atype != ATTRSTRING) {

                printf("capsType should be followed by a single string!\n");
                status = EGADS_ATTRERR;
                goto cleanup;
            }
        }

        // Look for name attribute and populate structure
        // capsGroup is a unique user defined string for a component
        // components can be of the same capsType with a different capsGroup

        // Retrieve the string following a capsGroup tag
        status =retrieve_CAPSGroupAttr(bodies[i], &surfaces[i].name);
        if (status != CAPS_SUCCESS) {
            if (status == EGADS_NOTFOUND) {
                printf("Warning: capsGroup not found on body %d!\n", i+1);
            } else {
                goto cleanup;
            }
        }

        // Identify lifting surfaces
        if (strcmp(surfaces[i].attribute,"Wing")    == 0 ||
            strcmp(surfaces[i].attribute,"Tail")    == 0 ||
            strcmp(surfaces[i].attribute,"HTail")   == 0 ||
            strcmp(surfaces[i].attribute,"VTail")   == 0 ||
            strcmp(surfaces[i].attribute,"Cannard") == 0 ||
            strcmp(surfaces[i].attribute,"Fin")     == 0   ) {

            // Extract data about the surface
            status = findSectionData(bodies[i], &nle, surfaces[i].xyz, &nte,
                                     xte, &surfaces[i].chordLength, &dtmp, &dtmp);
            if (status != CAPS_SUCCESS) {
                printf(" awaveAIM/aimPreAnalysis findSectionData = %d!\n",
                       status);
               goto cleanup;
            }

            // Determine the type of lifting surface (Wing, Fin, Cannard)
            if (strcmp(surfaces[i].attribute,"Wing") == 0) {

                J1 = 1; // awave header variable indicating a wing will be input

                surfaces[i].ndiv = NWAFOR; // number of divisions for a wing
                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 1;
                    nAwaveType[0]++;
                    locWing[i] = 1;
                }

            } else if (strcmp(surfaces[i].attribute,"Tail") == 0) {

                J5 = 1; // awave header variable indicating a cannard (horizontal tail) will be input

                surfaces[i].ndiv = NCANOR; // number of divisions for a cannard

                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 5;
                    nAwaveType[4]++;
                }

            } else if (strcmp(surfaces[i].attribute,"HTail") == 0) {

                J5 = 1; // awave header variable indicating a cannard (horizontal tail) will be input

                surfaces[i].ndiv = NCANOR; // number of divisions for a cannard
                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 5;
                    nAwaveType[4]++;
                }

            } else if (strcmp(surfaces[i].attribute,"VTail") == 0) {

                J4 = 1; // awave header variable indicating a fin will be input

                surfaces[i].ndiv = NFINOR; // number of divisions for a fin
                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 4;
                    nAwaveType[3]++;
                }

            } else if (strcmp(surfaces[i].attribute,"Fin") == 0) {

                J4 = 1; // awave header variable indicating a fin will be input

                surfaces[i].ndiv = NFINOR; // number of divisions for a fin
                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 4;
                    nAwaveType[3]++;
                }

            } else if (strcmp(surfaces[i].attribute,"Cannard") == 0) {

                J5 = 1; // awave header variable indicating a cannard (horizontal tail) will be input

                surfaces[i].ndiv = NCANOR; // number of divisions for a cannard
                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 5;
                    nAwaveType[4]++;
                }
            }

            // based on the type of lifting surface allocation memory to store its shape
            surfaces[i].x         = (double *) EG_alloc(surfaces[i].ndiv*sizeof(double));
            surfaces[i].camber    = (double *) EG_alloc(surfaces[i].ndiv*sizeof(double));
            surfaces[i].halfThick = (double *) EG_alloc(surfaces[i].ndiv*sizeof(double));

            if ((surfaces[i].x         == NULL) ||
                (surfaces[i].camber    == NULL) ||
                (surfaces[i].halfThick == NULL)   ) {

                printf(" awaveAIM/aimPreAnalysis %d ndiv = %d Allocation!\n",
                       i+1, surfaces[i].ndiv);

                status = EGADS_MALLOC;
                goto cleanup;
            }

            // Populate the lifting surface shape
            status = defineAwaveAirfoil(bodies[i],
                                        surfaces[i].ndiv,
                                        surfaces[i].x,
                                        surfaces[i].camber,
                                        surfaces[i].halfThick);
            if (status != CAPS_SUCCESS) goto cleanup;

            // print the result to standard out
            printf("Lifting Surface Section, Body ID %d\n", i+1);
            printf("\tXLE:   %8.6f %8.6f %8.6f\n", surfaces[i].xyz[0],
                   surfaces[i].xyz[1], surfaces[i].xyz[2]);
            printf("\tCHORD: %8.6f\n", surfaces[i].chordLength);
            printf("\tATTRIB:%s\n", surfaces[i].attribute);
            printf("\tNAME  :%s\n", surfaces[i].name);
            printf("\tAwave TYPE: %d\n", awaveType[i]);

#ifdef DEBUG
            printAwaveArray(surfaces[i].ndiv, surfaces[i].x, "X", NULL);
            printAwaveArray(surfaces[i].ndiv, surfaces[i].camber, "CAMBER", NULL);
            printAwaveArray(surfaces[i].ndiv, surfaces[i].halfThick,
                            "HALF THICKNESS", NULL);
#endif

            printf("\n");

        // Identify bodies (Fuselage, Pod)
        }  else if (strcmp(surfaces[i].attribute,"Fuse")     == 0 ||
                    strcmp(surfaces[i].attribute,"Fuselage") == 0 ||
                    strcmp(surfaces[i].attribute,"Store")    == 0   ) {

            // Determine the type of body
            if (strcmp(surfaces[i].attribute,"Fuse")     == 0 ||
                strcmp(surfaces[i].attribute,"Fuselage") == 0    ) {

                J2 = -1; // awave header variable indicating a fuselage will be input

                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 2;
                    nAwaveType[1]++;
                    locFuse[i] = 1;
                }

            } else if (strcmp(surfaces[i].attribute,"Store") == 0) {

                J3 = 1; // awave header variable indicating a pod will be input

                if (surfaces[i].xyz[1] >= 0.0) {
                    awaveType[i] = 3;
                    nAwaveType[2]++;
                }
            }

            // Set type to body of revolution
            surfaces[i].type = 1;

            // Determine of the body is a NODE
            status = aim_isNodeBody(bodies[i], data);
            if (status < EGADS_SUCCESS) goto cleanup;

            // Handle the special case of a NODE body
            if (status == EGADS_SUCCESS) {

                surfaces[i].area = 0.0;
                surfaces[i].radius = 0.0;
                surfaces[i].xyz[0] = data[0];
                surfaces[i].xyz[1] = data[1];
                surfaces[i].xyz[2] = data[2];

            } else {

                // Determine the geometric properties of the body
                if (faces != NULL) EG_free(faces);

                status = EG_getBodyTopos(bodies[i], NULL, FACE, &nFace, &faces);
                if ((status != EGADS_SUCCESS) || (faces == NULL)) goto cleanup;

                if (nFace != 1) {
                    printf(" awaveAIM/aimPreAnalysis body %d with %d faces should only have one face!\n",
                           i+1, nFace);
                    status = CAPS_BADOBJECT;
                    goto cleanup;
                }

                status = EG_getMassProperties(faces[0], massData);
                if (status != CAPS_SUCCESS) goto cleanup;


                // volume, surface area (length), cg(3), iniria(9)
                // for 2d, volume is area
                surfaces[i].area = massData[1];
                surfaces[i].radius = sqrt(massData[1] / PI);
                surfaces[i].xyz[0] = massData[2];
                surfaces[i].xyz[1] = massData[3];
                surfaces[i].xyz[2] = massData[4];
            }

            // Print the information about the body to standard out
            printf("Body Section, Body ID %d\n", i+1);
            printf("\tXCG:   %8.6f %8.6f %8.6f\n", surfaces[i].xyz[0],
                   surfaces[i].xyz[1], surfaces[i].xyz[2]);
            printf("\tAREA:  %8.6f\n", surfaces[i].area);
            printf("\tRAD:   %8.6f\n", surfaces[i].radius);
            printf("\tATTRIB:%s\n", surfaces[i].attribute);
            printf("\tNAME  :%s\n", surfaces[i].name);
            printf("\tAwave TYPE: %d\n\n", awaveType[i]);
        }
    }

    // Print the number of each Awave type to the screen
    // nAwaveType[0] = J1 etc ...
    printf("AwaveTypes: %d %d %d %d %d\n", nAwaveType[0],
                                           nAwaveType[1],
                                           nAwaveType[2],
                                           nAwaveType[3],
                                           nAwaveType[4]);

    // Determine input values for awave header values updating the many Prelimary inputs defined above
    NWAF = nAwaveType[0]; // number of wing airfoil sections

    // Fuselage input is different, up to 4 segments with 30 stations each
    NFUS = 1 + nAwaveType[1] / 30;
    if (nAwaveType[1] == 60) {
        NFUS = 3;
    } else if (nAwaveType[1] == 89) {
        NFUS = 4;
    } else if (nAwaveType[1] > 117) {
        NFUS = 5; // cause CAPS_IOERR in switch statement
    }

    // Must repeat final - initial section
    // 1 2 3 4 ... 29 30
    // 30 31   ... 58 59
    // 59 60   ... 87 88
    // 88 89   ... 116 117
    switch (NFUS) {
        case 1:
            NFORX[0] = nAwaveType[1];
            break;
        case 2:
            NFORX[0] = 30;
            NFORX[1] = nAwaveType[1] - 29;
            break;
        case 3:
            NFORX[0] = 30;
            NFORX[1] = 30;
            NFORX[2] = nAwaveType[1] - 58;
            break;
        case 4:
            NFORX[0] = 30;
            NFORX[1] = 30;
            NFORX[2] = 30;
            NFORX[3] = nAwaveType[1] - 87;
            break;
        default:
            printf("Number of fuselage stations %d is too large for Awave\n",
                   nAwaveType[1]);
            status = CAPS_RANGEERR;
            goto cleanup;
    }

    // xFuse etc is a single array the contains all fuselage information
    // populate
    j = 0;
    for (i = 0; i < nbody; i++) {
        if (locFuse[i]) {
            xFuse[i] = surfaces[i].xyz[0];
            zFuse[i] = surfaces[i].xyz[2];
            aFuse[i] = surfaces[i].area;
            iloc[j] = i;
            j++; // total # of fuselage sections
        }
    }

    // Awave can only accept 4, 30 section fuselage inputs at a time
    // break xFuse, etc into individual arrays of size 30
    n = 0;
    for (i = 0; i < j; i++) {
        if (i < 30) {
            xFuse1[n] = xFuse[iloc[i]];
            zFuse1[n] = zFuse[iloc[i]];
            aFuse1[n] = aFuse[iloc[i]];
            n++;
        } else if (i == 30) {
            xFuse1[n] = xFuse[iloc[i]];
            zFuse1[n] = zFuse[iloc[i]];
            aFuse1[n] = aFuse[iloc[i]];
            xFuse2[0] = xFuse[iloc[i]];
            zFuse2[0] = zFuse[iloc[i]];
            aFuse2[0] = aFuse[iloc[i]];
            n = 1;
        } else if (i < 59) {
            xFuse2[n] = xFuse[iloc[i]];
            zFuse2[n] = zFuse[iloc[i]];
            aFuse2[n] = aFuse[iloc[i]];
            n++;
        } else if (i == 59) {
            xFuse2[n] = xFuse[iloc[i]];
            zFuse2[n] = zFuse[iloc[i]];
            aFuse2[n] = aFuse[iloc[i]];
            xFuse3[0] = xFuse[iloc[i]];
            zFuse3[0] = zFuse[iloc[i]];
            aFuse3[n] = aFuse[iloc[i]];
            n = 1;
        } else if (i < 88) {
            xFuse3[n] = xFuse[iloc[i]];
            zFuse3[n] = zFuse[iloc[i]];
            aFuse3[n] = aFuse[iloc[i]];
            n++;
        } else if (i == 88) {
            xFuse3[n] = xFuse[iloc[i]];
            zFuse3[n] = zFuse[iloc[i]];
            aFuse3[n] = aFuse[iloc[i]];
            xFuse4[0] = xFuse[iloc[i]];
            zFuse4[0] = zFuse[iloc[i]];
            aFuse4[0] = aFuse[iloc[i]];
            n = 1;
        } else if (i < 117) {
            xFuse4[n] = xFuse[iloc[i]];
            zFuse4[n] = zFuse[iloc[i]];
            aFuse4[n] = aFuse[iloc[i]];
            n++;
        }
    }

    // Report the number of sections in each of the 4 fuslage entries to standard out
    printf("Points for all four fuselage sections :: %d %d %d %d\n", NFORX[0],
                                                                     NFORX[1],
                                                                     NFORX[2],
                                                                     NFORX[3]);

    // Loop over the surfaces[i] structures and determine the number of
    // pods (J3), fins (J4) and cannards (J5)
    // these are in logical "capsType" attribute containers
    // each as a unique "capsGroup" attribute

    // NP number of pods (max 9) capsTypes Store
    // NF number of fins (max 6) capsTypes VTail, Fin
    // NCAN number of cannards (max 2) capsTypes Cannard, Tail, HTail
    // all are set to 0 above

    nPod = 0;
    nFin = 0;
    nCan = 0;

    for (i = 0; i < nbody; i++) {
        if (surfaces[i].xyz[1] >= 0.0) {
            if (strcmp(surfaces[i].attribute,"Store") == 0 ) {

                if (nPod == 0 | i == 0) {
                    nPod++;
                    locPod[i] = nPod;

                } else if (strcmp(surfaces[i].name,surfaces[i-1].name) == 0) {

                    // If the input is a pod type with the same name as before then
                    locPod[i] = nPod; // this is a section of the same pod

                } else {

                    // If the input is a pod type with a different name as before then
                    nPod++;
                    locPod[i] = nPod; // this is a section of a new pod
                }
            } else if (strcmp(surfaces[i].attribute,"VTail") == 0 ||
                       strcmp(surfaces[i].attribute,"Fin")   == 0   ) {

                if (nFin == 0 | i == 0) {
                    nFin++;
                    locFin[i] = nFin;

                } else if (strcmp(surfaces[i].name,surfaces[i-1].name) == 0) {

                    // If the input is a pod type with the same name as before then
                    locFin[i] = nFin; // this is a section of the same pod

                } else {

                    // If the input is a pod type with a different name as before then
                    nFin++;
                    locFin[i] = nFin; // this is a section of a new pod
                }

            } else if (strcmp(surfaces[i].attribute,"Cannard") == 0 ||
                       strcmp(surfaces[i].attribute,"Tail")    == 0 ||
                       strcmp(surfaces[i].attribute,"HTail")   == 0   ) {

                if (nCan == 0 | i == 0) {
                    nCan++;
                    locCan[i] = nCan;

                } else if (strcmp(surfaces[i].name,surfaces[i-1].name) == 0) {

                    // If the input is a pod type with the same name as before then
                    locCan[i] = nCan; // this is a section of the same pod

                } else {

                    // If the input is a pod type with a different name as before then
                    nCan++;
                    locCan[i] = nCan; // this is a section of a new pod
                }
            }
        }
    }

    NP = nPod;
    NF = nFin;
    NCAN = nCan;

    // Check for awave maximums
    if (NP > 9) {
        printf("Error: Awave can only handle 9 pods, pods entered :: %d\n", NP);
        status = CAPS_RANGEERR;
        goto cleanup;
    }

    if (NF > 6) {
        printf("Error: Awave can only handle 6 fins, fins entered :: %d\n", NF);
        status = CAPS_RANGEERR;
        goto cleanup;
    }

    if (NCAN > 2) {
        printf("Error: Awave can only handle 2 cannards, cannards entered :: %d\n",
               NCAN);
        status = CAPS_RANGEERR;
        goto cleanup;
    }

    // ****************************************************************** //
    // ****************************************************************** //
    // ANALYSIS CODE INPUT FILE GENERATION //
    // ****************************************************************** //
    // ****************************************************************** //
    // Write Awave HEADER //
    tmpChar = convert_integerToString(J0,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J1,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J2,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J3,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J4,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J5,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(J6,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NWAF,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NWAFOR,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NFUS,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    for (i = 0; i < 4; i++) {

        tmpChar = convert_integerToString(NRADX[i],3,0); fprintf(fp,"%s",tmpChar);
        if (tmpChar != NULL) EG_free(tmpChar);

        tmpChar = convert_integerToString(NFORX[i],3,0); fprintf(fp,"%s",tmpChar);
        if (tmpChar != NULL) EG_free(tmpChar);
    }

    tmpChar = convert_integerToString(NP,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NPODOR,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NF,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NFINOR,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NCAN,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    tmpChar = convert_integerToString(NCANOR,3,0); fprintf(fp,"%s",tmpChar);
    if (tmpChar != NULL) EG_free(tmpChar);

    fprintf(fp,"\n");
    // END WRITE HEADER

    // REF AREA
    if (J0 == 1) {

        tmpChar = convert_doubleToString(Sref,7,0);
        fprintf(fp,"%s                                                                 ",
                tmpChar);
        if (tmpChar != NULL) EG_free(tmpChar);

        // In Awave colums 73+ on each line are used for comment space
        fprintf(fp,"REF AREA\n");
    }

    // WING DATA
    if (J1) {

        // PRINT CHORD LOCATIONS (x/c) * 100
        // This is the same for every section
        j = 1;
        for (i = 0; i < nbody; i++) {

            if (locWing[i]) {

                if (j) { // only print the x locations once

                    printAwaveArray(surfaces[i].ndiv, surfaces[i].x,
                                    "WING DATA 100(x/c)", fp);
                    j = 0;
                }

                // Print Leading Edge x,y,z then chord length for each section
                // did not use printAwaveArray because xyz and chordLength are not in the same array
                tmpChar = convert_doubleToString(surfaces[i].xyz[0],7,0);
                fprintf(fp,"%s",tmpChar);
                if (tmpChar != NULL) EG_free(tmpChar);

                tmpChar = convert_doubleToString(surfaces[i].xyz[1],7,0);
                fprintf(fp,"%s",tmpChar);
                if (tmpChar != NULL) EG_free(tmpChar);

                tmpChar = convert_doubleToString(surfaces[i].xyz[2],7,0);
                fprintf(fp,"%s",tmpChar);
                if (tmpChar != NULL) EG_free(tmpChar);

                tmpChar = convert_doubleToString(surfaces[i].chordLength,7,0);
                fprintf(fp,"%s",tmpChar);
                if (tmpChar != NULL) EG_free(tmpChar);
                fprintf(fp,"                                            CS%d X,Y,Z,CHORD\n",i);
            }
        }

        // Print camber and half thickness data for all wing sections
        for (i = 0; i < nbody; i++) {

            if (locWing[i]) {

/*@-bufferoverflowhigh@*/
                sprintf(tmpName,"CS%d CAMBER",i);
                printAwaveArray(surfaces[i].ndiv, surfaces[i].camber, tmpName,
                                fp);

                sprintf(tmpName,"CS%d HALF THICK",i);
                printAwaveArray(surfaces[i].ndiv, surfaces[i].halfThick, tmpName,
                                fp);
/*@+bufferoverflowhigh@*/
            }
        }
    }

    // PRINT FUSELAGE INPUTS
    if (J2) {

        if (NFORX[0]) {
            printAwaveArray(NFORX[0], xFuse1, "FUSE1 X" , fp);
            printAwaveArray(NFORX[0], zFuse1, "FUSE1 Z" , fp);
            printAwaveArray(NFORX[0], aFuse1, "FUSE1 AREA" , fp);
        }

        if (NFORX[1]) {
            printAwaveArray(NFORX[1], xFuse2, "FUSE2 X" , fp);
            printAwaveArray(NFORX[1], zFuse2, "FUSE2 Z" , fp);
            printAwaveArray(NFORX[1], aFuse2, "FUSE2 AREA" , fp);
        }

        if (NFORX[2]) {
            printAwaveArray(NFORX[2], xFuse3, "FUSE3 X" , fp);
            printAwaveArray(NFORX[2], zFuse3, "FUSE3 Z" , fp);
            printAwaveArray(NFORX[2], aFuse3, "FUSE3 AREA" , fp);
        }
        if (NFORX[3]) {
            printAwaveArray(NFORX[3], xFuse4, "FUSE4 X" , fp);
            printAwaveArray(NFORX[3], zFuse4, "FUSE4 Z" , fp);
            printAwaveArray(NFORX[3], aFuse4, "FUSE4 AREA" , fp);
        }

//        if (NFORX[0]) {
//            printAwaveArray(NFORX[0], xFuse1, "FUSE1 X" , fp);
//            printAwaveArray(NFORX[0], zFuse1, "FUSE1 Z" , fp);
//            printAwaveArray(NFORX[0], aFuse1, "FUSE1 AREA" , fp);
//        }
//
//        if (NFORX[1]) {
//            printAwaveArray(NFORX[1], xFuse1, "FUSE2 X" , fp);
//            printAwaveArray(NFORX[1], zFuse1, "FUSE2 Z" , fp);
//            printAwaveArray(NFORX[1], aFuse1, "FUSE2 AREA" , fp);
//        }
//
//        if (NFORX[2]) {
//            printAwaveArray(NFORX[2], xFuse1, "FUSE3 X" , fp);
//            printAwaveArray(NFORX[2], zFuse1, "FUSE3 Z" , fp);
//            printAwaveArray(NFORX[2], aFuse1, "FUSE3 AREA" , fp);
//        }
//        if (NFORX[3]) {
//            printAwaveArray(NFORX[3], xFuse1, "FUSE4 X" , fp);
//            printAwaveArray(NFORX[3], zFuse1, "FUSE4 Z" , fp);
//            printAwaveArray(NFORX[3], aFuse1, "FUSE4 AREA" , fp);
//        }
    }

    // PRINT POD DATA
    // NP - number of pods
    // NPODOR - number of points defining each pod (same for all pods)
    // locPod - location of Pods in nbody size array

    xPod = (double *) EG_alloc(NPODOR*sizeof(double));
    rPod = (double *) EG_alloc(NPODOR*sizeof(double));

    if ((xPod == NULL) || (rPod == NULL)) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Find the sections that define a pod
    for (n = 0; n < NP; n++) {

        for (i = 0; i < nbody; i++) {
            // Reuse locFuse because we are done with it
            locFuse[i] = 0;

            // Create a locPod the contains information for a single pod at a time
            if (locPod[i] == n+1) {
                locFuse[i] = 1;
            }
        }

        // Get awave values to define the pod
        status = defineAwavePod(nbody, surfaces, NPODOR, locFuse, xyzPod, xPod,
                                rPod);
        if (status != CAPS_SUCCESS) {
            printf(" awaveAIM/aimPreAnalysis defineAwavePod = %d\n", status);
        }

        /*@-bufferoverflowhigh@*/
        sprintf(tmpName,"POD %i",n);
        /*@+bufferoverflowhigh@*/

        printAwaveArray(3, xyzPod, tmpName , fp);
        printAwaveArray(NPODOR, xPod, "   X" , fp);
        printAwaveArray(NPODOR, rPod, "   RADIUS" , fp);
    }

    // PRINT FIN DATA
    if (J4) {

        // Only two sections can be printed for a fin
        // loop over all nbody and determine the first and last section location for each fin
        for (n = 0; n < NF; n++) {
            j = 0;
            for (i = 0; i < nbody; i++) {
                if (locFin[i] == n+1) {
                    iloc[j] = i;
                    j++;
                }
            }
            // Fin input is
            // x y z chord x y z chord for first and last section
            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[0],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[1],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[2],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].chordLength,7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[0],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[1],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[2],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].chordLength,7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            fprintf(fp,"                FIN %d \n",n);

            // ONLY shape of the first section is used to define the complete airfoil
            // no camber input is allowed
            printAwaveArray(surfaces[iloc[0]].ndiv, surfaces[iloc[0]].x,
                            "   CHORD 100(x/c)", fp);
            printAwaveArray(surfaces[iloc[0]].ndiv, surfaces[iloc[0]].halfThick,
                            "   HALF THICK", fp);
        }
    }

    // PRINT CANNARD DATA
    if (J5) {

        // Only two sections can be printed for a fin
        // loop over all nbody and determine the first and last section location for each fin
        for (n = 0; n < NCAN; n++) {
            j = 0;
            for (i = 0; i < nbody; i++) {
                if (locCan[i] == n+1) {
                    iloc[j] = i;
                    j++;
                }
            }
            // Fin input is
            // x y z chord x y z chord for first and last section
            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[0],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[1],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].xyz[2],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[0]].chordLength,7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[0],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[1],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].xyz[2],7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            tmpChar = convert_doubleToString(surfaces[iloc[j-1]].chordLength,7,0); fprintf(fp,"%s",tmpChar);
            if (tmpChar != NULL) EG_free(tmpChar);

            fprintf(fp,"                CANARD %d \n",n);

            // ONLY shape of the first section is used to define the complete airfoil
            // no camber input is allowed
            printAwaveArray(surfaces[iloc[0]].ndiv, surfaces[iloc[0]].x,
                            "   CHORD 100(x/c)", fp);
            printAwaveArray(surfaces[iloc[0]].ndiv, surfaces[iloc[0]].halfThick,
                            "   HALF THICK", fp);
        }
    }

    //  CASE CONTROL INPUT
    //  4-column right justified field

    /*
    Case Definition Input Data


    Case Control Card: Format(A4,11I4)  Integer input must be right
    justified in the indicated 4-column field.  As many Case Control
    Cards as desired may be input.  If NREST > 0, the Case Control
    Card is followed by a Restraint Card, but only one Restraint
    Card may be input per configuration.

    Variable   Columns    Description

    NCASE       1- 4     Case number or other 4 character descriptor

    MACH        5- 8     Mach number * 1000

    NX          9-12     Number of intervals on x axis

    NTHETA     13-16     Number of thetas

    NREST      17-20     Number of restraint points for drag minimization (maximum = 10)

    NCON       21-24     Configuration control
                         = 1, A new configuration follows this case.
                             (Title Card, Geometry Input Control Card,
                             etc.) This option is normally used when
                             only minor geometry changes are desired.
                         = 0, Otherwise (another case or nothing)

    ICYC       25-28     Number of Optimization cycles ( < 10)

    KKODE      29-32     Slope check control
                         = 0, Turn on slope checking
                         = 1, No slope checking

    JRST       33-36     Equivalent body data control
                         = 0, Compute equivalent body areas, drags etc.
                         = 1, Perform minimum calculations required for
                              wave drag (saves ~25% execution time)

     IALPH      37-40     Angle of attack x 100

     IUP1       41-44     Unit number for plot tape 1 which contains Mach
                          sliced area data for theta = -90 degrees
                          = 0, Plot tape not written
                          = Integer > 6, plot data on file PLOT1

     IUP2       45-48     Unit number for plot tape 2 which contains
                          equivalent body area data
                          = 0, Plot tape not written
                          = Integer > 6, plot data on file PLOT2


     Restraint Card: Format(10F7.0)  Input only if NREST > 0

     XREST(I)             X locations of fuselage restraint points.  The
                          fuselage cross-sectional area is held constant at
                          these locations during optimization.
     */

    // This AIM will only support analyis
    // Mulitple Mach and AoA entries can be made
    // A default Mach 1.0 design case is always created
    fprintf(fp,"DES 1000 100  32   0   0   0   0   0   0\n");

    printf("Number of Mach-Alpha cases = %d\n", inputs[inMach-1].length);

    if (inputs[0].length == 1) {

        // MACH
        Mach = 1000 * inputs[inMach-1].vals.real;
        fprintf(fp,"%4.0f", Mach); // case name
        fprintf(fp,"%4.0f", Mach); // Mach input

        // AOA
        AoA = 100 * inputs[inAlpha-1].vals.real;
        fprintf(fp," 100  32   0   0   0   0   1");
        fprintf(fp," %3.0f\n", AoA);

    } else {

        for (i = 0; i < inputs[0].length; i++) { // Multiple Mach, Angle of Attack pairs

            // MACH
            Mach = 1000 * inputs[inMach-1].vals.reals[i];
            fprintf(fp,"%4.0f", Mach); // case name
            fprintf(fp,"%4.0f", Mach); // Mach input

            // AOA
            AoA = 100 * inputs[inAlpha-1].vals.reals[i];
            fprintf(fp," 100  32   0   0   0   0   1");
            fprintf(fp," %3.0f\n", AoA);
        }
    }

    status = CAPS_SUCCESS;

// Free memory
cleanup:

    if (status != CAPS_SUCCESS)
      printf("Premature exit in awaveAIM preAnalysis status = %d\n", status);

    if (fp != NULL) fclose(fp);

    if (locPod  != NULL) EG_free(locPod);
    if (locCan  != NULL) EG_free(locCan);
    if (locFin  != NULL) EG_free(locFin);
    if (locFuse != NULL) EG_free(locFuse);
    if (locWing != NULL) EG_free(locWing);
    if (iloc    != NULL) EG_free(iloc);
    if (xFuse   != NULL) EG_free(xFuse);
    if (zFuse   != NULL) EG_free(zFuse);
    if (aFuse   != NULL) EG_free(aFuse);

    if (xPod != NULL) EG_free(xPod);
    if (rPod != NULL) EG_free(rPod);

    if (surfaces != NULL) {

        for (i = 0; i < nbody; i++) {
            (void) destroy_Section(&surfaces[i]);
        }

        EG_free(surfaces);
    }

    if (awaveType != NULL) EG_free(awaveType);

    if (nodes != NULL) EG_free(nodes);
    if (faces != NULL) EG_free(faces);

    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteAwave AIM Execution
   *
   * If auto execution is enabled when creating an awave AIM,
   * the AIM will execute awave just-in-time with the command line:
   *
   * \code{.sh}
   * awave awaveInput.txt > awaveOutput.txt
   * \endcode
   *
   * where preAnalysis generated the file "dealundoInput.txt" which contains the input information.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an awave AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or awave can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * awave.preAnalysis()
   *
   * print ("\n\nRunning......")
   * awave.system("awave awaveInput.txt > awaveOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * awave.postAnalysis()
   * \endcode
   */

  *state = 0;
  return aim_system(aimInfo, NULL, "awave awaveInput.txt > awaveOutput.txt");
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  // check the awave output file
  if (aim_isFile(aimInfo, "cdwave.txt") != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "awave execution did not produce cdwave.txt");
    return CAPS_EXECERR;
  }

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
#ifdef DEBUG
    printf(" awaveAIM/aimOutputs index = %d!\n", index);
#endif

    // Awave MAY NOT PRODUCE RESULTS FOR ALL INPUT MACH, AoA
    // ECHO THE SOVED INPUTS IN THE OUTPUT STRUCTURE

    /*! \page aimOutputsAwave AIM Outputs
	 * The main output for Awave is CDwave.  This reports wave drag coefficient with respect to the
	 * \ref aimInputsAwave given.  In addition, an echo of the Mach number and angle of attack inputs is
	 * provided.  This allows
	 * the user to ensure that the CDwave value matches the expected Mach, AoA input pair.  If a given pair
	 * does not execute then it will not appear in the results.  Thus, it is always good practice to do a sanity
	 * check using the echo of input values.
	 */

    if (index == outCDwave) {

        *aoname     = EG_strdup("CDwave");
        form->units = NULL;

        /*! \page aimOutputsAwave
         * - <B> CDwave = </B> Wave Drag Coefficient.
         */

    } else if (index == outMach) {
        *aoname     = EG_strdup("Mach");
        form->units = NULL;

        /*! \page aimOutputsAwave
         * - <B> MachOut = </B> Mach number.
         */

    } else if (index == outAlpha) {
        *aoname     = EG_strdup("Alpha");
        form->units = EG_strdup("degree");

        /*! \page aimOutputsAwave
         * - <B> Alpha = </B> Angle of attack (degree).
         */
    }

    form->type   = Double;
    form->lfixed = Change;
    form->sfixed = Fixed;
    form->dim    = Vector;
    form->nrow   = 1;
    form->ncol   = 1;
    form->vals.real = 0.0;
    form->vals.reals = NULL;

    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                  int index, capsValue *val)
{

    int status; // Function return status

    size_t linecap = 0;
    char   *str = "CASE", *valstr, *line = NULL;
    double tmp, Cd[100], Mach[100], AoA[100];
    int i, cnt = 0, DesCase = 0;

    FILE   *fp = NULL;

#ifdef DEBUG
    printf(" awaveAIM/aimCalcOutput instance = %d  index = %d!\n", inst, index);
#endif

    val->nrow = 1;
    val->ncol = 1;

    // ****************************************************************** //
    // PARSE OUTPUT DATA FILE //
    // ****************************************************************** //

    // Open the Awave output file
    fp = aim_fopen(aimInfo, "cdwave.txt", "r"); // EJA mod specific
    if (fp == NULL) {
        AIM_ERROR(aimInfo, " awaveAIM/aimCalcOutput Cannot open Output file!\n");
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Awave MAY NOT HAVE RAN ALL THE CASES
    // Comp MACH, AOA then assign Cd

    // Scan the file for the string
    valstr = NULL;

    while (getline(&line, &linecap, fp) >= 0) {

        if (line == NULL) continue;

        valstr = strstr(line, str);

        if (valstr != NULL) {

            // Skip the DES case Output information
            if (DesCase) {
#ifdef DEBUG
                printf("Awave OUTPUT FOUND :: %s",line);
#endif

                // Move down two lines MachOut
                (void) getline(&line, &linecap, fp);
                (void) getline(&line, &linecap, fp);
#ifdef DEBUG
                printf("\tLINE  %s", line);
#endif

                status = sscanf(line,"%*s %lf", &tmp);

                if (status == 0) {
                    AIM_ERROR(aimInfo, "Failed to parse: %s\n", line);
                    status = CAPS_IOERR;
                    goto cleanup;
                }

#ifdef DEBUG
                printf("\t\tMACH READ :: %8.6f\n", tmp);
#endif

                Mach[cnt] = tmp;

                // Move down 1 line AoAOut
                (void) getline(&line, &linecap, fp);
#ifdef DEBUG
                printf("\tLINE  %s", line);
#endif

                status = sscanf(line,"%*s %lf", &tmp);
                if (status == 0) {
                    AIM_ERROR(aimInfo, "Failed to parse: %s\n", line);
                    status = CAPS_IOERR;
                    goto cleanup;
                }

#ifdef DEBUG
                printf("\t\t AOA READ :: %8.6f\n", tmp);
#endif

                AoA[cnt] = tmp;

                // Move down 1 line CdWave
                (void) getline(&line, &linecap, fp);
#ifdef DEBUG
                printf("\tLINE  %s", line);
#endif

                status = sscanf(line,"%*s %lf", &tmp);
                if (status == 0) {
                    AIM_ERROR(aimInfo, "Failed to parse: %s\n", line);
                    status = CAPS_IOERR;
                    goto cleanup;
                }

#ifdef DEBUG
                printf("\t\t  CD READ :: %8.6f\n\n", tmp);
#endif

                Cd[cnt] = tmp;

                cnt++;
            }

            DesCase = 1;
        }
    }

    if (cnt < 1) {
        AIM_ERROR(aimInfo, "Cannot find '%s' in Output file!\n", str);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    val->nrow = cnt;
    val->ncol = 1;

    if (cnt == 1) {
        if (index == outCDwave) val->vals.real = Cd[0];
        if (index == outMach)   val->vals.real = Mach[0];
        if (index == outAlpha)  val->vals.real = AoA[0];

    } else if (cnt > 1) {

        AIM_ALLOC(val->vals.reals, val->nrow, double, aimInfo, status);

        for (i = 0; i < cnt; i++) {

            if (index == outCDwave) val->vals.reals[i] = Cd[i];
            if (index == outMach)   val->vals.reals[i] = Mach[i];
            if (index == outAlpha)  val->vals.reals[i] = AoA[i];

        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);
    if (line != NULL) free(line); // must use free

    return status;
}


// ********************** AIM Function Break *****************************
void aimCleanup(/*@unused@*/ void *aimStore)
{
#ifdef DEBUG
    printf(" awaveAIM/aimCleanup!\n");
#endif

    EG_free(aimStore);
}
