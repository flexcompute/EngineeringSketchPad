// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef __astrosCardS_H__
#define __astrosCardS_H__

#include "capsTypes.h"
#include "feaTypes.h" // Bring in FEA structures
#include "miscUtils.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * AEFACT
 * 
 * sid     Set identification number
 * D       Numbers (Real)
 * 
 * numD  Number of `D` values
 */
int astrosCard_aefact(FILE *fp, int *sid, int numD, double *D,
                       feaFileTypeEnum formatType);

/*
 * AERO
 * 
 * acsid   Aerodynamic coordinate system identification (Integer >= 0
 *         or Blank)
 * refc    Reference length (Real >= 0.0)
 * rhoref  Reference density (Real >= 0.0)
 */
int astrosCard_aero(FILE *fp, int *acsid, double *refc, double *rhoref,
                     feaFileTypeEnum formatType);

/*
 * AEROS
 * 
 * acsid   Aerodynamic coordinate system identification (Integer > 0,
 *         or blank)
 * rcsid   Reference coordinate system identification for rigid body motions
 *         (Integer > 0, or blank)
 * refc    Reference chord length (Real > 0.0, Default = 1.0)
 * refb    Reference span (Real > 0.0, Default = 1.0)
 * refs    Reference wing area (Real > 0.0, Default = 1.0)
 * gref    Reference grid point for stability derivative calculations
 *         (Integer > 0)
 * refd    Fuselage reference diameter (Real > 0.0, Default = 1.0)
 * refl    Fuselage reference length (Real > 0.0, Default = 1.0)
 */
int astrosCard_aeros(FILE *fp, int *acsid, int *rcsid, double *refc,
                      double *refb, double *refs, /*@null@*/ int *gref,
                     /*@null@*/ double *refd, /*@null@*/ double *refl,
                     feaFileTypeEnum formatType);

/*
 * AIRFOIL
 * 
 * acid    Associated aircraft component identification number referenced
 *         by a matching CAERO6 bulk data entry
 * cmpnt   Type of aircraft component (One of: "WING", "FIN", "CANARD")
 * cp      Coordinate system for airfoil (Integer > 0, or blank)
 * chord   Identification number of an AEFACT data entry containing a list
 *         of division points (in terms of percent chord) at which
 *         airfoil thickness and camber data are specified
 * usothk  Identification number of an AEFACT data entry defining either the
 *         upper surface ordinates in percent chord if `lso` is not blank,
 *         or the half thicknesses about the camber ordinates if `cam` is
 *         not blank (Integer > 0, or blank)
 * lso     Identification number of an AEFACT data entry defining the lower
 *         surface ordinates in percent chord (Integer > 0, or blank)
 * cam     Identification number of an AEFACT data entry defining the mean
 *         line (camber line) ordinates in percent chord (Integer)
 * radius  Radius of leading edge in percent chord (Real >= 0.0)
 * x1y1z1  Location of the airfoil leading edge in coordinate system `cp`
 *         (Real, y1 >= 0.0)
 * x12     Airfoil chord length in x-axis coordinate of system `cp`
 *         (Real > 0.0 or blank)
 * ipanel  Identification number of an AEFACT data containing a list of
 *         chord wise cuts in percent chord for wing paneling (Integer > 0,
 *         or blank)
 */
int astrosCard_airfoil(FILE *fp, int *acid, char *cmpnt, /*@null@*/ int *cp,
                       int *chord, int *usothk, int *lso, /*@null@*/ int *cam,
                       /*@null@*/ double *radius, double x1y1z1[3], double *x12,
                       /*@null@*/ int *ipanel, feaFileTypeEnum formatType);

/*
 * CAERO6
 * 
 * acid    Component identification number (Integer > 0)
 * cmpnt   Aircraft component (One of: "WING", "FIN", "CANARD")
 * cp      Coordinate system (Integer >= 0, or blank)
 * igrp    Group number for this component (Integer > 0)
 * lchord  Identification number of an AEFACT data entry containing a list
 *         of division points in percent chord for chord-wise boxes for
 *         the aerodynamic surface. If `lchord`
 * lspan
 */
int astrosCard_caero6(FILE *fp, int *acid, char *cmpnt, /*@null@*/ int *cp,
                      int *igrp, int *lchord, int *lspan,
                      feaFileTypeEnum formatType);

/*
 * CELAS2
 * 
 * eid     Element identification number (Integer > 0)
 * k       The value of the scalar spring (Real > 0.0)
 * G       Geometric grid point identification numbers (Integer >= 0)
 * C       Component numbers (6 >= Integer >= 0)
 * ge      Damping coefficient (Real >= 0.0)
 * s       Stress coefficient (Real >= 0.0)
 * tmin    Minimum value for design (Real)
 * tmax    Maximum value for design (Real)
 */
int astrosCard_celas2(FILE *fp, int *eid, double *k, int G[2], int C[2],
                      double *ge, double *s, /*@null@*/ double *tmin,
                      /*@null@*/ double *tmax, feaFileTypeEnum formatType);

/*
 * CIHEX1
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PIHEX property entry (Integer > 0,
 *         Default = `eid`)
 * G       Grid point identification numbers (Integer > 0)
 */
int astrosCard_cihex1(FILE *fp, int *eid, int *pid, int G[8],
                       feaFileTypeEnum formatType);

/*
 * CMASS2
 * 
 * eid     Element identification number (Integer > 0)
 * m       The value of the scalar mass (Real)
 * G       Geometric grid point identification number
 * C       Component number (6 >= Integer >= 0)
 * tmin    The minimum mass value for design (Real)
 * tmax    The maximum mass value for design (Real)
 */
int astrosCard_cmass2(FILE *fp, int *eid, double *m, int G[2], int C[2],
                      /*@null@*/ double *tmin, /*@null@*/ double *tmax,
                      feaFileTypeEnum formatType);
/*
 * CONM2
 * 
 * eid     Element identification number (Integer > 0)
 * g       Grid point identification number (Integer > 0)
 * cid     Coordinate system identification number (Integer >= -1)
 *         A `cid` of -1 allows the user to input `X` as the center
 *         of gravity location in the basic coordinate system. A `cid`
 *         of 0 implies the basic coordinate system
 * m       Mass value (Real)
 * X       Offset distances from the grid point to the center of gravity
 *         of the mass of the coordinate system defined by `cid`, unless
 *         `cid` = -1, in which case `X` values are the coordinates of the
 *         center of gravity of the mass in the basic coordinate system
 *         (Real, or blank)
 * I       Mass moments of inertia measured at the mass e.g., in
 *         coordinate system defined by `cid` (Real, or blank)
 * tmin    The minimum mass values for design (Real, or blank)
 * tmax    The maximum mass values for design (Real, or blank)
 */
int astrosCard_conm2(FILE *fp, int *eid, int *g, int *cid, double *m,
                     double X[3], double I[6], /*@null@*/ double *tmin,
                     /*@null@*/double *tmax, feaFileTypeEnum formatType);

/*
 * CQDMEM1
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PQDMEM1 property entry (Integer > 0,
 *         Default = `eid`)
 * G       Grid point identification numbers of connection points 
 *         (Integer > 0)
 * tm      Material property orientation specification (Real or Integer)
 * tmax    Maximum allowable element thickness in design (Real > 0.0, or blank)
 * 
 * tmType  The type of `tm` entry (1: Integer or 2: Double)
 */
int astrosCard_cqdmem1(FILE *fp, int *eid, int *pid, int G[4], int tmType,
                       /*@null@*/ void *tm, /*@null@*/ double *tmax,
                       feaFileTypeEnum formatType);

/*
 * CQUAD4
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PSHELL or PCOMPi entry (Integer > 0,
 *         Default = `eid`)
 * G       Grid point identification numbers of connection points
 *         (Integer > 0)
 * tm      Material property orientation specification (Real or blank,
 *         or 0 <= Integer < 1,000,000)
 * zoff    Offset of the element reference plane from the plane of grid
 *         points (Real, or blank)
 * tmax    Maximum allowable element thickness in design (Real > 0.0)
 * T       Membrane thicknesses of elements at grid points `G` (Real,
 *         or blank)
 * 
 * tmType  The type of `tm` entry (1: Integer or 2: Double)
 */
int astrosCard_cquad4(FILE *fp, int *eid, int *pid, int G[4], int tmType,
                      /*@null@*/ void *tm, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmax, /*@null@*/ double *T,
                      feaFileTypeEnum formatType);

/*
 * CROD
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PROD property entry (Integer > 0,
 *         Default = `eid`)
 * G       Grid point identification numbers of connection points 
 *         (Integer > 0)
 * tmax    Maximum allowable rod area in design (Real > 0.0, or blank)
 */
int astrosCard_crod(FILE *fp, int *eid, int *pid, int G[2],
                    /*@null@*/ double *tmax, feaFileTypeEnum formatType);

/*
 * CSHEAR
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PSHEAR property entry (Integer > 0,
 *         Default = `eid`)
 * G       Grid point identification numbers of connection points 
 *         (Integer > 0)
 * tmax    Maximum allowable thickness in design (Real > 0.0, or blank)
 */
int astrosCard_cshear(FILE *fp, int *eid, int *pid, int G[4],
                      /*@null@*/ double *tmax, feaFileTypeEnum formatType);

/*
 * CTRIA3
 * 
 * eid     Element identification number (Integer > 0)
 * pid     Identification number of a PSHELL or PCOMPi propetry entry
 *         (Integer > 0, Default = `eid`)
 * G       Grid point identification numbers of connection points
 *         (Integer > 0)
 * tm      Material property orientation specification (Real or blank,
 *         or 0 <= Integer < 1,000,000)
 * zoff    Offset of the element reference plane from the plane of grid
 *         points (Real, or blank)
 * tmax    Maximum allowable element thickness in design (Real > 0.0)
 * T       Membrane thicknesses of elements at grid points `G` (Real,
 *         or blank)
 * 
 * tmType  The type of `tm` entry (1: Integer or 2: Double)
 */
int astrosCard_ctria3(FILE *fp, int *eid, int *pid, int G[3], int tmType,
                      /*@null@*/ void *tm, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmax, /*@null@*/ double *T,
                      feaFileTypeEnum formatType);

/* 
 * DCONFLT
 * 
 * sid     Constraint set identification (Integer > 0)
 * vtype   Nature of the velocity referred to in the table 
 *         (Either "TRUE" or "EQUIV", Default = "TRUE")
 * gfact   Constraint scaling flag (Real > 0.0, Default = 0.10)
 * v       Velocity values (Real >= 0.0)
 * gam     Required damping value (Real)
 * 
 * numVal  Number of `v` and `gam` values
 */
int astrosCard_dconflt(FILE *fp, int *sid, char *vtype, double *gfact,
                       int numVal, double *v, double *gam, 
                       feaFileTypeEnum formatType);
/* 
 * DCONFRQ
 * 
 * sid     Constraint set identification (Integer > 0)
 * mode    Modal number of the frequency to be constrained (Integer > 0)
 * ctype   Constraint type. (Either "UPPER" or "LOWER", Default = "LOWER")
 * frqall  Frequency constraint in Hz (Real)
 */
int astrosCard_dconfrq(FILE *fp, int *sid, int *mode, char *ctype, 
                       double *frqall, feaFileTypeEnum formatType);

/*
 * DCONTWP
 * 
 * sid     Stress constraint set identification (Integer > 0)
 * xt      Tensile stress limit in the longitudinal direction (Real > 0.0)
 * xc      Compressive stress limit in the longitudinal direction
 *         (Real > 0.0)
 * yt      Tensile stress limit in the transverse direction (Real > 0.0)
 * yc      Compressive stress limit in the transverse direction
 *         (Real > 0.0)
 * ss      Shear stress limit for in-plane stress (Real > 0.0)
 * f12     Tsai-Wu interaction term (Real)
 * ptype   Property type (One of: "PQDMEM1", "PTRMEM", "PSHELL", "PCOMP",
 *         "PCOMP1", "PCOMP2")
 * layrnum The layer number of a composite element (Integer > 0 or blank)
 * PID     Property identification numbers (Integer > 0)
 * 
 * numPID  Number of `PID` values
 */
int astrosCard_dcontwp(FILE *fp, int *sid, double *xt, double *xc,
                        double *yt, double *yc, double *ss, double *f12, 
                        char *ptype, int *layrnum, int numPID, int *PID,
                        feaFileTypeEnum formatType);

/*
 * DCONVMP
 * 
 * sid     Stress constraint set identification (Integer > 0)
 * st      Tensile stress limit (Real > 0.0 or blank)
 * sc      Compressive stress limit (Real, Default = `st`)
 * ss      Shear stress limit (Real > 0.0, or blank)
 * ptype   Property type (One of: "PBAR", "PROD", "PSHEAR", "PQDMEM1",
 *         "PTRMEM", "PSHELL", "PCOMP", "PCOMP1", "PCOMP2")
 * layrnum The layer number of a composite element (Integer > 0 or blank)
 * PID     Property identification numbers (Integer > 0)
 * 
 * numPID  Number of `PID` values
 */
int astrosCard_dconvmp(FILE *fp, int *sid, double *st, /*@null@*/ double *sc,
                       /*@null@*/ double *ss, char *ptype,
                       /*@null@*/ int *layrnum, int numPID, int *PID,
                       feaFileTypeEnum formatType);

/*
 * DESVARP
 * 
 * dvid    Design variable identification (Integer > 0)
 * linkid  Link identification number referring to ELIST, ELISTM and/or
 *         PLIST, PLISTM entries (Integer > 0, Default = `dvid`)
 * vmin    Minimum allowable value of the design variable (Real >= 0.0,
 *         Default = 0.001)
 * vmax    Maximum allowable value of the design variable (Real >= 0.0,
 *         Default = 1000.0) 
 * vinit   Initial value of the design variable (Real, Default = 1.0)
 * layrnum Layer number if referencing a single layer of composite element(s)
 *         (Integer > 0 or blank)
 * layrlst Set identification number of PLYLIST entries specifying a set of
 *         composite layers to be linked (Integer > 0 or blank)
 * label   Optional user supplied label to define the design variable
 */
int astrosCard_desvarp(FILE *fp, int *dvid, int *linkid, double *vmin,
                       double *vmax, double *vinit, /*@null@*/ int *layrnum,
                       /*@null@*/ int *layrlst, /*@null@*/ char *label,
                       feaFileTypeEnum formatType);

// int astrosCard_dvgrid(FILE *fp, feaFileTypeEnum formatType);

/*
 * EIGR
 * 
 * sid     Set identification number (Integer > 0)
 * method  Method of eigenvalue extraction. Available options:
 *           - [ASTROS 11] INV, GIV
 *           - [ASTROS 12] SINV, GIV, MGIV, FEER
 *           - [ASTROS 20] GIV, MGIV
 * f1, f2  Frequency range of interest
 * ne      Estimate of number of roots in range (only used if `method` = sinv)
 * nd      Desired number of roots or eigenvectors/eigenvalues
 * e       [ASTROS 11, 12] Convergence test (Real, Default = 1e-6)
 *         [ASTROS 20] Mass orthogonality test parameter. 
 *         (Real > 0.0, Default = 1e-10) 
 * norm    Method for eigenvector normalization (one of "MASS", "MAX", or "POINT")
 * gid     Grid or scalar point identification number (Integer > 0)
 * dof     Component number (Integer 1-6)
 */
int astrosCard_eigr(FILE *fp, int *sid, char *method, double *f1,
                    double *f2, int *ne, int *nd, /*@null@*/ double *e,
                    char *norm, int *gid, int *dof, feaFileTypeEnum formatType);

/*
 * FLFACT
 * 
 * sid     Set identification number (Integer > 0)
 * F       Aerdynamic factor values (Real)
 * 
 * numF    Number of `F` values
 */
int astrosCard_flfact(FILE *fp, int *sid, int numF, double *F,
                       feaFileTypeEnum formatType);

/*
 * FLUTTER
 * 
 * sid     Set identification number (Integer > 0)
 * method  FLUTTER analysis method (Either "PK" or "PKIT", Default = "PK")
 * dens    Identification number of an FLFACT set specifying density ratios
 *         to be used in FLUTTER analysis (Integer > 0)
 * mach    Mach number to be used in the FLUTTER analysis (Real >= 0.0)
 * vel     Identification number of an FLFACT set specifying velocities
 *         to used in the FLUTTER analysis (Integer > 0)
 * mlist   Identification number of a SET1 set specifying a list of normal
 *         modes to be omitted from the FLUTTER analysis (Integer > 0,
 *         or blank)
 * klist   Identification number of an FLFACT set specifying a list of hard
 *         point reduced frequencies for the given Mach number for use in
 *         the FLUTTER analysis (Integer >= 0, or blank)
 * effid   Identification number of a CONEFFF set specifying control surface
 *         effectiveness values (Integer >= 0, or blank)
 * symxz   XZ-symmetry flag associated with the aerodynamics (Integer, 
 *         one of: +1, 0 or blank, -1)
 * symxy   XY-symmetry flag associated with the aerodynamics (Integer, 
 *         one of: +1, 0 or blank, -1)
 * eps     Convergence parameter for FLUTTER eigenvalue (Real, Default = 10-5)
 * curfit  Type of curve fit to be used in the PK FLUTTER analysis. (One of:
 *         "LINEAR", "QUAD", "CUBIC", "ORIG", Default = "LINEAR")
 * nroot   Requests that only the first `nroot` eigenvalues be found
 *         (Integer or blank)
 * vtype   Input velocities are in units of "TRUE" or "EQUIV" speed
 * gflut   The damping a mode must exceed to be considered a flutter
 *         crossing (Real >= 0.0, Default = 0.0)
 * gfilter The damping a mode must attain to be considered stable
 *         before a flutter crossing (Real, Default = 0.0)
 */
int astrosCard_flutter(FILE *fp, int *sid, char *method, int *dens,
                       double *mach, int *vel, /*@null@*/ int *mlist,
                       /*@null@*/ int *klist, /*@null@*/ int *effid,
                       int *symxz, int *symxy, /*@null@*/ double *eps,
                       /*@null@*/ char *curfit, /*@null@*/ int *nroot,
                       /*@null@*/ char *vtype, /*@null@*/ double *gflut,
                       /*@null@*/ double *gfilter, feaFileTypeEnum formatType);

/*
 * FORCE
 * 
 * sid     Load set identification number (Integer > 0)
 * g       Grid point identification number (Integer > 0)
 * cid     Coordinate system identification number (Integer >= 0, Default = 0)
 * f       Scale factor (Real)
 * N       Components of a vector measured in the coordinate system
 *         defined by `cid` (Real; must have at least one nonzero component)
 */
int astrosCard_force(FILE *fp, int *sid, int *g, int *cid, double *f,
                      double N[3], feaFileTypeEnum formatType);

/*
 * GRAV
 * 
 * sid     Set identification number (Integer > 0)
 * cid     Coordinate system identification number (Integer >= 0)
 * g       Gravity vector scale factor (Real != 0.0)
 * N       Gravity vector components (Real, at least one nonzero)
 */
int astrosCard_grav(FILE *fp, int *sid, int *cid, double *g, double N[3],
                     feaFileTypeEnum formatType);

/*
 * GRID
 * 
 * id      Grid point identification number (Integer > 0)
 * cp      Identification number of coordinate system in which the
 *         location of the grid points is define (Integer > 0, or blank)
 * xyz     Location of the grid point in the coordinate system (Real)
 * cd      Identification number of coordinate system in which displacements,
 *         degrees of freedom, constraints, and solution vectors are defined
 *         at the grid point (Integer > 0 or blank)
 * ps      Permanent single-point constraints associated with grid point
 *         (any of digits 1-6 with no embedded blanks) (Integer > 0 or blank)
 */
int astrosCard_grid(FILE *fp, int *id, /*@null@*/ int *cp,
                     double xyz[3], /*@null@*/ int *cd, /*@null@*/ int *ps,
                     feaFileTypeEnum formatType);

/*
 * MKAERO1
 * 
 * symxz   Symmetry flag (Integer)
 * symxy   Symmetry flag (Integer)
 * M       List of from 1 to 6 Mach numbers (Real >= 0.0 or blank)
 * K       List of from 1 to 8 reducted frequencies (Real >= 0.0 or blank)
 * 
 * numM    Number of Mach numbers
 * numK    Number of reduced frequencies
 */
int astrosCard_mkaero1(FILE *fp, int *symxz, int *symxy, int numM, double *,
                        int numK, double *K, feaFileTypeEnum formatType);

/*
 * MOMENT
 * 
 * sid     Load set identification number (Integer > 0)
 * g       Grid point identification number (Integer > 0)
 * cid     Coordinate system identification number (Integer >= 0)
 * m       Scale factor (Real)
 * N       Components of vector measured in coordinate system defined by CID (Real, 
 *         least one nonzero component)
 */
int astrosCard_moment(FILE *fp, int *sid, int *g, int *cid, double *m,
                       double N[3], feaFileTypeEnum formatType);

/*
 * PBAR
 * 
 * pid     Property identification number (Integer > 0)
 * mid     Material identification number (Integer > 0)
 * a       Area of bar cross-section (Real >= 0.0)
 * i1,i2   Area moments of inerta (Real >= 0)
 * j       Torsional constant (Real >= 0)
 * nsm     Nonstructural mass per unit length (Real >= 0.0)
 * tmin    The minimum cross-sectional area in design (Real, Default = 0.0001)
 * k1,k2   Area factor for shear (Real)
 * C,D,E,F Stress recovery coefficient pairs (Real or blank)
 * r12     Inertia linking terms for design (Real)
 * r22     ...
 * alpha   ...
 */
int astrosCard_pbar(FILE *fp, int *pid, int *mid, double *a, double *i1,
                    double *i2, double *j, double *nsm, /*@null@*/ double *tmin,
                    /*@null@*/ double *k1, /*@null@*/ double *k2,
                    /*@null@*/ double *C, /*@null@*/ double *D,
                    /*@null@*/ double *E, /*@null@*/ double *F,
                    /*@null@*/ double *r12, /*@null@*/ double *r22,
                    /*@null@*/ double *alpha, feaFileTypeEnum formatType);

/*
 * PCOMP
 * 
 * pid     Property identification number (Integer > 0)
 * z0      Offset of laminate lower surface from the element place (Real or blank)
 * nsm     Nonstructural mass per unit area (Real >= 0.0)
 * sbond   Allowable shear stress of the bonding material (Real >= 0.0)
 * ft      Failure theory ("HILL", "HOFF", "TSAI", "STRESS", or "STRAIN")
 * tmin    Minimum ply thickness for design (Real > 0.0 or blank)
 * lopt    Lamination generation option ("MEM" or blank)
 * MID     Material identification number of the i-th layer (Integer > 0)
 * T       Thickness of the i-th layer (Real > 0.0)
 * TH      Angle between the longitudinal direction of the fibers of the i-th
 *         layer and the material X-axis (Real)
 * SOUT    Stress output request for i-th layer, ("YES" or "NO", Default = "NO")
 * 
 * numLayers  Number of material layers
 * symmetricLaminate  If true, write symmetrical material layer fields
 */
int astrosCard_pcomp(FILE *fp, int *pid, /*@null@*/ double *z0,
                     /*@null@*/ double *nsm, double *sbond, /*@null@*/ char *ft,
                     /*@null@*/ double *tmin, /*@null@*/ char *lopt, int numLayers,
                     int *MID, double *T, double *TH, /*@null@*/ char **SOUT,
                     int symmetricLaminate, feaFileTypeEnum formatType);

/*
 * PIHEX
 * 
 * pid     Property identification number (Integer > 0)
 * mid     Material identification number (Integer > 0)
 * cid     Identification number of the coordinate system in which the material
 *         referenced by `mid` is defined
 * nip     Number of integration points along each edge of the element
 *         (Integer = 2, 3, 4, or blank)
 * ar      Maximum aspect ratio of the element (Real > 1.0 or blank)
 * alpha   Maximum angle in degrees between the normals of two subtriangles
 *         comprising a quadrilateral face (Real, 0.0 < alpha < 180.0, 
 *         or blank, Default = 45.0)
 * beta    Maximum angle in degrees between the vector connecting a corner
 *         point to an adjacent midside point and the vector connection that
 *         midside point and the other midside point or corner 
 *         (Real, 0.0 < beta < 180.0, or blank, Default = 45.0)
 */
int astrosCard_pihex(FILE *fp, int *pid, int *mid, /*@null@*/ int *cid,
                     /*@null@*/ int *nip, /*@null@*/ double *ar,
                     /*@null@*/ double *alpha, /*@null@*/ double *beta,
                     feaFileTypeEnum formatType);

/*
 * PLIST
 * 
 * linkid  Property list identifier (Integer > 0)
 * ptype   Property type associated with this list (e.g. PROD)
 * PID     Property entry identifications (Integer > 0, blank)
 */
int astrosCard_plist(FILE *fp, int *linkid, /*@null@*/ char *ptype, int numPID,
                     int *PID, feaFileTypeEnum formatType);

/*
 * PLOAD
 * 
 * sid    Load set identification number (Integer > 0)
 * p      Pressure (Real)
 * G      Grid point identification numbers (Integer > 0, G4 may be 0 or blank)
 * 
 * numG   Number of grid points defining the surface (3 or 4)
 */
int astrosCard_pload(FILE *fp, int *sid, double *P, int numG, int *G,
                      feaFileTypeEnum formatType);

/*
 * PLOAD4
 * 
 * lid    Load set identification number (Integer > 0)
 * eid    Element identification number (Integer > 0)
 * P      Pressure at the grid points defining the element surface (Real)
 * cid    Coordinate system identification number (Integer > 0, or blank)
 * V      Components of a vector in system `cid` that defines the direction
 *        of the grid point loads generated by the pressure (Real)
 * 
 * numP   Number of P values specified (1, 3, or 4)
 */
int astrosCard_pload4(FILE *fp, int *lid, int *eid, int numP, double *P,
                      /*@null@*/ int *cid, /*@null@*/ double *V,
                      feaFileTypeEnum formatType);

/*
 * PLYLIST
 * 
 * sid    Set identification number (Integer > 0)
 * P      List of ply numbers (Integer > 0)
 * 
 * numP   Number of plys
 */
int astrosCard_plylist(FILE *fp, int *sid, int numP, int *P,
                        feaFileTypeEnum formatType);

/*
 * PROD
 * 
 * pid    Property identification number (Integer > 0)
 * mid    Material identification number (Integer > 0)
 * a      Area of rod (Real >= 0, or blank)
 * j      Torsional constant (Real >= 0, or blank)
 * c      Coefficient to determine torsional stress (Real >= 0.0, or blank)
 * nsm    Nonstructural mass per unit length (Real >= 0, or blank)
 * tmin   Minimum rod area for design (Real > 0.0, or blank, Default = 0.0001)
 */
int astrosCard_prod(FILE *fp, int *pid, int *mid, double *a, double *j,
                     double *c, double *nsm, /*@null@*/ double *tmin,
                     feaFileTypeEnum formatType);

/*
 * PSHEAR
 */
int astrosCard_pshear(FILE *fp, int *pid, int *mid, double *t, 
                      /*@null@*/ double *nsm, /*@null@*/ double *tmin,
                      feaFileTypeEnum formatType);

/*
 * PQDMEM1
 * 
 * pid    Property identification number (Integer > 0)
 * mid    Material identification number (Integer > 0)
 * t      Thickness of membrane (Real >= 0, or blank)
 * nsm    Nonstructural mass per unit length (Real >= 0, or blank)
 * tmin   Minimum thickness for design (Real > 0.0, or blank, Default = 0.0001)
 */
int astrosCard_pqdmem1(FILE *fp, int *pid, int *mid, double *t, 
                       /*@null@*/ double *nsm, /*@null@*/ double *tmin,
                       feaFileTypeEnum formatType);

/*
 * PSHELL
 * 
 * pid    Property identification number (Integer > 0)
 * mid1   Material identification number for membrane (Integer > 0 or blank)
 * t      Default value for membrane thickness (Real > 0.0, or blank)
 * mid2   Material identification number for bending (Integer > 0 or blank)
 * i12t3  Bending stiffness parameter  (Real > 0.0, or blank, Default = 1.0)
 * mid3   Material identification number for transverse shear (Integer > 0,
 *        or blank, must be blank unless `mid2` > 0)
 * tst    Transverse shear thickness divided by membrane thickness
 *        (Real > 0.0 or blank, default = 0.833333)
 * nsm    Nonstructural mass per unit length (Real >= 0, or blank)
 * z1,z2  Fiber distances for stress computation (Real or blank)
 * mid4   Material identification number for membrane-bending coupling
 *        (Integer > 0 or blank, must blank unless `mid1` > 0 and `mid2` > 0,
 *        may not equal `mid1` or `mid2`)
 * mcsid  Identification number of material coordinate system 
 *        (Real or blank, or Integer >= 0)
 * scsid  Identification number of stress coordinate system 
 *        (Real or blank, or Integer >= 0)
 * zoff   Offset of the element reference plane from the plane of grid
 *        points (Real or blank, default = 0.0)
 * tmin   Minimum thickness for design (Real > 0.0 or blank, Default = 0.0001)
 * 
 * mcsidType  The type of `mcsid` entry (1: Integer or 2: Double)
 * scsidType  The type of `scsid` entry (1: Integer or 2: Double)
 */
int astrosCard_pshell(FILE *fp, int *pid, int *mid1, double *t,
                      /*@null@*/ int *mid2, /*@null@*/ double *i12t3,
                      /*@null@*/ int *mid3, /*@null@*/ double *tst,
                      /*@null@*/ double *nsm, /*@null@*/ double *z1,
                      /*@null@*/ double *z2, /*@null@*/ int *mid4, int mcsidType,
                      /*@null@*/ void *mcsid, int scsidType,
                      /*@null@*/ void *scsid, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmin, feaFileTypeEnum formatType);

/*
 * RBE2
 * 
 * setid   Multipoint constraint set identification number specifief in
 *         Solution Control (Integer > 0)
 * eid     Rigid body element identification number (Integer > 0)
 * gn      Grid point identification number at which all 6 independent
 *         degrees of freedom are assigned (Integer > 0)
 * cm      Component numbers of dependent degrees of freedom in the global
 *         coordinate system assigned by the elemnt at grid points `GM`
 *         (Integer > 0 or blank)
 * GM      Grid point identification numbers at which dependent degrees of
 *         freedom are assigned (Integer > 0)
 * 
 * numGM   Number of `GM` grid points
 */
int astrosCard_rbe2(FILE *fp, int *setid, int *eid, int *gn, int *cm,
                     int numGM, int *GM, feaFileTypeEnum formatType);

/*
 * RBE2
 * 
 * setid   Multipoint constraint set identification number specifief in
 *         Solution Control (Integer > 0)
 * eid     Rigid body element identification number (Integer > 0)
 * refg    Reference grid point identification number (Integer > 0)
 * refc    Component numbers of degrees of freedom in the global coordinate
 *         system that will be computed at `refg` (Integer > 0)
 * wt      Weighting factors (Real)
 * c       Component numbers of degrees of freedom in the global
 *         coordinate system assigned by the element at grid points `g`
 *         (Integer > 0 or blank)
 * g       Grid point identification number whose components `c` have 
 *         weighting factors `g` (Integer > 0)
 * gm      Grid point identification numbers with components in the m-set 
 *         (Integer > 0)
 * cm      Component numbers in the global coordinate system at 
 *         grid points `gm` which are placed in the m-set(Integer > 0)
 * 
 * numG    Number of `g`. `wt`, `c` values
 * numM    Number of `gm`, `cm` values
 */
int astrosCard_rbe3(FILE *fp, int *setid, int *eid, int *refg, int *refc,
                     int numG, double *wt, int *c, int *g,
                     int numM, /*@null@*/ int *gm, /*@null@*/ int *cm,
                     feaFileTypeEnum formatType);

// int astrosCard_rforce(FILE *fp, feaFileTypeEnum formatType);

/*
 * SPC
 * 
 * sid     Identification number of single-point constraint set
 * G       Grid or scalar point identification numbers (Integer > 0)
 * C       Component numbers of global coordinates (6 >= Integer >= 0;
 *         up to 6 unique digits may be placed in the field with no
 *         no embedded blanks)
 * D       Values of enforced displacements for all coordinates designed
 *         by G and C (Real)
 * 
 * numSPC  Number of `G`, `C`, `D` values
 */
int astrosCard_spc(FILE *fp, int *sid, int numSPC, int *G, int *C, double *D,
                    feaFileTypeEnum formatType);

/*
 * SPC1
 * 
 * sid     Identification number of single-point constraint set (Integer > 0)
 * c       Component number of global coordinate (6 >= Integer >= 0;
 *         up to 6 unique digits may be placed in the field with no
 *         no embedded blanks when grid points, must be null if scalar)
 * G       Grid or scalar point identification numbers (Integer > 0)
 * 
 * numSPC  Number of `G` values 
 */
int astrosCard_spc1(FILE *fp, int *sid, int *c, int numSPC, int *G,
                     feaFileTypeEnum formatType);

// int astrosCard_snorm(FILE *fp, feaFileTypeEnum formatType);

// int astrosCard_snormdt(FILE *fp, feaFileTypeEnum formatType);

/*
 * SPLINE1
 * 
 * eid     Element identification number (Integer > 0)
 * cp      Coordinate system defining the spline plane (Integer >= 0,
 *         or blank)
 * macroid Identification number of a CAEROi entry which defines the
 *         plane of spline (Integer > 0)
 * box1    First box whose motion is interpolated using this spline 
 *         (Integer > 0)
 * box2    Last box whose motion is interpolated using this spline
 *         (Integer > 0)
 * setg    Refers to a SETi entry which lists the structural grid points
 *         to which the spline is attached (Integer > 0)
 * dz      Linear attachment flexibility (Real >= 0)
 */
int astrosCard_spline1(FILE *fp, int *eid, /*@null@*/ int *cp, int *macroid,
                        int *box1, int *box2, int *setg, /*@null@*/ double *dz,
                        feaFileTypeEnum formatType);

/*
 * SUPORT
 * 
 * setid   Solution control SUPPORT set identification number (Integer > 0)
 * ID      Grid or scalar point identification numbers (Integer > 0)
 * C       Component numbers (zero or blank for scalar points; any unique
 *         combination of the digits 1 through 6 for grid points)
 * 
 * numSUP  Number of `ID`, `C` values
 */
int astrosCard_suport(FILE *fp, int *setid, int numSUP, int *ID, int *C,
                       feaFileTypeEnum formatType);

/*
 * TEMP
 * 
 * sid     Temperature set identification number (Integer > 0)
 * G       Grid point identification number (Integer > 0)
 * T       Temperature (Real)
 * 
 * numT    Number of `G`, `T` values
 */
int astrosCard_temp(FILE *fp, int *sid, int numT, int *G, double *T,
                     feaFileTypeEnum formatType);

/*
 * TEMPD
 * 
 * SID     Temperature set identification numbers (Integer > 0)
 * T       Temperatures (Real)
 * 
 * numT    Number of `SID`, `T` values
 */
int astrosCard_tempd(FILE *fp, int numT, int *SID, double *T,
                      feaFileTypeEnum formatType);

/*
 * TRIM parameter struct required to use `astrosCard_trim` function
 */
typedef struct {

    // Label defining aerodynamic trim parameter
    char *label;

    // Magnitude of the specified trim parameter
    double value;

    // Where trim parameter value is FREE
    int isFree;

} astrosCardTrimParamStruct;

// initiate astrosCardTrimParamStruct
int astrosCard_initiateTrimParam(astrosCardTrimParamStruct *param);

// destroy astrosCardTrimParamStruct
int astrosCard_destroyTrimParam(astrosCardTrimParamStruct *param);

/*
 * TRIM
 * 
 * trimid  Trim set identification number (Integer > 0)
 * mach    Mach number (Real >= 0.0)
 * qdp     Dynamic pressure (Real > 0.0)
 * trmtyp  Type of trim required (Character or blank)
 *           - blank   SUPORT controlled trim
 *           - ROLL    Axisymmetric roll trim (1 DOF)
 *           - LIFT    Symmetric trim of lift forces (1 DOF)
 *           - PITCH   Symmetric trim of lift and pitching moment (2 DOF)
 * effid   Identification number of CONEFFS entries which modify control
 *         surface effectiveness values (Integer >= 0, or blank)
 * vo      True velocity (Real > 0.0, or blank) 
 * 
 * numParam  Number of trim parameters
 * param     Trim parameters
 */
int astrosCard_trim(FILE *fp, int *trimid, double *mach, double *qdp,
                    /*@null@*/ char *trmtyp, /*@null@*/ int *effid,
                    /*@null@*/ double *vo, int numParam,
                    astrosCardTrimParamStruct *param, feaFileTypeEnum formatType);


#ifdef __cplusplus
}
#endif

#endif // __astrosCardS_H__
