/*
 ************************************************************************
 *                                                                      *
 * udpParsec -- udp file to generate a parsec airfoil (w/sensitivities) *
 *                                                                      *
 *            Written by Marshall Galbraith @ MIT                       *
 *            Patterned after .c code written by John Dannenhoffer      *
 *                                                @ Syracuse University *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2020
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

/* set to 1 for arc-lenght knots, -1 for equally spaced knots */
#define KNOTS -1

/* Surrealized spline fits */
#include "egadsSplineFit.h"
#include "egads_dot.h"

extern "C" {
#define NUMUDPARGS 5
#include "udpUtilities.h"
} /* extern "C" */

/* shorthands for accessing argument values and velocities */
#define YTE_VAL(   IUDP  ) ((double *) (udps[IUDP].arg[0].val))[0]
#define YTE_DOT(   IUDP  ) ((double *) (udps[IUDP].arg[0].dot))[0]
#define POLY_VAL(  IUDP,I) ((double *) (udps[IUDP].arg[1].val))[I]
#define POLY_DOT(  IUDP,I) ((double *) (udps[IUDP].arg[1].dot))[I]
#define PARAM_VAL( IUDP,I) ((double *) (udps[IUDP].arg[2].val))[I]
#define PARAM_DOT( IUDP,I) ((double *) (udps[IUDP].arg[2].dot))[I]
#define MEANLINE(  IUDP  ) ((int    *) (udps[IUDP].arg[3].val))[0]
#define ZTAIL_VAL( IUDP,I) (udps[IUDP].arg[4].size == 1 ? 0 : ((double *) (udps[IUDP].arg[4].val))[I])
#define ZTAIL_DOT( IUDP,I) (udps[IUDP].arg[4].size == 1 ? 0 : ((double *) (udps[IUDP].arg[4].dot))[I])

#define YTE_SURREAL(   IUDP  ) SurrealS<1>( YTE_VAL(  IUDP  ), YTE_DOT(  IUDP  ) )
#define POLY_SURREAL(  IUDP,I) SurrealS<1>( POLY_VAL( IUDP,I), POLY_DOT( IUDP,I) )
#define PARAM_SURREAL( IUDP,I) SurrealS<1>( PARAM_VAL(IUDP,I), PARAM_DOT(IUDP,I) )
#define ZTAIL_SURREAL( IUDP,I) SurrealS<1>( ZTAIL_VAL(IUDP,I), ZTAIL_DOT(IUDP,I) )

/* data about possible arguments */
static const char  *argNames[NUMUDPARGS] = {"yte",      "poly",       "param",     "meanline", "ztail",     };
static       int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, ATTRINT,    ATTRREALSEN, };
static       int    argIdefs[NUMUDPARGS] = {0,           0,           0,           0,          0,           };
static       double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,          0.,         0.,          };

struct udpDotCache_T
{
  int    sharpte;
  double yte_dot;
  double *poly_dot;
  double *param_dot;
  double ztail_dot[2];
} ;

#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data)
{
  udpDotCache_T *cache = (udpDotCache_T*)data;
  EG_free(cache->poly_dot);
  EG_free(cache->param_dot);
  EG_free(cache);
  return EGADS_SUCCESS;
}

extern "C" {
/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
} /* extern "C" */

/* functions for accessing argument values and surreal velocities
 * the compiler should inline these and they will go away in the shared library
 *
 * the function prototypes are declared, and then specializations are used to
 * return the appropriate macro
 */
template<class T> T    YTE(int iudp);
template<> double      YTE< double      >( int iudp ) { return YTE_VAL(iudp); }
template<> SurrealS<1> YTE< SurrealS<1> >( int iudp ) { return YTE_SURREAL(iudp); }

template<class T> T    POLY(int iudp, int i);
template<> double      POLY< double      >( int iudp, int i ) { return POLY_VAL(iudp,i); }
template<> SurrealS<1> POLY< SurrealS<1> >( int iudp, int i ) { return POLY_SURREAL(iudp,i); }

template<class T> T    PARAM(int iudp, int i);
template<> double      PARAM< double      >( int iudp, int i ) { return PARAM_VAL(iudp,i); }
template<> SurrealS<1> PARAM< SurrealS<1> >( int iudp, int i ) { return PARAM_SURREAL(iudp,i); }

template<class T> T    ZTAIL(int iudp, int i);
template<> double      ZTAIL< double      >( int iudp, int i ) { return ZTAIL_VAL(iudp,i); }
template<> SurrealS<1> ZTAIL< SurrealS<1> >( int iudp, int i ) { return ZTAIL_SURREAL(iudp,i); }

/* generating a FaceBody airfoil */
static int
buildFaceBodyAirfoil(ego context, int iudp, ego  *ebody);

/* computing the FaceBody airfoil sensitivity */
static int
sensFaceBodyAirfoil(ego ebody, int iudp, int npts, int *header, SurrealS<1> *pts, SurrealS<1> *rdata);

/* generating a WireBody meanline airfoil */
static int
buildWireBodyMeanline(ego context, int iudp, ego  *ebody);

/* computing the WireBody meanline sensitivity */
static int
sensWireBodyMeanline(ego ebody, int iudp, int npts, int *header, SurrealS<1> *pts, SurrealS<1> *rdata);

/* constructs a B-spline of the airfoil (or meanline) */
template<class T>
static int
parsecSplineFit(int iudp, int& npts, T **pts, int* header, T **rdata);

/* constructs the points used to for the airfoil (or meanline) B-spline */
template<class T>
static int
parsecPoints(int iudp, int& npts, T **pts);

/* constructs the top and bottom parsec polynomials */
template<class T>
static int
parsecPolyCoeff(int iudp, int& npoly, T **polyTop, T **polyBot);

/* function for computing polynomial fit */
template<class T>
static int
polyfit(const T& lefac, const T& yte, const T& x, const T& y, const T& c, const T& t, T poly[]);

/* function for evaluating the parsec polynomial */
template<class T>
static T
parsec(const double x, const T poly[], const int npoly);

/* matrix solve support routine for polyfit */
template<class T>
static int
matsol(T A[], T b[], int n, T x[]);

#define EPS12 1.0e-12  /* used to detect a singular matrix */
#define EPS06 1.0e-06

/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

extern "C" int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     i, iudp, npoly, nparam;
    double  rle, xtop, xbot;
    udpDotCache_T *cache;

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    npoly  = udps[0].arg[1].size;
    nparam = udps[0].arg[2].size;

    if (udps[0].arg[0].size > 1) {
        printf(" udpParsec.udpExecute: yte should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (npoly > 1 && nparam > 1) {
        printf(" udpParsec.udpExecute: poly and param cannot both be set\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (nparam > 1 && nparam != 9) {
        printf(" udpParsec.udpExecute: there should be 9 elements in Param\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (npoly > 1 && npoly%2 != 0) {
        printf(" udpParsec.udpExecute: there should be an even number of elements in Poly\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (npoly <= 1 && nparam <= 1) {
        printf(" udpParsec.udpExecute: neither poly nor param was set\n");
        status  = EGADS_NODATA;
        goto cleanup;

    } else if (!((udps[0].arg[4].size == 1 && ZTAIL_VAL(0,0) == 0) || udps[0].arg[4].size == 2)) {
        printf(" udpExecute: ztail should contain 0 or 2 values (upper,lower)\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* more input checks */
    if (nparam > 1) {
        /* don't use iudp here because it has not been cached yet. */
        rle  = PARAM_VAL(0,0);

        xtop = PARAM_VAL(0,1);
        xbot = PARAM_VAL(0,5);

        if (rle <= 0) {
            printf(" udpParsec.udpExecute: rle = %f <= 0\n", rle);
            status  = EGADS_RANGERR;
            goto cleanup;
        } else if (xtop <= 0) {
            printf(" udpParsec.udpExecute: xtop = %f <= 0\n", xtop);
            status  = EGADS_RANGERR;
            goto cleanup;
        } else if (xbot <= 0) {
            printf(" udpParsec.udpExecute: xbot = %f <= 0\n", xbot);
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }


    /* Cache copy of arguments for future use.
     * This also increments numUdp and copies udps[0] to udps[numUdp].
     * Caching should only be performed after checking for valid inputs.
     */
    status = cacheUdp(NULL);
    if (status < 0) {
        printf(" udpParsec.udpExecute: problem caching arguments\n");
        goto cleanup;
    }

    /* Set iudp to be consistent with udpSensitivity */
    iudp = numUdp;


    /* setup the cache for tracking dot changes */
    cache = (udpDotCache_T*)EG_alloc(sizeof(udpDotCache_T));
    udps[iudp].data = (void*)cache;

    cache->sharpte   = 0;
    cache->yte_dot   = 0;
    cache->poly_dot  = NULL;
    cache->param_dot = NULL;

    cache->poly_dot = (double*)EG_alloc(udps[iudp].arg[1].size*sizeof(double));
    if (cache->poly_dot == NULL && udps[iudp].arg[1].size > 0) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < udps[iudp].arg[1].size; i++) {
      cache->poly_dot[i] = 0;
    }

    cache->param_dot = (double*)EG_alloc(udps[iudp].arg[2].size*sizeof(double));
    if (cache->param_dot == NULL && udps[iudp].arg[2].size > 0) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < udps[iudp].arg[2].size; i++) {
      cache->param_dot[i] = PARAM_DOT(iudp,i);
    }
    for (i = 0; i < 2; i++) {
      cache->ztail_dot[i] = ZTAIL_DOT(iudp,i);
    }

#ifdef DEBUG
    printf("udpParsec.udpExecute\n");
    printf("yte(%d)      = %f\n", numUdp, YTE_VAL(     numUdp));
    for (int n = 0; n < npoly; n++) {
        printf("poly(%d,%d)  = %f\n", numUdp, n, POLY_VAL(numUdp,n) );
    }
    for (int n = 0; n < nparam; n++) {
        printf("param(%d,%d) = %f\n", numUdp, n, PARAM_VAL(numUdp,n) );
    }
    printf("meanline(%d) = %d\n", numUdp, MEANLINE(numUdp));
    printf("\n");
#endif

    /* if meanline[iudp] is zero, then create a FaceBody (the default) */
    if (MEANLINE(iudp) == 0) {

        status = buildFaceBodyAirfoil(context, iudp, ebody);
        if (status != EGADS_SUCCESS) goto cleanup;

    /* otherwise create a WireBody of the meanline */
    } else {

        status = buildWireBodyMeanline(context, iudp, ebody);
        if (status != EGADS_SUCCESS) goto cleanup;

    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }

    return status;
}


/*
 *****************************************************************************
 *                                                                           *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                           *
 *****************************************************************************
 */

extern "C" int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;
    int    npts, iudp, judp, i, ipnt, idotChange, nchild, stride, header[4];
    double point[18], point_dot[18];
    ego    eent, *echildren;
    SurrealS<1> *pts=NULL, *rdata=NULL;
    udpDotCache_T *cache;

    /* check that ebody matches one of the ebodys
     * this seems to be done for every udpSensitivity,
     * why isn't iudp just an argument to udpSensitivity?
     */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

#ifdef DEBUG
    int n, npoly, nparam;

    /* get arguments */
    npoly  = udps[iudp].arg[1].size;
    nparam = udps[iudp].arg[2].size;

    printf("udpParsec.udpSensitivity\n");
    printf("YTE_VAL(%d) = %f\n", iudp, YTE_VAL(iudp) );
    printf("YTE_DOT(%d) = %f\n", iudp, YTE_DOT(iudp) );
    for (n = 0; n < npoly; n++) {
        printf("POLY_VAL(%d,%d) = %f\n", iudp, n, POLY_VAL(iudp,n) );
    }
    for (n = 0; n < npoly; n++) {
        printf("POLY_DOT(%d,%d) = %f\n", iudp, n, POLY_DOT(iudp,n) );
    }
    for (n = 0; n < nparam; n++) {
        printf("PARAM_VAL(%d,%d) = %f\n", iudp, n, PARAM_VAL(iudp,n) );
    }
    for (n = 0; n < nparam; n++) {
        printf("PARAM_DOT(%d,%d) = %f\n", iudp, n, PARAM_DOT(iudp,n) );
    }
    printf("\n");
#endif


    /* check if velocities have changed */
    cache = (udpDotCache_T*)udps[iudp].data;
    idotChange = 0;
    if ( cache->yte_dot != YTE_DOT(iudp) ) {
        idotChange = 1;
    }
    for (i = 0; i < udps[iudp].arg[1].size; i++)
          if ( cache->poly_dot[i] != POLY_DOT(iudp,i) ) {
              idotChange = 1;
              break;
          }
    for (i = 0; i < udps[iudp].arg[2].size; i++)
        if ( cache->param_dot[i] != PARAM_DOT(iudp,i) ) {
            idotChange = 1;
            break;
        }
    for (i = 0; i < 2; i++)
        if ( cache->ztail_dot[i] != ZTAIL_DOT(iudp,i) ) {
            idotChange = 1;
            break;
        }

    /* build the sensitvity if needed */
    if (idotChange == 1) {

        cache->yte_dot = YTE_DOT(iudp);
        for (i = 0; i < udps[iudp].arg[1].size; i++) {
          cache->poly_dot[i] = POLY_DOT(iudp,i);
        }
        for (i = 0; i < udps[iudp].arg[2].size; i++) {
          cache->param_dot[i] = PARAM_DOT(iudp,i);
        }
        for (i = 0; i < 2; i++) {
          cache->ztail_dot[i] = ZTAIL_DOT(iudp,i);
        }

        /* get airfoil B-spline with it's sensitivities */
        status = parsecSplineFit(iudp, npts, &pts, header, &rdata);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* if meanline[iudp] is zero, then this is a FaceBody (the default) */
        if (MEANLINE(iudp) == 0) {

            status = sensFaceBodyAirfoil(ebody, iudp, npts, header, pts, rdata);
            if (status != EGADS_SUCCESS) goto cleanup;

        /* otherwise this is a WireBody of the meanline */
        } else {

            status = sensWireBodyMeanline(ebody, iudp, npts, header, pts, rdata);
            if (status != EGADS_SUCCESS) goto cleanup;

        }
    }


    /* find the ego entity */
    if (entType == OCSM_NODE) {
        status = EG_getBodyTopos(ebody, NULL, NODE, &nchild, &echildren);
        if (status != EGADS_SUCCESS) goto cleanup;

        stride = 0;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else if (entType == OCSM_EDGE) {
        status = EG_getBodyTopos(ebody, NULL, EDGE, &nchild, &echildren);
        if (status != EGADS_SUCCESS) goto cleanup;

        stride = 1;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else if (entType == OCSM_FACE) {
        status = EG_getBodyTopos(ebody, NULL, FACE, &nchild, &echildren);
        if (status != EGADS_SUCCESS) goto cleanup;

        stride = 2;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else {
        printf("udpSensitivity: bad entType=%d\n", entType);
        status = EGADS_GEOMERR;
        goto cleanup;
    }

    /* get the velocities from the entity */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        status = EG_evaluate_dot(eent, &(uvs[stride*ipnt]), NULL, point, point_dot);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* return the point velocity */
        vels[3*ipnt  ] = point_dot[0];
        vels[3*ipnt+1] = point_dot[1];
        vels[3*ipnt+2] = point_dot[2];
    }

cleanup:
    delete [] pts;
    EG_free( rdata );

    return status;
}


/*
 ***********************************************************************************
 *                                                                                 *
 *   buildFaceBodyAirfoil - builds a FaceBody of the airfoil geometry with EGADS   *
 *                                                                                 *
 *                          The B-spline airfoil approximation is constructed with *
 *                          parsecSplineFit, and the remainder of the function     *
 *                          constructs 2 node and 2 edge topologies as illustrated *
 *                          below:                                                 *
 *                                                                                 *
 *                                 edge 1                                          *
 *                              //======\\                                         *
 *     y                     //           \\                                       *
 *     ^           node 2  *       face     \\                                     *
 *     |                     \\============\\ *   node 1                           *
 *     |                          edge 2                                           *
 *     +-----> x                                                                   *
 *                                                                                 *
 *                          or with a blunt trailing edge with the topology:       *
 *                                                                                 *
 *                                 edge 1                                          *
 *                              //======\\                                         *
 *     y                     //            \\ *   node 1                           *
 *     ^           node 2  *       face       |   edge 3                           *
 *     |                     \\============\\ *   node 3                           *
 *     |                          edge 2                                           *
 *     +-----> x                                                                   *
 *                                                                                 *
 ***********************************************************************************
 */

static int
buildFaceBodyAirfoil(ego  context,      /* (in)  the EGADS context */
                     int  iudp,         /* (in)  udp index */
                     ego  *ebody)       /* (out) Body pointer */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int    i, sense[3], npts, nedge;
    double tle, data[18], tdata[2];
    double *rdata=NULL, *pts=NULL;
    ego    ecurve, enodes[2], eedges[3], enode1, enode2, enode3, eline, eloop, eplane, eface;
    int header[4];
    udpDotCache_T *cache;

    cache = (udpDotCache_T*)udps[iudp].data;

    /* get the B-spline and points used to construct the spline */
    status = parsecSplineFit(iudp, npts, &pts, header, &rdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* convert the B-spline to a CURVE */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rdata, &ecurve);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create Node at upper trailing edge (node 1)*/
    i = 0;
    data[0] = pts[3*i  ];
    data[1] = pts[3*i+1];
    data[2] = pts[3*i+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enode1));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* node at leading edge as a function of the spline */
    i = (npts - 1) / 2 + 3; /* leading edge index, with knot offset of 3 (cubic)*/
    tle = rdata[i];         /* leading edge t-value (should be very close to (0,0,0) */

    status = EG_evaluate(ecurve, &tle, data);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create node at leading edge (node 2)*/
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enode2));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make Edge for upper surface (edge 1)*/
    tdata[0] = 0;
    tdata[1] = tle;

    enodes[0] = enode1;
    enodes[1] = enode2;

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             tdata, 2, enodes, NULL, &(eedges[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create line segment at trailing edge */
    i = npts-1;
    data[0] = pts[3*i  ];
    data[1] = pts[3*i+1];
    data[2] = pts[3*i+2];
    data[3] = pts[0] - data[0];
    data[4] = pts[1] - data[1];
    data[5] = pts[2] - data[2];

    if (fabs(data[3]) > EPS06 || fabs(data[4]) > EPS06 || fabs(data[5]) > EPS06) {
        nedge = 3;
        cache->sharpte = 0;

        /* make the 3rd node on lower trailing edge */
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enode3));
        if (status != EGADS_SUCCESS) goto cleanup;

        /* make Edge for lower surface (edge 2)*/
        tdata[0] = tle;
        tdata[1] = 1;

        enodes[0] = enode2;
        enodes[1] = enode3;

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, enodes, NULL, &(eedges[1]));
        if (status != EGADS_SUCCESS) goto cleanup;

        /* make the trailing edge Line */
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* make Edge for the line */
        tdata[0] = 0;
        tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

        enodes[0] = enode3;
        enodes[1] = enode1;

        status = EG_makeTopology(context, eline, EDGE, TWONODE,
                                 tdata, 2, enodes, NULL, &(eedges[2]));
        if (status != EGADS_SUCCESS) goto cleanup;
    } else {
        nedge = 2;
        cache->sharpte = 1;

        /* make Edge for lower surface (edge 2)*/
        tdata[0] = tle;
        tdata[1] = 1;

        enodes[0] = enode2;
        enodes[1] = enode1;

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, enodes, NULL, &(eedges[1]));
        if (status != EGADS_SUCCESS) goto cleanup;
    }

    /* create loop of the two Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, nedge, eedges, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create a plane for the loop */
    data[0] = 0.;
    data[1] = 0.;
    data[2] = 0.;
    data[3] = 1.; data[4] = 0.; data[5] = 0.;
    data[6] = 0.; data[7] = 1.; data[8] = 0.;

    status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, data, &eplane);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create the face from the plane and the loop */
    status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                             NULL, 1, &eloop, sense, &eface);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create the face body */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
    if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
    delete [] pts;
    EG_free(rdata);

    return status;
}


/*
 ***********************************************************************************
 *                                                                                 *
 *   sensFaceBodyAirfoil - sensitivity of a FaceBody airfoil                       *
 *                                                                                 *
 ***********************************************************************************
 */

static int
sensFaceBodyAirfoil(ego ebody,          /* (in)  Body pointer */
                    int iudp,           /* (in)  udp index */
                    int npts,           /* (in)  number of points in the spline fit */
                    int *header,        /* (in)  spline header */
                    SurrealS<1> *pts,   /* (in)  spline points */
                    SurrealS<1> *rdata) /* (in)  spline data */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int         i, ipnt, nchild, nedge, nnode, oclass, mtype, *senses;
    double      data[18];
    SurrealS<1> tle, sdata[18];
    ego         eface, eref, eplane, eloop, *eedges, ecurve, eline, *enodes, *echildren;
    udpDotCache_T *cache;

    cache = (udpDotCache_T*)udps[iudp].data;

    /* get the face from the FACEBODY */
    status = EG_getTopology(ebody, &eref, &oclass, &mtype, data, &nchild, &echildren,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    eface = echildren[0];

    /* get the plane and the loop */
    status = EG_getTopology(eface, &eplane, &oclass, &mtype, data, &nchild, &echildren,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    eloop = echildren[0];

    /* get the edges from the loop */
    status = EG_getTopology(eloop, &eref, &oclass, &mtype, data, &nedge, &eedges,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the nodes and the curve from the first edge */
    status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, data, &nnode, &enodes,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the Curve */
    status = EG_setGeometry_dot(ecurve, CURVE, BSPLINE, header, rdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the Node at upper trailing edge */
    ipnt = 0;
    status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, &(pts[3*ipnt]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* node at leading edge as a function of the spline */
    i = (npts - 1) / 2 + 3; /* leading edge index, with knot offset of 3 (cubic)*/
    tle = rdata[i];         /* leading edge t-value (should be very close to (0,0,0) */

    status = EG_evaluate(ecurve, &tle, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the Node at leading edge */
    status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set Edge t-range sensitivity for upper surface */
    sdata[0] = 0;
    sdata[1] = tle;

    status = EG_setRange_dot(eedges[0], EDGE, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set Edge t-range sensitivity for lower surface */
    sdata[0] = sdata[1];
    sdata[1] = 1;

    status = EG_setRange_dot(eedges[1], EDGE, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;


    if (cache->sharpte == 0) {
        /* get trailing edge line and the lower trailing edge node from the 3rd edge */
        status = EG_getTopology(eedges[2], &eline, &oclass, &mtype, data, &nchild, &enodes,
                                &senses);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity of the Node at lower trailing edge */
        ipnt = npts - 1;
        status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, &(pts[3*ipnt]));
        if (status != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity of the line segment at trailing edge */
        sdata[0] = pts[3*ipnt  ];
        sdata[1] = pts[3*ipnt+1];
        sdata[2] = pts[3*ipnt+2];
        sdata[3] = pts[0] - sdata[0];
        sdata[4] = pts[1] - sdata[1];
        sdata[5] = pts[2] - sdata[2];

        status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, sdata);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* set Edge t-range sensitivity */
        sdata[0] = 0;
        sdata[1] = sqrt(sdata[3]*sdata[3] + sdata[4]*sdata[4] + sdata[5]*sdata[5]);

        status = EG_setRange_dot(eedges[2], EDGE, sdata);
        if (status != EGADS_SUCCESS) goto cleanup;
    }


    /* plane data */
    sdata[0] = 0.;
    sdata[1] = 0.;
    sdata[2] = 0.;
    sdata[3] = 1.; sdata[4] = 0.; sdata[5] = 0.;
    sdata[6] = 0.; sdata[7] = 1.; sdata[8] = 0.;

    /* set the sensitivity of the plane */
    status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;

cleanup:

    return status;
}


/*
 ***********************************************************************************
 *                                                                                 *
 *   buildWireBodyMeanline - builds a WireBody of the meanline geometry with EGADS *
 *                                                                                 *
 *                          The B-spline airfoil approximation is constructed with *
 *                          parsecSplineFit, and the remainder of the function     *
 *                          constructs 2 node and 1 edge topologies as illustrated *
 *                          below:                                                 *
 *                                                                                 *
 *                                 edge 1                                          *
 *                              //========\\                                       *
 *     y                     //             \\                                     *
 *     ^           node 1  *                  \\                                   *
 *     |                                       *  node 2                           *
 *     |                                                                           *
 *     +-----> x                                                                   *
 *                                                                                 *
 ***********************************************************************************
 */

static int
buildWireBodyMeanline(ego context,      /* (in)  the EGADS context */
                      int iudp,         /* (in)  iudp index */
                      ego *ebody)       /* (out) Body pointer */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int    ile, ite, sense[1], npts;
    double data[18], tdata[2], *pts=NULL, *rdata=NULL;
    ego    ecurve, enodes[2], eedges[1], eloop;

    /* create spline curve from LE to TE */
    int header[4];

    /* get the B-spline and points used to construct the spline */
    status = parsecSplineFit(iudp, npts, &pts, header, &rdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* the leading and trailing edge indexes (must be consistent with parsecPoints) */
    ile = (npts - 1) / 2;
    ite =  npts - 1;

    /* convert the B-spline to a CURVE */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rdata, &ecurve);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create Node at leading edge (node 1)*/
    data[0] = pts[3*ile  ];
    data[1] = pts[3*ile+1];
    data[2] = pts[3*ile+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create Node at trailing edge (node 2)*/
    data[0] = pts[3*ite  ];
    data[1] = pts[3*ite+1];
    data[2] = pts[3*ite+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* t-values at the leading and trailing edge nodes */
    tdata[0] = 0;
    tdata[1] = 1;

    /* make Edge for camberline (edge 1)*/
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create loop of the Edges */
    sense[0] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, OPEN,
                             NULL, 1, eedges, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create the WireBody (which will be returned) */
    status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                             NULL, 1, &eloop, NULL, ebody);
    if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
    delete [] pts;
    EG_free(rdata);

    return status;
}


/*
 ***********************************************************************************
 *                                                                                 *
 *   sensWireBodyMeanline - sensitivity of WireBody meanline                       *
 *                                                                                 *
 ***********************************************************************************
 */

static int
sensWireBodyMeanline(ego ebody,          /* (in)  Body pointer */
                     int iudp,           /* (in)  udp index */
                     int npts,           /* (in)  number of points in the spline fit */
                     int *header,
                     SurrealS<1> *pts,   /* (in)  spline points */
                     SurrealS<1> *rdata) /* (in)  spline data */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int         ile, ite, nchild, nedge, nnode, oclass, mtype, *senses;
    double      data[18];
    SurrealS<1> sdata[2];
    ego         eref, eloop, *eedges, ecurve, *enodes, *echildren;

    /* get the loop from the body */
    status = EG_getTopology(ebody, &eref, &oclass, &mtype, data, &nchild, &echildren,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    eloop = echildren[0];

    /* get the edges from the loop */
    status = EG_getTopology(eloop, &eref, &oclass, &mtype, data, &nedge, &eedges,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the nodes and the curve from the first edge */
    status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, data, &nnode, &enodes,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the Curve */
    status = EG_setGeometry_dot(ecurve, CURVE, BSPLINE, header, rdata);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* the leading and trailing edge indexes (must be consistent with parsecPoints) */
    ile = (npts - 1) / 2;
    ite =  npts - 1;

    /* set the sensitivity of the Node at trailing edge */
    status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, &(pts[3*ile]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the Node at leading edge */
    status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, &(pts[3*ite]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set Edge t-range sensitivity */
    sdata[0] = 0;
    sdata[1] = 1;

    status = EG_setRange_dot(eedges[0], EDGE, sdata);
    if (status != EGADS_SUCCESS) goto cleanup;

cleanup:

    return status;
}


/*
 **************************************************************************
 *                                                                        *
 *   parsecSplineFit - constructs a B-spline of the airfoil (or meanline) *
 *                                                                        *
 **************************************************************************
 */

template<class T>
static int
parsecSplineFit(int    iudp,         /* (in)  udp index */
                int&   npts,         /* (out) number of points in the spline fit */
                T      **pts,        /* (out) points used in the spline fit */
                int*   header,       /* (out) spline header data */
                T      **rdata)      /* (out) spline data */
{
    int status = EGADS_SUCCESS;      /* (out) return status */

    int     ile, nbspts;
    double  dxytol=1.0e-6; /* EG_spline1dFit tolerance */

    /* get airfoil (or meanline) points */
    status = parsecPoints(iudp, npts, pts);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* if meanline[iudp] is zero, then create a FaceBody (the default) */
    if (MEANLINE(iudp) == 0) {

        /* compute the spline fit
         * finite difference must use knots equally spaced (-npts in the argument below) in order for
         * the leading edge node sensitivity to be correct.
         * the t-value moves with arc-length based knots which causes problems
         */
        status = EG_spline1dFit<T>(0, KNOTS*npts, *pts, NULL, dxytol, header, rdata);
        if (status != EGADS_SUCCESS) goto cleanup;

    /* otherwise create a WireBody of the meanline */
    } else {

        /* modify the node count to be consistent with the meanline points */
        ile = (npts - 1) / 2;

        nbspts = ile + 1;

        /* compute the spline fit */
        status = EG_spline1dFit<T>(0, KNOTS*nbspts, *pts + 3*ile, NULL, dxytol, header, rdata);
        if (status != EGADS_SUCCESS) goto cleanup;
    }

cleanup:

    return status;
}


/*
 *****************************************************************************************
 *                                                                                       *
 *   parsecPoints - constructs the points used to for the airfoil (or meanline) B-spline *
 *                                                                                       *
 *****************************************************************************************
 */

template<class T>
static int
parsecPoints(int    iudp,         /* (in)  udp index */
             int&   npts,         /* (out) number of points */
             T      **ppts)       /* (out) airfoil (or meanline) points */
{
    int status = EGADS_SUCCESS;   /* (out) return status */

    int    i, j, npoly, ile, ite;
    double zeta, xx;
    T      yy, *polyTop=NULL, *polyBot=NULL, *pts=NULL;

    /* construct the parsec polynomial coefficients */
    status = parsecPolyCoeff(iudp, npoly, &polyTop, &polyBot);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* allocation required by Windows compiler */
    npts = 101;
    pts = new T[3*npts];
    if (pts == NULL) {
        status  = EGADS_MALLOC;
        goto cleanup;
    }

    /* points around airfoil (upper and then lower)
     */
    for (i = 0; i < npts; i++) {
        zeta = TWOPI * i / (npts-1);
        xx   = (1 + cos(zeta)) / 2;

        if (i < npts/2) {
            yy = parsec(xx, polyTop, npoly);

            pts[3*i  ] = xx;
            pts[3*i+1] = yy + ZTAIL<T>(iudp,0)*xx;
            pts[3*i+2] = 0;
        } else if (i == (npts - 1) / 2) {
            pts[3*i  ] = 0;
            pts[3*i+1] = 0;
            pts[3*i+2] = 0;
        } else {
            yy = parsec(xx, polyBot, npoly);

            pts[3*i  ] = xx;
            pts[3*i+1] = yy + ZTAIL<T>(iudp,1)*xx;
            pts[3*i+2] = 0;
        }
    }

    /* if meanline[iudp] is zero, then create a FaceBody (the default) */
    if (MEANLINE(iudp) == 0) {

        /* use the points as is */

    /* otherwise create a WireBody of the meanline */
    } else {

        /* collapse the upper and lower surfaces into the meanline
              (over-writing  the old lower surface) */
        ile = (npts - 1) / 2;
        ite =  npts - 1;
        for (i = ile+1, j = ile-1; i < ite; i++, j--) {
            pts[3*i  ] = (pts[3*i  ] + pts[3*j  ]) / 2;
            pts[3*i+1] = (pts[3*i+1] + pts[3*j+1]) / 2;
            pts[3*i+2] = (pts[3*i+2] + pts[3*j+2]) / 2;
        }
    }

    /* set the ouput */
    *ppts = pts;

cleanup:
    delete [] polyTop;
    delete [] polyBot;

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   parsecPoly - constructs the top and bottom parsec polynomials      *
 *                                                                      *
 ************************************************************************
 */

template<class T>
static int
parsecPolyCoeff(int    iudp,            /* (in)  udp index */
                int&   npoly,           /* (out) number of polynomial coefficients */
                T      **polyTop,       /* (out) polynomial coefficients for top */
                T      **polyBot)       /* (out) polynomial coefficients for bottom */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int nparam, ipoly;
    T   yte, rle, lefac, xtop, ytop, ctop, ttop, xbot, ybot, cbot, tbot;

    /* get arguments */
    npoly  = udps[iudp].arg[1].size;
    nparam = udps[iudp].arg[2].size;

    /* create the polynomial coefficients if Nparam > 0 */
    if (nparam > 1) {
        npoly = 6;

        yte  = YTE<T>(iudp);

        rle  = PARAM<T>(iudp,0);

        xtop = PARAM<T>(iudp,1);
        ytop = PARAM<T>(iudp,2);
        ctop = PARAM<T>(iudp,3);
        ttop = PARAM<T>(iudp,4);

        xbot = PARAM<T>(iudp,5);
        ybot = PARAM<T>(iudp,6);
        cbot = PARAM<T>(iudp,7);
        tbot = PARAM<T>(iudp,8);

        /* create top polynomial */
        *polyTop = new T[npoly];
        if (*polyTop == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        lefac = sqrt(2*rle);

        status = polyfit( lefac, yte, xtop, ytop, ctop, ttop, *polyTop );
        if (status != EGADS_SUCCESS) {
            printf(" udpParsec.parsecPolyCoeff: top matrix is singular\n");
            goto cleanup;
        }

        /* create bottom polynomial */
        *polyBot = new T[npoly];
        if (*polyBot == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        lefac = -sqrt(2*rle);

        status = polyfit( lefac, yte, xbot, ybot, cbot, tbot, *polyBot );
        if (status != EGADS_SUCCESS) {
            printf(" udpParsec.parsecPolyCoeff: bot matrix is singular\n");
            goto cleanup;
        }

    /* use given polynomial coefficients because Npoly > 0 */
    } else {
        npoly = 1 + npoly / 2;

        *polyTop = new T[npoly+1];
        if (*polyTop == NULL) {
            status  = EGADS_MALLOC;
            goto cleanup;
        }

        *polyBot = new T[npoly+1];
        if (*polyBot == NULL) {
            status  = EGADS_MALLOC;
            goto cleanup;
        }

        /* copy polynomial coefficients and close airfoil */
        (*polyTop)[npoly-1] = YTE<T>(iudp);
        for (ipoly = 0; ipoly < npoly-1; ipoly++) {
            (*polyTop)[ipoly  ]  = POLY<T>(iudp,ipoly);
            (*polyTop)[npoly-1] -= POLY<T>(iudp,ipoly);
        }

        (*polyBot)[npoly-1] = YTE<T>(iudp);
        for (ipoly = 0; ipoly < npoly-1; ipoly++) {
            (*polyBot)[ipoly  ]  = POLY<T>(iudp,ipoly+npoly-1);
            (*polyBot)[npoly-1] -= POLY<T>(iudp,ipoly+npoly-1);
        }
    }

cleanup:

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   parsec - Evaluates parsec polynomial at an x-coordinate            *
 *                                                                      *
 ************************************************************************
 */

template<class T>
static T
parsec(const double x,                  /* (in)  independent variable */
       const T      poly[],             /* (in)  polynomial coefficients */
       const int    npoly)              /* (in)  number of polynomial coefficients */
{
    T y = 0;                            /* (out) dependent variable */

    for (int ipoly = 0; ipoly < npoly; ipoly++) {
        y += poly[ipoly] * pow(x, ipoly+0.5);
    }

    return y;
}


/*
 ************************************************************************
 *                                                                      *
 *   polyfit - Compute the polynomial fit for parsec                    *
 *                                                                      *
 ************************************************************************
 */

template<class T>
static int
polyfit(const T&   lefac,               /* (in)  leading edge factor */
        const T&   yte,                 /* (in)  height of trailing edge */
        const T&   x,                   /* (in)  x-location at bot */
        const T&   y,                   /* (in)  y-location at bot */
        const T&   c,                   /* (in)  d2x/dy2 at bot */
        const T&   t,                   /* (in)  bot theta at trailing edge (deg) */
        T          poly[])              /* (out) polynomial coeffs fit to input data */
{
    T Amat[25], Rhs[5];

    Amat[ 0] = 1;
    Amat[ 1] = 1;
    Amat[ 2] = 1;
    Amat[ 3] = 1;
    Amat[ 4] = 1;
    Rhs[  0] = yte - lefac;

    Amat[ 5] =      pow(x, 1.5);
    Amat[ 6] =      pow(x, 2.5);
    Amat[ 7] =      pow(x, 3.5);
    Amat[ 8] =      pow(x, 4.5);
    Amat[ 9] =      pow(x, 5.5);
    Rhs[  1] = y - lefac * sqrt(x);

    Amat[10] =  1.5;
    Amat[11] =  2.5;
    Amat[12] =  3.5;
    Amat[13] =  4.5;
    Amat[14] =  5.5;
    Rhs[  2] = tan(t*PI/180) - 0.5 * lefac;

    Amat[15] =  1.5  * pow(x, 0.5);
    Amat[16] =  2.5  * pow(x, 1.5);
    Amat[17] =  3.5  * pow(x, 2.5);
    Amat[18] =  4.5  * pow(x, 3.5);
    Amat[19] =  5.5  * pow(x, 4.5);
    Rhs[  3] =  -.5  * lefac / sqrt(x);

    Amat[20] =  0.75 * pow(x, -.5);
    Amat[21] =  3.75 * pow(x, 0.5);
    Amat[22] =  8.75 * pow(x, 1.5);
    Amat[23] = 15.75 * pow(x, 2.5);
    Amat[24] = 24.75 * pow(x, 3.5);
    Rhs[  4] = c + 0.25 * lefac / pow(x, 1.5);

    poly[0] = lefac;

    return matsol(Amat, Rhs, 5, &poly[1]);
}


/*
 ************************************************************************
 *                                                                      *
 *   matsol - Gaussian elimination with partial pivoting                *
 *                                                                      *
 ************************************************************************
 */

template<class T>
int
matsol(T         A[],                   /* (in)  matrix to be solved (stored rowwise) */
                                        /* (out) upper-triangular form of matrix */
       T         b[],                   /* (in)  right hand side */
                                        /* (out) right-hand side after swapping */
       int       n,                     /* (in)  size of matrix */
       T         x[])                   /* (out) solution of A*x=b */
{
    int       status = EGADS_SUCCESS;   /* (out) return status */

    int   ir, jc, kc, imax;
    T     amax, swap, fact;

    /* --------------------------------------------------------------- */

    /* reduce each column of A */
    for (kc = 0; kc < n; kc++) {

        /* find pivot element */
        imax = kc;
        amax = fabs(A[kc*n+kc]);

        for (ir = kc+1; ir < n; ir++) {
            if (fabs(A[ir*n+kc]) > amax) {
                imax = ir;
                amax = fabs(A[ir*n+kc]);
            }
        }

        /* check for possibly-singular matrix (ie, near-zero pivot) */
        if (amax < EPS12) {
            status  = EGADS_DEGEN;
            goto cleanup;
        }

        /* if diagonal is not pivot, swap rows in A and b */
        if (imax != kc) {
            for (jc = 0; jc < n; jc++) {
                swap         = A[kc  *n+jc];
                A[kc  *n+jc] = A[imax*n+jc];
                A[imax*n+jc] = swap;
            }

            swap    = b[kc  ];
            b[kc  ] = b[imax];
            b[imax] = swap;
        }

        /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
        for (ir = kc+1; ir < n; ir++) {
            fact = A[ir*n+kc] / A[kc*n+kc];

            for (jc = kc+1; jc < n; jc++) {
                A[ir*n+jc] -= fact * A[kc*n+jc];
            }

            b[ir] -= fact * b[kc];

            A[ir*n+kc] = 0;
        }
    }

    /* back-substitution pass */
    x[n-1] = b[n-1] / A[(n-1)*n+(n-1)];

    for (jc = n-2; jc >= 0; jc--) {
        x[jc] = b[jc];
        for (kc = jc+1; kc < n; kc++) {
            x[jc] -= A[jc*n+kc] * x[kc];
        }
        x[jc] /= A[jc*n+jc];
    }

cleanup:
    return status;
}
