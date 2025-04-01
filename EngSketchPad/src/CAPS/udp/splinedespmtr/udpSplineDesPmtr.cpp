/*
 *************************************************************************
 *                                                                       *
 * udpSplineDesPmtr -- udf for b-spline design parameter                 *
 *                                                                       *
 *            Written by Marshall Galbraith @ MIT                        *
 *            Patterned after .c code written by John Dannenhoffer       *
 *                                                @ Syracuse University  *
 *                                                                       *
 *************************************************************************
 */

/*
 * Copyright (C) 2011/2025 Marshall Galbraith @ MIT
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

extern "C" {
#define NUMUDPARGS 5
#include "udpUtilities.h"
}

#define PARAM_ARG(UDP   )                    (UDP)->arg[0]
#define PARAM_VAL(UDP, I) ((const double *) ((UDP)->arg[0].val))[I]
#define PARAM_DOT(UDP, I) ((const double *) ((UDP)->arg[0].dot))[I]
#define VALUE_ARG(UDP   )                    (UDP)->arg[1]
#define VALUE_VAL(UDP, I) ((const double *) ((UDP)->arg[1].val))[I]
#define VALUE_DOT(UDP, I) ((const double *) ((UDP)->arg[1].dot))[I]
#define CONT_ARG( UDP   )                    (UDP)->arg[2]
#define CONT_VAL( UDP, I) ((const int *)    ((UDP)->arg[2].val))[I]
#define TS_ARG(   UDP   )                    (UDP)->arg[3]
#define TS_VAL(   UDP, I) ((const double *) ((UDP)->arg[3].val))[I]
#define TS_DOT(   UDP, I) ((const double *) ((UDP)->arg[3].dot))[I]
#define INTRP_ARG(UDP   )                    (UDP)->arg[4]
#define INTRP_VAL(UDP, I) ((double *)       ((UDP)->arg[4].val))[I]
#define INTRP_DOT(UDP, I) ((double *)       ((UDP)->arg[4].dot))[I]

#define PARAM_SURREAL(UDP, I) SurrealS<1>( PARAM_VAL(UDP, I), PARAM_DOT(UDP, I))
#define VALUE_SURREAL(UDP, I) SurrealS<1>( VALUE_VAL(UDP, I), VALUE_DOT(UDP, I))
#define TS_SURREAL(   UDP, I) SurrealS<1>( TS_VAL(   UDP, I), TS_DOT(   UDP, I))
#define INTRP_SURREAL(UDP, I) SurrealS<1>( INTRP_VAL(UDP, I), INTRP_DOT(UDP, I))

/* data about possible arguments */
static const char  *argNames[NUMUDPARGS] = {"param",     "value",     "cont",  "ts",        "interp"    };
static int          argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRINT, ATTRREALSEN, -ATTRREALSEN};
static int          argIdefs[NUMUDPARGS] = {0 ,          0 ,          2 ,      0,           0,          };
static double       argDdefs[NUMUDPARGS] = {0.,          0.,          0.,      0.,          0.,         };


#include "egadsSplineFit.h" /* EG_spline1dFit2 */
#include "egads_dot.h"

template<class T> T    PARAM(udp_T *udp, int i);
template<> double      PARAM< double      >(udp_T *udp, int i ) { return PARAM_VAL(udp,i); }
template<> SurrealS<1> PARAM< SurrealS<1> >(udp_T *udp, int i ) { return PARAM_SURREAL(udp,i); }

template<class T> T    VALUE(udp_T *udp, int i);
template<> double      VALUE< double      >(udp_T *udp, int i ) { return VALUE_VAL(udp,i); }
template<> SurrealS<1> VALUE< SurrealS<1> >(udp_T *udp, int i ) { return VALUE_SURREAL(udp,i); }

template<class T> T    TS(udp_T *udp, int i);
template<> double      TS< double      >( udp_T *udp, int i       ) { return TS_VAL(udp, i); }
template<> SurrealS<1> TS< SurrealS<1> >( udp_T *udp, int i       ) { return TS_SURREAL(udp, i); }

template<class T> T    INTRP(udp_T *udp, int i );
template<> double      INTRP< double      >( udp_T *udp, int i ) { return INTRP_VAL(udp, i); }
template<> SurrealS<1> INTRP< SurrealS<1> >( udp_T *udp, int i ) { return INTRP_SURREAL(udp, i); }

template class SurrealS<1>;
typedef SurrealS<1> SurrealS1;

#include <vector>

#define LENMSG 512

extern "C" {
/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
}

// Spline fitting tolerance
#define TOL 1e-6
#define NEWTONMAX 21

/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */
extern "C"
int
udpExecute(ego context,                 /* (in)  context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
  int    status = EGADS_SUCCESS;
  ego    ecurve, spline_body = NULL;
  double data[18], tdata[2], t, x, *rdata=NULL;
  int    header[4], sense[1], end1c, endnc;
  ego    eedge = NULL, enodes[2], eloop;
  int    i, j, k, nkn;
//  udpDotCache_T *cache;
  char          *message=NULL;
  void    *realloc_temp=NULL;
  udp_T   *udps = *Udps;
  udp_T   *udp = &udps[0];

  std::vector<int> cn;
  std::vector<double> kn;
  std::vector<double> xyz;

  ROUTINE(udpExecute);

  /* --------------------------------------------------------------- */

#ifdef DEBUG
  printf("udpSplineDesPmtr.udpExecute\n");
#endif

  /* default return values */
  *ebody  = NULL;
  *nMesh  = 0;
  *string = NULL;

  MALLOC(message, char, LENMSG);
  message[0] = '\0';

  if (PARAM_ARG(udp).ncol != 1) {
    snprintf(message, LENMSG, "Expecting 1 'param' columns. There are %d columns \n", PARAM_ARG(udp).ncol);
    status = EGADS_RANGERR;
    goto cleanup;
  }
  if (PARAM_ARG(udp).nrow < 2) {
    snprintf(message, LENMSG, "Expecting at least 2 'param' rows. There are %d rows \n", PARAM_ARG(udp).nrow);
    status = EGADS_RANGERR;
    goto cleanup;
  }
  if (VALUE_ARG(udp).ncol != 1) {
    snprintf(message, LENMSG, "Expecting 1 'value' columns. There are %d columns \n", VALUE_ARG(udp).ncol);
    status = EGADS_RANGERR;
    goto cleanup;
  }
  if (VALUE_ARG(udp).nrow < 2) {
    snprintf(message, LENMSG, "Expecting at least 2 'value' rows. There are %d rows \n", VALUE_ARG(udp).nrow);
    status = EGADS_RANGERR;
    goto cleanup;
  }
  if (PARAM_ARG(udp).nrow != VALUE_ARG(udp).nrow) {
    snprintf(message, LENMSG, "Expecting 'param' rows (%d) and 'value' rows (%d) to be equal\n", PARAM_ARG(udp).nrow, VALUE_ARG(udp).nrow);
    status = EGADS_RANGERR;
    goto cleanup;
  }
  if (CONT_ARG(udp).nrow > 1) {
    if (CONT_ARG(udp).ncol != 1) {
      snprintf(message, LENMSG, "Expecting 1 'cont' columns. There are %d columns \n", CONT_ARG(udp).ncol);
      status = EGADS_RANGERR;
      goto cleanup;
    }
    if (PARAM_ARG(udp).nrow != CONT_ARG(udp).nrow) {
      snprintf(message, LENMSG, "Expecting 'param' rows (%d) and 'cont' rows (%d) to be equal\n", PARAM_ARG(udp).nrow, CONT_ARG(udp).nrow);
      status = EGADS_RANGERR;
      goto cleanup;
    }
    for (i = 0; i < CONT_ARG(udp).nrow; i++) {
      if (CONT_VAL(udp, i) < 0 || CONT_VAL(udp, i) > 2) {
        snprintf(message, LENMSG, "Expecting 0 <= 'cont[%d]' = %d <= 2\n", i+1, CONT_VAL(udp, i));
        status = EGADS_RANGERR;
        goto cleanup;
      }
    }
  }
  if (TS_ARG(udp).size < 1) {
    snprintf(message, LENMSG, "Expecting at least 1 'ts'. There are %d \n", TS_ARG(udp).size);
    status = EGADS_RANGERR;
    goto cleanup;
  }

  /* Cache copy of arguments for future use.
   * This also increments numUdp and copies udps[numUdp] to udps[numUdp].
   * Caching should only be performed after checking for valid inputs.
   */
  status = cacheUdp(NULL);
  CHECK_STATUS(cacheUdp);
  udp = &udps[numUdp];

  /* make enough room for the interpolated output */
  INTRP_ARG(udp).size = TS_ARG(udp).nrow;
  INTRP_ARG(udp).nrow = TS_ARG(udp).nrow;
  INTRP_ARG(udp).ncol = 1;

  RALLOC(INTRP_ARG(udp).val, double, INTRP_ARG(udp).size);
  RALLOC(INTRP_ARG(udp).dot, double, INTRP_ARG(udp).size);

  cn.resize(PARAM_ARG(udp).nrow, 0);
  if (CONT_ARG(udp).nrow > 1) {
    for (i = 0; i < CONT_ARG(udp).nrow; i++) {
      cn[i] = 2-CONT_VAL(udp,i);
    }
  }

  end1c = cn.front();
  endnc = cn.back();

  /* use finite difference for linear */
  if (PARAM_ARG(udp).nrow == 2) end1c = endnc = 1;

  cn.front() = 0;
  cn.back() = 0;

  nkn = PARAM_ARG(udp).nrow;
  for (i = 0; i < (int)cn.size(); i++) {
    nkn += cn[i];
  }

  /* get the knot values */
  kn.resize(nkn);
  for (k = i = 0; i < PARAM_ARG(udp).nrow; i++) {
    for (j = 0; j <= cn[i]; j++) {
      kn[k++] = PARAM_VAL(udp, i);
    }
  }

  xyz.resize(3*nkn);
  for (k = i = 0; i < PARAM_ARG(udp).nrow; i++) {
    for (j = 0; j <= cn[i]; j++) {
      xyz[3*k+0] = PARAM_VAL(udp, i);
      xyz[3*k+1] = VALUE_VAL(udp, i);
      xyz[3*k+2] = 0;
      k++;
    }
  }

  status = EG_spline1dFit2c( end1c, endnc, kn.size(), xyz.data(), kn.data(), TOL, header, &rdata );
  CHECK_STATUS(EG_spline1dFit2);

  status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rdata, &ecurve);
  CHECK_STATUS(EG_makeGeometry);
  FREE(rdata);

  /* create Topology */

  /* first Node */
  i = 0;
  data[0] = PARAM_VAL(udp, i);
  data[1] = VALUE_VAL(udp, i);
  data[2] = 0;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &(enodes[0]));
  CHECK_STATUS(EG_makeTopology);

  /* second Node */
  i = PARAM_ARG(udp).nrow-1;
  data[0] = PARAM_VAL(udp, i);
  data[1] = VALUE_VAL(udp, i);
  data[2] = 0;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &(enodes[1]));
  CHECK_STATUS(EG_makeTopology);

  /* Edge */
  tdata[0] = kn.front();
  tdata[1] = kn.back();
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &(eedge));
  CHECK_STATUS(EG_makeTopology);

  /* create Loop of the Edge */
  sense[0] = SFORWARD;
  status = EG_makeTopology(context, NULL, LOOP, OPEN,
                           NULL, 1, &eedge, sense, &eloop);
  CHECK_STATUS(EG_makeTopology);

  /* create the WireBody */
  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, &spline_body);
  CHECK_STATUS(EG_makeTopology);

  /* Evaluate the spline */
  for (i = 0; i < TS_ARG(udp).size; i++) {
    t = x = TS_VAL(udp, i);
    for (j = 0; j < NEWTONMAX; j++) {
      status = EG_evaluate(ecurve, &t, data);
      CHECK_STATUS(EG_evaluate);

      if (fabs(x - data[0]) < TOL) break;

      // Newton update the t-value until x and data[0] match
      double R   = x - data[0];
      double R_t =   - data[3];
      double dt  = R/R_t;
      t -= dt;
    }
    INTRP_VAL(udp, i) = data[1];
    //std::cout << "interp val: " << i << " : " << t << " " << data[0] << " : " << data[1] << std::endl;
  }

  /* cleanup in reverse order */
  status = EG_deleteObject(eloop);
  CHECK_STATUS(EG_deleteObject);
  status = EG_deleteObject(eedge);
  CHECK_STATUS(EG_deleteObject);
  status = EG_deleteObject(ecurve);
  CHECK_STATUS(EG_deleteObject);
  status = EG_deleteObject(enodes[1]);
  CHECK_STATUS(EG_deleteObject);
  status = EG_deleteObject(enodes[0]);
  CHECK_STATUS(EG_deleteObject);

  /* remember this model (body) */
  udp->ebody = *ebody = spline_body;

#ifdef DEBUG
  printf("udpExecute -> *ebody=%lx\n", (long)(*ebody));
#endif

cleanup:
  if (strlen(message) > 0) {
    *string = message;
    printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
    printf("ERROR: status = %d \n", status);
    // freePrivateData(&ffdBox);
    FREE(message);
    *string = udpErrorStr(status);
  } else {
    FREE(message);
  }

  return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                      *
 ************************************************************************
 */
extern "C"
int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
  int    status = EGADS_SUCCESS;
  int    iudp = 0, judp, stride, ipnt, header[4], *senses;
  int    nchild, nedge, nnode;
  int    oclass, mtype, end1c, endnc;
  double point[18], point_dot[18], data[4], x;
  ego    eent, eref, eloop, *eedges, ecurve, *enodes, *echildren;
  SurrealS<1> t, pnt[3], tdata[2], sdata[18], *rdata=NULL;
  udp_T   *udp = NULL;

  int    i, j, k;
  char   *message=NULL;

  ROUTINE(udpSensitivity);

  /* --------------------------------------------------------------- */

  /* check that ebody matches one of the ebodys */
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

  udp = &udps[iudp];

#ifdef DEBUG
  printf("udpSplineDesPmtr.udpSensitivity iudp = %d\n", iudp);

  for (i = 0; i < PARAM_ARG(udp).nrow; i++) {
    printf("PARAM_VAL(%d) = %f %f\n", i, PARAM_VAL(udp,i) );
  }
  for (i = 0; i < VALUE_ARG(udp).nrow; i++) {
    printf("VALUE_VAL(%d) = %f %f\n", i, VALUE_VAL(udp,i) );
  }
  for (i = 0; i < CONT_ARG(udp).nrow; i++) {
    printf("CONT_VAL(%d) = %f %f\n", i, CONT_VAL(udp,i) );
  }
  for (i = 0; i < TS_ARG(udp).nrow; i++) {
      printf("TS_VAL(%d) = %f\n", i, TS_VAL(udp,i) );
  }
  for (i = 0; i < INTRP_ARG(udp).nrow; i++) {
      printf("INTRP_VAL(%d) = %f\n", i, INTRP_VAL(udp,i) );
  }
  printf("\n");
#endif

  MALLOC(message, char, LENMSG);
  message[0] = '\0';

  /* get the loop from the body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype, data, &nchild, &echildren,
                          &senses);
  CHECK_STATUS(EG_getTopology);
  eloop = echildren[0];

  /* get the edges from the loop */
  status = EG_getTopology(eloop, &eref, &oclass, &mtype, data, &nedge, &eedges,
                          &senses);
  CHECK_STATUS(EG_getTopology);

  /* get the nodes and the curve from the first edge */
  status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, data, &nnode, &enodes,
                          &senses);
  CHECK_STATUS(EG_getTopology);


  /* build the sensitivity if needed */
  if (udp->ndotchg > 0 || EG_hasGeometry_dot(ebody) != EGADS_SUCCESS) {

    std::vector<int> cn(PARAM_ARG(udp).nrow, 0);
    if (CONT_ARG(udp).nrow > 1) {
      for (i = 0; i < CONT_ARG(udp).nrow; i++) {
        cn[i] = 2-CONT_VAL(udp,i);
      }
    }

    end1c = cn.front();
    endnc = cn.back();

    /* use finite difference for linear */
    if (PARAM_ARG(udp).nrow == 2) end1c = endnc = 1;

    cn.front() = 0;
    cn.back() = 0;

    int nkn = PARAM_ARG(udp).nrow;
    for (i = 0; i < (int)cn.size(); i++) {
      nkn += cn[i];
    }

    std::vector<SurrealS1> kn(nkn);
    for (k = i = 0; i < PARAM_ARG(udp).nrow; i++) {
      for (j = 0; j <= cn[i]; j++) {
        kn[k++] = PARAM_SURREAL(udp, i);
      }
    }

    std::vector<SurrealS1> xyz(3*nkn);
    for (k = i = 0; i < PARAM_ARG(udp).nrow; i++) {
      for (j = 0; j <= cn[i]; j++) {
        xyz[3*k+0] = PARAM_SURREAL(udp, i);
        xyz[3*k+1] = VALUE_SURREAL(udp, i);
        xyz[3*k+2] = 0;
        k++;
      }
    }

    status = EG_spline1dFit2c( end1c, endnc, kn.size(), xyz.data(), kn.data(), TOL, header, &rdata );
    CHECK_STATUS(EG_spline1dFit2);

    /* set the sensitivity of the Curve */
    status = EG_setGeometry_dot(ecurve, CURVE, BSPLINE, header, rdata);
    CHECK_STATUS(EG_setGeometry_dot);
    FREE(rdata);

    /* set the sensitivity of the Node at trailing edge */
    i = 0;
    pnt[0] = PARAM_SURREAL(udp, i);
    pnt[1] = VALUE_SURREAL(udp, i);
    pnt[2] = 0;
    status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, pnt);
    CHECK_STATUS(EG_setGeometry_dot);

    /* set the sensitivity of the Node at leading edge */
    i = PARAM_ARG(udp).nrow-1;
    pnt[0] = PARAM_SURREAL(udp, i);
    pnt[1] = VALUE_SURREAL(udp, i);
    pnt[2] = 0;
    status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, pnt);
    CHECK_STATUS(EG_setGeometry_dot);

    /* set Edge t-range sensitivity */
    tdata[0] = kn.front();
    tdata[1] = kn.back();
    status = EG_setRange_dot(eedges[0], EDGE, tdata);
    CHECK_STATUS(EG_setRange_dot);

    /* update sensitivities for the interpolated data points*/
    for (i = 0; i < TS_ARG(udp).size; i++) {
      t = TS_SURREAL(udp, i);
      x = t.value();
      for (j = 0; j < NEWTONMAX; j++) {
        status = EG_evaluate(ecurve, &t, sdata);
        CHECK_STATUS(EG_evaluate);

        if (fabs(x - sdata[0]) < TOL) break;

        // Newton update the t-value until x and data[0] match
        SurrealS<1> R   = x - sdata[0];
        SurrealS<1> R_t =   - sdata[3];
        SurrealS<1> dt  = R/R_t;
        t -= dt;
      }
      INTRP_SURREAL(udp, i) = sdata[1];
      //std::cout << "interp sens: " << i << " : " << t << " " << sdata[0] << " : " << sdata[1] << std::endl;
    }
  }


  /* find the ego entity */
  if (entType == OCSM_NODE) {
    status = EG_getBodyTopos(ebody, NULL, NODE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 0;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else if (entType == OCSM_EDGE) {
    status = EG_getBodyTopos(ebody, NULL, EDGE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 1;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else if (entType == OCSM_FACE) {
    status = EG_getBodyTopos(ebody, NULL, FACE, &nchild, &echildren);
    CHECK_STATUS(EG_getBodyTopos);

    stride = 2;
    eent = echildren[entIndex-1];

    FREE(echildren);
  } else {
    printf("udpSensitivity: bad entType=%d\n", entType);
    status = EGADS_GEOMERR;
    goto cleanup;
  }


  /* get the velocities from the entity */
  for (ipnt = 0; ipnt < npnt; ipnt++) {
    status = EG_evaluate_dot(eent, &(uvs[stride*ipnt]), NULL, point, point_dot);
    CHECK_STATUS(EG_evaluate_dot);

    /* return the point velocity */
    vels[3*ipnt  ] = point_dot[0];
    vels[3*ipnt+1] = point_dot[1];
    vels[3*ipnt+2] = point_dot[2];
  }


cleanup:
  if (strlen(message) > 0) {
   printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
   printf("ERROR: status = %d \n", status);
  }

  FREE( rdata );
  FREE(message);

  return status;

}
