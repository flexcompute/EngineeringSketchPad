/*
 ************************************************************************
 *                                                                      *
 * udpSupell -- udp file to generate a super-ellipse                    *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 19
#include "udpUtilities.h"

#undef  OFFSET        /* needed to remove duplicate macro definitions */

/* shorthands for accessing argument values and velocities */
#define RX(      IUDP)  ((double *) (udps[IUDP].arg[ 0].val))[0]
#define RX_DOT(  IUDP)  ((double *) (udps[IUDP].arg[ 0].dot))[0]
#define RX_W(    IUDP)  ((double *) (udps[IUDP].arg[ 1].val))[0]
#define RX_W_DOT(IUDP)  ((double *) (udps[IUDP].arg[ 1].dot))[0]
#define RX_E(    IUDP)  ((double *) (udps[IUDP].arg[ 2].val))[0]
#define RX_E_DOT(IUDP)  ((double *) (udps[IUDP].arg[ 2].dot))[0]
#define RY(      IUDP)  ((double *) (udps[IUDP].arg[ 3].val))[0]
#define RY_DOT(  IUDP)  ((double *) (udps[IUDP].arg[ 3].dot))[0]
#define RY_S(    IUDP)  ((double *) (udps[IUDP].arg[ 4].val))[0]
#define RY_S_DOT(IUDP)  ((double *) (udps[IUDP].arg[ 4].dot))[0]
#define RY_N(    IUDP)  ((double *) (udps[IUDP].arg[ 5].val))[0]
#define RY_N_DOT(IUDP)  ((double *) (udps[IUDP].arg[ 5].dot))[0]
#define N(       IUDP)  ((double *) (udps[IUDP].arg[ 6].val))[0]
#define N_W(     IUDP)  ((double *) (udps[IUDP].arg[ 7].val))[0]
#define N_E(     IUDP)  ((double *) (udps[IUDP].arg[ 8].val))[0]
#define N_S(     IUDP)  ((double *) (udps[IUDP].arg[ 9].val))[0]
#define N_N(     IUDP)  ((double *) (udps[IUDP].arg[10].val))[0]
#define N_SW(    IUDP)  ((double *) (udps[IUDP].arg[11].val))[0]
#define N_SE(    IUDP)  ((double *) (udps[IUDP].arg[12].val))[0]
#define N_NW(    IUDP)  ((double *) (udps[IUDP].arg[13].val))[0]
#define N_NE(    IUDP)  ((double *) (udps[IUDP].arg[14].val))[0]
#define OFFSET(  IUDP)  ((double *) (udps[IUDP].arg[15].val))[0]
#define NQUAD(   IUDP)  ((int    *) (udps[IUDP].arg[16].val))[0]
#define NUMPNTS( IUDP)  ((int    *) (udps[IUDP].arg[17].val))[0]
#define SLPFACT( IUDP)  ((double *) (udps[IUDP].arg[18].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"rx",        "rx_w",      "rx_e",
                                      "ry",        "ry_s",      "ry_n",
                                      "n",         "n_w",       "n_e",       "n_s",    "n_n",
                                      "n_sw",      "n_se",      "n_nw",      "n_ne",
                                      "offset",    "nquad",     "numpnts",   "slpfact",     };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREALSEN,
                                      ATTRREALSEN, ATTRREALSEN, ATTRREALSEN,
                                      ATTRREAL,    ATTRREAL,    ATTRREAL,    ATTRREAL, ATTRREAL,
                                      ATTRREAL,    ATTRREAL,    ATTRREAL,    ATTRREAL,
                                      ATTRREAL,    ATTRINT,     ATTRINT,     ATTRREAL,         };
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,
                                      0,           0,           0,
                                      0,           0,           0,           0,        0,
                                      0,           0,           0,           0,
                                      0,           4,           11,          0,                };
static double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,
                                      0.,          0.,          0.,
                                      2.,          0.,          0.,          0.,       0.,
                                      0.,          0.,          0.,          0.,
                                      0.,          0.,          0.,          0.,               };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

#ifdef GRAFIC
   #include "grafic.h"
#endif


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     npnt, sense[4], sizes[2], periodic, i, add=1;
    double  rx_w, rx_e, ry_s, ry_n, n_sw, n_se, n_nw, n_ne, theta;
    double  *pnt_ne=NULL, *pnt_nw=NULL, *pnt_sw=NULL, *pnt_se=NULL, *pnt_save=NULL;
    double  data[18], tdata[2], result[3], range[4], eval[18], norm[3], dx, dy, ds;
    double  dxytol = 1.0e-6;
    char    *message=NULL;
    ego     enodes[5], eedges[4], ecurve[4], eloop, eface, enew;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("rx(     0) = %f\n", RX(     0));
    printf("rx_w(   0) = %f\n", RX_W(   0));
    printf("rx_e(   0) = %f\n", RX_E(   0));
    printf("ry(     0) = %f\n", RY(     0));
    printf("ry_s(   0) = %f\n", RY_S(   0));
    printf("ry_n(   0) = %f\n", RY_N(   0));
    printf("n(      0) = %f\n", N(      0));
    printf("n_w(    0) = %f\n", N_W(    0));
    printf("n_e(    0) = %f\n", N_E(    0));
    printf("n_s(    0) = %f\n", N_S(    0));
    printf("n_n(    0) = %f\n", N_N(    0));
    printf("n_sw(   0) = %f\n", N_SW(   0));
    printf("n_se(   0) = %f\n", N_SE(   0));
    printf("n_nw(   0) = %f\n", N_NW(   0));
    printf("n_ne(   0) = %f\n", N_NE(   0));
    printf("offset( 0) = %f\n", OFFSET( 0));
    printf("nquad(  0) = %d\n", NQUAD(  0));
    printf("numpnts(0) = %d\n", NUMPNTS(0));
    printf("slpfact(0) = %f\n", SLPFACT(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    for (i = 0; i < NUMUDPARGS; i++) {
        if (udps[0].arg[i].size > 1) {
            snprintf(message, 100, "SUPELL: all arguments should be a scalar");
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }

    if (NQUAD(0) != 1 && NQUAD(0) != 2 && NQUAD(0) != 4) {
        snprintf(message, 100, "SUPELL: nquad (%d) should be 1, 2, or 4", NQUAD(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* set up the parameters in each of the quadrants */
    rx_w = rx_e = 0;
    ry_s = ry_n = 0;
    n_sw = n_se = n_nw = n_ne = 2;

    if (RX(  0) > 0)  rx_w = rx_e = RX(  0);
    if (RX_W(0) > 0)  rx_w        = RX_W(0);
    if (RX_E(0) > 0)         rx_e = RX_E(0);

    if (RY(  0) > 0)  ry_s = ry_n = RY(  0);
    if (RY_S(0) > 0)  ry_s        = RY_S(0);
    if (RY_N(0) > 0)         ry_n = RY_N(0);

    if (N(   0) > 0)  n_sw = n_se = n_nw = n_ne = N(   0);
    if (N_W( 0) > 0)  n_sw         = n_nw       = N_W( 0);
    if (N_E( 0) > 0)         n_se        = n_ne = N_E( 0);
    if (N_S( 0) > 0)  n_sw = n_se               = N_S( 0);
    if (N_N( 0) > 0)                n_nw = n_ne = N_N( 0);
    if (N_SW(0) > 0)  n_sw                      = N_SW(0);
    if (N_SE(0) > 0)         n_se               = N_SE(0);
    if (N_NW(0) > 0)                n_nw        = N_NW(0);
    if (N_NE(0) > 0)                       n_ne = N_NE(0);

    if        (rx_w <= 0) {
        snprintf(message, 100, "SUPELL: rx_w should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (rx_e <= 0) {
        snprintf(message, 100, "SUPELL: rx_e should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (ry_s <= 0) {
        snprintf(message, 100, "SUPELL: ry_s should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (ry_n <= 0) {
        snprintf(message, 100, "SUPELL: ry_n should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (n_sw <= 0) {
        snprintf(message, 100, "SUPELL: n_sw should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (n_se <= 0) {
        snprintf(message, 100, "SUPELL: n_se should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (n_nw <= 0) {
        snprintf(message, 100, "SUPELL: n_nw should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (n_ne <= 0) {
        snprintf(message, 100, "SUPELL: n_ne should be positive");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (NUMPNTS(0) < 11) {
        snprintf(message, 100, "SUPELL: numpnts must be at least 11");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (SLPFACT(0) > 0 && OFFSET(0) != 0) {
        snprintf(message, 100, "SUPELL: both offset and slpfact cannot be set");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (SLPFACT(0) > 0 && (n_sw < 2 || n_se < 2 || n_nw < 2 || n_ne < 2)) {
        snprintf(message, 100, "SUPELL: slpfact cannot be set if n < 2");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("rx_w    = %f\n", rx_w);
    printf("rx_e    = %f\n", rx_e);
    printf("ry_s    = %f\n", ry_s);
    printf("ry_n    = %f\n", ry_n);
    printf("n_sw    = %f\n", n_sw);
    printf("n_se    = %f\n", n_se);
    printf("n_nw    = %f\n", n_nw);
    printf("n_ne    = %f\n", n_ne);
    printf("nquad   = %d\n", NQUAD(  numUdp));
    printf("numpnts = %d\n", NUMPNTS(numUdp));
    printf("slpfact = %f\n", SLPFACT(numUdp));
#endif

    /* mallocs required by Windows compiler */
    MALLOC(pnt_ne,   double, 9*NUMPNTS(0)+6);
    MALLOC(pnt_nw,   double, 9*NUMPNTS(0)+6);
    MALLOC(pnt_sw,   double, 9*NUMPNTS(0)+6);
    MALLOC(pnt_se,   double, 9*NUMPNTS(0)+6);
    MALLOC(pnt_save, double, 9*NUMPNTS(0)+6);

    /* create Nodes at the four cardinal directions */
    data[0] = +rx_e + OFFSET(0);
    data[1] =  0;
    data[2] =  0;
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    data[0] =  0;
    data[1] = +ry_n + OFFSET(0);
    data[2] =  0;
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    data[0] = -rx_w - OFFSET(0);
    data[1] =  0;
    data[2] =  0;
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    data[0] =  0;
    data[1] = -ry_s - OFFSET(0);
    data[2] =  0;
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[3]));
    CHECK_STATUS(EG_makeTopology);

    enodes[4] = enodes[0];

#ifdef DEBUG
    printf("enodes[0]=%llx\n", (long long)(enodes[0]));
    printf("enodes[1]=%llx\n", (long long)(enodes[1]));
    printf("enodes[2]=%llx\n", (long long)(enodes[2]));
    printf("enodes[3]=%llx\n", (long long)(enodes[3]));
#endif

    /* ne quadrant -- generate coordinates */
    npnt = 0;

    if (SLPFACT(0) < 0 && n_ne >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(NUMPNTS(0)-i) / (double)(NUMPNTS(0)) * PI/2;

            pnt_ne[3*npnt  ] = +rx_e * pow(cos(theta), 2.0/n_ne);
            pnt_ne[3*npnt+1] = -ry_n * pow(sin(theta), 2.0/n_ne);
            pnt_ne[3*npnt+2] = 0;
            npnt++;
        }
    }
    
    pnt_ne[3*npnt  ] = +rx_e;
    pnt_ne[3*npnt+1] =  0;
    pnt_ne[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) > 0) {
        pnt_ne[3*npnt  ] = +rx_e;
        pnt_ne[3*npnt+1] = +ry_n * SLPFACT(0);
        pnt_ne[3*npnt+2] = 0;
        npnt++;
    }

    for (i = 1; i < NUMPNTS(0)-1; i++) {
        theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

        pnt_ne[3*npnt  ] = +rx_e * pow(cos(theta), 2.0/n_ne);
        pnt_ne[3*npnt+1] = +ry_n * pow(sin(theta), 2.0/n_ne);
        pnt_ne[3*npnt+2] =  0;
        npnt++;
    }

    if (SLPFACT(0) > 0) {
        pnt_ne[3*npnt  ] = +rx_e * SLPFACT(0);
        pnt_ne[3*npnt+1] = +ry_n;
        pnt_ne[3*npnt+2] = 0;
        npnt++;
    }

    pnt_ne[3*npnt  ] =  0;
    pnt_ne[3*npnt+1] = +ry_n;
    pnt_ne[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) < 0 && n_ne >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(NUMPNTS(0)-1-i) / (double)(NUMPNTS(0)) * PI/2;

            pnt_ne[3*npnt  ] = -rx_e * pow(cos(theta), 2.0/n_ne);
            pnt_ne[3*npnt+1] = +ry_n * pow(sin(theta), 2.0/n_ne);
            pnt_ne[3*npnt+2] = 0;
            npnt++;
        }
    }

    /* create offset if required */
    if (OFFSET(0) != 0) {
        for (i = 0; i < npnt; i++) {
            pnt_save[3*i  ] = pnt_ne[3*i  ];
            pnt_save[3*i+1] = pnt_ne[3*i+1];
        }

        pnt_ne[       0] += OFFSET(0);
        pnt_ne[3*npnt-2] += OFFSET(0);

        for (i = 1; i < npnt-1; i++) {
            dx = pnt_save[3*i+3] - pnt_save[3*i-3];
            dy = pnt_save[3*i+4] - pnt_save[3*i-2];
            ds = sqrt(dx * dx + dy * dy);

            pnt_ne[3*i  ] += OFFSET(0) * dy / ds;
            pnt_ne[3*i+1] -= OFFSET(0) * dx / ds;
        }

#ifdef DEBUG
        for (i = 0; i < npnt; i++) {
            printf("ne %3d:  old: %10.4f %10.4f  new: %10.4f %10.4f\n",
                   i, pnt_save[3*i], pnt_save[3*i+1], pnt_ne[3*i], pnt_ne[3*i+1]);
        }
#endif
    }

    /* create spline curve and Edge (use OpenCASCADE's approximate so
       that subsequent HOLLOWs work) */
    sizes[0] = npnt;
    sizes[1] = 0;
    status = EG_approximate(context, 1, dxytol, sizes, pnt_ne, &(ecurve[0]));
    CHECK_STATUS(EG_approximate);

    data[0] = +rx_e + OFFSET(0);
    data[1] = 0;
    data[2] = 0;
    status = EG_invEvaluate(ecurve[0], data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = 0;
    data[1] = +ry_n + OFFSET(0);
    data[2] = 0;
    status = EG_invEvaluate(ecurve[0], data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

#ifdef GRAFIC
    {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    ilin[3], isym[3], nper[3], nline=0, nplot=0, ipnt;
        int    oclass, mtype, *header;
        float  xplot[2000], yplot[2000];
        double tt, *gdata;
        ego    eref;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt_ne[3*ipnt  ];
            yplot[nplot] = pnt_ne[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = npnt;
        nline++;

        for (ipnt = 0; ipnt < 1001; ipnt++) {
            tt = tdata[0] + (tdata[1]-tdata[0]) * (double)(ipnt) / 1000.;
            status = EG_evaluate(ecurve[0], &tt, data);
            CHECK_STATUS(EG_evaluate);

            xplot[nplot] = data[0];
            yplot[nplot] = data[1];
            nplot++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] = 1001;
        nline++;

        status = EG_getGeometry(ecurve[0], &oclass, &mtype, &eref, &header, &gdata);
        CHECK_STATUS(EG_getGeometry);

        for (ipnt = 0; ipnt < header[2]; ipnt++) {
            xplot[nplot] = gdata[header[3]+3*ipnt  ];
            yplot[nplot] = gdata[header[3]+3*ipnt+1];
            nplot++;
        }
        ilin[nline] = +GR_DOTTED;
        isym[nline] = +GR_SQUARE;
        nper[nline] = header[2];
        nline++;

        EG_free(header);
        EG_free(gdata );

        grinit_(&io_kbd, &io_scr, "northeast", strlen("northeast"));
        grline_(ilin, isym, &nline, "~x~y~O=in, S=cp", &indgr, xplot, yplot, nper, strlen("~x~y~O=in, S=cp"));
    }
#endif

    status = EG_makeTopology(context, ecurve[0], EDGE, TWONODE,
                             tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
    printf("eedges[0]=%llx\n", (long long)(eedges[0]));
#endif

    /* return WireBody if NQUAD==1 */
    if (NQUAD(numUdp) == 1) {
        sense[0] = SFORWARD;
        status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                 NULL, 1, &(eedges[0]), sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

        goto cleanup;
    }

    /* nw quadrant -- generate coordinates */
    npnt = 0;

    if (SLPFACT(0) < 0 && n_nw >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(i) / (double)(NUMPNTS(0)) * PI/2;

            pnt_nw[3*npnt  ] = +rx_w * pow(cos(theta), 2.0/n_nw);
            pnt_nw[3*npnt+1] = +ry_n * pow(sin(theta), 2.0/n_nw);
            pnt_nw[3*npnt+2] = 0;
            npnt++;
        }
    }
    
    pnt_nw[3*npnt  ] =  0;
    pnt_nw[3*npnt+1] = +ry_n;
    pnt_nw[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) > 0) {
        pnt_nw[3*npnt  ] = -rx_w * SLPFACT(0);
        pnt_nw[3*npnt+1] = +ry_n;
        pnt_nw[3*npnt+2] = 0;
        npnt++;
    }

    for (i = 1; i < NUMPNTS(0)-1; i++) {
        theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

        pnt_nw[3*npnt  ] = -rx_w * pow(sin(theta), 2.0/n_nw);
        pnt_nw[3*npnt+1] = +ry_n * pow(cos(theta), 2.0/n_nw);
        pnt_nw[3*npnt+2] =  0;
        npnt++;
    }

    if (SLPFACT(0) > 0) {
        pnt_nw[3*npnt  ] = -rx_w;
        pnt_nw[3*npnt+1] = +ry_n * SLPFACT(0);
        pnt_nw[3*npnt+2] = 0;
        npnt++;
    }

    pnt_nw[3*npnt  ] = -rx_w;
    pnt_nw[3*npnt+1] =  0;
    pnt_nw[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) < 0 && n_nw >= 2) {
        for (i = 1; i < NUMPNTS(0); i++) {
            theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

            pnt_nw[3*npnt  ] = -rx_w * pow(cos(theta), 2.0/n_nw);
            pnt_nw[3*npnt+1] = -ry_n * pow(sin(theta), 2.0/n_nw);
            pnt_nw[3*npnt+2] = 0;
            npnt++;
        }
    }

    /* create offset if required */
    if (OFFSET(0) != 0) {
        for (i = 0; i < npnt; i++) {
            pnt_save[3*i  ] = pnt_nw[3*i  ];
            pnt_save[3*i+1] = pnt_nw[3*i+1];
        }

        pnt_nw[       1] += OFFSET(0);
        pnt_nw[3*npnt-3] -= OFFSET(0);

        for (i = 1; i < npnt-1; i++) {
            dx = pnt_save[3*i+3] - pnt_save[3*i-3];
            dy = pnt_save[3*i+4] - pnt_save[3*i-2];
            ds = sqrt(dx * dx + dy * dy);

            pnt_nw[3*i  ] += OFFSET(0) * dy / ds;
            pnt_nw[3*i+1] -= OFFSET(0) * dx / ds;
        }

#ifdef DEBUG
        for (i = 0; i < npnt; i++) {
            printf("nw %3d:  old: %10.4f %10.4f  new: %10.4f %10.4f\n",
                   i, pnt_save[3*i], pnt_save[3*i+1], pnt_nw[3*i], pnt_nw[3*i+1]);
        }
#endif
    }

    /* create spline curve and Edge */
    sizes[0] = npnt;
    sizes[1] = 0;
    status = EG_approximate(context, 1, dxytol, sizes, pnt_nw, &(ecurve[1]));
    CHECK_STATUS(EG_approximate);

    data[0] = 0;
    data[1] = +ry_n + OFFSET(0);
    data[2] = 0;
    status = EG_invEvaluate(ecurve[1], data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = -rx_w - OFFSET(0);
    data[1] = 0;
    data[2] = 0;
    status = EG_invEvaluate(ecurve[1], data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

#ifdef GRAFIC
    {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    ilin[3], isym[3], nper[3], nline=0, nplot=0, ipnt;
        int    oclass, mtype, *header;
        float  xplot[2000], yplot[2000];
        double tt, *gdata;
        ego    eref;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt_nw[3*ipnt  ];
            yplot[nplot] = pnt_nw[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = npnt;
        nline++;

        for (ipnt = 0; ipnt < 1001; ipnt++) {
            tt = tdata[0] + (tdata[1]-tdata[0]) * (double)(ipnt) / 1000.;
            status = EG_evaluate(ecurve[1], &tt, data);
            CHECK_STATUS(EG_evaluate);

            xplot[nplot] = data[0];
            yplot[nplot] = data[1];
            nplot++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] = 1001;
        nline++;

        status = EG_getGeometry(ecurve[1], &oclass, &mtype, &eref, &header, &gdata);
        CHECK_STATUS(EG_getGeometry);

        for (ipnt = 0; ipnt < header[2]; ipnt++) {
            xplot[nplot] = gdata[header[3]+3*ipnt  ];
            yplot[nplot] = gdata[header[3]+3*ipnt+1];
            nplot++;
        }
        ilin[nline] = +GR_DOTTED;
        isym[nline] = +GR_SQUARE;
        nper[nline] = header[2];
        nline++;

        EG_free(header);
        EG_free(gdata );

        grinit_(&io_kbd, &io_scr, "northwest", strlen("northwest"));
        grline_(ilin, isym, &nline, "~x~y~O=in, S=cp", &indgr, xplot, yplot, nper, strlen("~x~y~O=in, S=cp"));
    }
#endif

    status = EG_makeTopology(context, ecurve[1], EDGE, TWONODE,
                             tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
    printf("eedges[1]=%llx\n", (long long)(eedges[1]));
#endif

    /* return WireBody if NQUAD==2 */
    if (NQUAD(numUdp) == 2) {
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                 NULL, 2, &(eedges[0]), sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

        goto cleanup;
    }

    /* sw quadrant -- generate coordinates */
    npnt = 0;

    if (SLPFACT(0) < 0 && n_sw >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(NUMPNTS(0)-i) / (double)(NUMPNTS(0)) * PI/2;

            pnt_sw[3*npnt  ] = -rx_w * pow(cos(theta), 2.0/n_sw);
            pnt_sw[3*npnt+1] = +ry_s * pow(sin(theta), 2.0/n_sw);
            pnt_sw[3*npnt+2] = 0;
            npnt++;
        }
    }
    
    pnt_sw[3*npnt  ] = -rx_w;
    pnt_sw[3*npnt+1] =  0;
    pnt_sw[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) > 0) {
        pnt_sw[3*npnt  ] = -rx_w;
        pnt_sw[3*npnt+1] = -ry_s * SLPFACT(0);
        pnt_sw[3*npnt+2] = 0;
        npnt++;
    }

    for (i = 1; i < NUMPNTS(0)-1; i++) {
        theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

        pnt_sw[3*npnt  ] = -rx_w * pow(cos(theta), 2.0/n_sw);
        pnt_sw[3*npnt+1] = -ry_s * pow(sin(theta), 2.0/n_sw);
        pnt_sw[3*npnt+2] =  0;
        npnt++;
    }

    if (SLPFACT(0) > 0) {
        pnt_sw[3*npnt  ] = -rx_w * SLPFACT(0);
        pnt_sw[3*npnt+1] = -ry_s;
        pnt_sw[3*npnt+2] = 0;
        npnt++;
    }

    pnt_sw[3*npnt  ] =  0;
    pnt_sw[3*npnt+1] = -ry_s;
    pnt_sw[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) < 0 && n_sw >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(NUMPNTS(0)-1-i) / (double)(NUMPNTS(0)) * PI/2;

            pnt_sw[3*npnt  ] = +rx_w * pow(cos(theta), 2.0/n_sw);
            pnt_sw[3*npnt+1] = -ry_s * pow(sin(theta), 2.0/n_sw);
            pnt_sw[3*npnt+2] = 0;
            npnt++;
        }
    }

    /* create offset if required */
    if (OFFSET(0) != 0) {
        for (i = 0; i < npnt; i++) {
            pnt_save[3*i  ] = pnt_sw[3*i  ];
            pnt_save[3*i+1] = pnt_sw[3*i+1];
        }

        pnt_sw[       0] -= OFFSET(0);
        pnt_sw[3*npnt-2] -= OFFSET(0);

        for (i = 1; i < npnt-1; i++) {
            dx = pnt_save[3*i+3] - pnt_save[3*i-3];
            dy = pnt_save[3*i+4] - pnt_save[3*i-2];
            ds = sqrt(dx * dx + dy * dy);

            pnt_sw[3*i  ] += OFFSET(0) * dy / ds;
            pnt_sw[3*i+1] -= OFFSET(0) * dx / ds;
        }

#ifdef DEBUG
        for (i = 0; i < npnt; i++) {
            printf("sw %3d:  old: %10.4f %10.4f  new: %10.4f %10.4f\n",
                   i, pnt_save[3*i], pnt_save[3*i+1], pnt_sw[3*i], pnt_sw[3*i+1]);
        }
#endif
    }

    /* create spline curve and Edge */
    sizes[0] = npnt;
    sizes[1] = 0;
    status = EG_approximate(context, 1, dxytol, sizes, pnt_sw, &(ecurve[2]));
    CHECK_STATUS(EG_approximate);

    data[0] = -rx_w - OFFSET(0);
    data[1] = 0;
    data[2] = 0;
    status = EG_invEvaluate(ecurve[2], data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = 0;
    data[1] = -ry_s - OFFSET(0);
    data[2] = 0;
    status = EG_invEvaluate(ecurve[2], data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

#ifdef GRAFIC
    {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    ilin[3], isym[3], nper[3], nline=0, nplot=0, ipnt;
        int    oclass, mtype, *header;
        float  xplot[2000], yplot[2000];
        double tt, *gdata;
        ego    eref;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt_sw[3*ipnt  ];
            yplot[nplot] = pnt_sw[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = npnt;
        nline++;

        for (ipnt = 0; ipnt < 1001; ipnt++) {
            tt = tdata[0] + (tdata[1]-tdata[0]) * (double)(ipnt) / 1000.;
            status = EG_evaluate(ecurve[2], &tt, data);
            CHECK_STATUS(EG_evaluate);

            xplot[nplot] = data[0];
            yplot[nplot] = data[1];
            nplot++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] = 1001;
        nline++;

        status = EG_getGeometry(ecurve[2], &oclass, &mtype, &eref, &header, &gdata);
        CHECK_STATUS(EG_getGeometry);

        for (ipnt = 0; ipnt < header[2]; ipnt++) {
            xplot[nplot] = gdata[header[3]+3*ipnt  ];
            yplot[nplot] = gdata[header[3]+3*ipnt+1];
            nplot++;
        }
        ilin[nline] = +GR_DOTTED;
        isym[nline] = +GR_SQUARE;
        nper[nline] = header[2];
        nline++;

        EG_free(header);
        EG_free(gdata );

        grinit_(&io_kbd, &io_scr, "southwest", strlen("southwwest"));
        grline_(ilin, isym, &nline, "~x~y~O=in, S=cp", &indgr, xplot, yplot, nper, strlen("~x~y~O=in, S=cp"));
    }
#endif

    status = EG_makeTopology(context, ecurve[2], EDGE, TWONODE,
                             tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
    CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
    printf("eedges[2]=%llx\n", (long long)(eedges[2]));
#endif

    /* se quadrant -- generate coordinates */
    npnt = 0;

    if (SLPFACT(0) < 0 && n_se >= 2) {
        for (i = 0; i < NUMPNTS(0)-1; i++) {
            theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

            pnt_se[3*npnt  ] = -rx_e * pow(cos(theta), 2.0/n_se);
            pnt_se[3*npnt+1] = -ry_s * pow(sin(theta), 2.0/n_se);
            pnt_se[3*npnt+2] = 0;
            npnt++;
        }
    }
    
    pnt_se[3*npnt  ] =  0;
    pnt_se[3*npnt+1] = -ry_s;
    pnt_se[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) > 0) {
        pnt_se[3*npnt  ] = +rx_e * SLPFACT(0);
        pnt_se[3*npnt+1] = -ry_s;
        pnt_se[3*npnt+2] = 0;
        npnt++;
    }

    for (i = 1; i < NUMPNTS(0)-1; i++) {
        theta = (double)(i) / (double)(NUMPNTS(0)-1) * PI/2;

        pnt_se[3*npnt  ] = +rx_e * pow(sin(theta), 2.0/n_se);
        pnt_se[3*npnt+1] = -ry_s * pow(cos(theta), 2.0/n_se);
        pnt_se[3*npnt+2] =  0;
        npnt++;
    }

    if (SLPFACT(0) > 0) {
        pnt_se[3*npnt  ] = +rx_e;
        pnt_se[3*npnt+1] = -ry_s * SLPFACT(0);
        pnt_se[3*npnt+2] = 0;
        npnt++;
    }

    pnt_se[3*npnt  ] = +rx_e;
    pnt_se[3*npnt+1] =  0;
    pnt_se[3*npnt+2] =  0;
    npnt++;

    if (SLPFACT(0) < 0 && n_se >= 2) {
        for (i = 0; i < NUMPNTS(0); i++) {
            theta = (double)(i+1) / (double)(NUMPNTS(0)) * PI/2;

            pnt_se[3*npnt  ] = +rx_e * pow(cos(theta), 2.0/n_se);
            pnt_se[3*npnt+1] = +ry_s * pow(sin(theta), 2.0/n_se);
            pnt_se[3*npnt+2] = 0;
            npnt++;
        }
    }

    /* create offset if required */
    if (OFFSET(0) != 0) {
        for (i = 0; i < npnt; i++) {
            pnt_save[3*i  ] = pnt_se[3*i  ];
            pnt_save[3*i+1] = pnt_se[3*i+1];
        }

        pnt_se[       1] -= OFFSET(0);
        pnt_se[3*npnt-3] += OFFSET(0);

        for (i = 1; i < npnt-1; i++) {
            dx = pnt_save[3*i+3] - pnt_save[3*i-3];
            dy = pnt_save[3*i+4] - pnt_save[3*i-2];
            ds = sqrt(dx * dx + dy * dy);

            pnt_se[3*i  ] += OFFSET(0) * dy / ds;
            pnt_se[3*i+1] -= OFFSET(0) * dx / ds;
        }

#ifdef DEBUG
        for (i = 0; i < npnt; i++) {
            printf("se %3d:  old: %10.4f %10.4f  new: %10.4f %10.4f\n",
                   i, pnt_save[3*i], pnt_save[3*i+1], pnt_se[3*i], pnt_se[3*i+1]);
        }
#endif
    }

    /* create spline curve and Edge */
    sizes[0] = npnt;
    sizes[1] = 0;
    status = EG_approximate(context, 1, dxytol, sizes, pnt_se, &(ecurve[3]));
    CHECK_STATUS(EG_approximate);

    data[0] = 0;
    data[1] = -ry_s - OFFSET(0);
    data[2] = 0;
    status = EG_invEvaluate(ecurve[3], data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = +rx_e + OFFSET(0);
    data[1] = 0;
    data[2] = 0;
    status = EG_invEvaluate(ecurve[3], data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

#ifdef GRAFIC
    {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    ilin[3], isym[3], nper[3], nline=0, nplot=0, ipnt;
        int    oclass, mtype, *header;
        float  xplot[2000], yplot[2000];
        double tt, *gdata;
        ego    eref;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt_se[3*ipnt  ];
            yplot[nplot] = pnt_se[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = npnt;
        nline++;

        for (ipnt = 0; ipnt < 1001; ipnt++) {
            tt = tdata[0] + (tdata[1]-tdata[0]) * (double)(ipnt) / 1000.;
            status = EG_evaluate(ecurve[3], &tt, data);
            CHECK_STATUS(EG_evaluate);

            xplot[nplot] = data[0];
            yplot[nplot] = data[1];
            nplot++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] = 1001;
        nline++;

        status = EG_getGeometry(ecurve[3], &oclass, &mtype, &eref, &header, &gdata);
        CHECK_STATUS(EG_getGeometry);

        for (ipnt = 0; ipnt < header[2]; ipnt++) {
            xplot[nplot] = gdata[header[3]+3*ipnt  ];
            yplot[nplot] = gdata[header[3]+3*ipnt+1];
            nplot++;
        }
        ilin[nline] = +GR_DOTTED;
        isym[nline] = +GR_SQUARE;
        nper[nline] = header[2];
        nline++;

        EG_free(header);
        EG_free(gdata );

        grinit_(&io_kbd, &io_scr, "southeast", strlen("southeast"));
        grline_(ilin, isym, &nline, "~x~y~O=in, S=cp", &indgr, xplot, yplot, nper, strlen("~x~y~O=in, S=cp"));
    }
#endif

    status = EG_makeTopology(context, ecurve[3], EDGE, TWONODE,
                             tdata, 2, &(enodes[3]), NULL, &(eedges[3]));
    CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
    printf("eedges[3]=%llx\n", (long long)(eedges[3]));
#endif

    /* create loop of the four Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;
    sense[3] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, 4, eedges, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    /* make Face from the loop */
    status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
    CHECK_STATUS(EG_makeFace);

    /* since this will make a PLANE, we need to add an Attribute
       to tell OpenCSM to scale the UVs when computing sensitivities */
    status = EG_attributeAdd(eface, "_scaleuv", ATTRINT, 1,
                             &add, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* find the direction of the Face normal */
    status = EG_getRange(eface, range, &periodic);
    CHECK_STATUS(EG_getRange);

    range[0] = (range[0] + range[1]) / 2;
    range[1] = (range[2] + range[3]) / 2;

    status = EG_evaluate(eface, range, eval);
    CHECK_STATUS(EG_evaluate);

    norm[0] = eval[4] * eval[8] - eval[5] * eval[7];
    norm[1] = eval[5] * eval[6] - eval[3] * eval[8];
    norm[2] = eval[3] * eval[7] - eval[4] * eval[6];

    /* if the normal is not positive, flip the Face */
    if (norm[2] < 0) {
        status = EG_flipObject(eface, &enew);
        CHECK_STATUS(EG_flipObject);

        eface = enew;
    }

    /* create the FaceBody (which will be returned) */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                             NULL, 1, &eface, sense, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(pnt_ne  );
    FREE(pnt_nw  );
    FREE(pnt_sw  );
    FREE(pnt_se  );
    FREE(pnt_save);

    if (strlen(message) > 0) {
        *string = message;
        printf("%s\n", message);
    } else if (status != EGADS_SUCCESS) {
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

int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt, nnode, nedge, nface, nchild, oclass, mtype, *senses;
    double data[18];
    double rx_w,     rx_e,     ry_s,     ry_n;
    double rx_w_dot, rx_e_dot, ry_s_dot, ry_n_dot;
    ego    eref, *echilds, *enodes, *eedges, *efaces, eent;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
#endif

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    /* set up the parameters in each of the quadrants */
    rx_w     = rx_e     = 0;
    rx_w_dot = rx_e_dot = 0;
    ry_s     = ry_n     = 0;
    ry_s_dot = ry_n_dot = 0;

    if (RX(  iudp) > 0)  {
        rx_w     = rx_e     = RX(      iudp);
        rx_w_dot = rx_e_dot = RX_DOT(  iudp);
    }
    if (RX_W(iudp) > 0) {
        rx_w                = RX_W(    iudp);
        rx_w_dot            = RX_W_DOT(iudp);
    }
    if (RX_E(iudp) > 0) {
        rx_e                = RX_E(    iudp);
        rx_e_dot            = RX_E_DOT(iudp);
    }

    if (RY(  iudp) > 0) {
        ry_s     = ry_n     = RY(      iudp);
        ry_s_dot = ry_n_dot = RY_DOT(  iudp);
    }
    if (RY_S(iudp) > 0) {
        ry_s                = RY_S(    iudp);
        ry_s_dot            = RY_S_DOT(iudp);
    }
    if (RY_N(iudp) > 0) {
        ry_n                = RY_N(    iudp);
        ry_n_dot            = RY_N_DOT(iudp);
    }

    /* find the ego entity */
    if (entType == OCSM_NODE) {
        status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        eent = enodes[entIndex-1];

        EG_free(enodes);
    } else if (entType == OCSM_EDGE) {
        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        CHECK_STATUS(EG_getBodyTopos);

        eent = eedges[entIndex-1];

        EG_free(eedges);
    } else if (entType == OCSM_FACE) {
        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        eent = efaces[entIndex-1];

        EG_free(efaces);
    } else {
        printf("udpSensitivity: bad entType=%d\n", entType);
        status = EGADS_ATTRERR;
        goto cleanup;
    }

    /* loop through the points */
    for (ipnt = 0; ipnt < npnt; ipnt++) {

        /* find the physical coordinates */
        if        (entType == OCSM_NODE) {
            status = EG_getTopology(eent, &eref, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);
        } else if (entType == OCSM_EDGE) {
            status = EG_evaluate(eent, &(uvs[ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        } else if (entType == OCSM_FACE) {
            status = EG_evaluate(eent, &(uvs[2*ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        }

        /* compute the sensitivity */
        if (data[0] >= 0 && data[1] >= 0) {
            vels[3*ipnt  ] = data[0] / rx_e * rx_e_dot;
            vels[3*ipnt+1] = data[1] / ry_n * ry_n_dot;
            vels[3*ipnt+2] = 0;
        } else if (data[1] >= 0) {
            vels[3*ipnt  ] = data[0] / rx_w * rx_w_dot;
            vels[3*ipnt+1] = data[1] / ry_n * ry_n_dot;
            vels[3*ipnt+2] = 0;
        } else if (data[0] <= 0) {
            vels[3*ipnt  ] = data[0] / rx_w * rx_w_dot;
            vels[3*ipnt+1] = data[1] / ry_s * ry_s_dot;
            vels[3*ipnt+2] = 0;
        } else {
            vels[3*ipnt  ] = data[0] / rx_e * rx_e_dot;
            vels[3*ipnt+1] = data[1] / ry_s * ry_s_dot;
            vels[3*ipnt+2] = 0;
        }
    }

    status = EGADS_SUCCESS;

    /* velocities for Nodes */
//$$$    if (entType == OCSM_NODE) {
//$$$        if        (entIndex == 1) {
//$$$            vels[0] = +rx_e_dot;
//$$$            vels[1] =  0;
//$$$            vels[2] =  0;
//$$$        } else if (entIndex == 2) {
//$$$            vels[0] =  0;
//$$$            vels[1] = +ry_n_dot;
//$$$            vels[2] =  0;
//$$$        } else if (entIndex == 3) {
//$$$            vels[0] = -rx_w_dot;
//$$$            vels[1] =  0;
//$$$            vels[2] =  0;
//$$$        } else if (entIndex == 4) {
//$$$            vels[0] =  0;
//$$$            vels[1] = -ry_s_dot;
//$$$            vels[2] =  0;
//$$$        } else {
//$$$            printf(" udpSensitivity: bad node index\n");
//$$$            return EGADS_INDEXERR;
//$$$        }
//$$$
    /* velocities for Edges */
//$$$    } else if (entType == OCSM_EDGE) {
//$$$        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
//$$$        if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$        if (entIndex >=1 && entIndex <= 4) {
//$$$            eedge = eedges[entIndex-1];
//$$$        } else {
//$$$            EG_free(eedges);
//$$$            printf(" udpSensitivity: bad edge index\n");
//$$$            return EGADS_INDEXERR;
//$$$        }
//$$$
//$$$        EG_free(eedges);
//$$$        eedges = NULL;
//$$$
//$$$        if        (entIndex == 1) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedge, &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                vels[3*ipnt  ] = rx_e_dot * xyz[0] / rx_e;
//$$$                vels[3*ipnt+1] = ry_n_dot * xyz[1] / ry_n;
//$$$                vels[3*ipnt+2] = 0;
//$$$
//$$$                if (fabs(n_ne_dot) > EPS06) {
//$$$                    theta = atan(pow((+rx_e*xyz[1])/(+ry_n*xyz[0]), n_ne/2));
//$$$                    if (theta > EPS06 && theta < PIo2-EPS06) {
//$$$                        vels[3*ipnt  ] -= n_ne_dot * 2 * xyz[0] * log(cos(theta)) / (n_ne*n_ne);
//$$$                        vels[3*ipnt+1] -= n_ne_dot * 2 * xyz[1] * log(sin(theta)) / (n_ne*n_ne);
//$$$                    }
//$$$                }
//$$$            }
//$$$        } else if (entIndex == 2) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedge, &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                vels[3*ipnt  ] = rx_w_dot * xyz[0] / rx_w;
//$$$                vels[3*ipnt+1] = ry_n_dot * xyz[1] / ry_n;
//$$$                vels[3*ipnt+2] = 0;
//$$$
//$$$                if (fabs(n_nw_dot) >EPS06) {
//$$$                    theta = atan(pow((-rx_w*xyz[1])/(+ry_n*xyz[0]), n_nw/2));
//$$$                    if (theta > EPS06 && theta < PIo2-EPS06) {
//$$$                        vels[3*ipnt  ] -= n_nw_dot * 2 * xyz[0] * log(cos(theta)) / (n_nw*n_nw);
//$$$                        vels[3*ipnt+1] -= n_nw_dot * 2 * xyz[1] * log(sin(theta)) / (n_nw*n_nw);
//$$$                    }
//$$$                }
//$$$            }
//$$$        } else if (entIndex == 3) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedge, &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                vels[3*ipnt  ] = rx_w_dot * xyz[0] / rx_w;
//$$$                vels[3*ipnt+1] = ry_s_dot * xyz[1] / ry_s;
//$$$                vels[3*ipnt+2] = 0;
//$$$
//$$$                if (fabs(n_sw_dot) > EPS06) {
//$$$                    theta = atan(pow((-rx_w*xyz[1])/(-ry_s*xyz[0]), n_sw/2));
//$$$                    if (theta > EPS06 && theta < PIo2-EPS06) {
//$$$                        vels[3*ipnt  ] -= n_sw_dot * 2 * xyz[0] * log(cos(theta)) / (n_sw*n_sw);
//$$$                        vels[3*ipnt+1] -= n_sw_dot * 2 * xyz[1] * log(sin(theta)) / (n_sw*n_sw);
//$$$                    }
//$$$                }
//$$$            }
//$$$        } else if (entIndex == 4) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedge, &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                vels[3*ipnt  ] = rx_e_dot * xyz[0] / rx_e;
//$$$                vels[3*ipnt+1] = ry_s_dot * xyz[1] / ry_s;
//$$$                vels[3*ipnt+2] = 0;
//$$$
//$$$                if (fabs(n_se_dot) > EPS06) {
//$$$                    theta = atan(pow((+rx_e*xyz[1])/(-ry_s*xyz[0]), n_se/2));
//$$$                    if (theta > EPS06 && theta < PIo2-EPS06) {
//$$$                        vels[3*ipnt  ] -= n_se_dot * 2 * xyz[0] * log(cos(theta)) / (n_se*n_se);
//$$$                        vels[3*ipnt+1] -= n_se_dot * 2 * xyz[1] * log(sin(theta)) / (n_se*n_se);
//$$$                    }
//$$$                }
//$$$            }
//$$$        } else {
//$$$            printf(" udpSensitivity: bad index\n");
//$$$            return EGADS_INDEXERR;
//$$$        }
//$$$
    /* velocities for Face */
//$$$    } else if (entType == OCSM_FACE) {
//$$$        if (entIndex == 1) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                vels[3*ipnt  ] = 0;
//$$$                vels[3*ipnt+1] = 0;
//$$$                vels[3*ipnt+2] = 0;
//$$$            }
//$$$        } else {
//$$$            printf(" udpSensitivity: bad Face index\n");
//$$$            return EGADS_INDEXERR;
//$$$        }
//$$$
//$$$    } else {
//$$$        printf(" udpSensitivity: bad entity type\n");
//$$$        return EGADS_GEOMERR;
//$$$    }
//$$$
//$$$    status = EGADS_SUCCESS;

cleanup:
    return status;
}
