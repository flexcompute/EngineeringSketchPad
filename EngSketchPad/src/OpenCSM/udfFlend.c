/*
 ************************************************************************
 *                                                                      *
 * udfFlend -- create a flended Body                                    *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2025  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPINPUTBODYS -2
#define NUMUDPARGS        7

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define SLOPEA(IUDP)   ((double *) (udps[IUDP].arg[0].val))[0]
#define SLOPEB(IUDP)   ((double *) (udps[IUDP].arg[1].val))[0]
#define TOLER( IUDP)   ((double *) (udps[IUDP].arg[2].val))[0]
#define EQUIS( IUDP)   ((int    *) (udps[IUDP].arg[3].val))[0]
#define NPNT(  IUDP)   ((int    *) (udps[IUDP].arg[4].val))[0]
#define PLOT(  IUDP)   ((int    *) (udps[IUDP].arg[5].val))[0]
#define METHOD(IUDP)   ((int    *) (udps[IUDP].arg[6].val))[0]

static char  *argNames[NUMUDPARGS] = {"slopea",  "slopeb",  "toler",  "equis", "npnt",  "plot",  "method", };
static int    argTypes[NUMUDPARGS] = {ATTRREAL,  ATTRREAL,  ATTRREAL, ATTRINT, ATTRINT, ATTRINT, ATTRINT,  };
static int    argIdefs[NUMUDPARGS] = {0,         0,         0,        0,       33,      0,       1,        };
static double argDdefs[NUMUDPARGS] = {1.00,      1.00,      1.0e-6,   0.,      0.,      0.,      0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

/* functions defined below */
static int exposedLoops(ego ebody, int *nloop, ego eloops[]);
static int fillPointsFromEdge(ego eedgeA, int senseA,
                              ego eedgeB, int senseB,
                              int npnt, double tA[], double pntA[], double tB[], double pntB[], int *NumUdp, udp_T *udps);
static int fillSlopesFromEdge(ego eedge, int sense, ego eface,
                              int npnt, double t[], double slp[]);
static int reorderLoop(ego eloopA, ego *eloopB);
static int slopeAtNode(ego eedge, int sense, ego eface1, ego eface2, ego ebody, double newSlope[]);

/* undocumented EGADS function that is used in BLEND */
extern int EG_spline2dAppx(      ego    context,      // EGADS context
                                 int    endc,         // =0 natural, =1 extrap, =2 dont use
                /*@null@*/ const double *uknot,       // knots in u direction     (imax+6)
                /*@null@*/ const double *vknot,       // knots in v direction     (jmax+6)
                /*@null@*/ const int    *vdata,       // used to set C0 and C1
                /*@null@*/ const double *wesT,        // tangency at imin         (3*jmax)
                /*@null@*/ const double *easT,        // tangency at imax         (3*jmax)
                /*@null@*/ const double *south,       // data for nose treatment
                /*@null@*/       double *snor,        // tangency at jmin         (3*imax)
                /*@null@*/ const double *north,       // data for tail treatment
                /*@null@*/       double *nnor,        // tangency at jmax         (3*imax)
                                 int    imax,
                                 int    jmax,
                                 const  double *xyz,  // x00, y00, z00, x10, ...  (3*imax*jmax)
                                 double tol,          // tolerance (or 0)
                                 ego    *esurf);


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  emodel,                 /* (in)  input model */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;     /* (in)  return status */

    ego     context, *ebodys;

    int     nedgeA, nfaceA, ifaceA, nloopA, *sensesA;
    int     nedgeB, nfaceB, ifaceB, nloopB, *sensesB;
    int     ipnt, nfacelist, brchattr[2], bodyattr[2];
    int     oclass, mtype, nchild, *senses, ntemp, atype, alen, itemp, iedge, jedge, npnt;
    int     senseA, senseB, nloop, iloop, indx, ij3;
    int     ismth, nsmth=100;
    CINT    *tempIlist;
    double  data[18], value, dot, fraci, dist, uvFace[2];
    double  oldSlope[3], newSlope[4], toler=1.0e-8;
    double  normAi[3], normAj[3], tangA[3], normBi[3], normBj[3], tangB[3];
    double  det, t1, t2, rhsA, rhsB, dt, dtmax;
    double  *spln=NULL, *west=NULL, *east=NULL;
    double  *tA=NULL, *pntA=NULL, *slpA=NULL, *tB=NULL, *pntB=NULL, *slpB=NULL;
    double  *pntC=NULL, *pntD=NULL, *dCdt=NULL, *tD=NULL;
    CDOUBLE *tempRlist;
    char    temp[1], *message=NULL;
    CCHAR   *tempClist;
    ego     ebodyA, eloopsA[2], *eedgesA, *efacesA=NULL, *efacesA2=NULL;
    ego     ebodyB, eloopsB[2], *eedgesB, *efacesB=NULL, *efacesB2=NULL;
    ego     eref, esurf, *etemps, *echilds, *eloops, emodel2, *efacelist=NULL;
    void    *modl;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("slopea(0) = %f\n", SLOPEA(0));
    printf("slopeb(0) = %f\n", SLOPEB(0));
    printf("toler( 0) = %f\n", TOLER( 0));
    printf("equis( 0) = %d\n", EQUIS( 0));
    printf("npnt ( 0) = %d\n", NPNT(  0));
    printf("plot(  0) = %d\n", PLOT(  0));
    printf("method(0) = %d\n", METHOD(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check/process arguments */
    if (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "slopea should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "slopeb should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "toler should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (TOLER(0) < 0) {
        snprintf(message, 100, "toler = %f <= 0", TOLER(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "\"plot\" should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[6].size > 1) {
        snprintf(message, 100, "\"method\" should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (NPNT(0) < 5) {
        snprintf(message, 100, "npnt = %d < 5", NPNT(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    }

#ifdef DEBUG
    printf("emodel:\n");
    ocsmPrintEgo(emodel);
#endif

    /* check that Model that was input that contains one or two Bodys */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1 && nchild != 2) {
        snprintf(message, 100, "expecting Model to contain one or two Bodys (not %d)", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("slopea(%d) = %f\n", numUdp, SLOPEA(numUdp));
    printf("slopeb(%d) = %f\n", numUdp, SLOPEB(numUdp));
    printf("toler( %d) = %f\n", numUdp, TOLER( numUdp));
    printf("equis( %d) = %d\n", numUdp, EQUIS( numUdp));
    printf("npnt(  %d) = %d\n", numUdp, NPNT(  numUdp));
    printf("plot(  %d) = %d\n", numUdp, PLOT(  numUdp));
    printf("method(%d) = %d\n", numUdp, METHOD(numUdp));
#endif

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* get pointer to the OpenCSM MODL */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    /* set up if we only have one Body */
    if (nchild == 1) {
        status = EG_copyObject(ebodys[0], NULL, &ebodyA);
        CHECK_STATUS(EG_copyObject);

        ebodyB = ebodyA;

        status = exposedLoops(ebodyA, &nloopA, eloopsA);
        CHECK_STATUS(exposedLoops);

        if (nloopA == 2) {
            eloopsB[0] = eloopsA[1];
        } else {
            snprintf(message, 100, "FLEND found BodyA contains %d Loops (expecting 2)", nloopA);
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

    /* set up if we have two Bodys */
    } else if (nchild == 2) {
        status = EG_copyObject(ebodys[0], NULL, &ebodyA);
        CHECK_STATUS(EG_copyObject);

        status = EG_copyObject(ebodys[1], NULL, &ebodyB);
        CHECK_STATUS(EG_copyObject);

        /* make a list of the exposed Loops in BodyA (when Faces with the
           _flend=remove attribute are removed) */
        status = exposedLoops(ebodyA, &nloopA, eloopsA);
        CHECK_STATUS(exposedLoops);

        if (nloopA != 1) {
            snprintf(message, 100, "FLEND found BodyA contains %d Loops (expecting 1)", nloopA);
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

        /* make a list of the exposed Loops in BodyB (when Faces with the
           _flend=remove attribute are removed) */
        status = exposedLoops(ebodyB, &nloopB, eloopsB);
        CHECK_STATUS(exposedLoops);

        if (nloopB != 1) {
            snprintf(message, 100, "FLEND found BodyB contains %d Loops (expecting 1)", nloopB);
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

    /* error */
    } else {
        snprintf(message, 100, "FLEND found Model contains %d Bodys (expecting 1 or 2)", nchild);
        status = OCSM_UDP_ERROR1;
        goto cleanup;
    }

    /* make sure both Loops have the same number of Edges */
    status = EG_getTopology(eloopsA[0], &eref, &oclass, &mtype,
                            data, &nedgeA, &eedgesA, &sensesA);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(eloopsB[0], &eref, &oclass, &mtype,
                            data, &nedgeB, &eedgesB, &sensesB);
    CHECK_STATUS(EG_getTopology);

    if (nedgeA != nedgeB) {
        printf("eloopsA[0]\n");
        ocsmPrintEgo(eloopsA[0]);
        printf("eloopsB[0]\n");
        ocsmPrintEgo(eloopsB[0]);

        snprintf(message, 100, "nedgeA=%d does not match nedgeB=%d", nedgeA, nedgeB);
        status = OCSM_UDP_ERROR2;
        goto cleanup;
    }

    /* reorder eloopsB[0] to minimize twist from eloopsA[0] */
    status = reorderLoop(eloopsA[0], eloopsB);
    CHECK_STATUS(reorderLoop);

#ifdef DEBUG
    printf("eloopsA[0] (after reordering)\n");
    ocsmPrintEgo(eloopsA[0]);
    printf("eloopsB[0] (after reordering)\n");
    ocsmPrintEgo(eloopsB[0]);
#endif

    /* get Edges for bounding loops (now that they have possibly been reordered) */
    status = EG_getTopology(eloopsA[0], &eref, &oclass, &mtype, data, &nedgeA, &eedgesA, &sensesA);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(eloopsB[0], &eref, &oclass, &mtype, data, &nedgeB, &eedgesB, &sensesB);
    CHECK_STATUS(EG_getTopology);

    /* get the Faces adjacent to the Edges in eedgesA */
    MALLOC(efacesA, ego, nedgeA);
    MALLOC(efacesB, ego, nedgeB);

    for (iedge = 0; iedge < nedgeA; iedge++) {
        efacesA[iedge] = NULL;

        status = EG_getBodyTopos(ebodyA, eedgesA[iedge], FACE, &ntemp, &etemps);
        CHECK_STATUS(EG_getBodyTopos);

        for (itemp = 0; itemp < ntemp; itemp++) {
            status = EG_attributeRet(etemps[itemp], "_flend", &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != EGADS_SUCCESS) {
                efacesA[iedge] = etemps[itemp];
                break;
            }
        }

        EG_free(etemps);

        if (efacesA[iedge] == NULL) {
            snprintf(message, 100, "eedgesA[%d] is not adjecent to one Face", iedge);
            status = OCSM_UDP_ERROR3;
            goto cleanup;
        }
    }

    /* get the Faces adjacent to the Edges in eedgesB */
    for (iedge = 0; iedge < nedgeB; iedge++) {
        efacesB[iedge] = NULL;

        status = EG_getBodyTopos(ebodyB, eedgesB[iedge], FACE, &ntemp, &etemps);
        CHECK_STATUS(EG_getBodyTopos);

        for (itemp = 0; itemp < ntemp; itemp++) {
            status = EG_attributeRet(etemps[itemp], "_flend", &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != EGADS_SUCCESS) {
                efacesB[iedge] = etemps[itemp];
                break;
            }
        }

        EG_free(etemps);

        if (efacesB[iedge] == NULL) {
            snprintf(message, 100, "eedgesB[%d] is not adjacent to one Face", iedge);
            status = OCSM_UDP_ERROR4;
            goto cleanup;
        }
    }

    /* get arrays to hold the west and east points and slopes */
    npnt = NPNT(0);

    MALLOC(tA,   double,   npnt*nedgeA);
    MALLOC(pntA, double, 3*npnt*nedgeA);
    MALLOC(slpA, double, 3*npnt*nedgeA);
    MALLOC(tB,   double,   npnt*nedgeA);
    MALLOC(pntB, double, 3*npnt*nedgeA);
    MALLOC(slpB, double, 3*npnt*nedgeA);
    MALLOC(pntC, double, 3*npnt*nedgeA);
    MALLOC(pntD, double, 3*npnt*nedgeA);
    MALLOC(tD,   double,   npnt       );
    MALLOC(dCdt, double, 3*npnt*nedgeA);

    /* get the (initial) Points from the Loops */
    for (iedge = 0; iedge < nedgeA; iedge++) {
        senseA = 0;
        senseB = 0;

        /* determine the sense of the Edge relative to efaces[iedge] */
        status = EG_getTopology(efacesA[iedge], &eref, &oclass, &mtype,
                                data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        for (iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            for (itemp = 0; itemp < ntemp; itemp++) {
                if (eedgesA[iedge] == etemps[itemp]) {
                    senseA = senses[itemp];
                    break;
                }
            }
            if (senseA != 0) break;
        }

        status = EG_getTopology(efacesB[iedge], &eref, &oclass, &mtype,
                                data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        for (iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            for (itemp = 0; itemp < ntemp; itemp++) {
                if (eedgesB[iedge] == etemps[itemp]) {
                    senseB = senses[itemp];
                    break;
                }
            }
            if (senseB != 0) break;
        }

        /* get the (initial) Points from the Loops */
        status = fillPointsFromEdge(eedgesA[iedge], sensesA[iedge],
                                    eedgesB[iedge], sensesB[iedge], npnt,
                                    &tA[npnt*iedge], &pntA[3*npnt*iedge],
                                    &tB[npnt*iedge], &pntB[3*npnt*iedge], NumUdp, udps);
        CHECK_STATUS(fillPointsFromEdge);

        /* get the (initial) slopes from the Loops */
#ifdef DEBUG
        printf("BodyA, iedge=%d (senseA=%d, iedge=%d, iface=%d)\n", iedge, senseA,
               EG_indexBodyTopo(ebodyA, eedgesA[iedge]), EG_indexBodyTopo(ebodyA, efacesA[iedge]));
#endif

        status = fillSlopesFromEdge(eedgesA[iedge], senseA, efacesA[iedge],
                                    npnt, &tA[npnt*iedge], &slpA[3*npnt*iedge]);
        CHECK_STATUS(fillSlopesFromEdge);

#ifdef DEBUG
        printf("BodyB, iedge=%d (senseB=%d, iedge=%d, iface=%d)\n", iedge, senseB,
               EG_indexBodyTopo(ebodyB, eedgesB[iedge]), EG_indexBodyTopo(ebodyB, efacesB[iedge]));
#endif

        status = fillSlopesFromEdge(eedgesB[iedge], senseB, efacesB[iedge],
                                    npnt, &tB[npnt*iedge], &slpB[3*npnt*iedge]);
        CHECK_STATUS(fillSlopesFromEdge);
    }

    /* find the slope at the Nodes and adjust the slopes for the adjacent Edges */
    if (METHOD(0) == 1) {
        for (iedge = 0; iedge < nedgeA; iedge++) {
            jedge = (nedgeA + iedge - 1) % nedgeA;

            /* find A slope at corner */
            if (efacesA[iedge] != efacesA[jedge]) {
                status = slopeAtNode(eedgesA[iedge], sensesA[iedge], efacesA[iedge], efacesA[jedge],
                                     ebodyA, newSlope);
                CHECK_STATUS(slopeAtNode);
            } else {
                newSlope[0] = slpA[3*(       iedge*npnt)  ] + slpA[3*(npnt-1+jedge*npnt)  ];
                newSlope[1] = slpA[3*(       iedge*npnt)+1] + slpA[3*(npnt-1+jedge*npnt)+1];
                newSlope[2] = slpA[3*(       iedge*npnt)+2] + slpA[3*(npnt-1+jedge*npnt)+2];

                newSlope[3] = sqrt(newSlope[0]*newSlope[0] + newSlope[1]*newSlope[1] + newSlope[2]*newSlope[2]);
                newSlope[0] /= newSlope[3];
                newSlope[1] /= newSlope[3];
                newSlope[2] /= newSlope[3];
            }

            /* adjust A slopes on iedge */
            oldSlope[0] = slpA[3*(iedge*npnt)  ];
            oldSlope[1] = slpA[3*(iedge*npnt)+1];
            oldSlope[2] = slpA[3*(iedge*npnt)+2];

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = (double)(ipnt) / (double)(npnt-1);

                slpA[3*(ipnt+iedge*npnt)  ] += (1-fraci) * (newSlope[0] - oldSlope[0]);
                slpA[3*(ipnt+iedge*npnt)+1] += (1-fraci) * (newSlope[1] - oldSlope[1]);
                slpA[3*(ipnt+iedge*npnt)+2] += (1-fraci) * (newSlope[2] - oldSlope[2]);
            }

            /* adjust A slopes on jedge */
            oldSlope[0] = slpA[3*(npnt-1+jedge*npnt)  ];
            oldSlope[1] = slpA[3*(npnt-1+jedge*npnt)+1];
            oldSlope[2] = slpA[3*(npnt-1+jedge*npnt)+2];

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = (double)(ipnt) / (double)(npnt-1);

                slpA[3*(ipnt+jedge*npnt)  ] += (  fraci) * (newSlope[0] - oldSlope[0]);
                slpA[3*(ipnt+jedge*npnt)+1] += (  fraci) * (newSlope[1] - oldSlope[1]);
                slpA[3*(ipnt+jedge*npnt)+2] += (  fraci) * (newSlope[2] - oldSlope[2]);
            }

            /* find B slope at corner */
            if (efacesB[iedge] != efacesB[jedge]) {
                status = slopeAtNode(eedgesB[iedge], sensesB[iedge], efacesB[iedge], efacesB[jedge],
                                     ebodyB, newSlope);
                CHECK_STATUS(slopeAtNode);
            } else {
                newSlope[0] = slpB[3*(       iedge*npnt)  ] + slpB[3*(npnt-1+jedge*npnt)  ];
                newSlope[1] = slpB[3*(       iedge*npnt)+1] + slpB[3*(npnt-1+jedge*npnt)+1];
                newSlope[2] = slpB[3*(       iedge*npnt)+2] + slpB[3*(npnt-1+jedge*npnt)+2];

                newSlope[3] = sqrt(newSlope[0]*newSlope[0] + newSlope[1]*newSlope[1] + newSlope[2]*newSlope[2]);
                newSlope[0] /= newSlope[3];
                newSlope[1] /= newSlope[3];
                newSlope[2] /= newSlope[3];
            }

            /* adjust B slopes on iedge */
            oldSlope[0] = slpB[3*(iedge*npnt)  ];
            oldSlope[1] = slpB[3*(iedge*npnt)+1];
            oldSlope[2] = slpB[3*(iedge*npnt)+2];

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = (double)(ipnt) / (double)(npnt-1);

                slpB[3*(ipnt+iedge*npnt)  ] += (1-fraci) * (newSlope[0] - oldSlope[0]);
                slpB[3*(ipnt+iedge*npnt)+1] += (1-fraci) * (newSlope[1] - oldSlope[1]);
                slpB[3*(ipnt+iedge*npnt)+2] += (1-fraci) * (newSlope[2] - oldSlope[2]);
            }

            /* adjust B slopes on jedge */
            oldSlope[0] = slpB[3*(npnt-1+jedge*npnt)  ];
            oldSlope[1] = slpB[3*(npnt-1+jedge*npnt)+1];
            oldSlope[2] = slpB[3*(npnt-1+jedge*npnt)+2];

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = (double)(ipnt) / (double)(npnt-1);

                slpB[3*(ipnt+jedge*npnt)  ] += (  fraci) * (newSlope[0] - oldSlope[0]);
                slpB[3*(ipnt+jedge*npnt)+1] += (  fraci) * (newSlope[1] - oldSlope[1]);
                slpB[3*(ipnt+jedge*npnt)+2] += (  fraci) * (newSlope[2] - oldSlope[2]);
            }
        }

        /* modify the slopes by the distances across the flend */
        for (iedge = 0; iedge < nedgeA; iedge++) {
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                indx = 3 * (ipnt + iedge * npnt);

                dist = sqrt((pntA[indx  ]-pntB[indx  ]) * (pntA[indx  ]-pntB[indx  ])
                           +(pntA[indx+1]-pntB[indx+1]) * (pntA[indx+1]-pntB[indx+1])
                           +(pntA[indx+2]-pntB[indx+2]) * (pntA[indx+2]-pntB[indx+2]));

                slpA[indx  ] *= dist;
                slpA[indx+1] *= dist;
                slpA[indx+2] *= dist;

                slpB[indx  ] *= dist;
                slpB[indx+1] *= dist;
                slpB[indx+2] *= dist;
            }
        }

    } else if (METHOD(0) == 2) {

        /* pntA is at tangent to faceA */
        /* pntB is at tangent to faceB */
        /* pntC is arbitrary point on intersection of faceA and faceB */
        /* dCdt is vector through C that points in intersection direction */
        /* pntD is evenly spaced points along intersection of faceA and faceB */

        /* first set up points on Nodes (at the end of iedge) */
        for (iedge = 0; iedge < nedgeA; iedge++) {
            jedge = (iedge + 1) % nedgeA;
            ipnt  = npnt - 1;
            ij3   = 3 * (ipnt + iedge * npnt);

#ifdef DEBUG
            printf("iedge=%d, jedge=%d\n", iedge, jedge);
#endif

            /* normalized normal at pntA for iedge */
            status = EG_getEdgeUV(efacesA[iedge], eedgesA[iedge], 0, tA[ipnt+iedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesA[iedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            normAi[0]  = data[4] * data[8] - data[5] * data[7];
            normAi[1]  = data[5] * data[6] - data[3] * data[8];
            normAi[2]  = data[3] * data[7] - data[4] * data[6];
            det        = sqrt(normAi[0]*normAi[0] + normAi[1]*normAi[1] + normAi[2]*normAi[2]);
            normAi[0] /= det;
            normAi[1] /= det;
            normAi[2] /= det;

            /* iedge and jedge are associated with the same Face on the A side */
            if (efacesA[iedge] == efacesA[jedge]) {
                normAj[0] = normAi[0];
                normAj[1] = normAi[1];
                normAj[2] = normAi[2];
                tangA[ 0] = 0;
                tangA[ 1] = 0;
                tangA[ 2] = 0;

            } else {

                /* normalized normal at pntA for jedge */
                status = EG_getEdgeUV(efacesA[jedge], eedgesA[jedge], 0, tA[jedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesA[jedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                normAj[0]  = data[4] * data[8] - data[5] * data[7];
                normAj[1]  = data[5] * data[6] - data[3] * data[8];
                normAj[2]  = data[3] * data[7] - data[4] * data[6];
                det        = sqrt(normAj[0]*normAj[0] + normAj[1]*normAj[1] + normAj[2]*normAj[2]);
                normAj[0] /= det;
                normAj[1] /= det;
                normAj[2] /= det;

                /* if nearly parallel, use average norm */
                if (fabs(normAi[0]*normAj[0]+normAi[1]*normAj[1]+normAi[2]*normAj[2]-1) < EPS06) {
                    normAi[0] = normAj[0] = (normAi[0] + normAj[0]) / 2;
                    normAi[1] = normAj[1] = (normAi[1] + normAj[1]) / 2;
                    normAi[2] = normAj[2] = (normAi[2] + normAj[2]) / 2;
                    tangA[ 0] = 0;
                    tangA[ 1] = 0;
                    tangA[ 2] = 0;

                /* otherwise set up tangent that is perpendicular to normAi and normAj */
                } else {
                    tangA[0]  = normAi[1] * normAj[2] - normAj[1] * normAi[2];
                    tangA[1]  = normAi[2] * normAj[0] - normAj[2] * normAi[0];
                    tangA[2]  = normAi[0] * normAj[1] - normAj[0] * normAi[1];
                    det       = sqrt(tangA[0]*tangA[0] + tangA[1]*tangA[1] + tangA[2]*tangA[2]);
                    tangA[0] /= det;
                    tangA[1] /= det;
                    tangA[2] /= det;
                }
            }

#ifdef DEBUG
            printf("   pntA  =%10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
            printf("   normAi=%10.5f %10.5f %10.5f\n", normAi[0], normAi[1], normAi[2]);
            printf("   normAj=%10.5f %10.5f %10.5f\n", normAj[0], normAj[1], normAj[2]);
            printf("   tangA =%10.5f %10.5f %10.5f\n", tangA[ 0], tangA[ 1], tangA[ 2]);
#endif

            /* normalized normal at pntB for iedge */
            status = EG_getEdgeUV(efacesB[iedge], eedgesB[iedge], 0, tB[ipnt+iedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesB[iedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            normBi[0]  = data[4] * data[8] - data[5] * data[7];
            normBi[1]  = data[5] * data[6] - data[3] * data[8];
            normBi[2]  = data[3] * data[7] - data[4] * data[6];
            det        = sqrt(normBi[0]*normBi[0] + normBi[1]*normBi[1] + normBi[2]*normBi[2]);
            normBi[0] /= det;
            normBi[1] /= det;
            normBi[2] /= det;

            /* iedge and jedge are associated with the same Face on the B side */
            if (efacesB[iedge] == efacesB[jedge]) {
                normBj[0] = normBi[0];
                normBj[1] = normBi[1];
                normBj[2] = normBi[2];
                tangB[ 0] = 0;
                tangB[ 1] = 0;
                tangB[ 2] = 0;

            } else {

                /* normalized normal at pntB for jedge */
                status = EG_getEdgeUV(efacesB[jedge], eedgesB[jedge], 0, tB[jedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesB[jedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                normBj[0]  = data[4] * data[8] - data[5] * data[7];
                normBj[1]  = data[5] * data[6] - data[3] * data[8];
                normBj[2]  = data[3] * data[7] - data[4] * data[6];
                det        = sqrt(normBj[0]*normBj[0] + normBj[1]*normBj[1] + normBj[2]*normBj[2]);
                normBj[0] /= det;
                normBj[1] /= det;
                normBj[2] /= det;

                /* if nearly parallel, use average norm */
                if (fabs(normBi[0]*normBj[0]+normBi[1]*normBj[1]+normBi[2]*normBj[2]-1) < EPS06) {
                    normBi[0] = normBj[0] = (normBi[0] + normBj[0]) / 2;
                    normBi[1] = normBj[1] = (normBi[1] + normBj[1]) / 2;
                    normBi[2] = normBj[2] = (normBi[2] + normBj[2]) / 2;
                    tangB[ 0] = 0;
                    tangB[ 1] = 0;
                    tangB[ 2] = 0;

                /* otherwise set up tangent that is perpendicular to normBi and normBj */
                } else {
                    tangB[0]  = normBi[1] * normBj[2] - normBj[1] * normBi[2];
                    tangB[1]  = normBi[2] * normBj[0] - normBj[2] * normBi[0];
                    tangB[2]  = normBi[0] * normBj[1] - normBj[0] * normBi[1];
                    det       = sqrt(tangB[0]*tangB[0] + tangB[1]*tangB[1] + tangB[2]*tangB[2]);
                    tangB[0] /= det;
                    tangB[1] /= det;
                    tangB[2] /= det;
                }
            }

#ifdef DEBUG
            printf("   pntB  =%10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
            printf("   normBi=%10.5f %10.5f %10.5f\n", normBi[0], normBi[1], normBi[2]);
            printf("   normBj=%10.5f %10.5f %10.5f\n", normBj[0], normBj[1], normBj[2]);
            printf("   tangB =%10.5f %10.5f %10.5f\n", tangB[ 0], tangB[ 1], tangB[ 2]);
#endif

            /* set up en conditions (arbirary --- not used) */
            dCdt[ij3  ] = 1;
            dCdt[ij3+1] = 0;
            dCdt[ij3+2] = 0;

            /* both A-side and B-side have corners */
            if ((tangA[0] != 0 || tangA[1] != 0 || tangA[2] != 0) &&
                (tangB[0] != 0 || tangB[1] != 0 || tangB[2] != 0)   ) {
                /* average of A nd B projections (see below) */
                t1 = (normBi[0] * (pntB[ij3  ]-pntA[ij3  ]) + normBi[1] * (pntB[ij3+1]-pntA[ij3+1]) + normBi[2] * (pntB[ij3+2]-pntA[ij3+2]))
                   / (normBi[0] * (           tangA[    0]) + normBi[1] * (           tangA[    1]) + normBi[2] * (           tangA[    2]));
                t2 = (normAi[0] * (pntA[ij3  ]-pntB[ij3  ]) + normAi[1] * (pntA[ij3+1]-pntB[ij3+1]) + normAi[2] * (pntA[ij3+2]-pntB[ij3+2]))
                   / (normAi[0] * (           tangB[    0]) + normAi[1] * (           tangB[    1]) + normAi[2] * (           tangB[    2]));

                pntC[ij3  ] = (pntA[ij3  ] + t1 * tangA[0] + pntB[ij3  ] + t2 * tangB[0]) / 2;
                pntC[ij3+1] = (pntA[ij3+1] + t1 * tangA[1] + pntB[ij3+1] + t2 * tangB[1]) / 2;
                pntC[ij3+2] = (pntA[ij3+2] + t1 * tangA[2] + pntB[ij3+2] + t2 * tangB[2]) / 2;

            /* A-side has a corner */
            } else if (tangA[0] != 0 || tangA[1] != 0 || tangA[2] != 0) {
                /* projection of pntA (in dirn tangA) onto plane B */
                t1 = (normBi[0] * (pntB[ij3  ]-pntA[ij3  ]) + normBi[1] * (pntB[ij3+1]-pntA[ij3+1]) + normBi[2] * (pntB[ij3+2]-pntA[ij3+2]))
                   / (normBi[0] * (           tangA[    0]) + normBi[1] * (           tangA[    1]) + normBi[2] * (           tangA[    2]));

                pntC[ij3  ] = pntA[ij3  ] + t1 * tangA[0];
                pntC[ij3+1] = pntA[ij3+1] + t1 * tangA[1];
                pntC[ij3+2] = pntA[ij3+2] + t1 * tangA[2];

            /* B-side has a corner */
            } else if (tangB[0] != 0 || tangB[1] != 0 || tangB[2] != 0) {
                /* projection of pntB (in dirn tangB) onto plane A */
                t2 = (normAi[0] * (pntA[ij3  ]-pntB[ij3  ]) + normAi[1] * (pntA[ij3+1]-pntB[ij3+1]) + normAi[2] * (pntA[ij3+2]-pntB[ij3+2]))
                   / (normAi[0] * (           tangB[    0]) + normAi[1] * (           tangB[    1]) + normAi[2] * (           tangB[    2]));

                pntC[ij3  ] = pntB[ij3  ] + t2 * tangB[0];
                pntC[ij3+1] = pntB[ij3+1] + t2 * tangB[1];
                pntC[ij3+2] = pntB[ij3+2] + t2 * tangB[2];

            /* neither side has a corner */
            } else {
                /* pntC is on plane intersection; dCdt is tangent along intersection */
                if        (fabs(det = normAi[1]*normBi[2]-normAi[2]*normBi[1]) > EPS03) {
#ifdef DEBUG
                    printf("   case X, det=%12.4e\n", det);
#endif
                    pntC[ij3  ] = (pntA[ij3  ] + pntB[ij3  ]) / 2;
                    rhsA        = normAi[0] * (pntA[ij3  ] - pntC[ij3  ]) + normAi[1] * pntA[ij3+1] + normAi[2] * pntA[ij3+2];
                    rhsB        = normBi[0] * (pntB[ij3  ] - pntC[ij3  ]) + normBi[1] * pntB[ij3+1] + normBi[2] * pntB[ij3+2];
                    pntC[ij3+1] = (rhsA     * normBi[2] - normAi[2] * rhsB    ) / det;
                    pntC[ij3+2] = (normAi[1] * rhsB     - rhsA     * normBi[1]) / det;

                    dCdt[ij3  ] = 1;
                    rhsA       -= normAi[0] + normAi[1] * pntC[ij3+1] + normAi[2] * pntC[ij3+2];
                    rhsB       -= normBi[0] + normBi[1] * pntC[ij3+1] + normBi[2] * pntC[ij3+2];
                    dCdt[ij3+1] = (rhsA     * normBi[2] - normAi[2] * rhsB    ) / det;
                    dCdt[ij3+2] = (normAi[1] * rhsB     - rhsA     * normBi[1]) / det;

                } else if (fabs(det = normAi[2]*normBi[0]-normAi[0]*normBi[2]) > EPS03) {
#ifdef DEBUG
                    printf("   case Y, det=%12.4e\n", det);
#endif
                    pntC[ij3+1] = (pntA[ij3+1] + pntB[ij3+1]) / 2;
                    rhsA        = normAi[1] * (pntA[ij3+1] - pntC[ij3+1]) + normAi[2] * pntA[ij3+2] + normAi[0] * pntA[ij3  ];
                    rhsB        = normBi[1] * (pntB[ij3+1] - pntC[ij3+1]) + normBi[2] * pntB[ij3+2] + normBi[0] * pntB[ij3  ];
                    pntC[ij3+2] = (rhsA     * normBi[0] - normAi[0] * rhsB    ) / det;
                    pntC[ij3  ] = (normAi[2] * rhsB     - rhsA     * normBi[2]) / det;

                    dCdt[ij3+1] = 1;
                    rhsA       -= normAi[1] + normAi[2] * pntC[ij3+2] + normAi[0] * pntC[ij3  ];
                    rhsB       -= normBi[1] + normBi[2] * pntC[ij3+2] + normBi[0] * pntC[ij3  ];
                    dCdt[ij3+2] = (rhsA     * normBi[0] - normAi[0] * rhsB    ) / det;
                    dCdt[ij3  ] = (normAi[2] * rhsB     - rhsA     * normBi[2]) / det;

                } else if (fabs(det = normAi[0]*normBi[1]-normAi[1]*normBi[0]) > EPS03) {
#ifdef DEBUG
                    printf("   case Z, det=%12.4e\n", det);
#endif
                    pntC[ij3+2] = (pntA[ij3+2] + pntB[ij3+2]) / 2;
                    rhsA        = normAi[2] * (pntA[ij3+2] - pntC[ij3+2]) + normAi[0] * pntA[ij3  ] + normAi[1] * pntA[ij3+1];
                    rhsB        = normBi[2] * (pntB[ij3+2] - pntC[ij3+2]) + normBi[0] * pntB[ij3  ] + normBi[1] * pntB[ij3+1];
                    pntC[ij3  ] = (rhsA     * normBi[1] - normAi[1] * rhsB    ) / det;
                    pntC[ij3+1] = (normAi[0] * rhsB     - rhsA     * normBi[0]) / det;

                    dCdt[ij3+2] = 1;
                    rhsA       -= normAi[2] + normAi[0] * pntC[ij3  ] + normAi[1] * pntC[ij3+1];
                    rhsB       -= normBi[2] + normBi[0] * pntC[ij3  ] + normBi[1] * pntC[ij3+1];
                    dCdt[ij3  ] = (rhsA     * normBi[1] - normAi[1] * rhsB    ) / det;
                    dCdt[ij3+1] = (normAi[0] * rhsB     - rhsA     * normBi[0]) / det;

                } else {
                    printf("cannot find plane-plane intersection\n");
                    exit(0);
                }

                /* adjust pntC so that it is at the average of the shortest distance
                   to both pntA and pntB */
                t1 = (dCdt[ij3] * (pntA[ij3]-pntC[ij3]) + dCdt[ij3+1] * (pntA[ij3+1]-pntC[ij3+1]) + dCdt[ij3+2] * (pntA[ij3+2]-pntC[ij3+2]))
                   / (dCdt[ij3] * (          dCdt[ij3]) + dCdt[ij3+1] * (            dCdt[ij3+1]) + dCdt[ij3+2] * (            dCdt[ij3+2]));
                t2 = (dCdt[ij3] * (pntB[ij3]-pntC[ij3]) + dCdt[ij3+1] * (pntB[ij3+1]-pntC[ij3+1]) + dCdt[ij3+2] * (pntB[ij3+2]-pntC[ij3+2]))
                   / (dCdt[ij3] * (          dCdt[ij3]) + dCdt[ij3+1] * (            dCdt[ij3+1]) + dCdt[ij3+2] * (            dCdt[ij3+2]));

                pntC[ij3  ] += (t1 + t2) / 2 * dCdt[ij3  ];
                pntC[ij3+1] += (t1 + t2) / 2 * dCdt[ij3+1];
                pntC[ij3+2] += (t1 + t2) / 2 * dCdt[ij3+2];
            }

#ifdef DEBUG
            printf("   pntC  =%10.5f %10.5f %10.5f\n", pntC[ij3], pntC[ij3+1], pntC[ij3+2]);
            printf("   dCdt  =%10.5f %10.5f %10.5f\n", dCdt[ij3], dCdt[ij3+1], dCdt[ij3+2]);
#endif

            /* copy pntC from end of iedge to beg of jedge */
            pntC[3*(jedge*npnt)  ] = pntC[ij3  ];
            pntC[3*(jedge*npnt)+1] = pntC[ij3+1];
            pntC[3*(jedge*npnt)+2] = pntC[ij3+2];

            dCdt[3*(jedge*npnt)  ] = dCdt[ij3  ];
            dCdt[3*(jedge*npnt)+1] = dCdt[ij3+1];
            dCdt[3*(jedge*npnt)+2] = dCdt[ij3+2];
        }

        /* now set up pntC and dCdt for all the intermediate points on iedge */
        for (iedge = 0; iedge < nedgeA; iedge++) {
            for (ipnt = 1; ipnt < npnt-1; ipnt++) {
                ij3 = 3 * (ipnt + iedge * npnt);

                /* normalized normal at pntA */
                status = EG_getEdgeUV(efacesA[iedge], eedgesA[iedge], 0, tA[ipnt+iedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesA[iedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                normAi[0]  = data[4] * data[8] - data[5] * data[7];
                normAi[1]  = data[5] * data[6] - data[3] * data[8];
                normAi[2]  = data[3] * data[7] - data[4] * data[6];
                det        = sqrt(normAi[0]*normAi[0] + normAi[1]*normAi[1] + normAi[2]*normAi[2]);
                normAi[0] /= det;
                normAi[1] /= det;
                normAi[2] /= det;

                /* normalized normal at pntB */
                status = EG_getEdgeUV(efacesB[iedge], eedgesB[iedge], 0, tB[ipnt+iedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesB[iedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                normBi[0]  = data[4] * data[8] - data[5] * data[7];
                normBi[1]  = data[5] * data[6] - data[3] * data[8];
                normBi[2]  = data[3] * data[7] - data[4] * data[6];
                det        = sqrt(normBi[0]*normBi[0] + normBi[1]*normBi[1] + normBi[2]*normBi[2]);
                normBi[0] /= det;
                normBi[1] /= det;
                normBi[2] /= det;

                /* parametric equation for the intersection of the normal plane at A
                   and normal plane at B */
                if        (fabs(det = normAi[1]*normBi[2]-normAi[2]*normBi[1]) > EPS03) {
                    pntC[ij3  ] = (pntA[ij3  ] + pntB[ij3  ]) / 2;
                    rhsA        = normAi[0] * (pntA[ij3  ] - pntC[ij3  ]) + normAi[1] * pntA[ij3+1] + normAi[2] * pntA[ij3+2];
                    rhsB        = normBi[0] * (pntB[ij3  ] - pntC[ij3  ]) + normBi[1] * pntB[ij3+1] + normBi[2] * pntB[ij3+2];
                    pntC[ij3+1] = (rhsA     * normBi[2] - normAi[2] * rhsB    ) / det;
                    pntC[ij3+2] = (normAi[1] * rhsB     - rhsA     * normBi[1]) / det;

                    dCdt[ij3  ] = 1;
                    rhsA       -= normAi[0] + normAi[1] * pntC[ij3+1] + normAi[2] * pntC[ij3+2];
                    rhsB       -= normBi[0] + normBi[1] * pntC[ij3+1] + normBi[2] * pntC[ij3+2];
                    dCdt[ij3+1] = (rhsA     * normBi[2] - normAi[2] * rhsB    ) / det;
                    dCdt[ij3+2] = (normAi[1] * rhsB     - rhsA     * normBi[1]) / det;

                } else if (fabs(det = normAi[2]*normBi[0]-normAi[0]*normBi[2]) > EPS03) {
                    pntC[ij3+1] = (pntA[ij3+1] + pntB[ij3+1]) / 2;
                    rhsA        = normAi[1] * (pntA[ij3+1] - pntC[ij3+1]) + normAi[2] * pntA[ij3+2] + normAi[0] * pntA[ij3  ];
                    rhsB        = normBi[1] * (pntB[ij3+1] - pntC[ij3+1]) + normBi[2] * pntB[ij3+2] + normBi[0] * pntB[ij3  ];
                    pntC[ij3+2] = (rhsA     * normBi[0] - normAi[0] * rhsB    ) / det;
                    pntC[ij3  ] = (normAi[2] * rhsB     - rhsA     * normBi[2]) / det;

                    dCdt[ij3+1] = 1;
                    rhsA       -= normAi[1] + normAi[2] * pntC[ij3+2] + normAi[0] * pntC[ij3  ];
                    rhsB       -= normBi[1] + normBi[2] * pntC[ij3+2] + normBi[0] * pntC[ij3  ];
                    dCdt[ij3+2] = (rhsA     * normBi[0] - normAi[0] * rhsB    ) / det;
                    dCdt[ij3  ] = (normAi[2] * rhsB     - rhsA     * normBi[2]) / det;

                } else if (fabs(det = normAi[0]*normBi[1]-normAi[1]*normBi[0]) > EPS03) {
                    pntC[ij3+2] = (pntA[ij3+2] + pntB[ij3+2]) / 2;
                    rhsA        = normAi[2] * (pntA[ij3+2] - pntC[ij3+2]) + normAi[0] * pntA[ij3  ] + normAi[1] * pntA[ij3+1];
                    rhsB        = normBi[2] * (pntB[ij3+2] - pntC[ij3+2]) + normBi[0] * pntB[ij3  ] + normBi[1] * pntB[ij3+1];
                    pntC[ij3  ] = (rhsA     * normBi[1] - normAi[1] * rhsB    ) / det;
                    pntC[ij3+1] = (normAi[0] * rhsB     - rhsA     * normBi[0]) / det;

                    dCdt[ij3+2] = 1;
                    rhsA       -= normAi[2] + normAi[0] * pntC[ij3  ] + normAi[1] * pntC[ij3+1];
                    rhsB       -= normBi[2] + normBi[0] * pntC[ij3  ] + normBi[1] * pntC[ij3+1];
                    dCdt[ij3  ] = (rhsA     * normBi[1] - normAi[1] * rhsB    ) / det;
                    dCdt[ij3+1] = (normAi[0] * rhsB     - rhsA     * normBi[0]) / det;

                } else {
                    printf("cannot find plane-plane intersection\n");
                    exit(0);
                }

                /* adjust pntC so that it is at the average of the shortest distance
                   to both pntA and pntB */
                t1 = (dCdt[ij3] * (pntA[ij3]-pntC[ij3]) + dCdt[ij3+1] * (pntA[ij3+1]-pntC[ij3+1]) + dCdt[ij3+2] * (pntA[ij3+2]-pntC[ij3+2]))
                   / (dCdt[ij3] * (          dCdt[ij3]) + dCdt[ij3+1] * (            dCdt[ij3+1]) + dCdt[ij3+2] * (            dCdt[ij3+2]));
                t2 = (dCdt[ij3] * (pntB[ij3]-pntC[ij3]) + dCdt[ij3+1] * (pntB[ij3+1]-pntC[ij3+1]) + dCdt[ij3+2] * (pntB[ij3+2]-pntC[ij3+2]))
                   / (dCdt[ij3] * (          dCdt[ij3]) + dCdt[ij3+1] * (            dCdt[ij3+1]) + dCdt[ij3+2] * (            dCdt[ij3+2]));

                pntC[ij3  ] += (t1 + t2) / 2 * dCdt[ij3  ];
                pntC[ij3+1] += (t1 + t2) / 2 * dCdt[ij3+1];
                pntC[ij3+2] += (t1 + t2) / 2 * dCdt[ij3+2];
            }

            /* find tD such that pntD are evenly spaced */
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                tD[ipnt] = 0;
            }

            for (ismth = 0; ismth < nsmth; ismth++) {
                dtmax = 0;

                /* update pntD */
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    ij3 = 3 * (ipnt + iedge * npnt);

                    pntD[ij3  ] = pntC[ij3  ] + tD[ipnt] * dCdt[ij3  ];
                    pntD[ij3+1] = pntC[ij3+1] + tD[ipnt] * dCdt[ij3+1];
                    pntD[ij3+2] = pntC[ij3+2] + tD[ipnt] * dCdt[ij3+2];
                }

                /* adjust t such the projections before and after i (in the dCdt direction)
                   are the same */
                for (ipnt = 1; ipnt < npnt-1; ipnt++) {
                    ij3 = 3 * (ipnt + iedge * npnt);

                    dt = (dCdt[ij3  ] * (pntD[ij3-3] - 2*pntD[ij3  ] + pntD[ij3+3])
                        + dCdt[ij3+1] * (pntD[ij3-2] - 2*pntD[ij3+1] + pntD[ij3+4])
                        + dCdt[ij3+2] * (pntD[ij3-1] - 2*pntD[ij3+2] + pntD[ij3+5]))
                       / 2 / (dCdt[ij3]*dCdt[ij3] + dCdt[ij3+1]*dCdt[ij3+1] + dCdt[ij3+2]*dCdt[ij3+2]);
                    if (fabs(dt) > dtmax) {
                        dtmax = fabs(dt);
                    }

                    tD[ipnt] += 0.5 * dt;
                }

                if (dtmax < EPS06) break;
            }

            /* slopes point to pntD */
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ij3 = 3 * (ipnt + iedge * npnt);

                slpA[ij3  ] = pntD[ij3  ] - pntA[ij3  ];
                slpA[ij3+1] = pntD[ij3+1] - pntA[ij3+1];
                slpA[ij3+2] = pntD[ij3+2] - pntA[ij3+2];

                slpB[ij3  ] = pntD[ij3  ] - pntB[ij3  ];
                slpB[ij3+1] = pntD[ij3+1] - pntB[ij3+1];
                slpB[ij3+2] = pntD[ij3+2] - pntB[ij3+2];
            }
        }

    /* compute slpA and slpB by subtracting the local normal components of
       the lines that connect pntA and pntB */
    } else if (METHOD(0) == 3) {
        double  norm[4], dxyz[3], A, B, C, D, E, F, factA, factB, proj;
        double  dslpi[3], dslpj[3];

        for (iedge = 0; iedge < nedgeA; iedge++) {
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ij3 = 3 * (ipnt + iedge * npnt);

                /* line between pntA and pntB */
                dxyz[0] = pntB[ij3  ] - pntA[ij3  ];
                dxyz[1] = pntB[ij3+1] - pntA[ij3+1];
                dxyz[2] = pntB[ij3+2] - pntA[ij3+2];

                /* subtract out normal for faceA */
                status = EG_getEdgeUV(efacesA[iedge], eedgesA[iedge], 0, tA[ipnt+iedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesA[iedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                norm[0] = data[4] * data[8] - data[5] * data[7];
                norm[1] = data[5] * data[6] - data[3] * data[8];
                norm[2] = data[3] * data[7] - data[4] * data[6];
                norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

                proj = dxyz[0] * norm[0] + dxyz[1] * norm[1] + dxyz[2] * norm[2];

                slpA[ij3  ] = dxyz[0] - proj * norm[0] / norm[3];
                slpA[ij3+1] = dxyz[1] - proj * norm[1] / norm[3];
                slpA[ij3+2] = dxyz[2] - proj * norm[2] / norm[3];

                /* subtract out normal for faceB */
                status = EG_getEdgeUV(efacesB[iedge], eedgesB[iedge], 0, tB[ipnt+iedge*npnt], uvFace);
                CHECK_STATUS(EG_getEdgeUV);

                status = EG_evaluate(efacesB[iedge], uvFace, data);
                CHECK_STATUS(EG_evaluate);

                norm[0] = data[4] * data[8] - data[5] * data[7];
                norm[1] = data[5] * data[6] - data[3] * data[8];
                norm[2] = data[3] * data[7] - data[4] * data[6];
                norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

                proj = dxyz[0] * norm[0] + dxyz[1] * norm[1] + dxyz[2] * norm[2];

                slpB[ij3  ] = -(dxyz[0] - proj * norm[0] / norm[3]);
                slpB[ij3+1] = -(dxyz[1] - proj * norm[1] / norm[3]);
                slpB[ij3+2] = -(dxyz[2] - proj * norm[2] / norm[3]);

                /* find factA and factB such that the distance between (pntA+factA*slpA) and
                   (pntB+factB*slpb) is minimized */
                A =   slpA[ij3  ]*slpA[ij3  ] + slpA[ij3+1]*slpA[ij3+1] + slpA[ij3+2]*slpA[ij3+2];
                B = - slpA[ij3  ]*slpB[ij3  ] - slpA[ij3+1]*slpB[ij3+1] - slpA[ij3+2]*slpB[ij3+2];
                C = - slpB[ij3  ]*slpA[ij3  ] - slpB[ij3+1]*slpA[ij3+1] - slpB[ij3+2]*slpA[ij3+2];
                D =   slpB[ij3  ]*slpB[ij3  ] + slpB[ij3+1]*slpB[ij3+1] + slpB[ij3+2]*slpB[ij3+2];

                det = A * D - B * C;

                if (fabs(det) > EPS20) {
                    E = + slpA[ij3  ] * (pntB[ij3  ] - pntA[ij3  ])
                        + slpA[ij3+1] * (pntB[ij3+1] - pntA[ij3+1])
                        + slpA[ij3+2] * (pntB[ij3+2] - pntA[ij3+2]);
                    F = + slpB[ij3  ] * (pntA[ij3  ] - pntB[ij3  ])
                        + slpB[ij3+1] * (pntA[ij3+1] - pntB[ij3+1])
                        + slpB[ij3+2] * (pntA[ij3+2] - pntB[ij3+2]);

                    factA = (E * D - B * F) / det;
                    factB = (A * F - E * C) / det;

                    /* modify the slopes by factA and factB */
                    slpA[ij3  ] *= factA;
                    slpA[ij3+1] *= factA;
                    slpA[ij3+2] *= factA;

                    slpB[ij3  ] *= factB;
                    slpB[ij3+1] *= factB;
                    slpB[ij3+2] *= factB;
                } else {
                    printf("WARNING in udfFlend: slopes could not be modified for (iedge=%d, ipnt=%d)\n", iedge, ipnt);
                }
            }
        }

        /* if the slopes at the end of one segment and the beginning of the next disagree,
           apply a correction that get extended into the adjoining Edge */
        for (iedge = 0; iedge < nedgeA; iedge++) {
            jedge = (iedge + 1) % nedgeA;

            /* remove normal component associated with beg of jedge from iedge's last slope */
            status = EG_getEdgeUV(efacesA[jedge], eedgesA[jedge], 0, tA[0+jedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesA[jedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            norm[0] = data[4] * data[8] - data[5] * data[7];
            norm[1] = data[5] * data[6] - data[3] * data[8];
            norm[2] = data[3] * data[7] - data[4] * data[6];
            norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

            proj = slpA[3*(npnt-1+iedge*npnt)  ] * norm[0]
                 + slpA[3*(npnt-1+iedge*npnt)+1] * norm[1]
                 + slpA[3*(npnt-1+iedge*npnt)+2] * norm[2];

            newSlope[0] = slpA[3*(npnt-1+iedge*npnt)  ] - proj * norm[0] / norm[3];
            newSlope[1] = slpA[3*(npnt-1+iedge*npnt)+1] - proj * norm[1] / norm[3];
            newSlope[2] = slpA[3*(npnt-1+iedge*npnt)+2] - proj * norm[2] / norm[3];

            /* remove normal component assoaiated with end of iedge from jedge's first slope */
            status = EG_getEdgeUV(efacesA[iedge], eedgesA[iedge], 0, tA[npnt-1+iedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesA[iedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            norm[0] = data[4] * data[8] - data[5] * data[7];
            norm[1] = data[5] * data[6] - data[3] * data[8];
            norm[2] = data[3] * data[7] - data[4] * data[6];
            norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

            proj = slpA[3*(0+jedge*npnt)  ] * norm[0]
                 + slpA[3*(0+jedge*npnt)+1] * norm[1]
                 + slpA[3*(0+jedge*npnt)+2] * norm[2];

            newSlope[0] = (newSlope[0] + slpA[3*(0+jedge*npnt)  ] - proj * norm[0] / norm[3]) / 2;
            newSlope[1] = (newSlope[1] + slpA[3*(0+jedge*npnt)+1] - proj * norm[1] / norm[3]) / 2;
            newSlope[2] = (newSlope[2] + slpA[3*(0+jedge*npnt)+2] - proj * norm[2] / norm[3]) / 2;

            dslpi[0] = newSlope[0] - slpA[3*(npnt-1+iedge*npnt)  ];
            dslpi[1] = newSlope[1] - slpA[3*(npnt-1+iedge*npnt)+1];
            dslpi[2] = newSlope[2] - slpA[3*(npnt-1+iedge*npnt)+2];

            dslpj[0] = newSlope[0] - slpA[3*(0     +jedge*npnt)  ];
            dslpj[1] = newSlope[1] - slpA[3*(0     +jedge*npnt)+1];
            dslpj[2] = newSlope[2] - slpA[3*(0     +jedge*npnt)+2];

            /* adjust end of iedge and beginning of jedge */
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = 1 - (double)(ipnt) / (double)(npnt-1);

                slpA[3*(npnt-1-ipnt+iedge*npnt)  ] += fraci * dslpi[0];
                slpA[3*(npnt-1-ipnt+iedge*npnt)+1] += fraci * dslpi[1];
                slpA[3*(npnt-1-ipnt+iedge*npnt)+2] += fraci * dslpi[2];

                slpA[3*(       ipnt+jedge*npnt)  ] += fraci * dslpj[0];
                slpA[3*(       ipnt+jedge*npnt)+1] += fraci * dslpj[1];
                slpA[3*(       ipnt+jedge*npnt)+2] += fraci * dslpj[2];
            }
        }

        for (iedge = 0; iedge < nedgeB; iedge++) {
            jedge = (iedge + 1) % nedgeB;

            /* remove normal component associated with beg of jedge from iedge's last slope */
            status = EG_getEdgeUV(efacesB[jedge], eedgesB[jedge], 0, tB[0+jedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesB[jedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            norm[0] = data[4] * data[8] - data[5] * data[7];
            norm[1] = data[5] * data[6] - data[3] * data[8];
            norm[2] = data[3] * data[7] - data[4] * data[6];
            norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

            proj = slpB[3*(npnt-1+iedge*npnt)  ] * norm[0]
                 + slpB[3*(npnt-1+iedge*npnt)+1] * norm[1]
                 + slpB[3*(npnt-1+iedge*npnt)+2] * norm[2];

            newSlope[0] = slpB[3*(npnt-1+iedge*npnt)  ] - proj * norm[0] / norm[3];
            newSlope[1] = slpB[3*(npnt-1+iedge*npnt)+1] - proj * norm[1] / norm[3];
            newSlope[2] = slpB[3*(npnt-1+iedge*npnt)+2] - proj * norm[2] / norm[3];

            /* remove normal component assoaiated with end of iedge from jedge's first slope */
            status = EG_getEdgeUV(efacesB[iedge], eedgesB[iedge], 0, tB[npnt-1+iedge*npnt], uvFace);
            CHECK_STATUS(EG_getEdgeUV);

            status = EG_evaluate(efacesB[iedge], uvFace, data);
            CHECK_STATUS(EG_evaluate);

            norm[0] = data[4] * data[8] - data[5] * data[7];
            norm[1] = data[5] * data[6] - data[3] * data[8];
            norm[2] = data[3] * data[7] - data[4] * data[6];
            norm[3] = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

            proj = slpB[3*(0+jedge*npnt)  ] * norm[0]
                 + slpB[3*(0+jedge*npnt)+1] * norm[1]
                 + slpB[3*(0+jedge*npnt)+2] * norm[2];

            newSlope[0] = (newSlope[0] + slpB[3*(0+jedge*npnt)  ] - proj * norm[0] / norm[3]) / 2;
            newSlope[1] = (newSlope[1] + slpB[3*(0+jedge*npnt)+1] - proj * norm[1] / norm[3]) / 2;
            newSlope[2] = (newSlope[2] + slpB[3*(0+jedge*npnt)+2] - proj * norm[2] / norm[3]) / 2;

            dslpi[0] = newSlope[0] - slpB[3*(npnt-1+iedge*npnt)  ];
            dslpi[1] = newSlope[1] - slpB[3*(npnt-1+iedge*npnt)+1];
            dslpi[2] = newSlope[2] - slpB[3*(npnt-1+iedge*npnt)+2];

            dslpj[0] = newSlope[0] - slpB[3*(0     +jedge*npnt)  ];
            dslpj[1] = newSlope[1] - slpB[3*(0     +jedge*npnt)+1];
            dslpj[2] = newSlope[2] - slpB[3*(0     +jedge*npnt)+2];

            /* adjust end of iedge and beginning of jedge */
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                fraci = 1 - (double)(ipnt) / (double)(npnt-1);

                slpB[3*(npnt-1-ipnt+iedge*npnt)  ] += fraci * dslpi[0];
                slpB[3*(npnt-1-ipnt+iedge*npnt)+1] += fraci * dslpi[1];
                slpB[3*(npnt-1-ipnt+iedge*npnt)+2] += fraci * dslpi[2];

                slpB[3*(       ipnt+jedge*npnt)  ] += fraci * dslpj[0];
                slpB[3*(       ipnt+jedge*npnt)+1] += fraci * dslpj[1];
                slpB[3*(       ipnt+jedge*npnt)+2] += fraci * dslpj[2];
            }
        }

    } else {
        printf("bad method=%d\n", METHOD(0));
        exit(0);
    }

    /* plot A-B connections and slopes */
    if (PLOT(0) == 1) {
        FILE *fp;

        fp = fopen("flend.plot", "w");
        if (fp != NULL) {

            for (iedge = 0; iedge < nedgeA; iedge++) {
                fprintf(fp, "%5d %5d flend_AB_%d|y\n", npnt, -1, iedge);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    fprintf(fp, " %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
                            pntA[3*(ipnt+iedge*npnt)  ],
                            pntA[3*(ipnt+iedge*npnt)+1],
                            pntA[3*(ipnt+iedge*npnt)+2],
                            pntB[3*(ipnt+iedge*npnt)  ],
                            pntB[3*(ipnt+iedge*npnt)+1],
                            pntB[3*(ipnt+iedge*npnt)+2]);
                }
            }

            for (iedge = 0; iedge < nedgeA; iedge++) {
                fprintf(fp, "%5d %5d flend_slpA_%d|m\n", npnt, -1, iedge);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    fprintf(fp, " %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
                            pntA[3*(ipnt+iedge*npnt)  ],
                            pntA[3*(ipnt+iedge*npnt)+1],
                            pntA[3*(ipnt+iedge*npnt)+2],
                            pntA[3*(ipnt+iedge*npnt)  ]+slpA[3*(ipnt+iedge*npnt)  ],
                            pntA[3*(ipnt+iedge*npnt)+1]+slpA[3*(ipnt+iedge*npnt)+1],
                            pntA[3*(ipnt+iedge*npnt)+2]+slpA[3*(ipnt+iedge*npnt)+2]);
                }
            }

            for (iedge = 0; iedge < nedgeB; iedge++) {
                fprintf(fp, "%5d %5d flend_slpB_%d|c\n", npnt, -1, iedge);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    fprintf(fp, " %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
                            pntB[3*(ipnt+iedge*npnt)  ],
                            pntB[3*(ipnt+iedge*npnt)+1],
                            pntB[3*(ipnt+iedge*npnt)+2],
                            pntB[3*(ipnt+iedge*npnt)  ]+slpB[3*(ipnt+iedge*npnt)  ],
                            pntB[3*(ipnt+iedge*npnt)+1]+slpB[3*(ipnt+iedge*npnt)+1],
                            pntB[3*(ipnt+iedge*npnt)+2]+slpB[3*(ipnt+iedge*npnt)+2]);
                }
            }
            
            fprintf(fp, "    0    0 end\n");
            fclose(fp);
        }
    }

    /* get list of Faces in BodyA and BodyB */
    status = EG_getBodyTopos(ebodyA, NULL, FACE, &nfaceA, &efacesA2);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebodyB, NULL, FACE, &nfaceB, &efacesB2);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(efacesA2);
    SPLINT_CHECK_FOR_NULL(efacesB2);

    /* get a list to hold all Faces */
    MALLOC(efacelist, ego, (nedgeA+nfaceA+nfaceB));

    status = ocsmEvalExpr(modl, "@nbody", &value, &dot, temp);
    CHECK_STATUS(ocsmEvalExpr);

    brchattr[0] = -1;                   // fixed in buildPrimitive because _markFaces_ is not set
    bodyattr[0] = value + 1.01;

    MALLOC(spln, double, 6*npnt);
    MALLOC(west, double, 3*npnt);
    MALLOC(east, double, 3*npnt);

    /* make a Face associated with each Edge for the flend */
    for (iedge = 0; iedge < nedgeA; iedge++) {
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            spln[6*ipnt  ] = pntA[3*(ipnt+iedge*npnt)  ];
            spln[6*ipnt+1] = pntA[3*(ipnt+iedge*npnt)+1];
            spln[6*ipnt+2] = pntA[3*(ipnt+iedge*npnt)+2];

            spln[6*ipnt+3] = pntB[3*(ipnt+iedge*npnt)  ];
            spln[6*ipnt+4] = pntB[3*(ipnt+iedge*npnt)+1];
            spln[6*ipnt+5] = pntB[3*(ipnt+iedge*npnt)+2];

            west[3*ipnt  ] = SLOPEA(0) * slpA[3*(ipnt+iedge*npnt)  ];
            west[3*ipnt+1] = SLOPEA(0) * slpA[3*(ipnt+iedge*npnt)+1];
            west[3*ipnt+2] = SLOPEA(0) * slpA[3*(ipnt+iedge*npnt)+2];

            east[3*ipnt  ] = SLOPEB(0) * slpB[3*(ipnt+iedge*npnt)  ];
            east[3*ipnt+1] = SLOPEB(0) * slpB[3*(ipnt+iedge*npnt)+1];
            east[3*ipnt+2] = SLOPEB(0) * slpB[3*(ipnt+iedge*npnt)+2];
        }

        status = EG_spline2dAppx(context, 0, NULL, NULL, NULL, west, east,
                                 NULL, NULL, NULL, NULL, 2, npnt, spln, toler, &esurf);
        CHECK_STATUS(EG_spline2dAppx);

        data[0] = 0;
        data[1] = 1;
        data[2] = 0;
        data[3] = 1;
        status = EG_makeFace(esurf, SFORWARD, data, &efacelist[iedge]);
        CHECK_STATUS(EG_makeFace);

        /* set _brch and _body Attributes on new Face */
        brchattr[1] = iedge + 1;
        bodyattr[1] = iedge + 1;

        status = EG_attributeAdd(efacelist[iedge], "_brch", ATTRINT, 2,
                                 brchattr, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_attributeAdd(efacelist[iedge], "_body", ATTRINT, 2,
                                 bodyattr, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);
    }

    /* add to efacelist the unmarked Faces in BodyA and BodyB */
    nfacelist = nedgeA;

    for (ifaceA = 0; ifaceA < nfaceA; ifaceA++) {
        status = EG_attributeRet(efacesA2[ifaceA], "_flend", &atype, &alen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status != EGADS_SUCCESS || atype != ATTRSTRING || strcmp(tempClist, "remove") != 0) {
            efacelist[nfacelist++] = efacesA2[ifaceA];
        }
    }

    if (ebodyA != ebodyB) {
        for (ifaceB = 0; ifaceB < nfaceB; ifaceB++) {
            status = EG_attributeRet(efacesB2[ifaceB], "_flend", &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != EGADS_SUCCESS || atype != ATTRSTRING || strcmp(tempClist, "remove") != 0) {
                efacelist[nfacelist++] = efacesB2[ifaceB];
            }
        }
    }

    /* sew the Faces into a single (output) Body */
    status = EG_sewFaces(nfacelist, efacelist, TOLER(0), 0, &emodel2);
    CHECK_STATUS(EG_sewFaces);

    status = EG_getTopology(emodel2, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (nchild != 1) {
        snprintf(message, 100, "expecting emodel to have only one child during SEW, but has %d... increase toler", nchild);
        status = OCSM_UDP_ERROR5;
        goto cleanup;
    }

    status = EG_copyObject(echilds[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    SPLINT_CHECK_FOR_NULL(*ebody);

    status = EG_getTopology(*ebody, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (mtype != SOLIDBODY) {
        snprintf(message, 100, "sewing resulted in mtype=%d (not %d)... increase toler", mtype, SOLIDBODY);
        status = OCSM_UDP_ERROR6;
        goto cleanup;
    }

    status = EG_deleteObject(emodel2);
    CHECK_STATUS(EG_deleteObject);

    /* clean up */
    status = EG_deleteObject(context);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(ebodyA);
    CHECK_STATUS(EG_deleteObject);

    if (ebodyB != ebodyA) {
        status = EG_deleteObject(ebodyB);
        CHECK_STATUS(EG_deleteObject);
    }

    /* return output value(s) */

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (efacesA2 != NULL) EG_free(efacesA2);
    if (efacesB2 != NULL) EG_free(efacesB2);

    FREE(efacelist);
    FREE(spln     );
    FREE(west     );
    FREE(east     );

    FREE(efacesA);
    FREE(tA     );
    FREE(pntA   );
    FREE(slpA   );
    FREE(efacesB);
    FREE(tB     );
    FREE(pntB   );
    FREE(slpB   );
    FREE(pntC   );
    FREE(pntD   );
    FREE(dCdt   );
    FREE(tD     );

    if (strlen(message) > 0) {
        *string = message;
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
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int iudp, judp;

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

    /* this routine is not written yet */
    return EGADS_NOLOAD;
}


/*
 ************************************************************************
 *                                                                      *
 *   exposedLoops - make 1 or 2 Loops from exposed Edges                *
 *                                                                      *
 ************************************************************************
 */

static int
exposedLoops(ego    ebody,              /* (in)  Body */
             int    *nloop,             /* (out) number of exposed Loops */
             ego    eloops[])           /* (out) up to 2 exposed Loops */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     nedge, iedge, nlist, ntemp, itemp, count;
    int     oclass, mtype, atype, alen;
    CINT    *tempIlist;
    CDOUBLE *tempRlist;
    CCHAR   *tempClist;
    ego     *eedges, topRef, prev, next, *etemps, edum;
    ego     *elist=NULL;

    ROUTINE(exposedLoops);

    /* --------------------------------------------------------------- */

    *nloop = 0;

    /* make a list of the Edges that bound only one Face (after the Faces
       with _flend=remove have been removed) */
    status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(eedges);

    nlist = 0;
    MALLOC(elist, ego, nedge);

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getInfo(eedges[iedge], &oclass, &mtype, &topRef, &prev, &next);
        CHECK_STATUS(EG_getInfo);

        if (mtype == DEGENERATE) continue;

        status = EG_getBodyTopos(ebody, eedges[iedge], FACE, &ntemp, &etemps);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(etemps);

        /* this Edge has only 1 adjacent Face */
        if (ntemp == 1) {
            count = 1;

        /* this Edge will be non-manifold when Faces with _flend-remove are removed */
        } else {
            count = 0;
            for (itemp = 0; itemp < ntemp; itemp++) {
                status = EG_attributeRet(etemps[itemp], "_flend", &atype, &alen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status == EGADS_SUCCESS) {
                    if (atype == ATTRSTRING && strcmp(tempClist, "remove") == 0) {
                        count++;
                    }
                }
            }
        }

        if (count == 1) {
            elist[nlist++] = eedges[iedge];
        }

        EG_free(etemps);
    }

#ifdef DEBUG
    printf("exposedLoops: nlist=%d\n", nlist);
    for (count = 0; count < nlist; count++) {
        printf("count=%d, iedge=%d\n", count, EG_indexBodyTopo(ebody, eedges[count]));
    }
#endif

    /* combine the Edges into Loops... first Loop */
    status = EG_makeLoop(nlist, elist, NULL, 0, &eloops[0]);
    CHECK_STATUS(EG_makeLoop);

    (*nloop)++;

    /* ...second Loop (if there are any unused Edges) */
    if (status > 0) {
        status = EG_makeLoop(nlist, elist, NULL, 0, &eloops[1]);
        CHECK_STATUS(EG_makeLoop);

        (*nloop)++;

        /* ... subsequent Loops (if there are any unused Edges)
               which are not returned */
        while (status > 0) {
            status = EG_makeLoop(nlist, elist, NULL, 0, &edum);
            CHECK_STATUS(EG_makeLoop);

            (*nloop)++;
        }
    }

    EG_free(eedges);

cleanup:
    FREE(elist);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fillPointsFromEdge - fill in points array associated with given Edge *
 *                                                                      *
 ************************************************************************
 */

static int
fillPointsFromEdge(ego    eedgeA,       /* (in)  Edge in BodyA */
                   int    senseA,       /* (in)  sense of edgeA */
                   ego    eedgeB,       /* (in)  Edge in BodyB */
                   int    senseB,       /* (in)  sense of edgeB */
                   int    npnt,         /* (in)  number of points points */
                   double tA[],         /* (out) t for Points along EdgeA */
                   double pntA[],       /* (out) Points along EdgeA */
                   double tB[],         /* (out) t for Points along EdgeB */
                   double pntB[],       /* (out) Points along EdgeB */
       /*@unused@*/int    *NumUdp,
                   udp_T  *udps)
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     periodic, ipnt;
    double  trangeA[2], trangeB[2], lenA, lenB, data[18];
    double  fraci, stgt, tleft, trite, tt, ss;

    ROUTINE(fillPointsFromEdge);

    /* --------------------------------------------------------------- */

    status = EG_getRange(eedgeA, trangeA, &periodic);
    CHECK_STATUS(EG_getRange);

    status = EG_arcLength(eedgeA, trangeA[0], trangeA[1], &lenA);
    CHECK_STATUS(EG_arcLength);

#ifdef DEBUG
    printf("sensesA=%d, trangeA=%10.5f %10.5f, lenA=%10.5f\n", senseA, trangeA[0], trangeA[1], lenA);
    (void) EG_evaluate(eedgeA, &trangeA[0], data);
    printf("    beg(A) = %10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
    (void) EG_evaluate(eedgeA, &trangeA[1], data);
    printf("    end(A) = %10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
#endif

    status = EG_getRange(eedgeB, trangeB, &periodic);
    CHECK_STATUS(EG_getRange);

    status = EG_arcLength(eedgeB, trangeB[0], trangeB[1], &lenB);
    CHECK_STATUS(EG_arcLength);

#ifdef DEBUG
    printf("sensesB=%d, trangeB=%10.5f %10.5f, lenB=%10.5f\n", senseB, trangeB[0], trangeB[1], lenB);
    (void) EG_evaluate(eedgeB, &trangeB[0], data);
    printf("    beg(B) = %10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
    (void) EG_evaluate(eedgeB, &trangeB[1], data);
    printf("    end(B) = %10.5f %10.5f %10.5f\n", data[0], data[1], data[2]);
#endif

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        fraci = (double)(ipnt) / (double)(npnt-1);

        /* eedgeA: equal arc-length spacing on both Edges */
        if (EQUIS(0) == 1) {
            if (senseA == SFORWARD) {
                stgt =    fraci  * lenA;
            } else {
                stgt = (1-fraci) * lenA;
            }

            tleft = trangeA[0];
            trite = trangeA[1];
            tt    = (tleft + trite) / 2;          // needed for scan-build
            while (trite-tleft > 1.0e-7) {
                tt = (tleft + trite) / 2;

                status = EG_arcLength(eedgeA, trangeA[0], tt, &ss);
                CHECK_STATUS(EG_arcLength);

                if (ss < stgt) {
                    tleft = tt;
                } else {
                    trite = tt;
                }
            }

        /* eedgeA: equal t spacing on Edge A and arc-length matching on Edge B */
        } else if (EQUIS(0) == 2) {
            if (senseA == SFORWARD) {
                tt = (1-fraci) * trangeA[0] + fraci * trangeA[1];
            } else {
                tt = (1-fraci) * trangeA[1] + fraci * trangeA[0];
            }

        /* eedgeA: equal t spacing on Edge B and arc-length matching on Edge A */
        } else if (EQUIS(0) == 3) {
            if (senseB == SFORWARD) {
                tt = (1-fraci) * trangeB[0] + fraci * trangeB[1];
            } else {
                tt = (1-fraci) * trangeB[1] + fraci * trangeB[0];
            }

            status = EG_arcLength(eedgeB, trangeB[0], tt, &ss);
            CHECK_STATUS(EG_arcLength);

            if (senseB == senseA) {
                stgt =         ss  / lenB * lenA;
            } else {
                stgt = (lenB - ss) / lenB * lenA;
            }

            tleft = trangeA[0];
            trite = trangeA[1];
            while (trite-tleft > 1.0e-7) {
                tt = (tleft + trite) / 2;

                status = EG_arcLength(eedgeA, trangeA[0], tt, &ss);
                CHECK_STATUS(EG_arcLength);

                if (ss < stgt) {
                    tleft = tt;
                } else {
                    trite = tt;
                }
            }

        /* eedgeA: equal t spacing on both Edges (the default) */
        } else {
            if (senseA == SFORWARD) {
                tt = (1-fraci) * trangeA[0] + fraci * trangeA[1];
            } else {
                tt = (1-fraci) * trangeA[1] + fraci * trangeA[0];
            }
        }

        status = EG_evaluate(eedgeA, &tt, data);
        CHECK_STATUS(EG_evaluate);

        tA[    ipnt  ] = tt;
        pntA[3*ipnt  ] = data[0];
        pntA[3*ipnt+1] = data[1];
        pntA[3*ipnt+2] = data[2];
#ifdef DEBUG
        printf("A: ipnt=%3d, fraci=%10.5f, tt=%10.5f, data=%10.5f %10.5f %10.5f\n",
               ipnt, fraci, tt, data[0], data[1], data[2]);
#endif

        /* eedgeB: equal arc-length spacing on both Edges */
        if (EQUIS(0) == 1) {
            if (senseB == SFORWARD) {
                stgt =    fraci  * lenB;
            } else {
                stgt = (1-fraci) * lenB;
            }

            tleft = trangeB[0];
            trite = trangeB[1];
            while (trite-tleft > 1.0e-7) {
                tt = (tleft + trite) / 2;

                status = EG_arcLength(eedgeB, trangeB[0], tt, &ss);
                CHECK_STATUS(EG_arcLength);

                if (ss < stgt) {
                    tleft = tt;
                } else {
                    trite = tt;
                }
            }

        /* eedgeB: equal t spacing on Edge A and arc-length matching on Edge B */
        } else if (EQUIS(0) == 2) {
            if (senseA == SFORWARD) {
                tt = (1-fraci) * trangeA[0] + fraci * trangeA[1];
            } else {
                tt = (1-fraci) * trangeA[1] + fraci * trangeA[0];
            }

            status = EG_arcLength(eedgeA, trangeA[0], tt, &ss);
            CHECK_STATUS(EG_arcLength);

            if (senseA == senseB) {
                stgt =         ss  / lenA * lenB;
            } else {
                stgt = (lenA - ss) / lenA * lenB;
            }

            tleft = trangeB[0];
            trite = trangeB[1];
            while (trite-tleft > 1.0e-7) {
                tt = (tleft + trite) / 2;

                status = EG_arcLength(eedgeB, trangeB[0], tt, &ss);
                CHECK_STATUS(EG_arcLength);

                if (ss < stgt) {
                    tleft = tt;
                } else {
                    trite = tt;
                }
            }

        /* eedgeB: equal t spacing on Edge B and arc-length matching on Edge A */
        } else if (EQUIS(0) == 3) {
            if (senseB == SFORWARD) {
                tt = (1-fraci) * trangeB[0] + fraci * trangeB[1];
            } else {
                tt = (1-fraci) * trangeB[1] + fraci * trangeB[0];
            }

        /* eedgeB: equal t spacing on both Edges (the default) */
        } else {
            if (senseB == SFORWARD) {
                tt = (1-fraci) * trangeB[0] + fraci * trangeB[1];
            } else {
                tt = (1-fraci) * trangeB[1] + fraci * trangeB[0];
            }
        }

        status = EG_evaluate(eedgeB, &tt, data);
        CHECK_STATUS(EG_evaluate);

        tB[    ipnt  ] = tt;
        pntB[3*ipnt  ] = data[0];
        pntB[3*ipnt+1] = data[1];
        pntB[3*ipnt+2] = data[2];
#ifdef DEBUG
        printf("B: ipnt=%3d, fraci=%10.5f, tt=%10.5f, data=%10.5f %10.5f %10.5f\n",
               ipnt, fraci, tt, data[0], data[1], data[2]);
#endif
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fillSlopesFromEdge - fill in slopes array associated with given Edge *
 *                                                                      *
 ************************************************************************
 */

static int
fillSlopesFromEdge(ego    eedge,
                   int    sense,
                   ego    eface,
                   int    npnt,
                   double t[],
                   double slp[])
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     ipnt, oclass, mtype;
    double  edata[18], uv[2], fdata[18], vec[4];
    ego     topRef, prev, next;

    ROUTINE(fillSlopesFromEdge);

    /* --------------------------------------------------------------- */

    status = EG_getInfo(eface, &oclass, &mtype, &topRef, &prev, &next);
    CHECK_STATUS(EG_getInfo);

#ifdef DEBUG
    printf("sense=%d, mtype=%d\n", sense, mtype);
#endif

    for (ipnt = 0; ipnt < npnt; ipnt++) {

        /* get the data associated with the Edge */
        status = EG_evaluate(eedge, &t[ipnt], edata);
        CHECK_STATUS(EG_evaluate);

        /* get the data associated with the Face */
        status = EG_getEdgeUV(eface, eedge, 0, t[ipnt], uv);
        CHECK_STATUS(EG_getEdgeUV);

        status = EG_evaluate(eface, uv, fdata);
        CHECK_STATUS(EG_evaluate);

        /* find the cross product of the Face normal and the Edge tangent vector */
        vec[0] = edata[4] * (fdata[3] * fdata[7] - fdata[4] * fdata[6])
               - edata[5] * (fdata[5] * fdata[6] - fdata[3] * fdata[8]);
        vec[1] = edata[5] * (fdata[4] * fdata[8] - fdata[5] * fdata[7])
               - edata[3] * (fdata[3] * fdata[7] - fdata[4] * fdata[6]);
        vec[2] = edata[3] * (fdata[5] * fdata[6] - fdata[3] * fdata[8])
               - edata[4] * (fdata[4] * fdata[8] - fdata[5] * fdata[7]);

        vec[3] = sense * mtype * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        slp[3*ipnt  ] = vec[0] / vec[3];
        slp[3*ipnt+1] = vec[1] / vec[3];
        slp[3*ipnt+2] = vec[2] / vec[3];

#ifdef DEBUG
        if (ipnt == 0 || ipnt == npnt-1) {
            printf("ipnt=%3d, xyz=%10.5f %10.5f %10.5f\n", ipnt, edata[0], edata[1], edata[2]);
            printf("          t  =%10.5f %10.5f %10.5f\n",       edata[3], edata[4], edata[5]);
            printf("          u  =%10.5f %10.5f %10.5f\n",       fdata[3], fdata[4], fdata[5]);
            printf("          v  =%10.5f %10.5f %10.5f\n",       fdata[6], fdata[7], fdata[8]);
            printf("          t*n=%10.5f %10.5f %10.5f\n",       slp[3*ipnt], slp[3*ipnt+1], slp[3*ipnt+2]);
        }
#endif
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   reorderLoop - reorder eloopsB to give minimum twist from eloopA    *
 *                                                                      *
 ************************************************************************
 */

static int
reorderLoop(ego    eloopA,            /* (in)  base Loop */
            ego    *eloopB)           /* (both) modified Loop */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int      iedge, jedge, itest, ishift;
    int      oclassi, oclassj, oclassk, oclassg, mtypei, mtypej, mtypek, mtypeg, nedgei, nedgej, nedgek;
    int      *sensesi, *sensesj, *sensesk, nnode, *sensesnew=NULL;
    double   uvlimitsi[4], uvlimitsj[4], uvlimitsk[4], data[18], *xyzi=NULL, *xyzj=NULL;
    double   ltest, lshift;
    ego      erefi, erefj, erefk, erefg, *eedgesi, *eedgesj, *eedgesk;
    ego      *enodes, *elist, *eedgesnew=NULL, context;

    ROUTINE(reorderLoop);

    /* --------------------------------------------------------------- */

    /* get info on eloopA */
    status = EG_getTopology(eloopA, &erefi,
                            &oclassi, &mtypei, uvlimitsi, &nedgei, &eedgesi, &sensesi);
    CHECK_STATUS(EG_getTopology);

    /* get info on eloopB */
    status = EG_getTopology(*eloopB, &erefj,
                            &oclassj, &mtypej, uvlimitsj, &nedgej, &eedgesj, &sensesj);
    CHECK_STATUS(EG_getTopology);

    /* set up the coordinates at the Nodes and midpoints in the two Loops */
    MALLOC(xyzi, double, 6*nedgei);
    MALLOC(xyzj, double, 6*nedgei);

    for (iedge = 0; iedge < nedgei; iedge++) {
        status = EG_getTopology(eedgesi[iedge], &erefk,
                                &oclassk, &mtypeg, uvlimitsk, &nnode, &enodes, &sensesk);
        CHECK_STATUS(EG_getTopology);

        if (sensesi[iedge] > 0) {
            status = EG_getTopology(enodes[0], &erefk,
                                    &oclassk, &mtypek, &(xyzi[6*iedge]), &nnode, &elist, &sensesk);
            CHECK_STATUS(EG_getTopology);
        } else {
            status = EG_getTopology(enodes[1], &erefk,
                                    &oclassk, &mtypek, &(xyzi[6*iedge]), &nnode, &elist, &sensesk);
            CHECK_STATUS(EG_getTopology);
        }

        if (mtypeg != DEGENERATE) {
            uvlimitsk[2] = (uvlimitsk[0] + uvlimitsk[1]) / 2;

            status = EG_evaluate(eedgesi[iedge], &uvlimitsk[2], data);
            CHECK_STATUS(EG_evaluate);

            xyzi[6*iedge+3] = data[0];
            xyzi[6*iedge+4] = data[1];
            xyzi[6*iedge+5] = data[2];
        } else {
            xyzi[6*iedge+3] = xyzi[3*iedge  ];
            xyzi[6*iedge+4] = xyzi[3*iedge+1];
            xyzi[6*iedge+5] = xyzi[3*iedge+2];
        }
    }

    for (iedge = 0; iedge < nedgej; iedge++) {
        status = EG_getTopology(eedgesj[iedge], &erefk,
                                &oclassk, &mtypeg, uvlimitsk, &nnode, &enodes, &sensesk);
        CHECK_STATUS(EG_getTopology);

        if (sensesj[iedge] > 0) {
            status = EG_getTopology(enodes[0], &erefk,
                                    &oclassk, &mtypek, &(xyzj[6*iedge]), &nnode, &elist, &sensesk);
            CHECK_STATUS(EG_getTopology);
        } else {
            status = EG_getTopology(enodes[1], &erefk,
                                    &oclassk, &mtypek, &(xyzj[6*iedge]), &nnode, &elist, &sensesk);
            CHECK_STATUS(EG_getTopology);
        }

        if (mtypeg != DEGENERATE) {
            uvlimitsk[2] = (uvlimitsk[0] + uvlimitsk[1]) / 2;

            status = EG_evaluate(eedgesj[iedge], &uvlimitsk[2], data);
            CHECK_STATUS(EG_evaluate);

            xyzj[6*iedge+3] = data[0];
            xyzj[6*iedge+4] = data[1];
            xyzj[6*iedge+5] = data[2];
        } else {
            xyzj[6*iedge+3] = xyzj[3*iedge  ];
            xyzj[6*iedge+4] = xyzj[3*iedge+1];
            xyzj[6*iedge+5] = xyzj[3*iedge+2];
        }
    }

    /* find the shift of eloopB that minimizes the distance between the Nodes of eloopA
       and the Nodes of eloopB */
    ishift = 0;
    lshift = HUGEQ;

    /* compute the sum of the distances for each candidate shift */
    for (itest = 0; itest < 2*nedgei; itest+=2) {

        /* simple shifts */
        ltest = 0;
        for (iedge = 0; iedge < 2*nedgei; iedge++) {
            jedge = (iedge + itest) % (2*nedgei);

            ltest += sqrt((xyzi[3*iedge  ]-xyzj[3*jedge  ])*(xyzi[3*iedge  ]-xyzj[3*jedge  ])
                         +(xyzi[3*iedge+1]-xyzj[3*jedge+1])*(xyzi[3*iedge+1]-xyzj[3*jedge+1])
                         +(xyzi[3*iedge+2]-xyzj[3*jedge+2])*(xyzi[3*iedge+2]-xyzj[3*jedge+2]));
        }
        if (ltest < lshift) {
            ishift = itest/2;
            lshift = ltest;
        }

        /* flipped shifts */
        ltest = 0;
        for (iedge = 0; iedge < 2*nedgei; iedge++) {
            jedge = (2 * nedgei - iedge + itest) % (2*nedgei);

            ltest += sqrt((xyzi[3*iedge  ]-xyzj[3*jedge  ])*(xyzi[3*iedge  ]-xyzj[3*jedge  ])
                         +(xyzi[3*iedge+1]-xyzj[3*jedge+1])*(xyzi[3*iedge+1]-xyzj[3*jedge+1])
                         +(xyzi[3*iedge+2]-xyzj[3*jedge+2])*(xyzi[3*iedge+2]-xyzj[3*jedge+2]));
        }
        if (ltest < lshift) {
            ishift = - (itest/2 + 1);
            lshift = ltest;
        }
    }

    /* create the new rotated Loop */
    if (ishift != 0) {
        MALLOC(eedgesnew, ego, 2*nedgej);
        MALLOC(sensesnew, int, 2*nedgej);

        /* shift the Loop by ishift */
        if (ishift > 0) {
            for (iedge = 0; iedge < nedgei; iedge++) {
                jedge = (iedge + ishift) % nedgei;
                eedgesnew[iedge] =  eedgesj[jedge];
                sensesnew[iedge] = +sensesj[jedge];
            }
        } else {
            for (iedge = 0; iedge < nedgei; iedge++) {
                jedge = (2 * nedgei - iedge - ishift - 2) % nedgei;
                eedgesnew[iedge] =  eedgesj[jedge];
                sensesnew[iedge] = -sensesj[jedge];
            }
        }

        status = EG_getTopology(*eloopB, &erefk,
                                &oclassk, &mtypek, uvlimitsk, &nedgek, &eedgesk, &sensesk);
        CHECK_STATUS(EG_getTopology);

        if (oclassk == FACE) {
            status = EG_getGeometry(erefk, &oclassg, &mtypeg, &erefg, NULL, NULL);
            CHECK_STATUS(EG_getGeometry);

            if (mtypeg != PLANE) {
                if (ishift > 0) {
                    for (iedge = 0; iedge < nedgei; iedge++) {
                        jedge = (iedge + ishift) % nedgei;
                        eedgesnew[iedge+nedgei] =  eedgesj[jedge+nedgei];
                        sensesnew[iedge+nedgei] = +sensesj[jedge+nedgei];
                    }
                } else {
                    for (iedge = 0; iedge < nedgei; iedge++) {
                        jedge = (2 * nedgei - iedge - ishift - 2) % nedgei;
                        eedgesnew[iedge+nedgei] =  eedgesj[jedge+nedgei];
                        sensesnew[iedge+nedgei] = -sensesj[jedge+nedgei];
                    }
                }
            }
        }

        /* make new Loop */
        status = EG_getContext(eloopA, &context);
        CHECK_STATUS(EG_getContext);

        status = EG_makeTopology(context, erefk, LOOP, CLOSED,
                                 NULL, nedgej, eedgesnew, sensesnew, eloopB);
        CHECK_STATUS(EG_makeTopology);
    }

cleanup:
    FREE(eedgesnew);
    FREE(sensesnew);

    FREE(xyzi);
    FREE(xyzj);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   slopeAtNode - compute new slope at Node                            *
 *                                                                      *
 ************************************************************************
 */

static int
slopeAtNode(ego    eedge,               /* (in)  Edge starting at Node */
            int    sense,               /* (in)  sense of eedge */
            ego    eface1,              /* (in)  Face associated with eedge */
            ego    eface2,              /* (in)  other Face that is associated with beg Node */
            ego    ebody,               /* (in)  parent Body */
            double newSlope[])          /* (out) new slope at Node */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     periodic, oclass, mtype, nnode, *senses, nedge2, iedge2, nface2;
    double  tang[18], uv[2], trange[2], tt, data[18];
    double  tlen;
    ego     eref, enode, *enodes, eedge2, *eedges2=NULL, *eefaces2=NULL;

    ROUTINE(fixSlopeAtNode);

    /* --------------------------------------------------------------- */

    status = EG_getRange(eedge, trange, &periodic);
    CHECK_STATUS(EG_getRange);

    if (sense == SFORWARD) {
        tt = trange[0];
    } else {
        tt = trange[1];
    }

    status = EG_getEdgeUV(eface1, eedge, 0, tt, uv);
    CHECK_STATUS(EG_getEdgeUV);

    /* find Node that from which Edge will emanate */
    status = EG_getTopology(eedge, &eref, &oclass, &mtype, data, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);

    if (sense == SFORWARD) {
        enode = enodes[0];
    } else {
        enode = enodes[1];
    }

    /* find an Edge which emanates from enode and which is shared by
       both eface1 and eface2 */
    status = EG_getBodyTopos(ebody, enode, EDGE, &nedge2, &eedges2);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(eedges2);

    eedge2 = NULL;
    for (iedge2 = 0; iedge2 < nedge2; iedge2++) {
        status = EG_getBodyTopos(ebody, eedges2[iedge2], FACE, &nface2, &eefaces2);
        CHECK_STATUS(EG_getBodyTopo);

        SPLINT_CHECK_FOR_NULL(eefaces2);

        if (nface2 < 2) {
        } else if (eefaces2[0] == eface1 && eefaces2[1] == eface2) {
            eedge2 = eedges2[iedge2];
            EG_free(eefaces2);
            break;
        } else if (eefaces2[0] == eface2 && eefaces2[1] == eface1) {
            eedge2 = eedges2[iedge2];
            EG_free(eefaces2);
            break;
        }

        EG_free(eefaces2);
    }

    EG_free(eedges2);

    if (eedge2 == NULL) goto cleanup;      // needed for splint

    status = EG_getTopology(eedge2, &eref, &oclass, &mtype, data, &nnode, &enodes, &senses);
    CHECK_STATUS(EG_getTopology);

    if        (enodes[0] == enode) {
        status = EG_evaluate(eedge2, &data[0], tang);
        CHECK_STATUS(EG_evaluate);

        tlen = -sqrt(tang[3]*tang[3] + tang[4]*tang[4] + tang[5]*tang[5]);

        newSlope[0] = tang[3] / tlen;
        newSlope[1] = tang[4] / tlen;
        newSlope[2] = tang[5] / tlen;
    } else if (enodes[1] == enode) {
        status = EG_evaluate(eedge2, &data[1], tang);
        CHECK_STATUS(EG_evaluate);

        tlen = +sqrt(tang[3]*tang[3] + tang[4]*tang[4] + tang[5]*tang[5]);

        newSlope[0] = tang[3] / tlen;
        newSlope[1] = tang[4] / tlen;
        newSlope[2] = tang[5] / tlen;
    } else {
        printf("did not match Nodes\n");
        exit(0);
    }

cleanup:
    return status;
}
