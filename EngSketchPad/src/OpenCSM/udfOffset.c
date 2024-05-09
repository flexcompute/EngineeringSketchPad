 /*
 ************************************************************************
 *                                                                      *
 * udfOffset -- make an offset of a planar WireBody or SheetBody        *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2024  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 5
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NODELIST(IUDP,I)  ((int    *) (udps[IUDP].arg[0].val))[I]
#define NODEDIST(IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define EDGELIST(IUDP,I)  ((int    *) (udps[IUDP].arg[2].val))[I]
#define FACELIST(IUDP,I)  ((int    *) (udps[IUDP].arg[3].val))[I]
#define DIST(    IUDP)    ((double *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"nodelist", "nodedist", "edgelist", "facelist", "dist",   };
static int    argTypes[NUMUDPARGS] = {ATTRINT,    ATTRREAL,   ATTRINT,    ATTRINT,    ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,          0,          0,          0,        };
static double argDdefs[NUMUDPARGS] = {0.,         1.,         0.,         0.,         0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"
#include <assert.h>

#define EPS12     1.0e-12
#define EPS30     1.0e-30

#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) < (B)) ? (B) : (A))

#ifdef GRAFIC
   #include "grafic.h"
#endif

/* structures */
typedef struct {
    double  xyz[3];           /* X,Y,Z coordinates */
    double  dist;             /* distance from Wnode */
    int     nedge;            /* number of incident Wedges */
    ego     enode;            /* ego in brep */
} wnode_T;

typedef struct {
    int     tag;              /* 1=can be offset from */
    int     ibeg;             /* Wnode at beginning */
    int     iend;             /* Wnode at end */
    int     ileft;            /* Wface on left */
    int     irite;            /* Wface on rite */
    int     ibleft;           /* Wedge at beginning associated with ileft */
    int     ibrite;           /* Wedge at beginning associated with irite */
    int     ieleft;           /* Wedge at end       associated with ileft */
    int     ierite;           /* Wedge at end       associated with irite */
    double  tbeg;             /* t at ibeg */
    double  tend;             /* t at end */
    ego     ecurve;           /* curve */
    ego     epleft;           /* pcurve associated with ileft */
    ego     eprite;           /* pcurve associated with irite */
    ego     eedge;            /* ego in brep */
} wedge_T;

typedef struct {
    int     tag;              /* 1=can be split */
    int     new;              /* =1 if new scribed Face */
    int     mtype;            /* SFORWARD or SREVERSE */
    ego     esurface;         /* surface */
    ego     eface;            /* ego in brep */
} wface_T;

typedef struct {
    int     nnode;            /* number of Wnodes */
    int     nedge;            /* number of Wedges */
    int     nface;            /* number of Wfaces */
    int     nface_orig;       /* number of Wfaces after wrepInit */
    wnode_T *node;            /* array  of Wnodes */
    wedge_T *edge;            /* array  of Wedges */
    wface_T *face;            /* array  of Wfaces */
    ego     ebody;            /* ego in brep */
} wrep_T;

typedef struct {
    int     iface;            /* Wface associated with this offset */
    int     ibeg;             /* Wnode at the beginning of this offset */
    int     iend;             /* Wnode at the end       of this offset */
    double  tbeg;             /* t     at the beginning of this offset */
    double  tend;             /* t     at the end       of this offset */
    ego     ecurve;           /* EGADS curve associated with this offset */
} offset_T;

/* prototypes */
static int makePlanarOffset(ego ebodyIn, ego *ebodyOut, int *NumUdp, udp_T *udps);
static int intersectLines(int npnt1, double xyz1[], int npnt2, double xyz2[], int iplane, double *t1, double *t2);

static int makeSolidOffset(ego ebodyIn, ego *ebodyOut, int *NumUdp, udp_T *udps);
static int makeOffsetCurve(wrep_T *wrep, int iface, int iedge, offset_T *offset);
static int wrepInit(ego ebody, wrep_T **wrep);
static int wrepBreakEdge(wrep_T *wrep, int iedge, double t, int *jnode);
static int wrepMakeEdge(wrep_T *wrep, int iface, int ibeg, int iend, /*@null@*/ego einput);
static int wrepMakeNode(wrep_T *wrep, int iface, double xyz[]);
static int wrep2ego(wrep_T wrep, ego context, ego *ebody);
static int wrepPrint(wrep_T *wrep);
static int wrepFree(wrep_T *wrep);

/* message (to be shared by all functions) */
static char message[1024];


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  emodel,                 /* (in)  Model containing Body */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nchild, *senses;
    double  data[18];
    ego     eref, *ebodys, *echilds;
    udp_T   *udps = *Udps;

#ifdef DEBUG
    int     i;
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("nodelist(0) =");
    for (i = 0; i < udps[0].arg[0].size; i++) {
        printf(" %d", NODELIST(0,i));
    }
    printf("\n");
    printf("nodedist(0) =");
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf(" %f", NODEDIST(0,i));
    }
    printf("\n");
    printf("edgelist(0) =");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf(" %d", EDGELIST(0,i));
    }
    printf("\n");
    printf("facelist(0) =");
    for (i = 0; i < udps[0].arg[3].size; i++) {
        printf(" %d", FACELIST(0,i));
    }
    printf("\n");
    printf("dist(    0) = %f\n", DIST(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 1024, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 1024, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* make sure input is WireBody or SheetBody */
    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    /* if planar, make planar offset */
    if (mtype == WIREBODY || mtype == FACEBODY || mtype == SHEETBODY) {
        status = makePlanarOffset(ebodys[0], ebody, NumUdp, udps);
        CHECK_STATUS(makePlanarOffset);

    /* of a solid, make offsets from specified Edges on specified Faces */
    } else if (mtype == SOLIDBODY) {

#ifdef GRAFIC
        {
            int     nface, iface, nloop, iloop, nedge, iedge, i, *senses, periodic;
            double  data[18], trange[4], tt;
            ego     *efaces, *eloops, *eedges, esurf, eref, topRef, prev, next;

            void    *realloc_temp=NULL;

            int     io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
            int     *ilin=NULL, *isym=NULL, *nper=NULL, mline, nline, mdata, ndata;
            float   *udata=NULL, *vdata=NULL;
            char    pltitl[255];

            status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
            CHECK_STATUS(EG_getBodyTopos);

            mline = 100;
            MALLOC(ilin, int, mline);
            MALLOC(isym, int, mline);
            MALLOC(nper, int, mline);

            mdata = 1000;
            MALLOC(udata, float, mdata);
            MALLOC(vdata, float, mdata);

            for (iface = 0; iface < nface; iface++) {
                ndata = 0;
                nline = 0;

                status = EG_getTopology(efaces[iface], &esurf, &oclass, &mtype,
                                        data, &nloop, &eloops, &senses);
                CHECK_STATUS(EG_getTopology);

                status = EG_getInfo(esurf, &oclass, &mtype, &topRef, &prev, &next);
                CHECK_STATUS(EG_getInfo);

                if (mtype == PLANE) continue;

                for (iloop = 0; iloop < nloop; iloop++) {
                    status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype,
                                            data, &nedge, &eedges, &senses);
                    CHECK_STATUS(EG_getTopology);

                    for (iedge = 0; iedge < nedge; iedge++) {
                        status = EG_getInfo(eedges[iedge], &oclass, &mtype, &topRef, &prev, &next);
                        CHECK_STATUS(EG_getInfo);

                        if (mtype == DEGENERATE) continue;

                        status = EG_getRange(eedges[iedge], trange, &periodic);
                        CHECK_STATUS(EG_getRange);

                        for (i = 0; i < 33; i++) {
                            tt = trange[0] + (trange[1] - trange[0]) * (double)(i) / 32.0;

                            status = EG_evaluate(eedges[iedge+nedge], &tt, data);
                            CHECK_STATUS(EG_evaluate);

                            udata[ndata] = data[0];
                            vdata[ndata] = data[1];
                            ndata++;

                            if (ndata >= mdata) {
                                mdata += 100;
                                RALLOC(udata, float, mdata);
                                RALLOC(vdata, float, mdata);
                            }
                        }

                        ilin[nline] = 1 + nline % 7;
                        isym[nline] = 1 + nline % 4;
                        nper[nline] = 33;
                        nline++;

                        if (nline >= mline) {
                            mline += 1000;
                            RALLOC(ilin, int, mline);
                            RALLOC(isym, int, mline);
                            RALLOC(nper, int, mline);
                        }
                    }
                }

                snprintf(pltitl, 254, "~u~v~Face %d", iface+1);
                grinit_(&io_kbd, &io_scr, "udfOffset", strlen("udfOffset"));
                grline_(ilin, isym, &nline, pltitl, &indgr, udata, vdata, nper, strlen(pltitl));
            }

            EG_free(efaces);

            FREE(ilin);
            FREE(isym);
            FREE(nper);
            FREE(udata);
            FREE(vdata);
        }
#endif

        status = makeSolidOffset(ebodys[0], ebody, NumUdp, udps);
        CHECK_STATUS(makeSolidOffset);

    } else {
        snprintf(message, 1024, "not a WireBody, SheetBody, or SolidBody");
        status = EGADS_NODATA;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("nodelist(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf(" %d", NODELIST(numUdp,i));
    }
    printf("\n");
    printf("nodedist(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf(" %f", NODEDIST(numUdp,i));
    }
    printf("\n");
    printf("edgelist(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
        printf(" %d", EDGELIST(numUdp,i));
    }
    printf("\n");
    printf("facelist(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[3].size; i++) {
        printf(" %d", FACELIST(numUdp,i));
    }
    printf("\n");
    printf("dist(    %d) = %f\n", numUdp, DIST(numUdp));
#endif

    /* set the output value */

    /* return the new Body */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (strlen(message) > 0) {
        MALLOC(*string, char, 1024);
        strncpy(*string, message, 1023);
    } else if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
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
 *   makePlanarOffset - build offset Body for a planar WireBody or SheetBody *
 *                                                                      *
 ************************************************************************
 */

static int
makePlanarOffset(ego    ebodyIn,        /* (in)  input Body */
                 ego    *ebodyOut,      /* (out) output Body */
     /*@unused@*/int    *NumUdp,
                 udp_T  *udps)
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nchild, *senses, iplane, nloop, nedge, medge, iedge, jedge, ipnt, nchange, isave;
    int     periodic, sizes[2], ntemp1, ntemp2, *ibeg=NULL, *iend=NULL, *active=NULL, *newSenses=NULL;
    double  data[18], bbox[6], xplane=0, yplane=0, zplane=0, tt, trange[4], alen, fact, frac;
    double  **xyz=NULL, *tbeg=NULL, *tend=NULL;
    ego     context, eref, *echilds, *eedges=NULL, *eloops=NULL;
    ego     *newCurves=NULL, *newNodes=NULL, *newEdges=NULL;
    ego     newLoop, newFace, newShell, topRef, eprev, enext;

    int     npnt_min=5, npnt_max=201;

#ifdef GRAFIC
    float   *xplot=NULL, *yplot=NULL;
    int     io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int     nline, nplot, *ilin=NULL, *isym=NULL, *nper=NULL;
    char    pltitl[80];
#endif

    ROUTINE(makePlanarOffset);

    /* --------------------------------------------------------------- */

    /* make sure DIST is only one non-zero scalar */
    if (udps[0].arg[4].size != 1) {
        snprintf(message, 1024, "\"dist\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (fabs(DIST(0)) < EPS06) {
        snprintf(message, 1024, "\"dist\" must be non-zero");
        status= EGADS_RANGERR;
        goto cleanup;
    }

    /* make sure body only contains one loop */
    status = EG_getBodyTopos(ebodyIn, NULL, LOOP, &nloop, &eloops);
    CHECK_STATUS(EG_getBodyTopos);

    if (nloop != 1) {
        snprintf(message, 1024, "input Body has %d Loops, but was expecting one\n", nloop);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getTopology(ebodyIn, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    /* if a WireBody, make sure it is not non-manifold */
    if (mtype == WIREBODY) {
        status = EG_getBodyTopos(ebodyIn, NULL, EDGE, &ntemp1, NULL);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(eloops);

        status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                                data, &ntemp2, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (ntemp1 != ntemp2) {
            snprintf(message, 1024, "Input WireBody must be manifold\n");
            status = EGADS_NOTBODY;
            goto cleanup;
        }
    }

    /* make sure input is aligned with one of the coordinate planes */
    status = EG_getBoundingBox(ebodyIn, bbox);
    CHECK_STATUS(EG_getBoundingBox);

    if        (fabs(bbox[5]-bbox[2]) < EPS06) {
        iplane = 3;
        zplane = (bbox[5] + bbox[2]) / 2;
    } else if (fabs(bbox[4]-bbox[1]) < EPS06) {
        iplane = 2;
        yplane = (bbox[4] + bbox[1]) / 2;
    } else if (fabs(bbox[3]-bbox[0]) < EPS06) {
        iplane = 1;
        xplane = (bbox[3] + bbox[0]) / 2;
    } else {
        snprintf(message, 1024, "input Body must be aligned with a coordinate plane: dx=%f dy=%f dz=%f\n",
               bbox[3]-bbox[0], bbox[4]-bbox[1], bbox[5]-bbox[2]);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getContext(ebodyIn, &context);
    CHECK_STATUS(EG_getContext);

    /* get the Edges associated with the input Body */
    SPLINT_CHECK_FOR_NULL(eloops);

    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                            data, &nedge, &eedges, &senses);
    CHECK_STATUS(EG_getTopology);

#ifdef DEBUG
    printf("Body has %d Edges\n", nedge);
#endif

    assert(nedge > 0);
    SPLINT_CHECK_FOR_NULL(eedges);

    /* allocate space for new Curves, Nodes, and Edges */
    MALLOC(ibeg,      int,     nedge  );
    MALLOC(iend,      int,     nedge  );
    MALLOC(tbeg,      double,  nedge  );
    MALLOC(tend,      double,  nedge  );
    MALLOC(active,    int,     nedge  );
    MALLOC(xyz,       double*, nedge  );
    MALLOC(newCurves, ego,     nedge  );
    MALLOC(newNodes,  ego,     nedge+1);
    MALLOC(newEdges,  ego,     nedge  );
    MALLOC(newSenses, int,     nedge  );

    for (iedge = 0; iedge < nedge; iedge++) {
        active[iedge] = 1;
        xyz[   iedge] = NULL;
    }

#ifdef GRAFIC
    MALLOC(ilin,  int,   4*nedge         );
    MALLOC(isym,  int,   4*nedge         );
    MALLOC(nper,  int,   4*nedge         );
    MALLOC(xplot, float, 4*nedge*npnt_max);
    MALLOC(yplot, float, 4*nedge*npnt_max);

    nplot = 0;
    nline = 0;
#endif

    /* make offset points for each of the Edges */
    for (iedge = 0; iedge < nedge; iedge++) {

        /* determine the number of offset points */
        status = EG_getRange(eedges[iedge], trange, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_arcLength(eedges[iedge], trange[0], trange[1], &alen);
        CHECK_STATUS(EG_arclength);

        ibeg[iedge] = 0;
        iend[iedge] = 10 * alen / fabs(DIST(0));
        if (iend[iedge] < npnt_min) iend[iedge] = npnt_min;
        if (iend[iedge] > npnt_max) iend[iedge] = npnt_max;

#ifdef DEBUG
        printf("ibeg[%d]=%d,  iend[%d]=%d\n", iedge, ibeg[iedge], iedge, iend[iedge]);
#endif

        /* generate the points */
        MALLOC(xyz[iedge], double, 3*iend[iedge]+3);

        for (ipnt = ibeg[iedge]; ipnt <= iend[iedge]; ipnt++) {
            if (senses[iedge] == SFORWARD) {
                tt = trange[0] + (double)(ipnt-ibeg[iedge])/(double)(iend[iedge]-ibeg[iedge]) * (trange[1]-trange[0]);
            } else {
                tt = trange[1] + (double)(ipnt-ibeg[iedge])/(double)(iend[iedge]-ibeg[iedge]) * (trange[0]-trange[1]);
            }

            status = EG_evaluate(eedges[iedge], &tt, data);
            CHECK_STATUS(EG_evaluate);

#ifdef GRAFIC
            /* original curve */
            if (iplane == 3) {
                xplot[nplot] = data[0];
                yplot[nplot] = data[1];
                nplot++;
            } else if (iplane == 2) {
                xplot[nplot] = data[2];
                yplot[nplot] = data[0];
                nplot++;
            } else {
                xplot[nplot] = data[1];
                yplot[nplot] = data[2];
                nplot++;
            }
#endif

            /* apply the offset */
            if (senses[iedge] == SFORWARD) {
                fact = +DIST(0) / sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);
            } else {
                fact = -DIST(0) / sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);
            }

            if (iplane == 3) {
                xyz[iedge][3*ipnt  ] = data[0] + fact * data[4];
                xyz[iedge][3*ipnt+1] = data[1] - fact * data[3];
                xyz[iedge][3*ipnt+2] = zplane;
            } else if (iplane == 2) {
                xyz[iedge][3*ipnt  ] = data[0] - fact * data[5];
                xyz[iedge][3*ipnt+1] = yplane;
                xyz[iedge][3*ipnt+2] = data[2] + fact * data[3];
            } else {
                xyz[iedge][3*ipnt  ] = xplane;
                xyz[iedge][3*ipnt+1] = data[1] + fact * data[5];
                xyz[iedge][3*ipnt+2] = data[2] - fact * data[4];
            }
        }

#ifdef GRAFIC
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_X;
        nper[nline] = iend[iedge] - ibeg[iedge] + 1;
        nline++;
#endif
    }

#ifdef GRAFIC
    for (iedge = 0; iedge < nedge; iedge++) {
        if (active[iedge] == 0) continue;

        /* offset curve */
        for (ipnt = ibeg[iedge]; ipnt <= iend[iedge]; ipnt++) {
            if (iplane == 3) {
                xplot[nplot] = xyz[iedge][3*ipnt  ];
                yplot[nplot] = xyz[iedge][3*ipnt+1];
                nplot++;
            } else if (iplane == 2) {
                xplot[nplot] = xyz[iedge][3*ipnt+2];
                yplot[nplot] = xyz[iedge][3*ipnt  ];
                nplot++;
            } else {
                xplot[nplot] = xyz[iedge][3*ipnt+1];
                yplot[nplot] = xyz[iedge][3*ipnt+2];
                nplot++;
            }
        }

        ilin[nline] = +GR_DOTTED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = iend[iedge] - ibeg[iedge] + 1;
        nline++;
    }
#endif

#ifdef DEBUG
    printf("before end adjustment\n");
    for (iedge = 0; iedge < nedge; iedge++) {
        printf("iedge=%d\n", iedge);
        for (ipnt = ibeg[iedge]; ipnt <= iend[iedge]; ipnt++) {
            printf("%3d %12.6f %12.6f %12.6f\n", ipnt, xyz[iedge][3*ipnt], xyz[iedge][3*ipnt+1], xyz[iedge][3*ipnt+2]);
        }
    }
#endif
#ifdef GRAFIC
    if (iplane == 3) {
        strcpy(pltitl, "~x~y~Original curves and offsets");
    } else if (iplane == 2) {
        strcpy(pltitl, "~z~x~Original curves and offsets");
    } else {
        strcpy(pltitl, "~y~z~Original curves and offsets");
    }

    grinit_(&io_kbd, &io_scr, "udfOffset", strlen("udfOffset"));
    grline_(ilin, isym, &nline, pltitl, &indgr, xplot, yplot, nper, strlen(pltitl));
#endif

    /* intersect each offset curve with its neighbors in order to
       find tbeg and tend of the intersections */
    status = EG_getInfo(eloops[0], &oclass, &mtype, &topRef, &eprev, &enext);
    CHECK_STATUS(EG_getInfo);

    medge = nedge;
    nchange = nedge - 1;
    while (nchange > 0) {
        nchange = 0;

        /* use whole range of each active Edge */
        for (iedge = 0; iedge < nedge; iedge++) {
            if (active[iedge] == 1) {
                tbeg[iedge] = 0;
                tend[iedge] = iend[iedge];
            }
        }

        /* find intersection for each pair of adjacent (active) Edges */
        for (iedge = 0; iedge < nedge; iedge++) {
            if (mtype == OPEN && iedge == medge-1) break;

            if (active[iedge] == 0) continue;     // ignore iedge if already inactive

            for (jedge = iedge+1; jedge != iedge; jedge++) {
                if (jedge >= nedge) jedge -= nedge;

                if (active[jedge] == 0) continue; // ignore jedge if already inactive

                status = intersectLines(iend[iedge]+1, xyz[iedge], iend[jedge]+1, xyz[jedge], iplane, &(tend[iedge]), &(tbeg[jedge]));
                if (status == OCSM_UDP_ERROR1) {
                    snprintf(message, 1024, "Edges %d and %d are parallel but not colinear\n", iedge+1, jedge+1);
                    goto cleanup;
                }
                CHECK_STATUS(intersectLines);

                break;
            }
        }

        /* inactivate any Edge in which t is not increasing */
        for (iedge = 0; iedge < nedge; iedge++) {
            if (active[iedge] == 1 && tend[iedge] <= tbeg[iedge]) {
                active[iedge] = 0;
                if (iedge == medge-1) medge--;
                nchange++;
            }
#ifdef DEBUG
            printf("    tbeg[%2d]=%12.4e, tend[%2d]=%12.4e,  active=%d\n", iedge, tbeg[iedge], iedge, tend[iedge], active[iedge]);
#endif
        }
    }

    /* reset ibeg and iend for each active Edge and adjust the first and last points
       so that the ends line up */
    for (iedge = 0; iedge < nedge; iedge++) {
        if (active[iedge] == 0) continue;

        if (nedge == 1) break;

        ibeg[iedge] = MAX(0, (int)tbeg[iedge]);
        frac        = tbeg[iedge] - ibeg[iedge];

        xyz[iedge][3*ibeg[iedge]  ] = (1-frac) * xyz[iedge][3*ibeg[iedge]  ] + frac * xyz[iedge][3*ibeg[iedge]+3];
        xyz[iedge][3*ibeg[iedge]+1] = (1-frac) * xyz[iedge][3*ibeg[iedge]+1] + frac * xyz[iedge][3*ibeg[iedge]+4];
        xyz[iedge][3*ibeg[iedge]+2] = (1-frac) * xyz[iedge][3*ibeg[iedge]+2] + frac * xyz[iedge][3*ibeg[iedge]+5];

        if (fabs(xyz[iedge][3*ibeg[iedge]  ]-xyz[iedge][3*ibeg[iedge]+3]) < EPS06 &&
            fabs(xyz[iedge][3*ibeg[iedge]+1]-xyz[iedge][3*ibeg[iedge]+4]) < EPS06 &&
            fabs(xyz[iedge][3*ibeg[iedge]+2]-xyz[iedge][3*ibeg[iedge]+5]) < EPS06   ) {
            ibeg[iedge]++;
        }

        iend[iedge] = MIN(iend[iedge], 1+(int)(tend[iedge]-EPS12));
        frac        = iend[iedge] - tend[iedge];

        xyz[iedge][3*iend[iedge]  ] = (1-frac) * xyz[iedge][3*iend[iedge]  ] + frac * xyz[iedge][3*iend[iedge]-3];
        xyz[iedge][3*iend[iedge]+1] = (1-frac) * xyz[iedge][3*iend[iedge]+1] + frac * xyz[iedge][3*iend[iedge]-2];
        xyz[iedge][3*iend[iedge]+2] = (1-frac) * xyz[iedge][3*iend[iedge]+2] + frac * xyz[iedge][3*iend[iedge]-1];

        if (fabs(xyz[iedge][3*iend[iedge]  ]-xyz[iedge][3*iend[iedge]-3]) < EPS06 &&
            fabs(xyz[iedge][3*iend[iedge]+1]-xyz[iedge][3*iend[iedge]-2]) < EPS06 &&
            fabs(xyz[iedge][3*iend[iedge]+2]-xyz[iedge][3*iend[iedge]-1]) < EPS06   ) {
            iend[iedge]--;
        }
    }

#ifdef DEBUG
    printf("after end adjustment\n");
    for (iedge = 0; iedge < nedge; iedge++) {
        printf("iedge=%d\n", iedge);
        for (ipnt = ibeg[iedge]; ipnt <= iend[iedge]; ipnt++) {
            printf("%3d %12.6f %12.6f %12.6f\n", ipnt, xyz[iedge][3*ipnt], xyz[iedge][3*ipnt+1], xyz[iedge][3*ipnt+2]);
        }
    }
#endif
    /* make the needed Nodes */
    isave = 0;
    jedge = 0;
    for (iedge = 0; iedge < nedge; iedge++) {
        if (active[iedge] == 0) continue;

#ifdef DEBUG
        printf("Node[%d] at (%10.5f %10.5f %10.5f\n", jedge, xyz[iedge][3*ibeg[iedge]], xyz[iedge][3*ibeg[iedge]+1], xyz[iedge][3*ibeg[iedge]+2]);
#endif

        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(xyz[iedge][3*ibeg[iedge]]), 0, NULL, NULL, &(newNodes[jedge]));
        CHECK_STATUS(EG_makeTopology);

        isave = iedge;

        jedge++;
    }

    if (mtype == CLOSED) {
        newNodes[jedge] = newNodes[0];
    } else {
#ifdef DEBUG
        printf("Node[%d] at (%10.5f %10.5f %10.5f\n", jedge, xyz[isave][3*ibeg[isave]], xyz[isave][3*ibeg[isave]+1], xyz[isave][3*ibeg[isave]+2]);
#endif

        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(xyz[isave][3*iend[isave]]), 0, NULL, NULL, &(newNodes[jedge]));
        CHECK_STATUS(EG_makeTopology);
    }

    /* make the Curves and Edges */
    jedge = 0;
    for (iedge = 0; iedge < nedge; iedge++) {
        if (active[iedge] == 0) continue;

        sizes[0] = iend[iedge] - ibeg[iedge] + 1;
        sizes[1] = 0;
        status = EG_approximate(context, 0, EPS06, sizes, &(xyz[iedge][3*ibeg[iedge]]), &newCurves[jedge]);
        CHECK_STATUS(EG_approximate);

        status = EG_getRange(newCurves[jedge], trange, &periodic);
        CHECK_STATUS(EG_getRange);

        newSenses[jedge] = SFORWARD;
        status = EG_makeTopology(context, newCurves[jedge], EDGE, TWONODE,
                                 trange, 2, &(newNodes[jedge]), &(newSenses[jedge]), &newEdges[jedge]);
        CHECK_STATUS(EG_makeTopology);

        jedge++;
    }

    /* make a Loop from the Edges */
    status = EG_makeTopology(context, NULL, LOOP, mtype,
                             NULL, jedge, newEdges, newSenses, &newLoop);
    CHECK_STATUS(EG_makeToplogy);

    /* get the mtype of the input (original) Body */
    status = EG_getInfo(ebodyIn, &oclass, &mtype, &topRef, &eprev, &enext);
    CHECK_STATUS(EG_getInfo);

    /* if input was a WireBody, we can make the WireBody directly from the Loop */
    if (mtype == WIREBODY) {
        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &newLoop, NULL, ebodyOut);
        CHECK_STATUS(EG_makeTopology);

    /* otherwise, make the Face, Sheet, Shell, and SheetBody */
    } else {
        status = EG_makeFace(newLoop, newSenses[0], NULL, &newFace);
        CHECK_STATUS(EG_makeFace);

        status = EG_makeTopology(context, NULL, SHELL, OPEN,
                                 NULL, 1, &newFace, NULL, &newShell);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                                 NULL, 1, &newShell, NULL, ebodyOut);
        CHECK_STATUS(EG_makeTopology);

        SPLINT_CHECK_FOR_NULL(*ebodyOut);

        /* tell OpenCSM to put _body, _brch, and Branch Attributes on the Faces */
        status = EG_attributeAdd(*ebodyOut, "__markFaces__", ATTRSTRING, 1, NULL, NULL, "true");
        CHECK_STATUS(EG_attributeAdd);
    }

cleanup:
    if (eloops != NULL) EG_free(eloops);

    if (xyz != NULL) {
        for (iedge = 0; iedge < nedge; iedge++) {
            FREE(xyz[iedge]);
        }
    }

    FREE(ibeg     );
    FREE(iend     );
    FREE(tbeg     );
    FREE(tend     );
    FREE(active   );
    FREE(xyz      );
    FREE(newCurves);
    FREE(newNodes );
    FREE(newEdges );
    FREE(newSenses);

#ifdef GRAFIC
    FREE(ilin );
    FREE(isym );
    FREE(nper );
    FREE(xplot);
    FREE(yplot);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   interesctLines - find intersection of two discrete lines           *
 *                                                                      *
 ************************************************************************
 */

static int
intersectLines(int    npnt1,            /* (in)  number of points in line 1 */
               double xyz1[],           /* (in)  xyz array       for line 1 */
               int    npnt2,            /* (in)  number of points in line 2 */
               double xyz2[],           /* (in)  xyz array       for line 2 */
               int    iplane,           /* (in)  1=yz, 2=zx, 3=xy */
               double *t1,              /* (out) parametric coordinate for line 1 */
               double *t2)              /* (out) parametric coordinate for line 2 */
{
    int status = EGADS_SUCCESS;         /* (out) return status */

    int    ipnt1, ipnt2, redo;
    double xA, yA, xB, yB, xC, yC, xD, yD, det, s1, s2;

    ROUTINE(intersect);

    /* --------------------------------------------------------------- */

    /* default returns */
    *t1 = -HUGEQ;
    *t2 = +HUGEQ;

    ipnt1 = npnt1 - 2;                  // start at end of line 1
    ipnt2 = 0;                          // start at beg of line 2

    redo = 1;
    while (redo != 0) {
        redo = 0;

        /* get coordinates */
        if (iplane == 1) {
            xA = xyz1[3*ipnt1+1];
            yA = xyz1[3*ipnt1+2];
            xB = xyz1[3*ipnt1+4];
            yB = xyz1[3*ipnt1+5];

            xC = xyz2[3*ipnt2+1];
            yC = xyz2[3*ipnt2+2];
            xD = xyz2[3*ipnt2+4];
            yD = xyz2[3*ipnt2+5];
        } else if (iplane == 2) {
            xA = xyz1[3*ipnt1+2];
            yA = xyz1[3*ipnt1  ];
            xB = xyz1[3*ipnt1+5];
            yB = xyz1[3*ipnt1+3];

            xC = xyz2[3*ipnt2+2];
            yC = xyz2[3*ipnt2  ];
            xD = xyz2[3*ipnt2+5];
            yD = xyz2[3*ipnt2+3];
        } else {
            xA = xyz1[3*ipnt1  ];
            yA = xyz1[3*ipnt1+1];
            xB = xyz1[3*ipnt1+3];
            yB = xyz1[3*ipnt1+4];

            xC = xyz2[3*ipnt2  ];
            yC = xyz2[3*ipnt2+1];
            xD = xyz2[3*ipnt2+3];
            yD = xyz2[3*ipnt2+4];
        }

        /* find the intersection */
        det = (xB-xA) * (yC-yD) - (xC-xD) * (yB-yA);
        if (fabs(det) < EPS30) {
            if (fabs(xB-xC) < EPS06 && fabs(yB-yC) < EPS06) {
                *t1 = npnt1 - 1;
                *t2 = 0;
            } else {
                snprintf(message, 1014, "no intersection");
                status = OCSM_UDP_ERROR1;
            }

            goto cleanup;
        } else {
            s1 = ((xC-xA) * (yC-yD) - (xC-xD) * (yC-yA)) / det;
            s2 = ((xB-xA) * (yC-yA) - (xC-xA) * (yB-yA)) / det;
        }

        /* if not in this segment in line 1, decrease ipnt1 if possible and redo */
        if (s1 < 0) {
            if (ipnt1 > 1) {
                ipnt1--;
                redo ++;
            } else {
                *t2 = ipnt2 + s2;
                goto cleanup;
            }
        }

        /* if not in this segment in line 2, increase ipnt2 if possible and redo */
        if (s2 > 1) {
            if (ipnt2 < npnt2-2) {
                ipnt2++;
                redo ++;
            } else {
                *t1 = ipnt1 + s1;
                goto cleanup;
            }
        }
    }

    *t1 = ipnt1 + s1;
    *t2 = ipnt2 + s2;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeSolidOffset - build offset Body for a SolidBody                *
 *                                                                      *
 ************************************************************************
 */

static int
makeSolidOffset(ego    ebodyIn,
                ego    *ebodyOut,
    /*@unused@*/int    *NumUdp,
                udp_T  *udps)
{
    int     status = EGADS_SUCCESS;

    int      i, iface, iedge, jedge, medge, inode, jnode;
    int      itry, ntry=10, iswap;
    int      *ibegList=NULL, *iendList=NULL;
    double   data1[18], data2[18], tbreak, arclength;
    double   xyz0[18], uv1[2], uv2[2], xyz_temp[18], fact1, fact2, dist, dtest, den, t1;
    ego      context, epcurvei, epcurvej, etemp;
    wrep_T   *wrep=NULL;
    offset_T *offset=NULL;

#ifdef GRAFIC
    int      io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int      nline, nplot, *ilin=NULL, *isym=NULL, *nper=NULL;
    int      ipnt, npnt=51;
    float    *uplot=NULL, *vplot=NULL;
    double   tt, xyz_out[3], uv_out[2], data[18];
    char     pltitl[80];
#endif

    ROUTINE(makeSolidOffset);

    /* --------------------------------------------------------------- */

    /* check inputs */
    if (udps[0].arg[0].size != udps[0].arg[1].size) {
        snprintf(message, 1024, "\"nodelist\" and \"nodedist must have the same number of entries");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    if (udps[0].arg[0].size > 1 || NODELIST(0,0) != 0) {
        for (i = 0; i < udps[0].arg[0].size; i++) {
            if (NODELIST(0,i) <= 0) {
                snprintf(message, 1024, "\"nodelist\" must contain all positive values");
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }

        for (i = 0; i < udps[0].arg[1].size; i++) {
            if (NODEDIST(0,i) <= 0) {
                snprintf(message, 1024, "\"nodedist\" must contain all positive values");
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }
    }

    for (i = 0; i < udps[0].arg[2].size; i++) {
        if (EDGELIST(0,i) <= 0) {
            snprintf(message, 1024, "\"edgelist\" must contain all positive values");
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }

    for (i = 0; i < udps[0].arg[3].size; i++) {
        if (FACELIST(0,i) <= 0) {
            snprintf(message, 1024, "\"facelist\" must contain all positive values");
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }

    if (udps[0].arg[4].size != 1) {
        snprintf(message, 1024, "\"dist\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (DIST(0) <= 0) {
        snprintf(message, 1024, "\"dist\" must be positive");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getContext(ebodyIn, &context);
    CHECK_STATUS(EG_getContext);

    /* make the winged-edge data structure */
    status = wrepInit(ebodyIn, &wrep);
    CHECK_STATUS(wrepInit);

    SPLINT_CHECK_FOR_NULL(wrep);

    /* mark the Nodes, Edges, and Faces that are in the nodeList, edgeList, and faceList */
    for (inode = 1; inode <= wrep->nnode; inode++) {
        wrep->node[inode].dist = DIST(0);
    }

    if (udps[0].arg[0].size > 1 || NODELIST(0,0) != 0) {
        for (i = 0; i < udps[0].arg[0].size; i++) {
            inode = NODELIST(0,i);
            if (inode >= 1 && inode <= wrep->nnode) {
                wrep->node[inode].dist = NODEDIST(0,i);
            } else {
                snprintf(message, 1024, "\"nodelist[%d]\" (%d) is out of bounds", i+1, NODELIST(0,i));
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }
    }

    for (i = 0; i < udps[0].arg[2].size; i++) {
        iedge = EDGELIST(0,i);
        if (iedge >= 1 && iedge <= wrep->nedge) {
            wrep->edge[iedge].tag = 1;
        } else {
            snprintf(message, 1024, "\"edgelist[%d]\" (%d) is out of bounds", i+1, EDGELIST(0,i));
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }

    for (i = 0; i < udps[0].arg[3].size; i++) {
        iface = FACELIST(0, i);
        if (iface >= 1 && iface <= wrep->nface) {
            wrep->face[iface].tag = 1;
        } else {
            snprintf(message, 1024, "\"facelist[%d]\" (%d) is out of bounds", i+1, FACELIST(0,i));
            status = EGADS_RANGERR;
            goto cleanup;
        }
    }

#ifdef DEBUG
    printf("initial wrep\n");
    status = wrepPrint(wrep);
    CHECK_STATUS(wrepPrint);
#endif

#ifdef GRAFIC
    MALLOC(ilin,  int,   5     *wrep->nedge);
    MALLOC(isym,  int,   5     *wrep->nedge);
    MALLOC(nper,  int,   5     *wrep->nedge);
    MALLOC(uplot, float, 5*npnt*wrep->nedge);
    MALLOC(vplot, float, 5*npnt*wrep->nedge);
#endif

    /* create array of possible offset Wedges for any Wface */
    medge = 3 * (wrep->nedge+1);
    MALLOC(offset, offset_T, medge);

    /* split any Wedge that is not tagged and which adjoins a
       Wedge that is */
    for (iedge = 1; iedge <= wrep->nedge; iedge++) {
        if (wrep->edge[iedge].tag != 0) continue;           // tagged     Wedge
        if (wrep->edge[iedge].ecurve == NULL) continue;     // degenerate Wedge

        /* (possibly) split at beginning */
        if ((wrep->face[wrep->edge[iedge].ileft ].tag == 1 &&
             wrep->edge[wrep->edge[iedge].ibleft].tag == 1   ) ||
            (wrep->face[wrep->edge[iedge].irite ].tag == 1 &&
             wrep->edge[wrep->edge[iedge].ibrite].tag == 1   )   ) {
            status = EG_arcLength(wrep->edge[iedge].ecurve,
                                  wrep->edge[iedge].tbeg, wrep->edge[iedge].tend, &arclength);
            CHECK_STATUS(EG_arclength);

            dist = wrep->node[wrep->edge[iedge].ibeg].dist;

            if (fabs(dist-arclength) < EPS06) {
                /* splitting is not required */
            } else if (dist < arclength) {
                tbreak = wrep->edge[iedge].tbeg
                      + (wrep->edge[iedge].tend - wrep->edge[iedge].tbeg) * dist / arclength;

                for (itry = 0; itry < ntry; itry++) {
                    status = EG_arcLength(wrep->edge[iedge].ecurve,
                                          wrep->edge[iedge].tbeg, tbreak, &arclength);
                    CHECK_STATUS(EG_arcLength);

                    if (fabs(dist-arclength) < EPS06) break;

                    tbreak = wrep->edge[iedge].tbeg
                        + (tbreak - wrep->edge[iedge].tbeg) * dist / arclength;
                }

                status = wrepBreakEdge(wrep, iedge, tbreak, &jnode);
                CHECK_STATUS(wrepBreakEdge);
            } else {
                snprintf(message, 1024, "Wedge %d is shorter than offset distance", iedge);
                status = OCSM_UDP_ERROR1;
                goto cleanup;
            }
        }

        /* split at end */
        if ((wrep->face[wrep->edge[iedge].ileft ].tag == 1 &&
             wrep->edge[wrep->edge[iedge].ieleft].tag == 1   ) ||
            (wrep->face[wrep->edge[iedge].irite ].tag == 1 &&
             wrep->edge[wrep->edge[iedge].ierite].tag == 1   )   ) {
            status = EG_arcLength(wrep->edge[iedge].ecurve,
                                  wrep->edge[iedge].tbeg, wrep->edge[iedge].tend, &arclength);
            CHECK_STATUS(EG_arclength);

            dist = wrep->node[wrep->edge[iedge].iend].dist;

            if (fabs(dist-arclength) < EPS06) {
                /* splitting is not required */
            } else if (dist < arclength) {
                tbreak = wrep->edge[iedge].tend
                      + (wrep->edge[iedge].tbeg - wrep->edge[iedge].tend) * dist / arclength;

                for (itry = 0; itry < ntry; itry++) {
                    status = EG_arcLength(wrep->edge[iedge].ecurve,
                                          tbreak, wrep->edge[iedge].tend, &arclength);
                    CHECK_STATUS(EG_arcLength);

                    if (fabs(dist-arclength) < EPS06) break;

                    tbreak = wrep->edge[iedge].tend
                        + (tbreak - wrep->edge[iedge].tend) * dist / arclength;
                }

                status = wrepBreakEdge(wrep, iedge, tbreak, &jnode);
                CHECK_STATUS(wrepBreakEdge);
            } else {
                snprintf(message, 1024, "Wedge %d is shorter than offset distance", iedge);
                status = OCSM_UDP_ERROR1;
                goto cleanup;
            }

        }
    }

#ifdef DEBUG
        printf("after splitting Wedges\n");
        status = wrepPrint(wrep);
        CHECK_STATUS(wrepPrint);
#endif

    /* process each tagged Wface */
    for (iface = 1; iface <= wrep->nface; iface++) {
        if (wrep->face[iface].tag == 0) continue;

#ifdef GRAFIC
        nplot = 0;
        nline = 0;

        /* draw outline */
        for (iedge = 1; iedge <= wrep->nedge; iedge++) {
            if (wrep->edge[iedge].ileft != iface &&
                wrep->edge[iedge].irite != iface   ) continue;

            jedge = iedge;
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                double uvout[2], xyzout[3];
                tt = wrep->edge[iedge].tbeg + (wrep->edge[iedge].tend - wrep->edge[iedge].tbeg) * (double)(ipnt) / (double)(npnt-1);

                status = EG_evaluate(wrep->edge[iedge].eedge, &tt, data);
                CHECK_STATUS(EG_evaluate);

                status = EG_invEvaluate(wrep->face[iface].esurface, data, uvout, xyzout);
                CHECK_STATUS(EG_invEvaluate);

                uplot[nplot] = uvout[0];
                vplot[nplot] = uvout[1];
                nplot++;
            }

            ilin[nline] = +GR_SOLID;
            isym[nline] = -GR_PLUS;
            nper[nline] = npnt;
            nline++;
        }
#endif

        /* initialize the offsets for this Wface */
        for (iedge = 0; iedge < medge; iedge++) {
            offset[iedge].iface  = 0;
            offset[iedge].ibeg   = 0;
            offset[iedge].iend   = 0;
            offset[iedge].tbeg   = 0;
            offset[iedge].tend   = 0;
            offset[iedge].ecurve = NULL;
        }

        /* loop through each Wedge that is tagged and that adjoins this Wface */
        for (iedge = wrep->nedge; iedge > 0; iedge--) {
            if (wrep->edge[iedge].tag == 0) continue;

            /* jedge is the previous Wedge */
            if        (wrep->edge[iedge].ileft == iface) {
                jedge = wrep->edge[iedge].ibleft;
                dist  = wrep->node[wrep->edge[iedge].ibeg].dist;
                fact1 = + dist * wrep->face[iface].mtype;
            } else if (wrep->edge[iedge].irite == iface) {
                jedge = wrep->edge[iedge].ibrite;
                dist  = wrep->node[wrep->edge[iedge].ibeg].dist;
                fact1 = - dist * wrep->face[iface].mtype;
            } else {
                continue;
            }

            /* no need to process if ibeg is already set */
            if (offset[iedge].ibeg != 0) {

            /* if the previous Wedge (jedge) is not tagged, set offset.ibeg to
               the end of jedge */
            } else if (wrep->edge[jedge].tag == 0) {
                if        (wrep->edge[jedge].ibeg == wrep->edge[iedge].ibeg) {
                    offset[iedge].ibeg = wrep->edge[jedge].iend;
                } else if (wrep->edge[jedge].iend == wrep->edge[iedge].ibeg) {
                    offset[iedge].ibeg = wrep->edge[jedge].ibeg;
                } else {
                    status = OCSM_INTERNAL_ERROR;
                    snprintf(message, 1024, "expecting Wedge %d to be attached to Wnode %d",
                             jedge, wrep->edge[iedge].iend);
                    goto cleanup;
                }

#ifdef GRAFIC
                status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[offset[iedge].ibeg].xyz, uv_out, xyz_out);
                CHECK_STATUS(EG_invEvaluate);

                uplot[nplot] = uv_out[0];
                vplot[nplot] = uv_out[1];
                nplot++;

                ilin[nline] = 0;
                isym[nline] = GR_STAR;
                nper[nline] = 1;
                nline++;
#endif

            /* if the previous Wedge (jedge) is tagged, find the intersection
               of offsets, create a Wnode and a Wedge */
            } else {
                status = EG_otherCurve(wrep->face[iface].esurface, wrep->edge[iedge].eedge, 0, &epcurvei);
                CHECK_STATUS(EG_otherCurve);

                status = EG_otherCurve(wrep->face[iface].esurface, wrep->edge[jedge].eedge, 0, &epcurvej);
                CHECK_STATUS(EG_otherCurve);

                /* beg of iedge, end of jedge */
                if ((wrep->edge[iedge].ileft == iface && wrep->edge[jedge].ileft == iface) ||
                    (wrep->edge[iedge].irite == iface && wrep->edge[jedge].irite == iface)   ) {
                    status = EG_evaluate(epcurvei, &wrep->edge[iedge].tbeg, data1);
                    CHECK_STATUS(EG_evaluate);

                    status = EG_evaluate(epcurvej, &wrep->edge[jedge].tend, data2);
                    CHECK_STATUS(EG_evaluate);

                    fact2 = +fact1;

                /* beg of iedge, beg of jedge */
                } else if ((wrep->edge[iedge].ileft == iface && wrep->edge[jedge].irite == iface) ||
                           (wrep->edge[iedge].irite == iface && wrep->edge[jedge].ileft == iface)   ) {
                    status = EG_evaluate(epcurvei, &wrep->edge[iedge].tbeg, data1);
                    CHECK_STATUS(EG_evaluate);

                    status = EG_evaluate(epcurvej, &wrep->edge[jedge].tbeg, data2);
                    CHECK_STATUS(EG_evaluate);

                    fact2 = -fact1;

                } else {
                    snprintf(message, 1024, "iedge=%d and jedge=%d do not have a common Wnode", iedge, jedge);
                    status = OCSM_INTERNAL_ERROR;
                    goto cleanup;
                }

                status = EG_evaluate(wrep->face[iface].esurface, data1, xyz0);
                CHECK_STATUS(EG_evaluate);

                for (itry = 0; itry < ntry; itry++) {
                    uv1[0] = data1[0] - fact1 * data1[3];
                    uv1[1] = data1[1] + fact1 * data1[2];

                    status = EG_evaluate(wrep->face[iface].esurface, uv1, xyz_temp);
                    CHECK_STATUS(EG_evaluate);

                    dtest = sqrt((xyz_temp[0]-xyz0[0]) * (xyz_temp[0]-xyz0[0])
                                +(xyz_temp[1]-xyz0[1]) * (xyz_temp[1]-xyz0[1])
                                +(xyz_temp[2]-xyz0[2]) * (xyz_temp[2]-xyz0[2]));

                    if (fabs(dtest-dist) < EPS06) break;
                    fact1 = fact1 * dist / dtest;
                }

                for (itry = 0; itry < ntry; itry++) {
                    uv2[0] = data2[0] - fact2 * data2[3];
                    uv2[1] = data2[1] + fact2 * data2[2];

                    status = EG_evaluate(wrep->face[iface].esurface, uv2, xyz_temp);
                    CHECK_STATUS(EG_evaluate);

                    dtest = sqrt((xyz_temp[0]-xyz0[0]) * (xyz_temp[0]-xyz0[0])
                                +(xyz_temp[1]-xyz0[1]) * (xyz_temp[1]-xyz0[1])
                                +(xyz_temp[2]-xyz0[2]) * (xyz_temp[2]-xyz0[2]));

                    if (fabs(dtest-dist) < EPS06) break;
                    fact2 = fact2 * dist / dtest;
                }

                den = (data1[3] * data2[2] - data1[2] * data2[3]);
                if (fabs(den) < EPS12) {
                    uv1[0] = (uv1[0] + uv2[0]) / 2;
                    uv1[1] = (uv1[1] + uv2[1]) / 2;
                } else {
                    t1 = ((uv2[1] - uv1[1]) * data2[2] - (uv2[0] - uv1[0]) * data2[3]) / den;
                    uv1[0] = uv1[0] + data1[2] * t1;
                    uv1[1] = uv1[1] + data1[3] * t1;
                }

#ifdef GRAFIC
                uplot[nplot] = data1[0];
                vplot[nplot] = data1[1];
                nplot++;

                uplot[nplot] = uv1[0];
                vplot[nplot] = uv1[1];
                nplot++;

                ilin[nline] = GR_DOTTED;
                isym[nline] = GR_CIRCLE;
                nper[nline] = 2;
                nline++;
#endif

                status = EG_evaluate(wrep->face[iface].esurface, uv1, data1);
                CHECK_STATUS(EG_evaluate);

                jnode = status = wrepMakeNode(wrep, iface, data1);
                CHECK_STATUS(wrepMakeNode);

                status = wrepMakeEdge(wrep, iface, wrep->edge[iedge].ibeg, jnode, NULL);
                CHECK_STATUS(wrepMakeEdge);

                offset[iedge].ibeg = jnode;

                if        (wrep->edge[jedge].ibeg == wrep->edge[iedge].ibeg) {
                    offset[jedge].ibeg = jnode;
                } else if (wrep->edge[jedge].iend == wrep->edge[iedge].ibeg) {
                    offset[jedge].iend = jnode;
                }
            }

            /* jedge is the next Wedge */
            if  (wrep->edge[iedge].ileft == iface) {
                jedge = wrep->edge[iedge].ieleft;
            } else {
                jedge = wrep->edge[iedge].ierite;
            }

            /* no need to process if iend is already set */
            if (offset[iedge].iend != 0) {

            /* if the next Wedge (jedge) is not tagged, set offset.iend to
               the end of jedge */
            } else if (wrep->edge[jedge].tag == 0) {
                if        (wrep->edge[jedge].ibeg == wrep->edge[iedge].iend) {
                    offset[iedge].iend = wrep->edge[jedge].iend;
                } else if (wrep->edge[jedge].iend == wrep->edge[iedge].iend) {
                    offset[iedge].iend = wrep->edge[jedge].ibeg;
                } else {
                    status = OCSM_INTERNAL_ERROR;
                    snprintf(message, 1024, "expecting Wedge %d to be attached to Wnode %d",
                             jedge, wrep->edge[iedge].iend);
                    goto cleanup;
                }

#ifdef GRAFIC
                status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[offset[iedge].iend].xyz, uv_out, xyz_out);
                CHECK_STATUS(EG_invEvaluate);

                uplot[nplot] = uv_out[0];
                vplot[nplot] = uv_out[1];
                nplot++;

                ilin[nline] = 0;
                isym[nline] = GR_STAR;
                nper[nline] = 1;
                nline++;
#endif

            /* if the next Wedge (jedge) is tagged, find the intersection
               of offsets, create a Wnode and a Wedge */
            } else {
                status = EG_otherCurve(wrep->face[iface].esurface, wrep->edge[iedge].eedge, 0, &epcurvei);
                CHECK_STATUS(EG_otherCurve);

                status = EG_otherCurve(wrep->face[iface].esurface, wrep->edge[jedge].eedge, 0, &epcurvej);
                CHECK_STATUS(EG_otherCurve);

                /* end of iedge, beg of jedge */
                if ((wrep->edge[iedge].ileft == iface && wrep->edge[jedge].ileft == iface) ||
                    (wrep->edge[iedge].irite == iface && wrep->edge[jedge].irite == iface)   ) {
                    status = EG_evaluate(epcurvei, &wrep->edge[iedge].tend, data1);
                    CHECK_STATUS(EG_evaluate);

                    status = EG_evaluate(epcurvej, &wrep->edge[jedge].tbeg, data2);
                    CHECK_STATUS(EG_evaluate);

                    fact2 = +fact1;

                /* end of iedge, end of jedge */
                } else if ((wrep->edge[iedge].ileft == iface && wrep->edge[jedge].irite == iface) ||
                           (wrep->edge[iedge].irite == iface && wrep->edge[jedge].ileft == iface)   ) {
                    status = EG_evaluate(epcurvei, &wrep->edge[iedge].tend, data1);
                    CHECK_STATUS(EG_evaluate);

                    status = EG_evaluate(epcurvej, &wrep->edge[jedge].tend, data2);
                    CHECK_STATUS(EG_evaluate);

                    fact2 = -fact1;

                } else {
                    snprintf(message, 1024, "iedge=%d and jedge=%d do not have a common Wnode", iedge, jedge);
                    status = OCSM_INTERNAL_ERROR;
                    goto cleanup;
                }

                status = EG_evaluate(wrep->face[iface].esurface, data1, xyz0);
                CHECK_STATUS(EG_evaluate);

                for (itry = 0; itry < ntry; itry++) {
                    uv1[0] = data1[0] - fact1 * data1[3];
                    uv1[1] = data1[1] + fact1 * data1[2];

                    status = EG_evaluate(wrep->face[iface].esurface, uv1, xyz_temp);
                    CHECK_STATUS(EG_evaluate);

                    dtest = sqrt((xyz_temp[0]-xyz0[0]) * (xyz_temp[0]-xyz0[0])
                                +(xyz_temp[1]-xyz0[1]) * (xyz_temp[1]-xyz0[1])
                                +(xyz_temp[2]-xyz0[2]) * (xyz_temp[2]-xyz0[2]));

                    if (fabs(dtest-dist) < EPS06) break;
                    fact1 = fact1 * dist / dtest;
                }

                for (itry = 0; itry < ntry; itry++) {
                    uv2[0] = data2[0] - fact2 * data2[3];
                    uv2[1] = data2[1] + fact2 * data2[2];

                    status = EG_evaluate(wrep->face[iface].esurface, uv2, xyz_temp);
                    CHECK_STATUS(EG_evaluate);

                    dtest = sqrt((xyz_temp[0]-xyz0[0]) * (xyz_temp[0]-xyz0[0])
                                +(xyz_temp[1]-xyz0[1]) * (xyz_temp[1]-xyz0[1])
                                +(xyz_temp[2]-xyz0[2]) * (xyz_temp[2]-xyz0[2]));

                    if (fabs(dtest-dist) < EPS06) break;
                    fact2 = fact2 * dist / dtest;
                }

                den = (data1[3] * data2[2] - data1[2] * data2[3]);
                if (fabs(den) < EPS12) {
                    uv1[0] = (uv1[0] + uv2[0]) / 2;
                    uv1[1] = (uv1[1] + uv2[1]) / 2;
                } else {
                    t1 = ((uv2[1] - uv1[1]) * data2[2] - (uv2[0] - uv1[0]) * data2[3]) / den;
                    uv1[0] = uv1[0] + data1[2] * t1;
                    uv1[1] = uv1[1] + data1[3] * t1;
                }

#ifdef GRAFIC
                uplot[nplot] = data1[0];
                vplot[nplot] = data1[1];
                nplot++;

                uplot[nplot] = uv1[0];
                vplot[nplot] = uv1[1];
                nplot++;

                ilin[nline] = GR_DOTTED;
                isym[nline] = GR_CIRCLE;
                nper[nline] = 2;
                nline++;
#endif

                status = EG_evaluate(wrep->face[iface].esurface, uv1, data1);
                CHECK_STATUS(EG_evaluate);

                jnode = status = wrepMakeNode(wrep, iface, data1);
                CHECK_STATUS(wrepMakeNode);

                status = wrepMakeEdge(wrep, iface, wrep->edge[iedge].iend, jnode, NULL);
                CHECK_STATUS(wrepMakeEdge);

                offset[iedge].iend = jnode;

                if        (wrep->edge[jedge].ibeg == wrep->edge[iedge].iend) {
                    offset[jedge].ibeg = jnode;
                } else if (wrep->edge[jedge].iend == wrep->edge[iedge].iend) {
                    offset[jedge].iend = jnode;
                }
            }

            /* create an offset curve and associated Wedge (which will make a Wface) */
            status = makeOffsetCurve(wrep, iface, iedge, &offset[iedge]);
            CHECK_STATUS(makeOffsetCurve);
        }

#ifdef DEBUG
        for (iedge = 1; iedge <= wrep->nedge; iedge++) {
            printf("iedge=%3d, iface=%3d, ibeg=%3d, iend=%3d, tbeg=%10.5f, tend=%10.5f, ecurve=%llx\n",
                   iedge, offset[iedge].iface, offset[iedge].ibeg, offset[iedge].iend, offset[iedge].tbeg, offset[iedge].tend, (long long)(offset[iedge].ecurve));
        }
#endif

        for (iedge = 1; iedge <= wrep->nedge; iedge++) {
            if (offset[iedge].ecurve == NULL) continue;

            /* orient the new Wedge so that its right side faces it Wedge */
            if        (wrep->edge[iedge].ileft == iface) {
                /* do nothing */
            } else if (wrep->edge[iedge].irite == iface) {
#ifdef DEBUG
                printf("swapping ends for iedge=%d\n", iedge);
#endif
                iswap              = offset[iedge].ibeg;
                offset[iedge].ibeg = offset[iedge].iend;
                offset[iedge].iend = iswap;

                status = EG_flipObject(offset[iedge].ecurve, &etemp);
                CHECK_STATUS(EG_flipObject);

                status = EG_deleteObject(offset[iedge].ecurve);
                CHECK_STATUS(EG_deleteObject);

                offset[iedge].ecurve = etemp;
            } else {
                snprintf(message, 1024, "neither ileft nor irite of iedge=%d is set to iface=%d", iedge, iface);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            }

            status = wrepMakeEdge(wrep, iface, offset[iedge].ibeg, offset[iedge].iend, offset[iedge].ecurve);
            CHECK_STATUS(wrepMakeEdge);

#ifdef GRAFIC
            jedge = status;
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                double uvout[2], xyzout[3];
                tt = wrep->edge[jedge].tbeg + (wrep->edge[jedge].tend - wrep->edge[jedge].tbeg) * (double)(ipnt) / (double)(npnt-1);

                status = EG_evaluate(wrep->edge[jedge].ecurve, &tt, data);
                CHECK_STATUS(EG_evaluate);

                status = EG_invEvaluate(wrep->face[iface].esurface, data, uvout, xyzout);
                CHECK_STATUS(EG_invEvaluate);

                uplot[nplot] = uvout[0];
                vplot[nplot] = uvout[1];
                nplot++;
            }

            ilin[nline] = +GR_DASHED;
            isym[nline] = -GR_PLUS;
            nper[nline] = npnt;
            nline++;
#endif
        }

#ifdef GRAFIC
        snprintf(pltitl, 79, "udfScribe --- offset ends");
        grinit_(&io_kbd, &io_scr, pltitl, STRLEN(pltitl));

        if (wrep->face[iface].mtype > 0) {
            snprintf(pltitl, 79, "~u~v~Face %d", iface);
        } else {
            for (i = 0; i < nplot; i++) {
                uplot[i] *= -1;
            }
            snprintf(pltitl, 79, "~-u~v~Face %d", iface);
        }
        grline_(ilin, isym, &nline, pltitl, &indgr, uplot, vplot, nper, STRLEN(pltitl));
#endif
    }

#ifdef DEBUG
    printf("final wrep\n");
    status = wrepPrint(wrep);
    CHECK_STATUS(wrepPrint);

    for (iface = 1; iface <= wrep->nface; iface++) {
        printf("Wface %3d: Wedges", iface);
        for (iedge = 1; iedge <= wrep->nedge; iedge++) {
            if (wrep->edge[iedge].ileft == iface) {
                printf(" %3d", -iedge);
            } else if (wrep->edge[iedge].irite == iface) {
                printf(" %3d", +iedge);
            }
        }
        printf("\n");
    }
#endif

    /* put the "scribeFace" attribute on Faces that are within
       the scribe */
    for (iface = 1; iface <= wrep->nface; iface++) {
        wrep->face[iface].new = 0;

        if (wrep->face[iface].tag != 0 || iface > wrep->nface_orig) {
            for (iedge = 1; iedge <= wrep->nedge; iedge++) {
                if (wrep->edge[iedge].tag == 0) continue;

                if (wrep->edge[iedge].ileft == iface ||
                    wrep->edge[iedge].irite == iface   ) {
                    wrep->face[iface].new = 1;
                }
            }
        }
    }

    FREE(offset);

#ifdef DEBUG
    printf("final wrep\n");
    status = wrepPrint(wrep);
    CHECK_STATUS(wrepPrint);
#endif

    /* try to make a new ebody */
    status = wrep2ego(*wrep, context, ebodyOut);
    CHECK_STATUS(wrep2ego);

    SPLINT_CHECK_FOR_NULL(*ebodyOut);

#ifdef DEBUG
    printf("ebody\n");
    ocsmPrintEgo(*ebodyOut);
#endif

    /* call to avoid compiler warning */
    if (wrep->nnode == -2 && wrep->nedge == -2 && wrep->nface == -2) {
        status = wrepPrint(wrep);
        CHECK_STATUS(wrepPrint);
    }

cleanup:
    /* free the WREP structure */
    wrepFree(wrep);
    FREE(wrep);

    FREE(offset);
    FREE(ibegList);
    FREE(iendList);

#ifdef GRAFIC
    FREE(ilin);
    FREE(isym);
    FREE(nper);
    FREE(uplot);
    FREE(vplot);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeOffsetCurve - make an (untrimmed) offset curve                 *
 *                                                                      *
 ************************************************************************
 */

static int
makeOffsetCurve(wrep_T    *wrep,        /* (in)  pointer to WREP structure */
                int       iface,        /* (in)  Wface index (bias-1) */
                int       iedge,        /* (in)  Wedge index (bias-1) */
                offset_T  *offset)      /* (in)  offset information */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     npnt, ipnt;
    int     oclass, mtype, sizes[2], i, *header;
    double  data1[18], tt, *gdata;
    double  tang[2], duv[2], duv_beg[2], duv_end[2];
    double  data[18], uv_out[2], xyz_out[3], frac;
    double  *uv1=NULL, *uv2=NULL, *xyz=NULL;
    double  toler=EPS06;
    ego     epcurve, topRef, prev, next, context, etemp, eref;

    ROUTINE(makeOffsetCurve);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter makeOffsetCurve(iface=%d, iedge=%d)\n", iface, iedge);
#endif

    npnt = 21;
    MALLOC(uv1, double, 2*npnt);        /* parametric coordinates of the Wedge */
    MALLOC(uv2, double, 3*npnt);        /* parametric coordinates of the offset */
    MALLOC(xyz, double, 3*npnt);        /* physical   coordinates of the offset */

    /* get the pcurve associated with this iedge,iface */
    status = EG_getInfo(wrep->face[iface].esurface, &oclass, &mtype, &topRef, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    if (mtype == PLANE) {
        status = EG_otherCurve(wrep->face[iface].esurface, wrep->edge[iedge].ecurve, 0.0, &epcurve);
        CHECK_STATUS(EG_otherCurve);
    } else if (wrep->edge[iedge].ileft == iface) {
        epcurve = wrep->edge[iedge].epleft;
    } else if (wrep->edge[iedge].irite == iface) {
        epcurve = wrep->edge[iedge].eprite;
    } else {
        snprintf(message, 1024, "Wface %d is not planar and does not point to Wedge %d", iface, iedge);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* find the uv along the Wedge (to be used below as the anchor
       for each tuft) */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        tt = wrep->edge[iedge].tbeg
           + (wrep->edge[iedge].tend - wrep->edge[iedge].tbeg) * (double)(ipnt) / (double)(npnt-1);

        status = EG_evaluate(epcurve, &tt, data1);
        CHECK_STATUS(EG_evaluate);

        uv1[2*ipnt  ] = data1[0];
        uv1[2*ipnt+1] = data1[1];
    }

    /* find the difference in UV between the ends of the Wedge and the ends of the offset */
    status = EG_evaluate(epcurve, &wrep->edge[iedge].tbeg, data);
    CHECK_STATUS(EG_evaluate);

    tang[0] = data[2] / sqrt(data[2] * data[2] + data[3] * data[3]);     // norm[1] = +tang[0]
    tang[1] = data[3] / sqrt(data[2] * data[2] + data[3] * data[3]);     // norm[0] = -tang[1]

    status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[offset->ibeg].xyz, uv_out, xyz_out);
    CHECK_STATUS(EG_invEvaluate);

    duv_beg[0] = +tang[0] * (uv_out[0] - uv1[0]) + tang[1] * (uv_out[1] - uv1[1]);
    duv_beg[1] = -tang[1] * (uv_out[0] - uv1[0]) + tang[0] * (uv_out[1] - uv1[1]);


    status = EG_evaluate(epcurve, &wrep->edge[iedge].tend, data);
    CHECK_STATUS(EG_evaluate);

    tang[0] = data[2] / sqrt(data[2] * data[2] + data[3] * data[3]);
    tang[1] = data[3] / sqrt(data[2] * data[2] + data[3] * data[3]);

    status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[offset->iend].xyz, uv_out, xyz_out);
    CHECK_STATUS(EG_invEvaluate);

    duv_end[0] = +tang[0] * (uv_out[0] - uv1[2*npnt-2]) + tang[1] * (uv_out[1] - uv1[2*npnt-1]);
    duv_end[1] = -tang[1] * (uv_out[0] - uv1[2*npnt-2]) + tang[0] * (uv_out[1] - uv1[2*npnt-1]);

    /* straight pcurve (for now) */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        frac   = (double)(ipnt) / (double)(npnt-1);
        tt     = (1-frac) * wrep->edge[iedge].tbeg + frac * wrep->edge[iedge].tend;
        duv[0] = (1-frac) * duv_beg[0]             + frac * duv_end[0];
        duv[1] = (1-frac) * duv_beg[1]             + frac * duv_end[1];

        status = EG_evaluate(epcurve, &tt, data);
        CHECK_STATUS(EG_evaluate);

        tang[0] = data[2] / sqrt(data[2] * data[2] + data[3] * data[3]);
        tang[1] = data[3] / sqrt(data[2] * data[2] + data[3] * data[3]);

        uv2[3*ipnt  ] = uv1[2*ipnt  ] + tang[0] * duv[0] - tang[1] * duv[1];
        uv2[3*ipnt+1] = uv1[2*ipnt+1] + tang[1] * duv[0] + tang[0] * duv[1];
        uv2[3*ipnt+2] = 0;
    }

#ifdef GRAFIC
    if (0) {
        int      io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int      nline=0, nplot=0, ilin[2], isym[2], nper[2];
        float    uplot[200], vplot[200];
        char     pltitl[80];

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            uplot[nplot] = uv1[2*ipnt  ];
            vplot[nplot] = uv1[2*ipnt+1];
            nplot++;
        }
        ilin[nline] = GR_SOLID;
        isym[nline] = GR_SQUARE;
        nper[nline] = npnt;
        nline++;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            uplot[nplot] = uv2[3*ipnt  ];
            vplot[nplot] = uv2[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = GR_SOLID;
        isym[nline] = GR_PLUS;
        nper[nline] = npnt;
        nline++;

        snprintf(pltitl, 79, "udfScribe");
        grinit_(&io_kbd, &io_scr, pltitl, STRLEN(pltitl));

        snprintf(pltitl, 79, "~u~v~makeOffsetCurve: iface=%d, iedge=%d", iface, iedge);
        grline_(ilin, isym, &nline, pltitl, &indgr, uplot, vplot, nper, STRLEN(pltitl));


    }
#endif

    /* generate the coordinates that are a given offset distance from the Wedge */
    for (ipnt = 0; ipnt < npnt; ipnt++) {

        status = EG_evaluate(wrep->face[iface].esurface, &uv2[3*ipnt], data1);
        CHECK_STATUS(EG_evaluate);

        /* this is needed to get tbeg and tend correct below */
        xyz[3*ipnt  ] = data1[0];
        xyz[3*ipnt+1] = data1[1];
        xyz[3*ipnt+2] = data1[2];
    }

    /* generate the EGADS curve */
    status = EG_getContext(wrep->face[iface].esurface, &context);
    CHECK_STATUS(EG_getContext);

    sizes[0] = npnt;
    sizes[1] = 0;
    status = EG_approximate(context, 0, EPS06, sizes, uv2, &etemp);
    CHECK_STATUS(EG_approximate);

    /* convert etemp (which is a BSPLINE Curve into a BSPLINE Pcurve) */
    status = EG_getGeometry(etemp, &oclass, &mtype, &eref, &header, &gdata);
    CHECK_STATUS(EG_getGeometry);

    if (oclass != CURVE || mtype != BSPLINE) {
        snprintf(message, 1024, "etemp: oclass=%d, mtype=%d", oclass, mtype);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(header);
    SPLINT_CHECK_FOR_NULL(gdata );

    status = EG_deleteObject(etemp);
    CHECK_STATUS(EG_deleteObject);

    for (i = 0; i < header[2]; i++) {
        gdata[header[3]+2*i  ] = gdata[header[3]+3*i  ];
        gdata[header[3]+2*i+1] = gdata[header[3]+3*i+1];
    }

    status = EG_makeGeometry(context, PCURVE, BSPLINE, eref, header, gdata, &epcurve);
    CHECK_STATUS(EG_makeGeometry);

    EG_free(header);
    EG_free(gdata );

    /* convert to EGADS curve */
    status = EG_otherCurve(wrep->face[iface].esurface, epcurve, toler, &(offset->ecurve));
    CHECK_STATUS(EG_otherCurve);

cleanup:
    FREE(uv1);
    FREE(uv2);
    FREE(xyz);

#ifdef DEBUG
    printf("exit makeOffsetCurve -> status=%d\n", status);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepInit - initialize a WREP structure from an EGADS Body          *
 *                                                                      *
 ************************************************************************
 */

static int
wrepInit(ego    ebody,                  /* (in)  Body to transform */
         wrep_T **wrep)                 /* (out) pointer to new WREP structure (freeable) */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     oclass, mtype, nchild, ichild, *senses;
    int     inode, iedge, jedge, kedge, iface, nloop, iloop;
    double  data[118];
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL, eref, esurface, *echilds;
    ego     *eloops;
    wrep_T  *newWrep=NULL;

    ROUTINE(wrepInit);

    /* --------------------------------------------------------------- */

    *wrep = NULL;

    /* set up initial data structure */
    MALLOC(newWrep, wrep_T, 1);

    newWrep->node = NULL;
    newWrep->edge = NULL;
    newWrep->face = NULL;
    newWrep->ebody = ebody;

    status = EG_getBodyTopos(ebody, NULL, NODE, &(newWrep->nnode), &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody, NULL, EDGE, &(newWrep->nedge), &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody, NULL, FACE, &(newWrep->nface), &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    MALLOC(newWrep->node, wnode_T, newWrep->nnode+1);
    MALLOC(newWrep->edge, wedge_T, newWrep->nedge+1);
    MALLOC(newWrep->face, wface_T, newWrep->nface+1);

    /* fill in the Wnode data structure */
    for (inode = 1; inode <= newWrep->nnode; inode++) {
        SPLINT_CHECK_FOR_NULL(enodes);

        status = EG_getTopology(enodes[inode-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        newWrep->node[inode].xyz[0] = data[0];
        newWrep->node[inode].xyz[1] = data[1];
        newWrep->node[inode].xyz[2] = data[2];
        newWrep->node[inode].nedge  = 0;
        newWrep->node[inode].enode  = enodes[inode-1];
    }

    /* fill in the Wedge data structure */
    for (iedge = 1; iedge <= newWrep->nedge; iedge++) {
        SPLINT_CHECK_FOR_NULL(eedges);

        status = EG_getTopology(eedges[iedge-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        newWrep->edge[iedge].tag    = 0;
        newWrep->edge[iedge].ibeg   = 0;
        newWrep->edge[iedge].iend   = 0;
        newWrep->edge[iedge].ileft  = 0;
        newWrep->edge[iedge].irite  = 0;
        newWrep->edge[iedge].ibleft = 0;
        newWrep->edge[iedge].ibrite = 0;
        newWrep->edge[iedge].ieleft = 0;
        newWrep->edge[iedge].ierite = 0;
        newWrep->edge[iedge].tbeg   = data[0];
        newWrep->edge[iedge].tend   = data[1];
        newWrep->edge[iedge].ecurve = eref;
        newWrep->edge[iedge].epleft = NULL;
        newWrep->edge[iedge].eprite = NULL;
        newWrep->edge[iedge].eedge  = eedges[iedge-1];

        for (inode = 1; inode <= newWrep->nnode; inode++) {
            if (echilds[0] == newWrep->node[inode].enode) {
                newWrep->edge[iedge].ibeg = inode;
                newWrep->node[inode].nedge++;
            }
            if (echilds[1] == newWrep->node[inode].enode) {
                newWrep->edge[iedge].iend = inode;
                newWrep->node[inode].nedge++;
            }
        }
    }

    /* fill in the Wface data structure */
    for (iface = 1; iface <= newWrep->nface; iface++) {
        SPLINT_CHECK_FOR_NULL(efaces);

        status = EG_getTopology(efaces[iface-1], &esurface, &oclass, &mtype,
                                data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        newWrep->face[iface].tag      = 0;
        newWrep->face[iface].mtype    = mtype;
        newWrep->face[iface].esurface = esurface;
        newWrep->face[iface].eface    = efaces[iface-1];

        for (iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            for (ichild = 0; ichild < nchild; ichild++) {
                iedge = status = EG_indexBodyTopo(ebody, echilds[ichild]);
                CHECK_STATUS(EG_indexBodyTopo);

                jedge = status = EG_indexBodyTopo(ebody, echilds[(ichild-1+nchild)%nchild]);
                CHECK_STATUS(EG_indexBodyTopo);

                kedge = status = EG_indexBodyTopo(ebody, echilds[(ichild+1+nchild)%nchild]);
                CHECK_STATUS(EG_indexBodyTopo);

                if (senses[ichild] == SFORWARD) {
                    newWrep->edge[iedge].ileft  = iface;
                    newWrep->edge[iedge].ibleft = jedge;
                    newWrep->edge[iedge].ieleft = kedge;
                    if (eref != NULL) {
                        newWrep->edge[iedge].epleft = echilds[ichild+nchild];
                    }
                } else {
                    newWrep->edge[iedge].irite  = iface;
                    newWrep->edge[iedge].ibrite = kedge;
                    newWrep->edge[iedge].ierite = jedge;
                    if (eref != NULL) {
                        newWrep->edge[iedge].eprite = echilds[ichild+nchild];
                    }
                }

            }
        }
    }

    newWrep->nface_orig = newWrep->nface;

    *wrep = newWrep;

cleanup:
    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepBreakEdge - break a given Wedge at the given t                 *
 *                                                                      *
 ************************************************************************
 */

static int
wrepBreakEdge(wrep_T  *wrep,            /* (in)  pointer to a WREP structure */
              int     iedge,            /* (in)  Wedge to break */
              double  t,                /* (in)  parametric coordinate at split */
              int     *jnode)           /* (out) Wnode at the break */
{
    int     status = EGADS_SUCCESS;     /* (out) <0 return status */
                                        /* (out) >0 is Wedge index */

    int     jedge, kedge, oclass, mtype;
    double  data[18];
    ego     context, enodes[2], topRef, prev, next;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(wrepBreakEdge);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter wrepBreakEdge(iedge=%d, t=%10.5f)\n", iedge, t);
#endif

    /* try to adjust t to put it between tbeg and tend if a CIRCLE */
    status = EG_getInfo(wrep->edge[iedge].ecurve, &oclass, &mtype, &topRef, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    if (mtype == CIRCLE) {
        if        (t-TWOPI >= wrep->edge[iedge].tbeg-EPS06 && t-TWOPI <= wrep->edge[iedge].tend+EPS06) {
            t -= TWOPI;
        } else if (t+TWOPI >= wrep->edge[iedge].tbeg-EPS06 && t+TWOPI <= wrep->edge[iedge].tend+EPS06) {
            t += TWOPI;
        }
    }

    /* return immediately if break is at end of Wedge */
    if        (fabs(t-wrep->edge[iedge].tbeg) < EPS06) {
        *jnode = wrep->edge[iedge].ibeg;
        status = iedge;
        goto cleanup;
    } else if (fabs(t-wrep->edge[iedge].tend) < EPS06) {
        *jnode = wrep->edge[iedge].iend;
        status = iedge;
        goto cleanup;
    }

    /* if t < tbeg, see if there is an adjoining Wedge that shares a curve with iedge */
    if (t < wrep->edge[iedge].tbeg) {
        if        (wrep->edge[wrep->edge[iedge].ibleft].ecurve == wrep->edge[iedge].ecurve) {
            iedge = wrep->edge[iedge].ibleft;
        } else if (wrep->edge[wrep->edge[iedge].ibrite].ecurve == wrep->edge[iedge].ecurve) {
            iedge = wrep->edge[iedge].ibrite;
        } else {
            snprintf(message, 1024, "cannot find adjoining Wedge that shares a curve with iedge=%d", iedge);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }
    }

    /* if t > tend, see if there is an adjoinig Wedge that shares a curve with iedge */
    if (t > wrep->edge[iedge].tend) {
        if        (wrep->edge[wrep->edge[iedge].ieleft].ecurve == wrep->edge[iedge].ecurve) {
            iedge = wrep->edge[iedge].ieleft;
        } else if (wrep->edge[wrep->edge[iedge].ierite].ecurve == wrep->edge[iedge].ecurve) {
            iedge = wrep->edge[iedge].ierite;
        } else {
            snprintf(message, 1024, "cannot find adjoining Wedge that shares a curve with iedge=%d", iedge);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }
    }

    status = EG_getContext(wrep->edge[iedge].eedge, &context);
    CHECK_STATUS(EG_getContext);

    /* get the coordinates of the new Wnode */
    status = EG_evaluate(wrep->edge[iedge].eedge, &t, data);
    CHECK_STATUS(EG_evaluate);

    /* add a Wnode */
    RALLOC(wrep->node, wnode_T, wrep->nnode+2);

    *jnode = wrep->nnode + 1;

    wrep->node[*jnode].xyz[0] = data[0];
    wrep->node[*jnode].xyz[1] = data[1];
    wrep->node[*jnode].xyz[2] = data[2];
    wrep->node[*jnode].nedge  = 2;
    wrep->node[*jnode].enode  = NULL;

    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(wrep->node[*jnode].enode));
    CHECK_STATUS(EG_makeTopology);

    (wrep->nnode)++;

    /* add a Wedge (the second part of iedge) */
    RALLOC(wrep->edge, wedge_T, wrep->nedge+2);

    jedge = wrep->nedge + 1;

    wrep->edge[jedge].tag    = wrep->edge[iedge].tag;
    wrep->edge[jedge].ibeg   = *jnode;
    wrep->edge[jedge].iend   = wrep->edge[iedge].iend;
    wrep->edge[jedge].ileft  = wrep->edge[iedge].ileft;
    wrep->edge[jedge].irite  = wrep->edge[iedge].irite;
    wrep->edge[jedge].ibleft = iedge;
    wrep->edge[jedge].ibrite = iedge;
    wrep->edge[jedge].ieleft = wrep->edge[iedge].ieleft;
    wrep->edge[jedge].ierite = wrep->edge[iedge].ierite;
    wrep->edge[jedge].tbeg   = t;
    wrep->edge[jedge].tend   = wrep->edge[iedge].tend;
    wrep->edge[jedge].ecurve = wrep->edge[iedge].ecurve;
    wrep->edge[jedge].epleft = wrep->edge[iedge].epleft;
    wrep->edge[jedge].eprite = wrep->edge[iedge].eprite;
    wrep->edge[jedge].eedge  = wrep->edge[iedge].eedge;

    data[0] = wrep->edge[jedge].tbeg;
    data[1] = wrep->edge[jedge].tend;

    enodes[0] = wrep->node[wrep->edge[jedge].ibeg].enode;
    enodes[1] = wrep->node[wrep->edge[jedge].iend].enode;

    status = EG_makeTopology(context, wrep->edge[jedge].ecurve, EDGE, TWONODE,
                             data, 2, enodes, NULL, &(wrep->edge[jedge].eedge));
    CHECK_STATUS(EG_makeTopology);

    (wrep->nedge)++;

    /* modify iedge (to be the first part of the original iedge) */
    wrep->edge[iedge].iend   = *jnode;
    wrep->edge[iedge].ieleft = jedge;
    wrep->edge[iedge].ierite = jedge;
    wrep->edge[iedge].tend   = t;

    data[0] = wrep->edge[iedge].tbeg;
    data[1] = wrep->edge[iedge].tend;

    enodes[0] = wrep->node[wrep->edge[iedge].ibeg].enode;
    enodes[1] = wrep->node[wrep->edge[iedge].iend].enode;

    status = EG_makeTopology(context, wrep->edge[iedge].ecurve, EDGE, TWONODE,
                             data, 2, enodes, NULL, &(wrep->edge[iedge].eedge));
    CHECK_STATUS(EG_makeTopology);

    /* modify Wedges that were incident on the end of iedge */
    kedge = wrep->edge[jedge].ieleft;
    if (wrep->edge[kedge].ibleft == iedge) {
        wrep->edge[kedge].ibleft =  jedge;
    }
    if (wrep->edge[kedge].ierite == iedge) {
        wrep->edge[kedge].ierite =  jedge;
    }

    kedge = wrep->edge[jedge].ierite;
    if (wrep->edge[kedge].ieleft == iedge) {
        wrep->edge[kedge].ieleft =  jedge;
    }
    if (wrep->edge[kedge].ibrite == iedge) {
        wrep->edge[kedge].ibrite =  jedge;
    }

    status = jedge;

cleanup:
#ifdef DEBUG
    printf("exit wrepBreakEdge -> status=%d, jnode=%d\n", status, *jnode);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepMakeEdge - build a new Wedge between Wnodes                    *
 *                                                                      *
 ************************************************************************
 */

static int
wrepMakeEdge(wrep_T  *wrep,             /* (in)  pointer to a WREP structure */
             int     iface,             /* (in)  Wface on which Wedge will reside */
             int     ibeg,              /* (in)  Wnode at begining of Wedge */
             int     iend,              /* (in)  Wnode at end      of Wedge */
   /*@null@*/ego     ecurvei)           /* (in)  ecurve for new Wedge */
{
    int     status = EGADS_SUCCESS;     /* (out) <0 return status */
                                        /* (out) >0 is Wedge index */

    int     oclass, mtype, knode, jedge, kedge, jface, kface, used;
    double  data[18], tout[2], uvout[4], xyzout[3];
    ego     ecurve, context, topRef, prev, next, epline, epleft, eprite;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(wrepMakeEdge);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter wrepMakeEdge(iface=%d, ibeg=%d, iend=%d, ecurvei=%llx)\n", iface, ibeg, iend, (long long)ecurvei);
#endif

    /* if ecurve is not given, create a new curve (straight trimmed pcurve) between ibeg and iend */
    if (ecurvei == NULL) {
        status = EG_getContext(wrep->face[iface].esurface, &context);
        CHECK_STATUS(EG_getContext);

        status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[ibeg].xyz, &uvout[0], xyzout);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(wrep->face[iface].esurface, wrep->node[iend].xyz, &uvout[2], xyzout);
        CHECK_STATUS(EG_invEvaluate);

        uvout[2] -= uvout[0];
        uvout[3] -= uvout[1];

        status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, uvout, &epline);
        CHECK_STATUS(EG_makeGeometry);

        data[0] = uvout[0];
        data[1] = uvout[1];
        status = EG_invEvaluate(epline, data, &tout[0], xyzout);
        CHECK_STATUS(EG_invEvaluate);

        data[0] = uvout[2] + uvout[0];
        data[1] = uvout[3] + uvout[1];
        status = EG_invEvaluate(epline, data, &tout[1], xyzout);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeGeometry(context, PCURVE, TRIMMED, epline, NULL, tout, &epleft);
        CHECK_STATUS(EG_makeGeometry);

        status = EG_otherCurve(wrep->face[iface].esurface, epleft, 0.0, &ecurve);
        CHECK_STATUS(EG_otherCurve);

        status = EG_deleteObject(epleft);
        CHECK_STATUS(EG_deleteObject);
    } else {
        ecurve = ecurvei;
    }

    status = EG_getInfo(wrep->face[iface].esurface, &oclass, &mtype, &topRef, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    SPLINT_CHECK_FOR_NULL(ecurve);

    if (mtype == PLANE) {
        epleft = NULL;
        eprite = NULL;
    } else {
        status = EG_otherCurve(wrep->face[iface].esurface, ecurve, 0.0, &epleft);
        if (status == EGADS_CONSTERR) {
            snprintf(message, 1024, "perhaps the offset is not contained to the Faces in facelist");
            status = OCSM_UDP_ERROR2;
            goto cleanup;
        }
        CHECK_STATUS(EG_otherCurve);

        eprite = epleft;
    }

    /* create a new Wedge between ibeg and iend */
    RALLOC(wrep->edge, wedge_T, wrep->nedge+2);

    jedge = wrep->nedge + 1;

    wrep->edge[jedge].tag    = 0;
    wrep->edge[jedge].ibeg   = ibeg;
    wrep->edge[jedge].iend   = iend;
    wrep->edge[jedge].ileft  = 0;                 // corrected below
    wrep->edge[jedge].irite  = 0;                 // corrected below
    wrep->edge[jedge].ibleft = 0;                 // corrected below
    wrep->edge[jedge].ibrite = 0;                 // corrected below
    wrep->edge[jedge].ieleft = 0;                 // corrected below
    wrep->edge[jedge].ierite = 0;                 // corrected below
    wrep->edge[jedge].tbeg   = 0;                 // corrected below
    wrep->edge[jedge].tend   = 0;                 // corrected below
    wrep->edge[jedge].ecurve = ecurve;
    /*@-kepttrans@*/
    wrep->edge[jedge].epleft = epleft;
    wrep->edge[jedge].eprite = eprite;
    /*@+kepttrans@*/
    wrep->edge[jedge].eedge  = NULL;

    status = EG_invEvaluate(ecurve, wrep->node[ibeg].xyz, tout, xyzout);
    CHECK_STATUS(EG_invEvaluate);
    wrep->edge[jedge].tbeg = tout[0];

    status = EG_invEvaluate(ecurve, wrep->node[iend].xyz, tout, xyzout);
    CHECK_STATUS(EG_invEvaluate);
    wrep->edge[jedge].tend = tout[0];

    (wrep->nedge)++;

    /* attaching to an isolated Wnode */
    if        (wrep->node[ibeg].nedge < 1) {
        wrep->edge[jedge].ileft = iface;
        wrep->edge[jedge].irite = iface;
    } else if (wrep->node[iend].nedge < 1) {
        wrep->edge[jedge].ileft = iface;
        wrep->edge[jedge].irite = iface;

    /* if ibeg and iend have a common Wface, set the ileft and irite of the new Wedge.
       note that iface is favored */
    } else {
        for (kface = 1; kface <= wrep->nface; kface++) {
            used = 0;
            for (kedge = 1; kedge < jedge; kedge++) {
                if (wrep->edge[kedge].ibeg == ibeg || wrep->edge[kedge].iend == ibeg) {
                    if (wrep->edge[kedge].ileft == kface || wrep->edge[kedge].irite == kface) {
                        used = 1;
                        break;
                    }
                }
            }
            if (used == 0) continue;

            used = 0;
            for (kedge = 1; kedge < jedge; kedge++) {
                if (wrep->edge[kedge].ibeg == iend || wrep->edge[kedge].iend == iend) {
                    if (wrep->edge[kedge].ileft == kface || wrep->edge[kedge].irite == kface) {
                        used = 1;
                        break;
                    }
                }
            }
            if (used == 0) continue;

            if (kface == iface) {
                wrep->edge[jedge].ileft = kface;
                wrep->edge[jedge].irite = kface;
                break;
            } else if (wrep->edge[jedge].ileft == 0 && wrep->edge[jedge].irite == 0) {
                wrep->edge[jedge].ileft = kface;
                wrep->edge[jedge].irite = kface;
            }
        }
    }

    if (wrep->edge[jedge].ileft == 0 || wrep->edge[jedge].irite == 0) {
        wrepPrint(wrep);
        snprintf(message, 1024, "eithe rileft or irite of jedge=%d is zero", jedge);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* determine the ileft and irite for the new Wedge and change
       the links at ibeg and iend to point to the new Wedge */
    for (kedge = 1; kedge < jedge; kedge++) {

        /* attaching to kedge, which has a Wnode with .nedge==1

           A <--             --> C
                 ---------->
           B -->    jedge    <-- D

         */

        /* case A */
        if (wrep->edge[kedge].ibeg == ibeg && wrep->edge[kedge].ibleft == 0
                                           && wrep->edge[kedge].ibrite == 0   ) {
#ifdef DEBUG
            printf("case A: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ibleft = jedge;
            wrep->edge[kedge].ibrite = jedge;
            wrep->edge[jedge].ibleft = kedge;
            wrep->edge[jedge].ibrite = kedge;
        }

        /* case B */
        if (wrep->edge[kedge].iend == ibeg && wrep->edge[kedge].ieleft == 0
                                           && wrep->edge[kedge].ierite == 0   ) {
#ifdef DEBUG
            printf("case B: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ieleft = jedge;
            wrep->edge[kedge].ierite = jedge;
            wrep->edge[jedge].ibleft = kedge;
            wrep->edge[jedge].ibrite = kedge;
        }

        /* case C */
        if (wrep->edge[kedge].ibeg == iend && wrep->edge[kedge].ibleft == 0
                                           && wrep->edge[kedge].ibrite == 0   ) {
#ifdef DEBUG
            printf("case C: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ibleft = jedge;
            wrep->edge[kedge].ibrite = jedge;
            wrep->edge[jedge].ieleft = kedge;
            wrep->edge[jedge].ierite = kedge;
        }

        /* case D */
        if (wrep->edge[kedge].iend == iend && wrep->edge[kedge].ieleft == 0
                                           && wrep->edge[kedge].ierite == 0   ) {
#ifdef DEBUG
            printf("case D: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ieleft = jedge;
            wrep->edge[kedge].ierite = jedge;
            wrep->edge[jedge].ieleft = kedge;
            wrep->edge[jedge].ierite = kedge;
        }

        /* attaching to kedge, which has a Wnode with .nedge>1

                 G E                 I K
                  \ ^               ^ /
                   v \    jedge    / v
                      ------------>
                   / ^             ^ \
                  v /               \ v
                 F H                 L J
        */

        /* case E */
        if (wrep->edge[kedge].ibeg == ibeg && wrep->edge[kedge].irite == wrep->edge[jedge].ileft) {
#ifdef DEBUG
            printf("case E: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].irite;

            wrep->edge[kedge].ibrite = jedge;
            wrep->edge[jedge].ibleft = kedge;
        }

        /* case F */
        if (wrep->edge[kedge].ibeg == ibeg && wrep->edge[kedge].ileft == wrep->edge[jedge].irite) {
#ifdef DEBUG
            printf("case F: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].irite = wrep->edge[kedge].ileft;

            wrep->edge[kedge].ibleft = jedge;
            wrep->edge[jedge].ibrite = kedge;
        }

        /* case G */
        if (wrep->edge[kedge].iend == ibeg && wrep->edge[kedge].ileft == wrep->edge[jedge].ileft) {
#ifdef DEBUG
            printf("case G: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;

            wrep->edge[kedge].ieleft = jedge;
            wrep->edge[jedge].ibleft = kedge;
        }

        /* case H */
        if (wrep->edge[kedge].iend == ibeg && wrep->edge[kedge].irite == wrep->edge[jedge].irite) {
#ifdef DEBUG
            printf("case H: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ierite = jedge;
            wrep->edge[jedge].ibrite = kedge;
        }

        /* case I */
        if (wrep->edge[kedge].ibeg == iend && wrep->edge[kedge].ileft == wrep->edge[jedge].ileft) {
#ifdef DEBUG
            printf("case I: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].ileft;

            wrep->edge[kedge].ibleft = jedge;
            wrep->edge[jedge].ieleft = kedge;
        }

        /* case J */
        if (wrep->edge[kedge].ibeg == iend && wrep->edge[kedge].irite == wrep->edge[jedge].irite) {
#ifdef DEBUG
            printf("case J: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].irite = wrep->edge[kedge].irite;

            wrep->edge[kedge].ibrite = jedge;
            wrep->edge[jedge].ierite = kedge;
        }

        /* case K */
        if (wrep->edge[kedge].iend == iend && wrep->edge[kedge].irite == wrep->edge[jedge].ileft) {
#ifdef DEBUG
            printf("case K: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].ileft = wrep->edge[kedge].irite;

            wrep->edge[kedge].ierite = jedge;
            wrep->edge[jedge].ieleft = kedge;
        }

        /* case L */
        if (wrep->edge[kedge].iend == iend && wrep->edge[kedge].ileft == wrep->edge[jedge].irite) {
#ifdef DEBUG
            printf("case L: jedge=%d, kedge=%d\n", jedge, kedge);
#endif
            wrep->edge[jedge].irite = wrep->edge[kedge].ileft;

            wrep->edge[kedge].ieleft = jedge;
            wrep->edge[jedge].ierite = kedge;
        }
    }

    /* if the number of incident Wedges on both ibeg and iend are greater than
       zero, this means we have closed a loop, so create a new Wface on the right of the new Wedge */
    if (wrep->node[ibeg].nedge > 0 && wrep->node[iend].nedge > 0) {
        RALLOC(wrep->face, wface_T, wrep->nface+2);

        jface = wrep->nface + 1;
#ifdef DEBUG
        printf("creating Wface %d\n", jface);
#endif

        wrep->face[jface].tag      = 0;
        wrep->face[jface].mtype    = wrep->face[iface].mtype;
        wrep->face[jface].esurface = wrep->face[iface].esurface;
        wrep->face[jface].eface    = wrep->face[iface].eface;

        (wrep->nface)++;

        /* change the ileft and irite pointers on all Wedges associated with the new Wface */
        wrep->edge[jedge].irite = jface;

        knode = wrep->edge[jedge].iend;
        kedge = wrep->edge[jedge].ierite;
        while (kedge != jedge) {
            if (kedge == 0) {
                wrepPrint(wrep);
                snprintf(message, 1024, "trouble setting ileft/irite for kedge=%d", kedge);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            }

            if        (wrep->edge[kedge].ieleft == 0) {
                wrep->edge[kedge].ileft = jface;
                wrep->edge[kedge].irite = jface;
                kedge = wrep->edge[kedge].ibleft;
            } else if (wrep->edge[kedge].ibeg == knode) {
                wrep->edge[kedge].irite = jface;
                knode = wrep->edge[kedge].iend;
                kedge = wrep->edge[kedge].ierite;
            } else if (wrep->edge[kedge].iend == knode) {
                wrep->edge[kedge].ileft = jface;
                knode = wrep->edge[kedge].ibeg;
                kedge = wrep->edge[kedge].ibleft;
            } else {
                wrepPrint(wrep);
                snprintf(message, 1024, "trouble setting ileft/irite for kedge=%d", kedge);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            }
        }
    }

    /* increment the number of Wedges incident on the two end-points */
    wrep->node[ibeg].nedge++;
    wrep->node[iend].nedge++;

    /* make sure all info is updated */
    if (wrep->edge[jedge].ileft  == 0 ||
        wrep->edge[jedge].irite  == 0 ||
        wrep->edge[jedge].ibleft == 0 ||
        wrep->edge[jedge].ibrite == 0   ) {
        printf("jedge=%d, ileft=%d, irite=%d, ibleft=%d, ibrite=%d\n", jedge,
               wrep->edge[jedge].ileft,
               wrep->edge[jedge].irite,
               wrep->edge[jedge].ibleft,
               wrep->edge[jedge].ibrite);
        snprintf(message, 1024, "info not updated for jedge=%d", jedge);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    status = jedge;

cleanup:
#ifdef DEBUG
    wrepPrint(wrep);
    printf("exit wrepMakeEdge -> status=%d\n", status);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepMakeNode - build a new Wnode                                   *
 *                                                                      *
 ************************************************************************
 */

static int
wrepMakeNode(wrep_T  *wrep,             /* (in)  pointer to a WREP structure */
             int     iface,             /* (in)  Wface on which Wnode should lie */
             double  xyz[])             /* (in)  coordinate of new Wnode */
{
    int     status = EGADS_SUCCESS;     /* (out) <0 return status */
                                        /* (out) >0 is Wnode index */

    int     jnode;
    double  uv_close[2], xyz_close[3];
    ego     context;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(wrepMakeNode);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter wrepMakeNode(iface=%d, xyz=%10.5f, %10.5f %10.5f)\n",
           iface, xyz[0], xyz[1], xyz[2]);
#endif

    /* check if Wnode already exists */
    for (jnode = 1; jnode <= wrep->nnode; jnode++) {
        if (fabs(wrep->node[jnode].xyz[0] - xyz[0]) < EPS06 &&
            fabs(wrep->node[jnode].xyz[1] - xyz[1]) < EPS06 &&
            fabs(wrep->node[jnode].xyz[2] - xyz[2]) < EPS06   ) {
            status = jnode;
            goto cleanup;
        }
    }

    status = EG_getContext(wrep->face[iface].esurface, &context);
    CHECK_STATUS(EG_getContext);

    /* adjust xyz so that it actually lies on the esurface */
    status = EG_invEvaluate(wrep->face[iface].esurface, xyz, uv_close, xyz_close);
    CHECK_STATUS(EG_invEvaluate);

    /* add a Wnode */
    RALLOC(wrep->node, wnode_T, wrep->nnode+2);

    jnode = wrep->nnode + 1;

    wrep->node[jnode].xyz[0] = xyz_close[0];
    wrep->node[jnode].xyz[1] = xyz_close[1];
    wrep->node[jnode].xyz[2] = xyz_close[2];
    wrep->node[jnode].nedge  = 0;
    wrep->node[jnode].enode  = NULL;

    status = EG_makeTopology(context, NULL, NODE, 0,
                             xyz_close, 0, NULL, NULL, &(wrep->node[jnode].enode));
    CHECK_STATUS(EG_makeTopology);

    (wrep->nnode)++;

    status = jnode;

cleanup:
#ifdef DEBUG
    printf("exit wrepMakeNode -> status=%d\n", status);
#endif

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrep2ego - convert a wrep into an EGADS Body                       *
 *                                                                      *
 ************************************************************************
 */

static int
wrep2ego(wrep_T  wrep,                  /* (in)  pointer to a WREP structure */
         ego     context,               /* (in)  EGADS context for Body */
         ego     *ebody)                /* (out) EGADS Body that was created */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     inode, iedge, nedge, iface, nloop, iloop, jloop, ibeg, oclass, mtype;
    int     *senses=NULL, *iedges=NULL, *iused=NULL;
    double  data[18];
    ego     enodes[2], *eedges=NULL, *eloops=NULL, *efaces=NULL, eshell;
    ego     oldEdge, topRef, prev, next;

    ROUTINE(wrep2ego);

    /* --------------------------------------------------------------- */

    *ebody = NULL;

    /* build the Nodes */
    for (inode = 1; inode <= wrep.nnode; inode++) {
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 wrep.node[inode].xyz, 0, NULL, NULL, &(wrep.node[inode].enode));
        CHECK_STATUS(EG_makeTopology);
    }

    /* build the Edges */
    for (iedge = 1; iedge <= wrep.nedge; iedge++) {
#ifdef DEBUG
        printf("building Edge %3d\n", iedge);
#endif
        oldEdge = wrep.edge[iedge].eedge;

        /* now make the Edge */
        data[0] = wrep.edge[iedge].tbeg;
        data[1] = wrep.edge[iedge].tend;

        enodes[0] = wrep.node[wrep.edge[iedge].ibeg].enode;
        enodes[1] = wrep.node[wrep.edge[iedge].iend].enode;

        if (wrep.edge[iedge].ecurve != NULL) {
            status = EG_makeTopology(context, wrep.edge[iedge].ecurve, EDGE, TWONODE,
                                     data, 2, enodes, NULL, &(wrep.edge[iedge].eedge));
            CHECK_STATUS(EG_makeTopology);
        } else {
            status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                     data, 1, enodes, NULL, &(wrep.edge[iedge].eedge));
            CHECK_STATUS(EG_makeTopology);
        }

        if (oldEdge != NULL) {
            status = EG_attributeDup(oldEdge, wrep.edge[iedge].eedge);
            CHECK_STATUS(EG_attributeDup);
        }
    }

    /* allocate temporary memory (to hold egos) */
    MALLOC(eloops, ego,   wrep.nedge+1);     // perhaps too big
    MALLOC(eedges, ego, 2*wrep.nedge+1);     // perhaps too big
    MALLOC(senses, int, 2*wrep.nedge+1);     // perhaps too big
    MALLOC(iedges, int,   wrep.nedge+1);     // perhaps too big
    MALLOC(efaces, ego,   wrep.nface+1);
    MALLOC(iused,  int,   wrep.nedge+1);

    /* build the Faces */
    for (iface = 1; iface <= wrep.nface; iface++) {
#ifdef DEBUG
        printf("building Face %3d\n", iface);
#endif
        nloop = 0;

        for (iedge = 1; iedge <= wrep.nedge; iedge++) {
            iused[iedge] = 0;
        }

        /* build the Loops in this Face */
        for (iloop = 0; iloop < wrep.nedge; iloop++) {
#ifdef DEBUG
            printf("    adding Loop %3d\n", iloop+1);
#endif
            ibeg      = 0;
            nedge     = 0;
            senses[0] = 0;

            /* find a non-degenerate Edge associated with this Face */
            for (iedge = 1; iedge <= wrep.nedge; iedge++) {
                if (iused[iedge] != 0) continue;
                if (wrep.edge[iedge].ibeg == wrep.edge[iedge].iend) continue;

                if        (wrep.edge[iedge].ileft == iface) {
#ifdef DEBUG
                    printf("        adding Edge %3d\n", iedge);
#endif
                    iedges[nedge] = iedge;
                    eedges[nedge] = wrep.edge[iedge].eedge;
                    senses[nedge] = SFORWARD;
                    nedge++;

                    iused[iedge] = 1;

                    ibeg  = iedge;
                    iedge = wrep.edge[iedge].ieleft;
                    break;
                } else if (wrep.edge[iedge].irite == iface) {
#ifdef DEBUG
                    printf("        adding Edge %3d\n", iedge);
#endif
                    iedges[nedge] = iedge;
                    eedges[nedge] = wrep.edge[iedge].eedge;
                    senses[nedge] = SREVERSE;
                    nedge++;

                    iused[iedge] = 1;

                    ibeg  = iedge;
                    iedge = wrep.edge[iedge].ibrite;
                    break;
                }
            }

            if (senses[0] == 0) break;            /* no Edge found */
            if (ibeg      == 0) break;

            /* keep adding Edges to this Loop until we get back to the beginning */
            while (iedge != ibeg) {
                if (wrep.edge[iedge].ileft == iface) {
#ifdef DEBUG
                    printf("        adding Edge %3d\n", iedge);
#endif
                    iedges[nedge] = iedge;
                    eedges[nedge] = wrep.edge[iedge].eedge;
                    senses[nedge] = SFORWARD;
                    nedge++;

                    iused[iedge] = 1;

                    iedge = wrep.edge[iedge].ieleft;
                } else if (wrep.edge[iedge].irite == iface) {
#ifdef DEBUG
                    printf("        adding Edge %3d\n", iedge);
#endif
                    iedges[nedge] = iedge;
                    eedges[nedge] = wrep.edge[iedge].eedge;
                    senses[nedge] = SREVERSE;
                    nedge++;

                    iused[iedge] = 1;

                    iedge = wrep.edge[iedge].ibrite;
                } else {
                    snprintf(message, 1024, "having trouble traversing the Loop");
                    status = OCSM_INTERNAL_ERROR;
                    goto cleanup;
                }
            }

            /* we should add pcurves here */
            status = EG_getInfo(wrep.face[iface].esurface, &oclass, &mtype, &topRef, &prev, &next);
            CHECK_STATUS(EG_getInfo);

            if (mtype != PLANE) {
                for (iedge = 0; iedge < nedge; iedge++) {
                    if (senses[iedge] == SFORWARD) {
                        eedges[iedge+nedge] = wrep.edge[iedges[iedge]].epleft;
                    } else {
                        eedges[iedge+nedge] = wrep.edge[iedges[iedge]].eprite;
                    }

                    if (eedges[iedge+nedge] == NULL) {
                        status = EG_otherCurve(wrep.face[iface].esurface, eedges[iedge], 0.0, &eedges[iedge+nedge]);
                        CHECK_STATUS(EG_otherCurve);
                    }
                }
            }

            /* make the Loop */
            if (mtype != PLANE) {
                status = EG_makeTopology(context, wrep.face[iface].esurface, LOOP, CLOSED,
                                         NULL, nedge, eedges, senses, &eloops[nloop]);
                CHECK_STATUS(EG_makeTopology);
            } else {
                status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                         NULL, nedge, eedges, senses, &eloops[nloop]);
                CHECK_STATUS(EG_makeTopology);
            }
            nloop++;
        }

        /* make the Face */
        senses[0] = SFORWARD;
        for (iloop = 1; iloop < nloop; iloop++) {
            senses[iloop] = SREVERSE;
        }
        status = EG_makeTopology(context, wrep.face[iface].esurface, FACE, SFORWARD,
                                 NULL, nloop, eloops, senses, &efaces[iface-1]);

        /* if we got a construction error, it might be because we have the wrong
           Loop identified as the outer Loop */
        for (iloop = 1; iloop < nloop; iloop++) {
            if (status != EGADS_CONSTERR) break;

            prev = eloops[0];
            for (jloop = 0; jloop < nloop-1; jloop++) {
                eloops[jloop] = eloops[jloop+1];
            }
            eloops[nloop-1] = prev;

            status = EG_makeTopology(context, wrep.face[iface].esurface, FACE, SFORWARD,
                                     NULL, nloop, eloops, senses, &efaces[iface-1]);
        }
        CHECK_STATUS(EG_makeTopology);

        /* copy the Attributes from the original Faces onto the new Face */
        status = EG_attributeDup(wrep.face[iface].eface, efaces[iface-1]);
        CHECK_STATUS(EG_attributeDup);

        /* put the "scribeFace" attribute on Faces that are within
           the scribe */
        if (wrep.face[iface].new == 1) {
            status = EG_attributeAdd(efaces[iface-1], "__offsetFace__", ATTRINT, 1,
                                     &iface, NULL, NULL);
            CHECK_STATUS(EG_attributeAdd);
        }
    }

    /* build a Shell */
#ifdef DEBUG
    printf("building Shell\n");
#endif
    status = EG_makeTopology(context, NULL, SHELL, CLOSED,
                             NULL, wrep.nface, efaces, NULL, &eshell);
    CHECK_STATUS(EG_makeTopology);

    /* build the new Body */
#ifdef DEBUG
    printf("building SolidBody\n");
#endif
    status = EG_makeTopology(context, NULL, BODY, SOLIDBODY,
                             NULL, 1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    FREE(eloops);
    FREE(eedges);
    FREE(senses);
    FREE(iedges);
    FREE(efaces);
    FREE(iused );

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepPrint - print a WREP structure                                 *
 *                                                                      *
 ************************************************************************
 */

static int
wrepPrint(wrep_T *wrep)                 /* (in)  pointer to a WREP structure */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    int     inode, iedge, iface;

    ROUTINE(wrepPrint);

    /* --------------------------------------------------------------- */

    printf("wrep associated with ebody=%12llx\n", (long long)wrep->ebody);

    printf("inode        x               y               z         nedg    dist           enode\n");
    for (inode = 1; inode <= wrep->nnode; inode++) {
        printf("%5d %15.8f %15.8f %15.8f %5d %10.7f %12llx\n", inode,
               wrep->node[inode].xyz[0],
               wrep->node[inode].xyz[1],
               wrep->node[inode].xyz[2],
               wrep->node[inode].nedge,
               wrep->node[inode].dist,
    (long long)wrep->node[inode].enode);
    }

    printf("iedge   tag  ibeg  iend ileft irite iblft ibrit ielft ierit        tbeg            tend           ecurve       epleft       eprite        eedge\n");
    for (iedge = 1; iedge <= wrep->nedge; iedge++) {
        printf("%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %15.8f %15.8f %12llx %12llx %12llx %12llx\n", iedge,
               wrep->edge[iedge].tag,
               wrep->edge[iedge].ibeg,
               wrep->edge[iedge].iend,
               wrep->edge[iedge].ileft,
               wrep->edge[iedge].irite,
               wrep->edge[iedge].ibleft,
               wrep->edge[iedge].ibrite,
               wrep->edge[iedge].ieleft,
               wrep->edge[iedge].ierite,
               wrep->edge[iedge].tbeg,
               wrep->edge[iedge].tend,
    (long long)wrep->edge[iedge].ecurve,
    (long long)wrep->edge[iedge].epleft,
    (long long)wrep->edge[iedge].eprite,
    (long long)wrep->edge[iedge].eedge);
    }

    printf("iface   tag mtype    esurface        eface\n");
    for (iface = 1; iface <= wrep->nface; iface++) {
        printf("%5d %5d %5d %12llx %12llx\n", iface,
               wrep->face[iface].tag,
               wrep->face[iface].mtype,
    (long long)wrep->face[iface].esurface,
    (long long)wrep->face[iface].eface);
    }

//cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   wrepFree - free a WREP structure                                   *
 *                                                                      *
 ************************************************************************
 */

static int
wrepFree(wrep_T   *wrep)                /* (in)  pointer to a WREP structure */
{
    int     status = EGADS_SUCCESS;     /* (out) return status */

    ROUTINE(wrepFree);

    /* --------------------------------------------------------------- */

    if (wrep != NULL) {
        FREE(wrep->node);
        FREE(wrep->edge);
        FREE(wrep->face);

        wrep->nnode = 0;
        wrep->nedge = 0;
        wrep->nface = 0;
    }

//cleanup:
    return status;
}
