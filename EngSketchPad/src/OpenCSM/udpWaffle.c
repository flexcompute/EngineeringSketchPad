/*
 ************************************************************************
 *                                                                      *
 * udpWaffle -- udp file to generate 2D waffle                          *
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

#define NUMUDPARGS 6
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define DEPTH(       IUDP)    ((double *) (udps[IUDP].arg[0].val))[0]
#define DEPTH_DOT(   IUDP)    ((double *) (udps[IUDP].arg[0].dot))[0]
#define SEGMENTS(    IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define SEGMENTS_DOT(IUDP,I)  ((double *) (udps[IUDP].arg[1].dot))[I]
#define FILENAME(    IUDP)    ((char   *) (udps[IUDP].arg[2].val))
#define PROGRESS(    IUDP)    ((int    *) (udps[IUDP].arg[3].val))[0]
#define LAYOUT(      IUDP)    ((int    *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"depth",     "segments",  "filename", "progress", "layout",  "rebuild",   };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRFILE,   ATTRINT,    ATTRINT,   ATTRREBUILD, };
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,          0,          0,         0,           };
static double argDdefs[NUMUDPARGS] = {1.,          0.,          0.,         0.,         0.,        0.,          };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#ifdef GRAFIC
    #include "grafic.h"
    #define DEBUG
#endif

/* include OpenCSM to make use of its evaluator */
#include "OpenCSM.h"

typedef struct {
    int    type;                        /* =0 for cpoint, =1 for point */
    double x;                           /*             X-coordinate */
    double x_dot;                       /* velocity of X-coordinate */
    double y;                           /*             Y-coordinate */
    double y_dot;                       /* velocity of Y-coordinate */
    char   name[80];                    /* name */
} pnt_T;

typedef struct {
    int    type;                        /* =0 for cline, =1 for line */
    int    ibeg;                        /* Point at beginning */
    int    iend;                        /* Point at end */
    int    num;                         /* number */
    int    idx;                         /* index */
    char   name[80];                    /* name */
    int    nattr;                       /* number of Attributes */
    char   **aname;                     /* array  of Attribute names */
    char   **avalu;                     /* array  of Attribute values */
} seg_T;

/* prototype for function defined below */
static int processSegments(int iudp,                int *npnt, pnt_T *pnt_p[], int *nseg, seg_T *seg_p[], int *NumUdp, udp_T *udps);
static int processFile(ego context, char message[], int *npnt, pnt_T *pnt_p[], int *nseg, seg_T *seg_p[], int *NumUdp, udp_T *udps);
static int modifySegments(                          int *npnt, pnt_T *pnt_p[], int *nseg, seg_T *seg_p[], char message[]);
static int getToken(char *text, int nskip, char sep, int maxtok, char *token);

#ifdef GRAFIC
static int     plotWaffle(int npnt, pnt_T pnt[], int nseg, seg_T seg[]);
#endif

#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))


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

    int     ipnt, npnt, nnode, iseg, jseg, nseg, jedge, jface, senses[4];
    int     ibeg, iend, seginfo[2], iattr, nodes[2];
    double  xyz[20], data[6], xyz_out[20];
    pnt_T   *pnt=NULL;
    seg_T   *seg=NULL;
    char    *message=NULL;
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL, ecurve, echild[4], eloop, eshell;
    ego     *ewires=NULL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("depth(   0) = %f\n", DEPTH(     0));
    printf("segments(0) = %f",   SEGMENTS(0,0));
    for (iseg = 1; iseg < udps[0].arg[1].size; iseg++) {
        printf(" %f", SEGMENTS(0,iseg));
    }
    printf("\n");
    printf("filename(0) = %s\n", FILENAME(  0));
    printf("progress(0) = %d\n", PROGRESS(  0));
    printf("layout(  0) = %d\n", LAYOUT(    0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 1024);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[1].size == 1 && STRLEN(FILENAME(0)) == 0) {
        snprintf(message, 1024, "must specify segments or filename");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1 && STRLEN(FILENAME(0)) > 0) {
        snprintf(message, 1024, "must specify segments or filename");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[0].size > 1) {
        snprintf(message, 1024, "depth should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (DEPTH(0) <= 0) {
        snprintf(message, 1024, "depth = %f <= 0", DEPTH(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (STRLEN(FILENAME(0)) == 0 && udps[0].arg[1].size%4 != 0) {
        snprintf(message, 1024, "segments must be divisible by 4");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("depth(   %d) = %f\n", numUdp, DEPTH(   numUdp));
    printf("segments(%d) = %f",   numUdp, SEGMENTS(numUdp,0));
    for (iseg = 1; iseg < udps[0].arg[1].size; iseg++) {
        printf(" %f", SEGMENTS(numUdp,iseg));
    }
    printf("\n");
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    printf("progress(%d) = %d\n", numUdp, PROGRESS(numUdp));
    printf("layout(  %d) = %d\n", numUdp, LAYOUT(  numUdp));
#endif

    /* if filename is given, process the file */
    if (STRLEN(FILENAME(numUdp)) > 0) {
        status = processFile(context, message, &npnt, &pnt, &nseg, &seg, NumUdp, udps);
        CHECK_STATUS(processFile);

    /* otherwise, process the Segments */
    } else {
        status = processSegments(0, &npnt, &pnt, &nseg, &seg, NumUdp, udps);
        CHECK_STATUS(processSegments);
    }

    SPLINT_CHECK_FOR_NULL(pnt);
    SPLINT_CHECK_FOR_NULL(seg);

    /* if layout was selected, generate a MODEL with a WireBody for each Segment */
    if (LAYOUT(0) != 0) {
        MALLOC(ewires, ego, nseg);
        MALLOC(enodes, ego, 2   );
        MALLOC(eedges, ego, 1   );

        for (iseg = 0; iseg < nseg; iseg++) {
            ibeg = seg[iseg].ibeg;
            iend = seg[iseg].iend;

            data[0] = pnt[ibeg].x;
            data[1] = pnt[ibeg].y;
            data[2] = 0;
            data[3] = pnt[iend].x;
            data[4] = pnt[iend].y;
            data[5] = 0;

            status = EG_makeTopology(context, NULL, NODE, 0,
                                     &(data[0]), 0, NULL, NULL, &(enodes[0]));
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, NULL, NODE, 0,
                                     &(data[3]), 0, NULL, NULL, &(enodes[1]));
            CHECK_STATUS(EG_makeTopology);

            data[3] -= data[0];
            data[4] -= data[1];

            status = EG_makeGeometry(context, CURVE, LINE, NULL,
                                     NULL, data, &ecurve);
            CHECK_STATUS(EG_makeGeometry);

            data[0] = 0;
            data[1] = sqrt(data[3]*data[3] + data[4]*data[4]);

            status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                     data, 2, enodes, NULL, eedges);
            CHECK_STATUS(EG_makeTopology);

            seginfo[0] = SFORWARD;
            status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                     NULL, 1, eedges, seginfo, &eloop);
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                     NULL, 1, &eloop, NULL, &(ewires[iseg]));
            CHECK_STATUS(EG_makeTopology);

            for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
                status = EG_attributeAdd(ewires[iseg], seg[iseg].aname[iattr], ATTRSTRING,
                                         1, NULL, NULL, seg[iseg].avalu[iattr]);
                CHECK_STATUS(EG_attributeAdd);
            }
        }

        /* make a Model of the WIreBodys */
        status = EG_makeTopology(context, NULL, MODEL, 0,
                                 NULL, nseg, ewires, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* set the output value(s) */

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

        goto cleanup;
    }

    /* modify Segments so that they are connected end to end */
    status = modifySegments(&npnt, &pnt, &nseg, &seg, message);
    CHECK_STATUS(modifySegments);

    /* show Points and Segments after intersections */
    if (PROGRESS(numUdp) != 0) {
        printf("after intersections\n");
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            printf("        Pnt %3d: %-20s %1d %10.5f %10.5f\n", ipnt,
                   pnt[ipnt].name, pnt[ipnt].type, pnt[ipnt].x, pnt[ipnt].y);
        }
        for (iseg = 0; iseg < nseg; iseg++) {
            printf("        Seg %3d: %-20s %1d %5d %5d\n", iseg,
                   seg[iseg].name, seg[iseg].type, seg[iseg].ibeg, seg[iseg].iend);
        }
    }

    /* allocate space for the egos */
    MALLOC(enodes, ego, (2*npnt       ));
    MALLOC(eedges, ego, (  npnt+2*nseg));
    MALLOC(efaces, ego, (         nseg));

    /* set up enodes at Z=0 and at Z=Depth */
    nnode = 0;

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        xyz[0] = pnt[ipnt].x;
        xyz[1] = pnt[ipnt].y;
        xyz[2] = 0;

        status = EG_makeTopology(context, NULL, NODE, 0,
                                 xyz, 0, NULL, NULL, &(enodes[nnode]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeAdd(enodes[nnode], "__nodes__", ATTRINT, 1,
                                 &nnode, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        nnode++;
    }

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        xyz[0] = pnt[ipnt].x;
        xyz[1] = pnt[ipnt].y;
        xyz[2] = DEPTH(0);

        status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0,
                                 NULL, NULL, &(enodes[nnode]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeAdd(enodes[nnode], "__nodes__", ATTRINT, 1,
                                 &nnode, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        nnode++;
    }

    /* set up the edges on the z=0 plane */
    jedge = 0;

    for (iseg = 0; iseg < nseg; iseg++) {
        xyz[0] = pnt[seg[iseg].ibeg].x;
        xyz[1] = pnt[seg[iseg].ibeg].y;
        xyz[2] = 0;
        xyz[3] = pnt[seg[iseg].iend].x - pnt[seg[iseg].ibeg].x;
        xyz[4] = pnt[seg[iseg].iend].y - pnt[seg[iseg].ibeg].y;
        xyz[5] = 0;

        status = EG_makeGeometry(context, CURVE, LINE, NULL,
                                 NULL, xyz, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        echild[0] = enodes[seg[iseg].ibeg];
        echild[1] = enodes[seg[iseg].iend];

        status = EG_invEvaluate(ecurve, xyz, &(data[0]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        xyz[0] = pnt[seg[iseg].iend].x;
        xyz[1] = pnt[seg[iseg].iend].y;
        xyz[2] = 0;

        status = EG_invEvaluate(ecurve, xyz, &(data[1]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 data, 2, echild, NULL, &(eedges[jedge]));
        CHECK_STATUS(EG_makeTopology);

        nodes[0] = seg[iseg].ibeg;
        nodes[1] = seg[iseg].iend;

        status = EG_attributeAdd(eedges[jedge], "__nodes__", ATTRINT, 2,
                                 nodes, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        jedge++;
    }

    /* set up the edges on the z=Depth plane */
    for (iseg = 0; iseg < nseg; iseg++) {
        xyz[0] = pnt[seg[iseg].ibeg].x;
        xyz[1] = pnt[seg[iseg].ibeg].y;
        xyz[2] = DEPTH(0);
        xyz[3] = pnt[seg[iseg].iend].x - pnt[seg[iseg].ibeg].x;
        xyz[4] = pnt[seg[iseg].iend].y - pnt[seg[iseg].ibeg].y;
        xyz[5] = 0;

        status = EG_makeGeometry(context, CURVE, LINE, NULL,
                                 NULL, xyz, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        echild[0] = enodes[seg[iseg].ibeg+nnode/2];
        echild[1] = enodes[seg[iseg].iend+nnode/2];

        status = EG_invEvaluate(ecurve, xyz, &(data[0]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        xyz[0] = pnt[seg[iseg].iend].x;
        xyz[1] = pnt[seg[iseg].iend].y;
        xyz[2] = DEPTH(0);

        status = EG_invEvaluate(ecurve, xyz, &(data[1]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 data, 2, echild, NULL, &(eedges[jedge]));
        CHECK_STATUS(EG_makeTopology);

        nodes[0] = seg[iseg].ibeg + npnt;
        nodes[1] = seg[iseg].iend + npnt;

        status = EG_attributeAdd(eedges[jedge], "__nodes__", ATTRINT, 2,
                                 nodes, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        jedge++;
    }

    /* set up the edges betwe-en z=0 and Depth */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        xyz[0] = pnt[ipnt].x;
        xyz[1] = pnt[ipnt].y;
        xyz[2] = 0;
        xyz[3] = 0;
        xyz[4] = 0;
        xyz[5] = 1;

        status = EG_makeGeometry(context, CURVE, LINE, NULL,
                                 NULL, xyz, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        echild[0] = enodes[ipnt     ];
        echild[1] = enodes[ipnt+npnt];

        status = EG_invEvaluate(ecurve, xyz, &(data[0]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        xyz[0] = pnt[ipnt].x;
        xyz[1] = pnt[ipnt].y;
        xyz[2] = DEPTH(0);

        status = EG_invEvaluate(ecurve, xyz, &(data[1]), xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 data, 2, echild, NULL, &(eedges[jedge]));
        CHECK_STATUS(EG_makeTopology);

        nodes[0] = ipnt;
        nodes[1] = ipnt + npnt;

        status = EG_attributeAdd(eedges[jedge], "__nodes__", ATTRINT, 2,
                                 nodes, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        jedge++;
    }

    /* set up the Faces */
    jface = 0;

    for (iseg = 0; iseg < nseg; iseg++) {
        echild[0] = eedges[iseg];
        echild[1] = eedges[2*nseg+seg[iseg].iend];
        echild[2] = eedges[  nseg+    iseg      ];
        echild[3] = eedges[2*nseg+seg[iseg].ibeg];

        senses[0] = SFORWARD;
        senses[1] = SFORWARD;
        senses[2] = SREVERSE;
        senses[3] = SREVERSE;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                 NULL, 4, echild, senses, &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[jface]));
        CHECK_STATUS(EG_makeFace);

        jseg = iseg + 1;
        status = EG_attributeAdd(efaces[jface], "segment", ATTRINT,
                                 1, &jseg, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        seginfo[0] = seg[iseg].num;
        seginfo[1] = seg[iseg].idx;
        status = EG_attributeAdd(efaces[jface], "waffleseg", ATTRINT,
                                 2, seginfo, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
            status = EG_attributeAdd(efaces[jface], seg[iseg].aname[iattr], ATTRSTRING,
                                     1, NULL, NULL, seg[iseg].avalu[iattr]);
            CHECK_STATUS(EG_attributeAdd);
        }

        status = EG_attributeAdd(efaces[jface], "__jface__", ATTRINT, 1,
                                 &jface, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        jface++;
    }

    /* make the sheet body */
    status = EG_makeTopology(context, NULL, SHELL, OPEN,
                             NULL, jface, efaces, NULL, &eshell);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                             NULL, 1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(efaces);
    FREE(eedges);
    FREE(enodes);
    FREE(ewires);

    FREE(pnt);

    if (seg != NULL) {
        for (iseg = 0; iseg < nseg; iseg++) {
            if (seg[iseg].aname != NULL) {
                for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
                    FREE(seg[iseg].aname[iattr]);
                }
                FREE(seg[iseg].aname);
            }
            if (seg[iseg].avalu != NULL) {
                for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
                    FREE(seg[iseg].avalu[iattr]);
                }
                FREE(seg[iseg].avalu);
            }
        }
        FREE(seg);
    }

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
 *   processSegments - segments was specified                           *
 *                                                                      *
 ************************************************************************
 */

static int
processSegments(int    iudp,            /* (in)  udp index */
                int    *npnt,           /* (out) number of Points */
                pnt_T  *pnt_p[],        /* (out) array  of Points */
                int    *nseg,           /* (both)number of Segments */
                seg_T  *seg_p[],        /* (out) array  of Segments */
    /*@unused@*/int    *NumUdp,
                udp_T  *udps)
{
    int     status = EGADS_SUCCESS;

    int     mpnt, jpnt, mseg, iseg;
    pnt_T   *pnt=NULL;
    seg_T   *seg=NULL;

    ROUTINE(processSegments);

    /* --------------------------------------------------------------- */

    /* get number of Segments from size of SEGMENTS() */
    *npnt = 0;
    *nseg = udps[iudp].arg[1].size / 4;

    /* default (empty) Point and Segment tables */
    mpnt = *nseg * 2;   /* perhaps too big */
    mseg = *nseg;

    MALLOC(pnt, pnt_T, mpnt);
    MALLOC(seg, seg_T, mseg);

    /* create Segments and unique Points */
    for (iseg = 0; iseg < *nseg; iseg++) {

        /* check if beginning of Segment matches any Points so far */
        seg[iseg].ibeg = -1;
        for (jpnt = 0; jpnt < *npnt; jpnt++) {
            if (fabs(SEGMENTS(iudp,4*iseg  )-pnt[jpnt].x) < EPS06 &&
                fabs(SEGMENTS(iudp,4*iseg+1)-pnt[jpnt].y) < EPS06   ) {
                seg[iseg].ibeg = jpnt;
                break;
            }
        }

        if (seg[iseg].ibeg < 0) {
            seg[iseg].ibeg = *npnt;

            pnt[*npnt].type  = 1;
            pnt[*npnt].x     = SEGMENTS(    iudp,4*iseg  );
            pnt[*npnt].x_dot = SEGMENTS_DOT(iudp,4*iseg  );
            pnt[*npnt].y     = SEGMENTS(    iudp,4*iseg+1);
            pnt[*npnt].y_dot = SEGMENTS_DOT(iudp,4*iseg+1);
            pnt[*npnt].name[0] = '\0';
            (*npnt)++;
        }

        /* check if end of Segment matches any Points so far */
        seg[iseg].iend = -1;
        for (jpnt = 0; jpnt < *npnt; jpnt++) {
            if (fabs(SEGMENTS(iudp,4*iseg+2)-pnt[jpnt].x) < EPS06 &&
                fabs(SEGMENTS(iudp,4*iseg+3)-pnt[jpnt].y) < EPS06   ) {
                seg[iseg].iend = jpnt;
                break;
            }
        }

        if (seg[iseg].iend < 0) {
            seg[iseg].iend = *npnt;

            pnt[*npnt].type  = 1;
            pnt[*npnt].x     = SEGMENTS(    iudp,4*iseg+2);
            pnt[*npnt].x_dot = SEGMENTS_DOT(iudp,4*iseg+2);
            pnt[*npnt].y     = SEGMENTS(    iudp,4*iseg+3);
            pnt[*npnt].y_dot = SEGMENTS_DOT(iudp,4*iseg+3);
            pnt[*npnt].name[0] = '\0';
            (*npnt)++;
        }

        seg[iseg].type  = 1;
        seg[iseg].num   = iseg + 1;
        seg[iseg].idx   = 1;
        seg[iseg].nattr = 0;
        seg[iseg].aname = NULL;
        seg[iseg].avalu = NULL;
    }

cleanup:
    /* return pointers to freeable memory */
    *pnt_p = pnt;
    *seg_p = seg;

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   processFile - filename was specified                               *
 *                                                                      *
 ************************************************************************
 */

/*
    notes:  keywords can either be lowercase or UPPERCASE (not mixedCase)
            keywords are shown here in UPPERCASE to distinguish them from variables
            segments are only created for LINEs
            segments are broken at intersection of other LINEs and at POINTs

    POINT  pointname AT  xloc           yloc            creates point at <xloc,yloc>
                     AT  x@pointname+dx y@pointname+dy  creates point <dx,dy> from named point
                     AT  xloc           y@pointname     creates point at same y as named point and at given xloc
                     AT  x@pointname    yloc            creates point at same x as named point and at given yloc
                     ON  linename FRAC  fractDist       creates point on line at given fractional distance
                     ON  linename XLOC  xloc            creates point on line at given xloc
                     ON  linename YLOC  yloc            creates point on line at given yloc
                     ON  linename PERP  pointname       creates point on line that is closest to point
                     ON  linename XSECT linename        creates point at intersection of two lines
                     OFF linename dist  pointname       creates point dist to the left of linename at pointname
    CPOINT -------------- same as POINT --------------- creates construction point ...
    LINE   linename  pointname pointname [attrName1=attrValue1 [...]]
                                                        creates line between points with given attributes
    CLINE  -------------- same as LINE ---------------- creates construction line ...
    SET    varname expression                           sets varname to expression
    PATBEG varname ncopy                                loops ncopy times with varname=1,...,ncopy
    PATBREAK val1 op val2                               op is LT LE EQ GE GT or NE
                                                        breaks if (val1 op val2) > 0
    PATEND
    IFTHEN val1 op val2                                 op is LT LE EQ GE GT or NE
    ENDIF
*/

static int
processFile(ego    context,             /* (in)  EGADS context */
            char   message[],           /* (out) error message */
            int    *npnt,               /* (out) number of Points */
            pnt_T  *pnt_p[],            /* (out) array  of Points */
            int    *nseg,               /* (both)number of Segments */
            seg_T  *seg_p[],            /* (out) array  of Segments */
            int    *NumUdp,
            udp_T  *udps)
{
    int    status = EGADS_SUCCESS;      /* (out) default return */

    int    nbrch, npmtr, npmtr_save=999999, ipmtr, nbody, type, nrow, ncol, itoken, ichar;
    int    mpnt=0, ipnt, jpnt, mseg=0, iseg, jseg, ibeg, iend, jbeg, jend, itype;
    int    mline, nline=0, mset=0, nset=0, i, *set_pmtr=NULL;
    int    *begline=NULL, *endline=NULL;
    int    pat_pmtr[ ]={ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int    pat_value[]={ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int    pat_end[  ]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    long   pat_seek[ ]={ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int    istream, ieof, iskip, ifthen, npat=0, patbreak, outLevel;
    double xloc,     yloc,     xvalue,     yvalue,     value,     val1,     val2;
    double xloc_dot, yloc_dot, xvalue_dot, yvalue_dot, value_dot, val1_dot, val2_dot;
    double frac,     N,     D,     s,     dist,     dx,     dy,     alen;
    double frac_dot, N_dot, D_dot, s_dot, dist_dot, dx_dot, dy_dot, alen_dot;
    pnt_T  *pnt=NULL;
    seg_T  *seg=NULL;
    char   templine[256], token[256], pname1[256], pname2[256], lname1[256], lname2[256];
    char   *filename, name[MAX_NAME_LEN], str[256];
    modl_T *MODL=NULL;
    void   *modl;
    void   *realloc_temp = NULL;            /* used by RALLOC macro */
    FILE   *fp=NULL;

    ROUTINE(processFile);

    /* --------------------------------------------------------------- */

    /* initially there are no Points or Segments */
    *npnt = 0;
    *nseg = 0;

    /* get pointer to model */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;
    if (MODL->perturb != NULL) {
        MODL = MODL->perturb;
    }

    /* get the outLevel from OpenCSM */
    outLevel = ocsmSetOutLevel(-1);

    /* make sure that there are no Parameters whose names start with
       x@ or y@ */
    status = ocsmInfo(MODL, &nbrch, &npmtr_save, &nbody);
    CHECK_STATUS(ocsmInfo);

    for (ipmtr = 1; ipmtr <= npmtr_save; ipmtr++) {
        status = ocsmGetPmtr(MODL, ipmtr, &type, &nrow, &ncol, name);
        CHECK_STATUS(ocsmGetPmtr);

        if (strncmp(name, "x@", 2) == 0 ||
            strncmp(name, "y@", 2) == 0   ) {
            snprintf(message, 1024, "cannot use waffle if you already have a Parameter named \"%s\"", name);
            status = EGADS_NODATA;
            goto cleanup;
        }
    }

    /* determine if filename actually contains the name of a file
       or a copy of a virtual file */
    filename = FILENAME(numUdp);

    if (strncmp(filename, "<<\n", 3) == 0) {
        istream = 1;
    } else {
        istream = 0;
    }

    /* if a file, open it now and remember if we are at the end of file */
    if (istream == 0) {
        fp = fopen(filename, "rb");
        if (fp == NULL) {
            snprintf(message, 1024, "processFile: could not open file \"%s\"", filename);
            status = EGADS_NOTFOUND;
            goto cleanup;
        }

        ieof = feof(fp);

    /* otherwise, find extents of each line of the input */
    } else {
        mline = 25;
        nline = 0;
        iend  = 0;

        MALLOC(begline, int, mline);
        MALLOC(endline, int, mline);

        for (i = 3; i < STRLEN(filename); i++) {
            if (iend == 0) {
                if (filename[i] == ' '  || filename[i] == '\t' ||
                    filename[i] == '\r' || filename[i] == '\n'   ) continue;

                iend = i;

                if (nline >= mline) {
                    mline += 25;
                    RALLOC(begline, int, mline);
                    RALLOC(endline, int, mline);
                }

                begline[nline] = i;
                endline[nline] = i + 1;
                nline++;
            } else if (filename[i] == '\r') {
                filename[i] = ' ';
            } else if (filename[i] == '\n') {
                endline[nline-1] = iend + 1;
                iend = 0;
            } else if (filename[i] != ' ' && filename[i] != '\t') {
                iend = i;
            }
        }

        /* remember if we are at the end of the stream */
        if (strlen(filename) <= 2) {
            ieof = 1;
        } else {
            ieof = 0;
        }
    }

    /* print copy of input file */
//$$$    for (i = 0; i < nline; i++) {
//$$$        memcpy(templine, filename+begline[i], endline[i]-begline[i]);
//$$$        templine[endline[i]-begline[i]] = '\0';
//$$$
//$$$        printf("%5d :%s:\n", i, templine);
//$$$    }

    /* default (empty) Point and Segment tables */
    mpnt = 10;
    mseg = 10;

    MALLOC(pnt, pnt_T, mpnt);
    MALLOC(seg, seg_T, mseg);

    /* by default, we are not skipping */
    iskip = 0;

    /* we start without any pending ifthens */
    ifthen = 0;

    /* read the file */
    while (ieof == 0) {

        /* read the next line */
        if (fp != NULL) {
            (void) fgets(templine, 255, fp);
            if (feof(fp) > 0) break;

            /* overwite the \n and \r at the end */
            if (STRLEN(templine) > 0 && templine[STRLEN(templine)-1] == '\n') {
                templine[STRLEN(templine)-1] = '\0';
            }
            if (STRLEN(templine) > 0 && templine[STRLEN(templine)-1] == '\r') {
                templine[STRLEN(templine)-1] = '\0';
            }
        } else {
            if (istream-1 >= nline) {
                break;
            }

            SPLINT_CHECK_FOR_NULL(begline);
            SPLINT_CHECK_FOR_NULL(endline);

            memcpy(templine, filename+begline[istream-1], endline[istream-1]-begline[istream-1]);
            templine[endline[istream-1]-begline[istream-1]] = '\0';

            istream++;
        }

        if (outLevel >= 1) {
            printf("    processing: %s", templine);
            if (fp == NULL) printf("\n");
        }

        /* get and process the first token */
        status = getToken(templine, 0, ' ', 255, token);
        if (status < SUCCESS) {
            snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
            goto cleanup;
        }

        /* skip blank line */
        if (strcmp(token, "") == 0) {
            continue;

        /* skip comment */
        } else if (token[0] == '#') {
            continue;

        /* skip comment */
        } else if (strcmp(token, "#") == 0) {
            continue;

        /* ENDIF */
        } else if (strcmp(token, "endif") == 0 ||
                   strcmp(token, "ENDIF") == 0   ) {
            if (ifthen > 0) {
                ifthen--;
            }

        /* skip line if there is an active IFTHEN */
        } else if (ifthen > 0) {
            if (outLevel >= 1) {
                printf("    ...skipping\n");
            }
            continue;

        /* POINT and CPOINT */
        } else if (strcmp(token,  "point") == 0 ||
                   strcmp(token,  "POINT") == 0 ||
                   strcmp(token, "cpoint") == 0 ||
                   strcmp(token, "CPOINT") == 0   ) {

            if (iskip > 0) continue;

            if (strcmp(token, "point") == 0 ||
                strcmp(token, "POINT") == 0   ) {
                itype = 1;
            } else {
                itype = 0;
            }

            /* pname1 */
            status = getToken(templine, 1, ' ', 255, pname1);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find pname1\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = getToken(templine, 2, ' ', 255, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            /* POINT pname1 AT xloc yloc */
            if        (strcmp(token, "at") == 0 ||
                       strcmp(token, "AT") == 0   ) {

                /* get xloc */
                status = getToken(templine, 3, ' ', 255, token);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find third token\nwhile processing: %s", templine);
                    goto cleanup;
                }

                status = ocsmEvalExpr(MODL, token, &xloc, &xloc_dot, str);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                    goto cleanup;
                }

                if (STRLEN(str) > 0) {
                    snprintf(message, 1024, "xvalue must be a number\nwhile processing: %s", templine);
                    status = EGADS_NODATA;
                    goto cleanup;
                }

                /* get yloc */
                status = getToken(templine, 4, ' ', 255, token);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find fourth token\nwhile processing: %s", templine);
                    goto cleanup;
                }

                status = ocsmEvalExpr(MODL, token, &yloc, &yloc_dot, str);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                    goto cleanup;
                } else if (STRLEN(str) > 0) {
                    snprintf(message, 1024, "yvalue must be a number\nwhile processing: %s", templine);
                    status = EGADS_NODATA;
                    goto cleanup;
                }

                /* compute (xvalue,yvalue) */
                xvalue     = xloc;
                xvalue_dot = xloc_dot;
                yvalue     = yloc;
                yvalue_dot = yloc_dot;

            /* POINT pname1 ON lname1 ... */
            } else if (strcmp(token, "on") == 0 ||
                       strcmp(token, "ON") == 0   ) {

                /* find lname1 */
                status = getToken(templine, 3, ' ', 255, lname1);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find lname1\nwhile processing: %s", templine);
                    goto cleanup;
                }

                iseg = -1;
                for (jseg = 0; jseg < *nseg; jseg++) {
                    if (strcmp(lname1, seg[jseg].name) == 0) {
                        iseg = jseg;
                    }
                }
                if (iseg < 0) {
                    snprintf(message, 1024, "line \"%s\" could not be found\nwhile processing: %s", lname1, templine);
                    status = EGADS_NOTFOUND;
                    goto cleanup;
                }

                ibeg = seg[iseg].ibeg;
                iend = seg[iseg].iend;

                /* determine the sub-operator */
                status = getToken(templine, 4, ' ', 255, token);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot determine sub-operator\nwhile processing: %s", templine);
                    goto cleanup;
                }

                /* POINT pname1 ON lname1 FRAC fracDist */
                if        (strcmp(token, "frac") == 0 ||
                           strcmp(token, "FRAC") == 0   ) {

                    /* get fracDist */
                    status = getToken(templine, 5, ' ', 255, token);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot find fifth token\nwhile processing: %s", templine);
                        goto cleanup;
                    }

                    status = ocsmEvalExpr(MODL, token, &frac, &frac_dot, str);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                        goto cleanup;
                    }

                    /* compute (xvalue,yvalue) */
                    xvalue     = (1-frac) * pnt[ibeg].x     + frac * pnt[iend].x;
                    xvalue_dot = (1-frac) * pnt[ibeg].x_dot + frac * pnt[iend].x_dot + frac_dot * (pnt[iend].x - pnt[ibeg].x);
                    yvalue     = (1-frac) * pnt[ibeg].y     + frac * pnt[iend].y;
                    yvalue_dot = (1-frac) * pnt[ibeg].y_dot + frac * pnt[iend].y_dot + frac_dot * (pnt[iend].y - pnt[ibeg].y);

                /* POINT pname1 ON lname1 XLOC xloc */
                } else if (strcmp(token, "xloc") == 0 ||
                           strcmp(token, "XLOC") == 0   ) {

                    if(fabs(pnt[ibeg].x-pnt[iend].x) < EPS06) {
                        snprintf(message, 1024, "cannot specify XLOC on a constant X line\nwhile processing: %s", templine);
                        status = EGADS_RANGERR;
                        goto cleanup;
                    }

                    /* get xloc */
                    status = getToken(templine, 5, ' ', 255, token);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot find fifth token\nwhile processing: %s", templine);
                        goto cleanup;
                    }

                    status = ocsmEvalExpr(MODL, token, &xloc, &xloc_dot, str);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                        goto cleanup;
                    }


                    /* compute (xvalue,yvalue) */
                    frac     = (xloc - pnt[ibeg].x) / (pnt[iend].x - pnt[ibeg].x);
                    frac_dot = ((pnt[iend].x - pnt[ibeg].x) * (xloc_dot        - pnt[ibeg].x_dot)
                               -(xloc        - pnt[ibeg].x) * (pnt[iend].x_dot - pnt[ibeg].x_dot))
                             /  (pnt[iend].x - pnt[ibeg].x) / (pnt[iend].x     - pnt[ibeg].x);

                    xvalue     = xloc;
                    xvalue_dot = xloc_dot;
                    yvalue     = (1-frac) * pnt[ibeg].y     + frac * pnt[iend].y;
                    yvalue_dot = (1-frac) * pnt[ibeg].y_dot + frac * pnt[iend].y_dot + frac_dot * (pnt[iend].y - pnt[ibeg].y);

                /* POINT pname1 ON lname1 YLOC yloc */
                } else if (strcmp(token, "yloc") == 0 ||
                           strcmp(token, "YLOC") == 0   ) {

                    if(fabs(pnt[ibeg].y-pnt[iend].y) < EPS06) {
                        snprintf(message, 1024, "cannot specify YLOC on a constant Y line\nwhile processing: %s", templine);
                        status = EGADS_RANGERR;
                        goto cleanup;
                    }

                    /* get yloc */
                    status = getToken(templine, 5, ' ', 255, token);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot find fifth token\nwhile processing: %s", templine);
                        goto cleanup;
                    }

                    status = ocsmEvalExpr(MODL, token, &yloc, &yloc_dot, str);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                        goto cleanup;
                    }

                    /* compute (xvalue,yvalue) */
                    frac     = (yloc - pnt[ibeg].y) / (pnt[iend].y - pnt[ibeg].y);
                    frac_dot = ((pnt[iend].y - pnt[ibeg].y) * (yloc_dot        - pnt[ibeg].y_dot)
                               -(yloc        - pnt[ibeg].y) * (pnt[iend].y_dot - pnt[ibeg].y_dot))
                             /  (pnt[iend].y - pnt[ibeg].y) / (pnt[iend].y     - pnt[ibeg].y);

                    xvalue     = (1-frac) * pnt[ibeg].x     + frac * pnt[iend].x;
                    xvalue_dot = (1-frac) * pnt[ibeg].x_dot + frac * pnt[iend].x_dot + frac_dot * (pnt[iend].x - pnt[ibeg].x);
                    yvalue     = yloc;
                    yvalue_dot = yloc_dot;

                /* POINT pname1 ON lname1 PERP pname2 */
                } else if (strcmp(token, "perp") == 0 ||
                           strcmp(token, "PERP") == 0   ) {

                    /* get pname2 */
                    status = getToken(templine, 5, ' ', 255, pname2);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot find fofth token\nwhile processing: %s", templine);
                        goto cleanup;
                    }

                    ipnt = -1;
                    for (jpnt = 0; jpnt < *npnt; jpnt++) {
                        if (strcmp(pname2, pnt[jpnt].name) == 0) {
                            ipnt = jpnt;
                        }
                    }
                    if (ipnt < 0) {
                        snprintf(message, 1024, "point \"%s\" could not be found\nwhile processing: %s", pname2, templine);
                        status = EGADS_NOTFOUND;
                        goto cleanup;
                    }

                    /* compute (xvalue,yvalue) */
                    N     = (pnt[ipnt].x     - pnt[ibeg].x    ) * (pnt[iend].x     - pnt[ibeg].x    )
                          + (pnt[ipnt].y     - pnt[ibeg].y    ) * (pnt[iend].y     - pnt[ibeg].y    );
                    N_dot = (pnt[ipnt].x_dot - pnt[ibeg].x_dot) * (pnt[iend].x     - pnt[ibeg].x    )
                          + (pnt[ipnt].x     - pnt[ibeg].x    ) * (pnt[iend].x_dot - pnt[ibeg].x_dot)
                          + (pnt[ipnt].y_dot - pnt[ibeg].y_dot) * (pnt[iend].y     - pnt[ibeg].y    )
                          + (pnt[ipnt].y     - pnt[ibeg].y    ) * (pnt[iend].y_dot - pnt[ibeg].y_dot);
                    D     = (pnt[iend].x     - pnt[ibeg].x    ) * (pnt[iend].x     - pnt[ibeg].x    )
                          + (pnt[iend].y     - pnt[ibeg].y    ) * (pnt[iend].y     - pnt[ibeg].y    );
                    D_dot = (pnt[iend].x_dot - pnt[ibeg].x_dot) * (pnt[iend].x     - pnt[ibeg].x    )
                          + (pnt[iend].x     - pnt[ibeg].x    ) * (pnt[iend].x_dot - pnt[ibeg].x_dot)
                          + (pnt[iend].y_dot - pnt[ibeg].y_dot) * (pnt[iend].y     - pnt[ibeg].y    )
                          + (pnt[iend].y     - pnt[ibeg].y    ) * (pnt[iend].y_dot - pnt[ibeg].y_dot);

                    frac     = N / D;
                    frac_dot = (D * N_dot - N * D_dot) / D / D;

                    xvalue     = (1-frac) * pnt[ibeg].x     + frac * pnt[iend].x;
                    xvalue_dot = (1-frac) * pnt[ibeg].x_dot + frac * pnt[iend].x_dot + frac_dot * (pnt[iend].x - pnt[ibeg].x);
                    yvalue     = (1-frac) * pnt[ibeg].y     + frac * pnt[iend].y;
                    yvalue_dot = (1-frac) * pnt[ibeg].y_dot + frac * pnt[iend].y_dot + frac_dot * (pnt[iend].y - pnt[ibeg].y);

                /* POINT pname1 ON lname1 XSECT lname2 */
                } else if (strcmp(token, "xsect") == 0 ||
                           strcmp(token, "XSECT") == 0   ) {

                    /* find lname2 */
                    status = getToken(templine, 5, ' ', 255, lname2);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot find fifth token\nwhile processing: %s", templine);
                        goto cleanup;
                    }

                    iseg = -1;
                    for (jseg = 0; jseg < *nseg; jseg++) {
                        if (strcmp(lname2, seg[jseg].name) == 0) {
                            iseg = jseg;
                        }
                    }
                    if (iseg < 0) {
                        snprintf(message, 1024, "line \"%s\" could not be found\nwhile processing: %s", lname2, templine);
                        status = EGADS_NOTFOUND;
                        goto cleanup;
                    }

                    jbeg = seg[iseg].ibeg;
                    jend = seg[iseg].iend;

                    /* compute (xvalue,yvalue) */
                    D     = (pnt[iend].x     - pnt[ibeg].x    ) * (pnt[jbeg].y     - pnt[jend].y    )
                          - (pnt[jbeg].x     - pnt[jend].x    ) * (pnt[iend].y     - pnt[ibeg].y    );
                    D_dot = (pnt[iend].x_dot - pnt[ibeg].x_dot) * (pnt[jbeg].y     - pnt[jend].y    )
                          + (pnt[iend].x     - pnt[ibeg].x    ) * (pnt[jbeg].y_dot - pnt[jend].y_dot)
                          - (pnt[jbeg].x_dot - pnt[jend].x_dot) * (pnt[iend].y     - pnt[ibeg].y    )
                          - (pnt[jbeg].x     - pnt[jend].x    ) * (pnt[iend].y_dot - pnt[ibeg].y_dot);

                    if (fabs(D) > EPS06) {
                        N     = (pnt[jbeg].x     - pnt[ibeg].x    ) * (pnt[jbeg].y     - pnt[jend].y    )
                              - (pnt[jbeg].x     - pnt[jend].x    ) * (pnt[jbeg].y     - pnt[ibeg].y    );
                        N_dot = (pnt[jbeg].x_dot - pnt[ibeg].x_dot) * (pnt[jbeg].y     - pnt[jend].y    )
                              + (pnt[jbeg].x     - pnt[ibeg].x    ) * (pnt[jbeg].y_dot - pnt[jend].y_dot)
                              - (pnt[jbeg].x_dot - pnt[jend].x_dot) * (pnt[jbeg].y     - pnt[ibeg].y    )
                              - (pnt[jbeg].x     - pnt[jend].x    ) * (pnt[jbeg].y_dot - pnt[ibeg].y_dot);

                        s     = N / D;
                        s_dot = (D * N_dot - N * D_dot) / D / D;

                        xvalue     = (1-s) * pnt[ibeg].x     + s * pnt[iend].x;
                        xvalue_dot = (1-s) * pnt[ibeg].x_dot + s * pnt[iend].x_dot + s_dot * (pnt[iend].x - pnt[ibeg].x);
                        yvalue     = (1-s) * pnt[ibeg].y     + s * pnt[iend].y;
                        yvalue_dot = (1-s) * pnt[ibeg].y_dot + s * pnt[iend].y_dot + s_dot * (pnt[iend].y - pnt[ibeg].y);
                    } else {
                        snprintf(message, 1024, "segments do not intersect\nwhile processing: %s", templine);
                        status = EGADS_NOTFOUND;
                        goto cleanup;
                    }

                } else {
                    snprintf(message, 1024, "fifth token should be FRAC, PERP, XLOC, YLOC, SAMEX, SAMEY, or XSECT\nwhile processing: %s", templine);
                    status = EGADS_RANGERR;
                    goto cleanup;
                }

            /* POINT pname1 OFF lname1 dist pname2 */
            } else if (strcmp(token, "off") == 0 ||
                       strcmp(token, "OFF") == 0   ) {

                /* find lname1 */
                status = getToken(templine, 3, ' ', 255, lname1);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find third token\nwhile processing: %s", templine);
                    goto cleanup;
                }

                iseg = -1;
                for (jseg = 0; jseg < *nseg; jseg++) {
                    if (strcmp(lname1, seg[jseg].name) == 0) {
                        iseg = jseg;
                    }
                }
                if (iseg < 0) {
                    snprintf(message, 1024, "line \"%s\" could not be found\nwhile processing: %s", lname1, templine);
                    status = EGADS_NOTFOUND;
                    goto cleanup;
                }

                ibeg = seg[iseg].ibeg;
                iend = seg[iseg].iend;

                /* get dist */
                status = getToken(templine, 4, ' ', 255, token);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find fourth token\nwhile processing: %s", templine);
                    goto cleanup;
                }

                status = ocsmEvalExpr(MODL, token, &dist, &dist_dot, str);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                    goto cleanup;
                }

                /* get pname2 */
                status = getToken(templine, 5, ' ', 255, pname2);
                if (status < SUCCESS) {
                    snprintf(message, 1024, "cannot find fifth token\nwhile processing: %s", templine);
                    goto cleanup;
                }

                ipnt = -1;
                for (jpnt = 0; jpnt < *npnt; jpnt++) {
                    if (strcmp(pname2, pnt[jpnt].name) == 0) {
                        ipnt = jpnt;
                    }
                }
                if (ipnt < 0) {
                    snprintf(message, 1024, "point \"%s\" could not be found\nwhile processing: %s", pname2, templine);
                    status = EGADS_NOTFOUND;
                    goto cleanup;
                }

                /* compute (xvalue,yvalue) */
                dx       = pnt[iend].x     - pnt[ibeg].x;
                dx_dot   = pnt[iend].x_dot - pnt[ibeg].x_dot;
                dy       = pnt[iend].y     - pnt[ibeg].y;
                dy_dot   = pnt[iend].y_dot - pnt[ibeg].y_dot;

                alen     = sqrt(dx * dx + dy * dy);
                alen_dot = (dx * dx_dot + dy * dy_dot) / alen;

                xvalue     = pnt[ipnt].x     - dist     * dy / alen;
                xvalue_dot = pnt[ipnt].x_dot - dist_dot * dy / alen - dy_dot * dist / alen + alen_dot * dist * dy / alen / alen;
                yvalue     = pnt[ipnt].y     + dist     * dx / alen;
                yvalue_dot = pnt[ipnt].y_dot + dist_dot * dx / alen + dx_dot * dist / alen - alen_dot * dist * dx / alen / alen;

            } else {
                snprintf(message, 1024, "third token should be AT, ON, or OFF\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

            /* make room for (possible) new Point */
            if (*npnt+1 >= mpnt) {
                mpnt += 10;

                RALLOC(pnt, pnt_T, mpnt);
            }

            /* if another Point has same name, remove the name */
            for (jpnt = 0; jpnt < *npnt; jpnt++) {
                if (strcmp(pnt[jpnt].name, pname1) == 0) {
                    strcpy(pnt[jpnt].name, "");
                    break;
                }
            }

            /* see if Point already exists at these coordintes  */
            ipnt = -1;
            for (jpnt = 0; jpnt < *npnt; jpnt++) {
                if (fabs(xvalue-pnt[jpnt].x) < EPS06 &&
                    fabs(yvalue-pnt[jpnt].y) < EPS06   ) {
                    ipnt = jpnt;
                    break;
                }
            }

            /* if Point already exists, just update its name */
            if (ipnt >= 0) {
                strcpy(pnt[ipnt].name, pname1);

                /* otherwise, create the new Point */
            } else {
                strcpy(pnt[*npnt].name, pname1);

                pnt[*npnt].type  = itype;
                pnt[*npnt].x     = xvalue;
                pnt[*npnt].x_dot = xvalue_dot;
                pnt[*npnt].y     = yvalue;
                pnt[*npnt].y_dot = yvalue_dot;
                (*npnt)++;

            }

            /* store coordinates in a local variable */
            (void) strcpy(str, "x@");
            (void) strcat(str, pname1);

            status = ocsmFindPmtr(MODL, str, OCSM_LOCALVAR, 1, 1, &ipmtr);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find \"%s\"\nwhile processing: %s", str, templine);
                goto cleanup;
            }

            /* note that we cannot call ocsmSetVeluD and ocsmSetVelD because
               they remove a perturbation and remove the velocities */
            MODL->pmtr[ipmtr].value[0] = xvalue;
            MODL->pmtr[ipmtr].dot[  0] = xvalue_dot;

            (void) strcpy(str, "y@");
            (void) strcat(str, pname1);

            status = ocsmFindPmtr(MODL, str, OCSM_LOCALVAR, 1, 1, &ipmtr);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find \"%s\"\nwhile processing: %s", str, templine);
                goto cleanup;
            }

            MODL->pmtr[ipmtr].value[0] = yvalue;
            MODL->pmtr[ipmtr].dot[  0] = yvalue_dot;

        } else if (strcmp(token,  "line") == 0 ||
                   strcmp(token,  "LINE") == 0 ||
                   strcmp(token, "cline") == 0 ||
                   strcmp(token, "CLINE") == 0 ||
                   strcmp(token,  "path") == 0 ||       // for backward compatability
                   strcmp(token,  "PATH") == 0   ) {    // for backward compatability

            if (iskip > 0) continue;

            status = getToken(templine, 1, ' ', 255, lname1);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = getToken(templine, 2, ' ', 255, pname1);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = getToken(templine, 3, ' ', 255, pname2);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find third token\nwhile processing: %s", templine);
                goto cleanup;
            }

            /* make room in the Segment table */
            if (*nseg+1 >= mseg) {
                mseg += 10;

                RALLOC(seg, seg_T, mseg);
            }

            /* set the Segment type (0  for cline, 1 for line) */
            if (strcmp(token, "line") == 0 ||
                strcmp(token, "LINE") == 0   ) {
                seg[*nseg].type = 1;
            } else {
                seg[*nseg].type = 0;
            }

            /* use last Point named pname1 (which takes care of redefinitions) */
            seg[*nseg].ibeg = -1;
            for (ipnt = *npnt-1; ipnt >= 0; ipnt--) {
                if (strcmp(pname1, pnt[ipnt].name) == 0) {
                    seg[*nseg].ibeg = ipnt;
                    break;
                }
            }
            if (seg[*nseg].ibeg < 0) {
                snprintf(message, 1024, "\"%s\" not found\nwhile processing: %s", pname1, templine);
                status = EGADS_NODATA;
                goto cleanup;
            }

            /* use last Point named pname2 (which takes care of redefinitions) */
            seg[*nseg].iend = -1;
            for (ipnt = *npnt-1; ipnt >= 0; ipnt--) {
                if (strcmp(pname2, pnt[ipnt].name) == 0) {
                    seg[*nseg].iend = ipnt;
                    break;
                }
            }
            if (seg[*nseg].iend < 0) {
                snprintf(message, 1024, "\"%s\" not found\nwhile processing: %s", pname2, templine);
                status = EGADS_NODATA;
                goto cleanup;
            }

            /* create Segment */
            seg[*nseg].num     = *nseg + 1;
            seg[*nseg].idx     = 1;
            seg[*nseg].name[0] = '\0';
            seg[*nseg].nattr   = 0;
            seg[*nseg].aname   = NULL;
            seg[*nseg].avalu   = NULL;

            strncpy(seg[*nseg].name, lname1, 80);

            /* process Attribute pairs */
            for (itoken = 4; itoken < 100; itoken++) {
                (void) getToken(templine, itoken, ' ', 255, token);

                if (STRLEN(token) == 0) break;

                if (strstr(token, "=") == NULL) {
                    snprintf(message, 1024, "attribute pair must contain = sign\nwhile processing: %s", templine);
                    status = EGADS_RANGERR;
                    goto cleanup;
                } else {
                    if (seg[*nseg].nattr == 0) {
                        MALLOC(seg[*nseg].aname, char*, 1);
                        MALLOC(seg[*nseg].avalu, char*, 1);
                    } else {
                        RALLOC(seg[*nseg].aname, char*, (seg[*nseg].nattr+1));
                        RALLOC(seg[*nseg].avalu, char*, (seg[*nseg].nattr+1));
                    }

                    seg[*nseg].aname[seg[*nseg].nattr] = NULL;
                    seg[*nseg].avalu[seg[*nseg].nattr] = NULL;

                    MALLOC(seg[*nseg].aname[seg[*nseg].nattr], char, 80);
                    MALLOC(seg[*nseg].avalu[seg[*nseg].nattr], char, 80);

                    for (ichar = 0; ichar < STRLEN(token); ichar++) {
                        if (token[ichar] == '=') {
                            token[ichar] = '\0';
                            break;
                        }
                    }

                    strncpy(seg[*nseg].aname[seg[*nseg].nattr], token, 80);

                    token[ichar] = '$';
                    status = ocsmEvalExpr(MODL, &(token[ichar]), &value, &value_dot, str);
                    if (status < SUCCESS) {
                        snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", &(token[ichar]), templine);
                        goto cleanup;
                    }

                    if (STRLEN(str) == 0) {
                        snprintf(message, 1024, "attribute value must be a string\nwhile processing: %s", templine);
                        status = EGADS_NODATA;
                        goto cleanup;
                    }

                    strncpy(seg[*nseg].avalu[seg[*nseg].nattr], str, 80);

                    (seg[*nseg].nattr)++;
                }
            }

            (*nseg)++;

        } else if (strcmp(token, "set") == 0 ||
                   strcmp(token, "SET") == 0   ) {

            if (nset >= mset) {
                mset += 10;
                RALLOC(set_pmtr, int, mset);
            }

            SPLINT_CHECK_FOR_NULL(set_pmtr);

            /* get the expression */
            status = getToken(templine, 2, ' ', 255, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &value, &value_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            /* set up Parameter */
            status = getToken(templine, 1, ' ', 255, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmFindPmtr(MODL, token, OCSM_LOCALVAR, 1, 1, &set_pmtr[nset]);
            CHECK_STATUS(ocsm_localvar);

            if (set_pmtr[nset] <= npmtr_save) {
                snprintf(message, 1024, "cannot use \"%s\" since it was previously defined\nwhile processing: %s", token, templine);
                status = EGADS_NONAME;
                goto cleanup;
            }

            MODL->pmtr[set_pmtr[nset]].value[0] = value;

        } else if (strcmp(token, "patbeg") == 0 ||
                   strcmp(token, "PATBEG") == 0   ) {

            if (npat < 9) {
                npat++;
            } else {
                snprintf(message, 1024, "PATBEGs nested too deeply\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

            if (iskip > 0) {
                iskip++;
                continue;
            }

            /* remember where we in the file are so that we can get back here */
            if (fp != NULL) {
                pat_seek[npat] = ftell(fp);
            } else {
                pat_seek[npat] = istream;
            }

            if (iskip > 0) {
                pat_end[npat] = -1;
                iskip++;
                continue;
            }

            /* get the number of replicates */
            status = getToken(templine, 2, ' ', 255, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &value, &value_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            pat_end[npat] = NINT(value);

            if (pat_end[npat] <= 0) {
                iskip++;
                continue;
            } else {
                pat_value[npat] = 1;
            }

            /* set up Parameter to hold pattern index */
            status = getToken(templine, 1, ' ', 255, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmFindPmtr(MODL, token, OCSM_LOCALVAR, 1, 1, &pat_pmtr[npat]);
            CHECK_STATUS(ocsm_localvar);

            if (pat_pmtr[npat] <= npmtr_save) {
                snprintf(message, 1024, "cannot use \"%s\" as pattern variable since it was previously defined\nwhile processing: %s", token, templine);
                status = EGADS_NONAME;
                goto cleanup;
            }

            MODL->pmtr[pat_pmtr[npat]].value[0] = pat_value[npat];

        } else if (strcmp(token, "patbreak") == 0 ||
                   strcmp(token, "PATBREAK") == 0   ) {

            if (pat_end[npat] < 0) {
                snprintf(message, 1024, "PATBREAK without PATBEG\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

            if (iskip > 0) {
                iskip--;
                continue;
            }

            /* if the expression is not positive, get the next line */
            /* get val1, op, and val2 */
            status = getToken(templine, 1, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &val1, &val1_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            status = getToken(templine, 3, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find third token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &val2, &val2_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            status = getToken(templine, 2, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            patbreak = 0;
            if (strcmp(token, "lt") == 0 || strcmp(token, "LT") == 0) {
                if (val1 < val2) {
                    patbreak = 1;
                }
            } else if (strcmp(token, "le") == 0 || strcmp(token, "LE") == 0) {
                if (val1 <= val2) {
                    patbreak = 1;
                }
            } else if (strcmp(token, "eq") == 0 || strcmp(token, "EQ") == 0) {
                if (val1 == val2) {
                    patbreak = 1;
                }
            } else if (strcmp(token, "ge") == 0 || strcmp(token, "GE") == 0) {
                if (val1 >= val2) {
                    patbreak = 1;
                }
            } else if (strcmp(token, "gt") == 0 || strcmp(token, "GT") == 0) {
                if (val1 > val2) {
                    patbreak = 1;
                }
            } else if (strcmp(token, "ne") == 0 || strcmp(token, "NE") == 0) {
                if (val1 != val2) {
                    patbreak = 1;
                }
            } else {
                snprintf(message, 1024, "op must be LT LE EQ GE GT or NE\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

            if (patbreak == 1) {
                iskip++;
            }

        } else if (strcmp(token, "patend") == 0 ||
                   strcmp(token, "PATEND") == 0   ) {

            if (pat_end[npat] < 0) {
                snprintf(message, 1024, "PATEND without PATBEG\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

            if (iskip > 0) {
                iskip--;
                continue;
            }

            /* for another replicate, increment the value and place file pointer
               just after the PATBEG statement */
            if (pat_value[npat] < pat_end[npat]) {
                (pat_value[npat])++;

                MODL->pmtr[pat_pmtr[npat]].value[0] = pat_value[npat];

                if (pat_seek[npat] != 0) {
                    if (fp != NULL) {
                        fseek(fp, pat_seek[npat], SEEK_SET);
                    } else {
                        istream = pat_seek[npat];
                    }
                }

                /* otherwise, reset the pattern variables */
            } else {

                pat_pmtr[npat] = -1;
                pat_end[ npat] = -1;

                npat--;
            }

        } else if (strcmp(token, "ifthen") == 0 ||
                   strcmp(token, "IFTHEN") == 0   ) {

            /* get val1, op, and val2 */
            status = getToken(templine, 1, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find first token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &val1, &val1_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            status = getToken(templine, 3, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find third token\nwhile processing: %s", templine);
                goto cleanup;
            }

            status = ocsmEvalExpr(MODL, token, &val2, &val2_dot, str);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot evaluate \"%s\"\nwhile processing: %s", token, templine);
                goto cleanup;
            }

            status = getToken(templine, 2, ' ', 256, token);
            if (status < SUCCESS) {
                snprintf(message, 1024, "cannot find second token\nwhile processing: %s", templine);
                goto cleanup;
            }

            if (ifthen > 0) {
                ifthen++;
            } else if (strcmp(token, "lt") == 0 || strcmp(token, "LT") == 0) {
                if (val1 >= val2) {
                    ifthen++;
                }
            } else if (strcmp(token, "le") == 0 || strcmp(token, "LE") == 0) {
                if (val1 > val2) {
                    ifthen++;
                }
            } else if (strcmp(token, "eq") == 0 || strcmp(token, "EQ") == 0) {
                if (val1 != val2) {
                    ifthen++;
                }
            } else if (strcmp(token, "ge") == 0 || strcmp(token, "GE") == 0) {
                if (val1 < val2) {
                    ifthen++;
                }
            } else if (strcmp(token, "gt") == 0 || strcmp(token, "GT") == 0) {
                if (val1 <= val2) {
                    ifthen++;
                }
            } else if (strcmp(token, "ne") == 0 || strcmp(token, "NE") == 0) {
                if (val1 == val2) {
                    ifthen++;
                }
            } else {
                snprintf(message, 1024, "op must be LT LE EQ GE GT or NE\nwhile processing: %s", templine);
                status = EGADS_RANGERR;
                goto cleanup;
            }

        } else {
            snprintf(message, 1024, "input should start with CPOINT, POINT, LINE, CLINE, PATBEG, or PATEND\nwhile processing: %s", templine);
            status = EGADS_RANGERR;
            goto cleanup;
        }

        /* show Points and Segments after processing this statement */
        if (PROGRESS(numUdp) != 0) {
            for (ipnt = 0; ipnt < *npnt; ipnt++) {
                printf("        Pnt %3d: %-20s %1d %10.5f %10.5f\n", ipnt,
                       pnt[ipnt].name, pnt[ipnt].type, pnt[ipnt].x, pnt[ipnt].y);
            }
            for (iseg = 0; iseg < *nseg; iseg++) {
                printf("        Seg %3d: %-20s %1d %5d %5d\n", iseg,
                       seg[iseg].name, seg[iseg].type, seg[iseg].ibeg, seg[iseg].iend);
            }
        }

#ifdef GRAFIC
        status = plotWaffle(*npnt, pnt, *nseg, seg);
        CHECK_STATUS(plotWaffle);
#endif
    }

cleanup:

    /* delete any Parameters that were added */
    if (MODL != NULL) {
        (void) ocsmInfo(MODL, &nbrch, &npmtr, &nbody);

        for (ipmtr = npmtr; ipmtr > npmtr_save; ipmtr--) {
            (void) ocsmDelPmtr(MODL, ipmtr);
        }
    }

    /*@ignore@*/
    FREE(begline);
    FREE(endline);
    FREE(set_pmtr);
    /*@end@*/

    /* return pointers to freeable memory */
    *pnt_p = pnt;
    *seg_p = seg;

    /* close the file */
    if (fp != NULL) {
        fclose(fp);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   modifySegments - make Segments link end to end                     *
 *                                                                      *
 ************************************************************************
 */

static int
modifySegments(int     *npnt,           /* (both) number of points */
               pnt_T   *pnt[],          /* (both) array  of points */
               int     *nseg,           /* (both) number of segments */
               seg_T   *seg[],          /* (both) array  of segments */
               char    message[])       /* (out)  error message */
{
    int     status = EGADS_SUCCESS;     /* (out)  return status */

    int     ipnt, iseg, jseg;
    int     ibeg, iend, jbeg, jend, jpnt, mpnt, mseg, iter, nchange, iattr;
    double  N,     D,     s,     t,     xx,     yy,     frac, dist;
    double  N_dot, D_dot, s_dot,        xx_dot, yy_dot;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(modifySegments);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("entering modifySegments\n");
    for (ipnt = 0; ipnt < *npnt; ipnt++) {
        printf("Pnt %3d: type=%d, x=%10.5f,%10.5f, y=%10.5f,%10.5f, name=%s\n", ipnt,
               (*pnt)[ipnt].type,
               (*pnt)[ipnt].x,
               (*pnt)[ipnt].x_dot,
               (*pnt)[ipnt].y,
               (*pnt)[ipnt].y_dot,
               (*pnt)[ipnt].name);
    }
    for (iseg = 0; iseg < *nseg; iseg++) {
        printf("Seg %3d: type=%d, ibeg=%3d, iend=%3d, num=%3d, idx=%3d, name=%s\n", iseg,
               (*seg)[iseg].type,
               (*seg)[iseg].ibeg,
               (*seg)[iseg].iend,
               (*seg)[iseg].num,
               (*seg)[iseg].idx,
               (*seg)[iseg].name);
    }
#endif

    /* current size of arrays */
    mpnt = *npnt;
    mseg = *nseg;

    /* check for intersections of Lines only */
    for (jseg = 0; jseg < *nseg; jseg++) {
        if ((*seg)[jseg].type == 0) continue;

        for (iseg = jseg+1; iseg < *nseg; iseg++) {
            if ((*seg)[iseg].type == 0) continue;

            ibeg = (*seg)[iseg].ibeg;
            iend = (*seg)[iseg].iend;
            jbeg = (*seg)[jseg].ibeg;
            jend = (*seg)[jseg].iend;

            D     = ((*pnt)[iend].x     - (*pnt)[ibeg].x    ) * ((*pnt)[jbeg].y     - (*pnt)[jend].y    )
                  - ((*pnt)[jbeg].x     - (*pnt)[jend].x    ) * ((*pnt)[iend].y     - (*pnt)[ibeg].y    );
            D_dot = ((*pnt)[iend].x_dot - (*pnt)[ibeg].x_dot) * ((*pnt)[jbeg].y     - (*pnt)[jend].y    )
                  - ((*pnt)[jbeg].x_dot - (*pnt)[jend].x_dot) * ((*pnt)[iend].y     - (*pnt)[ibeg].y    )
                  + ((*pnt)[iend].x     - (*pnt)[ibeg].x    ) * ((*pnt)[jbeg].y_dot - (*pnt)[jend].y_dot)
                  - ((*pnt)[jbeg].x     - (*pnt)[jend].x    ) * ((*pnt)[iend].y_dot - (*pnt)[ibeg].y_dot);
            if (fabs(D) > EPS06) {
                N      = ((*pnt)[jbeg].x     - (*pnt)[ibeg].x    ) * ((*pnt)[jbeg].y     - (*pnt)[jend].y    )
                       - ((*pnt)[jbeg].x     - (*pnt)[jend].x    ) * ((*pnt)[jbeg].y     - (*pnt)[ibeg].y    );
                N_dot  = ((*pnt)[jbeg].x_dot - (*pnt)[ibeg].x_dot) * ((*pnt)[jbeg].y     - (*pnt)[jend].y    )
                       - ((*pnt)[jbeg].x_dot - (*pnt)[jend].x_dot) * ((*pnt)[jbeg].y     - (*pnt)[ibeg].y    )
                       + ((*pnt)[jbeg].x     - (*pnt)[ibeg].x    ) * ((*pnt)[jbeg].y_dot - (*pnt)[jend].y_dot)
                       - ((*pnt)[jbeg].x     - (*pnt)[jend].x    ) * ((*pnt)[jbeg].y_dot - (*pnt)[ibeg].y_dot);

                s     = N / D;
                s_dot = (D * N_dot - N * D_dot) / D / D;

                N     = ((*pnt)[iend].x     - (*pnt)[ibeg].x    ) * ((*pnt)[jbeg].y     - (*pnt)[ibeg].y    )
                      - ((*pnt)[jbeg].x     - (*pnt)[ibeg].x    ) * ((*pnt)[iend].y     - (*pnt)[ibeg].y    );

                t     = N / D;

                if (s > -EPS06 && s < 1+EPS06 &&
                    t > -EPS06 && t < 1+EPS06   ) {
                    xx     = (1 - s) * (*pnt)[ibeg].x     + s * (*pnt)[iend].x;
                    xx_dot = (1 - s) * (*pnt)[ibeg].x_dot + s * (*pnt)[iend].x_dot + s_dot * ((*pnt)[iend].x - (*pnt)[ibeg].x);
                    yy     = (1 - s) * (*pnt)[ibeg].y     + s * (*pnt)[iend].y;
                    yy_dot = (1 - s) * (*pnt)[ibeg].y_dot + s * (*pnt)[iend].y_dot + s_dot * ((*pnt)[iend].y - (*pnt)[ibeg].y);

                    ipnt = -1;
                    for (jpnt = 0; jpnt < *npnt; jpnt++) {
                        if (fabs(xx-(*pnt)[jpnt].x) < EPS06 &&
                            fabs(yy-(*pnt)[jpnt].y) < EPS06   ) {
                            ipnt = jpnt;
                            break;
                        }
                    }

                    if (ipnt < 0) {
                        if (*npnt+1 >= mpnt) {
                            mpnt += 10;

                            RALLOC(*pnt, pnt_T, mpnt);
                        }

                        ipnt = *npnt;

                        (*pnt)[*npnt].type  = 1;
                        (*pnt)[*npnt].x     = xx;
                        (*pnt)[*npnt].x_dot = xx_dot;
                        (*pnt)[*npnt].y     = yy;
                        (*pnt)[*npnt].y_dot = yy_dot;
                        (*pnt)[*npnt].name[0] = '\0';
                        (*npnt)++;
                    }

                    if ((ibeg != ipnt && iend != ipnt) ||
                        (jbeg != ipnt && jend != ipnt)   ) {
                        if (*nseg+2 >= mseg) {
                            mseg += 10;

                            RALLOC(*seg, seg_T, mseg);
                        }
                    }

                    if (ibeg != ipnt && iend != ipnt) {
                        (*seg)[*nseg].type    = (*seg)[iseg].type;
                        (*seg)[*nseg].ibeg    = ipnt;
                        (*seg)[*nseg].iend    = iend;
                        (*seg)[*nseg].num     = (*seg)[iseg].num;
                        (*seg)[*nseg].idx     = 0;
                        (*seg)[*nseg].name[0] = '\0';

                        (*seg)[*nseg].nattr = (*seg)[iseg].nattr;
                        (*seg)[*nseg].aname = NULL;
                        (*seg)[*nseg].avalu = NULL;

                        if ((*seg)[*nseg].nattr > 0) {
                            MALLOC((*seg)[*nseg].aname, char*, (*seg)[*nseg].nattr);
                            MALLOC((*seg)[*nseg].avalu, char*, (*seg)[*nseg].nattr);
                        }

                        for (iattr = 0; iattr < (*seg)[*nseg].nattr; iattr++) {
                            (*seg)[*nseg].aname[iattr] = NULL;
                            (*seg)[*nseg].avalu[iattr] = NULL;

                            MALLOC((*seg)[*nseg].aname[iattr], char, 80);
                            MALLOC((*seg)[*nseg].avalu[iattr], char, 80);

                            strncpy((*seg)[*nseg].aname[iattr], (*seg)[iseg].aname[iattr], 80);
                            strncpy((*seg)[*nseg].avalu[iattr], (*seg)[iseg].avalu[iattr], 80);
                        }

                        (*nseg)++;

                        (*seg)[iseg].iend = ipnt;
                    }

                    if (jbeg != ipnt && jend != ipnt) {
                        (*seg)[*nseg].type    = (*seg)[jseg].type;
                        (*seg)[*nseg].ibeg    = ipnt;
                        (*seg)[*nseg].iend    = jend;
                        (*seg)[*nseg].num     = (*seg)[jseg].num;
                        (*seg)[*nseg].idx     = 0;
                        (*seg)[*nseg].name[0] = '\0';

                        (*seg)[*nseg].nattr = (*seg)[jseg].nattr;
                        (*seg)[*nseg].aname = NULL;
                        (*seg)[*nseg].avalu = NULL;

                        if ((*seg)[*nseg].nattr > 0) {
                            MALLOC((*seg)[*nseg].aname, char*, (*seg)[*nseg].nattr);
                            MALLOC((*seg)[*nseg].avalu, char*, (*seg)[*nseg].nattr);
                        }

                        for (iattr = 0; iattr < (*seg)[*nseg].nattr; iattr++) {
                            (*seg)[*nseg].aname[iattr] = NULL;
                            (*seg)[*nseg].avalu[iattr] = NULL;

                            MALLOC((*seg)[*nseg].aname[iattr], char, 80);
                            MALLOC((*seg)[*nseg].avalu[iattr], char, 80);

                            strncpy((*seg)[*nseg].aname[iattr], (*seg)[jseg].aname[iattr], 80);
                            strncpy((*seg)[*nseg].avalu[iattr], (*seg)[jseg].avalu[iattr], 80);
                        }

                        (*nseg)++;

                        (*seg)[jseg].iend = ipnt;
                    }
                }
            }
        }
    }

    /* break Lines at Points */
    for (iseg = 0; iseg < *nseg; iseg++) {
        if ((*seg)[iseg].type == 0) continue;

        ibeg = (*seg)[iseg].ibeg;
        iend = (*seg)[iseg].iend;

        for (ipnt = 0; ipnt < *npnt; ipnt++) {
            if ((*pnt)[ipnt].type == 0) continue;

            /* distance from Point to line */
            frac = ( ((*pnt)[ipnt].x - (*pnt)[ibeg].x) * ((*pnt)[iend].x - (*pnt)[ibeg].x)
                   + ((*pnt)[ipnt].y - (*pnt)[ibeg].y) * ((*pnt)[iend].y - (*pnt)[ibeg].y))
                 / ( ((*pnt)[iend].x - (*pnt)[ibeg].x) * ((*pnt)[iend].x - (*pnt)[ibeg].x)
                   + ((*pnt)[iend].y - (*pnt)[ibeg].y) * ((*pnt)[iend].y - (*pnt)[ibeg].y));

            if (frac < EPS06 || frac > 1-EPS06) continue;

            xx = (1-frac) * (*pnt)[ibeg].x + frac * (*pnt)[iend].x;
            yy = (1-frac) * (*pnt)[ibeg].y + frac * (*pnt)[iend].y;

            dist = sqrt( (xx - (*pnt)[ipnt].x) * (xx - (*pnt)[ipnt].x)
                       + (yy - (*pnt)[ipnt].y) * (yy - (*pnt)[ipnt].y));


            if (dist < EPS06) {

                /* make room for new Segment */
                if (*nseg+1 >= mseg) {
                    mseg += 10;

                    RALLOC(*seg, seg_T, mseg);
                }

                /* make second half */
                (*seg)[*nseg].type    = (*seg)[iseg].type;
                (*seg)[*nseg].ibeg    = ipnt;
                (*seg)[*nseg].iend    = (*seg)[iseg].iend;
                (*seg)[*nseg].num     = (*seg)[iseg].num;
                (*seg)[*nseg].idx     = 0;
                (*seg)[*nseg].name[0] = '\0';

                (*seg)[*nseg].nattr = (*seg)[iseg].nattr;
                (*seg)[*nseg].aname = NULL;
                (*seg)[*nseg].avalu = NULL;

                if ((*seg)[*nseg].nattr > 0) {
                    MALLOC((*seg)[*nseg].aname, char*, (*seg)[*nseg].nattr);
                    MALLOC((*seg)[*nseg].avalu, char*, (*seg)[*nseg].nattr);
                }

                for (iattr = 0; iattr < (*seg)[*nseg].nattr; iattr++) {
                    (*seg)[*nseg].aname[iattr] = NULL;
                    (*seg)[*nseg].avalu[iattr] = NULL;

                    MALLOC((*seg)[*nseg].aname[iattr], char, 80);
                    MALLOC((*seg)[*nseg].avalu[iattr], char, 80);

                    strncpy((*seg)[*nseg].aname[iattr], (*seg)[iseg].aname[iattr], 80);
                    strncpy((*seg)[*nseg].avalu[iattr], (*seg)[iseg].avalu[iattr], 80);
                }

                /* revise first half */
                (*seg)[iseg].iend = ipnt;

                (*nseg)++;

                /* start again at this Segment */
                iseg--;
                break;
            }
        }
    }


    /* remove the "cline" Segments */
    for (iseg = 0; iseg < *nseg; iseg++) {
        if ((*seg)[iseg].type == 0) {
            FREE((*seg)[iseg].aname);
            FREE((*seg)[iseg].avalu);

            for (jseg = iseg; jseg < *nseg-1; jseg++) {
                (*seg)[jseg].type  = (*seg)[jseg+1].type;
                (*seg)[jseg].ibeg  = (*seg)[jseg+1].ibeg;
                (*seg)[jseg].iend  = (*seg)[jseg+1].iend;
                (*seg)[jseg].num   = (*seg)[jseg+1].num;
                (*seg)[jseg].idx   = (*seg)[jseg+1].idx;
                (*seg)[jseg].nattr = (*seg)[jseg+1].nattr;
                (*seg)[jseg].aname = (*seg)[jseg+1].aname;
                (*seg)[jseg].avalu = (*seg)[jseg+1].avalu;

                strcpy((*seg)[jseg].name, (*seg)[jseg+1].name);
            }

            iseg--;
            (*nseg)--;
        }
    }

    /* check for degenerate Segments */
    for (jseg = 0; jseg < *nseg; jseg++) {
        if ((*seg)[jseg].ibeg == (*seg)[jseg].iend) {
            snprintf(message, 1024, "Segment %d is degenerate", iseg);
            status = EGADS_DEGEN;
            goto cleanup;
        }
    }

    /* assign the Segment indices */
    for (iter = 0; iter < *nseg; iter++) {
        nchange = 0;

        for (iseg = 0; iseg < *nseg; iseg++) {
            if ((*seg)[iseg].idx > 0) continue;

            for (jseg = 0; jseg < *nseg; jseg++) {
                if (jseg == iseg || (*seg)[jseg].idx <= 0) continue;

                if ((*seg)[iseg].num  == (*seg)[jseg].num  &&
                    (*seg)[iseg].ibeg == (*seg)[jseg].iend   ) {
                    (*seg)[iseg].idx = (*seg)[jseg].idx + 1;
                    nchange++;
                    break;
                }
            }
        }

        if (nchange == 0) break;
    }

#ifdef DEBUG
    printf("exiting modifySegments\n");
    for (ipnt = 0; ipnt < *npnt; ipnt++) {
        printf("Pnt %3d: type=%d, x=%10.5f,%10.5f, y=%10.5f,%10.5f, name=%s\n", ipnt,
               (*pnt)[ipnt].type,
               (*pnt)[ipnt].x,
               (*pnt)[ipnt].x_dot,
               (*pnt)[ipnt].y,
               (*pnt)[ipnt].y_dot,
               (*pnt)[ipnt].name);
    }
    for (iseg = 0; iseg < *nseg; iseg++) {
        printf("Seg %3d: type=%d, ibeg=%3d, iend=%3d, num=%3d, idx=%3d, name=%s\n", iseg,
               (*seg)[iseg].type,
               (*seg)[iseg].ibeg,
               (*seg)[iseg].iend,
               (*seg)[iseg].num,
               (*seg)[iseg].idx,
               (*seg)[iseg].name);
    }
#endif

cleanup:
return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   getToken - get a token from a string                               *
 *                                                                      *
 ************************************************************************
 */

static int
getToken(char   *text,                  /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   sep,                    /* (in)  seperator character */
         int    maxtok,                 /* (in)  size of token */
         char   *token)                 /* (out) token */
{
    int  lentok = 0;                    /* (out) length of token (or -1 if not found) */

    int     ibeg, i, j, count, iskip;

    /* --------------------------------------------------------------- */

    token[0] = '\0';
    lentok   = 0;

    /* convert tabs to spaces */
    for (i = 0; i < STRLEN(text); i++) {
        if (text[i] == '\t') {
            text[i] = ' ';
        }
    }

    /* skip past white spaces at beginning of text */
    ibeg = 0;
    while (ibeg < STRLEN(text)) {
        if (text[ibeg] != ' ' && text[ibeg] != '\t' && text[ibeg] != '\r') break;
        ibeg++;
    }

    /* count the number of separators */
    count = 0;
    for (i = ibeg; i < STRLEN(text); i++) {
        if (text[i] == sep) {
            count++;
            for (j = i+1; j < STRLEN(text); j++) {
                if (text[j] != sep) break;
                i++;
            }
        }
    }

    if (count < nskip) {
        lentok = -1;
        goto cleanup;
    }

    /* skip over nskip tokens */
    i = ibeg;
    for (iskip = 0; iskip < nskip; iskip++) {
        while (text[i] != sep) {
            i++;
        }
        while (text[i] == sep) {
            i++;
        }
    }

    /* extract the token we are looking for */
    while (text[i] != sep && i < STRLEN(text)) {
        token[lentok++] = text[i++];
        token[lentok  ] = '\0';

        if (lentok >= maxtok-1) {
            lentok = -1;
            goto cleanup;
        }
    }

    lentok = STRLEN(token);

cleanup:
    return lentok;
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
               int    npnt_in,          /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int     status = EGADS_SUCCESS;

    int     iudp, judp;

    int     ipnt_in, ipnt, npnt, iseg, nseg, ibeg, iend, outLevel_save;
    int     nnode, nedge, nface, periodic, attrType, attrLen, iattr;
    CINT    *tempIlist;
    double  frac, fracu, fracv, z_pnt_dot, z_beg_dot, z_end_dot, data[18], uvrange[4];
    CDOUBLE *tempRlist;
    char    message[80];
    CCHAR   *tempClist;
    pnt_T   *pnt=NULL;
    seg_T   *seg=NULL;
    ego     context, *enodes, *eedges, *efaces, eface;

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

    status = EG_getContext(ebody, &context);
    CHECK_STATUS(EG_getContext);

    /* re-process the Segments (to propagate the velocities) */
    if (STRLEN(FILENAME(iudp)) > 0) {
        outLevel_save = ocsmSetOutLevel(0);
        status = processFile(context, message, &npnt, &pnt, &nseg, &seg, NumUdp, udps);
        CHECK_STATUS(processFile);
        (void) ocsmSetOutLevel(outLevel_save);

        status = modifySegments(&npnt, &pnt, &nseg, &seg, message);
        CHECK_STATUS(modifySegments);
    } else {
        status = processSegments(iudp, &npnt, &pnt, &nseg, &seg, NumUdp, udps);
        CHECK_STATUS(processSegments);

        status = modifySegments(&npnt, &pnt, &nseg, &seg, message);
        CHECK_STATUS(modifySegments);
    }

    SPLINT_CHECK_FOR_NULL(pnt);
    SPLINT_CHECK_FOR_NULL(seg);

    /* get the velocity for a Node */
    if (entType == OCSM_NODE) {
        if (npnt_in > 1) {
            printf("for OCSM_NODE, npnt_in should be 1 (not %d)\n", npnt_in);
            status = EGADS_INDEXERR;
            goto cleanup;
        } else if (entIndex < 1 || entIndex > 2*npnt) {
            printf("entIndex=%d is not between 1 and %d\n", entIndex, 2*npnt);
            status = EGADS_INDEXERR;
            goto cleanup;
        }

        status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_attributeRet(enodes[entIndex-1], "__nodes__", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        CHECK_STATUS(EG_attributeRet);

        EG_free(enodes);

        if (tempIlist[0] < npnt) {
            ipnt      = tempIlist[0];
            z_pnt_dot = 0;
        } else {
            ipnt      = tempIlist[0] - npnt;
            z_pnt_dot = DEPTH_DOT(iudp);
        }

        vels[0] = pnt[ipnt].x_dot;
        vels[1] = pnt[ipnt].y_dot;
        vels[2] = z_pnt_dot;

    } else if (entType == OCSM_EDGE) {
        if (entIndex < 1 || entIndex > 2*nseg+npnt) {
            printf("entIndex=%d is not between 1 and %d\n", entIndex, 2*nseg+npnt);
            status = EGADS_INDEXERR;
            goto cleanup;
        }

        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_attributeRet(eedges[entIndex-1], "__nodes__", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        CHECK_STATUS(EG_attributeRet);

        status = EG_getRange(eedges[entIndex-1], uvrange, &periodic);
        CHECK_STATUS(EG_getRange);

        EG_free(eedges);

        if (tempIlist[0] < npnt) {
            ibeg      = tempIlist[0];
            z_beg_dot = 0;
        } else {
            ibeg      = tempIlist[0] - npnt;
            z_beg_dot = DEPTH_DOT(iudp);
        }
        if (tempIlist[1] < npnt) {
            iend      = tempIlist[1];
            z_end_dot = 0;
        } else {
            iend      = tempIlist[1] - npnt;
            z_end_dot = DEPTH_DOT(iudp);
        }

        for (ipnt_in = 0; ipnt_in < npnt_in; ipnt_in++) {
            frac = (uvs[ipnt_in] - uvrange[0]) / (uvrange[1] - uvrange[0]);
            if (frac < 0 || frac > 1) {
                printf("frac out of range: entIndex=%d, uvs[%d]=%12.6f, frac=%12.6f\n", entIndex, ipnt_in, uvs[ipnt_in], frac);
                exit(0);
            }

            vels[3*ipnt_in  ] = (1-frac) * pnt[ibeg].x_dot + frac * pnt[iend].x_dot;
            vels[3*ipnt_in+1] = (1-frac) * pnt[ibeg].y_dot + frac * pnt[iend].y_dot;
            vels[3*ipnt_in+2] = (1-frac) * z_beg_dot       + frac * z_end_dot;
        }

    } else if (entType == OCSM_FACE) {
        if (entIndex < 1 || entIndex > nseg) {
            printf("entIndex=%d is not between 1 and %d\n", entIndex, nseg);
            status = EGADS_INDEXERR;
            goto cleanup;
        }

        iseg = entIndex - 1;
        ibeg = seg[iseg].ibeg;
        iend = seg[iseg].iend;

        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        eface = efaces[entIndex-1];

        EG_free(efaces);

        for (ipnt_in = 0; ipnt_in < npnt_in; ipnt_in++) {
            status = EG_evaluate(eface, &(uvs[2*ipnt_in]), data);
            CHECK_STATUS(EG_evaluate);

            fracu = ((data[0]    -pnt[ibeg].x) * (pnt[iend].x-pnt[ibeg].x)
                    +(data[1]    -pnt[ibeg].y) * (pnt[iend].y-pnt[ibeg].y))
                  / ((pnt[iend].x-pnt[ibeg].x) * (pnt[iend].x-pnt[ibeg].x)
                    +(pnt[iend].y-pnt[ibeg].y) * (pnt[iend].y-pnt[ibeg].y));
            fracv = data[2] / DEPTH(iudp);

            if        (fracu < -EPS06 || fracu > 1+EPS06) {
                printf("frac out of range: iseg=%d, fracu=%12.6f\n", iseg, fracu);
                printf("data    :%12.6f,%12.6f,%12.6f\n", data[0], data[1], data[2]);
                printf("ibeg=%3d:%12.6f,%12.6f,%12.6f\n", ibeg, pnt[ibeg].x, pnt[ibeg].y, 0.);
                printf("iend=%3d:%12.6f,%12.6f,%12.6f\n", iend, pnt[iend].x, pnt[iend].y, DEPTH(iudp));
                exit(0);
            } else if (fracv < -EPS06 || fracv > 1+EPS06) {
                printf("frac out of range: iseg=%d, fracv=%12.6f\n", iseg, fracv);
                printf("data    :%12.6f,%12.6f,%12.6f\n", data[0], data[1], data[2]);
                printf("ibeg=%3d:%12.6f,%12.6f,%12.6f\n", ibeg, pnt[ibeg].x, pnt[ibeg].y, 0.);
                printf("iend=%3d:%12.6f,%12.6f,%12.6f\n", iend, pnt[iend].x, pnt[iend].y, DEPTH(iudp));
                exit(0);
            }

            vels[3*ipnt_in  ] = (1-fracu) * pnt[ibeg].x_dot + fracu * pnt[iend].x_dot;
            vels[3*ipnt_in+1] = (1-fracu) * pnt[ibeg].y_dot + fracu * pnt[iend].y_dot;
            vels[3*ipnt_in+2] =                               fracv * DEPTH_DOT(iudp);
        }

    } else {
        status = EGADS_ATTRERR;
        goto cleanup;
    }

    status = EGADS_SUCCESS;

cleanup:
    FREE(pnt);

    if (seg != NULL) {
        for (iseg = 0; iseg < nseg; iseg++) {
            if (seg[iseg].aname != NULL) {
                for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
                    FREE(seg[iseg].aname[iattr]);
                }
                FREE(seg[iseg].aname);
            }
            if (seg[iseg].avalu != NULL) {
                for (iattr = 0; iattr < seg[iseg].nattr; iattr++) {
                    FREE(seg[iseg].avalu[iattr]);
                }
                FREE(seg[iseg].avalu);
            }
        }
        FREE(seg);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotWaffle - plot Points and Segments                              *
 *                                                                      *
 ************************************************************************
 */

#ifdef GRAFIC
static int
plotWaffle(int    npnt,                 /* (in)  number of Points */
           pnt_T  pnt[],                /* (in)  array  of Points */
           int    nseg,                 /* (in)  number of Segments */
           seg_T  seg[])                /* (in)  array  of Segments */
{
    int    status = EGADS_SUCCESS;

    int    io_kbd=5, io_scr=6, nline=0, nplot=0;
    int    indgr=1+2+4+16+64, ipnt, iseg;
    int    *ilin=NULL, *isym=NULL, *nper=NULL;
    float  *xplot=NULL, *yplot=NULL;

    ROUTINE(plotWaffle);

    /* --------------------------------------------------------------- */

    if (npnt == 0 && nseg == 0) return status;

    MALLOC(xplot, float, (2*nseg+npnt));
    MALLOC(yplot, float, (2*nseg+npnt));
    MALLOC(ilin,  int,   (  nseg+npnt));
    MALLOC(isym,  int,   (  nseg+npnt));
    MALLOC(nper,  int,   (  nseg+npnt));

    nplot = 0;
    nline = 0;

    /* build plot arrays for Points */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        xplot[nplot] = pnt[ipnt].x;
        yplot[nplot] = pnt[ipnt].y;
        nplot++;

        if (pnt[ipnt].type == 0) {
            ilin[nline] = 0;
            isym[nline] = GR_X;
            nper[nline] = 1;
        } else {
            ilin[nline] = 0;
            isym[nline] = GR_CIRCLE;
            nper[nline] = 1;
        }
        nline++;
    }

    /* build plot arrays for Segments */
    for (iseg = 0; iseg < nseg; iseg++) {
        xplot[nplot] = pnt[seg[iseg].ibeg].x;
        yplot[nplot] = pnt[seg[iseg].ibeg].y;
        nplot++;

        xplot[nplot] = pnt[seg[iseg].iend].x;
        yplot[nplot] = pnt[seg[iseg].iend].y;
        nplot++;

        if (seg[iseg].type == 0) {
            ilin[nline] = GR_DOTTED;
            isym[nline] = 0;
            nper[nline] = 2;
        } else {
            ilin[nline] = GR_SOLID;
            isym[nline] = 0;
            nper[nline] = 2;
        }
        nline++;
    }

    /* generate plot */
    grinit_(&io_kbd, &io_scr, "udpWaffle", strlen("udpWaffle"));

    grline_(ilin, isym, &nline,                "~x~y~ ",
            &indgr, xplot, yplot, nper, strlen("~x~y~ "));

cleanup:
    FREE(xplot);
    FREE(yplot);
    FREE(ilin );
    FREE(isym );
    FREE(nper );

    return status;
}
#endif
