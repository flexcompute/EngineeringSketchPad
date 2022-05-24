/*
 ************************************************************************
 *                                                                      *
 * udfCompare -- compares points in tessfile and Body                   *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPINPUTBODYS 1
#define NUMUDPARGS       4

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define TESSFILE(IUDP)   ((char   *) (udps[IUDP].arg[0].val))
#define HISTFILE(IUDP)   ((char   *) (udps[IUDP].arg[1].val))
#define PLOTFILE(IUDP)   ((char   *) (udps[IUDP].arg[2].val))
#define TOLER(   IUDP)   ((double *) (udps[IUDP].arg[3].val))[0]

static char  *argNames[NUMUDPARGS] = {"tessfile", "histfile",  "plotfile",  "toler",  };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRSTRING,  ATTRSTRING,  ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,           0,           0,        };
static double argDdefs[NUMUDPARGS] = {0.,         0.,          0.,          1.0e-6,   };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* structure to hold Node info */
typedef struct {
    ego      enode;
    double   x;
    double   y;
    double   z;
} node_TT;

/* structure to hold Edge info */
typedef struct {
    ego      eedge;
    double   xmin;
    double   ymin;
    double   zmin;
    double   xmax;
    double   ymax;
    double   zmax;
} edge_TT;

/* structure to hold Face info */
typedef struct {
    ego      eface;
    double   xmin;
    double   ymin;
    double   zmin;
    double   xmax;
    double   ymax;
    double   zmax;
} face_TT;

/* declaration for gmgwCSM routines defined below */
static int      pointToBrepDist(double xyz[], int nnode, node_TT nodes[],
                                              int nedge, edge_TT edges[],
                                              int nface, face_TT faces[], double *dbest);
static int      addToHistogram(double entry, int nhist, double dhist[], int hist[]);
static int      printHistogram(FILE *fp, int nhist, double dhist[], int hist[]);

#define  EPS03     1.0e-03
#define  MAX(A,B)  (((A) < (B)) ? (B) : (A))

static void *realloc_temp=NULL;              /* used by RALLOC macro */


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
    int     status = EGADS_SUCCESS;
    int     oclass, mtype, nchild, *senses;
    int     tess_nnode, tess_nedge, tess_nface, nface, iface, nedge, iedge, nnode, inode;
    int     npnt, ipnt, ntri, itri, ndist, nbad, mbad, ibad, count, decade;
    int     idum1, idum2, idum3, idum4, idum5, idum6;
    double  data[4], dumu, dumv;
    double  xyz_in[3], dbest;
    double  distmax, distavg, distrms, tol;
    double  *bad=NULL;
    FILE    *fp_tess, *fp_hist, *fp_plot;
    ego     *enodes, *eedges, *efaces;
    node_TT *nodes=NULL;
    edge_TT *edges=NULL;
    face_TT *faces=NULL;
    ego     context, eref, topref, prev, next, *ebodys;

    int     ihist, nhist=28, hist[28];
    double  dhist[] = {1e-8, 2e-8, 5e-8,
                       1e-7, 2e-7, 5e-7,
                       1e-6, 2e-6, 5e-6,
                       1e-5, 2e-5, 5e-5,
                       1e-4, 2e-4, 5e-4,
                       1e-3, 2e-3, 5e-3,
                       1e-2, 2e-2, 5e-2,
                       1e-1, 2e-1, 5e-1,
                       1e+0, 2e+0, 5e+0,
                       1e+1};

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("tessfile(0) = %s\n", TESSFILE(0));
    printf("histfile(0) = %s\n", HISTFILE(0));
    printf("plotfile(0) = %s\n", PLOTFILE(0));
    printf("toler(   0) = %f\n", TOLER(   0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if (strlen(TESSFILE(0)) == 0) {
        printf(" udpExecute: tessfile must be specified\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (TOLER(0) < 0) {
        printf(" udpExecute: toler = %f < 0\n", TOLER(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        printf(" udpExecute: toler should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    }

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        printf(" udpExecute: expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("tessfile(%d) = %s\n", numUdp, TESSFILE(numUdp));
    printf("histfile(%d) = %s\n", numUdp, HISTFILE(numUdp));
    printf("plotfile(%d) = %s\n", numUdp, PLOTFILE(numUdp));
    printf("toler(   %d) = %f\n", numUdp, TOLER(   numUdp));
#endif

    /* make a copy of the Body (so that it does not get removed
       when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* open tessellation file */
    fp_tess = fopen(TESSFILE(numUdp), "r");
    if (fp_tess == NULL) {
        printf(" udpExecute: file \"%s\" not found\n", TESSFILE(numUdp));
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    fscanf(fp_tess, "%d %d %d", &tess_nnode, &tess_nedge, &tess_nface);

    printf("File \"%s\" contains %d Nodes, %d Edges, and %d Faces\n", TESSFILE(numUdp),
           tess_nnode, tess_nedge, tess_nface);

    /* open the other output files if specified */
    if (strlen(HISTFILE(numUdp)) > 0) {
        fp_hist = fopen(HISTFILE(numUdp), "w");
        if (fp_hist == NULL) {
            printf(" udpExecute: file \"%s\" could not be created\n", HISTFILE(numUdp));
            status = EGADS_RANGERR;
            goto cleanup;
        }
    } else {
        fp_hist = stdout;
    }

    if (strlen(PLOTFILE(numUdp)) > 0) {
        fp_plot = fopen(PLOTFILE(numUdp), "w");
        if (fp_plot == NULL) {
            printf(" udpExecute: file \"%s\" could not be created\n", PLOTFILE(numUdp));
            status = EGADS_RANGERR;
            goto cleanup;
        }
    } else {
        fp_plot = NULL;
    }

    /* make an array of the Node ego and bounding boxes */
    status = EG_getBodyTopos(*ebody, NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    MALLOC(nodes, node_TT, nnode);

    for (inode = 0; inode < nnode; inode++) {
        status = EG_getTopology(enodes[inode], &eref, &oclass, &mtype,
                                data, &nchild, &ebodys, &senses);
        CHECK_STATUS(EG_getTopology);

        nodes[inode].enode = enodes[inode];
        nodes[inode].x     = data[0];
        nodes[inode].y     = data[1];
        nodes[inode].z     = data[2];
    }

    EG_free(enodes);

    /* make an array of the Edge ego and bounding boxes */
    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    MALLOC(edges, edge_TT, nedge);

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getInfo(eedges[iedge], &oclass, &mtype, &topref, &prev, &next);
        CHECK_STATUS(EG_getInfo);

        if (mtype != DEGENERATE) {
            edges[iedge].eedge = eedges[iedge];

            status = EG_getBoundingBox(eedges[iedge], &(edges[iedge].xmin));
            CHECK_STATUS(EG_getBoundingBox);
        } else {
            edges[iedge].eedge = NULL;
        }
    }

    EG_free(eedges);

    /* make an array of the Face ego and bounding boxes */
    status = EG_getBodyTopos(*ebody, NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    MALLOC(faces, face_TT, nface);

    for (iface = 0; iface < nface; iface++) {
        faces[iface].eface = efaces[iface];

        status = EG_getBoundingBox(efaces[iface], &(faces[iface].xmin));
        CHECK_STATUS(EG_getBoundingBox);
    }

    EG_free(efaces);

    /* initialize the histogram */
    for (ihist = 0; ihist < nhist; ihist++) {
        hist[ihist] = 0;
    }

    distmax = 0;
    distavg = 0;
    distrms = 0;
    ndist   = 0;

    /* initialize the list of bad points */
    nbad = 0;
    mbad = 1000;

    MALLOC(bad, double, 4*mbad);

    printf("processing Nodes\n");

    /* process each point in TESSFILE */
    for (inode = 0; inode < tess_nnode; inode++) {
        fscanf(fp_tess, "%lf %lf %lf", &(xyz_in[0]), &(xyz_in[1]), &(xyz_in[2]));

        /* find closest location in Brep */
        status = pointToBrepDist(xyz_in, nnode, nodes, nedge, edges, nface, faces, &dbest);
        CHECK_STATUS(pointToBrepDist);

        distmax = MAX(distmax,  dbest);
        distavg =     distavg + dbest;
        distrms =     distrms + dbest * dbest;
        ndist++;

        /* add to histogram */
        status = addToHistogram(dbest, nhist, dhist, hist);
        CHECK_STATUS(addToHistogram);

        /* add to PLOTFILE (if it exists) */
        if (fp_plot != NULL && dbest > TOLER(numUdp)) {
            if (nbad >= mbad) {
                mbad += 1000;
                RALLOC(bad, double, 4*mbad);
            }

            bad[4*nbad  ] = xyz_in[0];
            bad[4*nbad+1] = xyz_in[1];
            bad[4*nbad+2] = xyz_in[2];
            bad[4*nbad+3] = dbest;
            nbad++;
        }
    }

    for (iedge = 0; iedge < tess_nedge; iedge++) {
        fscanf(fp_tess, "%d", &npnt);

        printf("processing Edge %4d (of %4d) with %5d points\n", iedge+1, tess_nedge, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fscanf(fp_tess, "%lf %lf %lf %lf", &(xyz_in[0]), &(xyz_in[1]), &(xyz_in[2]), &dumu);

            /* find closest location in Brep */
            status = pointToBrepDist(xyz_in, nnode, nodes, nedge, edges, nface, faces, &dbest);
            CHECK_STATUS(pointToBrepDist);

            distmax = MAX(distmax,  dbest);
            distavg =     distavg + dbest;
            distrms =     distrms + dbest * dbest;
            ndist++;

            /* add to histogram */
            status = addToHistogram(dbest, nhist, dhist, hist);
            CHECK_STATUS(addToHistogram);

            /* add to PLOTFILE (if it exists) */
            if (fp_plot != NULL && dbest > TOLER(numUdp)) {
                if (nbad >= mbad) {
                    mbad += 1000;
                    RALLOC(bad, double, 4*mbad);
                }

                bad[4*nbad  ] = xyz_in[0];
                bad[4*nbad+1] = xyz_in[1];
                bad[4*nbad+2] = xyz_in[2];
                bad[4*nbad+3] = dbest;
                nbad++;
            }
        }
    }

    for (iface = 0; iface < tess_nface; iface++) {
        fscanf(fp_tess, "%d %d", &npnt, &ntri);

        printf("processing Face %4d (of %4d) with %5d points\n", iface+1, tess_nface, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fscanf(fp_tess, "%lf %lf %lf %lf %lf %d %d", &(xyz_in[0]), &(xyz_in[1]), &(xyz_in[2]), &dumu, &dumv, &idum1, &idum2);

            /* find closest location in Brep */
            status = pointToBrepDist(xyz_in, nnode, nodes, nedge, edges, nface, faces, &dbest);
            CHECK_STATUS(pointToBrepDist);

            distmax = MAX(distmax,  dbest);
            distavg =     distavg + dbest;
            distrms =     distrms + dbest * dbest;
            ndist++;

            /* add to histogram */
            status = addToHistogram(dbest, nhist, dhist, hist);
            CHECK_STATUS(addToHistogram);

            /* add to PLOTFILE (if it exists) */
            if (fp_plot != NULL && dbest > TOLER(numUdp)) {
                if (nbad >= mbad) {
                    mbad += 1000;
                    RALLOC(bad, double, 4*mbad);
                }

                bad[4*nbad  ] = xyz_in[0];
                bad[4*nbad+1] = xyz_in[1];
                bad[4*nbad+2] = xyz_in[2];
                bad[4*nbad+3] = dbest;
                nbad++;
            }
        }

        for (itri = 0; itri < ntri; itri++) {
            fscanf(fp_tess, "%d %d %d %d %d %d", &idum1, &idum2, &idum3, &idum4, &idum5, &idum6);
        }
    }

    /* write the bad points */
    if (fp_plot != NULL) {
        for (decade = 0; decade < 10; decade++) {
            tol = pow(10, -decade);

            count = 0;
            for (ibad = 0; ibad < nbad; ibad++) {
                if (bad[4*ibad+3] >= tol) count++;
            }

            if (count > 0) {
                if        (decade == 0) {
                    fprintf(fp_plot, "%5d %5d Err>10^0|r\n",  count, 0);
                } else if (decade == 1) {
                    fprintf(fp_plot, "%5d %5d Err>10^-1|g\n", count, 0);
                } else if (decade == 2) {
                    fprintf(fp_plot, "%5d %5d Err>10^-2|b\n", count, 0);
                } else if (decade == 3) {
                    fprintf(fp_plot, "%5d %5d Err>10^-3|c\n", count, 0);
                } else if (decade == 4) {
                    fprintf(fp_plot, "%5d %5d Err>10^-4|m\n", count, 0);
                } else if (decade == 5) {
                    fprintf(fp_plot, "%5d %5d Err>10^-5|y\n", count, 0);
                } else {
                    fprintf(fp_plot, "%5d %5d Err>10^-%d\n",  count, 0, decade);
                }

                for (ibad = 0; ibad < nbad; ibad++) {
                    if (bad[4*ibad+3] >= tol) {
                        fprintf(fp_plot, "%15.8f %15.8f %15.8f\n",
                                bad[4*ibad], bad[4*ibad+1], bad[4*ibad+2]);
                        bad[4*ibad+3] = -1;
                    }
                }
            }
        }
    }

    /* close the files */
    if (fp_tess != NULL) fclose(fp_tess);
    if (fp_plot != NULL) fclose(fp_plot);

    /* write the HISTFILE */
    if (fp_hist != NULL) {
        printf("writing histfile \"%s\"\n", HISTFILE(numUdp));

        status = printHistogram(fp_hist, nhist, dhist, hist);

        fprintf(fp_hist, "    max dist = %12.4e\n", distmax);
        fprintf(fp_hist, "    avg dist = %12.4e\n", distavg/ndist);
        fprintf(fp_hist, "    rms dist = %12.4e\n", sqrt(distrms/ndist));

        /* close HISTFILE if not stdout */
        if (fp_hist != stdout) {
            fclose(fp_hist);
        }
    }

    /* no output value(s) */

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;
    goto cleanup;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(nodes);
    FREE(edges);
    FREE(faces);
    FREE(bad  );

    if (status != EGADS_SUCCESS) {
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


/***********************************************************************/
/*                                                                     */
/*   pointToBrepDist - find minimum distance of point to the Brep      */
/*                                                                     */
/***********************************************************************/

static int
pointToBrepDist(double  xyz[],          /* (in)  point to check */
                int     nnode,          /* (in)  number of Node info */
                node_TT nodes[],        /* (in)  array  of Node info */
                int     nedge,          /* (in)  number of Edge info */
                edge_TT edges[],        /* (in)  array  of Edge info */
                int     nface,          /* (in)  number of Face info */
                face_TT faces[],        /* (in)  array  of Face info */
                double  *dbest)         /* (out) maximum distance */
{
    int    status = EGADS_SUCCESS;

    int    inode, iedge, iface, ifbest=-1;
    int    oclass, mtype, nloop, *senses;
    double dtest, uv_out[2], xyz_out[3], data[18];
    ego    esurf, *eloops;

    ROUTINE(pointToBrepDist);

    /* --------------------------------------------------------------- */

    /* initialize */
    *dbest = 1e+20;

    /* do an exhaustive search through all Nodes whose bounding box
       contains xyz */
    for (inode = 0; inode < nnode; inode++) {

        /* skip Node if xyz is not in Node's bounding box */
        if (xyz[0] < nodes[inode].x-(*dbest) ||
            xyz[0] > nodes[inode].x+(*dbest) ||
            xyz[1] < nodes[inode].y-(*dbest) ||
            xyz[1] > nodes[inode].y+(*dbest) ||
            xyz[2] < nodes[inode].z-(*dbest) ||
            xyz[2] > nodes[inode].z+(*dbest)   ) continue;

        /* do inverse evaluation */
        xyz_out[0] = nodes[inode].x;
        xyz_out[1] = nodes[inode].y;
        xyz_out[2] = nodes[inode].z;

        /* compute distance */
        dtest = sqrt((xyz_out[0]-xyz[0]) * (xyz_out[0]-xyz[0])
                    +(xyz_out[1]-xyz[1]) * (xyz_out[1]-xyz[1])
                    +(xyz_out[2]-xyz[2]) * (xyz_out[2]-xyz[2]));

        /* save point if best so far */
        if (dtest < *dbest) {
            *dbest = dtest;
        }
    }

    /* do an exhaustive search through all Edges whose bounding box
       contains xyz */
    for (iedge = 0; iedge < nedge; iedge++) {
        if (edges[iedge].eedge == NULL) continue;

        /* skip Edge if xyz is not in Edge's bounding box */
        if (xyz[0] < edges[iedge].xmin-(*dbest) ||
            xyz[0] > edges[iedge].xmax+(*dbest) ||
            xyz[1] < edges[iedge].ymin-(*dbest) ||
            xyz[1] > edges[iedge].ymax+(*dbest) ||
            xyz[2] < edges[iedge].zmin-(*dbest) ||
            xyz[2] > edges[iedge].zmax+(*dbest)   ) continue;

        /* do inverse evaluation */
        status = EG_invEvaluate(edges[iedge].eedge,
                                xyz, uv_out, xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        /* compute distance */
        dtest = sqrt((xyz_out[0]-xyz[0]) * (xyz_out[0]-xyz[0])
                    +(xyz_out[1]-xyz[1]) * (xyz_out[1]-xyz[1])
                    +(xyz_out[2]-xyz[2]) * (xyz_out[2]-xyz[2]));

        /* save point if best so far */
        if (dtest < *dbest) {
            *dbest = dtest;
        }
    }

    /* do an exhaustive search through all Faces whose bounding box
       contains xyz */
    for (iface = 0; iface < nface; iface++) {

        /* skip Face if xyz is not in Face's bounding box */
        if (xyz[0] < faces[iface].xmin-(*dbest) ||
            xyz[0] > faces[iface].xmax+(*dbest) ||
            xyz[1] < faces[iface].ymin-(*dbest) ||
            xyz[1] > faces[iface].ymax+(*dbest) ||
            xyz[2] < faces[iface].zmin-(*dbest) ||
            xyz[2] > faces[iface].zmax+(*dbest)   ) continue;

        /* do inverse evaluation */
        status = EG_invEvaluate(faces[iface].eface,
                                xyz, uv_out, xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        /* compute distance */
        dtest = sqrt((xyz_out[0]-xyz[0]) * (xyz_out[0]-xyz[0])
                    +(xyz_out[1]-xyz[1]) * (xyz_out[1]-xyz[1])
                    +(xyz_out[2]-xyz[2]) * (xyz_out[2]-xyz[2]));

        /* save point if best so far */
        if (dtest < *dbest) {
            *dbest = dtest;
            ifbest = iface;
        }
    }

    /* if *dbest is not too close, try doing an inverse evaluation
       on the surface underlying ifbest */
    if (*dbest > 1e-3 && ifbest >= 0) {
        status = EG_getTopology(faces[ifbest].eface, &esurf, &oclass, &mtype,
                                data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_invEvaluate(esurf, xyz, uv_out, xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        /* compute distance */
        dtest = sqrt((xyz_out[0]-xyz[0]) * (xyz_out[0]-xyz[0])
                    +(xyz_out[1]-xyz[1]) * (xyz_out[1]-xyz[1])
                    +(xyz_out[2]-xyz[2]) * (xyz_out[2]-xyz[2]));

        /* save point if better that Face evaluation */
        if (dtest < 0.99*(*dbest)) {
            printf("WARNING:: point at (%12.5f %12.5f %12.5f) was obtained from underlying surface (dtest=%10.3e, dbest=%10.3e\n",
                   xyz[0], xyz[1], xyz[2], dtest, *dbest);

            *dbest = dtest;
        }
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   addToHistogram - add entry to histogram                           */
/*                                                                     */
/***********************************************************************/

static int
addToHistogram(double entry,            /* (in)  entry to add */
               int    nhist,            /* (in)  number of histogram entries */
               double dhist[],          /* (in)  histogram entries */
               int    hist[])           /* (both)histogram counts */
{
    int    status = EGADS_SUCCESS;      /* return status */

    int    ileft, imidl, irite;

    ROUTINE(addToHistogram);

    /* --------------------------------------------------------------- */

    /* binary search for correct histogram bin */
    ileft = 0;
    irite = nhist-1;

    while (irite-ileft > 1) {
        imidl = (ileft + irite) / 2;
        if (entry > dhist[imidl]) {
            ileft = imidl;
        } else {
            irite = imidl;
        }
    }

    hist[ileft]++;

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   printHistogram - print a histogram                                */
/*                                                                     */
/***********************************************************************/

static int
printHistogram(FILE   *fp,
               int    nhist,            /* (in)  number of histogram entries */
               double dhist[],          /* (in)  histogram entries */
               int    hist[])           /* (in)  histogram counts */
{
    int    status = EGADS_SUCCESS;      /* return status */

    int    ihist, ntotal, ix;
    double percent;

    ROUTINE(printHistogram);

    /* --------------------------------------------------------------- */

    /* compute and print total entries in histogram */
    ntotal = 0;
    for (ihist = 0; ihist < nhist; ihist++) {
        ntotal += hist[ihist];
    }

    /* print basic histogram */
    percent = 100.0 * (double)(hist[0]) / (double)(ntotal);
    fprintf(fp, "    %9d (%5.1f%%)                    < %8.1e   |",
            hist[0], percent, dhist[1]);
    for (ix = 0; ix < 20; ix++) {
        if (5.0*ix >= percent) break;

        if (ix%5 == 4) {
            fprintf(fp, "+");
        } else {
            fprintf(fp, "-");
        }
    }
    fprintf(fp, "\n");

    for (ihist = 1; ihist < nhist-2; ihist++) {
        percent = 100.0 * (double)(hist[ihist]) / (double)(ntotal);
        fprintf(fp, "    %9d (%5.1f%%) between %8.1e and %8.1e   |",
                hist[ihist], percent, dhist[ihist], dhist[ihist+1]);
        for (ix = 0; ix < 20; ix++) {
            if (5.0*ix >= percent) break;

            if (ix%5 == 4) {
                fprintf(fp, "+");
            } else {
                fprintf(fp, "-");
            }
        }
        fprintf(fp, "\n");
    }

    percent = 100.0 * (double)(hist[nhist-2]) / (double)(ntotal);
    fprintf(fp, "    %9d (%5.1f%%)       > %8.1e                |",
            hist[nhist-2], percent, dhist[nhist-2]);
    for (ix = 0; ix < 20; ix++) {
        if (5.0*ix >= percent) break;

        if (ix%5 == 4) {
            fprintf(fp, "+");
        } else {
            fprintf(fp, "-");
        }
    }
    fprintf(fp, "\n");

    fprintf(fp, "    %9d total\n", ntotal);

//cleanup:
    return status;
}
