/*
 ************************************************************************
 *                                                                      *
 * timPlugs -- Tool Integration Module for PLUGS                        *
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "egads.h"
#include "OpenCSM.h"
#include "common.h"
#include "tim.h"
#include "wsserver.h"

#define CINT     const int
#define CDOUBLE  const double
#define CCHAR    const char

static int       outLevel   = 1;

static int       plugsPhase1(modl_T *modl, int ibody, int npmtr, int pmtrindx[], int ncloud, double cloud[],                                         double *rmsbest);
static int       plugsPhase2(modl_T *modl, int ibody, int npmtr, int pmtrindx[], int ncloud, double cloud[], int face[], int *unclass, int *reclass, double *rmsbest);
static int       matsol(double A[], double b[], int n, double x[]);
static int       solsvd(double A[], double b[], int mrow, int ncol, double W[], double x[]);
static int       tridiag(int n, double a[], double b[], double c[], double d[], double x[]);

static int       plotPointCloud(esp_T *ESP);

typedef struct {
    int    ncloud;            /* number of points in cloud */
    int    nclass;            /* number of classified points */
    double *cloud;            /* array  of points in cloud */
    int    *face;             /* array  of Faces associated with points */

    double RMS;               /* RMS distance to cloud */

    int    npmtr;             /* number of DESPMTRs */
    int    *pmtrindx;         /* array  of DESPMTR indicies */
    double *pmtrorig;         /* array  of DESPMTR initial values */
} plugs_T;


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        void  *cloudfile)               /* (in)  name of file containing cloud */
{
    int    status=0;                    /* (out) return status */

    int     icloud, jmax, ipmtr;
    char    templine[128];
    FILE    *plot_fp;

    plugs_T *plugs=NULL;

    ROUTINE(timLoad(plugs));

    /* --------------------------------------------------------------- */

    outLevel = ocsmSetOutLevel(-1);

    /* create the plugs_T structure and initialize it */
    MALLOC(ESP->udata, plugs_T, 1);

    plugs = (plugs_T *) ESP->udata;

    /* initialize it */
    plugs->ncloud   = 0;
    plugs->nclass   = 0;
    plugs->cloud    = NULL;
    plugs->face     = NULL;
    plugs->RMS      = 1e+6;
    plugs->npmtr    = 0;
    plugs->pmtrindx = NULL;
    plugs->pmtrorig = NULL;

    /* make sure that there is a plotfile */
    SPLINT_CHECK_FOR_NULL(cloudfile);

    plot_fp = fopen(cloudfile, "r");
    if (plot_fp == NULL) {
        status = EGADS_NOTFOUND;
        SPRINT1(0, "ERROR:: cloudfile \"%s\" does not exist", (char *)cloudfile);
        goto cleanup;
    } else {
        SPRINT1(1, "Reading \"%s\"", (char *)cloudfile);
    }

   /* read the header */
    if (fscanf(plot_fp, "%d %d %s", &plugs->ncloud, &jmax, templine) != 3) {
        SPRINT0(0, "ERROR:: problem reading plotfile header");
        goto cleanup;
    }

    /* read the points and close file */
    MALLOC(plugs->cloud, double, 3*plugs->ncloud);
    MALLOC(plugs->face,  int,      plugs->ncloud);

    for (icloud = 0; icloud < plugs->ncloud; icloud++) {
        fscanf(plot_fp, "%lf %lf %lf", &plugs->cloud[3*icloud], &plugs->cloud[3*icloud+1], &plugs->cloud[3*icloud+2]);
        plugs->face[icloud] = 0;
    }

    fclose(plot_fp);

    /* add the cloud points to te display */
    plotPointCloud(ESP);

    /* set up the table of DESPMTRs */
    MALLOC(plugs->pmtrindx, int,    ESP->MODL->npmtr);
    MALLOC(plugs->pmtrorig, double, ESP->MODL->npmtr);

    for (ipmtr = 1; ipmtr <= ESP->MODL->npmtr; ipmtr++) {
        if (ESP->MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
            if (ESP->MODL->pmtr[ipmtr].nrow != 1 ||
                ESP->MODL->pmtr[ipmtr].ncol != 1   ) {
                SPRINT3(0, "ERROR:: DESPMTR %s is (%d*%d) and must be a scalar",
                        ESP->MODL->pmtr[ipmtr].name, ESP->MODL->pmtr[ipmtr].nrow, ESP->MODL->pmtr[ipmtr].ncol);
                status = -999;
                goto cleanup;
            }
            plugs->pmtrindx[plugs->npmtr] = ipmtr;
            plugs->pmtrorig[plugs->npmtr] = ESP->MODL->pmtr[ipmtr].value[0];

            SPRINT3(1, "initial DESPMTR %3d: %20s = %10.5f",
                    plugs->npmtr, ESP->MODL->pmtr[plugs->pmtrindx[plugs->npmtr]].name, plugs->pmtrorig[plugs->npmtr]);

            plugs->npmtr++;
        }
    }

    /* unset the verification flag */
    ESP->MODL->verify = 0;

    /* hold the UI when executing */
    status = 1;

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timMesg - get command, process, and return response               */
/*                                                                     */
/***********************************************************************/

int
timMesg(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        char  command[])                /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int     ibody, unclass, reclass, ipmtr;
    char    response[MAX_EXPR_LEN];

    modl_T  *MODL  =             ESP->MODL;
    plugs_T *plugs = (plugs_T *)(ESP->udata);

    ROUTINE(timMesg(plugs));

    /* --------------------------------------------------------------- */

    /* "phase1" */
    if (strncmp(command, "phase1|", 7) == 0) {
        ibody = MODL->nbody;

        status = plugsPhase1(MODL, ibody,
                             plugs->npmtr, plugs->pmtrindx,
                             plugs->ncloud, plugs->cloud,
                             &(plugs->RMS));

        SPRINT1(1, "\nAt end of phase1: RMS = %12.4e", plugs->RMS);
        for (ipmtr = 0; ipmtr < plugs->npmtr; ipmtr++) {
            SPRINT5(1, "%2d %3d %20s %12.6f (%12.6f)", ipmtr, plugs->pmtrindx[ipmtr],
                    MODL->pmtr[plugs->pmtrindx[ipmtr]].name,
                    MODL->pmtr[plugs->pmtrindx[ipmtr]].value[0],
                    plugs->pmtrorig[ipmtr]);
        }

        if (status < EGADS_SUCCESS) {
            snprintf(response, MAX_EXPR_LEN, "timMesg|plugs|phase1|%10.4e|ERROR:: %d detected", plugs->RMS, status);
        } else {
            snprintf(response, MAX_EXPR_LEN, "timMesg|plugs|phase1|%10.4e|%d", plugs->RMS, status);
        }

        tim_bcst("plugs", response);

    /* "phase2" */
    } else if (strncmp(command, "phase2|", 7) == 0) {
        ibody = MODL->nbody;

        status = plugsPhase2(MODL, ibody,
                             plugs->npmtr, plugs->pmtrindx,
                             plugs->ncloud, plugs->cloud, plugs->face,
                             &unclass, &reclass,
                             &(plugs->RMS));

        SPRINT1(1, "\nAt end of phase2: RMS = %12.4e", plugs->RMS);
        for (ipmtr = 0; ipmtr < plugs->npmtr; ipmtr++) {
            SPRINT5(1, "%2d %3d %20s %12.6f (%12.6f)", ipmtr, plugs->pmtrindx[ipmtr],
                    MODL->pmtr[plugs->pmtrindx[ipmtr]].name,
                    MODL->pmtr[plugs->pmtrindx[ipmtr]].value[0],
                    plugs->pmtrorig[ipmtr]);
        }

        if (status < EGADS_SUCCESS) {
            snprintf(response, MAX_EXPR_LEN, "timMesg|plugs|phase2|%10.4e|%d|%dERROR:: %d detected", plugs->RMS, unclass, reclass, status);
        } else {
            snprintf(response, MAX_EXPR_LEN, "timMesg|plugs|phase2|%10.4e|%d|%d|%d", plugs->RMS, unclass, reclass, status);
        }

        tim_bcst("plugs", response);

    /* "draw|" */
    } else if (strncmp(command, "draw|", 5) == 0) {
        plotPointCloud(ESP);

        tim_bcst("plugs", "timMesg|plugs|draw");
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    plugs_T *plugs = (plugs_T *)(ESP->udata);

    ROUTINE(timSave(plugs));

    /* --------------------------------------------------------------- */

    if (plugs == NULL) {
        goto cleanup;
    }

    /* set the verification flag */
    ESP->MODL->verify = 1;

    /* remove the allocated data from the PLUGS structure */
    FREE(plugs->cloud   );
    FREE(plugs->face    );
    FREE(plugs->pmtrindx);
    FREE(plugs->pmtrorig);

    FREE(ESP->udata);

    tim_bcst("plugs", "timSave|plugs|");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(esp_T *ESP,                     /* (in)  pointer to ESP structure */
/*@unused@*/int   unload)               /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int     ipmtr;
    plugs_T *plugs = (plugs_T *)(ESP->udata);

    ROUTINE(timQuit(plugs));

    /* --------------------------------------------------------------- */

    if (plugs == NULL) {
        goto cleanup;
    }

    /* return all DESPMTRs to their original values */
    if (plugs->pmtrindx != NULL && plugs->pmtrorig != NULL) {
        for (ipmtr = 0; ipmtr < plugs->npmtr; ipmtr++) {
            if (ESP->MODL                                     != NULL &&
                ESP->MODL->pmtr                               != NULL &&
                ESP->MODL->pmtr[plugs->pmtrindx[ipmtr]].value != NULL   ) {
                ESP->MODL->pmtr[plugs->pmtrindx[ipmtr]].value[0] = plugs->pmtrorig[ipmtr];
            }
        }
    }

    /* remove the allocated data from the PLUGS structure */
    FREE(plugs->cloud   );
    FREE(plugs->face    );
    FREE(plugs->pmtrindx);
    FREE(plugs->pmtrorig);

    FREE(ESP->udata);

    tim_bcst("plugs", "timQuit|plugs|");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plugsPhase1 - phase 1 of PLUGS                                    */
/*                                                                     */
/***********************************************************************/

static int
plugsPhase1(modl_T *MODL,               /* (in)  pointer to MODL */
            int    ibody,               /* (in)  Body index (bias-1) */
            int    npmtr,               /* (in)  number  of DESPMTRs */
            int    pmtrindx[],          /* (in)  indices of DESPMTRs */
            int    ncloud,              /* (in)  number of points in cloud */
            double cloud[],             /* (in)  array  of points in cloud */
            double *rmsbest)            /* (out) best rms distance to cloud */
{
    int    status = SUCCESS;            /* return status */

    int     icloud, ipmtr, jpmtr, inode, builtTo, nbody, ierr, iter, niter=20, idone, oldOutLevel;
    double  bbox_cloud[6], bbox_modl[6], vel[3], qerr[6], rms, lambda, value, dot, lbound, ubound, dmax;
    double  *ajac=NULL, *ajtj=NULL, *ajtq=NULL, *delta=NULL, *pmtrbest=NULL, *W=NULL;
    clock_t old_time, new_time;

    ROUTINE(plugsPhase1);

    /* --------------------------------------------------------------- */

    old_time = clock();

    *rmsbest = 0;

    assert(ncloud>0);

    /* if there are no DESPMTRs, simply return */
    if (npmtr == 0) {
        SPRINT0(1, "Phase1 will be skipped because npmtr=0");
        goto cleanup;
    }

    /* in phase1, the DESPMTRs are changed (as little as possible) to match
       the bounding box of the cloud.  hence there are 6 errors that
       are being minimized (xmin, ymin, zmin, xnax, ymax, zmax) */

    /* get needed arrays */
    MALLOC(ajac,     double,     6*npmtr);    // jacobian
    MALLOC(ajtj,     double, npmtr*npmtr);    // Jt * J
    MALLOC(ajtq,     double,       npmtr);    // Jt * Q
    MALLOC(delta,    double,       npmtr);    // proposed change in DESPMTRs
    MALLOC(pmtrbest, double,       npmtr);    // DESPMTRs associated with best
    MALLOC(W,        double,       npmtr);    // singular values

    /* find the bounding box of the cloud */
    bbox_cloud[0] = cloud[0];
    bbox_cloud[1] = cloud[1];
    bbox_cloud[2] = cloud[2];
    bbox_cloud[3] = cloud[0];
    bbox_cloud[4] = cloud[1];
    bbox_cloud[5] = cloud[2];

    for (icloud = 0; icloud < ncloud; icloud++) {
        if (cloud[3*icloud  ] < bbox_cloud[0]) bbox_cloud[0] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] < bbox_cloud[1]) bbox_cloud[1] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] < bbox_cloud[2]) bbox_cloud[2] = cloud[3*icloud+2];
        if (cloud[3*icloud  ] > bbox_cloud[3]) bbox_cloud[3] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] > bbox_cloud[4]) bbox_cloud[4] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] > bbox_cloud[5]) bbox_cloud[5] = cloud[3*icloud+2];
    }

    SPRINT3(1, "bbox of cloud: %10.5f %10.5f %10.5f",   bbox_cloud[0], bbox_cloud[1], bbox_cloud[2]);
    SPRINT3(1, "               %10.5f %10.5f %10.5f\n", bbox_cloud[3], bbox_cloud[4], bbox_cloud[5]);

    /* compute initial objective function */
    bbox_modl[0] = +HUGEQ;
    bbox_modl[1] = +HUGEQ;
    bbox_modl[2] = +HUGEQ;
    bbox_modl[3] = -HUGEQ;
    bbox_modl[4] = -HUGEQ;
    bbox_modl[5] = -HUGEQ;

    for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
        if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
            bbox_modl[0] = MODL->body[ibody].node[inode].x;
        }
        if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
            bbox_modl[1] = MODL->body[ibody].node[inode].y;
        }
        if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
            bbox_modl[2] = MODL->body[ibody].node[inode].z;
        }
        if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
            bbox_modl[3] = MODL->body[ibody].node[inode].x;
        }
        if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
            bbox_modl[4] = MODL->body[ibody].node[inode].y;
        }
        if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
            bbox_modl[5] = MODL->body[ibody].node[inode].z;
        }
    }

    SPRINT3(1, "bbox of MODL:  %10.5f %10.5f %10.5f",   bbox_modl[0], bbox_modl[1], bbox_modl[2]);
    SPRINT3(1, "               %10.5f %10.5f %10.5f\n", bbox_modl[3], bbox_modl[4], bbox_modl[5]);

    /* compute the initial errors */
    rms = 0;
    for (ierr = 0; ierr < 6; ierr++) {
        qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
        rms       += qerr[ierr] * qerr[ierr];
    }
    rms = sqrt(rms / 6);

    /* print the initial rms and DESPMTRs */
    SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", -1, rms);
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValu);

        SPRINT1x(1, " %10.5f", pmtrbest[ipmtr]);
    }
    SPRINT0(1, " ");

    /* initialize for levenberg-marquardt */
    *rmsbest = rms;
    lambda   = 1;

    /* levenberg-marquardt iterations */
    for (iter = 0; iter < niter; iter++) {

        /* find the bounding box of the MODL and its velocities */
        bbox_modl[0] = +HUGEQ;
        bbox_modl[1] = +HUGEQ;
        bbox_modl[2] = +HUGEQ;
        bbox_modl[3] = -HUGEQ;
        bbox_modl[4] = -HUGEQ;
        bbox_modl[5] = -HUGEQ;

        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmSetVelD(MODL, 0,               0, 0, 0.0);
            CHECK_STATUS(ocsmSetVelD);

            status = ocsmSetVelD(MODL, pmtrindx[ipmtr], 1, 1, 1.0);
            CHECK_STATUS(ocsmSetVelD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);
            CHECK_STATUS(ocsmBuild);

            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
                    bbox_modl[0] = MODL->body[ibody].node[inode].x;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[        ipmtr] = vel[0];
                }
                if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
                    bbox_modl[1] = MODL->body[ibody].node[inode].y;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[  npmtr+ipmtr] = vel[1];
                }
                if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
                    bbox_modl[2] = MODL->body[ibody].node[inode].z;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[2*npmtr+ipmtr] = vel[2];
                }
                if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
                    bbox_modl[3] = MODL->body[ibody].node[inode].x;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[3*npmtr+ipmtr] = vel[0];
                }
                if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
                    bbox_modl[4] = MODL->body[ibody].node[inode].y;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[4*npmtr+ipmtr] = vel[1];
                }
                if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
                    bbox_modl[5] = MODL->body[ibody].node[inode].z;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[5*npmtr+ipmtr] = vel[2];
                }
            }
        }

        /* compute the errors */
        for (ierr = 0; ierr < 6; ierr++) {
            qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
        }

        /* compute Jt * J and Jt * Q */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                ajtj[ipmtr*npmtr+jpmtr] = 0;
                for (ierr = 0; ierr < 6; ierr++) {
                    ajtj[ipmtr*npmtr+jpmtr] += ajac[ierr*npmtr+ipmtr] * ajac[ierr*npmtr+jpmtr];
                }
            }
            ajtj[ipmtr*npmtr+ipmtr] *= (1.0 + lambda);
            ajtq[ipmtr] = 0;
            for (ierr = 0; ierr < 6; ierr++) {
                ajtq[ipmtr] -= qerr[ierr] * ajac[ierr*npmtr+ipmtr];
            }
        }

        /* solve for potential change.  this is done with SVD because the
           matrix will be singular for any DESPMTR that does not (currently)
           affect the errors */
        status = solsvd(ajtj, ajtq, npmtr, npmtr, W, delta);
        CHECK_STATUS(solsvd);

        /* if all the delta's are small, no sense in any more iterations */
        dmax = 0;
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            if (fabs(delta[ipmtr]) > dmax) dmax = fabs(delta[ipmtr]);
        }

        if (dmax < EPS06) {
            SPRINT0(1, "maximum delta is small, so no more iterations");
            break;
        }

        /* temporarily apply change */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            value = pmtrbest[ipmtr] + delta[ipmtr];

            status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
            CHECK_STATUS(ocsmGetBnds);

            if (value < lbound) value = lbound;
            if (value > ubound) value = ubound;

            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);
        }

        /* rebuild and evaluate the new objective function */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);
        if (status < SUCCESS) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                SPRINT3(1, "error  DESPMTR %3d: %20s = %10.5f",
                        ipmtr, MODL->pmtr[pmtrindx[ipmtr]].name, MODL->pmtr[pmtrindx[ipmtr]].value[0]);
            }
        }
        CHECK_STATUS(ocsmBuild);

        bbox_modl[0] = +HUGEQ;
        bbox_modl[1] = +HUGEQ;
        bbox_modl[2] = +HUGEQ;
        bbox_modl[3] = -HUGEQ;
        bbox_modl[4] = -HUGEQ;
        bbox_modl[5] = -HUGEQ;

        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
                bbox_modl[0] = MODL->body[ibody].node[inode].x;
            }
            if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
                bbox_modl[1] = MODL->body[ibody].node[inode].y;
            }
            if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
                bbox_modl[2] = MODL->body[ibody].node[inode].z;
            }
            if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
                bbox_modl[3] = MODL->body[ibody].node[inode].x;
            }
            if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
                bbox_modl[4] = MODL->body[ibody].node[inode].y;
            }
            if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
                bbox_modl[5] = MODL->body[ibody].node[inode].z;
            }
        }

        rms = 0;
        for (ierr = 0; ierr < 6; ierr++) {
            qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
            rms       += qerr[ierr] * qerr[ierr];
        }
        rms = sqrt(rms / 6);

        SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", iter, rms);
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
            CHECK_STATUS(ocsmGetValu);

            SPRINT1x(1, " %10.5f", value);
        }

        /* if this new solution is better, accept it and decrease lambda
           (making it more newton-like) */
        if (rms < *rmsbest) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
                CHECK_STATUS(ocsmGetPmtr);
            }
            *rmsbest = rms;
            lambda   = MAX(1.0e-10, lambda/2);
            SPRINT1(1, "  accepted: lambda=%10.3e", lambda);

            /* check for convergence */
            idone = 1;
            for (ierr = 0; ierr < 6; ierr++) {
                if (fabs(qerr[ierr]) > EPS06) {
                    idone = 0;
                    break;
                }
            }
            if (idone == 1) {
                SPRINT0(1, "Phase 1 converged");
                break;
            }

        /* otherwise reject the step and increase lambda
           (making it more steppest-descent-like) */
        } else {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                CHECK_STATUS(ocsmSetValuD);
            }
            lambda = MIN(1.0e+10, lambda*2);
            SPRINT1(1, "  rejected: lambda=%10.3e", lambda);
        }
    }

cleanup:
    /* clear the velocities */
    (void) ocsmSetVelD(MODL, 0, 0, 0, 0.0);

    /* return CPU time */
    new_time = clock();
    SPRINT1(1, "Phase 1 CPUtime=%9.3f sec", (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

    FREE(ajac);
    FREE(ajtj);
    FREE(ajtq);
    FREE(delta);
    FREE(pmtrbest);
    FREE(W);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plugsPhase2 - phase 2 of PLUGS                                    */
/*                                                                     */
/***********************************************************************/

static int
plugsPhase2(modl_T *MODL,               /* (in)  pointer to MODL */
            int    ibody,               /* (in)  Body index (bias-1) */
            int    npmtr,               /* (in)  number  of DESPMTRs */
            int    pmtrindx[],          /* (in)  indices of DESPMTRs */
            int    ncloud,              /* (in)  number of points in cloud */
            double cloud[],             /* (in)  array  of points in cloud */
            int    face[],              /* (in)  array  of associated Faces */
            int    *unclass,            /* (out) number of unclassified points */
            int    *reclass,            /* (out) number of reclassified points */
            double *rmsbest)            /* (out) best rms distance to cloud */
{
    int    status = SUCCESS;            /* return status */

    int     iface, icloud, ipmtr, jpmtr, npnt, ntri, itri, ip0, ip1, ip2;
    int     nerr, nvar, ivar, count, ireclass, iter, niter=50;
    int     nbody, builtTo, oldOutLevel, periodic, ibest, scaleDiag, naccept, iperturb, limit;
    int     *prevface=NULL;
    CINT    *ptype, *pindx, *tris, *tric;
    double  bbox_cloud[6], dtest, value, dot, lbound, ubound, rms, data[18], rmsperturb, valperturb;
    double  dmax, lambda, uvrange[4], dbest, scaleFact, uv_guess[2], xyz_guess[3];
    double  atotal, massprops[14];
    double  *dist=NULL, *uvface=NULL, *velface=NULL;
    double  *beta=NULL, *delta=NULL, *qerr=NULL, *qerrbest=NULL, *ajac=NULL;
    double  *atri=NULL, *btri=NULL, *ctri=NULL, *dtri=NULL, *xtri=NULL;
    double  *mat=NULL, *rhs=NULL, *xxx=NULL;
    double  *pmtrbest=NULL, *pmtrsave=NULL;
    CDOUBLE *xyz, *uv;
    clock_t old_time, new_time;

    ROUTINE(plugsPhase2);

#define JAC(I,J)  ajac[(I)*npmtr+(J)]
#define JtJ(I,J)  ajtj[(I)*nvar +(J)]
#define JtQ(I)    ajtq[(I)]

    /* --------------------------------------------------------------- */

    old_time = clock();

    nvar = 2 * ncloud + npmtr;
    nerr = 3 * ncloud;

    *rmsbest = 0;

    assert(ncloud>0);

    /* get necessary arrays */
    MALLOC(prevface, int,      ncloud    );     // Face associated on previous pass
    MALLOC(dist,     double,   ncloud    );     // min dist from each cloud to tessellation
    MALLOC(uvface,   double, 2*ncloud    );     // temp (u,v) for ocsmGetVel
    MALLOC(velface,  double, 3*ncloud    );     // temp vel   from ocsmGetVel

    MALLOC(beta,     double,   nvar      );     // optimization variables
    MALLOC(delta,    double,   nvar      );     // chg in optimization variables
    MALLOC(qerr,     double,   nerr      );     // error to be minimized
    MALLOC(qerrbest, double,   nerr      );     // error to be minimized for best so far
    MALLOC(ajac,     double,   nerr*npmtr);     // d(qerr)/d(DESPMTR)

    MALLOC(atri, double, 2*ncloud);
    MALLOC(btri, double, 2*ncloud);
    MALLOC(ctri, double, 2*ncloud);
    MALLOC(dtri, double, 2*ncloud);
    MALLOC(xtri, double, 2*ncloud);

    MALLOC(mat,  double, npmtr*npmtr);
    MALLOC(rhs,  double,       npmtr);
    MALLOC(xxx,  double,       npmtr);

    MALLOC(pmtrbest, double,  npmtr      );     // best DESPMTRs so far
    MALLOC(pmtrsave, double,  npmtr      );     // best DESPMTRS on last pass

    /* remember the DESPMTRs at the begiining of this pass */
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValuD);

        pmtrsave[ipmtr] = pmtrbest[ipmtr];
    }

    /* find the bounding box of the cloud (needed for classification) */
    bbox_cloud[0] = cloud[0];
    bbox_cloud[1] = cloud[1];
    bbox_cloud[2] = cloud[2];
    bbox_cloud[3] = cloud[0];
    bbox_cloud[4] = cloud[1];
    bbox_cloud[5] = cloud[2];

    for (icloud = 0; icloud < ncloud; icloud++) {
        if (cloud[3*icloud  ] < bbox_cloud[0]) bbox_cloud[0] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] < bbox_cloud[1]) bbox_cloud[1] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] < bbox_cloud[2]) bbox_cloud[2] = cloud[3*icloud+2];
        if (cloud[3*icloud  ] > bbox_cloud[3]) bbox_cloud[3] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] > bbox_cloud[4]) bbox_cloud[4] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] > bbox_cloud[5]) bbox_cloud[5] = cloud[3*icloud+2];
    }

    /* rebuild and evaluate the new objective function */
    nbody = 0;
    oldOutLevel = ocsmSetOutLevel(0);
    status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
    (void) ocsmSetOutLevel(oldOutLevel);
    CHECK_STATUS(ocsmBuild);

    /* remember Face associations from last time */
    for (icloud = 0; icloud < ncloud; icloud++) {
        prevface[icloud] = face[icloud];
    }

    /* only classify points that are within 0.25 of the bounding box size */
    dmax = 0.25 * MAX(MAX(bbox_cloud[3]-bbox_cloud[0],
                          bbox_cloud[4]-bbox_cloud[1]),
                          bbox_cloud[5]-bbox_cloud[2]);

    /* associate each cloud point with the closest tessellation point */
    for (icloud = 0; icloud < ncloud; icloud++) {
        face[  icloud  ] = 0;
        dist[  icloud  ] = dmax * dmax;
        beta[2*icloud  ] = 0;
        beta[2*icloud+1] = 0;
    }

    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
        status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getFaceTess);

        for (icloud = 0; icloud < ncloud; icloud++) {

            /* compare with centroid of each triangle */
            for (itri = 0; itri < ntri; itri++) {
                ip0 = tris[3*itri  ] - 1;
                ip1 = tris[3*itri+1] - 1;
                ip2 = tris[3*itri+2] - 1;

                data[0] = (xyz[3*ip0  ] + xyz[3*ip1  ] + xyz[3*ip2  ]) / 3;
                data[1] = (xyz[3*ip0+1] + xyz[3*ip1+1] + xyz[3*ip2+1]) / 3;
                data[2] = (xyz[3*ip0+2] + xyz[3*ip1+2] + xyz[3*ip2+2]) / 3;

                dtest = (cloud[3*icloud  ]-data[0]) * (cloud[3*icloud  ]-data[0])
                      + (cloud[3*icloud+1]-data[1]) * (cloud[3*icloud+1]-data[1])
                      + (cloud[3*icloud+2]-data[2]) * (cloud[3*icloud+2]-data[2]);
                if (dtest < dist[icloud]-EPS06*dmax) {
                    face[  icloud  ] = iface;
                    dist[  icloud  ] = dtest;
                    beta[2*icloud  ] = (uv[2*ip0  ] + uv[2*ip1  ] + uv[2*ip2  ]) / 3;
                    beta[2*icloud+1] = (uv[2*ip0+1] + uv[2*ip1+1] + uv[2*ip2+1]) / 3;
                }
            }
        }
    }

    /* each Face must have its share of cloud points, defined as the larger of
       (arbitrarily) a tenth of the fractional area (as compared with the whole Body's area
       or (arbitrarily) 5, whichever is larger.  for Faces with fewer points,
       reclassify the cloud points that are closest to the center of the Face */
    status = EG_getMassProperties(MODL->body[ibody].ebody, massprops);
    CHECK_STATUS(EG_getMassProperties);

    atotal = massprops[1];

    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
        count = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            if (face[icloud] == iface) count++;
        }

        status = EG_getMassProperties(MODL->body[ibody].face[iface].eface, massprops);
        CHECK_STATUS(EG_getMassProperties);

        limit = 0.10 * ncloud * massprops[1] / atotal;
        if (limit < 5) limit = 5;

        for (ireclass = count; ireclass < limit; ireclass++) {
            status = EG_getRange(MODL->body[ibody].face[iface].eface, uvrange, &periodic);
            CHECK_STATUS(EG_getRange);

            uvrange[0] = (uvrange[0] + uvrange[1]) / 2;
            uvrange[1] = (uvrange[2] + uvrange[3]) / 2;

            status = EG_evaluate(MODL->body[ibody].face[iface].eface, uvrange, data);
            CHECK_STATUS(EG_evaluate);

            ibest = -1;
            dbest = HUGEQ;

            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] == iface) continue;

                dtest = (cloud[3*icloud  ]-data[0]) * (cloud[3*icloud  ]-data[0])
                      + (cloud[3*icloud+1]-data[1]) * (cloud[3*icloud+1]-data[1])
                      + (cloud[3*icloud+2]-data[2]) * (cloud[3*icloud+2]-data[2]);
                if (dtest < dbest) {
                    ibest = icloud;
                    dbest = dtest;
                }
            }

            face[  ibest  ] = iface;
            dist[  ibest  ] = dbest;
            beta[2*ibest  ] = uvrange[0];
            beta[2*ibest+1] = uvrange[1];

            SPRINT2(1, "WARNING:: reclassifying cloud point %5d to be associated with Face %d", ibest, iface);
        }
    }

    /* report the number of cloud points associated with each Face */
    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
        count = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            if (face[icloud] == iface) count++;
        }
        SPRINT2(1, "Face %3d has %5d cloud points", iface, count);
    }

    *unclass = 0;
    for (icloud = 0; icloud < ncloud; icloud++) {
        if (face[icloud] <= 0) (*unclass)++;
    }
    SPRINT1(1, "Unclassified %5d cloud points", *unclass);

    for (icloud = 0; icloud < ncloud; icloud++) {
        dist[icloud] = sqrt(dist[icloud]);
    }

    /* if the faceID associated with all cloud points are the same as
       the previous pass, we are done */
    *reclass = 0;
    for (icloud = 0; icloud < ncloud; icloud++) {
        if (face[icloud] != prevface[icloud]) {
            if (*reclass < 20) {
                SPRINT6(1, "    cloud point %5d (%10.5f, %10.5f, %10.5f) was reclassified (%3d to %3d)",
                        icloud, cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2], prevface[icloud], face[icloud]);
            }
            (*reclass)++;
        }
    }
    if (*reclass >= 20) {
        SPRINT1(1, "    ... too many to list (%d total)", *reclass);
    }
    /* compute the errors and the rms */
    rms = 0;
    for (icloud = 0; icloud < ncloud; icloud++) {
        iface = face[icloud];

        if (iface <= 0) {
            qerr[3*icloud  ] = 0;
            qerr[3*icloud+1] = 0;
            qerr[3*icloud+2] = 0;
            continue;
        }

        status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
        CHECK_STATUS(EG_evaluate);

        qerr[3*icloud  ] = cloud[3*icloud  ] - data[0];
        qerr[3*icloud+1] = cloud[3*icloud+1] - data[1];
        qerr[3*icloud+2] = cloud[3*icloud+2] - data[2];

        rms += qerr[3*icloud  ] * qerr[3*icloud  ]
            +  qerr[3*icloud+1] * qerr[3*icloud+1]
            +  qerr[3*icloud+2] * qerr[3*icloud+2];
    }

    rms = sqrt(rms / (3 * ncloud));

    /* print the rms and DESPMTRs */
    SPRINT2x(1, "\niter=%3d, rms=%10.3e, DESPMTRs=", -1, rms);
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValu);

        SPRINT1x(1, " %10.5f", pmtrbest[ipmtr]);
    }
    SPRINT0(1, " ");

    if (*unclass == 0 && *reclass == 0) {
        SPRINT0(1, "\n    Phase2 passes converged because points are classified same as previous pass\n");
        status = 0;
        goto cleanup;
    }

    /* save initial DESPMTRs as best so far for this pass */
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(beta[2*ncloud+ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValu);
    }
    *rmsbest = rms;

    /* initialize for levenberg-marquardt */
    lambda    = 1;
    scaleDiag = 0;
    scaleFact = 1;
    naccept   = 0;

    /* levenberg-marquardt iterations */
    for (iter = 0; iter < niter; iter++) {
        if (scaleDiag == 0) {
            for (icloud = 0; icloud < ncloud; icloud++) {
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    JAC(3*icloud  ,ipmtr) = 0;
                    JAC(3*icloud+1,ipmtr) = 0;
                    JAC(3*icloud+2,ipmtr) = 0;
                }
            }

            /* find the sensitivities of the errors w.r.t. the DESPMTRs */
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmSetVelD(MODL, 0,               0, 0, 0.0);
                CHECK_STATUS(ocsmSetVelD);

                status = ocsmSetVelD(MODL, pmtrindx[ipmtr], 1, 1, 1.0);
                CHECK_STATUS(ocsmSetVelD);

                nbody = 0;
                oldOutLevel = ocsmSetOutLevel(0);
                status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
                (void) ocsmSetOutLevel(oldOutLevel);
                CHECK_STATUS(ocsmBuild);

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    count = 0;
                    for (icloud = 0; icloud < ncloud; icloud++) {
                        if (face[icloud] == iface) {
                            uvface[2*count  ] = beta[2*icloud  ];
                            uvface[2*count+1] = beta[2*icloud+1];
                            count++;
                        }
                    }

                    status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, count, uvface, velface);
                    CHECK_STATUS(ocsmGetVel);

                    count = 0;
                    for (icloud = 0; icloud < ncloud; icloud++) {
                        if (face[icloud] == iface) {
                            JAC(3*icloud  ,ipmtr) = velface[3*count  ];
                            JAC(3*icloud+1,ipmtr) = velface[3*count+1];
                            JAC(3*icloud+2,ipmtr) = velface[3*count+2];
                            count++;
                        }
                    }
                }
            }

            /* initialize matrices */
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                    mat[ipmtr*npmtr+jpmtr] = 0;
                }
                rhs[ipmtr] = 0;
            }

            /* fill in matricies */
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface = face[icloud];

                /* tridiagonal matrix for (U,V)s */
                if (iface <= 0) {
                    atri[2*icloud  ] = 0;
                    btri[2*icloud  ] = 1;
                    ctri[2*icloud  ] = 0;
                    dtri[2*icloud  ] = 0;

                    atri[2*icloud+1] = 0;
                    btri[2*icloud+1] = 1;
                    ctri[2*icloud+1] = 0;
                    dtri[2*icloud+1] = 0;
                } else {
                    status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
                    CHECK_STATUS(EG_evaluate);

                    atri[2*icloud  ] =  0;
                    btri[2*icloud  ] = (data[3] * data[3] + data[4] * data[4] + data[5] * data[5]) * (1 + lambda);
                    ctri[2*icloud  ] =  data[3] * data[6] + data[4] * data[7] + data[5] * data[8];
                    dtri[2*icloud  ] =  data[3] * qerr[3*icloud  ]
                                     +  data[4] * qerr[3*icloud+1]
                                     +  data[5] * qerr[3*icloud+2];

                    atri[2*icloud+1] =  data[3] * data[6] + data[4] * data[7] + data[5] * data[8];
                    btri[2*icloud+1] = (data[6] * data[6] + data[7] * data[7] + data[8] * data[8]) * (1 + lambda);
                    ctri[2*icloud+1] =  0;
                    dtri[2*icloud+1] =  data[6] * qerr[3*icloud  ]
                                     +  data[7] * qerr[3*icloud+1]
                                     +  data[8] * qerr[3*icloud+2];
                }

                /* full matrix for DESPMTRs */
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                        mat[ipmtr*npmtr+jpmtr] += JAC(3*icloud  ,ipmtr) * JAC(3*icloud  ,jpmtr)
                                                + JAC(3*icloud+1,ipmtr) * JAC(3*icloud+1,jpmtr)
                                                + JAC(3*icloud+2,ipmtr) * JAC(3*icloud+2,jpmtr);
                    }

                    rhs[ipmtr] += JAC(3*icloud  ,ipmtr) * qerr[3*icloud  ]
                                + JAC(3*icloud+1,ipmtr) * qerr[3*icloud+1]
                                + JAC(3*icloud+2,ipmtr) * qerr[3*icloud+2];
                }
            }

            /* modify diagonal of full matrix for Levenberg-Marquardt */
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                mat[ipmtr*npmtr+ipmtr] *= (1 + lambda);
            }

        /* we have just REJECTed a change, so only the diagonal must be scaled
           because of a new lambda */
        } else {
            for (icloud = 0; icloud < ncloud; icloud++) {
                btri[2*icloud  ] *= scaleFact;
                btri[2*icloud+1] *= scaleFact;
            }
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                mat[ipmtr*npmtr+ipmtr] *= scaleFact;
            }
        }

        /* solve tridiagonal system to update the (U,V)s */
        status = tridiag(2*ncloud, atri, btri, ctri, dtri, xtri);
        CHECK_STATUS(tridiag);

        /* solve the full system to update the DESPMTRs */
        status = matsol(mat, rhs, npmtr, xxx);
        CHECK_STATUS(matsol);

        /* assemble the solutions into the delta array */
        for (icloud = 0; icloud < ncloud; icloud++) {
            delta[2*icloud  ] = xtri[2*icloud  ];
            delta[2*icloud+1] = xtri[2*icloud+1];
        }
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            delta[2*ncloud+ipmtr] = xxx[ipmtr];
        }

        /* check for no change in DESPMTRs (since it is unlikely that rms will
           ever go to zero) */
        dmax = 0;
        for (ivar = 0; ivar < 2*ncloud+npmtr; ivar++) {
            if (fabs(delta[ivar]) > dmax) {
                dmax = fabs(delta[ivar]);
            }
        }
        if (dmax < EPS06 && lambda <= 1) {
            SPRINT4(1, "    Pass converged, dmax=%10.3e,     rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                    dmax, *rmsbest, *reclass, *unclass);
            status = 1;
            goto cleanup;
        }

        /* temporarily apply change */
        for (icloud = 0; icloud < ncloud; icloud++) {
            if (face[icloud] > 0) {
                beta[2*icloud  ] += delta[2*icloud  ];
                beta[2*icloud+1] += delta[2*icloud+1];
            }
        }
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            value = pmtrbest[ipmtr] + delta[2*ncloud+ipmtr];

            status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
            CHECK_STATUS(ocsmGetBnds);

            if (value < lbound) value = lbound;
            if (value > ubound) value = ubound;
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);
        }

        /* rebuild and evaluate the new objective function */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);

        /* compute the errors and the rms */
        if (status != SUCCESS) {
            rms = 1 + *rmsbest;
        } else {
            rms = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                qerrbest[3*icloud  ] = qerr[3*icloud  ];
                qerrbest[3*icloud+1] = qerr[3*icloud+1];
                qerrbest[3*icloud+2] = qerr[3*icloud+2];

                iface = face[icloud];

                if (iface <= 0) {
                    qerr[3*icloud  ] = 0;
                    qerr[3*icloud+1] = 0;
                    qerr[3*icloud+2] = 0;
                    continue;
                }

                status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
                CHECK_STATUS(EG_evaluate);

                qerr[3*icloud  ] = cloud[3*icloud  ] - data[0];
                qerr[3*icloud+1] = cloud[3*icloud+1] - data[1];
                qerr[3*icloud+2] = cloud[3*icloud+2] - data[2];

                rms += qerr[3*icloud  ] * qerr[3*icloud  ]
                    +  qerr[3*icloud+1] * qerr[3*icloud+1]
                    +  qerr[3*icloud+2] * qerr[3*icloud+2];
            }

            rms = sqrt(rms / (3 * ncloud));
        }

        /* print the rms and DESPMTRs */
        SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", iter, rms);
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
            CHECK_STATUS(ocsmGetValu);

            SPRINT1x(1, " %10.5f", value);
        }

        /* if this new solution is better, accept it and decrease lambda
           (making it more newton-like) */
        if (rms < *rmsbest) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
                CHECK_STATUS(ocsmGetPmtr);
            }
            *rmsbest  = rms;
            scaleDiag = 0;
            naccept++;
            lambda    = MAX(1.0e-10, lambda/2);
            SPRINT1(1, "  accepted: lambda=%10.3e", lambda);

        /* otherwise reject the step and increase lambda
           (making it more steepest-descent-like) */
        } else {
            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] > 0) {
                    beta[2*icloud  ] -= delta[2*icloud  ];
                    beta[2*icloud+1] -= delta[2*icloud+1];
                }

                qerr[3*icloud  ] = qerrbest[3*icloud  ];
                qerr[3*icloud+1] = qerrbest[3*icloud+1];
                qerr[3*icloud+2] = qerrbest[3*icloud+2];
            }
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                CHECK_STATUS(ocsmSetValuD);
            }
            scaleDiag  = 1;
            scaleFact  = 1 / (1 + lambda);
            lambda     = MIN(1.0e+10, lambda*2);
            scaleFact *= (1 + lambda);
            SPRINT1(1, "  rejected: lambda=%10.3e", lambda);
        }

        /* if lambda gets very big, stop iterations */
        if (lambda > 100) {
            SPRINT4(1, "    Pass has stalled, lambda=%10.3e, rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                    lambda, *rmsbest, *reclass, *unclass);
            if (naccept > 0) {
                status = 2;
                goto cleanup;
            } else {
                break;
            }
        }
    }

    if (iter >= niter-1) {
        SPRINT3(1, "    Pass ran out of iterations,          rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                *rmsbest, *reclass, *unclass);
        if (naccept > 0) {
            status = 3;
            goto cleanup;
        }
    }

    /* if the last build was a rejection, rebuild with the best */
    if (rms >= *rmsbest) {
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
            CHECK_STATUS(ocsmSetValuD);
        }

        status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
        CHECK_STATUS(ocsmSetVelD);

        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);
        if (status < SUCCESS) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                SPRINT3(1, "error  DESPMTR %3d: %20s = %10.5f",
                        ipmtr, MODL->pmtr[pmtrindx[ipmtr]].name, MODL->pmtr[pmtrindx[ipmtr]].value[0]);
            }
        }
        CHECK_STATUS(ocsmBuild);
    }

    /* if none of the LM itertions were accepted, check if perturbing a DESPMTR will allow us to restart */
    if (naccept == 0) {
        SPRINT0(1, "\n    ERROR:: no LM iterations were accepted, so checking if any perturbations will be better\n");

        /* get the baseline */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        CHECK_STATUS(ocsmBuild);
        (void) ocsmSetOutLevel(oldOutLevel);

        rms      = 0;
        *unclass = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            iface       = face[  icloud  ];
            uv_guess[0] = beta[2*icloud  ];
            uv_guess[1] = beta[2*icloud+1];

            if (iface > 0) {
                status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                CHECK_STATUS(EG_invEvaluateGuess);

                rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                     + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                     + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
            } else {
                (*unclass)++;
            }
        }
        rms = sqrt(rms / (3*(ncloud-(*unclass))));

        SPRINT1(1, "    baseline                     rms=%12.5e\n", rms);
        iperturb   = 0;
        rmsperturb = rms;
        valperturb = HUGEQ;

        /* now pertrub each of the DESPMTRs, one at a time, and compute the rms by inverse evaluations */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {

            status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
            CHECK_STATUS(ocsmGetBnds);

            /* decrease by 5% */
            if (pmtrbest[ipmtr] > 0) {
                value = MAX(pmtrbest[ipmtr]/1.05, lbound);
            } else {
                value = MIN(pmtrbest[ipmtr]/1.05, ubound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            CHECK_STATUS(ocsmBuild);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms      = 0;
            *unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    (*unclass)++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-(*unclass))));

            SPRINT4(1, "    ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/best=%10.5f", ipmtr, value, rms, rms/rmsperturb);

            if (rms < rmsperturb) {
                iperturb   = -ipmtr;
                rmsperturb = rms;
                valperturb = value;
            }

            /* increase by 5% */
            if (pmtrbest[ipmtr] > 0) {
                value = MIN(pmtrbest[ipmtr]*1.05, ubound);
            } else {
                value = MAX(pmtrbest[ipmtr]*1.05, lbound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            CHECK_STATUS(ocsmBuild);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms      = 0;
            *unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    (*unclass)++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-(*unclass))));

            SPRINT4(1, "    ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/best=%10.5f\n", ipmtr, value, rms, rms/rmsperturb);

            if (rms < rmsperturb) {
                iperturb   = +ipmtr;
                rmsperturb = rms;
                valperturb = value;
            }

            /* set back to nominal value */
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
            CHECK_STATUS(ocsmSetValuD);
        }

        if (iperturb > 0) {
            pmtrbest[+iperturb] = valperturb;

            status = ocsmSetValuD(MODL, pmtrindx[+iperturb], 1, 1, pmtrbest[+iperturb]);
            CHECK_STATUS(ocsmSetValuD);

            SPRINT2(1, "    restarting with ipmtr=%2d perturbed to %12.5f", +iperturb, pmtrbest[+iperturb]);
            status = 4;
            goto cleanup;
        } else if (iperturb < 0) {
            pmtrbest[-iperturb] = valperturb;

            status = ocsmSetValuD(MODL, pmtrindx[-iperturb], 1, 1, pmtrbest[-iperturb]);
            CHECK_STATUS(ocsmSetValuD);

            SPRINT2(1, "    restarting with ipmtr=%2d perturbed to %12.5f", -iperturb, pmtrbest[-iperturb]);
            status = 4;
            goto cleanup;
        } else {
            SPRINT0(1, "    no perturbation succeeded, so stopping");
            status = 5;
            goto cleanup;
        }

    } else {
        dmax = 0;
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            if (fabs(pmtrbest[ipmtr]-pmtrsave[ipmtr]) > dmax) {
                dmax = fabs(pmtrbest[ipmtr]-pmtrsave[ipmtr]);
            }
        }

        if (*unclass == 0 && dmax < 1.0e-5) {
            SPRINT3(1, "\n    Phase2 passes converged because maximum DESPMTR change is %10.3e, reclass=%5d, unclass=%5d\n",
                    dmax, *reclass, *unclass);
            status = 0;
            goto cleanup;
        }
    }

cleanup:
    /* clear the velocities */
    (void) ocsmSetVelD(MODL, 0, 0, 0, 0.0);

    /* return CPU time */
    new_time = clock();
    SPRINT1(1, "Phase 2 CPUtime=%9.3f sec",
            (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

    FREE(prevface);
    FREE(dist    );
    FREE(uvface  );
    FREE(velface );

    FREE(beta    );
    FREE(delta   );
    FREE(qerr    );
    FREE(qerrbest);
    FREE(ajac    );

    FREE(atri);
    FREE(btri);
    FREE(ctri);
    FREE(dtri);
    FREE(xtri);

    FREE(mat);
    FREE(rhs);
    FREE(xxx);

    FREE(pmtrbest);
    FREE(pmtrsave);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   matsol - Gaussian elimination with partial pivoting                *
 *                                                                      *
 ************************************************************************
 */

static int
matsol(double    A[],                   /* (in)  matrix to be solved (stored rowwise) */
                                        /* (out) upper-triangular form of matrix */
       double    b[],                   /* (in)  right hand side */
                                        /* (out) right-hand side after swapping */
       int       n,                     /* (in)  size of matrix */
       double    x[])                   /* (out) solution of A*x=b */
{
    int       status = SUCCESS;         /* (out) return status */

    int       ir, jc, kc, imax;
    double    amax, swap, fact;

    ROUTINE(matsol);

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
            status = OCSM_SINGULAR_MATRIX;
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


/*
 ************************************************************************
 *                                                                      *
 *   solsvd - solve  A * x = b  via singular value decomposition        *
 *                                                                      *
 ************************************************************************
 */

static int
solsvd(double A[],                      /* (in)  mrow*ncol matrix */
       double b[],                      /* (in)  mrow-vector (right-hand sides) */
       int    mrow,                     /* (in)  number of rows */
       int    ncol,                     /* (in)  number of columns */
       double W[],                      /* (out) ncol-vector of singular values */
       double x[])                      /* (out) ncol-vector (solution) */
{
    int    status = SUCCESS;            /* (out) return status */

    int    irow, jcol, i, j, k, flag, its, jj, ip1=0, nm=0;
    double *U=NULL, *V=NULL, *r=NULL, *t=NULL;
    double wmin, wmax, s, anorm, c, f, g, h, scale, xx, yy, zz;

    ROUTINE(solsvd);

    /* --------------------------------------------------------------- */

    /* this routine is an adaptation of svf.f found in the netlib
       repository.  it is a modification of a routine from the eispack
       collection, which in turn is a translation of the algol procedure svd,
       num. math. 14, 403-420(1970) by golub and reinsch.
       handbook for auto. comp., vol ii-linear algebra, 134-151(1971). */

    /* default return */
    for (jcol = 0; jcol < ncol; jcol++) {
        x[jcol] = 0;
    }

    /* check for legal size for A */
    if (ncol <= 0) {
        status = OCSM_ILLEGAL_VALUE;
        goto cleanup;
    } else if (mrow < ncol) {
        status = OCSM_ILLEGAL_VALUE;
        goto cleanup;
    }

    MALLOC(U, double, mrow*ncol);
    MALLOC(V, double, ncol*ncol);
    MALLOC(r, double, ncol     );
    MALLOC(t, double, ncol     );

    /* initializations needed to avoid clang warning */
    for (i = 0; i < ncol; i++) {
        W[i] = 0;
        r[i] = 0;
        t[i] = 0;
    }

    /* initialize U to the original A */
    for (irow = 0; irow < mrow; irow++) {
        for (jcol = 0; jcol < ncol; jcol++) {
            U[irow*ncol+jcol] = A[irow*ncol+jcol];
        }
    }

    /* decompose A into U*W*V' */
    g     = 0;
    scale = 0;
    anorm = 0;

    /* Householder reduction of U to bidiagonal form */
    for (i = 0; i < ncol; i++) {
        ip1   = i + 1;
        r[i]  = scale * g;
        g     = 0;
        s     = 0;
        scale = 0;
        if (i < mrow) {
            for (k = i; k < mrow; k++) {
                scale += fabs(U[k*ncol+i]);
            }
            if (scale != 0) {
                for (k = i; k < mrow; k++) {
                    U[k*ncol+i] /= scale;
                    s           += U[k*ncol+i] * U[k*ncol+i];
                }
                f           = U[i*ncol+i];
                g           = -FSIGN(sqrt(s), f);
                h           = f * g - s;
                U[i*ncol+i] = f - g;
                for (j = ip1; j < ncol; j++) {
                    s = 0;
                    for (k = i; k < mrow; k++) {
                        s += U[k*ncol+i] * U[k*ncol+j];
                    }
                    f = s / h;
                    for (k = i; k < mrow; k++) {
                        U[k*ncol+j] += f * U[k*ncol+i];
                    }
                }
                for (k = i; k < mrow; k++) {
                    U[k*ncol+i] *= scale;
                }
            }
        }
        W[i]  = scale  * g;
        g     = 0;
        s     = 0;
        scale = 0;
        if (i < mrow && i+1 != ncol) {
            for (k = ip1; k < ncol; k++) {
                scale += fabs(U[i*ncol+k]);
            }
            if (scale != 0) {
                for (k = ip1; k < ncol; k++) {
                    U[i*ncol+k] /= scale;
                    s           += U[i*ncol+k] * U[i*ncol+k];
                }
                f           = U[i*ncol+ip1];
                g           = -FSIGN(sqrt(s), f);
                h           = f * g - s;
                U[i*ncol+ip1] = f - g;
                for (k = ip1; k < ncol; k++) {
                    r[k] = U[i*ncol+k] / h;
                }
                for (j = ip1; j < mrow; j++) {
                    s = 0;
                    for (k = ip1; k < ncol; k++) {
                        s += U[j*ncol+k] * U[i*ncol+k];
                    }
                    for (k = ip1; k < ncol; k++) {
                        U[j*ncol+k] += s * r[k];
                    }
                }
                for (k = ip1; k < ncol; k++) {
                    U[i*ncol+k] *= scale;
                }
            }
        }
        anorm = MAX(anorm, (fabs(W[i]) + fabs(r[i])));
    }

    /* accumulation of right-hand transformations */
    for (i = ncol-1; i >= 0; i--) {
        if (i < ncol-1) {
            if (g != 0) {
                for (j = ip1; j < ncol; j++) {
                    V[j*ncol+i] = (U[i*ncol+j] / U[i*ncol+ip1]) / g; /* avoid possible underflow */
                }
                for (j = ip1; j < ncol; j++) {
                    s = 0;
                    for (k = ip1; k < ncol; k++) {
                        s += U[i*ncol+k] * V[k*ncol+j];
                    }
                    for (k = ip1; k < ncol; k++) {
                        V[k*ncol+j] += s * V[k*ncol+i];
                    }
                }
            }
            for (j = ip1; j < ncol; j++) {
                V[i*ncol+j] = 0;
                V[j*ncol+i] = 0;
            }
        }
        V[i*ncol+i] = 1;
        g           = r[i];
        ip1         = i;
    }

    /* accumulation of left-side transformations */
    for (i = MIN(mrow, ncol)-1; i >= 0; i--) {
        ip1 = i + 1;
        g = W[i];
        for (j = ip1; j < ncol; j++) {
            U[i*ncol+j] = 0;
        }
        if (g != 0) {
            g = 1 / g;
            for (j = ip1; j < ncol; j++) {
                s = 0;
                for (k = ip1; k < mrow; k++) {
                    s += U[k*ncol+i] * U[k*ncol+j];
                }
                f = (s / U[i*ncol+i]) * g;
                for (k = i; k < mrow; k++) {
                    U[k*ncol+j] += f * U[k*ncol+i];
                }
            }
            for (j = i; j < mrow; j++) {
                U[j*ncol+i] *= g;
            }
        } else {
            for (j = i; j < mrow; j++) {
                U[j*ncol+i] = 0;
            }
        }
        ++U[i*ncol+i];
    }

    /* diagonalization of the bidiagonal form */

    /* loop over singular values */
    for (k = ncol-1; k >= 0; k--) {

        /* loop over allowed iterations */
        for (its = 0; its < 30; its++) {

            /* test for splitting */
            flag = 1;
            for (ip1 = k; ip1 >= 0; ip1--) {
                nm = ip1 - 1;

                if ((double)(fabs(r[ip1]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }

                assert (nm >= 0);                 /* needed to avoid clang warning */
                assert (nm < ncol);               /* needed to avoid clang warning */

                if ((double)(fabs(W[nm]) + anorm) == anorm) break;
            }
            if (flag) {
                c = 0;
                s = 1;
                for (i = ip1; i < k+1; i++) {
                    f    = s * r[i];
                    r[i] = c * r[i];
                    if ((double)(fabs(f) + anorm) == anorm) break;
                    g    = W[i];
                    if (fabs(f) > fabs(g)) {
                        h = fabs(f) * sqrt(1 + (g/f) * (g/f));
                    } else if (fabs(g) == 0) {
                        h = 0;
                    } else {
                        h = fabs(g) * sqrt(1 + (f/g) * (f/g));
                    }
                    W[i] = h;
                    h    = 1 / h;
                    c    = g * h;
                    s    = -f * h;
                    for (j = 0; j < mrow; j++) {
                        yy           = U[j*ncol+nm];
                        zz           = U[j*ncol+i ];
                        U[j*ncol+nm] = yy * c + zz * s;
                        U[j*ncol+i ] = zz * c - yy * s;
                    }
                }
            }

            /* test for convergence */
            zz = W[k];
            if (ip1 == k) {

                /* make singular values non-negative */
                if (zz < 0) {
                    W[k] = -zz;
                    for (j = 0; j < ncol; j++) {
                        V[j*ncol+k] = -V[j*ncol+k];
                    }
                }
                break;
            }

            assert (ip1 >= 0);                    /* needed to avoid clang warning */
            assert (ip1 < ncol);                  /* needed to avoid clang warning */

            /* shift from bottom 2*2 minor */
            xx = W[ip1];
            nm = k - 1;
            yy = W[nm];
            g  = r[nm];
            h  = r[k];
            f  = ((yy - zz) * (yy + zz) + (g - h) * (g + h)) / (2 * h * yy);
            g  = sqrt(f * f + 1);
            f  = ((xx - zz) * (xx + zz) + h * ((yy / (f + FSIGN(g, f))) - h)) / xx;

            /* next QR transformation */
            c = 1;
            s = 1;
            for (j = ip1; j <= nm; j++) {
                i    = j + 1;
                g    = r[i];
                yy   = W[i];
                h    = s * g;
                g    = c * g;
                if (fabs(f) > fabs(h)) {
                    zz = fabs(f) * sqrt(1 + (h/f) * (h/f));
                } else if (fabs(h) == 0) {
                    zz = 0;
                } else {
                    zz = fabs(h) * sqrt(1 + (f/h) * (f/h));
                }
                r[j] = zz;
                c    = f / zz;
                s    = h / zz;
                f    = xx * c + g * s;
                g    = g * c - xx * s;
                h    = yy * s;
                yy  *= c;
                for (jj = 0; jj < ncol; jj++) {
                    xx           = V[jj*ncol+j];
                    zz           = V[jj*ncol+i];
                    V[jj*ncol+j] = xx * c + zz * s;
                    V[jj*ncol+i] = zz * c - xx * s;
                }
                if (fabs(f) > fabs(h)) {
                    zz = fabs(f) * sqrt(1 + (h/f) * (h/f));
                } else if (fabs(h) == 0) {
                    zz = 0;
                } else {
                    zz = fabs(h) * sqrt(1 + (f/h) * (f/h));
                }

                /* rotation can be arbitrary if zz=0 */
                W[j] = zz;
                if (zz != 0) {
                    zz = 1 / zz;
                    c  = f * zz;
                    s  = h * zz;
                }
                f  = c * g  + s * yy;
                xx = c * yy - s * g;
                for (jj = 0; jj < mrow; jj++) {
                    yy           = U[jj*ncol+j];
                    zz           = U[jj*ncol+i];
                    U[jj*ncol+j] = yy * c + zz * s;
                    U[jj*ncol+i] = zz * c - yy * s;
                }
            }
            r[ip1] = 0;
            r[k  ] = f;
            W[k  ] = xx;
        }
    }

    /* find the largest singular value (for scaling) */
    wmax = 0;
    for (jcol = 0; jcol < ncol; jcol++) {
        if (W[jcol] > wmax) {
            wmax = W[jcol];
        }
    }

    /* set all singular values less than wmin to zero */
    wmin = wmax * 1.0e-6;
    for (jcol = 0; jcol < ncol; jcol++) {
        if (W[jcol] < wmin) {
            W[jcol] = 0;
        }
    }

    /* perform the back-substitution */
    for (j = 0; j < ncol; j++) {
        s = 0;
        if (W[j] != 0) {
            for (i = 0; i < mrow; i++) {
                s += U[i*ncol+j] * b[i];
            }
            s /= W[j];
        }
        t[j] = s;
    }

    for (j = 0; j < ncol; j++) {
        s = 0;
        for (k = 0; k < ncol; k++) {
            s += V[j*ncol+k] * t[k];
        }
        x[j] = s;
    }

cleanup:
    FREE(t);
    FREE(r);
    FREE(V);
    FREE(U);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   tridiag - solve tridiaginal system                                 *
 *                                                                      *
 ************************************************************************
 */

static int
tridiag(int    n,                       /* (in)  size of system */
        double a[],                     /* (in)  sub-diagonals */
        double b[],                     /* (in)  diagonals */
        double c[],                     /* (in)  super-diagonals */
        double d[],                     /* (in)  right-hand side */
        double x[])                     /* (out) solution */
{
    int    status = SUCCESS;            /* (out) return status */

    int    i;
    double W, *p=NULL, *q=NULL;

    ROUTINE(tridiag);

    /* --------------------------------------------------------------- */

    MALLOC(p, double, n);
    MALLOC(q, double, n);

    /* forward elimination */
    p[0] = -c[0] / b[0];
    q[0] =  d[0] / b[0];

    for (i = 1; i < n; i++) {
        W     = b[i] + c[i] * p[i-1];
        p[i] = -a[i] / W;
        q[i] = (d[i] - c[i] * q[i-1]) / W;
    }

    /* final solution */
    x[n-1] = q[n-1];

    /* back substitution */
    for (i = n-2; i >= 0; i--) {
        x[i] = p[i] * x[i+1] + q[i];
    }

cleanup:
    FREE(p);
    FREE(q);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotPointCloud - add point cloud to scene graph                    *
 *                                                                      *
 ************************************************************************
 */

static int
plotPointCloud(esp_T *ESP)              /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    npnt, icloud, nitems, attrs, igprim;
    float  color[3], *verts=NULL;
    char   gpname[255];
    wvData items[20];

    wvContext *cntxt =             ESP->cntxt;
    plugs_T   *plugs = (plugs_T *)(ESP->udata);

    ROUTINE(plotPointCloud);

    /* --------------------------------------------------------------- */

    if (cntxt == NULL) {
        goto cleanup;
    }

    MALLOC(verts, float, 3*plugs->ncloud);

    /* unclassified points */
    nitems = 0;

    npnt = 0;
    for (icloud = 0; icloud < plugs->ncloud; icloud++) {
        if (plugs->face[icloud] <= 0) {
            verts[3*npnt  ] = (float) plugs->cloud[3*icloud  ];
            verts[3*npnt+1] = (float) plugs->cloud[3*icloud+1];
            verts[3*npnt+2] = (float) plugs->cloud[3*icloud+2];

            npnt++;
        }
    }

    /* vertices */
    if (npnt > 0) {
        status = wv_setData(WV_REAL32, npnt, (void*)verts, WV_VERTICES, &(items[nitems]));
        CHECK_STATUS(wv_setData);

        wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
        nitems++;

        /* point color */
        color[0] = 1;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
        CHECK_STATUS(wv_setData);

        nitems++;

        attrs = WV_ON;
    }

    /* get index of current graphic primitive (if it exists) */
    strcpy(gpname, "PlotPoints: unclassified");
    igprim = wv_indexGPrim(cntxt, gpname);

    /* make graphic primitive */
    if (igprim < 0 && npnt > 0) {
        igprim = status = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
        CHECK_STATUS(wv_addGPrim);

        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].pSize = 3.0;
    } else if (npnt > 0) {
        status = wv_modGPrim(cntxt, igprim, nitems, items);
        CHECK_STATUS(wv_modGPrim);
    } else {
        wv_removeGPrim(cntxt, igprim);
    }

    /* classified points */
    nitems = 0;

    npnt = 0;
    for (icloud = 0; icloud < plugs->ncloud; icloud++) {
        if (plugs->face[icloud] > 0) {
            verts[3*npnt  ] = (float) plugs->cloud[3*icloud  ];
            verts[3*npnt+1] = (float) plugs->cloud[3*icloud+1];
            verts[3*npnt+2] = (float) plugs->cloud[3*icloud+2];

            npnt++;
        }
    }

    /* vertices */
    if (npnt > 0) {
        status = wv_setData(WV_REAL32, npnt, (void*)verts, WV_VERTICES, &(items[nitems]));
        CHECK_STATUS(wv_setData);

        wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
        nitems++;

        /* point color */
        color[0] = 0;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
        CHECK_STATUS(wv_setData);

        nitems++;

        attrs = WV_ON;
    }

    /* get index of current graphic primitive (if it exists) */
    SPLINT_CHECK_FOR_NULL(cntxt);

    strcpy(gpname, "PlotPoints: classified");
    igprim = wv_indexGPrim(cntxt, gpname);

    /* make graphic primitive */
    if (igprim < 0 && npnt > 0) {
        igprim = status = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
        CHECK_STATUS(wv_addGPrim);

        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].pSize = 3.0;
    } else if (npnt > 0) {
        status = wv_modGPrim(cntxt, igprim, nitems, items);
        CHECK_STATUS(wv_modGPrim);
    } else {
        wv_removeGPrim(cntxt, igprim);
    }

cleanup:
    FREE(verts);

    return status;
}
