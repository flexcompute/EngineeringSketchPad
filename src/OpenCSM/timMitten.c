/*
 ************************************************************************
 *                                                                      *
 * timMitten -- Tool Integration Module for micro-gloves application    *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
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

#include "egads.h"
#include "common.h"
#include "OpenCSM.h"
#include "tim.h"
#include "wsserver.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

typedef struct {
    double    xcent;
    double    ycent;
    double    zcent;
    double    xsize;
    double    ysize;
    double    zsize;

    char      bodyName[MAX_NAME_LEN];
} mitten_T;

static int mittenBuildBox(tim_T *TIM, double rotate);


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(tim_T *TIM,                     /* (in)  pointer to TIM structure */
        void  *data)                    /* (in)  component name */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    mitten_T *mitn;

    ROUTINE(timLoad(mitten));

    /* --------------------------------------------------------------- */

    /* create the mitten_T structure */
    MALLOC(TIM->udata, mitten_T, 1);

    mitn = (mitten_T *) TIM->udata;

    /* initialize the structure */
    mitn->xcent = 0;
    mitn->ycent = 0;
    mitn->zcent = 0;
    mitn->xsize = 1;
    mitn->ysize = 1;
    mitn->zsize = 1;

    if (data == NULL) {
        mitn->bodyName[0] = '\0';
    } else {
        /*@ignore@*/
        strcpy(mitn->bodyName, (char *)data);
        /*@end@*/
    }

    /* build the initial body */
    status = mittenBuildBox(TIM, 0);
    CHECK_STATUS(mittenBuildBox);

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(tim_T *TIM)                     /* (in)  pointer to TI Mstructure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    iface, igprim;
    char   gpname[255];
    FILE   *fp=NULL;

    mitten_T *mitn = (mitten_T *)(TIM->udata);

    ROUTINE(timSave(mitten));

    /* --------------------------------------------------------------- */

    /* add the mitten Body to the bottm of the .csm file */
    fp = fopen(TIM->MODL->filename, "a");
    if (fp == NULL) {
        printf("the .csm file does not exist or could not be opened\n");
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    fprintf(fp, "\n### begin Body created by mitten\n\n");
    fprintf(fp, "DESPMTR   %s:xcent   %10.5f\n", mitn->bodyName, mitn->xcent);
    fprintf(fp, "DESPMTR   %s:ycent   %10.5f\n", mitn->bodyName, mitn->ycent);
    fprintf(fp, "DESPMTR   %s:zcent   %10.5f\n", mitn->bodyName, mitn->zcent);
    fprintf(fp, "DESPMTR   %s:xsize   %10.5f\n", mitn->bodyName, mitn->xsize);
    fprintf(fp, "DESPMTR   %s:ysize   %10.5f\n", mitn->bodyName, mitn->ysize);
    fprintf(fp, "DESPMTR   %s:zsize   %10.5f\n", mitn->bodyName, mitn->zsize);
    fprintf(fp, "BOX       %s:xcent-%s:xsize/2   %s:ycent-%s:ysize/2   %s:zcent-%s:zsize/2 \\\n",
            mitn->bodyName, mitn->bodyName, mitn->bodyName, mitn->bodyName, mitn->bodyName, mitn->bodyName);
    fprintf(fp, "          %s:xsize   %s:ysize   %s:zsize\n",
            mitn->bodyName, mitn->bodyName, mitn->bodyName);
    fprintf(fp, "ATTRIBUTE _name $%s\n",    mitn->bodyName);
    fprintf(fp, "\n### end Body created by mitten\n\n");

    fclose(fp);

    /* remove the Body from the scene graph */
    for (iface = 1; iface <= 6; iface++) {
        snprintf(gpname, 254, "Body %s Face %d", mitn->bodyName, iface);

        igprim = wv_indexGPrim(TIM->ctxt, gpname);

        if (igprim > 0) {
            wv_removeGPrim(TIM->ctxt, igprim);
        }
    }

    /* free up the mitten_T structure */
    FREE(TIM->udata);

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(tim_T *TIM)                     /* (in)  pointer to TIM structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int      iface, igprim;
    char     gpname[255];
    mitten_T *mitn = (mitten_T *)(TIM->udata);

    ROUTINE(timQuit(mitten));

    /* --------------------------------------------------------------- */

    /* remove the Body from the scene graph */
    for (iface = 1; iface <= 6; iface++) {
        snprintf(gpname, 254, "Body %s Face %d", mitn->bodyName, iface);

        igprim = wv_indexGPrim(TIM->ctxt, gpname);
        printf("iface=%d, igprim=%d\n", iface, igprim);

        if (igprim > 0) {
            wv_removeGPrim(TIM->ctxt, igprim);
        }
    }

    /* free up the mitten_T structure */
    FREE(TIM->udata);

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timMesg - get command, process, and return response               */
/*                                                                     */
/***********************************************************************/

int
timMesg(tim_T *TIM,                     /* (in)  pointer to TIM structure */
        char  command[],                /* (in)  command */
        int   max_resp_len,             /* (in)  length of response */
        char  response[])               /* (out) response */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    double delta=0, fact=1, angle=0;
    char   arg1[MAX_EXPR_LEN], *pEnd;
    mitten_T *mitn = (mitten_T *)(TIM->udata);

    ROUTINE(timMesg(mitten));

    /* --------------------------------------------------------------- */

    /* "xcent|delta|" */
    if        (strncmp(command, "xcent|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) delta = strtod(arg1, &pEnd);

        mitn->xcent += delta;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "xcent");

    /* "xsize|fact|" */
    } else if (strncmp(command, "xsize|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) fact = strtod(arg1, &pEnd);

        mitn->xsize *= fact;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "xsize");

    /* "ycent|delta|" */
    } else if (strncmp(command, "ycent|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) delta = strtod(arg1, &pEnd);

        mitn->ycent += delta;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "ycent");

    /* "ysize|fact|" */
    } else if (strncmp(command, "ysize|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) fact = strtod(arg1, &pEnd);

        mitn->ysize *= fact;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "ysize");

    /* "zcent|delta|" */
    } else if (strncmp(command, "zcent|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) delta = strtod(arg1, &pEnd);

        mitn->zcent += delta;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "zcent");

    /* "zsize|fact|" */
    } else if (strncmp(command, "zsize|", 6) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) fact = strtod(arg1, &pEnd);

        mitn->zsize *= fact;

        mittenBuildBox(TIM, 0);

        snprintf(response, max_resp_len, "zsize");

    /* "rotate|angle|" */
    } else if (strncmp(command, "rotate|", 7) == 0) {

        /* extract argument */
        if (GetToken(command, 1, '|', arg1)) angle = strtod(arg1, &pEnd);

        mittenBuildBox(TIM, angle);

        snprintf(response, max_resp_len, "rotate|%10.3f", angle);

    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   mittenBuildBox - build a box and put it on the screen             */
/*                                                                     */
/***********************************************************************/

static int
mittenBuildBox(tim_T    *TIM,           /* (in)  pointer to TIM structure */
               double   angle)          /* (in)  rotation angle */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int     iface, npnt, ntri, itri, k, nitems, igprim, nseg, *segs=NULL, attrs;
    CINT    *ptype, *pindx, *tris, *tric;
    float   color[3];
    double  data[6], params[3], sinz, cosz, matrix[12];
    CDOUBLE *xyz, *uv;
    char    gpname[255];
    ego     etemp, exform, ebody, etess;
    wvData  items[20];
    mitten_T *mitn = (mitten_T *)(TIM->udata);

    ROUTINE(mittenBuildBox);

    /* --------------------------------------------------------------- */

    /* build the new body to be visualized */
    data[0] = mitn->xcent - mitn->xsize / 2;
    data[1] = mitn->ycent - mitn->ysize / 2;
    data[2] = mitn->zcent - mitn->zsize / 2;
    data[3] = mitn->xsize;
    data[4] = mitn->ysize;
    data[5] = mitn->zsize;

    status = EG_makeSolidBody(TIM->MODL->context, BOX, data, &etemp);
    CHECK_STATUS(EG_makeSolidBody);

    /* rotate it */
    cosz = cos(angle * PI / 180);
    sinz = sin(angle * PI / 180);

    matrix[ 0] = cosz;  matrix[ 1] =-sinz;  matrix[ 2] = 0; matrix[ 3] = 0;
    matrix[ 4] = sinz;  matrix[ 5] = cosz;  matrix[ 6] = 0; matrix[ 7] = 0;
    matrix[ 8] = 0;     matrix[ 9] = 0;     matrix[10] = 1; matrix[11] = 0;

    status = EG_makeTransform(TIM->MODL->context, matrix, &exform);
    CHECK_STATUS(EG_makeTransform);

    status = EG_copyObject(etemp, exform, &ebody);
    CHECK_STATUS(EG_copyObject);

    status = EG_deleteObject(etemp);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(exform);
    CHECK_STATUS(EG_deleteObject);

    /* tessellate it */
    params[0] = 0.25;
    params[1] = 0.25;
    params[2] = 10;

    status = EG_makeTessBody(ebody, params, &etess);
    CHECK_STATUS(EG_makeTessBody);

    /* generate the scene graph info needed to visualize it */
    for (iface = 1; iface <= 6; iface++) {
        nitems = 0;

        status = EG_getTessFace(etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        /* vertices */
        status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
        CHECK_STATUS(wv_setData);

        wv_adjustVerts(&(items[nitems]), TIM->sgFocus);
        nitems++;

        /* loop through the triangles and build up the segment table (bias-1) */
        nseg = 0;
        for (itri = 0; itri < ntri; itri++) {
            for (k = 0; k < 3; k++) {
                if (tric[3*itri+k] < itri+1) {
                    nseg++;
                }
            }
        }

        MALLOC(segs, int, 2*nseg);

        nseg = 0;
        for (itri = 0; itri < ntri; itri++) {
            for (k = 0; k < 3; k++) {
                if (tric[3*itri+k] < itri+1) {
                    segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                    segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                    nseg++;
                }
            }
        }

        /* triangles */
        status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
        CHECK_STATUS(wv_setData);
        nitems++;

        color[0] = 1.00;
        color[1] = 0.75;
        color[2] = 0.75;

        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
        CHECK_STATUS(wv_setData);
        nitems++;

        /* triangle backface color */
        color[0] = 0.50;
        color[1] = 0.50;
        color[2] = 0.50;

        status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems]));
        CHECK_STATUS(wv_setData);
        nitems++;

        /* segment indices */
        status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[nitems]));
        CHECK_STATUS(wv_setData);
        nitems++;

        FREE(segs);

        /* segment colors */
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;

        status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems]));
        CHECK_STATUS(wv_setData);
        nitems++;

        /* name */
        snprintf(gpname, 254, "Body %s Face %d", mitn->bodyName, iface);

        /* get index of current graphic primitive (if it exists) */
        SPLINT_CHECK_FOR_NULL(TIM->ctxt);

        igprim = wv_indexGPrim(TIM->ctxt, gpname);

        /* make graphic primitive */
        if (igprim < 0) {
            attrs = WV_ON | WV_ORIENTATION;

            igprim = status = wv_addGPrim(TIM->ctxt, gpname, WV_TRIANGLE, attrs, nitems, items);
            CHECK_STATUS(wv_addGPrim);

            SPLINT_CHECK_FOR_NULL(TIM->ctxt->gPrims);

            TIM->ctxt->gPrims[igprim].lWidth = 1.0;
        } else {
            status = wv_modGPrim(TIM->ctxt, igprim, nitems, items);
            CHECK_STATUS(wv_modGPrim);
        }
    }

    /* delete the tessellation and EGADS objects */
    status = EG_deleteObject(etess);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(ebody);
    CHECK_STATUS(EG_deleteObject);

cleanup:
    FREE(segs);

    return status;
}
