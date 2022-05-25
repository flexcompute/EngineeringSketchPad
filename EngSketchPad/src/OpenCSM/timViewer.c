/*
 ************************************************************************
 *                                                                      *
 * timViewer -- Tool Integration Module for simple viewer               *
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

#include "egads.h"
#include "common.h"
#include "OpenCSM.h"
#include "caps.h"
#include "tim.h"
#include "wsserver.h"
#include "emp.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define TESS_PARAM_0      0.0250
#define TESS_PARAM_1      0.0075
#define TESS_PARAM_2      20.0

static void       addToSgMetaData(int *sgMetaDataUsed, int sgMetaDataSize, char sgMetaData[], char format[], ...);
static int        buildSceneGraphAIM(esp_T *ESP, char aimName[]);
static int        buildSceneGraphBOUND(esp_T *esp, char boundName[], char aimName[], char dataName[]);
static int        buildSceneGraphMODL(esp_T *ESP);
static void       spec_col(float scalar, float out[]);

/***********************************************************************/
/*                                                                     */
/* macros (including those that go along with common.h)                */
/*                                                                     */
/***********************************************************************/

static void *realloc_temp=NULL;            /* used by RALLOC macro */

#define  RED(COLOR)      (float)(COLOR / 0x10000        ) / (float)(255)
#define  GREEN(COLOR)    (float)(COLOR / 0x00100 % 0x100) / (float)(255)
#define  BLUE(COLOR)     (float)(COLOR           % 0x100) / (float)(255)

/***********************************************************************/
/*                                                                     */
/* structures                                                          */
/*                                                                     */
/***********************************************************************/

static float  lims[2]  = {-1.0, +1.0};

/* blue-white-red spectrum */
static float color_map[256*3] =
{ 0.0000, 0.0000, 1.0000,    0.0078, 0.0078, 1.0000,   0.0156, 0.0156, 1.0000,    0.0234, 0.0234, 1.0000,
  0.0312, 0.0312, 1.0000,    0.0391, 0.0391, 1.0000,   0.0469, 0.0469, 1.0000,    0.0547, 0.0547, 1.0000,
  0.0625, 0.0625, 1.0000,    0.0703, 0.0703, 1.0000,   0.0781, 0.0781, 1.0000,    0.0859, 0.0859, 1.0000,
  0.0938, 0.0938, 1.0000,    0.1016, 0.1016, 1.0000,   0.1094, 0.1094, 1.0000,    0.1172, 0.1172, 1.0000,
  0.1250, 0.1250, 1.0000,    0.1328, 0.1328, 1.0000,   0.1406, 0.1406, 1.0000,    0.1484, 0.1484, 1.0000,
  0.1562, 0.1562, 1.0000,    0.1641, 0.1641, 1.0000,   0.1719, 0.1719, 1.0000,    0.1797, 0.1797, 1.0000,
  0.1875, 0.1875, 1.0000,    0.1953, 0.1953, 1.0000,   0.2031, 0.2031, 1.0000,    0.2109, 0.2109, 1.0000,
  0.2188, 0.2188, 1.0000,    0.2266, 0.2266, 1.0000,   0.2344, 0.2344, 1.0000,    0.2422, 0.2422, 1.0000,
  0.2500, 0.2500, 1.0000,    0.2578, 0.2578, 1.0000,   0.2656, 0.2656, 1.0000,    0.2734, 0.2734, 1.0000,
  0.2812, 0.2812, 1.0000,    0.2891, 0.2891, 1.0000,   0.2969, 0.2969, 1.0000,    0.3047, 0.3047, 1.0000,
  0.3125, 0.3125, 1.0000,    0.3203, 0.3203, 1.0000,   0.3281, 0.3281, 1.0000,    0.3359, 0.3359, 1.0000,
  0.3438, 0.3438, 1.0000,    0.3516, 0.3516, 1.0000,   0.3594, 0.3594, 1.0000,    0.3672, 0.3672, 1.0000,
  0.3750, 0.3750, 1.0000,    0.3828, 0.3828, 1.0000,   0.3906, 0.3906, 1.0000,    0.3984, 0.3984, 1.0000,
  0.4062, 0.4062, 1.0000,    0.4141, 0.4141, 1.0000,   0.4219, 0.4219, 1.0000,    0.4297, 0.4297, 1.0000,
  0.4375, 0.4375, 1.0000,    0.4453, 0.4453, 1.0000,   0.4531, 0.4531, 1.0000,    0.4609, 0.4609, 1.0000,
  0.4688, 0.4688, 1.0000,    0.4766, 0.4766, 1.0000,   0.4844, 0.4844, 1.0000,    0.4922, 0.4922, 1.0000,
  0.5000, 0.5000, 1.0000,    0.5078, 0.5078, 1.0000,   0.5156, 0.5156, 1.0000,    0.5234, 0.5234, 1.0000,
  0.5312, 0.5312, 1.0000,    0.5391, 0.5391, 1.0000,   0.5469, 0.5469, 1.0000,    0.5547, 0.5547, 1.0000,
  0.5625, 0.5625, 1.0000,    0.5703, 0.5703, 1.0000,   0.5781, 0.5781, 1.0000,    0.5859, 0.5859, 1.0000,
  0.5938, 0.5938, 1.0000,    0.6016, 0.6016, 1.0000,   0.6094, 0.6094, 1.0000,    0.6172, 0.6172, 1.0000,
  0.6250, 0.6250, 1.0000,    0.6328, 0.6328, 1.0000,   0.6406, 0.6406, 1.0000,    0.6484, 0.6484, 1.0000,
  0.6562, 0.6562, 1.0000,    0.6641, 0.6641, 1.0000,   0.6719, 0.6719, 1.0000,    0.6797, 0.6797, 1.0000,
  0.6875, 0.6875, 1.0000,    0.6953, 0.6953, 1.0000,   0.7031, 0.7031, 1.0000,    0.7109, 0.7109, 1.0000,
  0.7188, 0.7188, 1.0000,    0.7266, 0.7266, 1.0000,   0.7344, 0.7344, 1.0000,    0.7422, 0.7422, 1.0000,
  0.7500, 0.7500, 1.0000,    0.7578, 0.7578, 1.0000,   0.7656, 0.7656, 1.0000,    0.7734, 0.7734, 1.0000,
  0.7812, 0.7812, 1.0000,    0.7891, 0.7891, 1.0000,   0.7969, 0.7969, 1.0000,    0.8047, 0.8047, 1.0000,
  0.8125, 0.8125, 1.0000,    0.8203, 0.8203, 1.0000,   0.8281, 0.8281, 1.0000,    0.8359, 0.8359, 1.0000,
  0.8438, 0.8438, 1.0000,    0.8516, 0.8516, 1.0000,   0.8594, 0.8594, 1.0000,    0.8672, 0.8672, 1.0000,
  0.8750, 0.8750, 1.0000,    0.8828, 0.8828, 1.0000,   0.8906, 0.8906, 1.0000,    0.8984, 0.8984, 1.0000,
  0.9062, 0.9062, 1.0000,    0.9141, 0.9141, 1.0000,   0.9219, 0.9219, 1.0000,    0.9297, 0.9297, 1.0000,
  0.9375, 0.9375, 1.0000,    0.9453, 0.9453, 1.0000,   0.9531, 0.9531, 1.0000,    0.9609, 0.9609, 1.0000,
  0.9688, 0.9688, 1.0000,    0.9766, 0.9766, 1.0000,   0.9844, 0.9844, 1.0000,    0.9922, 0.9922, 1.0000,
  1.0000, 1.0000, 1.0000,    1.0000, 0.9922, 0.9922,   1.0000, 0.9844, 0.9844,    1.0000, 0.9766, 0.9766,
  1.0000, 0.9688, 0.9688,    1.0000, 0.9609, 0.9609,   1.0000, 0.9531, 0.9531,    1.0000, 0.9453, 0.9453,
  1.0000, 0.9375, 0.9375,    1.0000, 0.9297, 0.9297,   1.0000, 0.9219, 0.9219,    1.0000, 0.9141, 0.9141,
  1.0000, 0.9062, 0.9062,    1.0000, 0.8984, 0.8984,   1.0000, 0.8906, 0.8906,    1.0000, 0.8828, 0.8828,
  1.0000, 0.8750, 0.8750,    1.0000, 0.8672, 0.8672,   1.0000, 0.8594, 0.8594,    1.0000, 0.8516, 0.8516,
  1.0000, 0.8438, 0.8438,    1.0000, 0.8359, 0.8359,   1.0000, 0.8281, 0.8281,    1.0000, 0.8203, 0.8203,
  1.0000, 0.8125, 0.8125,    1.0000, 0.8047, 0.8047,   1.0000, 0.7969, 0.7969,    1.0000, 0.7891, 0.7891,
  1.0000, 0.7812, 0.7812,    1.0000, 0.7734, 0.7734,   1.0000, 0.7656, 0.7656,    1.0000, 0.7578, 0.7578,
  1.0000, 0.7500, 0.7500,    1.0000, 0.7422, 0.7422,   1.0000, 0.7344, 0.7344,    1.0000, 0.7266, 0.7266,
  1.0000, 0.7188, 0.7188,    1.0000, 0.7109, 0.7109,   1.0000, 0.7031, 0.7031,    1.0000, 0.6953, 0.6953,
  1.0000, 0.6875, 0.6875,    1.0000, 0.6797, 0.6797,   1.0000, 0.6719, 0.6719,    1.0000, 0.6641, 0.6641,
  1.0000, 0.6562, 0.6562,    1.0000, 0.6484, 0.6484,   1.0000, 0.6406, 0.6406,    1.0000, 0.6328, 0.6328,
  1.0000, 0.6250, 0.6250,    1.0000, 0.6172, 0.6172,   1.0000, 0.6094, 0.6094,    1.0000, 0.6016, 0.6016,
  1.0000, 0.5938, 0.5938,    1.0000, 0.5859, 0.5859,   1.0000, 0.5781, 0.5781,    1.0000, 0.5703, 0.5703,
  1.0000, 0.5625, 0.5625,    1.0000, 0.5547, 0.5547,   1.0000, 0.5469, 0.5469,    1.0000, 0.5391, 0.5391,
  1.0000, 0.5312, 0.5312,    1.0000, 0.5234, 0.5234,   1.0000, 0.5156, 0.5156,    1.0000, 0.5078, 0.5078,
  1.0000, 0.5000, 0.5000,    1.0000, 0.4922, 0.4922,   1.0000, 0.4844, 0.4844,    1.0000, 0.4766, 0.4766,
  1.0000, 0.4688, 0.4688,    1.0000, 0.4609, 0.4609,   1.0000, 0.4531, 0.4531,    1.0000, 0.4453, 0.4453,
  1.0000, 0.4375, 0.4375,    1.0000, 0.4297, 0.4297,   1.0000, 0.4219, 0.4219,    1.0000, 0.4141, 0.4141,
  1.0000, 0.4062, 0.4062,    1.0000, 0.3984, 0.3984,   1.0000, 0.3906, 0.3906,    1.0000, 0.3828, 0.3828,
  1.0000, 0.3750, 0.3750,    1.0000, 0.3672, 0.3672,   1.0000, 0.3594, 0.3594,    1.0000, 0.3516, 0.3516,
  1.0000, 0.3438, 0.3438,    1.0000, 0.3359, 0.3359,   1.0000, 0.3281, 0.3281,    1.0000, 0.3203, 0.3203,
  1.0000, 0.3125, 0.3125,    1.0000, 0.3047, 0.3047,   1.0000, 0.2969, 0.2969,    1.0000, 0.2891, 0.2891,
  1.0000, 0.2812, 0.2812,    1.0000, 0.2734, 0.2734,   1.0000, 0.2656, 0.2656,    1.0000, 0.2578, 0.2578,
  1.0000, 0.2500, 0.2500,    1.0000, 0.2422, 0.2422,   1.0000, 0.2344, 0.2344,    1.0000, 0.2266, 0.2266,
  1.0000, 0.2188, 0.2188,    1.0000, 0.2109, 0.2109,   1.0000, 0.2031, 0.2031,    1.0000, 0.1953, 0.1953,
  1.0000, 0.1875, 0.1875,    1.0000, 0.1797, 0.1797,   1.0000, 0.1719, 0.1719,    1.0000, 0.1641, 0.1641,
  1.0000, 0.1562, 0.1562,    1.0000, 0.1484, 0.1484,   1.0000, 0.1406, 0.1406,    1.0000, 0.1328, 0.1328,
  1.0000, 0.1250, 0.1250,    1.0000, 0.1172, 0.1172,   1.0000, 0.1094, 0.1094,    1.0000, 0.1016, 0.1016,
  1.0000, 0.0938, 0.0938,    1.0000, 0.0859, 0.0859,   1.0000, 0.0781, 0.0781,    1.0000, 0.0703, 0.0703,
  1.0000, 0.0625, 0.0625,    1.0000, 0.0547, 0.0547,   1.0000, 0.0469, 0.0469,    1.0000, 0.0391, 0.0391,
  1.0000, 0.0312, 0.0312,    1.0000, 0.0234, 0.0234,   1.0000, 0.0156, 0.0156,    1.0000, 0.0078, 0.0078 };

static int outLevel = 1;


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(/*@unused@*/esp_T *ESP,                     /* (in)  pointer to ESP structure */
        /*@unused@*/void  *data)                    /* (in)  component name */
{
    int    status;                      /* (out) return status */

    ROUTINE(timLoad(mitten));

    /* --------------------------------------------------------------- */

    /* hold the UI when executing */
    status = 1;

//cleanup:
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

    int    ibody, iface;
    char   *arg1=NULL, *arg2=NULL, *arg3=NULL;
    modl_T *MODL;

    ROUTINE(timMesg(mitten));

    /* --------------------------------------------------------------- */

    /* "MODL|" */
    if        (strncmp(command, "MODL", 4) == 0) {

        buildSceneGraphMODL(ESP);

    /* "AIM|aimName|" */
    } else if (strncmp(command, "AIM|", 4) == 0) {

        /* extract argument */
        GetToken(command, 1, '|', &arg1);

        SPLINT_CHECK_FOR_NULL(arg1);

        buildSceneGraphAIM(ESP, arg1);

        FREE(arg1);

    /* "BOUND|boundName|aimName|dataName|" */
    } else if (strncmp(command, "BOUND|", 6) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        GetToken(command, 2, '|', &arg2);
        GetToken(command, 3, '|', &arg3);

        SPLINT_CHECK_FOR_NULL(arg1);
        SPLINT_CHECK_FOR_NULL(arg2);
        SPLINT_CHECK_FOR_NULL(arg3);

        buildSceneGraphBOUND(ESP, arg1, arg2, arg3);

        FREE(arg1);
        FREE(arg2);
        FREE(arg3);

    /* "red|" */
    } else if (strncmp(command, "red|", 4) == 0) {

        if (ESP == NULL) {
            SPRINT0(0, "WARNING:: not running via serveESP");
            goto cleanup;
        } else {
            MODL = ESP->MODL;

            for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                if (MODL->body[ibody].onstack != 1) continue;

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    MODL->body[ibody].face[iface].gratt.color = 0x00FF0000;
                }
            }
        }

        buildSceneGraphMODL(ESP);

    /* "green|" */
    } else if (strncmp(command, "green|", 6) == 0) {

        if (ESP == NULL) {
            SPRINT0(0, "WARNING:: not running via serveESP");
            goto cleanup;
        } else {
            MODL = ESP->MODL;

            for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                if (MODL->body[ibody].onstack != 1) continue;

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    MODL->body[ibody].face[iface].gratt.color = 0x0000FF00;
                }
            }
        }

        buildSceneGraphMODL(ESP);

    /* "blue|" */
    } else if (strncmp(command, "blue|", 5) == 0) {

        if (ESP == NULL) {
            SPRINT0(0, "WARNING:: not running via serveESP");
            goto cleanup;
        } else {
            MODL = ESP->MODL;

            for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                if (MODL->body[ibody].onstack != 1) continue;

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    MODL->body[ibody].face[iface].gratt.color = 0x000000FF;
                }
            }
        }

        buildSceneGraphMODL(ESP);
    }

cleanup:
    FREE(arg1);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(/*@unused@*/esp_T *ESP)                     /* (in)  pointer to TI Mstructure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timSave(mitten));

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(/*@unused@*/esp_T *ESP,                     /* (in)  pointer to ESP structure */
        /*@unused@*/int   unload)                   /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timQuit(mitten));

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   buildSceneGraphMODL - make a scene graph for the MODL             */
/*                                                                     */
/***********************************************************************/

static int
buildSceneGraphMODL(esp_T  *ESP)
{
    int       status = SUCCESS;         /* return status */

    int       ibody, jbody, iface, iedge, inode, iattr, nattr, atype, alen;
    int       npnt, ipnt, ntri, itri, igprim, nseg, i, k;
    int       attrs, head[3], nitems, nnode, nedge, nface;
    int       oclass, mtype, nchild, *senses;
    int       *segs=NULL, *ivrts=NULL;
    int        npatch2, ipatch, n1, n2, i1, i2, *Tris=NULL;
    CINT      *ptype, *pindx, *tris, *tric, *pvindex, *pbounds;
    float     color[18], *pcolors=NULL, *plotdata=NULL, *segments=NULL, *tuft=NULL;
    double    bigbox[6], box[6], size, xyz_dum[6], axis[18];
    CDOUBLE   *xyz, *uv, *t;
    char      gpname[MAX_STRVAL_LEN], bname[MAX_NAME_LEN];
    ego       ebody, etess, eface, eedge, enode, eref, *echilds;
    ego       *enodes=NULL, *eedges=NULL, *efaces=NULL;

    int       attrLen, attrType;
    CINT      *tempIlist, *nquad=NULL;
    CDOUBLE   *tempRlist;
    CCHAR     *attrName, *tempClist;

    modl_T    *MODL;
    wvContext *cntxt;
    wvData    items[6];

/* global variables associated with scene graph meta-data */
#define MAX_METADATA_CHUNK 32000
    int       sgMetaDataSize = 0;
    int       sgMetaDataUsed = 0;
    char      *sgMetaData    = NULL;
    char      *sgFocusData   = NULL;

    ROUTINE(buildSceneGraphMODL);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        SPRINT0(0, "WARNING:: not running via serveESP");
        goto cleanup;
    } else {
        MODL  = ESP->MODL;
        cntxt = ESP->cntxt;
    }

    sgMetaDataSize = MAX_METADATA_CHUNK;

    MALLOC(sgMetaData,  char, sgMetaDataSize);
    MALLOC(sgFocusData, char, MAX_STRVAL_LEN);

    /* set the mutex.  if the mutex wss already set, wait until released */
    EMP_LockSet(ESP->sgMutex);

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    if (MODL == NULL) {
        goto cleanup;
    }

    /* close the key (from a previous view?) */
    status = wv_setKey(cntxt, 0, NULL, lims[0], lims[1], NULL);
    if (status != SUCCESS) {
        SPRINT1(9, "ERROR:: wv_setKey -> status=%d", status);
    }
    wv_broadcastText("setWvKey|off|");

    /* find the values needed to adjust the vertices */
    bigbox[0] = bigbox[1] = bigbox[2] = +HUGEQ;
    bigbox[3] = bigbox[4] = bigbox[5] = -HUGEQ;

    /* use Bodys on stack */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        ebody = MODL->body[ibody].ebody;

        status = EG_getBoundingBox(ebody, box);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_getBoundingBox(%d) -> status=%d", ibody, status);
        }

        if (box[0] < bigbox[0]) bigbox[0] = box[0];
        if (box[1] < bigbox[1]) bigbox[1] = box[1];
        if (box[2] < bigbox[2]) bigbox[2] = box[2];
        if (box[3] > bigbox[3]) bigbox[3] = box[3];
        if (box[4] > bigbox[4]) bigbox[4] = box[4];
        if (box[5] > bigbox[5]) bigbox[5] = box[5];
    }

                                    size = bigbox[3] - bigbox[0];
    if (size < bigbox[4]-bigbox[1]) size = bigbox[4] - bigbox[1];
    if (size < bigbox[5]-bigbox[2]) size = bigbox[5] - bigbox[2];

    ESP->sgFocus[0] = (bigbox[0] + bigbox[3]) / 2;
    ESP->sgFocus[1] = (bigbox[1] + bigbox[4]) / 2;
    ESP->sgFocus[2] = (bigbox[2] + bigbox[5]) / 2;
    ESP->sgFocus[3] = size;

    /* generate the scene graph focus data */
    snprintf(sgFocusData, MAX_STRVAL_LEN-1, "sgFocus|[%20.12e,%20.12e,%20.12e,%20.12e]",
             ESP->sgFocus[0], ESP->sgFocus[1], ESP->sgFocus[2], ESP->sgFocus[3]);

    /* initialize the scene graph meta data */
    sgMetaData[ 0] = '\0';
    sgMetaDataUsed = 0;

    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "sgData|{");

    /* loop through the Bodys */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        if (MODL->body[ibody].eebody == NULL) {
            ebody = MODL->body[ibody].ebody;

            EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
            EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
            EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        } else {
            ebody = MODL->body[ibody].eebody;

            EG_getBodyTopos(ebody, NULL, NODE,  &nnode, &enodes);
            EG_getBodyTopos(ebody, NULL, EEDGE, &nedge, &eedges);
            EG_getBodyTopos(ebody, NULL, EFACE, &nface, &efaces);
        }

        /* set up Body name */
        snprintf(bname, MAX_NAME_LEN-1, "Body %d", ibody);

        status = EG_attributeRet(ebody, "_name", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS && attrType == ATTRSTRING) {
            snprintf(bname, MAX_NAME_LEN-1, "%s", tempClist);
        }

        /* check for duplicate Body names */
        for (jbody = 1; jbody < ibody; jbody++) {
            if (MODL->body[jbody].onstack != 1) continue;

            status = EG_attributeRet(MODL->body[jbody].ebody, "_name", &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status == SUCCESS && attrType == ATTRSTRING) {
                if (strcmp(tempClist, bname) == 0) {
                    SPRINT2(0, "WARNING:: duplicate Body name (%s) found; being changed to \"Body %d\"", bname, ibody);
                    snprintf(bname, MAX_NAME_LEN-1, "Body %d", ibody);
                }
            }
        }

        /* Body info for NodeBody, SheetBody or SolidBody */
        snprintf(gpname, MAX_STRVAL_LEN, "%s", bname);
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        /* add Node to meta data (if there is room) */
        if (nattr > 0) {
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[", gpname);
        } else {
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[\"body\",\"%d\"", gpname, ibody);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(ebody, iattr, &attrName, &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iattr, status);
            }

            if (attrType == ATTRCSYS) continue;

            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

            if        (attrType == ATTRINT) {
                for (i = 0; i < attrLen ; i++) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                }
            } else if (attrType == ATTRREAL) {
                for (i = 0; i < attrLen ; i++) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                }
            } else if (attrType == ATTRSTRING) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
            }

            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
        }
        sgMetaDataUsed--;
        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");

        if (MODL->body[ibody].eetess == NULL) {
            etess = MODL->body[ibody].etess;
        } else {
            etess = MODL->body[ibody].eetess;
        }

        /* loop through the Faces within the current Body */
        for (iface = 1; iface <= nface; iface++) {
            nitems = 0;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Face %d", bname, iface);
            attrs = WV_ON | WV_ORIENTATION;

            status = EG_getQuads(etess, iface,
                                 &npnt, &xyz, &uv, &ptype, &pindx, &npatch2);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getQuads(%d,%d) -> status=%d", ibody, iface, status);
            }

            status = EG_attributeRet(etess, ".tessType",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);

            /* new-style Quads */
            if (status == SUCCESS && atype == ATTRSTRING && (
                    strcmp(tempClist, "Quad") == 0 || strcmp(tempClist, "Mixed") == 0)) {
                status = EG_attributeRet(etess, ".mixed",
                                         &atype, &alen, &nquad, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeRet(%d,%d) -> status=%d", ibody, iface, status);
                }

                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
                nitems++;

                /* loop through the triangles and build up the segment table */
                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            nseg++;
                        }
                    }
                }

                MALLOC(segs, int, 2*nseg);

                /* create segments between Triangles (bias-1) */
                SPLINT_CHECK_FOR_NULL(nquad);

                nseg = 0;
                for (itri = 0; itri < ntri-2*nquad[iface-1]; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                            segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                            nseg++;
                        }
                    }
                }

                /* create segments between Triangles (but not within the pair) (bias-1) */
                for (; itri < ntri; itri++) {
                    if (tric[3*itri  ] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                    itri++;

                    if (tric[3*itri  ] < itri) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                }

                /* triangles */
                status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

            /* patches in old-style Quads */
            } else if (npatch2 > 0) {
                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
                nitems++;

                /* loop through the patches and build up the triangle and
                   segment tables (bias-1) */
                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    ntri += 2 * (n1-1) * (n2-1);
                    nseg += n1 * (n2-1) + n2 * (n1-1);
                }

                MALLOC(Tris, int, 3*ntri);
                MALLOC(segs, int, 2*nseg);

                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    for (i2 = 1; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            Tris[3*ntri  ] = pvindex[(i1-1)+n1*(i2-1)];
                            Tris[3*ntri+1] = pvindex[(i1  )+n1*(i2-1)];
                            Tris[3*ntri+2] = pvindex[(i1  )+n1*(i2  )];
                            ntri++;

                            Tris[3*ntri  ] = pvindex[(i1  )+n1*(i2  )];
                            Tris[3*ntri+1] = pvindex[(i1-1)+n1*(i2  )];
                            Tris[3*ntri+2] = pvindex[(i1-1)+n1*(i2-1)];
                            ntri++;
                        }
                    }

                    for (i2 = 0; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            segs[2*nseg  ] = pvindex[(i1-1)+n1*(i2)];
                            segs[2*nseg+1] = pvindex[(i1  )+n1*(i2)];
                            nseg++;
                        }
                    }

                    for (i1 = 0; i1 < n1; i1++) {
                        for (i2 = 1; i2 < n2; i2++) {
                            segs[2*nseg  ] = pvindex[(i1)+n1*(i2-1)];
                            segs[2*nseg+1] = pvindex[(i1)+n1*(i2  )];
                            nseg++;
                        }
                    }
                }

                /* two triangles for rendering each quadrilateral */
                status = wv_setData(WV_INT32, 3*ntri, (void*)Tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(Tris);

            /* Triangles */
            } else {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
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
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;
            }

            /* constant triangle colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].face[iface].gratt.color);
                color[1] = GREEN(MODL->body[ibody].face[iface].gratt.color);
                color[2] = BLUE( MODL->body[ibody].face[iface].gratt.color);
            } else {
                color[0] = 0.75;
                color[1] = 0.75;
                color[2] = 1.00;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* triangle backface color */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].face[iface].gratt.bcolor);
                color[1] = GREEN(MODL->body[ibody].face[iface].gratt.bcolor);
                color[2] = BLUE( MODL->body[ibody].face[iface].gratt.bcolor);
            } else {
                color[0] = 0.50;
                color[1] = 0.50;
                color[2] = 0.50;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* segment indices */
            status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            FREE(segs);

            /* segment colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].face[iface].gratt.mcolor);
                color[1] = GREEN(MODL->body[ibody].face[iface].gratt.mcolor);
                color[2] = BLUE( MODL->body[ibody].face[iface].gratt.mcolor);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* get Attributes for the Face */
            SPLINT_CHECK_FOR_NULL(efaces);
            eface  = efaces[iface-1];
            status = EG_attributeNum(eface, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iface, status);
            }

            /* add Face to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eface, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iface, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Edges within the current Body */
        for (iedge = 1; iedge <= nedge; iedge++) {
            nitems = 0;

            /* skip edge if this is a NodeBody */
            if (MODL->body[ibody].botype == OCSM_NODE_BODY) continue;

            status = EG_getTessEdge(etess, iedge, &npnt, &xyz, &t);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTessEdge(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Edge %d", bname, iedge);
            attrs = WV_ON;

            /* vertices */
            status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* segments (bias-1) */
            MALLOC(ivrts, int, 2*(npnt-1));

            for (ipnt = 0; ipnt < npnt-1; ipnt++) {
                ivrts[2*ipnt  ] = ipnt + 1;
                ivrts[2*ipnt+1] = ipnt + 2;
            }

            status = wv_setData(WV_INT32, 2*(npnt-1), (void*)ivrts, WV_INDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* line colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].edge[iedge].gratt.color);
                color[1] = GREEN(MODL->body[ibody].edge[iedge].gratt.color);
                color[2] = BLUE( MODL->body[ibody].edge[iedge].gratt.color);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* points */
            MALLOC(ivrts, int, npnt);

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ivrts[ipnt] = ipnt + 1;
            }

            status = wv_setData(WV_INT32, npnt, (void*)ivrts, WV_PINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* point colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].edge[iedge].gratt.mcolor);
                color[1] = GREEN(MODL->body[ibody].edge[iedge].gratt.mcolor);
                color[2] = BLUE( MODL->body[ibody].edge[iedge].gratt.mcolor);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2 ]= 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 2 (does not work for ANGLE) */
                cntxt->gPrims[igprim].lWidth = 2.0;

                /* make point size 5 */
                cntxt->gPrims[igprim].pSize  = 5.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = npnt - 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/ESP->sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_addArrowHeads(%d,%d) -> status=%d", ibody, iedge, status);
                }
            }

            SPLINT_CHECK_FOR_NULL(eedges);
            eedge  = eedges[iedge-1];
            status = EG_attributeNum(eedge, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* add Edge to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eedge, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Nodes within the current Body */
        for (inode = 1; inode <= nnode; inode++) {
            nitems = 0;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Node %d", bname, inode);

            /* default for NodeBodys is to turn the Node on */
            if (MODL->body[ibody].botype == OCSM_NODE_BODY) {
                attrs = WV_ON;
            } else {
                attrs = 0;
            }

            SPLINT_CHECK_FOR_NULL(enodes);
            enode  = enodes[inode-1];

            /* two copies of vertex (to get around possible wv limitation) */
            status = EG_getTopology(enode, &eref, &oclass, &mtype,
                                    xyz_dum, &nchild, &echilds, &senses);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTopology(%d,%d) -> status=%d", ibody, inode, status);
            }

            xyz_dum[3] = xyz_dum[0];
            xyz_dum[4] = xyz_dum[1];
            xyz_dum[5] = xyz_dum[2];

            status = wv_setData(WV_REAL64, 2, (void*)xyz_dum, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* node color (black) */
            color[0] = 0;   color[1] = 0;   color[2] = 0;

            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].node[inode].gratt.color);
                color[1] = GREEN(MODL->body[ibody].node[inode].gratt.color);
                color[2] = BLUE( MODL->body[ibody].node[inode].gratt.color);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].pSize = 6.0;
            }

            status = EG_attributeNum(enode, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* add Node to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(enode, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Csystems associated with the current Body */
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            nitems = 0;

            status = EG_attributeGet(ebody, iattr, &attrName, &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_attributeGet -> status=%d", status);
            }

            if (attrType != ATTRCSYS) continue;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Csys %s", bname, attrName);
            attrs = WV_ON | WV_SHADING | WV_ORIENTATION;

            /* vertices */
            axis[ 0] = tempRlist[attrLen  ];
            axis[ 1] = tempRlist[attrLen+1];
            axis[ 2] = tempRlist[attrLen+2];
            axis[ 3] = tempRlist[attrLen  ] + tempRlist[attrLen+ 3];
            axis[ 4] = tempRlist[attrLen+1] + tempRlist[attrLen+ 4];
            axis[ 5] = tempRlist[attrLen+2] + tempRlist[attrLen+ 5];
            axis[ 6] = tempRlist[attrLen  ];
            axis[ 7] = tempRlist[attrLen+1];
            axis[ 8] = tempRlist[attrLen+2];
            axis[ 9] = tempRlist[attrLen  ] + tempRlist[attrLen+ 6];
            axis[10] = tempRlist[attrLen+1] + tempRlist[attrLen+ 7];
            axis[11] = tempRlist[attrLen+2] + tempRlist[attrLen+ 8];
            axis[12] = tempRlist[attrLen  ];
            axis[13] = tempRlist[attrLen+1];
            axis[14] = tempRlist[attrLen+2];
            axis[15] = tempRlist[attrLen  ] + tempRlist[attrLen+ 9];
            axis[16] = tempRlist[attrLen+1] + tempRlist[attrLen+10];
            axis[17] = tempRlist[attrLen+2] + tempRlist[attrLen+11];

            status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* line colors (multiple colors requires that WV_SHADING be set above) */
            color[ 0] = 1;   color[ 1] = 0;   color[ 2] = 0;
            color[ 3] = 1;   color[ 4] = 0;   color[ 5] = 0;
            color[ 6] = 0;   color[ 7] = 1;   color[ 8] = 0;
            color[ 9] = 0;   color[10] = 1;   color[11] = 0;
            color[12] = 0;   color[13] = 0;   color[14] = 1;
            color[15] = 0;   color[16] = 0;   color[17] = 1;
            status = wv_setData(WV_REAL32, 6, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 1 */
                cntxt->gPrims[igprim].lWidth = 1.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/ESP->sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_addArrowHeads -> status=%d", status);
                }
            }

            /* add Csystem to meta data (if there is room) */
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[],", gpname);
        }
    }

    /* axes */
    nitems = 0;

    /* name and attributes */
    snprintf(gpname, MAX_STRVAL_LEN-1, "Axes");
    attrs = 0;

    /* vertices */
    axis[ 0] = MIN(2*bigbox[0]-bigbox[3], 0);   axis[ 1] = 0;   axis[ 2] = 0;
    axis[ 3] = MAX(2*bigbox[3]-bigbox[0], 0);   axis[ 4] = 0;   axis[ 5] = 0;

    axis[ 6] = 0;   axis[ 7] = MIN(2*bigbox[1]-bigbox[4], 0);   axis[ 8] = 0;
    axis[ 9] = 0;   axis[10] = MAX(2*bigbox[4]-bigbox[1], 0);   axis[11] = 0;

    axis[12] = 0;   axis[13] = 0;   axis[14] = MIN(2*bigbox[2]-bigbox[5], 0);
    axis[15] = 0;   axis[16] = 0;   axis[17] = MAX(2*bigbox[5]-bigbox[2], 0);
    status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
    }

    wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
    nitems++;

    /* line color */
    color[0] = 0.7;   color[1] = 0.7;   color[2] = 0.7;
    status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
    }
    nitems++;

    /* make graphic primitive */
    igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
    if (igprim < 0) {
        SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
    } else {
        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].lWidth = 1.0;
    }

    /* finish the scene graph meta data */
    sgMetaDataUsed--;
    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "}");

    /* send the scene graph meta data if it has not already been sent */
    if (STRLEN(sgMetaData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgData: %s", sgMetaData);
        wv_broadcastText(sgMetaData);

        /* nullify meta data so that it does not get sent again */
        sgMetaData[0]  = '\0';
        sgMetaDataUsed = 0;
    }

    if (STRLEN(sgFocusData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgFocus: %s", sgFocusData);
        wv_broadcastText(sgFocusData);
    }

cleanup:

    /* release the mutex so that another thread can access the scene graph */
    if (ESP != NULL) {
        EMP_LockRelease(ESP->sgMutex);
    }

    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    FREE(plotdata);
    FREE(pcolors );
    FREE(segments);
    FREE(segs    );
    FREE(ivrts   );
    FREE(tuft    );
    FREE(Tris    );

    FREE(sgFocusData);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   buildSceneGraphAIM - make a scene graph for wv                    */
/*                                                                     */
/***********************************************************************/

static int
buildSceneGraphAIM(esp_T  *ESP,
                   char   aimName[])
{
    int       status = SUCCESS;         /* return status */

    int       ibody, jbody, itess, iface, iedge, inode, iattr, nattr, atype, alen;
    int       npnt, ipnt, ntri, itri, igprim, nseg, i, k, istat;
    int       attrs, head[3], nitems, nnode, nedge, nface;
    int       oclass, mtype, nchild, *senses, *segs=NULL, *ivrts=NULL;
    int       npatch2, ipatch, n1, n2, i1, i2, *Tris=NULL;
    CINT      *ptype, *pindx, *tris, *tric, *pvindex, *pbounds;
    float     color[18];
    double    bigbox[6], bbox[6], size, xyz_dum[6], axis[18], params[3];
    CDOUBLE   *xyz, *uv, *t;
    char      gpname[MAX_STRVAL_LEN], bname[MAX_NAME_LEN];
    ego       ebody, etess, eface, eedge, enode, eref, *echilds, etemp;
    ego       *enodes=NULL, *eedges=NULL, *efaces=NULL, *etesss=NULL;

    int       attrLen, attrType, nbody, ntess, nerror;
    CINT      *tempIlist, *nquad=NULL;
    CDOUBLE   *tempRlist;
    char      *units;
    CCHAR     *attrName, *tempClist;

    modl_T    *MODL;
    capsObj   aimObj;
    capsErrs  *errors;
    wvContext *cntxt;
    wvData    items[6];

/* global variables associated with scene graph meta-data */
#define MAX_METADATA_CHUNK 32000
    int       sgMetaDataSize = 0;
    int       sgMetaDataUsed = 0;
    char      *sgMetaData    = NULL;
    char      *sgFocusData   = NULL;

    ROUTINE(buildSceneGraphAIM);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        SPRINT0(0, "WARNING:: not running via serveESP");
        goto cleanup;
    } else {
        MODL  = ESP->MODL;
        cntxt = ESP->cntxt;
    }

    /* find the AIM with the given name */
    status = caps_childByName(ESP->CAPS, ANALYSIS, NONE, aimName, &aimObj, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS) {
        SPRINT2(0, "ERROR:: caps_childByName -> status=%d, nerror=%d", status, nerror);
    }

    sgMetaDataSize = MAX_METADATA_CHUNK;

    MALLOC(sgMetaData,  char, sgMetaDataSize);
    MALLOC(sgFocusData, char, MAX_STRVAL_LEN);

    /* set the mutex.  if the mutex wss already set, wait until released */
    EMP_LockSet(ESP->sgMutex);

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    /* close the key (from a previous view?) */
    status = wv_setKey(cntxt, 0, NULL, lims[0], lims[1], NULL);
    if (status != SUCCESS) {
        SPRINT1(9, "ERROR:: wv_setKey -> status=%d", status);
    }
    wv_broadcastText("setWvKey|off|");

    /* get the number of Bodys associated with this AIM */
    status = caps_size(aimObj, BODIES, NONE, &nbody, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS) {
        SPRINT3(0, "ERROR:: caps_size -> status=%d, nbody=%d, nerror=%d", status, nbody, nerror);
    }

    /* get the list of Tessellations associated with this AIM */
    status = caps_getTessels(aimObj, &ntess, &etesss, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS) {
        SPRINT3(0, "ERROR:: caps_getTessels -> status=%d, ntess=%d, nerror=%d", status, ntess, nerror);
    }

    /* find the values needed to adjust the vertices */
    bigbox[0] = bigbox[1] = bigbox[2] = +HUGEQ;
    bigbox[3] = bigbox[4] = bigbox[5] = -HUGEQ;

    /* use Bodys on stack */
    ocsmSetOutLevel(2);
    for (ibody = 1; ibody <= nbody; ibody++) {
        status = caps_bodyByIndex(aimObj, ibody, &ebody, &units);
        if (status != SUCCESS) {
            SPRINT3(0, "ERROR:: caps_bodyByIndex(%d) -> status=%d, units=%s", ibody, status, units);
        }

        status = EG_getBoundingBox(ebody, bbox);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_getBoundingBox(%d) -> status=%d", ibody, status);
        }

        if (bbox[0] < bigbox[0]) bigbox[0] = bbox[0];
        if (bbox[1] < bigbox[1]) bigbox[1] = bbox[1];
        if (bbox[2] < bigbox[2]) bigbox[2] = bbox[2];
        if (bbox[3] > bigbox[3]) bigbox[3] = bbox[3];
        if (bbox[4] > bigbox[4]) bigbox[4] = bbox[4];
        if (bbox[5] > bigbox[5]) bigbox[5] = bbox[5];
    }
    ocsmSetOutLevel(1);

                                    size = bigbox[3] - bigbox[0];
    if (size < bigbox[4]-bigbox[1]) size = bigbox[4] - bigbox[1];
    if (size < bigbox[5]-bigbox[2]) size = bigbox[5] - bigbox[2];

    ESP->sgFocus[0] = (bigbox[0] + bigbox[3]) / 2;
    ESP->sgFocus[1] = (bigbox[1] + bigbox[4]) / 2;
    ESP->sgFocus[2] = (bigbox[2] + bigbox[5]) / 2;
    ESP->sgFocus[3] = size;

    /* generate the scene graph focus data */
    snprintf(sgFocusData, MAX_STRVAL_LEN-1, "sgFocus|[%20.12e,%20.12e,%20.12e,%20.12e]",
             ESP->sgFocus[0], ESP->sgFocus[1], ESP->sgFocus[2], ESP->sgFocus[3]);

    /* initialize the scene graph meta data */
    sgMetaData[ 0] = '\0';
    sgMetaDataUsed = 0;

    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "sgData|{");

    /* loop through the Bodys */
    for (ibody = 1; ibody <= nbody; ibody++) {
        status = caps_bodyByIndex(aimObj, ibody, &ebody, &units);
        if (status != SUCCESS) {
            SPRINT3(0, "ERROR:: caps_bodyByIndex(%d) -> status=%d, units=%s", ibody, status, units);
        }

        EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
        EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);

        /* set up Body name */
        snprintf(bname, MAX_NAME_LEN-1, "Body %d", ibody);

        status = EG_attributeRet(ebody, "_name", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS && attrType == ATTRSTRING) {
            snprintf(bname, MAX_NAME_LEN-1, "%s", tempClist);
        }

        /* Body info for NodeBody, SheetBody or SolidBody */
        snprintf(gpname, MAX_STRVAL_LEN, "%s", bname);
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        /* add Node to meta data (if there is room) */
        if (nattr > 0) {
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[", gpname);
        } else {
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[\"body\",\"%d\"", gpname, ibody);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(ebody, iattr, &attrName, &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iattr, status);
            }

            if (attrType == ATTRCSYS) continue;

            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

            if        (attrType == ATTRINT) {
                for (i = 0; i < attrLen ; i++) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                }
            } else if (attrType == ATTRREAL) {
                for (i = 0; i < attrLen ; i++) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                }
            } else if (attrType == ATTRSTRING) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
            }

            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
        }
        sgMetaDataUsed--;
        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");

        etess = NULL;

        /* see if there is a Tessellation associated with ebody */
        for (itess = 0; itess < ntess; itess++) {
            SPLINT_CHECK_FOR_NULL(etesss);

            status = EG_statusTessBody(etesss[itess], &etemp, &istat, &npnt);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_statusTessBody -> status=%d, istat=%d, npnt=%d", status, istat, npnt);
            }

            if (status == SUCCESS && istat == 1 && etemp == ebody) {
                etess = etesss[itess];
                break;
            }
        }

        /* if etess is still NULL, find the Body in the MODL that matches
           ebody and use its Tessellation */
        if (etess == NULL) {
            for (jbody = 1; jbody <= MODL->nbody; jbody++) {
                if (MODL->body[jbody].ebody == ebody) {
                    etess = MODL->body[jbody].etess;
                    break;
                }
            }
        }

        /* if etess is still NULL, make a tessellation */
        if (etess == NULL) {
            status = EG_getBoundingBox(ebody, bbox);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_getBoundingBox -> status=%d", status);
            }

            size = sqrt(SQR(bbox[3]-bbox[0]) + SQR(bbox[4]-bbox[1]) + SQR(bbox[5]-bbox[2]));

            params[0] = TESS_PARAM_0 * size;
            params[1] = TESS_PARAM_1 * size;
            params[2] = TESS_PARAM_2;

            status = EG_makeTessBody(ebody, params, &etess);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_makeTessBody -> status=%d", status);
            }

            for (jbody = 1; jbody <= MODL->nbody; jbody++) {
                if (MODL->body[jbody].ebody == ebody) {
                    MODL->body[jbody].etess = etess;
                    break;
                }
            }
        }

        SPLINT_CHECK_FOR_NULL(etess);

        /* loop through the Faces within the current Body */
        for (iface = 1; iface <= nface; iface++) {
            nitems = 0;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Face %d", bname, iface);
            attrs = WV_ON | WV_ORIENTATION;

            status = EG_getQuads(etess, iface,
                                 &npnt, &xyz, &uv, &ptype, &pindx, &npatch2);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getQuads(%d,%d) -> status=%d", ibody, iface, status);
            }

            status = EG_attributeRet(etess, ".tessType",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);

            /* new-style Quads */
            if (status == SUCCESS && atype == ATTRSTRING && (
                    strcmp(tempClist, "Quad") == 0 || strcmp(tempClist, "Mixed") == 0)) {
                status = EG_attributeRet(etess, ".mixed",
                                         &atype, &alen, &nquad, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeRet(%d,%d) -> status=%d", ibody, iface, status);
                }

                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
                nitems++;

                /* loop through the triangles and build up the segment table */
                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            nseg++;
                        }
                    }
                }

                MALLOC(segs, int, 2*nseg);

                /* create segments between Triangles (bias-1) */
                SPLINT_CHECK_FOR_NULL(nquad);

                nseg = 0;
                for (itri = 0; itri < ntri-2*nquad[iface-1]; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                            segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                            nseg++;
                        }
                    }
                }

                /* create segments between Triangles (but not within the pair) (bias-1) */
                for (; itri < ntri; itri++) {
                    if (tric[3*itri  ] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                    itri++;

                    if (tric[3*itri  ] < itri) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                }

                /* triangles */
                status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

            /* patches in old-style Quads */
            } else if (npatch2 > 0) {
                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
                nitems++;

                /* loop through the patches and build up the triangle and
                   segment tables (bias-1) */
                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    ntri += 2 * (n1-1) * (n2-1);
                    nseg += n1 * (n2-1) + n2 * (n1-1);
                }

                MALLOC(Tris, int, 3*ntri);
                MALLOC(segs, int, 2*nseg);

                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    for (i2 = 1; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            Tris[3*ntri  ] = pvindex[(i1-1)+n1*(i2-1)];
                            Tris[3*ntri+1] = pvindex[(i1  )+n1*(i2-1)];
                            Tris[3*ntri+2] = pvindex[(i1  )+n1*(i2  )];
                            ntri++;

                            Tris[3*ntri  ] = pvindex[(i1  )+n1*(i2  )];
                            Tris[3*ntri+1] = pvindex[(i1-1)+n1*(i2  )];
                            Tris[3*ntri+2] = pvindex[(i1-1)+n1*(i2-1)];
                            ntri++;
                        }
                    }

                    for (i2 = 0; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            segs[2*nseg  ] = pvindex[(i1-1)+n1*(i2)];
                            segs[2*nseg+1] = pvindex[(i1  )+n1*(i2)];
                            nseg++;
                        }
                    }

                    for (i1 = 0; i1 < n1; i1++) {
                        for (i2 = 1; i2 < n2; i2++) {
                            segs[2*nseg  ] = pvindex[(i1)+n1*(i2-1)];
                            segs[2*nseg+1] = pvindex[(i1)+n1*(i2  )];
                            nseg++;
                        }
                    }
                }

                /* two triangles for rendering each quadrilateral */
                status = wv_setData(WV_INT32, 3*ntri, (void*)Tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(Tris);

            /* Triangles */
            } else {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
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
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;
            }

            /* constant triangle colors */
            color[0] = 0.75;
            color[1] = 1.00;
            color[2] = 0.75;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* triangle backface color */
            color[0] = 0.50;
            color[1] = 0.50;
            color[2] = 0.50;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* segment indices */
            status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            FREE(segs);

            /* segment colors */
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* get Attributes for the Face */
            SPLINT_CHECK_FOR_NULL(efaces);
            eface  = efaces[iface-1];
            status = EG_attributeNum(eface, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iface, status);
            }

            /* add Face to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eface, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iface, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Edges within the current Body */
        for (iedge = 1; iedge <= nedge; iedge++) {
            nitems = 0;

            status = EG_getTessEdge(etess, iedge, &npnt, &xyz, &t);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTessEdge(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* skip if a degenerate Edge */
            if (npnt < 2) continue;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Edge %d", bname, iedge);
            attrs = WV_ON;

            /* vertices */
            status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* segments (bias-1) */
            MALLOC(ivrts, int, 2*(npnt-1));

            for (ipnt = 0; ipnt < npnt-1; ipnt++) {
                ivrts[2*ipnt  ] = ipnt + 1;
                ivrts[2*ipnt+1] = ipnt + 2;
            }

            status = wv_setData(WV_INT32, 2*(npnt-1), (void*)ivrts, WV_INDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* line colors */
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* points */
            MALLOC(ivrts, int, npnt);

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ivrts[ipnt] = ipnt + 1;
            }

            status = wv_setData(WV_INT32, npnt, (void*)ivrts, WV_PINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* point colors */
            color[0] = 0;
            color[1] = 0;
            color[2]= 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 2 (does not work for ANGLE) */
                cntxt->gPrims[igprim].lWidth = 2.0;

                /* make point size 5 */
                cntxt->gPrims[igprim].pSize  = 5.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = npnt - 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/ESP->sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_addArrowHeads(%d,%d) -> status=%d", ibody, iedge, status);
                }
            }

            SPLINT_CHECK_FOR_NULL(eedges);
            eedge  = eedges[iedge-1];
            status = EG_attributeNum(eedge, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* add Edge to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eedge, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Nodes within the current Body */
        for (inode = 1; inode <= nnode; inode++) {
            nitems = 0;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Node %d", bname, inode);

            /* default for NodeBodys is to turn the Node on */
            if (nedge == 1 && nnode == 1) {
                attrs = WV_ON;
            } else {
                attrs = 0;
            }

            SPLINT_CHECK_FOR_NULL(enodes);
            enode  = enodes[inode-1];

            /* two copies of vertex (to get around possible wv limitation) */
            status = EG_getTopology(enode, &eref, &oclass, &mtype,
                                    xyz_dum, &nchild, &echilds, &senses);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTopology(%d,%d) -> status=%d", ibody, inode, status);
            }

            xyz_dum[3] = xyz_dum[0];
            xyz_dum[4] = xyz_dum[1];
            xyz_dum[5] = xyz_dum[2];

            status = wv_setData(WV_REAL64, 2, (void*)xyz_dum, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* node color (black) */
            color[0] = 0;   color[1] = 0;   color[2] = 0;

            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].pSize = 6.0;
            }

            status = EG_attributeNum(enode, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* add Node to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[",  gpname);
            } else {
                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(enode, iattr, &attrName, &attrType, &attrLen,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (attrType == ATTRCSYS) continue;

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\",\"", attrName);

                if        (attrType == ATTRINT) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %d", tempIlist[i]);
                    }
                } else if (attrType == ATTRREAL) {
                    for (i = 0; i < attrLen ; i++) {
                        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %f", tempRlist[i]);
                    }
                } else if (attrType == ATTRSTRING) {
                    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, " %s ", tempClist);
                }

                addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");
        }

        /* loop through the Csystems associated with the current Body */
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            nitems = 0;

            status = EG_attributeGet(ebody, iattr, &attrName, &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_attributeGet -> status=%d", status);
            }

            if (attrType != ATTRCSYS) continue;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Csys %s", bname, attrName);
            attrs = WV_ON | WV_SHADING | WV_ORIENTATION;

            /* vertices */
            axis[ 0] = tempRlist[attrLen  ];
            axis[ 1] = tempRlist[attrLen+1];
            axis[ 2] = tempRlist[attrLen+2];
            axis[ 3] = tempRlist[attrLen  ] + tempRlist[attrLen+ 3];
            axis[ 4] = tempRlist[attrLen+1] + tempRlist[attrLen+ 4];
            axis[ 5] = tempRlist[attrLen+2] + tempRlist[attrLen+ 5];
            axis[ 6] = tempRlist[attrLen  ];
            axis[ 7] = tempRlist[attrLen+1];
            axis[ 8] = tempRlist[attrLen+2];
            axis[ 9] = tempRlist[attrLen  ] + tempRlist[attrLen+ 6];
            axis[10] = tempRlist[attrLen+1] + tempRlist[attrLen+ 7];
            axis[11] = tempRlist[attrLen+2] + tempRlist[attrLen+ 8];
            axis[12] = tempRlist[attrLen  ];
            axis[13] = tempRlist[attrLen+1];
            axis[14] = tempRlist[attrLen+2];
            axis[15] = tempRlist[attrLen  ] + tempRlist[attrLen+ 9];
            axis[16] = tempRlist[attrLen+1] + tempRlist[attrLen+10];
            axis[17] = tempRlist[attrLen+2] + tempRlist[attrLen+11];

            status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
            }

            wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
            nitems++;

            /* line colors (multiple colors requires that WV_SHADING be set above) */
            color[ 0] = 1;   color[ 1] = 0;   color[ 2] = 0;
            color[ 3] = 1;   color[ 4] = 0;   color[ 5] = 0;
            color[ 6] = 0;   color[ 7] = 1;   color[ 8] = 0;
            color[ 9] = 0;   color[10] = 1;   color[11] = 0;
            color[12] = 0;   color[13] = 0;   color[14] = 1;
            color[15] = 0;   color[16] = 0;   color[17] = 1;
            status = wv_setData(WV_REAL32, 6, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 1 */
                cntxt->gPrims[igprim].lWidth = 1.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/ESP->sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_addArrowHeads -> status=%d", status);
                }
            }

            /* add Csystem to meta data (if there is room) */
            addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "\"%s\":[],", gpname);
        }

        if (enodes != NULL) EG_free(enodes);
        if (eedges != NULL) EG_free(eedges);
        if (efaces != NULL) EG_free(efaces);

        enodes = NULL;
        eedges = NULL;
        efaces = NULL;
    }

    /* axes */
    nitems = 0;

    /* name and attributes */
    snprintf(gpname, MAX_STRVAL_LEN-1, "Axes");
    attrs = 0;

    /* vertices */
    axis[ 0] = MIN(2*bigbox[0]-bigbox[3], 0);   axis[ 1] = 0;   axis[ 2] = 0;
    axis[ 3] = MAX(2*bigbox[3]-bigbox[0], 0);   axis[ 4] = 0;   axis[ 5] = 0;

    axis[ 6] = 0;   axis[ 7] = MIN(2*bigbox[1]-bigbox[4], 0);   axis[ 8] = 0;
    axis[ 9] = 0;   axis[10] = MAX(2*bigbox[4]-bigbox[1], 0);   axis[11] = 0;

    axis[12] = 0;   axis[13] = 0;   axis[14] = MIN(2*bigbox[2]-bigbox[5], 0);
    axis[15] = 0;   axis[16] = 0;   axis[17] = MAX(2*bigbox[5]-bigbox[2], 0);
    status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
    }

    wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
    nitems++;

    /* line color */
    color[0] = 0.7;   color[1] = 0.7;   color[2] = 0.7;
    status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
    }
    nitems++;

    /* make graphic primitive */
    igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
    if (igprim < 0) {
        SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
    } else {
        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].lWidth = 1.0;
    }

    /* finish the scene graph meta data */
    sgMetaDataUsed--;
    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "}");

    /* send the scene graph meta data if it has not already been sent */
    if (STRLEN(sgMetaData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgData: %s", sgMetaData);
        wv_broadcastText(sgMetaData);

        /* nullify meta data so that it does not get sent again */
        sgMetaData[0]  = '\0';
        sgMetaDataUsed = 0;
    }

    if (STRLEN(sgFocusData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgFocus: %s", sgFocusData);
        wv_broadcastText(sgFocusData);
    }

cleanup:

    /* release the mutex so that another thread can access the scene graph */
    if (ESP != NULL) {
        EMP_LockRelease(ESP->sgMutex);
    }

    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    FREE(segs    );
    FREE(ivrts   );
    FREE(Tris    );

    FREE(sgFocusData);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   buildSceneGraphBOUND - make a scene graph for wv                  */
/*                                                                     */
/***********************************************************************/

static int
buildSceneGraphBOUND(esp_T  *ESP,
                     char   boundName[],
                     char   aimName[],
                     char   dataName[])
{
    int       status = SUCCESS;         /* return status */

    int       ibody, iface, nerror, i, rank, npnt, ipnt, igprim;
    int       attrs, nitems;
    int       ntri, *tris=NULL, nseg, *segs=NULL, ndtris, *dtris=NULL, ndsegs, *dsegs=NULL;
    float     value, color[18], *pcolors=NULL;
    double    bigbox[6], size, axis[18];
    double    *xyz, *data;
    char      gpname[MAX_STRVAL_LEN];
    char      *units;

    capsObj   boundObj, vsetObj, xyzObj, dsetObj, link, temp;
    capsErrs  *errors;
    enum capsoType   otype;
    enum capssType   stype;
    enum capsfType   ftype;
    enum capsdMethod dmethod;
    char      *name;
    capsOwn   last;
    wvContext *cntxt;
    wvData    items[6];

/* global variables associated with scene graph meta-data */
#define MAX_METADATA_CHUNK 32000
    int       sgMetaDataSize = 0;
    int       sgMetaDataUsed = 0;
    char      *sgMetaData    = NULL;
    char      *sgFocusData   = NULL;

    ROUTINE(buildSceneGraphBOUND);

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        SPRINT0(0, "WARNING:: not running via serveESP");
        goto cleanup;
    } else {
        cntxt = ESP->cntxt;
    }

    /* find the BOUND with the given name */
    status = caps_childByName(ESP->CAPS, BOUND, NONE, boundName, &boundObj, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0) {
        SPRINT3(0, "ERROR:: caps_childByName(%s) -> status=%d, nerror=%d", boundName, status, nerror);
        goto cleanup;
    }

    /* find the VERTEXSET with the given name */
    status = caps_childByName(boundObj, VERTEXSET, CONNECTED, aimName, &vsetObj, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0) {
        SPRINT3(0, "ERROR:: caps_childByName(%s) -> status=%d, nerror=%d", aimName, status, nerror);
        goto cleanup;
    }

    /* get the coordinate DATASET associated with this VERTEXSET */
    status = caps_childByName(vsetObj, DATASET, NONE, "xyz", &xyzObj, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0) {
        SPRINT3(0, "ERROR:: caps_childByName(%s) -> status=%d, nerror=%d", "xyz", status, nerror);
        goto cleanup;
    }

    /* find the DATASET with the given name */
    status = caps_childByName(vsetObj, DATASET, NONE, dataName, &dsetObj, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0) {
        SPRINT3(0, "ERROR:: caps_childByName(%s) -> status=%d, nerror=%d", dataName, status, nerror);
        goto cleanup;
    }

    sgMetaDataSize = MAX_METADATA_CHUNK;

    MALLOC(sgMetaData,  char, sgMetaDataSize);
    MALLOC(sgFocusData, char, MAX_STRVAL_LEN);

    /* set the mutex.  if the mutex wss already set, wait until released */
    EMP_LockSet(ESP->sgMutex);

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    /* get the coordinates associated with this VERTEXSET */
    status = caps_getData(xyzObj, &npnt, &rank, &xyz, &units, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0 || rank != 3) {
        SPRINT3(0, "ERROR:: caps_getData(xyz) -> status=%d, rank=%d, nerror=%d", status, rank, nerror);
    }

    status = caps_getData(dsetObj, &npnt, &rank, &data, &units, &nerror, &errors);
    (void) caps_printErrors(NULL, nerror, errors);
    if (status != SUCCESS || nerror > 0) {
        SPRINT2(0, "ERROR:: caps_getData(data) -> status=%d, nerror=%d", status, nerror);
    }

    /* open the key */
    lims[0] = +HUGEQ;
    lims[1] = -HUGEQ;
    
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (rank == 1) {
            value = data[ipnt];
        } else {
            value = 0;
            for (i = 0; i < rank; i++) {
                value += data[ipnt*rank+i] * data[ipnt*rank+i];
            }
            value = sqrt(value);
        }
        if (value < lims[0]) lims[0] = value;
        if (value > lims[1]) lims[1] = value;
    }

    status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], dataName);
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setKey -> status=%d", status);
    }
    wv_broadcastText("setWvKey|on|");

    /* get the triangulations associated with the VERTEXSET */
    tris  = NULL;
    segs  = NULL;
    dtris = NULL;
    dsegs = NULL;

    status = caps_getTriangles(vsetObj, &ntri, &tris, &nseg, &segs, &ndtris, &dtris, &ndsegs, &dsegs);
        SPRINT2(0, "ERROR:: caps_triangulate -> status=%d, nerror=%d", status, nerror);
    if (status != SUCCESS || nerror > 0) {
    }

    /* find the values needed to adjust the vertices */
    bigbox[0] = bigbox[1] = bigbox[2] = +HUGEQ;
    bigbox[3] = bigbox[4] = bigbox[5] = -HUGEQ;

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (xyz[3*ipnt  ] < bigbox[0]) bigbox[0] = xyz[3*ipnt  ];
        if (xyz[3*ipnt+1] < bigbox[1]) bigbox[1] = xyz[3*ipnt+1];
        if (xyz[3*ipnt+2] < bigbox[2]) bigbox[2] = xyz[3*ipnt+2];
        if (xyz[3*ipnt  ] > bigbox[3]) bigbox[3] = xyz[3*ipnt  ];
        if (xyz[3*ipnt+1] > bigbox[4]) bigbox[4] = xyz[3*ipnt+1];
        if (xyz[3*ipnt+2] > bigbox[5]) bigbox[5] = xyz[3*ipnt+2];
    }

                                    size = bigbox[3] - bigbox[0];
    if (size < bigbox[4]-bigbox[1]) size = bigbox[4] - bigbox[1];
    if (size < bigbox[5]-bigbox[2]) size = bigbox[5] - bigbox[2];

    ESP->sgFocus[0] = (bigbox[0] + bigbox[3]) / 2;
    ESP->sgFocus[1] = (bigbox[1] + bigbox[4]) / 2;
    ESP->sgFocus[2] = (bigbox[2] + bigbox[5]) / 2;
    ESP->sgFocus[3] = size;

    /* generate the scene graph focus data */
    snprintf(sgFocusData, MAX_STRVAL_LEN-1, "sgFocus|[%20.12e,%20.12e,%20.12e,%20.12e]",
             ESP->sgFocus[0], ESP->sgFocus[1], ESP->sgFocus[2], ESP->sgFocus[3]);

    /* initialize the scene graph meta data */
    sgMetaData[ 0] = '\0';
    sgMetaDataUsed = 0;

    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "sgData|{");

    /* loop through the Faces within the current Body */
    ibody = 0;
    iface = 0;
    while (1) {
        nitems = 0;

        /* name and attributes */
        if (ibody == 0) {
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Face %d", aimName, iface);
        } else {
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Face %d", name, iface);
        }
        attrs = WV_ON | WV_SHADING | WV_ORIENTATION;

        /* vertices */
        status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
        if (status != SUCCESS) {
            ibody=0;
            SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
        }

        wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
        nitems++;

        /* triangles */
        if (ntri > 0) {
            SPLINT_CHECK_FOR_NULL(tris);

            status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* smooth triangle colors */
            MALLOC(pcolors, float, 3*npnt);
            
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                if (rank == 1) {
                    value = data[ipnt];
                } else {
                    value = 0;
                    for (i = 0; i < rank; i++) {
                        value += data[ipnt*rank+i] * data[ipnt*rank+i];
                    }
                    value = sqrt(value);
                }

                spec_col((float)(value), &(pcolors[3*ipnt]));
            }

            status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            FREE(pcolors);

            /* triangle backface color */
            color[0] = 0.50;
            color[1] = 0.50;
            color[2] = 0.50;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;
        }

        /* segment indices */
        if (nseg > 0) {
            SPLINT_CHECK_FOR_NULL(segs);

            status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* segment colors */
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;
        }

        /* make graphic primitive */
        igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
        if (igprim < 0) {
            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
        } else {
            SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

            cntxt->gPrims[igprim].lWidth = 1.0;
        }

        sgMetaDataUsed--;
        addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "],");

        /* we can have at most one link */
        if (ibody == 1) break;

        /* determine if there is another data linked to this dataset */
        status = caps_dataSetInfo(dsetObj, &ftype, &link, &dmethod);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: caps_dataSetInfo(DATASET) -> status=%d", status);
        }

        if (link == NULL) {
            SPRINT0(2, "this dataset is not linked");
            break;
        } else {
            SPRINT0(2, "this dataset is linked");
            ibody++;
            dsetObj = link;
        }

        /* we have a linked dataset, so find its vertexset */
        status = caps_info(dsetObj, &name, &otype, &stype, &link, &vsetObj, &last);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: capsInfo(DATASET) -> status=%d", status);
        }

        /* get the coordinate DATASET associated with this VERTEXSET */
        status = caps_childByName(vsetObj, DATASET, NONE, "xyz", &xyzObj, &nerror, &errors);
        (void) caps_printErrors(NULL, nerror, errors);
        if (status != SUCCESS || nerror > 0) {
            SPRINT3(0, "ERROR:: caps_childByName(%s) -> status=%d, nerror=%d", "xyz", status, nerror);
            goto cleanup;
        }

        /* get the name of the new VERTEXSET */
        status = caps_info(vsetObj, &name, &otype, &stype, &link, &temp, &last);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: capsInfo(VERTEXSET) -> status=%d", status);
        }

        /* get the coordinates associated with this VERTEXSET */
        status = caps_getData(xyzObj, &npnt, &rank, &xyz, &units, &nerror, &errors);
        (void) caps_printErrors(NULL, nerror, errors);
        if (status != SUCCESS || nerror > 0 || rank != 3) {
            SPRINT3(0, "ERROR:: caps_getData(xyz) -> status=%d, rank=%d, nerror=%d", status, rank, nerror);
        }

        status = caps_getData(dsetObj, &npnt, &rank, &data, &units, &nerror, &errors);
        (void) caps_printErrors(NULL, nerror, errors);
        if (status != SUCCESS || nerror > 0) {
            SPRINT2(0, "ERROR:: caps_getData(data) -> status=%d, nerror=%d", status, nerror);
        }

        /* get the triangulations associated with the VERTEXSET */
        if (tris  != NULL) EG_free(tris );
        if (segs  != NULL) EG_free(segs );
        if (dtris != NULL) EG_free(dtris);
        if (dsegs != NULL) EG_free(dsegs);

        tris  = NULL;
        segs  = NULL;
        dtris = NULL;
        dsegs = NULL;

        status = caps_getTriangles(vsetObj, &ntri, &tris, &nseg, &segs, &ndtris, &dtris, &ndsegs, &dsegs);
        SPRINT2(0, "ERROR:: caps_triangulate -> status=%d, nerror=%d", status, nerror);
        if (status != SUCCESS || nerror > 0) {
        }
    }

    /* axes */
    nitems = 0;

    /* name and attributes */
    snprintf(gpname, MAX_STRVAL_LEN-1, "Axes");
    attrs = 0;

    /* vertices */
    axis[ 0] = MIN(2*bigbox[0]-bigbox[3], 0);   axis[ 1] = 0;   axis[ 2] = 0;
    axis[ 3] = MAX(2*bigbox[3]-bigbox[0], 0);   axis[ 4] = 0;   axis[ 5] = 0;

    axis[ 6] = 0;   axis[ 7] = MIN(2*bigbox[1]-bigbox[4], 0);   axis[ 8] = 0;
    axis[ 9] = 0;   axis[10] = MAX(2*bigbox[4]-bigbox[1], 0);   axis[11] = 0;

    axis[12] = 0;   axis[13] = 0;   axis[14] = MIN(2*bigbox[2]-bigbox[5], 0);
    axis[15] = 0;   axis[16] = 0;   axis[17] = MAX(2*bigbox[5]-bigbox[2], 0);
    status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
    }

    wv_adjustVerts(&(items[nitems]), ESP->sgFocus);
    nitems++;

    /* line color */
    color[0] = 0.7;   color[1] = 0.7;   color[2] = 0.7;
    status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
    }
    nitems++;

    /* make graphic primitive */
    igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
    if (igprim < 0) {
        SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
    } else {
        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].lWidth = 1.0;
    }

    /* finish the scene graph meta data */
    sgMetaDataUsed--;
    addToSgMetaData(&sgMetaDataUsed, sgMetaDataSize, sgMetaData, "}");

    /* send the scene graph meta data if it has not already been sent */
    if (STRLEN(sgMetaData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgData: %s", sgMetaData);
        wv_broadcastText(sgMetaData);

        /* nullify meta data so that it does not get sent again */
        sgMetaData[0]  = '\0';
        sgMetaDataUsed = 0;
    }

    if (STRLEN(sgFocusData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgFocus: %s", sgFocusData);
        wv_broadcastText(sgFocusData);
    }

cleanup:

    /* release the mutex so that another thread can access the scene graph */
    if (ESP != NULL) {
        EMP_LockRelease(ESP->sgMutex);
    }

    if (tris  != NULL) EG_free(tris );
    if (segs  != NULL) EG_free(segs );
    if (dtris != NULL) EG_free(dtris);
    if (dsegs != NULL) EG_free(dsegs);

    FREE(sgFocusData);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   addToSgMetaData - add text to sgMetaData (with buffer length protection) */
/*                                                                     */
/***********************************************************************/

static void
addToSgMetaData(int   *sgMetaDataUsed,  /* (both)SG data used so far */
                int    sgMetaDataSize,  /* (in)  size of SG meta data */
                char   sgMetaData[],    /* (in)  SG meta data */
                char   format[],        /* (in)  format specifier */
                 ...)                   /* (in)  variable arguments */
{
    int      status=SUCCESS, newchars;

    va_list  args;

    ROUTINE(addToSgMetaData);

    /* --------------------------------------------------------------- */

    /* make sure sgMetaData is big enough */
    if ((*sgMetaDataUsed)+1024 >= sgMetaDataSize) {
        sgMetaDataSize += 1024 + MAX_METADATA_CHUNK;
        RALLOC(sgMetaData, char, sgMetaDataSize);
    }

    /* set up the va structure */
    va_start(args, format);

    newchars = vsnprintf(sgMetaData+(*sgMetaDataUsed), 1024, format, args);
    (*sgMetaDataUsed) += newchars;

    /* clean up the va structure */
    va_end(args);

cleanup:
    if (status != SUCCESS) {
        SPRINT0(0, "ERROR:: sgMetaData could not be increased");
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   spec_col - return color for a given scalar value                  */
/*                                                                     */
/***********************************************************************/

static void
spec_col(float  scalar,
         float  color[])
{
    int   indx;
    float frac;

    /* --------------------------------------------------------------- */

    if (lims[0] == lims[1]) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
    } else if (scalar <= lims[0]) {
        color[0] = color_map[0];
        color[1] = color_map[1];
        color[2] = color_map[2];
    } else if (scalar >= lims[1]) {
        color[0] = color_map[3*255  ];
        color[1] = color_map[3*255+1];
        color[2] = color_map[3*255+2];
    } else {
        frac  = 255.0 * (scalar - lims[0]) / (lims[1] - lims[0]);
        if (frac < 0  ) frac = 0;
        if (frac > 255) frac = 255;
        indx  = frac;
        frac -= indx;
        if (indx == 255) {
            indx--;
            frac += 1.0;
        }

        color[0] = frac * color_map[3*(indx+1)  ] + (1.0-frac) * color_map[3*indx  ];
        color[1] = frac * color_map[3*(indx+1)+1] + (1.0-frac) * color_map[3*indx+1];
        color[2] = frac * color_map[3*(indx+1)+2] + (1.0-frac) * color_map[3*indx+2];
    }
}
