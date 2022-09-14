/*
 ************************************************************************
 *                                                                      *
 * udfPrintBbox -- print bounding box info for Body on top of the stack *
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

#define NUMUDPARGS 1
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"foo",   };
static int    argTypes[NUMUDPARGS] = {ATTRINT, };
static int    argIdefs[NUMUDPARGS] = {0,       };
static double argDdefs[NUMUDPARGS] = {0.,      };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */


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
    int     inode, nnode, iedge, nedge, iface, nface;
    double  data[18], bbox[6];
    ego     context, eref, *ebodys, *echilds, *enodes, *eedges, *efaces;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

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

    /* check arguments */

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
#endif

    /* make a copy of the Body (so that it does not get removed
     when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* Body bbox */
    status = EG_getBoundingBox(*ebody, bbox);
    CHECK_STATUS(EG_getBoundingBox);
    printf("    Body        x:%10.5f %10.5f   y:%10.5f %10.5f   z:%10.5f %10.5f\n\n",
               bbox[0], bbox[3], bbox[1], bbox[4], bbox[2], bbox[5]);
    

    /* Node bboxes */
    status = EG_getBodyTopos(*ebody, NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);
    
    for (inode = 0; inode < nnode; inode++) {
        status = EG_getTopology(enodes[inode], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        printf("    Node %4d   x:%10.5f              y:%10.5f              z:%10.5f\n",
               inode+1, data[0], data[1], data[2]);
    }
    printf("\n");
    EG_free(enodes);

    /* Edge bboxes */
    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);
    
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getBoundingBox(eedges[iedge], bbox);
        CHECK_STATUS(EG_getBoundingBox);

        printf("    Edge %4d   x:%10.5f %10.5f   y:%10.5f %10.5f   z:%10.5f %10.5f\n",
               iedge+1, bbox[0], bbox[3], bbox[1], bbox[4], bbox[2], bbox[5]);
    }
    printf("\n");
    EG_free(eedges);

    /* Face bboxes */
    status = EG_getBodyTopos(*ebody, NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);
    
    for (iface = 0; iface < nface; iface++) {
        status = EG_getBoundingBox(efaces[iface], bbox);
        CHECK_STATUS(EG_getBoundingBox);

        printf("    Face %4d   x:%10.5f %10.5f   y:%10.5f %10.5f   z:%10.5f %10.5f\n",
               iface+1, bbox[0], bbox[3], bbox[1], bbox[4], bbox[2], bbox[5]);
    }
    printf("\n");
    EG_free(efaces);

    /* set the output value */

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
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
