/*
 ************************************************************************************************
 *                                                                                              *
 * udfNacelle -- generate a customizable aircraft engine nacelle                                *
 *                                                                                              *
 *                       Written by Aidan Hoff                                                  *
 *                       Modified by John F. Dannenhoffer, III                                  *
 *                                                                                              *
 ************************************************************************************************
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
 *    MA  02110-1301  USA
 */

/* definitions for udpUtilities */
#define NUMUDPINPUTBODYS 1
#define NUMUDPARGS 7

#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FRAD(     IUDP, I) ((double *) (udps[IUDP].arg[0].val))[I] // E, N, W, S
#define FRAD_SIZE(IUDP   )              udps[IUDP].arg[0].size
#define ARAD(     IUDP, I) ((double *) (udps[IUDP].arg[1].val))[I] // E, N, W, S
#define ARAD_SIZE(IUDP   )              udps[IUDP].arg[1].size
#define FPOW(     IUDP, I) ((double *) (udps[IUDP].arg[2].val))[I] // NE, NW, SW, SE
#define FPOW_SIZE(IUDP   )              udps[IUDP].arg[2].size
#define APOW(     IUDP, I) ((double *) (udps[IUDP].arg[3].val))[I] // NE, NW, SW, SE
#define APOW_SIZE(IUDP   )              udps[IUDP].arg[3].size
#define LENGTH(   IUDP   ) ((double *) (udps[IUDP].arg[4].val))[0]
#define DELTAH(   IUDP   ) ((double *) (udps[IUDP].arg[5].val))[0]
#define RAKEANG(  IUDP   ) ((double *) (udps[IUDP].arg[6].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"f_rad",   "a_rad",   "f_pow",   "a_pow",   "length",  "deltah",  "rakeang", };
static int    argTypes[NUMUDPARGS] = { ATTRREAL,  ATTRREAL,  ATTRREAL,  ATTRREAL,  ATTRREAL,  ATTRREAL,  ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = { 0,         0,         0,         0,         0,         0,         0,        };
static double argDdefs[NUMUDPARGS] = { 0.,        0.,        0.,        0.,        0.,        0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

#define           MAXLOC          0.40
#define           SHARPTE         1
#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define           SIGN(A)         (((A) < 0) ? -1 : (((A) > 0) ? +1 : 0))


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

    double  data[18];

    int     oclass, mtype, nchild, *senses, isect;
    double   leadx,  leady,  leadz;
    double  trailx, traily, trailz;
    double  thetax, thetay, thetaz;
    double  phi, scale;
    double  dirBeg[4], dirEnd[4], matrix[12];
    ego     eref, context, *ebodys, *echilds=NULL, *echilds2;
    ego     etrans, ecross[37], equad[4], ehalf[2], ewhole;

    char    *message=NULL;

    ROUTINE(udpExecute);

    /* -------------------------------------------------------------- */

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

#ifdef DEBUG
    // arguments here
#endif

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 100, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* check arguments */
    if (udps[0].arg[0].size != 4) {
        snprintf(message, 100, "four values needed for foreward radii");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)FRAD(0,0) < 0) {
        snprintf(message, 100, "FRAD(0,0) = %f < 0", FRAD(0,0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)FRAD(0,1) < 0) {
        snprintf(message, 100, "FRAD(0,1) = %f < 0", FRAD(0,1));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)FRAD(0,2) < 0) {
        snprintf(message, 100, "FRAD(0,2) = %f < 0", FRAD(0,2));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)FRAD(0,3) < 0) {
        snprintf(message, 100, "FRAD(0,3) = %f < 0", FRAD(0,3));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size != 4) {
        snprintf(message, 100, "four values needed for aft radii");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)ARAD(0,0) < 0) {
        snprintf(message, 100, "ARAD(0,0) = %f < 0", ARAD(0,0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)ARAD(0,1) < 0) {
        snprintf(message, 100, "ARAD(0,1) = %f < 0", ARAD(0,1));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)ARAD(0,2) < 0) {
        snprintf(message, 100, "ARAD(0,2) = %f < 0", ARAD(0,2));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)ARAD(0,3) < 0) {
        snprintf(message, 100, "ARAD(0,3) = %f < 0", ARAD(0,3));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size != 4) {
        snprintf(message, 100, "four values needed for foreward powers");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size != 4) {
        snprintf(message, 100, "four values needed for aft powers");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        snprintf(message, 100, "length should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)LENGTH(0) < 0) {
        snprintf(message, 100, "length = %f < 0", LENGTH(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "delta_height should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[6].size > 1) {
        snprintf(message, 100, "rake_angle should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    // arguments here
#endif

    /* get the Face from the input Body */
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nchild, &echilds);
    CHECK_STATUS(EG_getBodyTopos);

    if (nchild != 1) {
        snprintf(message, 100, "input Body should contain 1 Face (not %d)", nchild);
        status = EGADS_NODATA;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(echilds);

    /* Creating airfoils at 10 degree intervals */
    for (isect = 0; isect < 25; isect++) {
        
        /* compute the angle ([0, 2*pi] radians) */
        phi = isect * PI / 12;

        /* calculate the position on the super-ellipses */
        if        (isect <= 6){
            leadz  =  FRAD(0,0)*(pow(cos(   phi), (2/FPOW(0,0))));
            leady  =  FRAD(0,1)*(pow(sin(   phi), (2/FPOW(0,0))));

            trailz =  ARAD(0,0)*(pow(cos(   phi), (2/APOW(0,0))));
            traily =  ARAD(0,1)*(pow(sin(   phi), (2/APOW(0,0))));
        } else if (isect <= 12){
            leadz  = -FRAD(0,2)*(pow(cos(PI-phi), (2/FPOW(0,1))));
            leady  =  FRAD(0,1)*(pow(sin(PI-phi), (2/FPOW(0,1))));

            trailz = -ARAD(0,2)*(pow(cos(PI-phi), (2/APOW(0,1))));
            traily =  ARAD(0,1)*(pow(sin(PI-phi), (2/APOW(0,1))));
        } else if (isect <= 18){
            leadz  = -FRAD(0,2)*(pow(cos(phi-PI), (2/FPOW(0,2))));
            leady  = -FRAD(0,3)*(pow(sin(phi-PI), (2/FPOW(0,2))));

            trailz = -ARAD(0,2)*(pow(cos(phi-PI), (2/APOW(0,2))));
            traily = -ARAD(0,3)*(pow(sin(phi-PI), (2/APOW(0,2))));
        } else {
            leadz  =  FRAD(0,0)*(pow(cos(  -phi), (2/FPOW(0,3))));
            leady  = -FRAD(0,3)*(pow(sin(  -phi), (2/FPOW(0,3))));

            trailz =  ARAD(0,0)*(pow(cos(  -phi), (2/APOW(0,3))));
            traily = -ARAD(0,3)*(pow(sin(  -phi), (2/APOW(0,3))));
        }

        /* account for height difference of super-ellipse centers */
        traily += DELTAH(0);

        /* account for the rake angle */
        leadx  = -leady * sin(RAKEANG(0) * PI / 180);
        leady *=          cos(RAKEANG(0) * PI / 180);

        /* account for the length of the nacelle */
        trailx = LENGTH(0);

        /* scale the airfoil */
        scale = sqrt((trailx-leadx) * (trailx-leadx) +
                     (traily-leady) * (traily-leady) +
                     (trailz-leadz) * (trailz-leadz) );

        /* determine how the airfoil should be rotated into place */
        thetaz =  atan((traily-leady)/(trailx-leadx));
        thetay = -atan((trailz-leadz)/(trailx-leadx));

        /* orientation for sections */
        thetax = (90 - 15 * isect) * PI / 180;

        /* compute the transformation matrix */
        matrix[ 0] =  scale *  cos(thetay) * cos(thetaz)                                           ;
        matrix[ 1] =  scale * (sin(thetax) * sin(thetay) - cos(thetax) * cos(thetay) * sin(thetaz));
        matrix[ 2] =  scale * (cos(thetax) * sin(thetay) + cos(thetay) * sin(thetax) * sin(thetaz));
        matrix[ 3] =  leadx                                                                        ;

        matrix[ 4] =  scale *  sin(thetaz)                                                         ;
        matrix[ 5] =  scale *  cos(thetax) * cos(thetaz)                                           ;
        matrix[ 6] = -scale *  cos(thetaz) * sin(thetax)                                           ;
        matrix[ 7] =  leady                                                                        ;

        matrix[ 8] = -scale *  cos(thetaz) * sin(thetay)                                           ;
        matrix[ 9] =  scale * (cos(thetay) * sin(thetax) + cos(thetax) * sin(thetay) * sin(thetaz));
        matrix[10] =  scale * (cos(thetax) * cos(thetay) - sin(thetax) * sin(thetay) * sin(thetaz));
        matrix[11] =  leadz                                                                        ;

        /* make the transform object */
        status = EG_makeTransform(context, matrix, &etrans);
        CHECK_STATUS(EG_makeTransform);

        status = EG_copyObject(echilds[0], etrans, &ecross[isect]);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(etrans);
        CHECK_STATUS(EG_deleteObject);
    }

    EG_free(echilds);

    /* make the blends */

    /* quadrant I (+y,+z)*/
    if (FPOW(0,0) < 2 || APOW(0,0) < 2) {
        status = EG_blend(7, &ecross[0], NULL, NULL, &equad[0]);
        CHECK_STATUS(EG_blend);
    } else {
        dirBeg[0] = 1; dirBeg[1] = 0; dirBeg[2] = +1; dirBeg[3] =  0;
        dirEnd[0] = 1; dirEnd[1] = 0; dirEnd[2] =  0; dirEnd[3] = +1;

        status = EG_blend(7, &ecross[0], dirBeg, dirEnd, &equad[0]);
        CHECK_STATUS(EG_blend);
    }

    /* quadrant II (+y,-z) */
    if (FPOW(0,1) < 2 || APOW(0,1) < 2) {
        status = EG_blend(7, &ecross[6], NULL, NULL, &equad[1]);
        CHECK_STATUS(EG_blend);
    } else {
        dirBeg[0] = 1; dirBeg[1] = 0; dirBeg[2] =  0; dirBeg[3] = -1;
        dirEnd[0] = 1; dirEnd[1] = 0; dirEnd[2] = +1; dirEnd[3] =  0;

        status = EG_blend(7, &ecross[6], dirBeg, dirEnd, &equad[1]);
        CHECK_STATUS(EG_blend);
    }

    /* quadrant III (-y,-z) */
    if (FPOW(0,2) < 2 || APOW(0,2) < 2) {
        status = EG_blend(7, &ecross[12], NULL, NULL, &equad[2]);
        CHECK_STATUS(EG_blend);
    } else {
        dirBeg[0] = 1; dirBeg[1] = 0; dirBeg[2] = -1; dirBeg[3] =  0;
        dirEnd[0] = 1; dirEnd[1] = 0; dirEnd[2] =  0; dirEnd[3] = -1;

        status = EG_blend(7, &ecross[12], dirBeg, dirEnd, &equad[2]);
        CHECK_STATUS(EG_blend);
    }

    /* quadrant IV (-y,+z) */
    if (FPOW(0,3) < 2 || APOW(0,3) < 2) {
        status = EG_blend(7, &ecross[18], NULL, NULL, &equad[3]);
        CHECK_STATUS(EG_blend);
    } else {
        dirBeg[0] = 1; dirBeg[1] = 0; dirBeg[2] =  0; dirBeg[3] = +1;
        dirEnd[0] = 1; dirEnd[1] = 0; dirEnd[2] = -1; dirEnd[3] =  0;

        status = EG_blend(7, &ecross[18], dirBeg, dirEnd, &equad[3]);
        CHECK_STATUS(EG_blend);
    }

    /* union the four quadrants into one body */
    status = EG_generalBoolean(equad[0], equad[1], FUSION, 0.0, &ehalf[0]);
    CHECK_STATUS(EG_generalBoolean);

    status = EG_deleteObject(equad[0]);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(equad[1]);
    CHECK_STATUS(EG_deleteObject);

    status = EG_generalBoolean(equad[2], equad[3], FUSION, 0.0, &ehalf[1]);
    CHECK_STATUS(EG_generalBoolean);

    status = EG_deleteObject(equad[2]);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(equad[3]);
    CHECK_STATUS(EG_deleteObject);

    status = EG_generalBoolean(ehalf[0], ehalf[1], FUSION, 0.0, &ewhole);
    CHECK_STATUS(EG_generalBoolean);

    status = EG_deleteObject(ehalf[0]);
    CHECK_STATUS(EG_deleteObject);

    status = EG_deleteObject(ehalf[1]);
    CHECK_STATUS(EG_deleteObject);

    /* extract the Body from the Model */
    status = EG_getTopology(ewhole, &eref, &oclass, &mtype, data, &nchild, &echilds2, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_copyObject(echilds2[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    SPLINT_CHECK_FOR_NULL(*ebody);
    
    status = EG_deleteObject(ewhole);
    CHECK_STATUS(EG_deleteObject);

    status = EG_attributeAdd(*ebody, "__markFaces__", ATTRSTRING, 1, NULL, NULL, "true");
    CHECK_STATUS(EG_attributeAdd);

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
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
 *   udpSensitivity - return sensitivity derivatives for the            *
 *                    "real" argument                                   *
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
    int    iudp, judp;

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
