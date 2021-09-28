/*
 ************************************************************************
 *                                                                      *
 * udpCsm -- recursive call to OpenCSM                                  *
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

/* the number of "input" Bodys */
#define NUMUDPINPUTBODYS 0

/* the number of arguments (specified below) */
#define NUMUDPARGS 4

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* include OpenCSM */
#include "OpenCSM.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME( IUDP  )  ((char   *) (udps[IUDP].arg[0].val))
#define PMTRNAME( IUDP  )  ((char   *) (udps[IUDP].arg[1].val))
#define PMTRVALUE(IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]
#define VOLUME(   IUDP  )  ((double *) (udps[IUDP].arg[3].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", "pmtrname", "pmtrvalue", "volume",  };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRSTRING, ATTRREAL,    -ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,          0,           0,         };
static double argDdefs[NUMUDPARGS] = {0.,         0.,         0.,          0.,        };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"


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

    void    *modl=NULL;
    int     buildTo, builtTo, nbody, count, i, j, k, ipmtr, found, myUdp;
    int     ibody, jbody, nnode, inode, nedge, iedge, nface, iface;
    int     attrtype, attrlen;
    int     outLevel=0;
    CINT    *tempIlist;
    double  mprop[14];
    char    temp[256];
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL;
    void    *save_modl;
    modl_T  *MODL;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename      = %s\n", FILENAME( 0));
    printf("pmtrname      = %s\n", PMTRNAME( 0));
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf("pmtrvalue[%3d]= %f\n", i, PMTRVALUE(0,i));
    }
#endif

    outLevel = EG_setOutLevel(context, 1       );
    (void)     EG_setOutLevel(context, outLevel);

    /* remember pointer to calling progrm's modl */
    status = EG_getUserPointer(context, (void**)(&(save_modl)));
    CHECK_STATUS(EG_getUserPointer);

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    count = 0;
    if (STRLEN(PMTRNAME(0)) > 0) {
        if (PMTRNAME(0)[0] == ';') {
            printf(" udpExecute: pmtrname cannot start with semi-colon\n");
            status = OCSM_ILLEGAL_VALUE;
            goto cleanup;
        }

        /* count number of names in pmtrname (exclude possible ; at beg and end) */
        count = 1;
        for (i = 1; i < STRLEN(PMTRNAME(0))-1; i++) {
            if (PMTRNAME(0)[i] == ';') count++;
        }

        if (udps[0].arg[2].size != count) {
            printf(" udpExecute: pmtrname and pmtrvalue should have same number of entries\n");
            status = OCSM_ILLEGAL_VALUE;
            goto cleanup;
        }
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

    myUdp = numUdp;

#ifdef DEBUG
    printf("filename[ %d]  = %s\n", myUdp, FILENAME( myUdp));
    printf("pmtrname[ %d]  = %s\n", myUdp, PMTRNAME( myUdp));
    for (i = 0; i < udps[myUdp].arg[2].size; i++) {
        printf("pmtrvalue[%3d]= %f\n", i, PMTRVALUE(myUdp,i));
    }
#endif

    if (outLevel > 0) {
        printf("\n>>> diverting to \"%s\"\n", FILENAME(myUdp));
    }

    /* load the .csm file */
    status = ocsmLoad(FILENAME(myUdp), &modl);
    CHECK_STATUS(ocsmLoad);

    SPLINT_CHECK_FOR_NULL(modl);

    MODL = (modl_T *)modl;

    /* make the new MODL use the same context as the caller */
    EG_deleteObject(MODL->context);
    MODL->context = context;

    /* adjust design parameters, which are specified in pairs in
     pmtrname and pmtrvalue */
    j = 0;
    for (i = 0; i < count; i++) {
        k       = 0;
        temp[0] = '\0';

        while (PMTRNAME(myUdp)[j] != ';' && j < STRLEN(PMTRNAME(myUdp))) {
            temp[k++] = PMTRNAME(myUdp)[j];
            temp[k  ] = '\0';

            j++;
        }
        j++;

        found = 0;
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (strcmp(MODL->pmtr[ipmtr].name, temp) == 0 && MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
                MODL->pmtr[ipmtr].value[0] = PMTRVALUE(myUdp,i);
                MODL->pmtr[ipmtr].dot[  0] = 0;
                found = 1;
                break;
            }
        }
        if (found == 0) {
            printf(" udpExecute: problem finding design parameter \"%s\"\n", temp);
            goto cleanup;
        } else if (outLevel > 0) {
            printf("--> changing \"%s\" to %10.5f\n", temp, PMTRVALUE(myUdp,i));
        }
    }

    /* build the MODL */
    buildTo = 0;
    nbody   = 0;
    status = ocsmBuild(modl, buildTo, &builtTo, &nbody, NULL);
    CHECK_STATUS(ocsmBuild);

    ibody = -1;
    for (jbody = MODL->nbody; jbody >= 1; jbody++) {
        if (MODL->body[jbody].onstack == 1) {
            ibody = jbody;
            break;
        }
    }

    if (ibody <= 0) {
        printf(" udpExecute: no Bodys were left on stack\n");
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* make a copy of the last Body on the stack */
    status = EG_copyObject(MODL->body[MODL->nbody].ebody, NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    SPLINT_CHECK_FOR_NULL(*ebody);

    /* remove _hist and __trace__ attributes (for now), which means
       that we will not be able to track sensitivities for the Body
       created by the .csm file */
    status = EG_getBodyTopos(*ebody, NULL, NODE, &nnode, &enodes);
    printf("nnode=%d\n", nnode);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(enodes);

    for (inode = 0; inode < nnode; inode++) {
        status = EG_attributeRet(enodes[inode], "__trace__",
                                 &attrtype, &attrlen, &tempIlist, NULL, NULL);
        if (status == SUCCESS) {
            status = EG_attributeDel(enodes[inode], "__trace__");
            CHECK_STATUS(EG_attributeDel);
        }
    }

    EG_free(enodes);   

    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
    printf("nedge=%d\n", nedge);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(eedges);

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_attributeRet(eedges[iedge], "__trace__",
                                 &attrtype, &attrlen, &tempIlist, NULL, NULL);
        if (status == SUCCESS) {
            status = EG_attributeDel(eedges[iedge], "__trace__");
            CHECK_STATUS(EG_attributeDel);
        }
    }

    EG_free(eedges);

    status = EG_getBodyTopos(*ebody, NULL, FACE, &nface, &efaces);
    printf("nface=%d\n", nface);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(efaces);

    for (iface = 0; iface < nface; iface++) {
        status = EG_attributeRet(efaces[iface], "_hist",
                                 &attrtype, &attrlen, &tempIlist, NULL, NULL);
        if (status == SUCCESS) {
            status = EG_attributeDel(efaces[iface], "_hist");
            CHECK_STATUS(EG_attributeDel);
        }
    }

    for (iface = 0; iface < nface; iface++) {
        status = EG_attributeRet(efaces[iface], "__trace__",
                                 &attrtype, &attrlen, &tempIlist, NULL, NULL);
        if (status == SUCCESS) {
            status = EG_attributeDel(efaces[iface], "__trace__");
            CHECK_STATUS(EG_attributeDel);
        }
    }

    EG_free(efaces);

    if (outLevel > 0) {
        printf("<<< returning from diversion to \"%s\"\n\n", FILENAME(myUdp));
    }

    /* restore user data to original modl */
    status = EG_setUserPointer(context, save_modl);
    CHECK_STATUS(EG_setUserPointer);

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, mprop);
    CHECK_STATUS(EG_getMassProperties);

    VOLUME(0) = mprop[0];

    /* remember this model (Body) */
    udps[myUdp].ebody = *ebody;
    goto cleanup;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (modl != NULL) {
        (void) ocsmFree(modl);
    }

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
    int    iudp, judp;

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
