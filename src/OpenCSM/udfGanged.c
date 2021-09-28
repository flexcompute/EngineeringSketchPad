/*
 ************************************************************************
 *                                                                      *
 * udfGanged -- test ganged Boolean operations                          *
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

#define NUMUDPINPUTBODYS -999
#define NUMUDPARGS          2

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define OP(   IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define TOLER(IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]

static char  *argNames[NUMUDPARGS] = {"op",       "toler",  };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,        };
static double argDdefs[NUMUDPARGS] = {0.,         0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

//#define PRINT_TIMES 1
#ifdef  PRINT_TIMES
#include <sys/time.h>
#include <sys/resource.h>
#endif


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

    int      oclass, mtype, nchild, *senses, i;
    double   data[18];
    ego      eref, *ebodysB=NULL, etools, eresult, *echilds;
    void     *modl;

#ifdef PRINT_TIMES
    struct rusage  beg_rusage, end_rusage;
    struct timeval beg_time,   end_time;
    double         user_time, sys_time, wall_time;
#endif

    ROUTINE(udpExecute);
    
#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("op(   0) = %s\n", OP(   0));
    printf("toler(0) = %f\n", TOLER(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check/process arguments */
    if (strcmp(OP(0), "SUBTRACT") != 0 && strcmp(OP(0), "UNION") != 0   ) {
        printf(" udpExecute: op should be SUBTRACT or UNION\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    }

    /* check that Model was input that contains two or more Bodys */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild < 2) {
        printf(" udpExecute: expecting Model to contain at least two Bodys (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

#ifdef DEBUG
    printf("emodel\n");
    ocsmPrintEgo(emodel);
#endif

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("op(   %d) = %s\n", numUdp, OP(   numUdp));
    printf("toler(%d) = %f\n", numUdp, TOLER(numUdp));
#endif

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* return the modfied Model that contains the two input Bodys */
//$$$    *ebody = emodel;

#ifdef DEBUG
    printf("*ebody\n");
    (void) ocsmPrintEgo(*ebody);
#endif

    /* get pointer to model */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    /* the left Body is ebodys[0] */
    
    /* make a Model of the rest of the Bodys in emodel */

    MALLOC(ebodysB, ego, (nchild-1));

    for (i = 1; i < nchild; i++) {
        status = EG_copyObject(ebodys[i], NULL, &(ebodysB[i-1]));
        CHECK_STATUS(EG_copyObject);
    }

    status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                             nchild-1, ebodysB, NULL, &etools);
    CHECK_STATUS(EG_makeTopology);

    /* perform the ganged boolean operation */

#ifdef PRINT_TIMES
    gettimeofday(&beg_time, NULL);
    getrusage(RUSAGE_SELF, &beg_rusage);
#endif
    if (strcmp(OP(0), "SUBTRACT") == 0) {
        status = EG_generalBoolean(ebodys[0], etools, 1, TOLER(0), &eresult);
        CHECK_STATUS(EG_generalBoolean);
    } else {
        status = EG_generalBoolean(ebodys[0], etools, 3, TOLER(0), &eresult);
        CHECK_STATUS(EG_generalBoolean);
    }
#ifdef PRINT_TIMES
    gettimeofday(&end_time, NULL);
    getrusage(RUSAGE_SELF, &end_rusage);
    user_time = (end_rusage.ru_utime.tv_sec+end_rusage.ru_utime.tv_usec*0.000001)
               -(beg_rusage.ru_utime.tv_sec+beg_rusage.ru_utime.tv_usec*0.000001);
    sys_time  = (end_rusage.ru_stime.tv_sec+end_rusage.ru_stime.tv_usec*0.000001)
               -(beg_rusage.ru_stime.tv_sec+beg_rusage.ru_stime.tv_usec*0.000001);
    wall_time =((long)end_time.tv_sec+(long)end_time.tv_usec*0.000001)
              -((long)beg_time.tv_sec+(long)beg_time.tv_usec*0.000001);

    printf("user:%.3f    sys:%.3f     tot:%.3f     wall=%.3f\n",
           user_time, sys_time, user_time+sys_time, wall_time);
#endif

    status = EG_deleteObject(etools);
    CHECK_STATUS(EG_deleteObject);

    status = EG_getTopology(eresult, &eref, &oclass, &mtype, 
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_copyObject(echilds[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    status = EG_deleteObject(eresult);
    CHECK_STATUS(EG_deleteObject);

    if (nchild != 1) {
        printf(" udpExecute: expecting 1 result Body, got %d\n", nchild);
        status = EGADS_TOPOERR;
        goto cleanup;
    }

    /* return output value(s) */

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;
    goto cleanup;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (ebodysB != NULL) EG_free(ebodysB);

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
