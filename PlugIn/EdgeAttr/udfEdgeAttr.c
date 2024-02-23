/*
 ************************************************************************
 *                                                                      *
 * udfEdgeAttr -- adds an Edge Attribute selected by the closest pts    *
 *                                                                      *
 ************************************************************************
 */

#define NUMUDPARGS 5
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values */
#define ATTRNAME( IUDP)   ((char   *) (udps[IUDP].arg[0].val))
#define ATTRSTR(  IUDP)   ((char   *) (udps[IUDP].arg[1].val))
#define ATTRINTS( IUDP,I) ((int    *) (udps[IUDP].arg[2].val))[I]
#define ATTRREALS(IUDP,I) ((double *) (udps[IUDP].arg[3].val))[I]
#define XYZS(     IUDP,I) ((double *) (udps[IUDP].arg[4].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"attrname",  "attrstr",   "attrints",
                                      "attrreals", "xyzs"   };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING,  ATTRSTRING,  ATTRINT,
                                      ATTRREAL,    ATTRREAL };
static int    argIdefs[NUMUDPARGS] = {0,           0,           -999999,
                                      0,           0        };
static double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,
                                      -999999.,    0.       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

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

    int     i, j, npts, oclass, mtype, nchild, aType, aLen, nedge, index, *senses;
    char    *message = NULL;
    double  data[4], uv[2], xyz[3], d, dist, mdist, *xyza;
    ego     eref, *ebodys, *edges = NULL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("attrname( 0) = %s\n",  ATTRNAME( 0));
    printf("attrstr(  0) = %s\n",  ATTRSTR(  0));
    printf("attrints( 0) = %d",    ATTRINTS( 0,0));
    for (i = 1; i < udps[0].arg[2].size; i++)
        printf(" %d", ATTRINTS(0,i));
    printf("\n");
    printf("attrreals(0) = %lf",   ATTRREALS(0,0));
    for (i = 1; i < udps[0].arg[3].size; i++)
        printf(" %lf", ATTRREALS(0,i));
    printf("\n");
    printf("xyzs(     0) = %lf",   XYZS(     0,0));
    for (i = 1; i < udps[0].arg[4].size; i++)
        printf(" %lf", XYZS(0,i));
    printf("\n");
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" EdgeAttr udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        printf(" EdgeAttr udpExecute: expecting Model to contain one Body (not %d)\n",
               nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* check arguments */
    if (STRLEN(ATTRNAME(0)) == 0) {
        snprintf(message, 100, "attrname must be set!");
        status = EGADS_ATTRERR;
        goto cleanup;
    }
    aLen = 0;
    if (STRLEN(ATTRSTR(0)) != 0) {
        aLen  = STRLEN(ATTRSTR(0));
        aType = ATTRSTRING;
    }
    if ((udps[0].arg[2].size != 1) || (ATTRINTS(0,0) != -999999)) {
        if (aLen != 0) {
            snprintf(message, 100,
                     "Only one of attrstr, attrints or attreals must be set!");
            status = EGADS_ATTRERR;
            goto cleanup;
        }
        aLen  = udps[0].arg[2].size;
        aType = ATTRINT;
    }
    if ((udps[0].arg[3].size != 1) || (ATTRREALS(0,0) != -999999.)) {
        if (aLen != 0) {
            snprintf(message, 100,
                     "Only one of attrstr, attrints or attreals must be set!");
            status = EGADS_ATTRERR;
            goto cleanup;
        }
        aLen  = udps[0].arg[3].size;
        aType = ATTRREAL;
    }
    if (aLen == 0) {
        snprintf(message, 100, "attrstr, attrints or attreals must be set!");
        status = EGADS_ATTRERR;
        goto cleanup;
    }
    npts = udps[0].arg[4].size/3;
    if ((npts == 0) || (npts*3 != udps[0].arg[4].size)) {
        snprintf(message, 100, "xyzs length must be a factor of 3 but is %d!",
                 udps[0].arg[4].size);
        status = EGADS_INDEXERR;
        goto cleanup;
    }
#ifdef DEBUG
    printf("attribute is %d with length of %d\n", aType, aLen);
#endif

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("attrname( %d) = %s\n",  numUdp, ATTRNAME( numUdp));
    printf("attrstr(  %d) = %s\n",  numUdp, ATTRSTR(  numUdp));
    printf("attrints( %d) = %d",    numUdp, ATTRINTS( numUdp,0));
    for (i = 1; i < udps[numUdp].arg[2].size; i++)
        printf(" %d", ATTRINTS(numUdp,i));
    printf("\n");
    printf("attrreals(%d) = %lf",   numUdp, ATTRREALS(numUdp,0));
    for (i = 1; i < udps[numUdp].arg[3].size; i++)
        printf(" %lf", ATTRREALS(numUdp,i));
    printf("\n");
    printf("xyzs(     %d) = %lf",   numUdp, XYZS(     numUdp,0));
    for (i = 1; i < udps[numUdp].arg[4].size; i++)
        printf(" %lf", XYZS(numUdp,i));
    printf("\n");
#endif

    /* make a copy of the Body (so that it does not get removed
       when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* add a special Attribute to the Body to tell OpenCSM that there
       is no topological change and hence it should not adjust the
       Attributes on the Body in finishBody() */
    status = EG_attributeAdd(*ebody, "__noTopoChange__", ATTRSTRING,
                             0, NULL, NULL, "udfEdgeAttr");
    CHECK_STATUS(EG_attributeAdd);
  
    /* find the closest Edge to the input Points */
    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &edges);
    CHECK_STATUS(EG_getBodyTopos);
  
    index = -1;
    dist  =  1.e200;
    xyza  = (double *) udps[numUdp].arg[4].val;
    for (i = 0; i < nedge; i++) {
        if (edges[i]->mtype == DEGENERATE) continue;
        mdist = 0.0;
        for (j = 0; j < npts; j++) {
            status = EG_invEvaluate(edges[i], &xyza[3*j], uv, xyz);
            CHECK_STATUS(EG_invEvaluate);
            d = sqrt((xyz[0]-xyza[3*j  ])*(xyz[0]-xyza[3*j  ]) +
                     (xyz[1]-xyza[3*j+1])*(xyz[1]-xyza[3*j+1]) +
                     (xyz[2]-xyza[3*j+2])*(xyz[2]-xyza[3*j+2]));
            if (d > mdist) mdist = d;
        }
        if (mdist < dist) {
            dist  = mdist;
            index = i;
        }
    }
    if (index == -1) {
        snprintf(message, 100, "No Edge Selected!");
        status = EGADS_EMPTY;
        goto cleanup;
    }
#ifdef DEBUG
    printf(" Edge %d: max distance = %le\n", index+1, dist);
#endif
    status = EG_attributeAdd(edges[index], ATTRNAME(numUdp), aType, aLen,
                             (int *)    udps[numUdp].arg[2].val,
                             (double *) udps[numUdp].arg[3].val,
                             (char *)   udps[numUdp].arg[1].val);
    CHECK_STATUS(EG_attributeAdd);

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

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
    if (edges != NULL) EG_free(edges);

    return status;
}


/*
 *****************************************************************************
 *                                                                           *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                           *
 *****************************************************************************
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

    /* this routine is not needed */
    return EGADS_NOLOAD;
}
