/*
 ************************************************************************
 *                                                                      *
 * udpUtilities -- utility routine for UDPs                             *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2025  John F. Dannenhoffer, III (Syracuse University)
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

#include "OpenCSM.h"
#include "common.h"

#define STRLEN(A)   (int)strlen(A)

static int  findInstance(ego ebody, int numudp, udp_T udps[]);

/*
  functions defined in the UDP/UDF:
    udpExecute      execute teh primitive
    udpSensitivity  return sensitivity
    udpMesh         return mesh associated with primitive (optional)
  functions defined below:
    udpErrorStr     print error string
    udpInitialize   initialize and get info about the list of areguments
    udpNumBodys     number of Bodys expected
    udpBodyList     list of Bodys input to a UDF
    udpReset        reset the arguments to their defaults
    udpSet          set an argument
    udpGet          return an output parameter
    udpVel          set velocity of an argument
    udpNumCals      number of calls to udpSet since last udpExecute
    udpNumVels      number of calls to udpVel since last udpSensitivity
    udpClean        clean the udp cache
    udpMesh         return mesh associated with primitive
    udpFree         free all memory associated with UDP/UDF
    CacheUdp        cache arguments of current instance
    findInstance    find UDP instance that matches ebody
    udpCleanupAll   free up UDP/UDF dispatch table
  macros that can be defined in UDP/UDP for the private data:
    FREEUDPDATA     free private data
    COPYUDPDATA     copy private data
*/


/*
 ************************************************************************
 *                                                                      *
 *   udpErrorStr - print error string                                   *
 *                                                                      *
 ************************************************************************
 */

/*@null@*/ char *
udpErrorStr(int stat)                   /* (in)  status number */
{
    char *string;                       /* (out) error message */

    /* --------------------------------------------------------------- */

    string = (char *) EG_alloc(25*sizeof(char));
    if (string == NULL) {
        return string;
    }
    snprintf(string, 25, "EGADS status = %d", stat);

    return string;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpInitialize - initialize and get info about the list of arguments*
 *                                                                      *
 ************************************************************************
 */

int
udpInitialize(int    *nArgs,            /* (out) number of arguments */
              char   **namex[],         /* (out) array  of argument names */
              int    *typex[],          /* (out) array  of argument types */
              int    *idefault[],       /* (out) array  of argument integer defaults */
              double *ddefault[],       /* (out) array  of argument double  defaults */
              udp_T  *Udps[])           /* (both)array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iarg, ival, *pint;
    double *p1dbl, *p2dbl;
    char   *pchar;
    udp_T  *udps;
    void   *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(udpInitialize);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpInitialize()\n");
#endif

    /* make the initial array that holds instances of this udp */
    if (*Udps == NULL) {

        MALLOC(*Udps, udp_T, 1);
        udps = *Udps;

        udps[0].narg = NUMUDPARGS;
        udps[0].arg  = NULL;
        RALLOC(udps[0].arg, udparg_T, NUMUDPARGS);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            udps[0].arg[iarg].type = argTypes[iarg];
            udps[0].arg[iarg].size = 1;
            udps[0].arg[iarg].nrow = 1;
            udps[0].arg[iarg].ncol = 1;
            udps[0].arg[iarg].val  = NULL;
            udps[0].arg[iarg].dot  = NULL;

            if        (argTypes[iarg] == +ATTRSTRING) {
                MALLOC(udps[0].arg[iarg].val, char, MAX_EXPR_LEN);
            } else if (argTypes[iarg] == +ATTRFILE) {
                MALLOC(udps[0].arg[iarg].val, char, MAX_EXPR_LEN);
            } else if (argTypes[iarg] == +ATTRINT ||
                       argTypes[iarg] == -ATTRINT   ) {
                MALLOC(udps[0].arg[iarg].val, int, 1);
            } else if (argTypes[iarg] == +ATTRREAL ||
                       argTypes[iarg] == -ATTRREAL   ) {
                MALLOC(udps[0].arg[iarg].val, double, 1);
            } else if (argTypes[iarg] == +ATTRREALSEN ||
                       argTypes[iarg] == -ATTRREALSEN   ) {
                MALLOC(udps[0].arg[iarg].val, double, 1);
                MALLOC(udps[0].arg[iarg].dot, double, 1);
            } else if (argTypes[iarg] == +ATTRREBUILD) {
            } else if (argTypes[iarg] == +ATTRRECYCLE) {
            } else if (argTypes[iarg] == 0) {
                MALLOC(udps[0].arg[iarg].val, double, 1);
            } else {
                printf("bad argType[%d]=%d in udpInitialize\n", iarg, argTypes[iarg]);
                exit(0);
            }
        }
        udps[0].data = NULL;
    }
    udps = *Udps;

    /* initialize the array elements that will hold the "current" settings */
    udps[0].ebody    = NULL;
    udps[0].ndotchg  = 1;
    udps[0].bodyList = NULL;

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if        (argTypes[iarg] == +ATTRSTRING) {
            pchar  = (char   *) (udps[0].arg[iarg].val);

            strcpy(pchar, "");
        } else if (argTypes[iarg] == +ATTRFILE) {
            pchar  = (char   *) (udps[0].arg[iarg].val);

            strcpy(pchar, "");
        } else if (argTypes[iarg] == +ATTRINT) {
            RALLOC(udps[0].arg[iarg].val, int, udps[0].arg[iarg].size);

            pint  = (int    *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                pint[ival] = argIdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREAL) {
            RALLOC(udps[0].arg[iarg].val, double, udps[0].arg[iarg].size);

            p1dbl  = (double *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p1dbl[ival] = argDdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREALSEN) {
            RALLOC(udps[0].arg[iarg].val, double, udps[0].arg[iarg].size);
            RALLOC(udps[0].arg[iarg].dot, double, udps[0].arg[iarg].size);

            p1dbl = (double *) (udps[0].arg[iarg].val);
            p2dbl = (double *) (udps[0].arg[iarg].dot);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p1dbl[ival] = argDdefs[iarg];
                p2dbl[ival] = 0;
            }
        } else if (argTypes[iarg] == +ATTRREBUILD) {
            /* do nothing */
        } else if (argTypes[iarg] == +ATTRRECYCLE) {
            /* do nothing */
        }
    }

    /* set up returns that describe the UDP */
    *nArgs    = NUMUDPARGS;
    *namex    = (char**)argNames;
    *typex    = argTypes;
    *idefault = argIdefs;
    *ddefault = argDdefs;

cleanup:
#ifdef DEBUG
    printf("exit  udpInitialize -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpNumBodys - number of Bodys expected in first arg to udpExecute  *
 *                                                                      *
 ************************************************************************
 */

int
udpNumBodys()                           /* (out) number of expected Bodys */
                                        /*       >0 number of Bodys */
                                        /*       <0 maximum number of Bodys */
{
    int    num=0;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpNumBodys()\n");
#endif

#ifdef NUMUDPINPUTBODYS
    num = NUMUDPINPUTBODYS;
#endif

//cleanup:
#ifdef DEBUG
    printf("exit  udpNumBodys -> num=%d\n", num);
#endif
    return num;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpBodyList - list of Bodys input to a UDF                         *
 *                                                                      *
 ************************************************************************
 */

int
udpBodyList(ego       ebody,            /* (in)  Body pointer */
            const int **bodyList,       /* (out) 0-terminated list of Bodys used by UDF */
            int       numudp,           /* (in)  number of instances */
            udp_T     udps[])           /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iudp;

    ROUTINE(udpBodyList);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpBodyList(ebody=%llx)\n", (long long)ebody);
#endif

    /* check that ebody matches one of the ebodys */
    iudp = status = findInstance(ebody, numudp, udps);
    CHECK_STATUS(findBody);

    /* return the bodyList */
    *bodyList = udps[iudp].bodyList;

cleanup:
#ifdef DEBUG
    printf("exit  udpBodyList -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpReset - reset the arguments to their defaults                   *
 *                                                                      *
 ************************************************************************
 */

int
udpReset(int   *NumUdp,                 /* (both)number of instances */
         udp_T udps[])                  /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iarg, ival, iudp;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpReset()\n");
#endif

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        udps[0].arg[iarg].size = 1;

        if        (argTypes[iarg] == +ATTRSTRING) {
            strcpy((char*)(udps[0].arg[iarg].val), "");
        } else if (argTypes[iarg] == +ATTRFILE) {
            strcpy((char*)(udps[0].arg[iarg].val), "");
        } else if (argTypes[iarg] == +ATTRINT) {
            int    *p  = (int    *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p[ival] = argIdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREAL) {
            double *p  = (double *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p[ival] = argDdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREALSEN) {
            double *p1 = (double *) (udps[0].arg[iarg].val);
            double *p2 = (double *) (udps[0].arg[iarg].dot);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p1[ival] = argDdefs[iarg];
                p2[ival] = 0;
            }
        } else if (argTypes[iarg] == +ATTRREBUILD) {
            /* do nothing */
        } else if (argTypes[iarg] == +ATTRRECYCLE) {
            /* do nothing */
        }
    }

    /* reset all the velocities */
    for (iudp = 0; iudp <= numUdp; iudp++) {
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == +ATTRREALSEN) {
                double *p2 = (double *) (udps[iudp].arg[iarg].dot);

                for (ival = 0; ival < udps[iudp].arg[iarg].size; ival++) {
                    p2[ival] = 0;
                }
            }
        }
    }

//cleanup:
#ifdef DEBUG
    printf("exit  udpReset -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpFree - free up all memory associated with UDP/UDF               *
 *                                                                      *
 ************************************************************************
 */

int
udpFree(int   numudp,                   /* (in)  number of instances */
        udp_T udps[])                   /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iarg, iudp;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpFree(numudp=%d, udps=%llx)\n", numudp, (long long)udps);
#endif

    if (udps == NULL) goto cleanup;

    for (iudp = 0; iudp <= numudp; iudp++) {

        /* arguments */
        for (iarg = 0; iarg < udps[iudp].narg; iarg++) {
            if (udps[iudp].arg[iarg].val != NULL) {
                EG_free(udps[iudp].arg[iarg].val);
            }

            if (udps[iudp].arg[iarg].dot != NULL) {
                EG_free(udps[iudp].arg[iarg].dot);
            }
        }
        EG_free(udps[iudp].arg);
        udps[iudp].arg = NULL;

        /* bodyList (but not 0 because it is same as another) */
        if (iudp > 0 && udps[iudp].bodyList != NULL) {
            EG_free(udps[iudp].bodyList);
            udps[iudp].bodyList = NULL;
        }

        /* private data */
        if (udps[iudp].data != NULL) {
#ifdef FREEUDPDATA
            FREEUDPDATA(udps[iudp].data);
#else
            EG_free(udps[iudp].data);
#endif
        }
        udps[iudp].data = NULL;
    }

    EG_free(udps);

cleanup:
#ifdef DEBUG
    printf("exit  udpFree -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSet - set an argument                                           *
 *                                                                      *
 ************************************************************************
 */

int
udpSet(char   name[],                   /* (in)  argument name */
       void   *value,                   /* (in)  pointer to values (can be char* or double*) */
       int    nrow,                     /* (in)  number of rows */
       int    ncol,                     /* (in)  number of columns */
       char   message[],                /* (out) error message (if any) */
       udp_T  udps[])                   /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    i, iarg, nvalue, ivalue;
    char   lowername[257];
    void   *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(udpSet);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpSet(name=%s, nrow=%d, ncol=%d)\n", name, nrow, ncol);
#endif

    if (message != NULL) {
        message[0] = '\0';
    }

    if        (name  == NULL) {
        status = EGADS_NONAME;
        goto cleanup;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        status = EGADS_NONAME;
        goto cleanup;
    } else if (value == NULL) {
        status = EGADS_NULLOBJ;
        goto cleanup;
    }

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    nvalue = nrow * MAX(ncol, 1);

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            udps[0].arg[iarg].nrow = nrow;
            udps[0].arg[iarg].ncol = ncol;
            udps[0].arg[iarg].size = nrow * MAX(ncol, 1);

            if        (argTypes[iarg] == +ATTRSTRING) {
                char   *p;
                char   *valueP = (char   *) value;

#ifdef DEBUG
                printf("   value=%s\n", valueP);
#endif

                RALLOC(udps[0].arg[iarg].val, char, nvalue+1);
                p = (char   *) (udps[0].arg[iarg].val);

                memcpy(p, valueP, (nvalue+1)*sizeof(char));
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRFILE) {
                char   *p;
                char   *valueP = (char   *) value;

#ifdef DEBUG
                printf("   value=%s\n", valueP);
#endif

                RALLOC(udps[0].arg[iarg].val, char, nvalue+1);
                p = (char   *) (udps[0].arg[iarg].val);

                memcpy(p, valueP, (nvalue+1)*sizeof(char));
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRINT) {
                int    *p;
                double *valueP = (double *) value;

#ifdef DEBUG
                printf("   value=%f", valueP[0]);
                for (ivalue = 1; ivalue < nvalue; ivalue++) {
                    printf(" %f", valueP[ivalue]);
                }
                printf("\n");
#endif

                RALLOC(udps[0].arg[iarg].val, int, nvalue);
                p = (int    *) (udps[0].arg[iarg].val);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p[ivalue] = NINT(valueP[ivalue]);
                }
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRREAL) {
                double *p;
                double *valueP = (double *) value;

#ifdef DEBUG
                printf("   value=%f", valueP[0]);
                for (ivalue = 1; ivalue < nvalue; ivalue++) {
                    printf(" %f", valueP[ivalue]);
                }
                printf("\n");
#endif

                RALLOC(udps[0].arg[iarg].val, double, nvalue);
                p = (double *) (udps[0].arg[iarg].val);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p[ivalue] = valueP[ivalue];
                }
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRREALSEN) {
                double *p1, *p2;
                double *valueP = (double *) value;

#ifdef DEBUG
                printf("   value=%f", valueP[0]);
                for (ivalue = 1; ivalue < nvalue; ivalue++) {
                    printf(" %f", valueP[ivalue]);
                }
                printf("\n");
#endif

                RALLOC(udps[0].arg[iarg].val, double, nvalue);
                RALLOC(udps[0].arg[iarg].dot, double, nvalue);
                p1 = (double *) (udps[0].arg[iarg].val);
                p2 = (double *) (udps[0].arg[iarg].dot);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p1[ivalue] = valueP[ivalue];
                    p2[ivalue] = 0;
                }
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRREBUILD) {
                /* do nothing */
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRRECYCLE) {
                /* do nothing */
                goto cleanup;
            }
        }
    }

    if (message != NULL) {
        snprintf(message, 256, "Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == +ATTRSTRING  ||
                argTypes[iarg] == +ATTRFILE    ||
                argTypes[iarg] == +ATTRINT     ||
                argTypes[iarg] == +ATTRREAL    ||
                argTypes[iarg] == +ATTRREALSEN ||
                argTypes[iarg] == +ATTRREBUILD ||
                argTypes[iarg] == +ATTRRECYCLE   ) {
                strncat(message, " ",            256);
                strncat(message, argNames[iarg], 256);
//$$$                snprintf(message, 256, "%s %s", message, argNames[iarg]);
            }
        }
    }

    status = EGADS_INDEXERR;

cleanup:
#ifdef DEBUG
    printf("exit  udpSet -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpGet - return an output parameter                                *
 *                                                                      *
 ************************************************************************
 */

int
udpGet(ego    ebody,                    /* (in)  Body pointer */
       char   name[],                   /* (in)  argument name */
       int    *nrow,                    /* (out) number of rows */
       int    *ncol,                    /* (out) number of columns */
       void   **val,                    /* (out) argument values (can be int* or double*) */
       void   **dot,                    /* (out) argument velocities (double*) */
       char   message[],                /* (out) error message */
       int    numudp,                   /* (in)  number of instances */
       udp_T  udps[])                   /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    i, iudp, iarg;
    char   lowername[257];
    void   *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(udpGet);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpGet(ebody=%llx, name=%s)\n", (long long)ebody, name);
#endif

    *nrow = 0;
    *ncol = 0;

    if (message != NULL) {
        message[0] = '\0';
    }

    if (name == NULL) {
        status = EGADS_NONAME;
        goto cleanup;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        status =  EGADS_NONAME;
        goto cleanup;
    }

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    /* check that ebody matches one of the ebodys */
    iudp = status = findInstance(ebody, numudp, udps);
    CHECK_STATUS(findBody);

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            if        (argTypes[iarg] == -ATTRINT) {
                int *ival, *valP=NULL;

                *nrow = udps[iudp].arg[iarg].nrow;
                *ncol = udps[iudp].arg[iarg].ncol;

                RALLOC(valP, int, (*nrow)*(*ncol));

                ival = (int*)(udps[iudp].arg[iarg].val);

                for (i = 0; i < (*nrow)*(*ncol); i++) {
                    valP[i] = ival[i];
#ifdef DEBUG
                    printf("   ival[%d]=%d\n", i, valP[i]);
#endif
                }

                *val = valP;
                *dot = NULL;
                goto cleanup;
            } else if (argTypes[iarg] == -ATTRREAL) {
                double *dval, *valP=NULL;

                *nrow = udps[iudp].arg[iarg].nrow;
                *ncol = udps[iudp].arg[iarg].ncol;

                RALLOC(valP, double, (*nrow)*(*ncol));

                dval = (double*)(udps[iudp].arg[iarg].val);

                for (i = 0; i < (*nrow)*(*ncol); i++) {
                    valP[i] = dval[i];
#ifdef DEBUG
                    printf("   dval[%d]=%f\n", i, valP[i]);
#endif
                }

                *val = valP;
                *dot = NULL;
                goto cleanup;
            } else if (argTypes[iarg] == -ATTRREALSEN) {
                double *dval, *ddot, *valP=NULL, *dotP=NULL;

                *nrow = udps[iudp].arg[iarg].nrow;
                *ncol = udps[iudp].arg[iarg].ncol;

                RALLOC(valP, double, (*nrow)*(*ncol));
                RALLOC(dotP, double, (*nrow)*(*ncol));

                dval = (double*)(udps[iudp].arg[iarg].val);
                ddot = (double*)(udps[iudp].arg[iarg].dot);

                for (i = 0; i < (*nrow)*(*ncol); i++) {
                    valP[i] = dval[i];
                    dotP[i] = ddot[i];
#ifdef DEBUG
                    printf("   dval[%d]=%f, ddot[%d]=%f\n", i, valP[i], i, dotP[i]);
#endif
                }

                *val = valP;
                *dot = dotP;
                goto cleanup;
            }
        }
    }

    if (message != NULL) {
        snprintf(message, 256, "Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == -ATTRINT    ||
                argTypes[iarg] == -ATTRREAL   ||
                argTypes[iarg] == -ATTRREALSEN  ) {
                strncat(message, " ",            256);
                strncat(message, argNames[iarg], 256);
//$$$                snprintf(message, 256, "%s %s", message, argNames[iarg]);
            }
        }
    }

    status = EGADS_INDEXERR;

cleanup:
#ifdef DEBUG
    printf("exit  udpGet -> status=%d, nrow=%d, ncol=%d\n", status, *nrow, *ncol);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpVel - set velocity of an argument                               *
 *                                                                      *
 ************************************************************************
 */

int
udpVel(ego    ebody,                    /* (in)  Body pointer */
       char   name[],                   /* (in)  argument name */
       double dot[],                    /* (in)  argument velocity(s) */
       int    ndot,                     /* (in)  number of velocities */
       int    numudp,                   /* (in)  number of instances */
       udp_T  udps[])                   /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    i, iudp, iarg, idot, hasdots;
    char   lowername[257];

    ROUTINE(udpVel);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpVel(ebody=%llx, name=%s, dot=%f", (long long)ebody, name, dot[0]);
    for (i = 1; i < ndot; i++) {
        printf(" %f", dot[i]);
    }
    printf(")\n");
#endif

    if        (name  == NULL) {
        status = EGADS_NONAME;
        goto cleanup;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        status = EGADS_NONAME;
        goto cleanup;
    } else if (dot == NULL) {
        status = EGADS_NULLOBJ;
        goto cleanup;
    }

    /* check that ebody matches one of the ebodys */
    iudp = status = findInstance(ebody, numudp, udps);
    if (status == EGADS_NOTMODEL) goto cleanup;
    CHECK_STATUS(findBody);

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            if (argTypes[iarg] == +ATTRREALSEN) {
                double *p = (double *) (udps[iudp].arg[iarg].dot);

                if (ndot != udps[iudp].arg[iarg].size) {
                    printf("size of value and dot do not agree for name=%s, ndot=%d, udps[%d].arg[%d].size=%d\n",
                           lowername, ndot, iudp, iarg, udps[iudp].arg[iarg].size);
                    status = EGADS_RANGERR;
                    goto cleanup;
                }

                for (idot = 0; idot < ndot; idot++) {
                    if (fabs(p[idot]-dot[idot]) > EPS12) {
                        udps[iudp].ndotchg++;
                        }

                    p[idot] = dot[idot];
                }
                status = EGADS_SUCCESS;
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRREBUILD) {
                /* do nothing */
                status = EGADS_SUCCESS;
                goto cleanup;
            } else if (argTypes[iarg] == +ATTRRECYCLE) {
                /* do nothing */
                status = EGADS_SUCCESS;
                goto cleanup;
            } else {
                hasdots = 0;
                for (idot = 0; idot < ndot; idot++) {
                    if (fabs(dot[idot]) > 1.0e-6) {
                        hasdots++;
                    }
                }
                if (hasdots == 0) {
                    status = EGADS_SUCCESS;
                    goto cleanup;
                } else {
                    status = EGADS_INDEXERR;
                    goto cleanup;
                }
            }
        }
    }

    status = EGADS_INDEXERR;

cleanup:
#ifdef DEBUG
    printf("exit  udpVel -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpPost - reset ndotchg flag                                       *
 *                                                                      *
 ************************************************************************
 */

int
udpPost(ego    ebody,                   /* (in)  Body pointer */
        int    numudp,                  /* (in)  number of instances */
        udp_T  udps[])                  /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iudp;

    ROUTINE(udpPost);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpPost(ebody=%llx)\n", (long long)ebody);
#endif

    /* check that ebody matches one of the ebodys */
    iudp = status = findInstance(ebody, numudp, udps);
    if (status == EGADS_NOTMODEL) goto cleanup;
    CHECK_STATUS(findBody);

    udps[iudp].ndotchg = 0;

cleanup:
#ifdef DEBUG
    printf("exit  udpPost -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpClean - clean the udp cache                                     *
 *                                                                      *
 ************************************************************************
 */

int
udpClean(int   *NumUdp,                 /* (both)number of instances */
         udp_T udps[])                  /* (in)  array  of instances */
{
    int    status = EGADS_SUCCESS;

    int    iarg;

    ROUTINE(udpClean);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpClean()\n");
#endif

    /* decrement number of udps based upon NULL ebody entries */
    while (numUdp > 0) {
        if (udps[numUdp].ebody == NULL) {

            /* arguments */
            for (iarg = 0; iarg < udps[numUdp].narg; iarg++) {
                if (udps[numUdp].arg[iarg].val != NULL) {
                    EG_free(udps[numUdp].arg[iarg].val);
                }

                if (udps[numUdp].arg[iarg].dot != NULL) {
                    EG_free(udps[numUdp].arg[iarg].dot);
                }
            }
            EG_free(udps[numUdp].arg);

            /* bodyList */
            if (udps[numUdp].bodyList != NULL) {
                EG_free(udps[numUdp].bodyList);
                udps[numUdp].bodyList = NULL;
            }

            /* private data */
            if (udps[numUdp].data != NULL) {
#ifdef FREEUDPDATA
                FREEUDPDATA(udps[numUdp].data);
#else
                EG_free(udps[numUdp].data);
#endif
            }
            udps[numUdp].data = NULL;

            numUdp--;
        } else {
            break;
        }
    }

//cleanup:
#ifdef DEBUG
    printf("exit  udpClean -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpMesh - return mesh associated with primitive                    *
 *                                                                      *
 ************************************************************************
 */

int
udpMesh(ego    ebody,                   /* (in)  Body pointer for mesh */
/*@unused@*/int    iMesh,               /* (in)  mesh index */
        int    *imax,                   /* (out) i-dimension of mesh */
        int    *jmax,                   /* (out) j-dimension of mesh */
        int    *kmax,                   /* (out) k-dimension of mesh */
        double *mesh[])                 /* (out) array of mesh points */
{
    int    status = EGADS_SUCCESS;

    ROUTINE(udpMesh);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter udpMesh(ebody=%llx)\n", (long long)ebody);
#endif

    *imax = 0;
    *jmax = 0;
    *kmax = 0;
    *mesh = NULL;

    /* check that ebody matches one of the ebodys */
    status = findInstance(ebody, numUdp, udps);
    CHECK_STATUS(findBody);

cleanup:
#ifdef DEBUG
    printf("exit  udpMesg -> status=%d, imax=%d, jmax=%d, kmax=%d\n", status, *imax, *jmax, *kmax);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   cacheUdp - cache arguments of current instance                     *
 *                                                                      *
 ************************************************************************
 */

static int
CacheUdp(/*@null@*/ego emodel,          /* (in)  Model with __bodyList__ */
         int   *NumUdp,                 /* (both)number of instances */
         udp_T *Udps[])                 /* (both)array  of instances */
{
    int     status = SUCCESS;

    int     iarg, isize, attrType, attrLen, i;
    CINT    *tempIlist;
    CDOUBLE *tempRlist;
    CCHAR   *tempClist;
    udp_T   *udps;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */

    ROUTINE(cacheUdp);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter cacheUdp(emodel=%llx)\n", (long long)emodel);
#endif

    udps = *Udps;

    /* create the BodyList (0-terminated) */
    if (emodel != NULL) {
        status = EG_attributeRet(emodel, "__bodyList__", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS && attrType == ATTRINT) {
            FREE(  udps[0].bodyList);
            MALLOC(udps[0].bodyList, int, attrLen+1);

            for (i = 0; i < attrLen; i++) {
                udps[0].bodyList[i] = tempIlist[i];
            }
            udps[0].bodyList[attrLen] = 0;
        }
    }

    /* increment number of UDPs in the cache */
    numUdp++;

#ifdef DEBUG
    printf("copying from udps[0] to udps[%d]\n", numUdp);
#endif

    /* increase the arrays to make room for the new UDP */
    RALLOC(*Udps, udp_T, numUdp+1);
    udps = *Udps;

    udps[numUdp].ebody = NULL;

    /* copy info from udps[0] into udps[numUdp] */
    udps[numUdp].ndotchg  = udps[0].ndotchg;
    udps[numUdp].bodyList = udps[0].bodyList;

#ifdef COPYUDPDATA
    COPYUDPDATA(udps[0].data, &(udps[numUdp].data));
#else
    udps[numUdp].data = NULL;
#endif

    udps[numUdp].narg = NUMUDPARGS;
    udps[numUdp].arg  = NULL;
    RALLOC(udps[numUdp].arg, udparg_T, NUMUDPARGS);
    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        udps[numUdp].arg[iarg].type = argTypes[iarg];
        udps[numUdp].arg[iarg].nrow = udps[0].arg[iarg].nrow;
        udps[numUdp].arg[iarg].ncol = udps[0].arg[iarg].ncol;
        udps[numUdp].arg[iarg].size = udps[0].arg[iarg].size;
        udps[numUdp].arg[iarg].val  = NULL;
        udps[numUdp].arg[iarg].dot  = NULL;
        isize = udps[numUdp].arg[iarg].size;

        if        (argTypes[iarg] == +ATTRSTRING) {
            MALLOC(udps[numUdp].arg[iarg].val, char, isize+1);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, (isize+1)*sizeof(char  ));
        } else if (argTypes[iarg] == +ATTRFILE) {
            MALLOC(udps[numUdp].arg[iarg].val, char, isize+1);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, (isize+1)*sizeof(char  ));
        } else if (argTypes[iarg] == +ATTRINT ||
                   argTypes[iarg] == -ATTRINT   ) {
            MALLOC(udps[numUdp].arg[iarg].val, int, isize);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(int   ));
        } else if (argTypes[iarg] == +ATTRREAL ||
                   argTypes[iarg] == -ATTRREAL   ) {
            MALLOC(udps[numUdp].arg[iarg].val, double, isize);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));
        } else if (argTypes[iarg] == +ATTRREALSEN ||
                   argTypes[iarg] == -ATTRREALSEN   ) {
            MALLOC(udps[numUdp].arg[iarg].val, double, isize);
            MALLOC(udps[numUdp].arg[iarg].dot, double, isize);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));
            memcpy(udps[numUdp].arg[iarg].dot, udps[0].arg[iarg].dot, isize*sizeof(double));
        } else if (argTypes[iarg] == +ATTRREBUILD) {
        } else if (argTypes[iarg] == +ATTRRECYCLE) {
        } else if (argTypes[iarg] == 0) {
            MALLOC(udps[numUdp].arg[iarg].val, double, 1);

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));
        } else {
            printf("bad argType[%d]=%d in CacheUdp\n", iarg, argTypes[iarg]);
            exit(0);
        }
    }

cleanup:
#ifdef DEBUG
    printf("exit  cacheUdp -> status=%d\n", status);
#endif
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   findInstance - find UDP instance that matches ebody                *
 *                                                                      *
 ************************************************************************
 */

static int
findInstance(ego    ebody,              /* (in)  Body to match */
             int    numudp,             /* (in)  number of instances */
             udp_T  udps[])             /* (in)  array  of instances */
{
    int    iudp, judp, status, oclass, mtype;

    ROUTINE(finBody);

    /* --------------------------------------------------------------- */

    iudp = EGADS_NOTMODEL;

    for (judp = 1; judp <= numudp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            goto cleanup;
        }
    }

    /* udps[].ebody will be a MODEL is the UDP/UDF returned more than
       one Body.  in this case, we should return iudp=1 so that checks
       in udpSensitivity and udp_postSens do not return an error */
    if (numudp > 0) {
        status = EG_getInfo(udps[1].ebody, &oclass, &mtype, NULL, NULL, NULL);
        if (status < EGADS_SUCCESS) {
            return status;
        } else if (oclass == MODEL) {
            iudp = 1;
        }
    }

#ifdef DEBUG
    printf("in findInstance(ebody=%llx)\n", (long long)ebody);
    for (judp = 1; judp <= numudp; judp++) {
        printf("    udps[%2d]=%llx\n", judp, (long long)(udps[judp].ebody));
    }
#endif

cleanup:
    return iudp;
}
