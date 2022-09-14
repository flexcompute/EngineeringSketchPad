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

#include "OpenCSM.h"

#define STRLEN(A)   (int)strlen(A)


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
              double *ddefault[])       /* (out) array  of argument double  defaults */
{
    int    iarg, ival, *pint;
    double *p1dbl, *p2dbl;
    char   *pchar;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpInitialize()\n");
#endif

    /* make the initial array that holds the udps */
    if (udps == NULL) {

        udps = (udp_T *) EG_alloc(sizeof(udp_T));
        if (udps == NULL) return EGADS_MALLOC;

        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            udps[0].arg[iarg].size = 1;

            if        (argTypes[iarg] == +ATTRSTRING) {
                udps[0].arg[iarg].val = (char   *) EG_alloc(MAX_EXPR_LEN*sizeof(char  ));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == +ATTRFILE) {
                udps[0].arg[iarg].val = (char   *) EG_alloc(MAX_EXPR_LEN*sizeof(char  ));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == +ATTRINT) {
                udps[0].arg[iarg].val = (int    *) EG_alloc(  sizeof(int   ));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == -ATTRINT) {
                udps[0].arg[iarg].val = (int    *) EG_alloc(  sizeof(int   ));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == +ATTRREAL) {
                udps[0].arg[iarg].val = (double *) EG_alloc(  sizeof(double));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == -ATTRREAL) {
                udps[0].arg[iarg].val = (double *) EG_alloc(  sizeof(double));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == +ATTRREALSEN) {
                udps[0].arg[iarg].val = (double *) EG_alloc(  sizeof(double));
                udps[0].arg[iarg].dot = (double *) EG_alloc(  sizeof(double));
                if (udps[0].arg[iarg].val == NULL || udps[0].arg[iarg].dot == NULL) return EGADS_MALLOC;
            } else if (argTypes[iarg] == 0) {
                udps[0].arg[iarg].val = (double *) EG_alloc(  sizeof(double));
                udps[0].arg[iarg].dot = NULL;
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
            } else {
                printf("bad argType[%d]=%d in udpInitialize\n", iarg, argTypes[iarg]);
                exit(0);
            }
        }
    }

    /* initialize the array elements that will hold the "current" settings */
    udps[0].ebody    = NULL;
    udps[0].bodyList = NULL;
    udps[0].data     = NULL;

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if        (argTypes[iarg] == +ATTRSTRING) {
            pchar  = (char   *) (udps[0].arg[iarg].val);

            strcpy(pchar, "");
        } else if (argTypes[iarg] == +ATTRFILE) {
            pchar  = (char   *) (udps[0].arg[iarg].val);

            strcpy(pchar, "");
        } else if (argTypes[iarg] == +ATTRINT) {
            udps[0].arg[iarg].val = (int *) EG_reall(udps[0].arg[iarg].val, udps[0].arg[iarg].size*sizeof(int));
            if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;

            pint  = (int    *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                pint[ival] = argIdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREAL) {
            udps[0].arg[iarg].val = (double *) EG_reall(udps[0].arg[iarg].val, udps[0].arg[iarg].size*sizeof(double));
            if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;

            p1dbl  = (double *) (udps[0].arg[iarg].val);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p1dbl[ival] = argDdefs[iarg];
            }
        } else if (argTypes[iarg] == +ATTRREALSEN) {
            udps[0].arg[iarg].val = (double *) EG_reall(udps[0].arg[iarg].val, udps[0].arg[iarg].size*sizeof(double));
            udps[0].arg[iarg].dot = (double *) EG_reall(udps[0].arg[iarg].dot, udps[0].arg[iarg].size*sizeof(double));
            if (udps[0].arg[iarg].val == NULL || udps[0].arg[iarg].dot == NULL) return EGADS_MALLOC;

            p1dbl = (double *) (udps[0].arg[iarg].val);
            p2dbl = (double *) (udps[0].arg[iarg].dot);

            for (ival = 0; ival < udps[0].arg[iarg].size; ival++) {
                p1dbl[ival] = argDdefs[iarg];
                p2dbl[ival] = 0;
            }
        }
    }

    /* set up returns that describe the UDP */
    *nArgs    = NUMUDPARGS;
    *namex    = (char**)argNames;
    *typex    = argTypes;
    *idefault = argIdefs;
    *ddefault = argDdefs;

    return EGADS_SUCCESS;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpNumBodys - number of Bodys expected in first arg to udpExecute  *
 *                                                                      *
 ************************************************************************
 */

int
udpNumBodys()
{

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpNumBodys()\n");
#endif

#ifdef NUMUDPINPUTBODYS
    return NUMUDPINPUTBODYS;
#else
    return 0;
#endif
}



/*
 ************************************************************************
 *                                                                      *
 *   udpBodyList - list of Bodys input to a UDF                         *
 *                                                                      *
 ************************************************************************
 */

int
udpBodyList(ego    ebody,               /* (in)  Body pointer */
            const int **bodyList)       /* (out) 0-terminated list of Bodys used by UDF */
{
    int     iudp, judp;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpBodyList(ebody=%llx)\n", (long long)ebody);
#endif

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

    /* return the bodyList */
    *bodyList = udps[iudp].bodyList;

    return EGADS_SUCCESS;
}



/*
 ************************************************************************
 *                                                                      *
 *   udpReset - reset the arguments to their defaults                   *
 *                                                                      *
 ************************************************************************
 */

int
udpReset(int flag)                      /* (in)  flag: 0=reset current, 1=close up */
{
    int  iarg, ival, iudp;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpReset(flag=%d)\n", flag);
#endif

    /* reset the "current" settings */
    if (flag == 0) {
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

    /* called when closing up */
    } else if (udps != NULL) {
        for (iudp = 0; iudp <= numUdp; iudp++) {

            /* arguments */
            for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
                if (udps[iudp].arg[iarg].val != NULL) {
                    EG_free(udps[iudp].arg[iarg].val);
                }

                if (udps[iudp].arg[iarg].dot != NULL) {
                    EG_free(udps[iudp].arg[iarg].dot);
                }
            }

            /* bodyList (but not 0 because it is same as another) */
            if (iudp > 0 && udps[iudp].bodyList != NULL) {
                EG_free(udps[iudp].bodyList);
                udps[iudp].bodyList = NULL;
            }

            /* private data (but not 0 because it is same as another) */
            if (iudp > 0 && udps[iudp].data != NULL) {
#ifdef FREEUDPDATA
                FREEUDPDATA(udps[iudp].data);
#else
                EG_free(udps[iudp].data);
#endif
                udps[iudp].data = NULL;
            }
        }

        EG_free(udps);
        udps   = NULL;
        numUdp = 0;
    }

    return EGADS_SUCCESS;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSet - set an argument                                           *
 *                                                                      *
 ************************************************************************
 */

int
udpSet(char name[],                     /* (in)  argument name */
       void *value,                     /* (in)  pointer to values (can be char* or double*) */
       int  nvalue,                     /* (in)  number  of values */
       char message[])                  /* (out) error message (if any) */
{
    int  i, iarg, ivalue;
    char lowername[257];

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSet(name=%s, nvalue=%d)\n", name, nvalue);
#endif

    if (message != NULL) {
        message[0] = '\0';
    }

    if        (name  == NULL) {
        return EGADS_NONAME;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        return EGADS_NONAME;
    } else if (value == NULL) {
        return EGADS_NULLOBJ;
    }

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            udps[0].arg[iarg].size = nvalue;

            if        (argTypes[iarg] == +ATTRSTRING) {
                char   *p;
                char   *valueP = (char   *) value;

#ifdef DEBUG
                printf("   value=%s\n", valueP);
#endif

                udps[0].arg[iarg].val = (char *) EG_reall(udps[0].arg[iarg].val, (nvalue+1)*sizeof(char  ));
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
                p = (char   *) (udps[0].arg[iarg].val);

                memcpy(p, valueP, (nvalue+1)*sizeof(char));
                return EGADS_SUCCESS;
            } else if (argTypes[iarg] == +ATTRFILE) {
                char   *p;
                char   *valueP = (char   *) value;

#ifdef DEBUG
                printf("   value=%s\n", valueP);
#endif

                udps[0].arg[iarg].val = (char *) EG_reall(udps[0].arg[iarg].val, (nvalue+1)*sizeof(char  ));
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
                p = (char   *) (udps[0].arg[iarg].val);

                memcpy(p, valueP, (nvalue+1)*sizeof(char));
                return EGADS_SUCCESS;
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

                udps[0].arg[iarg].val = (int   *) EG_reall(udps[0].arg[iarg].val, nvalue*sizeof(int   ));
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
                p = (int    *) (udps[0].arg[iarg].val);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p[ivalue] = NINT(valueP[ivalue]);
                }
                return EGADS_SUCCESS;
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

                udps[0].arg[iarg].val = (double *) EG_reall(udps[0].arg[iarg].val, nvalue*sizeof(double));
                if (udps[0].arg[iarg].val == NULL) return EGADS_MALLOC;
                p = (double *) (udps[0].arg[iarg].val);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p[ivalue] = valueP[ivalue];
                }
                return EGADS_SUCCESS;
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

                udps[0].arg[iarg].val = (double *) EG_reall(udps[0].arg[iarg].val, nvalue*sizeof(double));
                udps[0].arg[iarg].dot = (double *) EG_reall(udps[0].arg[iarg].dot, nvalue*sizeof(double));
                if (udps[0].arg[iarg].val == NULL || udps[0].arg[iarg].dot == NULL) return EGADS_MALLOC;
                p1 = (double *) (udps[0].arg[iarg].val);
                p2 = (double *) (udps[0].arg[iarg].dot);

                for (ivalue = 0; ivalue < nvalue; ivalue++) {
                    p1[ivalue] = valueP[ivalue];
                    p2[ivalue] = 0;
                }
                return EGADS_SUCCESS;
            }
        }
    }

    if (message != NULL) {
        snprintf(message, 256, "Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == +ATTRSTRING ||
                argTypes[iarg] == +ATTRFILE   ||
                argTypes[iarg] == +ATTRINT    ||
                argTypes[iarg] == +ATTRREAL   ||
                argTypes[iarg] == +ATTRREALSEN  ) {
                strncat(message, " ",            256);
                strncat(message, argNames[iarg], 256);
//$$$                snprintf(message, 256, "%s %s", message, argNames[iarg]);
            }
        }

        printf(" udpSet: %s\n", message);
    } else {
        printf(" udpSet: Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == +ATTRSTRING ||
                argTypes[iarg] == +ATTRFILE   ||
                argTypes[iarg] == +ATTRINT    ||
                argTypes[iarg] == +ATTRREAL   ||
                argTypes[iarg] == +ATTRREALSEN  ) {
                printf(" %s", argNames[iarg]);
            }
        }
        printf("\n");
    }

    return EGADS_INDEXERR;
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
       void   *value,                   /* (out) argument value (can be int* or double*) */
       char   message[])
{
    int  i, iudp, judp, iarg;
    char lowername[257];

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpGet(ebody=%llx, name=%s)\n", (long long)ebody, name);
#endif

    if (message != NULL) {
        message[0] = '\0';
    }

    if (name == NULL) {
        return EGADS_NONAME;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        return EGADS_NONAME;
    }

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
#ifdef DEBUG
        printf("udpGet: udps[%d].ebody=%lx\n", judp, (long)(udps[judp].ebody));
#endif
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            if        (argTypes[iarg] == -ATTRINT) {
                int    *p      = (int    *) (udps[0].arg[iarg].val);
                int    *valueP = (int    *) value;

                *valueP = p[0];
                return EGADS_SUCCESS;
            } else if (argTypes[iarg] == -ATTRREAL) {
                double *p      = (double *) (udps[0].arg[iarg].val);
                double *valueP = (double *) value;

                *valueP = p[0];
                return EGADS_SUCCESS;
            }
        }
    }

    if (message != NULL) {
        snprintf(message, 256, "Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == -ATTRINT ||
                argTypes[iarg] == -ATTRREAL  ) {
                strncat(message, " ",            256);
                strncat(message, argNames[iarg], 256);
//$$$                snprintf(message, 256, "%s %s", message, argNames[iarg]);
            }
        }

        printf(" udpGet: %s\n", message);
    } else {
        printf(" udpGet: Parameter \"%s\" not known.  should be one of:", name);
        for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
            if (argTypes[iarg] == -ATTRINT ||
                argTypes[iarg] == -ATTRREAL  ) {
                printf(" %s", argNames[iarg]);
            }
        }
        printf("\n");
    }

    return EGADS_INDEXERR;
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
       int    ndot)                     /* (in)  number of velocities */
{
    int  i, iudp, judp, iarg, idot, hasdots;
    char lowername[257];

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpVel(ebody=%llx, name=%s, dot=%f", (long long)ebody, name, dot[0]);
    for (i = 1; i < ndot; i++) {
        printf(" %f", dot[i]);
    }
    printf(")\n");
#endif

    if        (name  == NULL) {
        return EGADS_NONAME;
    } else if (STRLEN(name) == 0 || STRLEN(name) > 255) {
        return EGADS_NONAME;
    } else if (dot == NULL) {
        return EGADS_NULLOBJ;
    }

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
#ifdef DEBUG
        printf("udpVel: udps[%d].ebody=%lx\n", judp, (long)(udps[judp].ebody));
#endif
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    for (i = 0; i <= STRLEN(name); i++) {
        lowername[i] = tolower(name[i]);
    }

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        if (strcmp(lowername, argNames[iarg]) == 0) {
            if (argTypes[iarg] == +ATTRREALSEN) {
                double *p = (double *) (udps[iudp].arg[iarg].dot);

                if (ndot != udps[iudp].arg[iarg].size) {
                    printf("size of value and dot do not agree\n");
                    return EGADS_RANGERR;
                }

                for (idot = 0; idot < ndot; idot++) {
                    p[idot] = dot[idot];
                }
                return EGADS_SUCCESS;
//$$$            } else if (argTypes[iarg] == +ATTRREAL) {
//$$$                double *p = (double *) (udps[iudp].arg[iarg].dot);
//$$$
//$$$                if (ndot != udps[iudp].arg[iarg].size) {
//$$$                    printf("size of value and dot do not agree\n");
//$$$                    return EGADS_RANGERR;
//$$$                }
//$$$
//$$$                for (idot = 0; idot < ndot; idot++) {
//$$$                    p[idot] = dot[idot];
//$$$                }
//$$$                return EGADS_NOTFOUND;
            } else {
                hasdots = 0;
                for (idot = 0; idot < ndot; idot++) {
                    if (fabs(dot[idot]) > 1.0e-6) {
                        hasdots++;
                    }
                }
                if (hasdots == 0) {
                    return EGADS_SUCCESS;
                } else {
                    return EGADS_INDEXERR;
                }
            }
        }
    }

    return EGADS_INDEXERR;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpClean - clean the udp cache                                     *
 *                                                                      *
 ************************************************************************
 */

int
udpClean(ego ebody)                     /* (in)   Body pointer to clean from cache */
{
    int iudp, judp, iarg;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpClean(ebody=%llx)\n", (long long)ebody);
#endif

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
#ifdef DEBUG
        printf("udpClean: udps[%d].ebody=%lx\n", judp, (long)(udps[judp].ebody));
#endif
        if (ebody == udps[judp].ebody) {
            iudp = judp;

            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    /* remove ebody from the udp's cache */
    udps[iudp].ebody = NULL;

    /* remove any private data from udp's cache */
    if (udps[iudp].data != NULL) {
#ifdef FREEUDPDATA
        FREEUDPDATA(udps[iudp].data);
#else
        EG_free(udps[iudp].data);
#endif
        udps[iudp].data = NULL;
    }

    /* decrement number of udps based upon NULL ebody entries */
    while (numUdp > 0) {
        if (udps[numUdp].ebody == NULL) {

            /* arguments */
            for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
                if (udps[numUdp].arg[iarg].val != NULL) {
                    EG_free(udps[numUdp].arg[iarg].val);
                }

                if (udps[numUdp].arg[iarg].dot != NULL) {
                    EG_free(udps[numUdp].arg[iarg].dot);
                }
            }

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
                udps[numUdp].data = NULL;
            }

            numUdp--;
        } else {
            break;
        }
    }

    return EGADS_SUCCESS;
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
    int    iudp, judp;

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpMesh(ebody=%llx)\n", (long long)ebody);
#endif

    *imax = 0;
    *jmax = 0;
    *kmax = 0;
    *mesh = NULL;

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
#ifdef DEBUG
        printf("udpMesh: udps[%d].ebody=%lx\n", judp, (long)(udps[judp].ebody));
#endif
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    return EGADS_SUCCESS;
}


/*
 ************************************************************************
 *                                                                      *
 *   cacheUdp - cache arguments of current instance                     *
 *                                                                      *
 ************************************************************************
 */

static int
cacheUdp(/*@null@*/ego emodel)          /* (in)  Model with __bodyList__ */
{
    int     status = SUCCESS;
    
    int     iarg, isize, attrType, attrLen, i;
    CINT    *tempIlist;
    CDOUBLE *tempRlist;
    CCHAR   *tempClist;

    ROUTINE(cacheUdp);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("cacheUdp(emodel=%lx)\n", (long)emodel);
#endif

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

    /* increase the arrays to make room for the new UDP */
    udps = (udp_T *) EG_reall(udps, (numUdp+1)*sizeof(udp_T));
    if (udps == NULL) return EGADS_MALLOC;
    udps[numUdp].ebody = NULL;

    /* copy info from udps[0] into udps[numUdp] */
    udps[numUdp].bodyList = udps[0].bodyList;
    udps[numUdp].data     = udps[0].data;

    for (iarg = 0; iarg < NUMUDPARGS; iarg++) {
        isize = udps[numUdp].arg[iarg].size = udps[0].arg[iarg].size;

        if        (argTypes[iarg] == +ATTRSTRING) {
            udps[numUdp].arg[iarg].val = (char *) EG_alloc((isize+1)*sizeof(char  ));
            udps[numUdp].arg[iarg].dot = NULL;
            if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, (isize+1)*sizeof(char  ));
        } else if (argTypes[iarg] == +ATTRFILE) {
            udps[numUdp].arg[iarg].val = (char *) EG_alloc((isize+1)*sizeof(char  ));
            udps[numUdp].arg[iarg].dot = NULL;
            if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

            memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, (isize+1)*sizeof(char  ));
        } else if (argTypes[iarg] == +ATTRINT ||
                   argTypes[iarg] == -ATTRINT   ) {
            if (isize > 0) {
                udps[numUdp].arg[iarg].val = (int *) EG_alloc(isize*sizeof(int   ));
                udps[numUdp].arg[iarg].dot = NULL;
                if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

                memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(int   ));
            } else {
                udps[numUdp].arg[iarg].val = NULL;
                udps[numUdp].arg[iarg].dot = NULL;
            }
        } else if (argTypes[iarg] == +ATTRREAL ||
                   argTypes[iarg] == -ATTRREAL   ) {
            if (isize > 0) {
                udps[numUdp].arg[iarg].val = (double *) EG_alloc(isize*sizeof(double));
                udps[numUdp].arg[iarg].dot = NULL;
                if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

                memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));
            } else {
                udps[numUdp].arg[iarg].val = NULL;
                udps[numUdp].arg[iarg].dot = NULL;
            }
        } else if (argTypes[iarg] == +ATTRREALSEN) {
            if (isize > 0) {
                udps[numUdp].arg[iarg].val = (double *) EG_alloc(isize*sizeof(double));
                if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

                memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));

                udps[numUdp].arg[iarg].dot = (double *) EG_alloc(isize*sizeof(double));
                if (udps[numUdp].arg[iarg].dot == NULL) return EGADS_MALLOC;

                memcpy(udps[numUdp].arg[iarg].dot, udps[0].arg[iarg].dot, isize*sizeof(double));
            } else {
                udps[numUdp].arg[iarg].val = NULL;
                udps[numUdp].arg[iarg].dot = NULL;
            }
        } else if (argTypes[iarg] == 0) {
            if (isize > 0) {
                udps[numUdp].arg[iarg].val = (double *) EG_alloc(isize*sizeof(double));
                udps[numUdp].arg[iarg].dot = NULL;
                if (udps[numUdp].arg[iarg].val == NULL) return EGADS_MALLOC;

                memcpy(udps[numUdp].arg[iarg].val, udps[0].arg[iarg].val, isize*sizeof(double));
            } else {
                udps[numUdp].arg[iarg].val = NULL;
                udps[numUdp].arg[iarg].dot = NULL;
            }
        }
    }

    status = EGADS_SUCCESS;

cleanup:
    return status;
}
