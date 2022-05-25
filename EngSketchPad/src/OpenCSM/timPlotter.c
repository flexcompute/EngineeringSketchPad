/*
 ************************************************************************
 *                                                                      *
 * timPlotter -- Tool Integration Module for 2D plotter                 *
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

typedef struct {
    int    npnt;
    double *x;
    double *y;
    char   *style;
} line_T;

typedef struct {
    char   *title;
    char   *xlabel;
    char   *ylabel;
    int    nline;
    line_T *lines;
} plotter_T;

static void *realloc_temp=NULL;              /* used by RALLOC macro */


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
/*@unused@*/void  *udata)               /* (in)  not used */
{
    int    status=0;                    /* (out) return status */

    plotter_T *plotter;

    ROUTINE(timLoad(plotter));

    /* --------------------------------------------------------------- */

    outLevel = ocsmSetOutLevel(-1);

    /* create the plotter_T structure and initialize it */
    MALLOC(ESP->udata2, plotter_T, 1);

    plotter = (plotter_T *) ESP->udata2;

    /* initialize it */
    plotter->title  = NULL;
    plotter->xlabel = NULL;
    plotter->ylabel = NULL;
    plotter->nline  = 0;
    plotter->lines  = NULL;

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

    int    iline, i;
    char   *arg1=NULL, *arg2=NULL, *arg3=NULL, *xystr=NULL;
    char   *pEnd;
    int    nresponse=0;
    char   *response=NULL;

    plotter_T *plotter = (plotter_T *) ESP->udata2;

    ROUTINE(timMesg(plotter));

    /* --------------------------------------------------------------- */

    /* "new|title|xlabel|ylabel|" */
    if (strncmp(command, "new|", 4) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        SPLINT_CHECK_FOR_NULL(arg1);

        GetToken(command, 2, '|', &arg2);
        SPLINT_CHECK_FOR_NULL(arg2);

        GetToken(command, 3, '|', &arg3);
        SPLINT_CHECK_FOR_NULL(arg3);

        /* delete any previous plots */
        for (iline = 0; iline < plotter->nline; iline++) {
            FREE(plotter->lines[iline].x     );
            FREE(plotter->lines[iline].y     );
            FREE(plotter->lines[iline].style);
        }

        FREE(plotter->title );
        FREE(plotter->xlabel);
        FREE(plotter->ylabel);
        FREE(plotter->lines );

        plotter->nline = 0;

        /* store title, xlabel, and ylabel */
        MALLOC(plotter->title,  char, strlen(arg1)+1);
        MALLOC(plotter->xlabel, char, strlen(arg2)+1);
        MALLOC(plotter->ylabel, char, strlen(arg2)+1);

        strcpy(plotter->title,  arg1);
        strcpy(plotter->xlabel, arg2);
        strcpy(plotter->ylabel, arg3);

    /* "add|x1;x2;...|y1;y2;...|style|" */
    } else if (strncmp(command, "add|", 4) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        SPLINT_CHECK_FOR_NULL(arg1);

        GetToken(command, 2, '|', &arg2);
        SPLINT_CHECK_FOR_NULL(arg2);

        GetToken(command, 3, '|', &arg3);
        SPLINT_CHECK_FOR_NULL(arg3);

        /* initialize a new line */
        RALLOC(plotter->lines, line_T, plotter->nline+1);

        plotter = (plotter_T *) ESP->udata2;
        iline   = plotter->nline;

        plotter->lines[iline].x     = NULL;
        plotter->lines[iline].y     = NULL;
        plotter->lines[iline].style = NULL;

        plotter->lines[iline].npnt = 1;
        for (i = 1; i < strlen(arg1)-1; i++) {
            if (arg1[i] == ';') {
                (plotter->lines[iline].npnt)++;
            }
        }

        strcat(arg1, ";");    /* needed so last token can be retrieved by GetToken */
        strcat(arg2, ";");

        MALLOC(plotter->lines[iline].x,     double, plotter->lines[iline].npnt);
        MALLOC(plotter->lines[iline].y,     double, plotter->lines[iline].npnt);
        MALLOC(plotter->lines[iline].style, char,   strlen(arg3)+1);

        for (i = 0; i < plotter->lines[iline].npnt; i++) {
            GetToken(arg1, i, ';', &xystr);
            SPLINT_CHECK_FOR_NULL(xystr);
            plotter->lines[iline].x[i] = strtod(xystr, &pEnd);
            FREE(xystr);

            GetToken(arg2, i, ';', &xystr);
            SPLINT_CHECK_FOR_NULL(xystr);
            plotter->lines[iline].y[i] = strtod(xystr, &pEnd);
            FREE(xystr);
        }

        strcpy(plotter->lines[iline].style, arg3);

        (plotter->nline)++;

    /* "show" */
    } else if (strncmp(command, "show", 4) == 0) {

        if (plotter->nline <= 0) goto cleanup;

        /* tell ESP that we are starting an overlay */
        tim_bcst("plotter", "overlayBeg|pyscript|plotter|");

        /* build the json stream */
        nresponse = 10000;
        MALLOC(response, char, nresponse);

        snprintf(response, nresponse, "timMesg|plotter|show|{\"title\":\"%s\", \"xlabel\":\"%s\", \"ylabel\":\"%s\", \"lines\":[",
                 plotter->title, plotter->xlabel, plotter->ylabel);

        MALLOC(xystr, char, 25);

        for (iline = 0; iline < plotter->nline; iline++) {
            strcat(response, "{\"x\":[");
            for (i = 0; i < plotter->lines[iline].npnt; i++) {
                snprintf(xystr, 24, "%e", plotter->lines[iline].x[i]);
                strcat(response, xystr);
                if (i < plotter->lines[iline].npnt-1) {
                    strcat(response, ",");
                } else {
                    strcat(response, "],\"y\":[");
                }
            }
            for (i = 0; i < plotter->lines[iline].npnt; i++) {
                snprintf(xystr, 24, "%e", plotter->lines[iline].y[i]);
                strcat(response, xystr);
                if (i < plotter->lines[iline].npnt-1) {
                    strcat(response, ",");
                } else {
                    strcat(response, "],\"style\":\"");
                    strcat(response, plotter->lines[iline].style);
                    strcat(response, "\"}");
                }
            }
            if (iline < plotter->nline-1) {
                strcat(response, ",");
            } else {
                strcat(response, "]}|");
            }
        }

        tim_bcst("plotter", response);
    }

cleanup:
    FREE(arg1);
    FREE(arg2);
    FREE(arg3);
    FREE(xystr);
    FREE(response);

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

    ROUTINE(timSave(plotter));

    /* --------------------------------------------------------------- */

    /* same as quit */
    status = timQuit(ESP, 0);

//cleanup:
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

    int    iline;

    plotter_T *plotter = (plotter_T *)(ESP->udata2);

    ROUTINE(timQuit(plotter));

    /* --------------------------------------------------------------- */

    /* do nothing if we have already cleared plotetr data */
    if (plotter == NULL) goto cleanup;

    /* delete any previous plot data */
    for (iline = 0; iline < plotter->nline; iline++) {
        FREE(plotter->lines[iline].x    );
        FREE(plotter->lines[iline].y    );
        FREE(plotter->lines[iline].style);
    }

    FREE(plotter->title );
    FREE(plotter->xlabel);
    FREE(plotter->ylabel);
    FREE(plotter->lines );

    plotter->nline = 0;

    FREE(ESP->udata2);

    tim_bcst("plotter", "timQuit|plotter|");

cleanup:
    return status;
}
