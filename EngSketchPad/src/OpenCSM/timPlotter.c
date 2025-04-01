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
 * Copyright (C) 2013/2025  John F. Dannenhoffer, III (Syracuse University)
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

/*
   timLoad|plotter|

   timMesg|plotter|new|title|xlabel|ylabel|ylabel2|

   timMesg|plotter|add|xvalue1;xvalue2;...|yvalue1;yvalue2;...|type|

   type:
     r red      - solid     o circle       2 ylabel2
     g green    : dotted    x x-mark
     b blue     _ dashed    + plus
     c cyan     ; dot-dash  * star
     m magenta              s square
     y yellow               ^ triangle-up
     k black                v triangle-down
     w white

   timMesg|plotter|show|
   timMesg|plotter|show|nohold|
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "OpenCSM.h"
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
    char   *ylabel2;
    int    nline;
    line_T *lines;
} plotter_T;


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

    if (ESP == NULL) {
        printf("ERROR:: cannot run timPlotter without serveESP\n");
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    /* create the plotter_T structure */
    if (ESP->nudata >= MAX_TIM_NESTING) {
        printf("ERROR:: cannot nest more than %d TIMs\n", MAX_TIM_NESTING);
        exit(0);
    }

    ESP->nudata++;
    MALLOC(ESP->udata[ESP->nudata-1], plotter_T, 1);

    strcpy(ESP->timName[ESP->nudata-1], "plotter");

    plotter = (plotter_T *) (ESP->udata[ESP->nudata-1]);

    /* initialize the structure */
    plotter->title   = NULL;
    plotter->xlabel  = NULL;
    plotter->ylabel  = NULL;
    plotter->ylabel2 = NULL;
    plotter->nline   = 0;
    plotter->lines   = NULL;

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
    char   *arg1=NULL, *arg2=NULL, *arg3=NULL, *arg4=NULL, *xystr=NULL;
    char   *pEnd;
    int    nresponse=0;
    char   *response=NULL;
    void   *realloc_temp = NULL;            /* used by RALLOC macro */

    plotter_T *plotter = (plotter_T *) (ESP->udata[ESP->nudata-1]);

    ROUTINE(timMesg(plotter));

    /* --------------------------------------------------------------- */

    /* "new|title|xlabel|ylabel|ylabel2|" */
    if (strncmp(command, "new|", 4) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        SPLINT_CHECK_FOR_NULL(arg1);

        GetToken(command, 2, '|', &arg2);
        SPLINT_CHECK_FOR_NULL(arg2);

        GetToken(command, 3, '|', &arg3);
        SPLINT_CHECK_FOR_NULL(arg3);

        GetToken(command, 4, '|', &arg4);
        SPLINT_CHECK_FOR_NULL(arg4);

        /* delete any previous plots */
        for (iline = 0; iline < plotter->nline; iline++) {
            FREE(plotter->lines[iline].x     );
            FREE(plotter->lines[iline].y     );
            FREE(plotter->lines[iline].style);
        }

        FREE(plotter->title  );
        FREE(plotter->xlabel );
        FREE(plotter->ylabel );
        FREE(plotter->ylabel2);
        FREE(plotter->lines  );

        plotter->nline = 0;

        /* store title, xlabel, and ylabel */
        MALLOC(plotter->title,  char, strlen(arg1)+1);
        MALLOC(plotter->xlabel, char, strlen(arg2)+1);
        MALLOC(plotter->ylabel, char, strlen(arg3)+1);

        strcpy(plotter->title,  arg1);
        strcpy(plotter->xlabel, arg2);
        strcpy(plotter->ylabel, arg3);

        /* stroe ylabel if given */
        if (arg4 != NULL) {
            MALLOC(plotter->ylabel2, char, strlen(arg4)+1);

            strcpy(plotter->ylabel2, arg4);
        }

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

        plotter = (plotter_T *) (ESP->udata[ESP->nudata-1]);
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

        if (plotter->ylabel2 == NULL) {
            snprintf(response, nresponse, "timMesg|plotter|show|{\"title\":\"%s\", \"xlabel\":\"%s\", \"ylabel\":\"%s\", \"lines\":[",
                     plotter->title, plotter->xlabel, plotter->ylabel);
        } else {
            snprintf(response, nresponse, "timMesg|plotter|show|{\"title\":\"%s\", \"xlabel\":\"%s\", \"ylabel\":\"%s\", \"ylabel2\":\"%s\", \"lines\":[",
                     plotter->title, plotter->xlabel, plotter->ylabel, plotter->ylabel2);
        }

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

        /* automatically hold ESP if the "nohold" option is not given */
        if (strncmp(command, "show|nohold|", 12) != 0) {
            tim_hold("pyscript", "plotter");
        }
    }

cleanup:
    FREE(arg1);
    FREE(arg2);
    FREE(arg3);
    FREE(arg4);
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

    int    i, iline;

    plotter_T *plotter;

    ROUTINE(timQuit(plotter));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "plotter") != 0) {
        printf("WARNING:: TIM on top of stack is not \"plotter\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    } else {
        plotter = (plotter_T *)(ESP->udata[ESP->nudata-1]);
    }

    /* do nothing if we have already cleared plotter data */
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

    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    tim_bcst("plotter", "timQuit|plotter|");

cleanup:
    return status;
}
