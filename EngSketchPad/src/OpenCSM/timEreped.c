/*
 ************************************************************************
 *                                                                      *
 * timEreped -- Tool Integration Module for Erep editor                 *
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
#include "OpenCSM.h"
#include "common.h"
#include "tim.h"


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(/*@unused@*/esp_T *ESP,                     /* (in)  pointer to ESP structure */
        /*@unused@*/void  *data)                    /* (in)  user-supplied data */
{
    int    status;                      /* (out) return status */

    /* does nothing, since all done in ESP-ereped.js */

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
timMesg(esp_T  *ESP,                    /* (in)  pointer to ESP structure */
        char   command[])               /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    ibody=0;
    double dihedral=0;
    char   *arg1=NULL, *arg2=NULL, *arg3=NULL;
    char   response[MAX_EXPR_LEN];
    char   *pEnd;

    ROUTINE("timMesg(ereped)");

    /* --------------------------------------------------------------- */

    /* "makeEBody|ibody|dihedral|ents|" */
    if (strncmp(command, "makeEBody|", 10) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        SPLINT_CHECK_FOR_NULL(arg1);
        if (strlen(arg1) > 0) ibody    = strtol(arg1, &pEnd, 10);
        FREE(arg1);
        
        GetToken(command, 2, '|', &arg2);
        SPLINT_CHECK_FOR_NULL(arg2);
        if (strlen(arg2) > 0) dihedral = strtod(arg2, &pEnd);
        FREE(arg2);
        
        GetToken(command, 3, '|', &arg3);
        SPLINT_CHECK_FOR_NULL(arg3);

        /* make or delete the EBody (by modifying MODL) */
        if (strcmp(arg3, ".") == 0) {
            status = ocsmMakeEBody(ESP->MODL, ibody, dihedral, NULL);
        } else {
            status = ocsmMakeEBody(ESP->MODL, ibody, dihedral, arg3);
        }

        FREE(arg1);

        /* build the response */
        if (status == SUCCESS) {
            snprintf(response, MAX_EXPR_LEN, "timMesg|ereped|makeEBody|");
        } else {
            snprintf(response, MAX_EXPR_LEN, "timMesg|ereped|makeEBody|ERROR:: unable to make EBody(s) status=%d", status);
        }
        tim_bcst("ereped", response);
    }

cleanup:
    FREE(arg1);
    FREE(arg2);
    FREE(arg3);
    
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(/*@unused@*/esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* does nothing, since all done in ESP-ereped.js */

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

    /* does nothing, since all done in ESP-ereped.js */

//cleanup:
    return status;
}
