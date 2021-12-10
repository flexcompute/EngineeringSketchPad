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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "egads.h"
#include "OpenCSM.h"
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
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* does nothing, since all done in ESP-ereped.js */

//cleanup:
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
timQuit(/*@unused@*/esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* does nothing, since all done in ESP-ereped.js */

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
        char   command[],               /* (in)  command */
        int    *max_resp_len,           /* (in)  length of response */
        char   *response[])             /* (out) response */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    ibody;
    double dihedral;
    char   arg1[MAX_EXPR_LEN], arg2[MAX_EXPR_LEN], arg3[MAX_EXPR_LEN];
    char   *pEnd;

    /* --------------------------------------------------------------- */

    /* "makeEBody|ibody|dihedral|ents|" */
    if (strncmp(command, "makeEBody|", 10) == 0) {

        /* extract arguments */
        ibody    = 0;
        dihedral = 0;
        if (GetToken(command, 1, '|', arg1)) ibody    = strtol(arg1, &pEnd, 10);
        if (GetToken(command, 2, '|', arg2)) dihedral = strtod(arg2, &pEnd);
        GetToken(command, 3, '|', arg3);

        /* make or delete the EBody (by modifying MODL) */
        if (strcmp(arg3, ".") == 0) {
            status = ocsmMakeEBody(ESP->MODL, ibody, dihedral, NULL);
        } else {
            status = ocsmMakeEBody(ESP->MODL, ibody, dihedral, arg3);
        }

        /* build the response */
        if (status == SUCCESS) {
            snprintf(*response, *max_resp_len, "makeEBody|");
        } else {
            snprintf(*response, *max_resp_len, "makeEBody|ERROR:: unable to make EBody(s) status=%d", status);
        }
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timGetModl - get the active MODL                                  */
/*                                                                     */
/***********************************************************************/

int
timGetModl(/*@unused@*/void **myModl)               /* (out) pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSetModl - set the active MODL                                  */
/*                                                                     */
/***********************************************************************/

int
timSetModl(/*@unused@*/void *myModl)                /* (in)  pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timViewModl - view the active MODL in serveESP                    */
/*                                                                     */
/***********************************************************************/

int
timViewModl(/*@unused@*/void *myModl)               /* (in)  pointer to active MODL */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}
