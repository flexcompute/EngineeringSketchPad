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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

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

    int    ibody=0, nnode=0, nedge=0, nface=0, test;
    double dihedral=0;
    char   *arg1=NULL, *arg2=NULL, *arg3=NULL, *arg4=NULL;
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
        if (ESP->batch == 0) {
            tim_bcst("ereped", response);
        } else {
            printf("%s\n", response);
        }

    /* "checkEBody|ibody|nnode|nedge|nface|" */
    } else if (strncmp(command, "checkEBody|", 11) == 0) {

        /* extract arguments */
        GetToken(command, 1, '|', &arg1);
        SPLINT_CHECK_FOR_NULL(arg1);
        if (strlen(arg1) > 0) ibody = strtol(arg1, &pEnd, 10);
        FREE(arg1);
        
        GetToken(command, 2, '|', &arg2);
        SPLINT_CHECK_FOR_NULL(arg2);
        if (strlen(arg2) > 0) nnode = strtol(arg2, &pEnd, 10);
        FREE(arg2);
        
        GetToken(command, 3, '|', &arg3);
        SPLINT_CHECK_FOR_NULL(arg3);
        if (strlen(arg3) > 0) nedge = strtol(arg3, &pEnd, 10);
        FREE(arg3);
        
        GetToken(command, 4, '|', &arg4);
        SPLINT_CHECK_FOR_NULL(arg4);
        if (strlen(arg4) > 0) nface = strtol(arg4, &pEnd, 10);
        FREE(arg4);

        /* make sure valid ibody is given */
        if (ibody < 1 || ibody > ESP->MODL->nbody) {
            printf("Invalid Body index (%d)\n", ibody);
            status = OCSM_BODY_NOT_FOUND;
            goto cleanup;
        }

        /* make sure we have an EBody */
        if (ESP->MODL->body[ibody].eebody == NULL) {
            printf("Body %d does not have an EBody\n", ibody);
            status = OCSM_BODY_NOT_FOUND;
            goto cleanup;
        }

        /* check number of Nodes */
        status = EG_getBodyTopos(ESP->MODL->body[ibody].eebody, NULL, NODE, &test, NULL);
        CHECK_STATUS(EG_getBodyTopos);

        if (test != nnode) {
            printf("Have %d Nodes but was expecting %d\n", test, nnode);
            status = OCSM_NODE_NOT_FOUND;
            goto cleanup;
        }

        /* check number of EEdges */
        status = EG_getBodyTopos(ESP->MODL->body[ibody].eebody, NULL, EEDGE, &test, NULL);
        CHECK_STATUS(EG_getBodyTopos);

        if (test != nedge) {
            printf("Have %d Edges but was expecting %d\n", test, nedge);
            status = OCSM_EDGE_NOT_FOUND;
            goto cleanup;
        }
            
        /* check number of EFaces */
        status = EG_getBodyTopos(ESP->MODL->body[ibody].eebody, NULL, EFACE, &test, NULL);
        CHECK_STATUS(EG_getBodyTopos);

        if (test != nface) {
            printf("Have %d Faces but was expecting %d\n", test, nface);
            status = OCSM_FACE_NOT_FOUND;
            goto cleanup;
        }
    }

cleanup:
    FREE(arg1);
    FREE(arg2);
    FREE(arg3);
    FREE(arg4);
    
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
