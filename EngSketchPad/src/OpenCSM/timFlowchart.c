/*
 ************************************************************************
 *                                                                      *
 * timFlowchart -- Tool Integration Module for CAPS flowcharts          *
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
#include "common.h"
#include "OpenCSM.h"
#include "caps.h"
#include "tim.h"
#include "wsserver.h"
#include "emp.h"

extern int caps_outputObjects(capsObj probObj, char **stream);


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        /*@unused@*/void  *data)                    /* (in)  component name */
{
    int    status=SUCCESS;              /* (out) return status */

    char   *json_stream=NULL;
    FILE   *fp;
    
    ROUTINE(timLoad(flowchart));

    /* --------------------------------------------------------------- */

    if (ESP == NULL) {
        goto cleanup;
    }

    /* get the JSON string from CAPS that is needed by the flowcharter */
    status = caps_outputObjects(ESP->CAPS, &json_stream);
    if (status < SUCCESS) goto cleanup;

    SPLINT_CHECK_FOR_NULL(json_stream);

    fp = fopen("../ESP/ESP-flowchart-data.js", "w");
    if (fp != NULL) {
        fprintf(fp, "%s\n", json_stream);
        fclose(fp);
    }
    
    EG_free(json_stream);

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
timMesg(
/*@unused@*/esp_T *ESP,                 /* (in)  pointer to ESP structure */
        char  command[])                /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timMesg(flowchart));

    /* --------------------------------------------------------------- */

    /* "show|" */
    if        (strncmp(command, "show", 4) == 0) {

        /* tell ESP that we are starting an overlay */
        tim_bcst("flowchart", "overlayBeg|pyscript|flowchart|");

        tim_bcst("flowchart", "timMesg|flowchart|show");

        /* automatically hold ESP */
        tim_hold("pyscript", "flowchart");
    }

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(/*@unused@*/esp_T *ESP)                     /* (in)  pointer to TI Mstructure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    ROUTINE(timSave(flowchart));

    /* --------------------------------------------------------------- */

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

    ROUTINE(timQuit(flowchart));

    /* --------------------------------------------------------------- */

//cleanup:
    return status;
}
