/*
 ************************************************************************
 *                                                                      *
 * timSlugs -- Tool Integration Module for SLUGS                        *
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
#include "tim.h"


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(/*@unused@*/esp_T *TIM,                     /* (in)  pointer to TIM structure */
        /*@unused@*/void  *data)                    /* (in)  user-supplied data */
{
    int    status;                      /* (out) return status */

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
timMesg(/*@unused@*/esp_T *TIM,                     /* (in)  pointer to TIM structure */
        /*@unused@*/char  command[])                /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(/*@unused@*/esp_T *TIM)                     /* (in)  pointer to TIM structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(/*@unused@*/esp_T *TIM,                     /* (in)  pointer to TIM structure */
        /*@unused@*/int   unload)                   /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

//cleanup:
    return status;
}
