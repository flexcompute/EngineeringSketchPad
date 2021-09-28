/*
 ************************************************************************
 *                                                                      *
 * tim.h -- header for the Tool Integration Modules                     *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2021  John F. Dannenhoffer, III (Syracuse University)
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

#ifndef TIM_H
#define TIM_H

#include "wsserver.h"

/* define needed for WIN32 */
#ifdef WIN32
    #define snprintf _snprintf
#endif

/*
 ************************************************************************
 *                                                                      *
 * Prototypes for functions to be available to be called by serveESP    *
 *                                                                      *
 ************************************************************************
 */

typedef struct {
    modl_T    *MODL;                    /* pointer to OpenCSM MODL */
    wvContext *ctxt;                    /* WebViewer context */
    float     sgFocus[4];               /* scene graph focus */
    void      *udata;                   /* pointer to user data */
} tim_T;

/* open a tim instance */
extern int
tim_load(char   timName[],              /* (in)  name of tim */
         tim_T  *TIM,                   /* (in)  pointer to TIM structure */
/*@null@*/void  *data);                 /* (in)  user-provided data */

/* save tim data and close tim instance */
extern int
tim_save(char   timName[]);             /* (in)  name of tim */

/* close tim instance without saving */
extern int
tim_quit(char   timName[]);             /* (in)  name of tim */

/* get command, process, and return response */
extern int
tim_mesg(char   timName[],              /* (in)  name of tim */
         char   command[],              /* (in)  command to process */
         int    len_response,           /* (in)  length of response */
         char   response[]);            /* (out) response */

/* return pointer to TIM structure */
extern int
tim_data(char   timName[],
         tim_T  **TIM);

/* free up all tim data */
extern void
tim_free();

/*
 ************************************************************************
 *                                                                      *
 * Prototypes for functions to be provided by timXxx.c                  *
 *    (see above for argument descriptions)                             *
 *                                                                      *
 ************************************************************************
 */

int
timLoad(tim_T  *TIM,
/*@null@*/void *data);

int
timSave(tim_T  *TIM);

int
timQuit(tim_T  *TIM);

int
timMesg(tim_T  *TIM,
        char   command[],
        int    len_response,
        char   response[]);

/*
 ************************************************************************
 *                                                                      *
 * Prototypes for helper functions available to TIMs                    *
 *                                                                      *
 ************************************************************************
 */

extern int GetToken(char text[], int nskip, char sep, char token[]);

#endif
