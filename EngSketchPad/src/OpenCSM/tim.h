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
 * Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
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

#ifdef WIN32
    #define snprintf _snprintf
#endif

#include "esp.h"

/*
 ************************************************************************
 *                                                                      *
 * Prototypes for functions to be available to be called by serveESP    *
 *                                                                      *
 ************************************************************************
 */

/* open a tim instance */
extern int
tim_load(char   timName[],              /* (in)  name of TIM */
         esp_T  *ESP,                   /* (in)  pointer to ESP structure */
/*@null@*/void  *data);                 /* (in)  user-provided data */

/* get command, process, and return response */
/*@null@*/extern int
tim_mesg(char   timName[],              /* (in)  name of TIM */
         char   command[]);             /* (in)  command to process */

/* set a lock on all TIMs */
extern void
tim_lock();

/* hold a TIM until the lock is lifted */
extern int
tim_hold(char   timName[],              /* (in)  name of TIM to hold */
         char   overlay[]);             /* (in)  name of overlay to lift hold */

/* lift the hold on a TIM */
extern int
tim_lift(char   timName[]);             /* (in)  name of TIM to unlock */

/* broadcast a message to all browsers */
extern int
tim_bcst(char   timName[],              /* (in)  name of TIM */
         char   text[]);                /* (in)  text to broadcast */

/* save TIM data and close TIM instance */
extern int
tim_save(char   timName[]);             /* (in)  name of TIM */

/* close TIM instance without saving */
extern int
tim_quit(char   timName[]);             /* (in)  name of TIM */

/* free up all TIM data */
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
timLoad(esp_T  *TIM,                    /* (in)  pointer to ESP structure */
/*@null@*/void *data);                  /* (in)  TIM-specific structure */

int
timMesg(esp_T  *TIM,                    /* (in)  pointer to ESP structure */
        char   command[]);              /* (in)  caommand to execute */

int
timSave(esp_T  *TIM);                   /* (in)  pointer to ESP structure */

int
timQuit(esp_T  *TIM,                    /* (in)  pointer to ESP structure */
        int    unload);                 /* (in)  flag to unload */

/*
 ************************************************************************
 *                                                                      *
 * Prototypes for helper functions available to TIMs                    *
 *                                                                      *
 ************************************************************************
 */

extern int
GetToken(char text[],                   /* (in)  original string */
         int  nskip,                    /* (in)  number of separators to skip */
         char sep,                      /* (in)  separator character */
         char *token[]);                /* (in)  token (freeable) */

// values for timState[]
#define TIM_INACTIVE   0
#define TIM_LOADING    1
#define TIM_READY      2
#define TIM_EXECUTING  3
#define TIM_CLOSING    4

// value for TIM errors

#endif
