/*
 ************************************************************************
 *                                                                      *
 * esp.h -- structure for sharing top-level pointers                    *
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

#ifndef ESP_H
#define ESP_H

#include "wsserver.h"

#define MAX_TIM_NESTING  10

typedef struct {
    ego       EGADS;                    /* pointer to EGADS object */
    modl_T    *MODL;                    /* pointer to OpenCSM MODL */
    modl_T    *MODLorig;                /* pointer to OpenCSM MODL before CAPS */
    void      *CAPS;                    /* capsProject */
    wvContext *cntxt;                   /* WebViewer context */
    float     sgFocus[4];               /* scene graph focus */
    void      *sgMutex;                 /* mutex associated with scene graphs */
    int       curTim;                   /* ID of current TIM (or -1 if serveESP) */
    int       nudata;                   /* number of nested TIMs */
    void      *udata[MAX_TIM_NESTING];  /* array  of nested TIMs */
    char      timName[MAX_TIM_NESTING][10];  /* names of nested TIM */
} esp_T;

#endif
