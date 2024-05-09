/*
 ************************************************************************
 *                                                                      *
 * udpUtilities.h -- header for udp*.c                                  *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2024  John F. Dannenhoffer, III (Syracuse University)
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

#ifndef _UDPUTILITIES_H_
#define _UDPUTILITIES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "egads.h"
#include "udp.h"

//$$$#define DEBUG_UDP

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define  udpExecute(A,B,C,D)         UdpExecute(    A, B, C, D,       int *NumUdp, udp_T **Udps)
#define  udpSensitivity(A,B,C,D,E,F) UdpSensitivity(A, B, C, D, E, F, int *NumUdp, udp_T  *udps)
#define  udpMesh(A,B,C,D,E,F)        UdpMesh(       A, B, C, D, E, F, int *NumUdp, udp_T  *udps)
#define  cacheUdp(A)                 CacheUdp(      A,                     NumUdp, Udps); udps = *Udps
#define  numUdp  (*NumUdp)

#endif  /* _UDPUTILITIES_H_ */

