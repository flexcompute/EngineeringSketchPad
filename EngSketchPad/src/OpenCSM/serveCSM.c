/*
 ************************************************************************
 *                                                                      *
 * serveCSM.c -- server for driving OpenCSM                             *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2012/2022  John F. Dannenhoffer, III (Syracuse University)
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
#include <time.h>
#include <unistd.h>
#include <assert.h>

#define PLUGS_PRUNE              0      // if >0, prune points and write new.cloud
#define PLUGS_CREATE_CSM_FILES   0      // create plugs_pass_xx.csm files
#define PLUGS_CREATE_FINAL_PLOT  0      // create final.plot        file
#define PLUGS_CREATE_EGADS_FILES 0      // create plugs_x_y.egads   files (use data/plugs/showHistory.csm to see)
#define PLUGS_SHOW_SENSITIVITY   0      // show sensitivity in final results via inverse evaluations */
//#define MIXED_TESS               1      // display mixed tri/quads tessellations

#define TESS_PARAM_0             0.0250
#define TESS_PARAM_1             0.0075
#define TESS_PARAM_2             20.0

#ifdef WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <winsock2.h>
    #define snprintf     _snprintf
    #define SLEEP(msec)  Sleep(msec)
    #define SLASH '\\'
#else
    #include <unistd.h>
    #define SLEEP(msec)  usleep(1000*msec)
    #define SLASH '/'
#endif

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define STRNCPY(A, B, LEN) strncpy(A, B, LEN); A[LEN-1] = '\0';

#include "common.h"
#include "OpenCSM.h"

#include "udp.h"
#include "egg.h"
#include "emp.h"

#include "wsserver.h"

/***********************************************************************/
/*                                                                     */
/* macros (including those that go along with common.h)                */
/*                                                                     */
/***********************************************************************/

static void *realloc_temp=NULL;            /* used by RALLOC macro */

#define  RED(COLOR)      (float)(COLOR / 0x10000        ) / (float)(255)
#define  GREEN(COLOR)    (float)(COLOR / 0x00100 % 0x100) / (float)(255)
#define  BLUE(COLOR)     (float)(COLOR           % 0x100) / (float)(255)

/***********************************************************************/
/*                                                                     */
/* structures                                                          */
/*                                                                     */
/***********************************************************************/

#define BlueWhiteRed

/* blue-white-red spectrum */
#ifdef  BlueWhiteRed
static float color_map[256*3] =
{ 0.0000, 0.0000, 1.0000,    0.0078, 0.0078, 1.0000,   0.0156, 0.0156, 1.0000,    0.0234, 0.0234, 1.0000,
  0.0312, 0.0312, 1.0000,    0.0391, 0.0391, 1.0000,   0.0469, 0.0469, 1.0000,    0.0547, 0.0547, 1.0000,
  0.0625, 0.0625, 1.0000,    0.0703, 0.0703, 1.0000,   0.0781, 0.0781, 1.0000,    0.0859, 0.0859, 1.0000,
  0.0938, 0.0938, 1.0000,    0.1016, 0.1016, 1.0000,   0.1094, 0.1094, 1.0000,    0.1172, 0.1172, 1.0000,
  0.1250, 0.1250, 1.0000,    0.1328, 0.1328, 1.0000,   0.1406, 0.1406, 1.0000,    0.1484, 0.1484, 1.0000,
  0.1562, 0.1562, 1.0000,    0.1641, 0.1641, 1.0000,   0.1719, 0.1719, 1.0000,    0.1797, 0.1797, 1.0000,
  0.1875, 0.1875, 1.0000,    0.1953, 0.1953, 1.0000,   0.2031, 0.2031, 1.0000,    0.2109, 0.2109, 1.0000,
  0.2188, 0.2188, 1.0000,    0.2266, 0.2266, 1.0000,   0.2344, 0.2344, 1.0000,    0.2422, 0.2422, 1.0000,
  0.2500, 0.2500, 1.0000,    0.2578, 0.2578, 1.0000,   0.2656, 0.2656, 1.0000,    0.2734, 0.2734, 1.0000,
  0.2812, 0.2812, 1.0000,    0.2891, 0.2891, 1.0000,   0.2969, 0.2969, 1.0000,    0.3047, 0.3047, 1.0000,
  0.3125, 0.3125, 1.0000,    0.3203, 0.3203, 1.0000,   0.3281, 0.3281, 1.0000,    0.3359, 0.3359, 1.0000,
  0.3438, 0.3438, 1.0000,    0.3516, 0.3516, 1.0000,   0.3594, 0.3594, 1.0000,    0.3672, 0.3672, 1.0000,
  0.3750, 0.3750, 1.0000,    0.3828, 0.3828, 1.0000,   0.3906, 0.3906, 1.0000,    0.3984, 0.3984, 1.0000,
  0.4062, 0.4062, 1.0000,    0.4141, 0.4141, 1.0000,   0.4219, 0.4219, 1.0000,    0.4297, 0.4297, 1.0000,
  0.4375, 0.4375, 1.0000,    0.4453, 0.4453, 1.0000,   0.4531, 0.4531, 1.0000,    0.4609, 0.4609, 1.0000,
  0.4688, 0.4688, 1.0000,    0.4766, 0.4766, 1.0000,   0.4844, 0.4844, 1.0000,    0.4922, 0.4922, 1.0000,
  0.5000, 0.5000, 1.0000,    0.5078, 0.5078, 1.0000,   0.5156, 0.5156, 1.0000,    0.5234, 0.5234, 1.0000,
  0.5312, 0.5312, 1.0000,    0.5391, 0.5391, 1.0000,   0.5469, 0.5469, 1.0000,    0.5547, 0.5547, 1.0000,
  0.5625, 0.5625, 1.0000,    0.5703, 0.5703, 1.0000,   0.5781, 0.5781, 1.0000,    0.5859, 0.5859, 1.0000,
  0.5938, 0.5938, 1.0000,    0.6016, 0.6016, 1.0000,   0.6094, 0.6094, 1.0000,    0.6172, 0.6172, 1.0000,
  0.6250, 0.6250, 1.0000,    0.6328, 0.6328, 1.0000,   0.6406, 0.6406, 1.0000,    0.6484, 0.6484, 1.0000,
  0.6562, 0.6562, 1.0000,    0.6641, 0.6641, 1.0000,   0.6719, 0.6719, 1.0000,    0.6797, 0.6797, 1.0000,
  0.6875, 0.6875, 1.0000,    0.6953, 0.6953, 1.0000,   0.7031, 0.7031, 1.0000,    0.7109, 0.7109, 1.0000,
  0.7188, 0.7188, 1.0000,    0.7266, 0.7266, 1.0000,   0.7344, 0.7344, 1.0000,    0.7422, 0.7422, 1.0000,
  0.7500, 0.7500, 1.0000,    0.7578, 0.7578, 1.0000,   0.7656, 0.7656, 1.0000,    0.7734, 0.7734, 1.0000,
  0.7812, 0.7812, 1.0000,    0.7891, 0.7891, 1.0000,   0.7969, 0.7969, 1.0000,    0.8047, 0.8047, 1.0000,
  0.8125, 0.8125, 1.0000,    0.8203, 0.8203, 1.0000,   0.8281, 0.8281, 1.0000,    0.8359, 0.8359, 1.0000,
  0.8438, 0.8438, 1.0000,    0.8516, 0.8516, 1.0000,   0.8594, 0.8594, 1.0000,    0.8672, 0.8672, 1.0000,
  0.8750, 0.8750, 1.0000,    0.8828, 0.8828, 1.0000,   0.8906, 0.8906, 1.0000,    0.8984, 0.8984, 1.0000,
  0.9062, 0.9062, 1.0000,    0.9141, 0.9141, 1.0000,   0.9219, 0.9219, 1.0000,    0.9297, 0.9297, 1.0000,
  0.9375, 0.9375, 1.0000,    0.9453, 0.9453, 1.0000,   0.9531, 0.9531, 1.0000,    0.9609, 0.9609, 1.0000,
  0.9688, 0.9688, 1.0000,    0.9766, 0.9766, 1.0000,   0.9844, 0.9844, 1.0000,    0.9922, 0.9922, 1.0000,
  1.0000, 1.0000, 1.0000,    1.0000, 0.9922, 0.9922,   1.0000, 0.9844, 0.9844,    1.0000, 0.9766, 0.9766,
  1.0000, 0.9688, 0.9688,    1.0000, 0.9609, 0.9609,   1.0000, 0.9531, 0.9531,    1.0000, 0.9453, 0.9453,
  1.0000, 0.9375, 0.9375,    1.0000, 0.9297, 0.9297,   1.0000, 0.9219, 0.9219,    1.0000, 0.9141, 0.9141,
  1.0000, 0.9062, 0.9062,    1.0000, 0.8984, 0.8984,   1.0000, 0.8906, 0.8906,    1.0000, 0.8828, 0.8828,
  1.0000, 0.8750, 0.8750,    1.0000, 0.8672, 0.8672,   1.0000, 0.8594, 0.8594,    1.0000, 0.8516, 0.8516,
  1.0000, 0.8438, 0.8438,    1.0000, 0.8359, 0.8359,   1.0000, 0.8281, 0.8281,    1.0000, 0.8203, 0.8203,
  1.0000, 0.8125, 0.8125,    1.0000, 0.8047, 0.8047,   1.0000, 0.7969, 0.7969,    1.0000, 0.7891, 0.7891,
  1.0000, 0.7812, 0.7812,    1.0000, 0.7734, 0.7734,   1.0000, 0.7656, 0.7656,    1.0000, 0.7578, 0.7578,
  1.0000, 0.7500, 0.7500,    1.0000, 0.7422, 0.7422,   1.0000, 0.7344, 0.7344,    1.0000, 0.7266, 0.7266,
  1.0000, 0.7188, 0.7188,    1.0000, 0.7109, 0.7109,   1.0000, 0.7031, 0.7031,    1.0000, 0.6953, 0.6953,
  1.0000, 0.6875, 0.6875,    1.0000, 0.6797, 0.6797,   1.0000, 0.6719, 0.6719,    1.0000, 0.6641, 0.6641,
  1.0000, 0.6562, 0.6562,    1.0000, 0.6484, 0.6484,   1.0000, 0.6406, 0.6406,    1.0000, 0.6328, 0.6328,
  1.0000, 0.6250, 0.6250,    1.0000, 0.6172, 0.6172,   1.0000, 0.6094, 0.6094,    1.0000, 0.6016, 0.6016,
  1.0000, 0.5938, 0.5938,    1.0000, 0.5859, 0.5859,   1.0000, 0.5781, 0.5781,    1.0000, 0.5703, 0.5703,
  1.0000, 0.5625, 0.5625,    1.0000, 0.5547, 0.5547,   1.0000, 0.5469, 0.5469,    1.0000, 0.5391, 0.5391,
  1.0000, 0.5312, 0.5312,    1.0000, 0.5234, 0.5234,   1.0000, 0.5156, 0.5156,    1.0000, 0.5078, 0.5078,
  1.0000, 0.5000, 0.5000,    1.0000, 0.4922, 0.4922,   1.0000, 0.4844, 0.4844,    1.0000, 0.4766, 0.4766,
  1.0000, 0.4688, 0.4688,    1.0000, 0.4609, 0.4609,   1.0000, 0.4531, 0.4531,    1.0000, 0.4453, 0.4453,
  1.0000, 0.4375, 0.4375,    1.0000, 0.4297, 0.4297,   1.0000, 0.4219, 0.4219,    1.0000, 0.4141, 0.4141,
  1.0000, 0.4062, 0.4062,    1.0000, 0.3984, 0.3984,   1.0000, 0.3906, 0.3906,    1.0000, 0.3828, 0.3828,
  1.0000, 0.3750, 0.3750,    1.0000, 0.3672, 0.3672,   1.0000, 0.3594, 0.3594,    1.0000, 0.3516, 0.3516,
  1.0000, 0.3438, 0.3438,    1.0000, 0.3359, 0.3359,   1.0000, 0.3281, 0.3281,    1.0000, 0.3203, 0.3203,
  1.0000, 0.3125, 0.3125,    1.0000, 0.3047, 0.3047,   1.0000, 0.2969, 0.2969,    1.0000, 0.2891, 0.2891,
  1.0000, 0.2812, 0.2812,    1.0000, 0.2734, 0.2734,   1.0000, 0.2656, 0.2656,    1.0000, 0.2578, 0.2578,
  1.0000, 0.2500, 0.2500,    1.0000, 0.2422, 0.2422,   1.0000, 0.2344, 0.2344,    1.0000, 0.2266, 0.2266,
  1.0000, 0.2188, 0.2188,    1.0000, 0.2109, 0.2109,   1.0000, 0.2031, 0.2031,    1.0000, 0.1953, 0.1953,
  1.0000, 0.1875, 0.1875,    1.0000, 0.1797, 0.1797,   1.0000, 0.1719, 0.1719,    1.0000, 0.1641, 0.1641,
  1.0000, 0.1562, 0.1562,    1.0000, 0.1484, 0.1484,   1.0000, 0.1406, 0.1406,    1.0000, 0.1328, 0.1328,
  1.0000, 0.1250, 0.1250,    1.0000, 0.1172, 0.1172,   1.0000, 0.1094, 0.1094,    1.0000, 0.1016, 0.1016,
  1.0000, 0.0938, 0.0938,    1.0000, 0.0859, 0.0859,   1.0000, 0.0781, 0.0781,    1.0000, 0.0703, 0.0703,
  1.0000, 0.0625, 0.0625,    1.0000, 0.0547, 0.0547,   1.0000, 0.0469, 0.0469,    1.0000, 0.0391, 0.0391,
  1.0000, 0.0312, 0.0312,    1.0000, 0.0234, 0.0234,   1.0000, 0.0156, 0.0156,    1.0000, 0.0078, 0.0078 };

/* blue-green-red spectrum */
#else
static float color_map[256*3] =
{ 0.0000, 0.0000, 1.0000,    0.0000, 0.0157, 1.0000,   0.0000, 0.0314, 1.0000,    0.0000, 0.0471, 1.0000,
  0.0000, 0.0627, 1.0000,    0.0000, 0.0784, 1.0000,   0.0000, 0.0941, 1.0000,    0.0000, 0.1098, 1.0000,
  0.0000, 0.1255, 1.0000,    0.0000, 0.1412, 1.0000,   0.0000, 0.1569, 1.0000,    0.0000, 0.1725, 1.0000,
  0.0000, 0.1882, 1.0000,    0.0000, 0.2039, 1.0000,   0.0000, 0.2196, 1.0000,    0.0000, 0.2353, 1.0000,
  0.0000, 0.2510, 1.0000,    0.0000, 0.2667, 1.0000,   0.0000, 0.2824, 1.0000,    0.0000, 0.2980, 1.0000,
  0.0000, 0.3137, 1.0000,    0.0000, 0.3294, 1.0000,   0.0000, 0.3451, 1.0000,    0.0000, 0.3608, 1.0000,
  0.0000, 0.3765, 1.0000,    0.0000, 0.3922, 1.0000,   0.0000, 0.4078, 1.0000,    0.0000, 0.4235, 1.0000,
  0.0000, 0.4392, 1.0000,    0.0000, 0.4549, 1.0000,   0.0000, 0.4706, 1.0000,    0.0000, 0.4863, 1.0000,
  0.0000, 0.5020, 1.0000,    0.0000, 0.5176, 1.0000,   0.0000, 0.5333, 1.0000,    0.0000, 0.5490, 1.0000,
  0.0000, 0.5647, 1.0000,    0.0000, 0.5804, 1.0000,   0.0000, 0.5961, 1.0000,    0.0000, 0.6118, 1.0000,
  0.0000, 0.6275, 1.0000,    0.0000, 0.6431, 1.0000,   0.0000, 0.6588, 1.0000,    0.0000, 0.6745, 1.0000,
  0.0000, 0.6902, 1.0000,    0.0000, 0.7059, 1.0000,   0.0000, 0.7216, 1.0000,    0.0000, 0.7373, 1.0000,
  0.0000, 0.7529, 1.0000,    0.0000, 0.7686, 1.0000,   0.0000, 0.7843, 1.0000,    0.0000, 0.8000, 1.0000,
  0.0000, 0.8157, 1.0000,    0.0000, 0.8314, 1.0000,   0.0000, 0.8471, 1.0000,    0.0000, 0.8627, 1.0000,
  0.0000, 0.8784, 1.0000,    0.0000, 0.8941, 1.0000,   0.0000, 0.9098, 1.0000,    0.0000, 0.9255, 1.0000,
  0.0000, 0.9412, 1.0000,    0.0000, 0.9569, 1.0000,   0.0000, 0.9725, 1.0000,    0.0000, 0.9882, 1.0000,
  0.0000, 1.0000, 0.9961,    0.0000, 1.0000, 0.9804,   0.0000, 1.0000, 0.9647,    0.0000, 1.0000, 0.9490,
  0.0000, 1.0000, 0.9333,    0.0000, 1.0000, 0.9176,   0.0000, 1.0000, 0.9020,    0.0000, 1.0000, 0.8863,
  0.0000, 1.0000, 0.8706,    0.0000, 1.0000, 0.8549,   0.0000, 1.0000, 0.8392,    0.0000, 1.0000, 0.8235,
  0.0000, 1.0000, 0.8078,    0.0000, 1.0000, 0.7922,   0.0000, 1.0000, 0.7765,    0.0000, 1.0000, 0.7608,
  0.0000, 1.0000, 0.7451,    0.0000, 1.0000, 0.7294,   0.0000, 1.0000, 0.7137,    0.0000, 1.0000, 0.6980,
  0.0000, 1.0000, 0.6824,    0.0000, 1.0000, 0.6667,   0.0000, 1.0000, 0.6510,    0.0000, 1.0000, 0.6353,
  0.0000, 1.0000, 0.6196,    0.0000, 1.0000, 0.6039,   0.0000, 1.0000, 0.5882,    0.0000, 1.0000, 0.5725,
  0.0000, 1.0000, 0.5569,    0.0000, 1.0000, 0.5412,   0.0000, 1.0000, 0.5255,    0.0000, 1.0000, 0.5098,
  0.0000, 1.0000, 0.4941,    0.0000, 1.0000, 0.4784,   0.0000, 1.0000, 0.4627,    0.0000, 1.0000, 0.4471,
  0.0000, 1.0000, 0.4314,    0.0000, 1.0000, 0.4157,   0.0000, 1.0000, 0.4000,    0.0000, 1.0000, 0.3843,
  0.0000, 1.0000, 0.3686,    0.0000, 1.0000, 0.3529,   0.0000, 1.0000, 0.3373,    0.0000, 1.0000, 0.3216,
  0.0000, 1.0000, 0.3059,    0.0000, 1.0000, 0.2902,   0.0000, 1.0000, 0.2745,    0.0000, 1.0000, 0.2588,
  0.0000, 1.0000, 0.2431,    0.0000, 1.0000, 0.2275,   0.0000, 1.0000, 0.2118,    0.0000, 1.0000, 0.1961,
  0.0000, 1.0000, 0.1804,    0.0000, 1.0000, 0.1647,   0.0000, 1.0000, 0.1490,    0.0000, 1.0000, 0.1333,
  0.0000, 1.0000, 0.1176,    0.0000, 1.0000, 0.1020,   0.0000, 1.0000, 0.0863,    0.0000, 1.0000, 0.0706,
  0.0000, 1.0000, 0.0549,    0.0000, 1.0000, 0.0392,   0.0000, 1.0000, 0.0235,    0.0000, 1.0000, 0.0078,
  0.0078, 1.0000, 0.0000,    0.0235, 1.0000, 0.0000,   0.0392, 1.0000, 0.0000,    0.0549, 1.0000, 0.0000,
  0.0706, 1.0000, 0.0000,    0.0863, 1.0000, 0.0000,   0.1020, 1.0000, 0.0000,    0.1176, 1.0000, 0.0000,
  0.1333, 1.0000, 0.0000,    0.1490, 1.0000, 0.0000,   0.1647, 1.0000, 0.0000,    0.1804, 1.0000, 0.0000,
  0.1961, 1.0000, 0.0000,    0.2118, 1.0000, 0.0000,   0.2275, 1.0000, 0.0000,    0.2431, 1.0000, 0.0000,
  0.2588, 1.0000, 0.0000,    0.2745, 1.0000, 0.0000,   0.2902, 1.0000, 0.0000,    0.3059, 1.0000, 0.0000,
  0.3216, 1.0000, 0.0000,    0.3373, 1.0000, 0.0000,   0.3529, 1.0000, 0.0000,    0.3686, 1.0000, 0.0000,
  0.3843, 1.0000, 0.0000,    0.4000, 1.0000, 0.0000,   0.4157, 1.0000, 0.0000,    0.4314, 1.0000, 0.0000,
  0.4471, 1.0000, 0.0000,    0.4627, 1.0000, 0.0000,   0.4784, 1.0000, 0.0000,    0.4941, 1.0000, 0.0000,
  0.5098, 1.0000, 0.0000,    0.5255, 1.0000, 0.0000,   0.5412, 1.0000, 0.0000,    0.5569, 1.0000, 0.0000,
  0.5725, 1.0000, 0.0000,    0.5882, 1.0000, 0.0000,   0.6039, 1.0000, 0.0000,    0.6196, 1.0000, 0.0000,
  0.6353, 1.0000, 0.0000,    0.6510, 1.0000, 0.0000,   0.6667, 1.0000, 0.0000,    0.6824, 1.0000, 0.0000,
  0.6980, 1.0000, 0.0000,    0.7137, 1.0000, 0.0000,   0.7294, 1.0000, 0.0000,    0.7451, 1.0000, 0.0000,
  0.7608, 1.0000, 0.0000,    0.7765, 1.0000, 0.0000,   0.7922, 1.0000, 0.0000,    0.8078, 1.0000, 0.0000,
  0.8235, 1.0000, 0.0000,    0.8392, 1.0000, 0.0000,   0.8549, 1.0000, 0.0000,    0.8706, 1.0000, 0.0000,
  0.8863, 1.0000, 0.0000,    0.9020, 1.0000, 0.0000,   0.9176, 1.0000, 0.0000,    0.9333, 1.0000, 0.0000,
  0.9490, 1.0000, 0.0000,    0.9647, 1.0000, 0.0000,   0.9804, 1.0000, 0.0000,    0.9961, 1.0000, 0.0000,
  1.0000, 0.9882, 0.0000,    1.0000, 0.9725, 0.0000,   1.0000, 0.9569, 0.0000,    1.0000, 0.9412, 0.0000,
  1.0000, 0.9255, 0.0000,    1.0000, 0.9098, 0.0000,   1.0000, 0.8941, 0.0000,    1.0000, 0.8784, 0.0000,
  1.0000, 0.8627, 0.0000,    1.0000, 0.8471, 0.0000,   1.0000, 0.8314, 0.0000,    1.0000, 0.8157, 0.0000,
  1.0000, 0.8000, 0.0000,    1.0000, 0.7843, 0.0000,   1.0000, 0.7686, 0.0000,    1.0000, 0.7529, 0.0000,
  1.0000, 0.7373, 0.0000,    1.0000, 0.7216, 0.0000,   1.0000, 0.7059, 0.0000,    1.0000, 0.6902, 0.0000,
  1.0000, 0.6745, 0.0000,    1.0000, 0.6588, 0.0000,   1.0000, 0.6431, 0.0000,    1.0000, 0.6275, 0.0000,
  1.0000, 0.6118, 0.0000,    1.0000, 0.5961, 0.0000,   1.0000, 0.5804, 0.0000,    1.0000, 0.5647, 0.0000,
  1.0000, 0.5490, 0.0000,    1.0000, 0.5333, 0.0000,   1.0000, 0.5176, 0.0000,    1.0000, 0.5020, 0.0000,
  1.0000, 0.4863, 0.0000,    1.0000, 0.4706, 0.0000,   1.0000, 0.4549, 0.0000,    1.0000, 0.4392, 0.0000,
  1.0000, 0.4235, 0.0000,    1.0000, 0.4078, 0.0000,   1.0000, 0.3922, 0.0000,    1.0000, 0.3765, 0.0000,
  1.0000, 0.3608, 0.0000,    1.0000, 0.3451, 0.0000,   1.0000, 0.3294, 0.0000,    1.0000, 0.3137, 0.0000,
  1.0000, 0.2980, 0.0000,    1.0000, 0.2824, 0.0000,   1.0000, 0.2667, 0.0000,    1.0000, 0.2510, 0.0000,
  1.0000, 0.2353, 0.0000,    1.0000, 0.2196, 0.0000,   1.0000, 0.2039, 0.0000,    1.0000, 0.1882, 0.0000,
  1.0000, 0.1725, 0.0000,    1.0000, 0.1569, 0.0000,   1.0000, 0.1412, 0.0000,    1.0000, 0.1255, 0.0000,
  1.0000, 0.1098, 0.0000,    1.0000, 0.0941, 0.0000,   1.0000, 0.0784, 0.0000,    1.0000, 0.0627, 0.0000,
  1.0000, 0.0471, 0.0000,    1.0000, 0.0314, 0.0000,   1.0000, 0.0157, 0.0000,    1.0000, 0.0000, 0.0000 };
#endif

/***********************************************************************/
/*                                                                     */
/* global variables                                                    */
/*                                                                     */
/***********************************************************************/

/* global variable holding a MODL */
static void       *modl      = NULL;   /* pointer to current MODL */
static int        addVerify  = 0;      /* =1 to create .csm_verify file */
static int        allVels    = 0;      /* =1 to compute Node/Edge/Faces vels */
static int        batch      = 0;      /* =0 to enable visualization */
static int        checkMass  = 0;      /* =1 to check mass properties */
static int        checkPara  = 0;      /* =1 to check for parallelism */
static int        dumpEgads  = 0;      /* =1 to dump Bodys as they are created */
static double     histDist   = 0;      /* >0 to generate histograms of distance from points to brep */
static int        loadEgads  = 0;      /* =1 to read Bodys instead of creating */
static int        onormal    = 0;      /* =1 to use orthonormal (not perspective) */
static int        outLevel   = 1;      /* default output level */
static int        plotCP     = 0;      /* =1 to plot Bspline control polygons */
static int        plugsOn    = 0;      /* =1 to enable PLUGS */
static int        plugs      = 0;      /* =0  run phase1 only */
                                       /* <0  run abs(plugs)        passes of phase2 */
                                       /* >0  run phase 1 and plugs passes of phase2 */
static int        printStack = 0;      /* =1 to print stack after every command */
static int        skipBuild  = 0;      /* =1 to skip initial build */
static int        skipTess   = 0;      /* -1 to skip tessellation at end of build */
static int        tessel     = 0;      /* =1 for tessellation sensitivities */
static int        verify     = 0;      /* =1 to enable verification */
static char       *filename  = NULL;   /* name of .csm file */
static char       *vrfyname  = NULL;   /* name of .vfy file */
static char       *despname  = NULL;   /* name of DESPMTRs file */
static char       *dictname  = NULL;   /* name of dictionary file */
static char       *dxddname  = NULL;   /* name of despmtr for -dxdd */
static char       *ptrbname  = NULL;   /* name of perturbation file */
static char       *eggname   = NULL;   /* name of external grid generator */
static char       *plotfile  = NULL;   /* name of plotdata file */
static char       *tessfile  = NULL;   /* name of tessellation file */
static char       *BDFname   = NULL;   /* name of BDF file to plot */

/* global variables associated with graphical user interface (gui) */
#define MAX_CLIENTS  100
static wvContext  *cntxt     = NULL;   /* context for the WebViewer */
static int         port      = 7681;   /* port number */
static int         serverNum = -1;     /* server number */

/* global variables associated with multiple users */
static void        *mutex    = NULL;   /* mutex for browserMessage */
static char        usernames[1024];    /* |-separated list of users */
static int         hasBall   = 0;      /* index of user with ball */

/* global variables associated with undo */
#define MAX_UNDOS  100
static int         nundo       = 0;    /* number of undos */
static void        **undo_modl = NULL;
static char        **undo_text = NULL;

/* global variables associated with scene graph meta-data */
#define MAX_METADATA_CHUNK 32000
static int         sgMetaDataSize = 0;
static int         sgMetaDataUsed = 0;
static char        *sgMetaData    = NULL;
static char        *sgFocusData   = NULL;

/* global variables associated with updated filelist */
static int        updatedFilelist = 1;
static char       *filelist       = NULL;

/* global variables associated with pending errors */
static int        pendingError =  0;
static int        successBuild = -1;

/* global variables associated with sensitivities */
static int        plotType = 0;
static float      lims[2]  = {-1.0, +1.0};
static int        haveDots = 0, sensPost = 0;
static double     sensLo = 0.0, sensHi = 0.0;
static char       *dotName=NULL;

/* global variables associated with Sketches */
static char       *skbuff = NULL;

/* global variables associated with StepThru mode */
static int        curStep = 0;
static float      sgFocus[4];

/* global variables associated with the response buffer */
static int        max_resp_len = 0;
static int        response_len = 0;
static char       *response = NULL;
static int        max_mesg_len = 0;
static int        messages_len = 0;
static char       *messages = NULL;

/* global variables associated with journals */
static FILE       *jrnl_out = NULL;    /* output journal file */

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

/* declarations for high-level routines defined below */

/* declarations for routines defined below */
static void       addToResponse(char text[]);
static void       addToSgMetaData(char format[], ...);
static int        applyDisplacement(modl_T *MODL, int ipmtr);
       void       browserMessage(void *udata, void *wsi, char text[], int lena);
static int        buildBodys(int buildTo, int *builtTo, int *buildStatus, int *nwarn);
static int        buildSceneGraph();
static int        buildSceneGraphBody(int ibody);
static void       cleanupMemory(int quiet);
static int        getToken(char text[], int nskip, char sep, char token[]);
static int        maxDistance(modl_T *MODL1, modl_T *MODL2, int ibody, double *dist);
       void       mesgCallbackFromOpenCSM(char mesg[]);
static int        processBrowserToServer(char text[]);
       void       sizeCallbackFromOpenCSM(void *modl, int ipmtr, int nrow, int ncel);
static void       spec_col(float scalar, float out[]);
static int        storeUndo(char cmd[], char arg[]);
static int        updateModl(modl_T *src_MODL, modl_T *tgt_MODL);

static int        addToHistogram(double entry, int nhist, double dhist[], int hist[]);
static int        printHistogram(int nhist, double dhist[], int hist[]);

static int        computeMassProps(modl_T *MODL, int ibody, double props[]);
static int        checkForGanged(modl_T *MODL);
static int        checkParallel(modl_T *MODL);

static int        plugsMain(modl_T *MODL, int npass, int ncloud, double cloud[]);
static int        plugsPhase1(modl_T *modl,            int ibody, int npmtr, int pmtrindx[], int ncloud, double cloud[], double *rmsbest, /*@null@*/FILE *fp_plugs);
static int        plugsPhase2(modl_T *modl, int npass, int ibody, int npmtr, int pmtrindx[], int ncloud, double cloud[], double *rmsbest, /*@null@*/FILE *fp_plugs);
static int        matsol(double A[], double b[], int n, double x[]);
static int        solsvd(double A[], double b[], int mrow, int ncol, double W[], double x[]);
static int        tridiag(int n, double a[], double b[], double c[], double d[], double x[]);
static int        writeSensFile(modl_T *modl, int ibody, char filename[]);


/***********************************************************************/
/*                                                                     */
/*   main program                                                      */
/*                                                                     */
/***********************************************************************/

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array  of arguments */
{
    int       imajor, iminor, status, i, bias, showUsage=0;
    int       builtTo, buildStatus, nwarn=0, plugs_save;
    int       npmtrs, type, nrow, irow, ncol, icol, ii, *ipmtrs=NULL, *irows=NULL, *icols=NULL;
    int       ipmtr, jpmtr, iundo; //, inode, iedge, iface;
    float     fov, zNear, zFar;
    float     eye[3]    = {0.0, 0.0, 7.0};
    float     center[3] = {0.0, 0.0, 0.0};
    float     up[3]     = {0.0, 1.0, 0.0};
    double    data[18], bbox[6], value, dot, *values=NULL, dist, size, err, maxerr=0;
    double    *cloud=NULL;
    char      *casename=NULL, *jrnlname=NULL, *tempname=NULL, *pmtrname=NULL;
    char      *dirname=NULL, *basename=NULL;
    char      pname[MAX_NAME_LEN], strval[MAX_STRVAL_LEN];
    char      *text=NULL, *esp_start;
    CCHAR     *OCC_ver;
    FILE      *jrnl_in=NULL, *ptrb_fp, *vrfy_fp, *auto_fp, *plot_fp;
    modl_T    *MODL=NULL;

    int       ibody;

    clock_t   old_time, new_time, old_totaltime, new_totaltime;

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    old_totaltime = clock();

    /* dynamically allocated array so that everything is on heap (not stack) */
    sgMetaDataSize = MAX_METADATA_CHUNK;

    MALLOC(dotName,     char, MAX_STRVAL_LEN  );
    MALLOC(sgMetaData,  char, sgMetaDataSize  );
    MALLOC(sgFocusData, char, MAX_STRVAL_LEN  );
    MALLOC(casename,    char, MAX_FILENAME_LEN);
    MALLOC(jrnlname,    char, MAX_FILENAME_LEN);
    MALLOC(tempname,    char, MAX_FILENAME_LEN);
    MALLOC(filename,    char, MAX_FILENAME_LEN);
    MALLOC(vrfyname,    char, MAX_FILENAME_LEN);
    MALLOC(despname,    char, MAX_FILENAME_LEN);
    MALLOC(dictname,    char, MAX_FILENAME_LEN);
    MALLOC(dxddname,    char, MAX_FILENAME_LEN);
    MALLOC(ptrbname,    char, MAX_FILENAME_LEN);
    MALLOC(pmtrname,    char, MAX_FILENAME_LEN);
    MALLOC(dirname,     char, MAX_FILENAME_LEN);
    MALLOC(basename,    char, MAX_FILENAME_LEN);
    MALLOC(eggname,     char, MAX_FILENAME_LEN);
    MALLOC(plotfile,    char, MAX_FILENAME_LEN);
    MALLOC(tessfile,    char, MAX_FILENAME_LEN);
    MALLOC(BDFname,     char, MAX_FILENAME_LEN);
    MALLOC(text,        char, MAX_STR_LEN     );

    MALLOC(undo_modl, void*, MAX_UNDOS+1);
    MALLOC(undo_text, char*, MAX_UNDOS+1);
    for (i = 0; i <= MAX_UNDOS; i++) {
        undo_text[i] = NULL;
    }
    for (i = 0; i <= MAX_UNDOS; i++) {
        MALLOC(undo_text[i], char, MAX_NAME_LEN);
    }

    strcpy(usernames, "|");

    /* get the flags and casename(s) from the command line */
    casename[ 0] = '\0';
    jrnlname[ 0] = '\0';
    filename[ 0] = '\0';
    vrfyname[ 0] = '\0';
    despname[ 0] = '\0';
    dictname[ 0] = '\0';
    dxddname[ 0] = '\0';
    ptrbname[ 0] = '\0';
    pmtrname[ 0] = '\0';
    eggname[  0] = '\0';
    plotfile[ 0] = '\0';
    tessfile[ 0] = '\0';
    BDFname[  0] = '\0';

    for (i = 1; i < argc; i++) {
        if        (strcmp(argv[i], "--") == 0) {
            /* ignore (needed for gdb) */
        } else if (strcmp(argv[i], "-addVerify") == 0) {
            addVerify = 1;
        } else if (strcmp(argv[i], "-allVels") == 0) {
            allVels = 1;
        } else if (strcmp(argv[i], "-batch") == 0) {
            batch = 1;
        } else if (strcmp(argv[i], "-checkMass") == 0) {
            checkMass = 1;
        } else if (strcmp(argv[i], "-checkPara") == 0) {
            checkPara = 1;
        } else if (strcmp(argv[i], "-despmtrs") == 0) {
            if (i < argc-1) {
                STRNCPY(despname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-dict") == 0) {
            if ( i < argc-1) {
                STRNCPY(dictname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-dumpEgads") == 0) {
            dumpEgads = 1;
        } else if (strcmp(argv[i], "-dxdd") == 0) {
            if (i < argc-1) {
                STRNCPY(dxddname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-egg") == 0) {
            if (i < argc-1) {
                STRNCPY(eggname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-histDist") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%lf", &histDist);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-help") == 0 ||
                   strcmp(argv[i], "-h"   ) == 0   ) {
            showUsage = 1;
            break;
        } else if (strcmp(argv[i], "-jrnl") == 0) {
            if (i < argc-1) {
                STRNCPY(jrnlname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-loadEgads") == 0) {
            loadEgads = 1;
        } else if (strcmp(argv[i], "-onormal") == 0) {
            onormal = 1;
        } else if (strcmp(argv[i], "-outLevel") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &outLevel);
                if (outLevel < 0) outLevel = 0;
                if (outLevel > 3) outLevel = 3;
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-plot") == 0) {
            if (i < argc-1) {
                STRNCPY(plotfile, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-plotBDF") == 0) {
            if (i < argc-1) {
                STRNCPY(BDFname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-plotCP") == 0) {
            plotCP = 1;
        } else if (strcmp(argv[i], "-plugs") == 0) {
            if (i < argc-1) {
                plugsOn = 1;
                sscanf(argv[++i], "%d", &plugs);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-port") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &port);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-printStack") == 0) {
            printStack = 1;
        } else if (strcmp(argv[i], "-ptrb") == 0) {
            if ( i < argc-1) {
                STRNCPY(ptrbname, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-tessel") == 0 || strcmp(argv[i], "-sensTess") == 0) {
            tessel = 1;
        } else if (strcmp(argv[i], "-skipBuild") == 0) {
            skipBuild = 1;
        } else if (strcmp(argv[i], "-skipTess") == 0) {
            skipTess = 1;
        } else if (strcmp(argv[i], "-tess") == 0) {
            if (i < argc-1) {
                STRNCPY(tessfile, argv[++i], MAX_FILENAME_LEN);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-verify") == 0) {
            verify = 1;
        } else if (strcmp(argv[i], "--version") == 0 ||
                   strcmp(argv[i], "-version" ) == 0 ||
                   strcmp(argv[i], "-v"       ) == 0   ) {
            (void) ocsmVersion(&imajor, &iminor);
            SPRINT2(0, "OpenCSM version: %2d.%02d", imajor, iminor);
            EG_revision(&imajor, &iminor, &OCC_ver);
            SPRINT3(0, "EGADS   version: %2d.%02d (with %s)", imajor, iminor, OCC_ver);
            exit(EXIT_SUCCESS);
        } else if (STRLEN(casename) == 0) {
            STRNCPY(casename, argv[i], MAX_FILENAME_LEN);
        } else {
            SPRINT1(0, "two casenames given (%s)", argv[i]);
            showUsage = 1;
            break;
        }
    }

    (void) ocsmVersion(&imajor, &iminor);

    if (showUsage) {
        SPRINT2(0, "serveCSM version %2d.%02d\n", imajor, iminor);
        SPRINT0(0, "proper usage: 'serveCSM [casename[.csm]] [options...]'");
        SPRINT0(0, "   where [options...] = -addVerify");
        SPRINT0(0, "                        -allVels");
        SPRINT0(0, "                        -batch");
        SPRINT0(0, "                        -checkMass");
        SPRINT0(0, "                        -checkPara");
        SPRINT0(0, "                        -despmtrs despname");
        SPRINT0(0, "                        -dict dictname");
        SPRINT0(0, "                        -dumpEgads");
        SPRINT0(0, "                        -dxdd despmtr");
        SPRINT0(0, "                        -egg eggname");
        SPRINT0(0, "                        -help  -or-  -h");
        SPRINT0(0, "                        -histDist dist");
        SPRINT0(0, "                        -jrnl jrnlname");
        SPRINT0(0, "                        -loadEgads");
        SPRINT0(0, "                        -onormal");
        SPRINT0(0, "                        -outLevel X");
        SPRINT0(0, "                        -plot plotfile");
        SPRINT0(0, "                        -plotBDF BDFname");
        SPRINT0(0, "                        -plotCP");
        SPRINT0(0, "                        -plugs npass");
        SPRINT0(0, "                        -port X");
        SPRINT0(0, "                        -printStack");
        SPRINT0(0, "                        -ptrb ptrbname");
        SPRINT0(0, "                        -skipBuild");
        SPRINT0(0, "                        -skipTess");
        SPRINT0(0, "                        -tess tessfile");
        SPRINT0(0, "                        -verify");
        SPRINT0(0, "                        -version  -or-  -v  -or-  --version");
        SPRINT0(0, "STOPPING...\a");
        status = -998;
        goto cleanup;
    }

    /* if you specify -dxdd or -skipTess, then batch is automatically enabled */
    if (strlen(dxddname) > 0 || skipTess == 1) {
        batch = 1;
    }

    /* welcome banner */
    SPRINT0(1, "**********************************************************");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*                    Program serveCSM                    *");
    SPRINT2(1, "*                     version %2d.%02d                      *", imajor, iminor);
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*        written by John Dannenhoffer, 2010/2022         *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "**********************************************************\n");

    SPRINT1(1, "    casename    = %s", casename   );
    SPRINT1(1, "    addVerify   = %d", addVerify  );
    SPRINT1(1, "    allVels     = %d", allVels    );
    SPRINT1(1, "    batch       = %d", batch      );
    SPRINT1(1, "    checkMass   = %d", checkMass  );
    SPRINT1(1, "    checkPara   = %d", checkPara  );
    SPRINT1(1, "    despmtrs    = %s", despname   );
    SPRINT1(1, "    dictname    = %s", dictname   );
    SPRINT1(1, "    dxddname    = %s", dxddname   );
    SPRINT1(1, "    dumpEgads   = %d", dumpEgads  );
    SPRINT1(1, "    eggname     = %s", eggname    );
    SPRINT1(1, "    jrnl        = %s", jrnlname   );
    SPRINT1(1, "    loadEgads   = %d", loadEgads  );
    SPRINT1(1, "    onormal     = %d", onormal    );
    SPRINT1(1, "    outLevel    = %d", outLevel   );
    SPRINT1(1, "    plotfile    = %s", plotfile   );
    SPRINT1(1, "    plotBDF     = %s", BDFname    );
    SPRINT1(1, "    plugsOn     = %d", plugsOn    );
    SPRINT1(1, "    plugs       = %d", plugs      );
    SPRINT1(1, "    port        = %d", port       );
    SPRINT1(1, "    printStack  = %d", printStack );
    SPRINT1(1, "    ptrbname    = %s", ptrbname   );
    SPRINT1(1, "    skipBuild   = %d", skipBuild  );
    SPRINT1(1, "    skipTess    = %d", skipTess   );
    SPRINT1(1, "    tessfile    = %s", tessfile   );
    SPRINT1(1, "    verify      = %d", verify     );
    SPRINT1(1, "    ESP_ROOT    = %s", /*@ignore@*/getenv("ESP_ROOT")/*@end@*/);
    SPRINT0(1, " ");

    plugs_save = plugs;

    /* set OCSMs output level */
    (void) ocsmSetOutLevel(outLevel);

    /* allocate nominal response and messages buffers */
    max_resp_len = 4096;
    MALLOC(response, char, max_resp_len+1);
    response[0] = '\0';

    max_mesg_len = 4096;
    MALLOC(messages, char, max_mesg_len+1);
    messages[0] = '\0';

    /* allocate sketch buffer */
    MALLOC(skbuff, char, MAX_STR_LEN);

    /* add .csm extension if a valid extension is not given */
    if (STRLEN(casename) > 0) {
        STRNCPY(filename, casename, MAX_FILENAME_LEN);
        if        (strstr(casename, ".csm"  ) != NULL ||
                   strstr(casename, ".cpc"  ) != NULL   ) {
            /* valid extension given */
        } else if (strstr(casename, ".stp"  ) != NULL ||
                   strstr(casename, ".step" ) != NULL ||
                   strstr(casename, ".STP"  ) != NULL ||
                   strstr(casename, ".STEP" ) != NULL   ) {
            auto_fp = fopen("autoStep.csm", "w");
            if (auto_fp != NULL) {
                strcpy(filename, casename);
                for (i = 0; i < strlen(filename); i++) {
                    if (filename[i] == '\\') filename[i] = '/';
                }
                fprintf(auto_fp, "# autoStep.csm (automatically generated)\n");
                fprintf(auto_fp, "IMPORT  %s  -1\n", filename);
                fprintf(auto_fp, "END\n");
                fclose(auto_fp);

                strcpy(filename, "autoStep.csm");
                SPRINT1(0, "Generated \"%s\" input file", filename);
            }
        } else if (strstr(casename, ".igs"  ) != NULL ||
                   strstr(casename, ".iges" ) != NULL ||
                   strstr(casename, ".IGS"  ) != NULL ||
                   strstr(casename, ".IGES" ) != NULL   ) {
            auto_fp = fopen("autoIges.csm", "w");
            if (auto_fp != NULL) {
                strcpy(filename, casename);
                for (i = 0; i < strlen(filename); i++) {
                    if (filename[i] == '\\') filename[i] = '/';
                }
                fprintf(auto_fp, "# autoIges.csm (automatically generated)\n");
                fprintf(auto_fp, "IMPORT  %s  -1\n", filename);
                fprintf(auto_fp, "END\n");
                fclose(auto_fp);

                strcpy(filename, "autoIges.csm");
                SPRINT1(0, "Generated \"%s\" imput file", filename);
            }
        } else if (strstr(casename, ".egads") != NULL ||
                   strstr(casename, ".EGADS") != NULL   ) {
            auto_fp = fopen("autoEgads.csm", "w");
            if (auto_fp != NULL) {
                strcpy(filename, casename);
                for (i = 0; i < strlen(filename); i++) {
                    if (filename[i] == '\\') filename[i] = '/';
                }
                fprintf(auto_fp, "# autoEgads.csm (automatically generated)\n");
                fprintf(auto_fp, "IMPORT  %s  -1\n", filename);
                fprintf(auto_fp, "END\n");
                fclose(auto_fp);

                strcpy(filename, "autoEgads.csm");
                SPRINT1(0, "Generated \"%s\" input file", filename);
            }
        } else {
            /* append .csm extension */
            strcat(filename, ".csm");
        }
    } else {
        casename[0] = '\0';
        filename[0] = '\0';
    }

    /* create the verify filename */
    if (verify == 1 || addVerify == 1) {
        EG_revision(&imajor, &iminor, &OCC_ver);

        /* get basename and dirname */
        i = strlen(filename) - 1;
        while (i >= 0) {
            if (filename[i] == '/' || filename[i] == '\\') {
                i++;
                break;
            }
            i--;
        }
        if (i == -1) {
            strcpy(dirname, ".");
            strcpy(basename, filename);
        } else {
            strcpy(basename, &(filename[i]));
            strcpy(dirname, filename);
            dirname[i-1] = '\0';
        }

        /* remove .csm or .cpc extension */
        i = strlen(basename);
        basename[i-4] = '\0';

        /* create the full vrfyname */
        snprintf(vrfyname, MAX_FILENAME_LEN-1, "%s%cverify_%s%c%s.vfy", dirname, SLASH, &(OCC_ver[strlen(OCC_ver)-5]), SLASH, basename);
    } else {
        vrfyname[0] = '\0';
    }

    /* read the .csm file and create the MODL */
    old_time = clock();
    status   = ocsmLoad(filename, &modl);
    new_time = clock();
    SPRINT3(1, "--> ocsmLoad(%s) -> status=%d (%s)", filename, status, ocsmGetText(status));
    SPRINT1(1, "==> ocsmLoad CPUtime=%9.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    if (status < SUCCESS && batch == 1) {
        SPRINT0(0, "ERROR:: problem in ocsmLoad");
        status = -999;
        goto cleanup;
    } else if (status < SUCCESS) {
        SPRINT0(0, "ERROR:: problem in ocsmLoad\a");
        pendingError = 1;
    } else {
        MODL = (modl_T*)modl;
    }

    SPLINT_CHECK_FOR_NULL(modl);

    if (pendingError == 0) {
        status = ocsmLoadDict(modl, dictname);
        if (status < EGADS_SUCCESS) goto cleanup;
    }

    if (pendingError == 0) {
        status = ocsmRegMesgCB(modl, mesgCallbackFromOpenCSM);
        if (status < EGADS_SUCCESS) goto cleanup;

        status = ocsmRegSizeCB(modl, sizeCallbackFromOpenCSM);
        if (status < EGADS_SUCCESS) goto cleanup;
    }

    if (strlen(despname) > 0) {
        status = ocsmUpdateDespmtrs(modl, despname);
        if (status < EGADS_SUCCESS) goto cleanup;
    }

    if (pendingError == 0) {
        MODL = (modl_T *)modl;
        if (filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        updatedFilelist = 1;
    }

    /* if verify is on, add verification data from .vfy file to Branches */
    if (verify == 1 && pendingError == 0) {
        old_time = clock();
        status   = ocsmLoad(vrfyname, &modl);
        MODL = (modl_T*)modl;
        new_time = clock();
        SPRINT3(1, "--> ocsmLoad(%s) -> status=%d (%s)", vrfyname, status, ocsmGetText(status));
        SPRINT1(1, "==> ocsmLoad CPUtime=%9.3f sec",
                (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));
    }

    /* check that Branches are properly ordered */
    if (pendingError == 0) {
        old_time = clock();
        status   = ocsmCheck(MODL);
        new_time = clock();
        SPRINT2(1, "--> ocsmCheck() -> status=%d (%s)",
                status, ocsmGetText(status));
        SPRINT1(1, "==> ocsmCheck CPUtime=%10.3f sec",
                (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

        if (status < SUCCESS && batch == 1) {
            SPRINT0(0, "ERROR:: problem in ocsmCheck");
            status = -999;
            goto cleanup;
        } else if (status < SUCCESS) {
            SPRINT0(0, "ERROR:: problem in ocsmCheck\a");
            pendingError = 1;
        }
    }

    /* print out the global Attributes and Parameters */
    if (batch == 1 && pendingError == 0) {
        SPRINT0(1, "External Parameter(s):");
        if (outLevel > 0) {
            status = ocsmPrintPmtrs(MODL, "");
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmPrintPmtrs -> status=%d", status);
            }
        }

        SPRINT0(1, "Global Attribute(s):");
        if (outLevel > 0) {
            status = ocsmPrintAttrs(MODL, "");
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmPrintAttrs -> status=%d", status);
            }
        }
    }

    /* set the external grid generator */
    if (pendingError == 0) {
        status = ocsmSetEgg(MODL, eggname);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmSetEgg -> status=%d", status);
            status = -999;
            goto cleanup;
        }
    }

    /* open the output journal file */
    snprintf(tempname, MAX_FILENAME_LEN, "port%d.jrnl", port);
    jrnl_out = fopen(tempname, "w");

    /* initialize the scene graph meta data */
    if (batch == 0) {
        sgMetaData[ 0] = '\0';
        sgMetaDataUsed = 0;
        sgFocusData[0] = '\0';
    }

    /* create the WebViewer context */
    if (batch == 0) {
        bias  =  1;
        if (onormal == 0) {
            fov    = 30.0;
            zNear  =  1.0;
            zFar   = 10.0;
        } else {
            eye[2] = 200;
            fov    =   1;
            zNear  = 195;
            zFar   = 205;
        }
        cntxt = wv_createContext(bias, fov, zNear, zFar, eye, center, up);
        if (cntxt == NULL) {
            SPRINT0(0, "ERROR:: failed to create wvContext");
            status = -999;
            goto cleanup;
        }

        wv_setCallBack(cntxt, browserMessage);
    }

    /* build the Bodys from the MODL */
    if (pendingError == 0) {
        status = buildBodys(0, &builtTo, &buildStatus, &nwarn);

        if (status != SUCCESS || buildStatus != SUCCESS || builtTo < 0) {
            successBuild = -1;
        } else {
            successBuild = builtTo;
        }

        /* interactive mode */
        if (batch == 0) {

            /* uncaught signal */
            if (builtTo < 0) {
                SPRINT2(0, "build() detected \"%s\" at %s",
                        ocsmGetText(buildStatus), MODL->brch[1-builtTo].name);
                SPRINT0(0, "Configuration only built up to detected error\a");
                pendingError = -builtTo;

            /* error in ocsmBuild that does not cause a signal to be raised */
            } else if (buildStatus != SUCCESS) {
                SPRINT0(0, "build() detected an error that did not raise a signal");
                SPRINT0(0, "Configuration only built up to detected error\a");
                pendingError = 299;

            /* error in ocsmcheck */
            } else if (status != SUCCESS) {
                SPRINT2(0, "ERROR:: build() detected %d (%s)",
                        status, ocsmGetText(status));
            }

        /* batch mode */
        } else {

            /* uncaught signal */
            if (builtTo < 0) {
                status = -999;
                goto cleanup;

            /* error in ocsmBuild that does not cause a signal to be raised */
            } else if (buildStatus != SUCCESS) {
                SPRINT2(0, "ERROR:: build() detected %d (%s)",
                        buildStatus, ocsmGetText(buildStatus));
                status = -999;
                goto cleanup;

            /* error in ocsmCheck */
            } else if (buildStatus != SUCCESS) {
                goto cleanup;

            }
        }

        /* if the -dxdd option is set, process it now */
        if (strlen(dxddname) > 0) {
            int  nbody;
            char sensfilename[MAX_FILENAME_LEN+1], *beg, *mid, *end;

            strcpy(sensfilename, dxddname);
            strcat(sensfilename, ".sens");

            /* extract row and/or column is given */
            if (strlen(dxddname) > 0 && strchr(dxddname, '[') != NULL) {
                beg = strchr(dxddname, '[');
                mid = strchr(dxddname, ',');
                end = strchr(dxddname, ']');

                if (beg == NULL || mid == NULL || end == NULL) {
                    SPRINT0(0, "if -dxdd is given, dxddname must be of form \"name\" or \"name[irow,icol]\"\a");
                    SPRINT0(0, "STOPPING...\a");
                    return EXIT_FAILURE;
                }

                *end = '\0';
                icol = atoi(mid+1);

                *mid = '\0';
                irow = atoi(beg+1);

                *beg = '\0';
            } else {
                irow = 1;
                icol = 1;
            }

            /* set the velocity of the DESPMTR */
            ipmtr = -1;
            for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
                if (MODL->pmtr[jpmtr].type == OCSM_DESPMTR   &&
                    strcmp(MODL->pmtr[jpmtr].name, dxddname) == 0  ) {
                    ipmtr = jpmtr;
                    break;
                }
            }

            if (ipmtr < 0) {
                SPRINT1(0, "ERROR:: no DESPMTR named \"%s\" found", dxddname);
                status = -999;
                goto cleanup;
            } else if (irow < 1 || irow > MODL->pmtr[ipmtr].nrow) {
                SPRINT2(0, "ERROR:: irow=%d is not between 1 and %d\n", irow, MODL->pmtr[ipmtr].nrow);
                status = -999;
                goto cleanup;
            } else if (icol < 1 || icol > MODL->pmtr[ipmtr].ncol) {
                SPRINT2(0, "ERROR:: icol=%d is not between 1 and %d\n", icol, MODL->pmtr[ipmtr].ncol);
                status = -999;
                goto cleanup;
            }

            status = ocsmSetVelD(MODL, 0,     0,    0,    0.0);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmSetVelD(clear) -> status=%d\n", status);
                status = -999;
                goto cleanup;
            }

            status = ocsmSetVelD(MODL, ipmtr, irow, icol, 1.0);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmSetVelD(set) -> status=%d\n", status);
                status = -999;
                goto cleanup;
            }

            /* rebuild to compute the velocities */
            nbody = 0;
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmBuild -> status=%d\n", status);
                status = -999;
                goto cleanup;
            }

            /* write the .tess file and exit */
            status = writeSensFile(MODL, MODL->nbody, sensfilename);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: writeSensFile -> status=%d\n", status);
                status = -999;
                goto cleanup;
            }

            SPRINT1(0, "==> \"%s\" has been written", sensfilename);
            goto cleanup;
        }
    }

    /* if there is a perturbation file, read it and create
       the perturbed MODL */
    if (STRLEN(ptrbname) > 0 && pendingError == 0) {

        /* open the ptrb file */
        ptrb_fp = fopen(ptrbname, "r");
        if (ptrb_fp == NULL) {
            SPRINT1(0, "ERROR:: perturbation \"%s\" not found", ptrbname);
            status = -999;
            goto cleanup;
        }

        SPRINT1(0, "--> Opening perturbation \"%s\"", ptrbname);

        /* determine the size of the file and allocate necessary arrays */
        npmtrs = 0;
        while (1) {
            status = fscanf(ptrb_fp, "%s %d %d %lf\n", pmtrname, &irow, &icol, &value);
            if (status < 4) break;
            npmtrs++;
        }

        MALLOC(ipmtrs, int,    npmtrs);
        MALLOC(irows,  int,    npmtrs);
        MALLOC(icols,  int,    npmtrs);
        MALLOC(values, double, npmtrs);

        /* re-read file and set up arrays */
        rewind(ptrb_fp);

        for (ii = 0; ii < npmtrs; ii++) {
            status = fscanf(ptrb_fp, "%s %d %d %lf\n",
                            pmtrname, &(irows[ii]), &icols[ii], &values[ii]);
            if (status != 4) {
                SPRINT1(0, "ERROR:: wrong number of values in \"%s\"", ptrbname);
                status = -999;
                goto cleanup;
            }

            status = ocsmFindPmtr(MODL, pmtrname, OCSM_DESPMTR, irows[ii], icols[ii], &ipmtrs[ii]);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: ocsmFindPmtr(%s) detected %d (%s)",
                        pmtrname, status, ocsmGetText(status));
                status = -999;
                goto cleanup;
            }

            SPRINT4(0, "    %20s[%2d,%2d] = %12.6f", pmtrname, irows[ii], icols[ii], values[ii]);
        }

        fclose(ptrb_fp);

        /* create the perturbed MODL */
        status = ocsmPerturb(MODL, npmtrs, ipmtrs, irows, icols, values);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: ocsmPerturb detected %d (%s)",
                    status, ocsmGetText(status));
            status = -999;
            goto cleanup;
        }

        /* find the maximum difference between the tessellation points in the
           base and perturbed Bodys */
        status = maxDistance(MODL, MODL->perturb, MODL->nbody, &dist);
        CHECK_STATUS(maxDistance);
        SPRINT1(1, "    maxDistance=%e", dist);

        SPRINT0(0, "--> Closing perturbation");

        FREE(values);
        FREE(icols );
        FREE(irows );
        FREE(ipmtrs);
    }

    /* if there is a tessellation file, read it and overwrite the tessellation
       of the last Body */
    if (STRLEN(tessfile) > 0 && pendingError == 0) {

        status = ocsmUpdateTess(MODL, MODL->nbody, tessfile);
        if (status == SUCCESS) {
            SPRINT1(1, "--> tessellation updated using \"%s\"", tessfile);
        } else {
            SPRINT1(0, "ERROR:: error update tessellation using \"%s\"", tessfile);
        }

        if (batch == 0) {
            buildSceneGraph();
        }
    }

    /* if plugs is on, run it now */
    if (plugsOn > 0 && pendingError == 0) {
        int     npass, ncloud, icloud, jmax;
        char    templine[128];

        /* make sure that there is a plotfile */
        plot_fp = fopen(plotfile, "r");
        if (plot_fp == NULL) {
            SPRINT1(0, "ERROR:: plotfile \"%s\" does not exist", plotfile);
            goto cleanup;
        } else {
            SPRINT1(1, "Running PLUGS for points in \"%s\"", plotfile);
        }

        /* read the header */
        if (fscanf(plot_fp, "%d %d %s", &ncloud, &jmax, templine) != 3) {
            SPRINT0(0, "ERROR:: problem reading plotfile header");
            goto cleanup;
        }

        /* read the points and close file */
        MALLOC(cloud, double, 3*ncloud);

        for (icloud = 0; icloud < ncloud; icloud++) {
            fscanf(plot_fp, "%lf %lf %lf", &cloud[3*icloud], &cloud[3*icloud+1], &cloud[3*icloud+2]);
        }

        fclose(plot_fp);

        /* prune out points that are with XXX of each other */
        if (PLUGS_PRUNE > 0) {
            int    jcloud;

            for (icloud = 1; icloud < ncloud; icloud++) {
                for (jcloud = 0; jcloud < icloud; jcloud++) {
                    dist = (cloud[3*icloud  ] - cloud[3*jcloud  ]) * (cloud[3*icloud  ] - cloud[3*jcloud  ])
                         + (cloud[3*icloud+1] - cloud[3*jcloud+1]) * (cloud[3*icloud+1] - cloud[3*jcloud+1])
                         + (cloud[3*icloud+2] - cloud[3*jcloud+2]) * (cloud[3*icloud+2] - cloud[3*jcloud+2]);
                    if (dist < PLUGS_PRUNE) {
                        cloud[3*icloud  ] = cloud[3*ncloud-3];
                        cloud[3*icloud+1] = cloud[3*ncloud-2];
                        cloud[3*icloud+2] = cloud[3*ncloud-1];
                        icloud--;
                        ncloud--;
                        break;
                    }
                }
            }

            plot_fp = fopen("new.cloud", "w");
            SPLINT_CHECK_FOR_NULL(plot_fp);

            fprintf(plot_fp, "%5d    0 new_cloud\n", ncloud);
            for (icloud = 0; icloud < ncloud; icloud++) {
                fprintf(plot_fp, "%22.15e %22.15e %22.15e\n",
                        cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2]);
            }
            fprintf(plot_fp, "    0    0 end\n");
            fclose(plot_fp);

            exit(0);
        }

        npass = plugs;

        /* execute plugs */
        status = plugsMain(MODL, npass, ncloud, cloud);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: plugsMain detected %d (%s)",
                    status, ocsmGetText(status));
            status = -999;
            goto cleanup;
        }

        FREE(cloud);

        /* enable verification */
        MODL->verify = verify;
    }

    /* process the input journal file if jrnlname exists */
    if (STRLEN(jrnlname) > 0) {
        SPRINT1(0, "\n==> Opening input journal file \"%s\"\n", jrnlname);

        jrnl_in = fopen(jrnlname, "r");
        if (jrnl_in == NULL) {
            SPRINT0(0, "ERROR:: Journal file cannot be opened");
            status = -999;
            goto cleanup;
        } else {
            while (1) {
                if (fgets(text, MAX_STR_LEN, jrnl_in) == NULL) break;
                if (feof(jrnl_in)) break;

                /* if there is a \n, convert it to a \0 */
                text[MAX_STR_LEN-1] = '\0';

                for (i = 0; i < STRLEN(text); i++) {
                    if (text[i] == '\n') {
                        text[i] =  '\0';
                        break;
                    }
                }

                status = processBrowserToServer(text);
                MODL = (modl_T*)modl;

                if (status < SUCCESS) {
                    fclose(jrnl_in);
                    goto cleanup;
                }
            }

            fclose(jrnl_in);

            SPRINT0(0, "\n==> Closing input journal file\n");
        }
    }

    if (pendingError == 0 && plugs_save >= 0 && MODL->sigCode < SUCCESS) {
        SPRINT2(0, "ERROR:: build not completed because error %d (%s) was detected",
                MODL->sigCode, ocsmGetText(MODL->sigCode));
        status = MODL->sigCode;
        goto cleanup;
    }

    /* if discrete displacement surfaces are specified, apply them now */
    if (pendingError == 0) {
        ipmtr = -1;
        for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
            if (strcmp(MODL->pmtr[jpmtr].name, "dds_spec") == 0) {
                ipmtr = jpmtr;
                break;
            }
        }
        if (ipmtr > 0) {
            old_time = clock();
            status = applyDisplacement(MODL, ipmtr);
            new_time = clock();
            SPRINT3(0, "--> applyDisplacement(ipmtr=%d) -> status=%d (%s)",
                    ipmtr, status, ocsmGetText(status));
            SPRINT1(0, "==> applyDisplacement CPUtime=%10.3f sec",
                    (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));
        }
    }

    /* get the command to start the client (if any) */
    if (batch == 0) {
        esp_start = getenv("ESP_START");
    }

    /* start the server */
    if (batch == 0) {
        /* create a lock */
        mutex = EMP_LockCreate();
        if (mutex == NULL) {
            SPRINT0(0, "WARNING:: lock could not be created");
        }

        status = SUCCESS;
        serverNum = wv_startServer(port, NULL, NULL, NULL, 0, cntxt);
        if (serverNum == 0) {

            /* stay alive a long as we have a client */
            while (wv_statusServer(0)) {
                SLEEP(500);

                /* start the browser if the first time through this loop */
                if (status == SUCCESS) {
                    if (esp_start != NULL) {
                        status += system(esp_start);
                    }

                    status++;
                }
            }
        }

        /* destroy the lock */
        if (mutex != NULL) {
            EMP_LockDestroy(mutex);
        }
    }

    /* cleanup and exit */
    MODL = (modl_T*)modl;

    /* update the thread using the context */
    if (MODL->context != NULL) {
        status = EG_updateThread(MODL->context);
        CHECK_STATUS(EG_updateThread);
    }

    /* print mass properties for all Bodys on stack */
    SPRINT0(1, "Mass properties of Bodys on stack");
    SPRINT0(1, "ibody    volume       area;len      xcg          ycg          zcg            Ixx          Ixy          Ixz          Iyy          Iyz          Izz");
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        status = EG_getMassProperties(MODL->body[ibody].ebody, data);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_getMassProperties(%d) -> status=%d\n", ibody, status);
        }

        SPRINT12(1, "%5d %12.4e %12.4e  %12.4e %12.4e %12.4e   %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e", ibody,
                 data[ 0], data[ 1], data[ 2], data[ 3], data[ 4],
                 data[ 5], data[ 6], data[ 7], data[ 9], data[10], data[13]);

        /* check area and volume */
        if (checkMass == 1) {
            double   props[14];

            status = computeMassProps(MODL, ibody, props);
            CHECK_STATUS(computeMassProps);

            SPRINT11(1, "      %12.4e %12.4e  %12.4e %12.4e %12.4e   %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                    props[ 0], props[ 1], props[ 2], props[ 3], props[ 4],
                    props[ 5], props[ 6], props[ 7], props[ 9], props[10], props[13]);

            /* compute relative error, normlized to the size of the configuration */
            for (i = 0; i < 13; i++) {
                if        (MODL->body[ibody].botype == OCSM_WIRE_BODY) {
                    size = data[1];                    // length
                } else if (MODL->body[ibody].botype == OCSM_SHEET_BODY) {
                    size = sqrt(data[1]);              // sqrt(area)
                } else if (MODL->body[ibody].botype == OCSM_SOLID_BODY) {
                    size = pow(data[0], 0.33333);      // cuberoot(volume)
                } else {
                    size = 1;
                }

                err = fabs(props[i] - data[i]) / size;
                if (err > maxerr) maxerr = err;
            }
        }
    }

    if (checkMass == 1) {
        SPRINT1(1, "Maximum massprop error = %10.3e (generated by -checkMass option)", maxerr);
    }

    /* print values of any output Parameters */
    SPRINT0(1, "Output Parameters");
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        status = ocsmGetPmtr(MODL, ipmtr, &type, &nrow, &ncol, pname);
        CHECK_STATUS(ocsmGetPmtr);

        if (type == OCSM_OUTPMTR) {
            if (nrow == 0 && ncol == 0) {
                status = ocsmGetValuS(MODL, ipmtr, strval);
                CHECK_STATUS(ocsmGetValuS);

                SPRINT2(1, "    %-20s %s", pname, strval);
            } else if (nrow > 1 || ncol > 1) {
                SPRINT1(1, "    %-20s", pname);
                for (irow = 1; irow <= nrow; irow++) {
                    for (icol = 1; icol <= ncol; icol++) {
                        status = ocsmGetValu(MODL, ipmtr, irow, icol, &value, &dot);
                        CHECK_STATUS(ocsmGetValu);

                        SPRINT4(1, "               [%3d,%3d] %11.5f %11.5f", irow, icol, value, dot);
                    }
                }
            } else {
                status = ocsmGetValu(MODL, ipmtr, 1, 1, &value, &dot);
                CHECK_STATUS(ocsmGetValu);

                SPRINT3(1, "    %-20s %11.5f %11.5f", pname, value, dot);
            }
        }
    }

    /* special code to analyze the parallelizablility of the build */
    if (checkPara == 1) {
        status = checkParallel(MODL);
        CHECK_STATUS(checkParallel);
    }

    /* check for possible ganged Boolean operations */
    status = checkForGanged(MODL);
    CHECK_STATUS(checkForGanged);

    /* special code to automatically add solution verification
       (via assertions) for all Bodys on the stack.  note: this makes
       a copy of the "*.csm" file and calls it "*.csm_verify" */
    if (addVerify) {
        SPRINT1(0, "WARNING:: writing verification data to \"%s\"", vrfyname);

        vrfy_fp = fopen(vrfyname, "w");
        if (vrfy_fp == NULL) {
            SPRINT1(0, "ERROR:: \"%s\" could not be created", vrfyname);
            status = -999;
            goto cleanup;
        }

        fprintf(vrfy_fp, "#======================================#\n");
        fprintf(vrfy_fp, "# automatically generated verification #\n");
        fprintf(vrfy_fp, "# OpenCSM %2d.%02d      %s #\n", imajor, iminor, &(OCC_ver[strlen(OCC_ver)-17]));
        fprintf(vrfy_fp, "#======================================#\n");

        for (ibody = 1; ibody <= MODL->nbody; ibody++) {
            if (MODL->body[ibody].onstack != 1) continue;

            fprintf(vrfy_fp, "select    body %d\n", ibody);

            if        (MODL->body[ibody].botype == OCSM_NODE_BODY) {
                fprintf(vrfy_fp, "   assert  %8d      @itype       0  1\n", 0);
            } else if (MODL->body[ibody].botype == OCSM_WIRE_BODY) {
                fprintf(vrfy_fp, "   assert  %8d      @itype       0  1\n", 1);
            } else if (MODL->body[ibody].botype == OCSM_SHEET_BODY) {
                fprintf(vrfy_fp, "   assert  %8d      @itype       0  1\n", 2);
            } else if (MODL->body[ibody].botype == OCSM_SOLID_BODY) {
                fprintf(vrfy_fp, "   assert  %8d      @itype       0  1\n", 3);
            }

            status = EG_getBoundingBox(MODL->body[ibody].ebody, bbox);
            if (status != SUCCESS) {
                SPRINT2(0, "ERROR:: EG_getBoundingBox(%d) -> status=%d\n", ibody, status);
            }

            status = EG_getMassProperties(MODL->body[ibody].ebody, data);
            if (status != SUCCESS) {
                SPRINT2(0, "ERROR:: EG_getMassProperties(%d) -> status=%d\n", ibody, status);
            }

            if (MODL->body[ibody].botype == OCSM_SHEET_BODY ||
                MODL->body[ibody].botype == OCSM_SOLID_BODY   ) {
                if        (fabs(data[0]) > 0.001) {
                    fprintf(vrfy_fp, "   assert %15.7e  @volume  -.001  1\n", data[0]);
                } else if (fabs(data[0]) < 1e-10) {
                    fprintf(vrfy_fp, "   assert %15.7e  @volume  0.001  1\n", 0.0);
                } else {
                    fprintf(vrfy_fp, "   assert %15.7e  @volume  0.001  1\n", data[0]);
                }
                if        (fabs(data[1]) > 0.001) {
                    fprintf(vrfy_fp, "   assert %15.7e  @area    -.001  1\n", data[1]);
                } else if (fabs(data[1]) < 1e-10) {
                    fprintf(vrfy_fp, "   assert %15.7e  @area    0.001  1\n", 0.0);
                } else {
                    fprintf(vrfy_fp, "   assert %15.7e  @area    0.001  1\n", data[1]);
                }
            } else if (MODL->body[ibody].botype == OCSM_WIRE_BODY) {
                if        (fabs(data[1]) > 0.001) {
                    fprintf(vrfy_fp, "   assert %15.7e  @length  -.001  1\n", data[1]);
                } else if (fabs(data[1]) < 1e-10) {
                    fprintf(vrfy_fp, "   assert %15.7e  @length  0.001  1\n", 0.0);
                } else {
                    fprintf(vrfy_fp, "   assert %15.7e  @length  0.001  1\n", data[1]);
                }
            }

            if        (fabs(data[2]) < 1e-10) {
                fprintf(vrfy_fp, "   assert %15.7e  @xcg     0.001  1\n", 0.0);
            } else if (bbox[3]-bbox[0] < 0.001) {
                fprintf(vrfy_fp, "   assert %15.7e  @xcg     -.001  1\n", data[2]);
            } else {
                fprintf(vrfy_fp, "   assert %15.7e  @xcg    %15.7e  1\n", data[2], 0.001*(bbox[3]-bbox[0]));
            }

            if        (fabs(data[3]) < 1e-10) {
                fprintf(vrfy_fp, "   assert %15.7e  @ycg     0.001  1\n", 0.0);
            } else if (bbox[4]-bbox[1] < 0.001) {
                fprintf(vrfy_fp, "   assert %15.7e  @ycg     -.001  1\n", data[3]);
            } else {
                fprintf(vrfy_fp, "   assert %15.7e  @ycg    %15.7e  1\n", data[3], 0.001*(bbox[4]-bbox[1]));
            }

            if        (fabs(data[4]) < 1e-10) {
                fprintf(vrfy_fp, "   assert %15.7e  @zcg     0.001  1\n", 0.0);
            } else if (bbox[5]-bbox[2] < 0.001) {
                fprintf(vrfy_fp, "   assert %15.7e  @zcg     -.001  1\n", data[4]);
            } else {
                fprintf(vrfy_fp, "   assert %15.7e  @zcg    %15.7e  1\n", data[4], 0.001*(bbox[5]-bbox[2]));
            }

            fprintf(vrfy_fp, "   assert  %8d      @nnode       0  1\n", MODL->body[ibody].nnode);
            fprintf(vrfy_fp, "   assert  %8d      @nedge       0  1\n", MODL->body[ibody].nedge);
            fprintf(vrfy_fp, "   assert  %8d      @nface       0  1\n", MODL->body[ibody].nface);

            fprintf(vrfy_fp, "\n");
        }

        fprintf(vrfy_fp, "end\n");
        fclose(vrfy_fp);
    }

    /* generate a histogram of the distance of plot points to Brep */
    if (histDist > 0 && strlen(plotfile) == 0) {
        SPRINT0(0, "WARNING:: Cannot choose -histDist without -plot");
    } else if (histDist > 0) {
        int     imax, jmax, j, iface, atype, alen, count=0;
        int     ibest, jbest;
        CINT    *tempIlist;
        double  dtest, dbest=0, xbest=0, ybest=0, zbest=0, ubest=0, vbest=0;
        double  dworst, drms, dultim, xyz_in[3], uv_out[2], xyz_out[3];
        CDOUBLE *tempRlist;
        char    templine[128];
        CCHAR   *tempClist;
        FILE    *fp_plot, *fp_all, *fp_bad;

        /* initialize the histogram */
        int    nhist=28, hist[28];
        double dhist[] = {1e-8, 2e-8, 5e-8,
                          1e-7, 2e-7, 5e-7,
                          1e-6, 2e-6, 5e-6,
                          1e-5, 2e-5, 5e-5,
                          1e-4, 2e-4, 5e-4,
                          1e-3, 2e-3, 5e-3,
                          1e-2, 2e-2, 5e-2,
                          1e-1, 2e-1, 5e-1,
                          1e+0, 2e+0, 5e+0,
                          1e+1};
        for (i = 0; i < 28; i++) {
            hist[i] = 0;
        }

        /* put the bounding box info as an attribute on each Face */
        for (ibody = 1; ibody <= MODL->nbody; ibody++) {
            if (MODL->body[ibody].onstack != 1) continue;

            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                status = EG_getBoundingBox(MODL->body[ibody].face[iface].eface, bbox);
                CHECK_STATUS(EG_getBoundingBox);

                status = EG_attributeAdd(MODL->body[ibody].face[iface].eface,
                                         "..bbox..", ATTRREAL, 6, NULL, bbox, NULL);
                CHECK_STATUS(EG_attributeAdd);
            }
        }

        /* open the plotfile */
        fp_plot = fopen(plotfile, "r");
        if (fp_plot == NULL) {
            SPRINT1(0, "ERROR:: pntsfile \"%s\" does not exist", plotfile);
            goto cleanup;
        } else {
            SPRINT1(1, "Computing distances to \"%s\"", plotfile);
        }

        /* open the bad point file */
        fp_bad = fopen("bad.points", "w");
        if (fp_bad == NULL) {
            SPRINT0(0, "ERROR:: could not create \"bad.points\"");
            goto cleanup;
        }

        /* open the all point file */
        fp_all = fopen("all.points", "w");
        if (fp_all == NULL) {
            SPRINT0(0, "ERROR:: could not create \"all.points\"");
            goto cleanup;
        }

        /* process each point in the plotfile */
        old_time = clock();
        dultim = 0;
        while (1) {
            if (fscanf(fp_plot, "%d %d %s", &imax, &jmax, templine) != 3) break;
            if (jmax == 0) jmax = 1;

            SPRINT3x(1, "imax=%8d, jmax=%8d, %-32s", imax, jmax, templine); fflush(NULL);

            dworst = 0;
            dbest = HUGEQ;
            drms  = 0;
            xbest = 0;
            ybest = 0;
            zbest = 0;
            ubest = -10;
            vbest = -10;
            ibest = -1;
            jbest = -1;

            for (j = 0; j < jmax; j++) {
                for (i = 0; i < imax; i++) {

                    /* read the point */
                    fscanf(fp_plot, "%lf %lf %lf", &(xyz_in[0]), &(xyz_in[1]), &(xyz_in[2]));

                    /* first look at the best Faces from last time, since the new point
                       tends to be near the previous point (and we can take advantage
                       of the bounding box screening) */
                    if (ibest > 0) {
                        status = EG_invEvaluate(MODL->body[ibest].face[jbest].eface,xyz_in, uv_out, xyz_out);
                        if (status != EGADS_DEGEN) {      // this happens if we try to inv eval along a degenerate Edge
                            CHECK_STATUS(EG_invEvaluate);

                            dbest = sqrt(SQR(xyz_out[0]-xyz_in[0]) + SQR(xyz_out[1]-xyz_in[1]) + SQR(xyz_out[2]-xyz_in[2]));
                            xbest = xyz_out[0];
                            ybest = xyz_out[1];
                            zbest = xyz_out[2];
                            ubest = uv_out[ 0];
                        }
                    }

                    /* now look at all Faces */
                    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                        if (MODL->body[ibody].onstack != 1) continue;

                        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                            status = EG_attributeRet(MODL->body[ibody].face[iface].eface, "..bbox..",
                                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
                            CHECK_STATUS(EG_attributeRet);

                            if (xyz_in[0] > tempRlist[0]-dbest && xyz_in[0] < tempRlist[3]+dbest &&
                                xyz_in[1] > tempRlist[1]-dbest && xyz_in[1] < tempRlist[4]+dbest &&
                                xyz_in[2] > tempRlist[2]-dbest && xyz_in[2] < tempRlist[5]+dbest   ) {

                                status = EG_invEvaluate(MODL->body[ibody].face[iface].eface, xyz_in, uv_out, xyz_out);
                                if (status != EGADS_DEGEN) {      // this happens if we try to inv eval along a degenerate Edge
                                    CHECK_STATUS(EG_invEvaluate);

                                    dtest = sqrt(SQR(xyz_out[0]-xyz_in[0]) + SQR(xyz_out[1]-xyz_in[1]) + SQR(xyz_out[2]-xyz_in[2]));
                                    if (dtest < dbest) {
                                        dbest = dtest;
                                        xbest = xyz_out[0];
                                        ybest = xyz_out[1];
                                        zbest = xyz_out[2];
                                        ubest = uv_out[ 0];
                                        vbest = uv_out[ 1];
                                        ibest = ibody;
                                        jbest = iface;
                                    }
                                }
                            }
                        }
                    }
                    if (dbest > dworst) {
                        dworst = dbest;
                    }
                    drms += dbest;

                    fprintf(fp_all, "%20.12f %20.12f %20.12f %5d %5d %20.12f %20.12f %20.12f %12.3e\n",
                            xyz_in[0], xyz_in[1], xyz_in[2], ibest, jbest, xbest, ybest, zbest, dbest);

                    if (dbest > histDist) {
                        fprintf(fp_bad, "%5d%5d point_%d_%d_%d\n",   1, 0, count, ibest, jbest);
                        fprintf(fp_bad, "%20.12f %20.12f %20.12f\n", xyz_in[0], xyz_in[1], xyz_in[2]);
                        fprintf(fp_bad, "%5d%5d line_%d_%f_%f\n",    2, 1, count, ubest, vbest);
                        fprintf(fp_bad, "%20.12f %20.12f %20.12f\n", xyz_in[0], xyz_in[1], xyz_in[2]);
                        fprintf(fp_bad, "%20.12f %20.12f %20.12f\n", xbest,     ybest,     zbest    );
                        count++;
                    }

                    /* add to histogram */
                    status = addToHistogram(dbest, nhist, dhist, hist);
                    CHECK_STATUS(addToHistogram);
                }
            }
            SPRINT2(1, " dworst=%12.3e, drms=%12.3e", dworst, sqrt(drms/imax/jmax));
            if (dworst > dultim) {
                dultim = dworst;
            }
        }
        SPRINT1(1, "dultim=%12.3e", dultim);

        /* close the plotfile and points files */
        fclose(fp_plot);
        fclose(fp_bad );
        fclose(fp_all );

        /* print the histogram */
        new_time = clock();
        SPRINT1(0, "Distance of plot points from Bodys on stack\nCPUtime=%9.3f sec",
                (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));
        status = printHistogram(nhist, dhist, hist);
        CHECK_STATUS(printHistogram);
    }

    /* free up storage associated with OpenCSM and EGADS */
    cleanupMemory(0);

    /* free up stoage associated with GUI */
    wv_cleanupServers();

    /* free up undo storage */
    for (iundo = nundo-1; iundo >= 0; iundo--) {
        (void) ocsmFree(undo_modl[iundo]);
    }

    if (undo_text != NULL) {
        for (i = 0; i <= MAX_UNDOS; i++) {
            FREE(undo_text[i]);
        }
        FREE(undo_text);
    }

    FREE(undo_modl  );

    nundo = 0;

    /* summary data */
    new_totaltime = clock();

    SPRINT1(1, "    Total CPU time = %.3f sec",
            (double)(new_totaltime - old_totaltime) / (double)(CLOCKS_PER_SEC));

    if (strlen(vrfyname) > 0) {
        vrfy_fp = fopen(vrfyname, "r");
        if (vrfy_fp != NULL) {
            fclose(vrfy_fp);

            if (nwarn == 0) {
                SPRINT0(0, "==> serveCSM completed successfully");
            } else {
                SPRINT1(0, "==> serveCSM completed successfully with %d warnings", nwarn);
            }
        } else {
            if (nwarn == 0) {
                SPRINT0(0, "==> serveCSM completed successfully with no verification data");
            } else {
                SPRINT1(0, "==> serveCSM completed successfully with %d warnings and no verification data", nwarn);
            }
        }
    } else {
        if (nwarn == 0) {
            SPRINT0(0, "==> serveCSM completed successfully");
        } else {
            SPRINT1(0, "==> serveCSM completed successfully with %d warnings", nwarn);
        }
    }

    status = SUCCESS;

cleanup:
    FREE(filelist);

    if (jrnl_out != NULL) fclose(jrnl_out);
    jrnl_out = NULL;

    if (undo_modl != NULL) {
        for (iundo = nundo-1; iundo >= 0; iundo--) {
            (void) ocsmFree(undo_modl[iundo]);
        }
        FREE(undo_modl);
    }

    if (undo_text != NULL) {
        for (i = 0; i <= MAX_UNDOS; i++) {
            FREE(undo_text[i]);
        }
        FREE(undo_text);
    }
    FREE(undo_modl  );
    FREE(text       );
    FREE(BDFname    );
    FREE(tessfile   );
    FREE(plotfile   );
    FREE(eggname    );
    FREE(basename   );
    FREE(dirname    );
    FREE(pmtrname   );
    FREE(ptrbname   );
    FREE(dictname   );
    FREE(dxddname   );
    FREE(vrfyname   );
    FREE(despname   );
    FREE(filename   );
    FREE(tempname   );
    FREE(jrnlname   );
    FREE(casename   );
    FREE(sgFocusData);
    FREE(sgMetaData );
    FREE(dotName    );

    FREE(values);
    FREE(icols );
    FREE(irows );
    FREE(ipmtrs);

    FREE(cloud);

    FREE(response);
    FREE(messages);
    FREE(skbuff  );

    if (status == -998) {
        return(EXIT_FAILURE);
    } else if (status < 0) {
        cleanupMemory(1);
        return(EXIT_FAILURE);
    } else {
        return(EXIT_SUCCESS);
    }
}


/***********************************************************************/
/*                                                                     */
/*   addToResponse - add text to response (with buffer length protection) */
/*                                                                     */
/***********************************************************************/

static void
addToResponse(char   text[])            /* (in)  text to add */
{
    int    status=SUCCESS, text_len;

    ROUTINE(addToResponse);

    /* --------------------------------------------------------------- */

    text_len = STRLEN(text);

    while (response_len+text_len > max_resp_len-2) {
        max_resp_len += 4096;
        SPRINT1(2, "increasing max_resp_len=%d", max_resp_len);

        RALLOC(response, char, max_resp_len+1);
    }

    strcpy(&(response[response_len]), text);

    response_len += text_len;

cleanup:
    if (status != SUCCESS) {
        SPRINT0(0, "ERROR:: max_resp_len could not be increased");
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   addToSgMetaData - add text to sgMetaData (with buffer length protection) */
/*                                                                     */
/***********************************************************************/

static void
addToSgMetaData(char   format[],        /* (in)  format specifier */
                 ...)                   /* (in)  variable arguments */
{
    int      status=SUCCESS, newchars;

    va_list  args;

    ROUTINE(addToSgMetaData);

    /* --------------------------------------------------------------- */

    /* make sure sgMetaData is big enough */
    if (sgMetaDataUsed+1024 >= sgMetaDataSize) {
        sgMetaDataSize += 1024 + MAX_METADATA_CHUNK;
        RALLOC(sgMetaData, char, sgMetaDataSize);
    }

    /* set up the va structure */
    va_start(args, format);

    newchars = vsnprintf(sgMetaData+sgMetaDataUsed, 1024, format, args);
    sgMetaDataUsed += newchars;

    /* clean up the va structure */
    va_end(args);

cleanup:
    if (status != SUCCESS) {
        SPRINT0(0, "ERROR:: sgMetaData could not be increased");
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   applyDisplacement - apply discrete displacement surface           */
/*                                                                     */
/***********************************************************************/

static int
applyDisplacement(modl_T *MODL,
                  int    ipmtr)
{
    int       status = SUCCESS;         /* return status */

    double    fact, limsrc[4], limtgt[4], *xyz_new=NULL, *uv_new=NULL;
    double    value, dot, uv_in[2], xyz_out[18], anorm[4];
    CDOUBLE   *xyz, *uv;
    int       isrc, itgtb, itgtf, periodic, npnt, ipnt, ntri, itri;
    int       type, nrow, ncol, irow;
    int       npnt_new, ntri_new, *tris_new=NULL, state;
    CINT      *ptype, *pindx, *tris, *tric;
    char      name[MAX_NAME_LEN];
    ego       esrc, etgt, etess, ebody;

    ROUTINE(applyDisplacement);

    /* --------------------------------------------------------------- */

    /* get the dds_spec design parameter */
    status = ocsmGetPmtr(MODL, ipmtr, &type, &nrow, &ncol, name);
    CHECK_STATUS(ocsmGetPmtr);

    if (type != OCSM_DESPMTR) {
        SPRINT0(0, "ERROR:: dds_spec is not an EXTERNAL parameter");
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    } else if (ncol != 4) {
        SPRINT1(0, "ERROR:: ncol=%d (and not 4)", ncol);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* loop through the various specifications */
    for (irow = 0; irow < nrow; irow++) {
        status = ocsmGetValu(MODL, ipmtr, irow+1, 1, &value, &dot);
        CHECK_STATUS(ocsmGetValu);
        isrc = NINT(value);
        if (isrc < 1 || isrc > MODL->nbody) break;

        status = ocsmGetValu(MODL, ipmtr, irow+1, 2, &value, &dot);
        CHECK_STATUS(ocsmGetValu);
        itgtb = NINT(value);
        if (itgtb < 1 || itgtb > MODL->nbody) break;

        status = ocsmGetValu(MODL, ipmtr, irow+1, 3, &value, &dot);
        CHECK_STATUS(ocsmGetValu);
        itgtf = NINT(value);
        if (itgtf < 1 || itgtf > MODL->body[itgtb].nface) break;

        status = ocsmGetValu(MODL, ipmtr, irow+1, 4, &value, &dot);
        CHECK_STATUS(ocsmGetValu);
        fact = value;
        if (fact == 0) break;

        SPRINT4(1, "    displacing itgt=%d:%d with isrc=%d with fact=%f",
                itgtb, itgtf, isrc, fact);

        /* get the ego for the source Face and its limits */
        esrc = MODL->body[isrc].face[1].eface;

        status = EG_getRange(esrc, limsrc, &periodic);
        CHECK_STATUS(EG_getRange);

        /* get the ego for the target Face and its limits */
        etgt = MODL->body[itgtb].face[itgtf].eface;

        status = EG_getRange(etgt, limtgt, &periodic);
        CHECK_STATUS(EG_getRange);

        /* get the tessellation for the target Face */
        etess = MODL->body[itgtb].etess;

        status = EG_getTessFace(etess, itgtf, &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        npnt_new = npnt;
        MALLOC(xyz_new, double, 3*npnt_new);
        MALLOC(uv_new,  double, 2*npnt_new);

        for (ipnt = 0;  ipnt < npnt_new; ipnt++) {
            xyz_new[3*ipnt  ] = xyz[3*ipnt  ];
            xyz_new[3*ipnt+1] = xyz[3*ipnt+1];
            xyz_new[3*ipnt+2] = xyz[3*ipnt+2];
            uv_new[ 2*ipnt  ] = uv[ 2*ipnt  ];
            uv_new[ 2*ipnt+1] = uv[ 2*ipnt+1];
        }

        ntri_new =  ntri;
        MALLOC(tris_new, int, 3*ntri_new);

        for (itri = 0; itri < ntri_new; itri++) {
            tris_new[3*itri  ] = tris[3*itri  ];
            tris_new[3*itri+1] = tris[3*itri+1];
            tris_new[3*itri+2] = tris[3*itri+2];
        }

        /* uodate the interior points */
        for (ipnt = 0; ipnt < npnt_new; ipnt++) {
            if (ptype[ipnt] < 0) {
                /* normal in target */
                status = EG_evaluate(etgt, &(uv_new[2*ipnt]), xyz_out);
                CHECK_STATUS(EG_evaluate);

                anorm[0] = xyz_out[4] * xyz_out[8] - xyz_out[5] * xyz_out[7];
                anorm[1] = xyz_out[5] * xyz_out[6] - xyz_out[3] * xyz_out[8];
                anorm[2] = xyz_out[3] * xyz_out[7] - xyz_out[4] * xyz_out[6];
                anorm[3] = sqrt(anorm[0] * anorm[0] + anorm[1] * anorm[1] + anorm[2] * anorm[2]);

            /* evaluate z in source */
                uv_in[0] = limsrc[0] + (limsrc[1] - limsrc[0]) * (uv_new[2*ipnt  ] - limtgt[0])
                                     / (limtgt[1] - limtgt[0]);
                uv_in[1] = limsrc[2] + (limsrc[3] - limsrc[2]) * (uv_new[2*ipnt+1] - limtgt[2])
                                     / (limtgt[3] - limtgt[2]);

                status = EG_evaluate(esrc, uv_in, xyz_out);
                CHECK_STATUS(EG_evaluate);

                /* apply the displacement */
                xyz_new[3*ipnt  ] += fact * anorm[0] / anorm[3] * xyz_out[2];
                xyz_new[3*ipnt+1] += fact * anorm[1] / anorm[3] * xyz_out[2];
                xyz_new[3*ipnt+2] += fact * anorm[2] / anorm[3] * xyz_out[2];
            }
        }

        /* open the tessellation for editing */
        status = EG_openTessBody(etess);
        CHECK_STATUS(EG_openTessBody);

        /* update the tessellation */
        status = EG_setTessFace(etess, itgtf,
                                npnt_new, xyz_new, uv_new,
                                ntri_new, tris_new);
        CHECK_STATUS(EG_setTessFace);

        /* cheking the status of the etess closes it */
        status = EG_statusTessBody(etess, &ebody, &state, &npnt);
        CHECK_STATUS(EG_statusTessBody);

        FREE( xyz_new);
        FREE(  uv_new);
        FREE(tris_new);
    }

    /* update the scene graph */
    if (batch == 0) {
        buildSceneGraph();
    }

cleanup:
    FREE( xyz_new);
    FREE(  uv_new);
    FREE(tris_new);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   browserMessage - called when client sends a message to the server */
/*                                                                     */
/***********************************************************************/

void
browserMessage(
   /*@unused@*/void    *udata,
   /*@unused@*/void    *wsi,
               char    text[],
  /*@unused@*/ int     lena)
{
    int       status;

    int       sendKeyData, ibody, onstack;
    char      message2[MAX_LINE_LEN];

    modl_T    *MODL = (modl_T*)modl;

    ROUTINE(browserMessage);

    /* --------------------------------------------------------------- */

    /* update the thread using the context */
    if (MODL->context != NULL) {
        status = EG_updateThread(MODL->context);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: EG_updateThread -> status=%d", status);
        }
    }

    /* set (or wait and set) lock so that messages do not overlap */
    if (mutex != NULL) {
        EMP_LockSet(mutex);
    }

    /* process the Message */
    (void) processBrowserToServer(text);

    /* send the response */
    SPRINT1(2, "\n<<< server2browser: %s", response);
//    wv_sendText(wsi, response);
    wv_broadcastText(response);

    /* if the sensitivities were just computed, send a message to
       inform the user about the lowest and highest sensitivity values */
    if (sensPost > 0) {
        snprintf(message2, MAX_LINE_LEN-1, "Sensitivities are in the range between %f and %f", sensLo, sensHi);
//        wv_sendText(wsi, message2);
        wv_broadcastText(message2);

        sensPost = 0;
    }

    sendKeyData = 0;

    /* send filenames if they have been updated */
    if (updatedFilelist == 1) {
        MODL = (modl_T *) modl;
        if (filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        SPLINT_CHECK_FOR_NULL(filelist);

        snprintf(message2, MAX_LINE_LEN-1, "getFilenames|%s", filelist);
        SPRINT1(2, "\n<<< server2browser: getFilenames|%s", filelist);
        wv_broadcastText(message2);

        updatedFilelist = 0;
    }

    /* send the scene graph meta data if it has not already been sent */
    if (STRLEN(sgMetaData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgData: %s", sgMetaData);
        wv_broadcastText(sgMetaData);

        /* nullify meta data so that it does not get sent again */
        sgMetaData[0]  = '\0';
        sgMetaDataUsed = 0;
        sendKeyData    = 1;
    }

    if (STRLEN(sgFocusData) > 0) {
        SPRINT1(2, "\n<<< server2browser: sgFocus: %s", sgFocusData);
        wv_broadcastText(sgFocusData);

        sendKeyData    = 1;
    }

    /* either open or close the key */
    if (sendKeyData == 1) {
        if (haveDots >  1) {
            if (tessel == 0) {
                status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Geom: d(norm)/d(***)");
            } else {
                status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Tess: d(norm)/d(***)");
            }
            SPRINT0(2, "\n<<< server2browser: setWvKey|on|");
            wv_broadcastText("setWvKey|on|");
        } else if (haveDots == 1) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], dotName         );
            SPRINT0(2, "\n<<< server2browser: setWvKey|on|");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 1) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Normalized U");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 2) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Normalized V");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 3) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Minimum Curv");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 4) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Maximum Curv");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 5) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "Gaussian Curv");
            wv_broadcastText("setWvKey|on|");
        } else if (plotType == 6) {
            status = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], "normals");
            wv_broadcastText("setWvKey|on|");
        } else {
            status = wv_setKey(cntxt,   0, NULL,      lims[0], lims[1], NULL            );
            SPRINT0(2, "\n<<< server2browser: setWvKey|off|");
            wv_broadcastText("setWvKey|off|");
        }
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setKet -> status=%d", status);
        }
    }

    /* send an error message (and the messages buffer) if one is pending */
    if (pendingError > 0) {
        snprintf(response, max_resp_len, "%s|%s|",
                 MODL->sigMesg, messages);

        SPRINT1(2, "\n<<< server2browser: %s", response);
        wv_broadcastText(response);

        pendingError =  0;
        successBuild = -1;

    } else if (pendingError == -1) {
        snprintf(response, max_resp_len, "ERROR:: could not find Design Velocities; shown as zeros|%s|",
            messages);

        SPRINT1(2, "\n<<< server2browser: %s", response);
        wv_broadcastText(response);

        pendingError =  0;
        successBuild = -1;

    } else if (successBuild >= 0) {
        onstack = 0;
        for (ibody = 1; ibody <= MODL->nbody; ibody++) {
            onstack += MODL->body[ibody].onstack;
        }

        snprintf(response, max_resp_len, "build|%d|%d|%s|",
                 successBuild, onstack, messages);

        SPRINT1(2, "\n<< server2browser: %s", response);
        wv_broadcastText(response);

        pendingError =  0;
        successBuild = -1;
    }

    messages[0] = '\0';
    messages_len = 0;

cleanup:
    /* release the lock */
    if (mutex != NULL) {
        EMP_LockRelease(mutex);
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   buildBodys - build Bodys and update scene graph                   */
/*                                                                     */
/***********************************************************************/

static int
buildBodys(int     buildTo,             /* (in)  last Branch to execute */
           int     *builtTo,            /* (out) last Branch successfully executed */
           int     *buildStatus,        /* (out) status returned from build */
           int     *nwarn)              /* (out) number of warnings */
{
    int       status = SUCCESS;         /* return status */

    int       nbody;
    int       showUVXYZ=0;              /* set to 1 to see UVXYZ at all corners */
    modl_T    *MODL = (modl_T*)modl;

    int       status2, ipmtr;
    clock_t   old_time, new_time;

    ROUTINE(buildBodys);

    /* --------------------------------------------------------------- */

    /* default return */
    *builtTo     = 0;
    *buildStatus = SUCCESS;
    *nwarn       = 0;

    /* restart StepThru mode */
    curStep = 0;

    /* if there are no Branches, simply return */
    if (MODL == NULL) {
        SPRINT0(1, "--> no MODL, so skipping build");
    } else {
        /* check that Branches are properly ordered */
        old_time = clock();
        status   = ocsmCheck(MODL);
        new_time = clock();
        SPRINT2(1, "--> ocsmCheck() -> status=%d (%s)",
                status, ocsmGetText(status));
        SPRINT1(1, "==> ocsmCheck CPUtime=%10.3f sec",
                (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

        if (status < SUCCESS) {
            goto cleanup;
        }

        /* set the verification flag (but not for original build if Plugs is on) */
        if (plugsOn == 0) {
            MODL->verify = verify;
        }

        /* set the dumpEgads and loadEgads flags */
        MODL->dumpEgads = dumpEgads;
        MODL->loadEgads = loadEgads;

        /* set the printStack flag */
        MODL->printStack = printStack;

        /* set the skip tessellation flag */
        MODL->tessAtEnd = 1 - skipTess;

        /* set the build erep flag */
        if (plotType == 10) {
            MODL->erepAtEnd = 1;
        } else {
            MODL->erepAtEnd = 0;
        }

        /* build the Bodys */
        if (skipBuild == 1) {
            SPRINT0(1, "--> skipping initial build");
            skipBuild = 0;
        } else {
            messages[0] = '\0';
            messages_len = 0;

            nbody        = 0;
            old_time     = clock();
            *buildStatus = ocsmBuild(MODL, buildTo, builtTo, &nbody, NULL);
            new_time = clock();
            SPRINT5(1, "--> ocsmBuild(buildTo=%d) -> status=%d (%s), builtTo=%d, nbody=%d",
                    buildTo, *buildStatus, ocsmGetText(*buildStatus), *builtTo, nbody);
            SPRINT1(1, "==> ocsmBuild CPUtime=%10.3f sec",
                    (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

            *nwarn = 0;
            for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
                if (strcmp(MODL->pmtr[ipmtr].name, "@nwarn") == 0) {
                    *nwarn = NINT(MODL->pmtr[ipmtr].value[0]);
                    break;
                }
            }

            /* print out the Branches */
            if (outLevel > 0 && MODL->sigCode == 0) {
                status2 = ocsmPrintBrchs(MODL, "");
                if (status2 != SUCCESS) {
                    SPRINT1(0, "ERROR:: ocsmPrintBrchs -> status=%d", status2);
                }
            }

            /* print out the Bodys */
            if (outLevel > 0 && MODL->sigCode == 0) {
                status2 = ocsmPrintBodys(MODL, "");
                if (status2 != SUCCESS) {
                    SPRINT1(0, "ERROR:: ocsmPrintBodys -> status=%d", status2);
                }
            }

            /* print out UV and XYZ at each corner of each Face */
            if (showUVXYZ) {
                int    ibody, iedge, iface, periodic;
                double uvrange[4], uv[2], xyz[18];

                SPRINT0(1, "ibody iedge        T                           X            Y            Z");
                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    if (MODL->body[ibody].onstack == 1) {
                        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                            (void) EG_getRange(MODL->body[ibody].edge[iedge].eedge, uvrange, &periodic);

                            uv[0] = uvrange[0];
                            (void) EG_evaluate(MODL->body[ibody].edge[iedge].eedge, uv, xyz);
                            SPRINT6(1, "%5d %5d   %12.6f                %12.6f %12.6f %12.6f",
                                    ibody, iedge, uv[0], xyz[0], xyz[1], xyz[2]);

                            uv[0] = uvrange[1];
                            (void) EG_evaluate(MODL->body[ibody].edge[iedge].eedge, uv, xyz);
                            SPRINT6(1, "%5d %5d   %12.6f                %12.6f %12.6f %12.6f",
                                    ibody, iedge, uv[0], xyz[0], xyz[1], xyz[2]);
                        }
                    }
                }

                SPRINT0(1, "ibody iface        U            V              X            Y            Z");
                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    if (MODL->body[ibody].onstack == 1) {
                        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                            (void) EG_getRange(MODL->body[ibody].face[iface].eface, uvrange, &periodic);

                            uv[0] = uvrange[0];   uv[1] = uvrange[2];
                            (void) EG_evaluate(MODL->body[ibody].face[iface].eface, uv, xyz);
                            SPRINT7(1, "%5d %5d   %12.6f %12.6f   %12.6f %12.6f %12.6f",
                                    ibody, iface, uv[0], uv[1], xyz[0], xyz[1], xyz[2]);

                            uv[0] = uvrange[1];   uv[1] = uvrange[2];
                            (void) EG_evaluate(MODL->body[ibody].face[iface].eface, uv, xyz);
                            SPRINT7(1, "%5d %5d   %12.6f %12.6f   %12.6f %12.6f %12.6f",
                                    ibody, iface, uv[0], uv[1], xyz[0], xyz[1], xyz[2]);

                            uv[0] = uvrange[0];   uv[1] = uvrange[3];
                            (void) EG_evaluate(MODL->body[ibody].face[iface].eface, uv, xyz);
                            SPRINT7(1, "%5d %5d   %12.6f %12.6f   %12.6f %12.6f %12.6f",
                                    ibody, iface, uv[0], uv[1], xyz[0], xyz[1], xyz[2]);

                            uv[0] = uvrange[1];   uv[1] = uvrange[3];
                            (void) EG_evaluate(MODL->body[ibody].face[iface].eface, uv, xyz);
                            SPRINT7(1, "%5d %5d   %12.6f %12.6f   %12.6f %12.6f %12.6f",
                                    ibody, iface, uv[0], uv[1], xyz[0], xyz[1], xyz[2]);
                        }
                    }
                }
            }
        }

        if (batch == 1) {
            if (*buildStatus < SUCCESS) {
                SPRINT2(0, "ERROR:: build not completed because error %d (%s) was detected",
                        *buildStatus, ocsmGetText(*buildStatus));

                status = -999;
                goto cleanup;
            } else if (*buildStatus > SUCCESS) {
                SPRINT1(0, "ERROR:: build not completed because user-thrown signal %d was uncaught",
                        *buildStatus);

                status = -999;
                goto cleanup;
            }
        }
    }

    /* disable -loadEgads flag (so that this flag is only enabled
       for first build (at the most) */
    loadEgads = 0;

    // this is different from gmgwCSM
    if (batch == 0) {
        buildSceneGraph();
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   buildSceneGraph - make a scene graph for wv                       */
/*                                                                     */
/***********************************************************************/

static int
buildSceneGraph()
{
    int       status = SUCCESS;         /* return status */

    int       ibody, jbody, iface, iedge, inode, iattr, nattr, ipmtr, irc, newStyleQuads, atype, alen;
    int       npnt, ipnt, ntri, itri, igprim, nseg, i, j, k, ij, ij1, ij2, ngrid, ncrod, nctri3, ncquad4;
    int       imax, jmax, ibeg, iend, isw, ise, ine, inw;
    int       attrs, head[3], nitems, nnode, nedge, nface;
    int       oclass, mtype, nchild, *senses, *header;
    int       *segs=NULL, *ivrts=NULL;
    int        npatch2, ipatch, n1, n2, i1, i2, *Tris=NULL;
    CINT      *ptype, *pindx, *tris, *tric, *pvindex, *pbounds;
    float     color[18], *pcolors=NULL, *plotdata=NULL, *segments=NULL, *tuft=NULL;
    double    bigbox[6], box[6], size, size2=0, xyz_dum[6], *vel=NULL, velmag;
    double    uvlimits[4], ubar, vbar, rcurv, xplot, yplot, zplot, fplot;
    double    data[18], normx, normy, normz, norm, axis[18], *cp;
    CDOUBLE   *xyz, *uv, *t;
    char      gpname[MAX_STRVAL_LEN], bname[MAX_NAME_LEN], temp[MAX_FILENAME_LEN];
    char      text[1025], dum[81];
    FILE      *fp;
    ego       ebody, etess, eface, eedge, enode, eref, prev, next, *echilds, esurf;
    ego       *enodes=NULL, *eedges=NULL, *efaces=NULL;

    int       nlist, itype;
    CINT      *tempIlist;
#if MIXED_TESS
    CINT      *nquad=NULL;
#endif
    CDOUBLE   *tempRlist;
    CCHAR     *attrName, *tempClist;

    wvData    items[6];

    modl_T    *MODL = (modl_T*)modl;

    ROUTINE(buildSceneGraph);

    /* --------------------------------------------------------------- */

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    if (MODL == NULL) {
        goto cleanup;
    }

    /* find the values needed to adjust the vertices */
    bigbox[0] = bigbox[1] = bigbox[2] = +HUGEQ;
    bigbox[3] = bigbox[4] = bigbox[5] = -HUGEQ;

    /* use Bodys on stack */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        ebody = MODL->body[ibody].ebody;

        status = EG_getBoundingBox(ebody, box);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_getBoundingBox(%d) -> status=%d", ibody, status);
        }

        if (box[0] < bigbox[0]) bigbox[0] = box[0];
        if (box[1] < bigbox[1]) bigbox[1] = box[1];
        if (box[2] < bigbox[2]) bigbox[2] = box[2];
        if (box[3] > bigbox[3]) bigbox[3] = box[3];
        if (box[4] > bigbox[4]) bigbox[4] = box[4];
        if (box[5] > bigbox[5]) bigbox[5] = box[5];
    }

    /* use any plotdata if it exists */
    if (plotfile != NULL) {
        fp = fopen(plotfile, "r");
        if (fp != NULL) {
            while (1) {
                if (fscanf(fp, "%d %d %s", &imax, &jmax, temp) != 3) break;

                if        (imax > 0 && jmax ==  0) {
                    npnt = imax;
                } else if (imax > 0 && jmax ==  1) {
                    npnt = imax;
                } else if (imax > 0 && jmax == -1) {
                    npnt = 2 * imax;
                } else if (imax > 0 && jmax == -2) {
                    npnt = 3 * imax;
                } else if (imax > 0 && jmax == -3) {
                    npnt = 3 * imax;
                } else if (imax > 0 && jmax == -4) {
                    npnt = 4 * imax;
                } else if (imax > 0 && jmax >   0) {
                    npnt = imax * jmax;
                } else {
                    break;
                }

                for (i = 0; i < npnt; i++) {
                    if (jmax == -3 || jmax == -4) {
                        fscanf(fp, "%lf %lf %lf %lf", &xplot, &yplot, &zplot, &fplot);
                    } else {
                        fscanf(fp, "%lf %lf %lf",     &xplot, &yplot, &zplot        );
                    }

                    if (xplot < bigbox[0]) bigbox[0] = xplot;
                    if (yplot < bigbox[1]) bigbox[1] = yplot;
                    if (zplot < bigbox[2]) bigbox[2] = zplot;
                    if (xplot > bigbox[3]) bigbox[3] = xplot;
                    if (yplot > bigbox[4]) bigbox[4] = yplot;
                    if (zplot > bigbox[5]) bigbox[5] = zplot;
                }

                if (feof(fp)) break;
            }

            fclose(fp);
            fp = NULL;
        }
    }

                                    size = bigbox[3] - bigbox[0];
    if (size < bigbox[4]-bigbox[1]) size = bigbox[4] - bigbox[1];
    if (size < bigbox[5]-bigbox[2]) size = bigbox[5] - bigbox[2];

    sgFocus[0] = (bigbox[0] + bigbox[3]) / 2;
    sgFocus[1] = (bigbox[1] + bigbox[4]) / 2;
    sgFocus[2] = (bigbox[2] + bigbox[5]) / 2;
    sgFocus[3] = size;

    /* generate the scene graph focus data */
    snprintf(sgFocusData, MAX_STRVAL_LEN-1, "sgFocus|[%20.12e,%20.12e,%20.12e,%20.12e]",
             sgFocus[0], sgFocus[1], sgFocus[2], sgFocus[3]);

    /* keep track of minimum and maximum velocities */
    sensLo = +HUGEQ;
    sensHi = -HUGEQ;

    /* initialize the scene graph meta data */
    sgMetaData[ 0] = '\0';
    sgMetaDataUsed = 0;

    addToSgMetaData("sgData|{");

    /* loop through the Bodys */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        if (MODL->body[ibody].eebody == NULL) {
            ebody = MODL->body[ibody].ebody;

            EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
            EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
            EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        } else {
            ebody = MODL->body[ibody].eebody;

            EG_getBodyTopos(ebody, NULL, NODE,  &nnode, &enodes);
            EG_getBodyTopos(ebody, NULL, EEDGE, &nedge, &eedges);
            EG_getBodyTopos(ebody, NULL, EFACE, &nface, &efaces);
        }

        /* set up Body name */
        snprintf(bname, MAX_NAME_LEN-1, "Body %d", ibody);

        status = EG_attributeRet(ebody, "_name", &itype, &nlist,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS && itype == ATTRSTRING) {
            snprintf(bname, MAX_NAME_LEN-1, "%s", tempClist);
        }

        /* check for duplicate Body names */
        for (jbody = 1; jbody < ibody; jbody++) {
            if (MODL->body[jbody].onstack != 1) continue;

            status = EG_attributeRet(MODL->body[jbody].ebody, "_name", &itype, &nlist,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status == SUCCESS && itype == ATTRSTRING) {
                if (strcmp(tempClist, bname) == 0) {
                    SPRINT2(0, "WARNING:: duplicate Body name (%s) found; being changed to \"Body %d\"", bname, ibody);
                    snprintf(bname, MAX_NAME_LEN-1, "Body %d", ibody);
                }
            }
        }

        /* Body info for NodeBody, SheetBody or SolidBody */
        snprintf(gpname, MAX_STRVAL_LEN, "%s", bname);
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        /* add Node to meta data (if there is room) */
        if (nattr > 0) {
            addToSgMetaData("\"%s\":[", gpname);
        } else {
            addToSgMetaData("\"%s\":[\"body\",\"%d\"", gpname, ibody);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(ebody, iattr, &attrName, &itype, &nlist,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iattr, status);
            }

            if (itype == ATTRCSYS) continue;

            addToSgMetaData("\"%s\",\"", attrName);

            if        (itype == ATTRINT) {
                for (i = 0; i < nlist ; i++) {
                    addToSgMetaData(" %d", tempIlist[i]);
                }
            } else if (itype == ATTRREAL) {
                for (i = 0; i < nlist ; i++) {
                    addToSgMetaData(" %f", tempRlist[i]);
                }
            } else if (itype == ATTRSTRING) {
                addToSgMetaData(" %s ", tempClist);
            }

            addToSgMetaData("\",");
        }
        sgMetaDataUsed--;
        addToSgMetaData("],");

        if (MODL->body[ibody].eetess == NULL) {
            etess = MODL->body[ibody].etess;
        } else {
            etess = MODL->body[ibody].eetess;
        }

        /* determine if any of the external Parameters have a velocity */
        haveDots   = 0;
        dotName[0] = '\0';

        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
                for (irc = 0; irc < (MODL->pmtr[ipmtr].nrow)*(MODL->pmtr[ipmtr].ncol); irc++) {
                    if (MODL->pmtr[ipmtr].dot[irc] != 0) {
                        if (fabs(MODL->pmtr[ipmtr].dot[irc]-1) < EPS06) {
                            if (haveDots == 0) {
                                if (tessel == 0) {
                                    snprintf(dotName, MAX_STRVAL_LEN-1, "Geom: d(norm)/d(%s)", MODL->pmtr[ipmtr].name);
                                } else {
                                    snprintf(dotName, MAX_STRVAL_LEN-1, "Tess: d(norm)/d(%s)", MODL->pmtr[ipmtr].name);
                                }
                            } else {
                                if (tessel == 0) {
                                    snprintf(dotName, MAX_STRVAL_LEN-1, "Geom: d(norm)/d(***)");
                                } else {
                                    snprintf(dotName, MAX_STRVAL_LEN-1, "Tess: d(norm)/d(***)");
                                }
                            }
                            haveDots++;
                        } else if (MODL->pmtr[ipmtr].dot[irc] != 0) {
                            if (tessel == 0) {
                                snprintf(dotName, MAX_STRVAL_LEN-1, "Geom: d(norm)/d(***)");
                            } else {
                                snprintf(dotName, MAX_STRVAL_LEN-1, "Tess: d(norm)/d(***)");
                            }
                            haveDots++;
                        }
                    }
                }
            }
        }

        /* get bounding box info if non-zero plottype */
        if (plotType > 0) {
            box[0] = box[1] = box[2] = 0;
            box[3] = box[4] = box[5] = 1;
            status = EG_getBoundingBox(ebody, box);
            if (status != SUCCESS) {
                SPRINT2(0, "ERROR:: EG_getBoundingBox(%d) -> status=%d", ibody, status);
            }

            size2 = (box[3] - box[0]) * (box[3] - box[0])
                  + (box[4] - box[1]) * (box[4] - box[1])
                  + (box[5] - box[2]) * (box[5] - box[2]);
        }

        /* determine if new-style quadding */
        newStyleQuads = 0;
        status = EG_attributeRet(etess, ".tessType",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS) {
#if MIXED_TESS
            newStyleQuads = 1;
            status = EG_attributeRet(etess, ".mixed",
                                     &atype, &alen, &nquad, &tempRlist, &tempClist);
            if (status != SUCCESS) newStyleQuads = 0;
#else
            if (strcmp(tempClist, "Quad") == 0) {
                newStyleQuads = 1;
            }
#endif
        }

        /* loop through the Faces within the current Body */
        for (iface = 1; iface <= nface; iface++) {
            nitems = 0;

            status = EG_getQuads(etess, iface,
                                 &npnt, &xyz, &uv, &ptype, &pindx,
                                 &npatch2);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getQuads(%d,%d) -> status=%d", ibody, iface, status);
            }

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Face %d", bname, iface);
            if (haveDots >= 1 || plotType > 0) {
                attrs = WV_ON | WV_SHADING;
            } else {
                attrs = WV_ON | WV_ORIENTATION;
            }

            /* render the Quadrilaterals (if they exist) */
            if (npatch2 > 0) {
                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* loop through the patches and build up the triangle and
                   segment tables (bias-1) */
                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    ntri += 2 * (n1-1) * (n2-1);
                    nseg += n1 * (n2-1) + n2 * (n1-1);
                }

                MALLOC(Tris, int, 3*ntri);
                MALLOC(segs, int, 2*nseg);

                ntri = 0;
                nseg = 0;

                for (ipatch = 1; ipatch <= npatch2; ipatch++) {
                    status = EG_getPatch(etess, iface, ipatch, &n1, &n2, &pvindex, &pbounds);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_getPatch(%d,%d) -> status=%d\n", ibody, iface, status);
                    }

                    for (i2 = 1; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            Tris[3*ntri  ] = pvindex[(i1-1)+n1*(i2-1)];
                            Tris[3*ntri+1] = pvindex[(i1  )+n1*(i2-1)];
                            Tris[3*ntri+2] = pvindex[(i1  )+n1*(i2  )];
                            ntri++;

                            Tris[3*ntri  ] = pvindex[(i1  )+n1*(i2  )];
                            Tris[3*ntri+1] = pvindex[(i1-1)+n1*(i2  )];
                            Tris[3*ntri+2] = pvindex[(i1-1)+n1*(i2-1)];
                            ntri++;
                        }
                    }

                    for (i2 = 0; i2 < n2; i2++) {
                        for (i1 = 1; i1 < n1; i1++) {
                            segs[2*nseg  ] = pvindex[(i1-1)+n1*(i2)];
                            segs[2*nseg+1] = pvindex[(i1  )+n1*(i2)];
                            nseg++;
                        }
                    }

                    for (i1 = 0; i1 < n1; i1++) {
                        for (i2 = 1; i2 < n2; i2++) {
                            segs[2*nseg  ] = pvindex[(i1)+n1*(i2-1)];
                            segs[2*nseg+1] = pvindex[(i1)+n1*(i2  )];
                            nseg++;
                        }
                    }
                }

                /* two triangles for rendering each quadrilateral */
                status = wv_setData(WV_INT32, 3*ntri, (void*)Tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(Tris);

            /* render the Quadrilaterals if new-style quadding */
            } else if (newStyleQuads == 1) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* loop through the triangles and build up the segment table */
                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            nseg++;
                        }
                    }
                }

                MALLOC(segs, int, 2*nseg);

#if MIXED_TESS
                /* create segments between Triangles (bias-1) */
                nseg = 0;
                for (itri = 0; itri < ntri-2*nquad[iface-1]; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                            segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                            nseg++;
                        }
                    }
                }

                /* create segments between Triangles (but not within the pair) (bias-1) */
                for (; itri < ntri; itri++) {
                    if (tric[3*itri  ] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                    itri++;

                    if (tric[3*itri  ] < itri) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                }
#else
                /* create segments between Triangles (but not within the pair) (bias-1) */
                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    if (tric[3*itri  ] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri+2) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                    itri++;

                    if (tric[3*itri  ] < itri) {
                        segs[2*nseg  ] = tris[3*itri+1];
                        segs[2*nseg+1] = tris[3*itri+2];
                        nseg++;
                    }
                    if (tric[3*itri+1] < itri) {
                        segs[2*nseg  ] = tris[3*itri+2];
                        segs[2*nseg+1] = tris[3*itri  ];
                        nseg++;
                    }
                    if (tric[3*itri+2] < itri) {
                        segs[2*nseg  ] = tris[3*itri  ];
                        segs[2*nseg+1] = tris[3*itri+1];
                        nseg++;
                    }
                }
#endif

                /* triangles */
                status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

            /* render the Triangles (if quads do not exist) */
            } else {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                }

                /* skip if no Triangles either */
                if (ntri <= 0) continue;

                /* vertices */
                status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* loop through the triangles and build up the segment table (bias-1) */
                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            nseg++;
                        }
                    }
                }

                MALLOC(segs, int, 2*nseg);

                nseg = 0;
                for (itri = 0; itri < ntri; itri++) {
                    for (k = 0; k < 3; k++) {
                        if (tric[3*itri+k] < itri+1) {
                            segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                            segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                            nseg++;
                        }
                    }
                }

                /* triangles */
                status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;
            }

            /* smooth colors (sensitivities) */
            if (haveDots >= 1) {
                sensPost = 1;

                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                SPLINT_CHECK_FOR_NULL(efaces);
                status = EG_getInfo(efaces[iface-1],
                                    &oclass, &mtype, &eref, &prev, &next);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getInfo(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                if (tessel == 0) {
                    MALLOC(vel, double, 3*npnt);

                    if (MODL->body[ibody].eebody == NULL) {
                        status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, npnt, NULL, vel);
                    } else {
                        status = ocsmGetVel(MODL, ibody, OCSM_EFACE, iface, npnt, NULL, vel);
                    }
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: ocsmGetVel(%d,%d) -> status=%d", ibody, iface, status);
                        goto cleanup;
                    }
                } else {
//$$$
                    status = ocsmGetTessVel(MODL, ibody, OCSM_FACE, iface, (const double**)(&vel));
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: ocsmGetTessVel(%d,%d) -> status=%d", ibody, iface, status);
                        goto cleanup;
                    }
                }

                SPLINT_CHECK_FOR_NULL(vel);

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    /* correct for NaN */
                    if (vel[3*ipnt  ] != vel[3*ipnt  ] ||
                        vel[3*ipnt+1] != vel[3*ipnt+1] ||
                        vel[3*ipnt+2] != vel[3*ipnt+2]   ) {
                        SPRINT1(0, "WARNING:: vel[%d] = NaN (being changed to 0)", ipnt);
                        velmag = 0;

                    /* find signed velocity magnitude */
                    } else if (tessel == 0) {
                        SPLINT_CHECK_FOR_NULL(efaces);
                        status = EG_evaluate(efaces[iface-1], &(uv[2*ipnt]), data);
                        if (status != SUCCESS) {
                            SPRINT3(0, "ERROR:: EG_evaluate(%d,%d) -> srtatus=%d", ibody, iface, status);
                            goto cleanup;
                        }

                        normx  = data[4] * data[8] - data[5] * data[7];
                        normy  = data[5] * data[6] - data[3] * data[8];
                        normz  = data[3] * data[7] - data[4] * data[6];

                        velmag = mtype * ( vel[3*ipnt  ] * normx
                                          +vel[3*ipnt+1] * normy
                                          +vel[3*ipnt+2] * normz)
                               / sqrt(normx * normx + normy * normy + normz * normz);

                        if (velmag != velmag) {
                            SPRINT1(0, "WARNING:: vel[%d] = NaN (being changed to 0)", ipnt);
                            velmag = 0;
                        }

                    /* find unsigned velocity magnitude */
                    } else {
                        velmag = sqrt( vel[3*ipnt  ] * vel[3*ipnt  ]
                                     + vel[3*ipnt+1] * vel[3*ipnt+1]
                                     + vel[3*ipnt+2] * vel[3*ipnt+2]);
                    }

                    spec_col((float)(velmag), &(pcolors[3*ipnt]));

                    if (velmag < sensLo) sensLo = velmag;
                    if (velmag > sensHi) sensHi = velmag;
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

                if (tessel == 0) {
                    FREE(vel);
                }

            /* smooth colors (normalized U) */
            } else if (plotType == 1) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                SPLINT_CHECK_FOR_NULL(efaces);
                status = EG_getTopology(efaces[iface-1], &eref, &oclass, &mtype, uvlimits,
                                        &nchild, &echilds, &senses);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR::EG_getTopology(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    ubar = (uv[2*ipnt  ] - uvlimits[0]) / (uvlimits[1] - uvlimits[0]);
                    spec_col((float)(ubar), &(pcolors[3*ipnt]));
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

            /* smooth colors (normalized V) */
            } else if (plotType == 2) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                SPLINT_CHECK_FOR_NULL(efaces);
                status = EG_getTopology(efaces[iface-1], &eref, &oclass, &mtype, uvlimits,
                                        &nchild, &echilds, &senses);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTopology(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    vbar = (uv[2*ipnt+1] - uvlimits[2]) / (uvlimits[3] - uvlimits[2]);
                    spec_col((float)(vbar), &(pcolors[3*ipnt]));
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

            /* smooth colors (minimum Curv) */
            } else if (plotType == 3) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                SPLINT_CHECK_FOR_NULL(efaces);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    status = EG_curvature(efaces[iface-1], &(uv[2*ipnt]), data);
                    if (status != SUCCESS) {
                        rcurv = 0;
                    } else {
                        rcurv = MIN(data[0], data[4]) * sqrt(size2);
                    }
                    spec_col((float)(rcurv), &(pcolors[3*ipnt]));
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

            /* smooth colors (maximum Curv) */
            } else if (plotType == 4) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                SPLINT_CHECK_FOR_NULL(efaces);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    status = EG_curvature(efaces[iface-1], &(uv[2*ipnt]), data);
                    if (status != SUCCESS) {
                        rcurv = 0;
                    } else {
                        rcurv = MAX(data[0], data[4]) * sqrt(size2);
                    }
                    spec_col((float)(rcurv), &(pcolors[3*ipnt]));
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

            /* smooth colors (Gaussian Curv) */
            } else if (plotType == 5) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                MALLOC(pcolors, float, 3*npnt);

                SPLINT_CHECK_FOR_NULL(efaces);
                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    status = EG_curvature(efaces[iface-1], &(uv[2*ipnt]), data);
                    if (status != SUCCESS) {
                        rcurv = 0;
                    } else if (MIN(fabs(data[0]), fabs(data[4])) < 0.00001 * MAX(fabs(data[0]), fabs(data[4]))) {
                        rcurv = 0;
                    } else if (data[0]*data[4] > 0) {
                        rcurv = +pow(fabs(data[0]*data[4]*size2), 0.25);
                    } else if (data[0]*data[4] < 0) {
                        rcurv = -pow(fabs(data[0]*data[4]*size2), 0.25);
                    } else {
                        rcurv = 0;
                    }
                    spec_col((float)(rcurv), &(pcolors[3*ipnt]));
                }

                status = wv_setData(WV_REAL32, npnt, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                FREE(pcolors);

            /* constant triangle colors */
            } else {
                if (MODL->body[ibody].eebody == NULL) {
                    color[0] = RED(  MODL->body[ibody].face[iface].gratt.color);
                    color[1] = GREEN(MODL->body[ibody].face[iface].gratt.color);
                    color[2] = BLUE( MODL->body[ibody].face[iface].gratt.color);
                } else {
                    color[0] = 0.75;
                    color[1] = 0.75;
                    color[2] = 1.00;
                }
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;
            }

            /* triangle backface color */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].face[iface].gratt.bcolor);
                color[1] = GREEN(MODL->body[ibody].face[iface].gratt.bcolor);
                color[2] = BLUE( MODL->body[ibody].face[iface].gratt.bcolor);
            } else {
                color[0] = 0.50;
                color[1] = 0.50;
                color[2] = 0.50;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* segment indices */
            status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            FREE(segs);

            /* segment colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].face[iface].gratt.mcolor);
                color[1] = GREEN(MODL->body[ibody].face[iface].gratt.mcolor);
                color[2] = BLUE( MODL->body[ibody].face[iface].gratt.mcolor);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* if plotCP is set and the associated surface is a BSPLINE, plot
               the control polygon (note that the control polygon is not listed in
               ESP as part of the Body, but separately at the bottom) */
            if (plotCP == 1) {
                SPLINT_CHECK_FOR_NULL(efaces);
                status = EG_getTopology(efaces[iface-1], &esurf, &oclass, &mtype,
                                        data, &nchild, &echilds, &senses);

                if (status == SUCCESS) {
                    status = EG_getGeometry(esurf, &oclass, &mtype,
                                            &eref, &header, &cp);

                    if (status == SUCCESS && oclass == SURFACE && mtype == BSPLINE) {
                        nitems = 0;

                        snprintf(gpname, MAX_STRVAL_LEN-1, "PlotCP: %d:%d", ibody, iface);
                        attrs = WV_ON;

                        /* control points */
                        status = wv_setData(WV_REAL64, header[2]*header[5], &(cp[header[3]+header[6]]),
                                            WV_VERTICES, &(items[nitems]));
                        if (status != SUCCESS) {
                            SPRINT3(0, "ERROR:: wv_setdata(%d,%d) -> status=%d", ibody, iface, status);
                        }

                        wv_adjustVerts(&(items[nitems]), sgFocus);
                        nitems++;

                        MALLOC(segs, int, 4*header[2]*header[5]);
                        nseg = 0;

                        /* i=constant lines (bias-1) */
                        for (i = 0; i < header[2]; i++) {
                            for (j = 0; j < header[5]-1; j++) {
                                segs[2*nseg  ] = 1 + (i) + (j  ) * header[2];
                                segs[2*nseg+1] = 1 + (i) + (j+1) * header[2];
                                nseg++;
                            }
                        }

                        /* j=constant lines (bias-1) */
                        for (j = 0; j < header[5]; j++) {
                            for (i = 0; i < header[2]-1; i++) {
                                segs[2*nseg  ] = 1 + (i  ) + (j) * header[2];
                                segs[2*nseg+1] = 1 + (i+1) + (j) * header[2];
                                nseg++;
                            }
                        }

                        status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_INDICES, &(items[nitems]));
                        if (status != SUCCESS) {
                            SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                        }
                        nitems++;

                        FREE(segs);

                        /* color */
                        color[0] = 0;   color[1] = 0;   color[2] = 0;
                        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                        if (status != SUCCESS) {
                            SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                        }
                        nitems++;

                        /* make graphic primitive */
                        igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                        if (igprim < 0) {
                            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                        }
                    }
                }
            }

            /* if tessel is set and we are asking for smooth colors
               (sensitivities) is set, draw tufts, whose length can be changed
               by recomputing sensitivity with a different value for the velocity
               of the Design Parameter.  (note that the tufts cannot be toggled in ESP) */
            if (tessel == 1 && haveDots >= 1) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                nitems = 0;

                /* name and attributes */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: Face_%d:%d_tufts", ibody, iface);
                attrs = WV_ON;

//$$$
                status = ocsmGetTessVel(MODL, ibody, OCSM_FACE, iface, (const double**)(&vel));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: ocsmGetTessVel(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                /* create tufts */
                MALLOC(tuft, float, 6*npnt);

                SPLINT_CHECK_FOR_NULL(vel);

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    tuft[6*ipnt  ] = xyz[3*ipnt  ];
                    tuft[6*ipnt+1] = xyz[3*ipnt+1];
                    tuft[6*ipnt+2] = xyz[3*ipnt+2];
                    tuft[6*ipnt+3] = xyz[3*ipnt  ] + vel[3*ipnt  ];
                    tuft[6*ipnt+4] = xyz[3*ipnt+1] + vel[3*ipnt+1];
                    tuft[6*ipnt+5] = xyz[3*ipnt+2] + vel[3*ipnt+2];
                }

                status = wv_setData(WV_REAL32, 2*npnt, (void*)tuft, WV_VERTICES, &(items[nitems]));

                FREE(tuft);

                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* tuft color */
                color[0] = 0;   color[1] = 0;   color[2] = 1;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                /* make graphic primitive for tufts */
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                }
            }

            /* get Attributes for the Face */
            SPLINT_CHECK_FOR_NULL(efaces);
            eface  = efaces[iface-1];
            status = EG_attributeNum(eface, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iface, status);
            }

            /* add Face to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData("\"%s\":[",  gpname);
            } else {
                addToSgMetaData("\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eface, iattr, &attrName, &itype, &nlist,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iface, status);
                }

                if (itype == ATTRCSYS) continue;

                addToSgMetaData("\"%s\",\"", attrName);

                if        (itype == ATTRINT) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %d", tempIlist[i]);
                    }
                } else if (itype == ATTRREAL) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %f", tempRlist[i]);
                    }
                } else if (itype == ATTRSTRING) {
                    addToSgMetaData(" %s ", tempClist);
                }

                addToSgMetaData("\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData("],");

            /* surface Normals */
            if (plotType == 6) {
                status = EG_getTessFace(etess, iface,
                                        &npnt, &xyz, &uv, &ptype, &pindx,
                                        &ntri, &tris, &tric);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessFace(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                nitems = 0;

                /* name and attributes */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: Face_%d:%d_norms", ibody, iface);
                attrs = WV_ON;

                /* create tufts */
                MALLOC(tuft, float, 6*npnt);

                status = EG_getTopology(MODL->body[ibody].face[iface].eface, &eref, &oclass, &mtype,
                                        data, &nchild, &echilds, &senses);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTopology(%d,%d) -> status=%d", ibody, iface, status);
                    goto cleanup;
                }

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(uv[2*ipnt]), data);
                    if (status != SUCCESS) {
                        SPRINT3(0, "ERROR:: EG_evaluate(%d,%d) -> status=%d", ibody, iface, status);
                        goto cleanup;
                    }

                    normx = (data[4] * data[8] - data[5] * data[7]);
                    normy = (data[5] * data[6] - data[3] * data[8]);
                    normz = (data[6] * data[7] - data[4] * data[6]);
                    norm  = sqrt(normx*normx + normy*normy + normz*normz);

                    tuft[6*ipnt  ] = xyz[3*ipnt  ];
                    tuft[6*ipnt+1] = xyz[3*ipnt+1];
                    tuft[6*ipnt+2] = xyz[3*ipnt+2];
                    tuft[6*ipnt+3] = xyz[3*ipnt  ] + mtype * lims[1] * normx / norm;
                    tuft[6*ipnt+4] = xyz[3*ipnt+1] + mtype * lims[1] * normy / norm;
                    tuft[6*ipnt+5] = xyz[3*ipnt+2] + mtype * lims[1] * normz / norm;
                }

                status = wv_setData(WV_REAL32, 2*npnt, (void*)tuft, WV_VERTICES, &(items[nitems]));

                FREE(tuft);

                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* tuft color */
                color[0] = 0;   color[1] = 0;   color[2] = 1;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iface, status);
                }
                nitems++;

                /* make graphic primitive for tufts */
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                }

            }
        }

        /* loop through the Edges within the current Body */
        for (iedge = 1; iedge <= nedge; iedge++) {
            nitems = 0;

            /* skip edge if this is a NodeBody */
            if (MODL->body[ibody].botype == OCSM_NODE_BODY) continue;

            status = EG_getTessEdge(etess, iedge, &npnt, &xyz, &t);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTessEdge(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* get Edge velocities (for testing) */
            if (allVels == 1 && haveDots > 0) {
                MALLOC(vel, double, 3*npnt);

                if (MODL->body[ibody].eebody == NULL) {
                    status = ocsmGetVel(MODL, ibody, OCSM_EDGE, iedge, npnt, NULL, vel);
                } else {
                    status = ocsmGetVel(MODL, ibody, OCSM_EEDGE, iedge, npnt, NULL, vel);
                }
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: ocsmGetVel(ibody=%d, iedge=%d) -> status=%d", ibody, iedge, status);
                }

                FREE(vel);
            }

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Edge %d", bname, iedge);
            attrs = WV_ON;

            /* vertices */
            status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* segments (bias-1) */
            MALLOC(ivrts, int, 2*(npnt-1));

            for (ipnt = 0; ipnt < npnt-1; ipnt++) {
                ivrts[2*ipnt  ] = ipnt + 1;
                ivrts[2*ipnt+1] = ipnt + 2;
            }

            status = wv_setData(WV_INT32, 2*(npnt-1), (void*)ivrts, WV_INDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* line colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].edge[iedge].gratt.color);
                color[1] = GREEN(MODL->body[ibody].edge[iedge].gratt.color);
                color[2] = BLUE( MODL->body[ibody].edge[iedge].gratt.color);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* points */
            MALLOC(ivrts, int, npnt);

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ivrts[ipnt] = ipnt + 1;
            }

            status = wv_setData(WV_INT32, npnt, (void*)ivrts, WV_PINDICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            FREE(ivrts);

            /* point colors */
            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].edge[iedge].gratt.mcolor);
                color[1] = GREEN(MODL->body[ibody].edge[iedge].gratt.mcolor);
                color[2] = BLUE( MODL->body[ibody].edge[iedge].gratt.mcolor);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2 ]= 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 2 (does not work for ANGLE) */
                cntxt->gPrims[igprim].lWidth = 2.0;

                /* make point size 5 */
                cntxt->gPrims[igprim].pSize  = 5.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = npnt - 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_addArrowHeads(%d,%d) -> status=%d", ibody, iedge, status);
                }
            }

            SPLINT_CHECK_FOR_NULL(eedges);
            eedge  = eedges[iedge-1];
            status = EG_attributeNum(eedge, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* if tessel is set and we are asking for smooth colors
               (sensitivities) is set, draw tufts, whose length can be changed
               by recomputing sensitivity with a different value for the velocity
               of the Design Parameter.  (note that the tufts cannot be toggled in ESP) */
            if (tessel == 1 && haveDots >= 1) {
                status = EG_getTessEdge(etess, iedge,
                                        &npnt, &xyz, &uv);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_getTessEdge(%d,%d) -> status=%d", ibody, iedge, status);
                    goto cleanup;
                }

                nitems = 0;

                /* name and attributes */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: Edge_%d:%d_tufts", ibody, iedge);
                attrs = WV_ON;

//$$$
                status = ocsmGetTessVel(MODL, ibody, OCSM_EDGE, iedge, (const double**)(&vel));
                CHECK_STATUS(ocsmGetTessVel);

                SPLINT_CHECK_FOR_NULL(vel);

                /* create tufts */
                MALLOC(tuft, float, 6*npnt);

                for (ipnt = 0; ipnt < npnt; ipnt++) {
                    tuft[6*ipnt  ] = xyz[3*ipnt  ];
                    tuft[6*ipnt+1] = xyz[3*ipnt+1];
                    tuft[6*ipnt+2] = xyz[3*ipnt+2];
                    tuft[6*ipnt+3] = xyz[3*ipnt  ] + vel[3*ipnt  ];
                    tuft[6*ipnt+4] = xyz[3*ipnt+1] + vel[3*ipnt+1];
                    tuft[6*ipnt+5] = xyz[3*ipnt+2] + vel[3*ipnt+2];
                }

                status = wv_setData(WV_REAL32, 2*npnt, (void*)tuft, WV_VERTICES, &(items[nitems]));

                FREE(tuft);

                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* tuft color */
                color[0] = 1;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
                }
                nitems++;

                /* make graphic primitive for tufts */
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                }

            }

            /* add Edge to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData("\"%s\":[",  gpname);
            } else {
                addToSgMetaData("\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(eedge, iattr, &attrName, &itype, &nlist,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (itype == ATTRCSYS) continue;

                addToSgMetaData("\"%s\",\"", attrName);

                if        (itype == ATTRINT) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %d", tempIlist[i]);
                    }
                } else if (itype == ATTRREAL) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %f", tempRlist[i]);
                    }
                } else if (itype == ATTRSTRING) {
                    addToSgMetaData(" %s ", tempClist);
                }

                addToSgMetaData("\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData("],");
        }

        /* loop through the Nodes within the current Body */
        for (inode = 1; inode <= nnode; inode++) {
            nitems = 0;

            /* get Node velocities (for testing) */
            if (allVels == 1 && haveDots > 0) {
                npnt = 1;
                MALLOC(vel, double, 3*npnt);

                status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, npnt, NULL, vel);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: ocsmGetVel(ibody=%d, inode=%d) -> status=%d", ibody, inode, status);
                }

                FREE(vel);
            }

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Node %d", bname, inode);

            /* default for NodeBodys is to turn the Node on */
            if (MODL->body[ibody].botype == OCSM_NODE_BODY) {
                attrs = WV_ON;
            } else {
                attrs = 0;
            }

            SPLINT_CHECK_FOR_NULL(enodes);
            enode  = enodes[inode-1];

            /* two copies of vertex (to get around possible wv limitation) */
            status = EG_getTopology(enode, &eref, &oclass, &mtype,
                                    xyz_dum, &nchild, &echilds, &senses);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_getTopology(%d,%d) -> status=%d", ibody, inode, status);
            }

            xyz_dum[3] = xyz_dum[0];
            xyz_dum[4] = xyz_dum[1];
            xyz_dum[5] = xyz_dum[2];

            status = wv_setData(WV_REAL64, 2, (void*)xyz_dum, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* node color (black) */
            color[0] = 0;   color[1] = 0;   color[2] = 0;

            if (MODL->body[ibody].eebody == NULL) {
                color[0] = RED(  MODL->body[ibody].node[inode].gratt.color);
                color[1] = GREEN(MODL->body[ibody].node[inode].gratt.color);
                color[2] = BLUE( MODL->body[ibody].node[inode].gratt.color);
            } else {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: wv_setData(%d,%d) -> status=%d", ibody, iedge, status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].pSize = 6.0;
            }

            status = EG_attributeNum(enode, &nattr);
            if (status != SUCCESS) {
                SPRINT3(0, "ERROR:: EG_attributeNum(%d,%d) -> status=%d", ibody, iedge, status);
            }

            /* add Node to meta data (if there is room) */
            if (nattr > 0) {
                addToSgMetaData("\"%s\":[",  gpname);
            } else {
                addToSgMetaData("\"%s\":[]", gpname);
            }

            for (iattr = 1; iattr <= nattr; iattr++) {
                status = EG_attributeGet(enode, iattr, &attrName, &itype, &nlist,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status != SUCCESS) {
                    SPRINT3(0, "ERROR:: EG_attributeGet(%d,%d) -> status=%d", ibody, iedge, status);
                }

                if (itype == ATTRCSYS) continue;

                addToSgMetaData("\"%s\",\"", attrName);

                if        (itype == ATTRINT) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %d", tempIlist[i]);
                    }
                } else if (itype == ATTRREAL) {
                    for (i = 0; i < nlist ; i++) {
                        addToSgMetaData(" %f", tempRlist[i]);
                    }
                } else if (itype == ATTRSTRING) {
                    addToSgMetaData(" %s ", tempClist);
                }

                addToSgMetaData("\",");
            }
            sgMetaDataUsed--;
            addToSgMetaData("],");
        }

        /* loop through the Csystems associated with the current Body */
        status = EG_attributeNum(ebody, &nattr);
        if (status != SUCCESS) {
            SPRINT2(0, "ERROR:: EG_attributeNum(%d) -> status=%d", ibody, status);
        }

        for (iattr = 1; iattr <= nattr; iattr++) {
            nitems = 0;

            status = EG_attributeGet(ebody, iattr, &attrName, &itype, &nlist,
                                     &tempIlist, &tempRlist, &tempClist);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_attributeGet -> status=%d", status);
            }

            if (itype != ATTRCSYS) continue;

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "%s Csys %s", bname, attrName);
            attrs = WV_ON | WV_SHADING | WV_ORIENTATION;

            /* vertices */
            axis[ 0] = tempRlist[nlist  ];
            axis[ 1] = tempRlist[nlist+1];
            axis[ 2] = tempRlist[nlist+2];
            axis[ 3] = tempRlist[nlist  ] + tempRlist[nlist+ 3];
            axis[ 4] = tempRlist[nlist+1] + tempRlist[nlist+ 4];
            axis[ 5] = tempRlist[nlist+2] + tempRlist[nlist+ 5];
            axis[ 6] = tempRlist[nlist  ];
            axis[ 7] = tempRlist[nlist+1];
            axis[ 8] = tempRlist[nlist+2];
            axis[ 9] = tempRlist[nlist  ] + tempRlist[nlist+ 6];
            axis[10] = tempRlist[nlist+1] + tempRlist[nlist+ 7];
            axis[11] = tempRlist[nlist+2] + tempRlist[nlist+ 8];
            axis[12] = tempRlist[nlist  ];
            axis[13] = tempRlist[nlist+1];
            axis[14] = tempRlist[nlist+2];
            axis[15] = tempRlist[nlist  ] + tempRlist[nlist+ 9];
            axis[16] = tempRlist[nlist+1] + tempRlist[nlist+10];
            axis[17] = tempRlist[nlist+2] + tempRlist[nlist+11];

            status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
            }

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* line colors (multiple colors requires that WV_SHADING be set above) */
            color[ 0] = 1;   color[ 1] = 0;   color[ 2] = 0;
            color[ 3] = 1;   color[ 4] = 0;   color[ 5] = 0;
            color[ 6] = 0;   color[ 7] = 1;   color[ 8] = 0;
            color[ 9] = 0;   color[10] = 1;   color[11] = 0;
            color[12] = 0;   color[13] = 0;   color[14] = 1;
            color[15] = 0;   color[16] = 0;   color[17] = 1;
            status = wv_setData(WV_REAL32, 6, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                /* make line width 1 */
                cntxt->gPrims[igprim].lWidth = 1.0;

                /* add arrow heads (requires that WV_ORIENTATION be set above) */
                head[0] = 1;
                status = wv_addArrowHeads(cntxt, igprim, 0.10/sgFocus[3], 1, head);
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_addArrowHeads -> status=%d", status);
                }
            }

            /* add Csystem to meta data (if there is room) */
            addToSgMetaData("\"%s\":[],", gpname);
        }
    }

    /* axes */
    nitems = 0;

    /* name and attributes */
    snprintf(gpname, MAX_STRVAL_LEN-1, "Axes");
    attrs = 0;

    /* vertices */
    axis[ 0] = MIN(2*bigbox[0]-bigbox[3], 0);   axis[ 1] = 0;   axis[ 2] = 0;
    axis[ 3] = MAX(2*bigbox[3]-bigbox[0], 0);   axis[ 4] = 0;   axis[ 5] = 0;

    axis[ 6] = 0;   axis[ 7] = MIN(2*bigbox[1]-bigbox[4], 0);   axis[ 8] = 0;
    axis[ 9] = 0;   axis[10] = MAX(2*bigbox[4]-bigbox[1], 0);   axis[11] = 0;

    axis[12] = 0;   axis[13] = 0;   axis[14] = MIN(2*bigbox[2]-bigbox[5], 0);
    axis[15] = 0;   axis[16] = 0;   axis[17] = MAX(2*bigbox[5]-bigbox[2], 0);
    status = wv_setData(WV_REAL64, 6, (void*)axis, WV_VERTICES, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(axis) -> status=%d", status);
    }

    wv_adjustVerts(&(items[nitems]), sgFocus);
    nitems++;

    /* line color */
    color[0] = 0.7;   color[1] = 0.7;   color[2] = 0.7;
    status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
    if (status != SUCCESS) {
        SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
    }
    nitems++;

    /* make graphic primitive */
    igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
    if (igprim < 0) {
        SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
    } else {
        SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

        cntxt->gPrims[igprim].lWidth = 1.0;
    }

    /* extra plotting data */
    SPLINT_CHECK_FOR_NULL(plotfile);

    if (STRLEN(plotfile) > 0) {
        fp = fopen(plotfile, "r");
        if (fp == NULL) {
            SPRINT1(0, "ERROR:: plotfile \"%s\" does not exist", plotfile);
            goto cleanup;
        } else {
            SPRINT1(1, "Opening \"%s\"", plotfile);
        }

        /* read multiple data sets */
        while (1) {
            if (fscanf(fp, "%d %d %s", &imax, &jmax, temp) != 3) break;

            /* points (imax=npts, jmax=0) */
            if (imax > 0 && jmax == 0) {
                SPRINT2(1, "    plotting %d points (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotPoints: %.114s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 3*imax);

                for (i = 0; i < imax; i++) {
                    fscanf(fp, "%f %f %f", &plotdata[3*i  ],
                                           &plotdata[3*i+1],
                                           &plotdata[3*i+2]);
                }

                /* points in plotdata */
                status = wv_setData(WV_REAL32, imax, plotdata, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
                }

                FREE(plotdata);

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* point color */
                color[0] = 0;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }
                nitems++;

                /* make graphic primitive */
                attrs  = WV_ON;
                igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                } else {
                    SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                    cntxt->gPrims[igprim].pSize = 3.0;
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* polyline (imax=npnt, jmax=1) */
            } else if (imax > 1 && jmax == 1) {
                SPRINT2(1, "    plotting line with %d points (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: %.116s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 3*imax);

                for (i = 0; i < imax; i++) {
                    fscanf(fp, "%f %f %f", &plotdata[3*i  ],
                                           &plotdata[3*i+1],
                                           &plotdata[3*i+2]);
                }

                /* create the segments */
                MALLOC(segments, float, 6*(imax-1));
                nseg = 0;

                for (i = 0; i < imax-1; i++) {
                    segments[6*nseg  ] = plotdata[3*i  ];
                    segments[6*nseg+1] = plotdata[3*i+1];
                    segments[6*nseg+2] = plotdata[3*i+2];

                    segments[6*nseg+3] = plotdata[3*i+3];
                    segments[6*nseg+4] = plotdata[3*i+4];
                    segments[6*nseg+5] = plotdata[3*i+5];

                    nseg++;
                }

                FREE(plotdata);

                status = wv_setData(WV_REAL32, 2*nseg, (void*)segments, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(segments) -> status=%d", status);
                }

                FREE(segments);

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* line color */
                color[0] = 0;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }
                nitems++;

                /* make graphic primitive and set line width */
                attrs = WV_ON;
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                } else {
                    SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                    cntxt->gPrims[igprim].lWidth = 1.0;
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* many line segments (imax=nseg, jmax=-1) */
            } else if (imax > 0 && jmax == -1) {
                SPRINT2(1, "    plotting %d lines with 2 points each (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: %.116s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 6*imax);

                for (i = 0; i < imax; i++) {
                    fscanf(fp, "%f %f %f", &plotdata[6*i  ],
                                           &plotdata[6*i+1],
                                           &plotdata[6*i+2]);
                    fscanf(fp, "%f %f %f", &plotdata[6*i+3],
                                           &plotdata[6*i+4],
                                           &plotdata[6*i+5]);
                }

                status = wv_setData(WV_REAL32, 2*imax, (void*)plotdata, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
                }

                FREE(plotdata);

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* line color */
                color[0] = 0;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }
                nitems++;

                /* make graphic primitive and set line width */
                attrs = WV_ON;
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                } else {
                    SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                    cntxt->gPrims[igprim].lWidth = 1.0;
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* many triangles (imax=ntri, jmax=-2) */
            } else if (imax > 0 && jmax == -2) {
                SPRINT2(1, "    plotting %d triangles (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotTris: %.114s", temp);

                MALLOC(plotdata, float, 9*imax);

                for (ij = 0; ij < imax; ij++) {
                    fscanf(fp, "%f %f %f %f %f %f %f %f %f", &plotdata[9*ij  ],
                                                             &plotdata[9*ij+1],
                                                             &plotdata[9*ij+2],
                                                             &plotdata[9*ij+3],
                                                             &plotdata[9*ij+4],
                                                             &plotdata[9*ij+5],
                                                             &plotdata[9*ij+6],
                                                             &plotdata[9*ij+7],
                                                             &plotdata[9*ij+8]);
                }

                /* points in plotdata */
                status = wv_setData(WV_REAL32, 3*imax, plotdata, WV_VERTICES, &(items[nitems++]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
                }

                FREE(plotdata);

                wv_adjustVerts(&(items[nitems-1]), sgFocus);

                /* triangle colors */
                color[0] = 0;   color[1] = 1;   color[2] = 1;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems++]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }

                /* triangle back colors */
                color[0] = 0;   color[1] = 0.5;   color[2] = 0.5;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_BCOLOR, &(items[nitems++]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }

                /* triangle sides */
                MALLOC(segs, int, 6*imax);

                /* build up the sides (bias-1) */
                for (ij = 0; ij < imax; ij++) {
                    segs[6*ij  ] = 3 * ij + 1;
                    segs[6*ij+1] = 3 * ij + 2;
                    segs[6*ij+2] = 3 * ij + 2;
                    segs[6*ij+3] = 3 * ij + 3;
                    segs[6*ij+4] = 3 * ij + 3;
                    segs[6*ij+5] = 3 * ij + 1;
                }

                status = wv_setData(WV_INT32, 6*imax, (void*)segs, WV_LINDICES, &(items[nitems++]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(segs) -> status=%d", status);
                }

                FREE(segs);

                /* triangle side color */
                color[0] = 1;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[nitems++]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }

                /* make graphic primitive */
                attrs  = WV_ON | WV_LINES;
                igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
               }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* filled Triangles (imax=ntri, jmax=-3) */
            } else if (imax > 1 && jmax == -3) {
                SPRINT2(1, "   plotting %d filled triangles (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotTris: %.11s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 12*imax);

                for (i = 0; i < imax; i++) {
                    fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f",
                           &plotdata[9*i  ], &plotdata[9*i+1], &plotdata[9*i+2], &plotdata[9*imax+3*i  ],
                           &plotdata[9*i+3], &plotdata[9*i+4], &plotdata[9*i+5], &plotdata[9*imax+3*i+1],
                           &plotdata[9*i+6], &plotdata[9*i+7], &plotdata[9*i+8], &plotdata[9*imax+3*i+2]);
                }

                status = wv_setData(WV_REAL32, 3*imax, (void*)plotdata, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* function values */
                MALLOC(pcolors, float, 9*imax);

                lims[0] = -1;  lims[1] = +1;
                for (i = 0; i < imax; i++) {
                    spec_col(plotdata[9*imax+3*i  ], &(pcolors[9*i  ]));
                    spec_col(plotdata[9*imax+3*i+1], &(pcolors[9*i+3]));
                    spec_col(plotdata[9*imax+3*i+2], &(pcolors[9*i+6]));
                }

                status = wv_setData(WV_REAL32, 3*imax, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT2(0, "ERROR:: wv_setData(%s) -> status=%d", temp, status);
                }
                nitems++;

                FREE(pcolors);
                FREE(plotdata);

                /* make graphic primitive */
                attrs = WV_ON | WV_SHADING;
                igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* filled Quads (imax=nquad, jmax=-4) */
            } else if (imax > 1 && jmax == -4) {
                SPRINT2(1, "   plotting %d filled quads (%s)", imax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotTris: %.11s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 24*imax);

                for (i = 0; i < imax; i++) {
                    fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                           &plotdata[18*i   ], &plotdata[18*i+ 1], &plotdata[18*i+ 2], &plotdata[18*imax+6*i  ],
                           &plotdata[18*i+ 3], &plotdata[18*i+ 4], &plotdata[18*i+ 5], &plotdata[18*imax+6*i+1],
                           &plotdata[18*i+ 6], &plotdata[18*i+ 7], &plotdata[18*i+ 8], &plotdata[18*imax+6*i+2],
                           &plotdata[18*i+ 9], &plotdata[18*i+10], &plotdata[18*i+11], &plotdata[18*imax+6*i+3]);

                    plotdata[18*i+12] = plotdata[18*i  ];
                    plotdata[18*i+13] = plotdata[18*i+1];
                    plotdata[18*i+14] = plotdata[18*i+2];
                    plotdata[18*i+15] = plotdata[18*i+6];
                    plotdata[18*i+16] = plotdata[18*i+7];
                    plotdata[18*i+17] = plotdata[18*i+8];

                    plotdata[18*imax+6*i+4] = plotdata[18*imax+6*i  ];
                    plotdata[18*imax+6*i+5] = plotdata[18*imax+6*i+2];
                }

                status = wv_setData(WV_REAL32, 6*imax, (void*)plotdata, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
                }

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* function values */
                MALLOC(pcolors, float, 18*imax);

                lims[0] = -1;  lims[1] = +1;
                for (i = 0; i < imax; i++) {
                    spec_col(plotdata[18*imax+6*i  ], &(pcolors[18*i   ]));
                    spec_col(plotdata[18*imax+6*i+1], &(pcolors[18*i+ 3]));
                    spec_col(plotdata[18*imax+6*i+2], &(pcolors[18*i+ 6]));
                    spec_col(plotdata[18*imax+6*i+3], &(pcolors[18*i+ 9]));
                    spec_col(plotdata[18*imax+6*i+4], &(pcolors[18*i+12]));
                    spec_col(plotdata[18*imax+6*i+5], &(pcolors[18*i+15]));
                }

                status = wv_setData(WV_REAL32, 6*imax, (void*)pcolors, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT2(0, "ERROR:: wv_setData(%s) -> status=%d", temp, status);
                }
                nitems++;

                FREE(pcolors);
                FREE(plotdata);

                /* make graphic primitive */
                attrs = WV_ON | WV_SHADING;
                igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* grid (imax, jmax) */
            } else if (imax > 1 && jmax > 1) {
                SPRINT3(1, "    plotting grid with %dx%d points (%s)", imax, jmax, temp);
                nitems = 0;

                /* name */
                snprintf(gpname, MAX_STRVAL_LEN-1, "PlotGrid: %.116s", temp);

                /* read the plotdata */
                MALLOC(plotdata, float, 3*imax*jmax);

                for (ij = 0; ij < imax*jmax; ij++) {
                    fscanf(fp, "%f %f %f", &plotdata[3*ij  ],
                                           &plotdata[3*ij+1],
                                           &plotdata[3*ij+2]);
                }

                /* create the segments */
                nseg = imax * (jmax-1) + (imax-1) * jmax;
                MALLOC(segments, float, 6*nseg);

                nseg = 0;

#ifndef __clang_analyzer__
                for (j = 0; j < jmax; j++) {
                    for (i = 0; i < imax-1; i++) {
                        ij1 = (i  ) + j * imax;
                        ij2 = (i+1) + j * imax;

                        segments[6*nseg  ] = plotdata[3*ij1  ];
                        segments[6*nseg+1] = plotdata[3*ij1+1];
                        segments[6*nseg+2] = plotdata[3*ij1+2];

                        segments[6*nseg+3] = plotdata[3*ij2  ];
                        segments[6*nseg+4] = plotdata[3*ij2+1];
                        segments[6*nseg+5] = plotdata[3*ij2+2];

                        nseg++;
                    }
                }

                for (i = 0; i < imax; i++) {
                    for (j = 0; j < jmax-1; j++) {
                        ij1 = i + (j  ) * imax;
                        ij2 = i + (j+1) * imax;

                        segments[6*nseg  ] = plotdata[3*ij1  ];
                        segments[6*nseg+1] = plotdata[3*ij1+1];
                        segments[6*nseg+2] = plotdata[3*ij1+2];

                        segments[6*nseg+3] = plotdata[3*ij2  ];
                        segments[6*nseg+4] = plotdata[3*ij2+1];
                        segments[6*nseg+5] = plotdata[3*ij2+2];

                        nseg++;
                    }
                }
#endif

                FREE(plotdata);

                status = wv_setData(WV_REAL32, 2*nseg, (void*)segments, WV_VERTICES, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(segments) -> status=%d", status);
                }

                FREE(segments);

                wv_adjustVerts(&(items[nitems]), sgFocus);
                nitems++;

                /* grid color */
                color[0] = 0;   color[1] = 0;   color[2] = 0;
                status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
                }
                nitems++;

                /* make graphic primitive and set line width */
                attrs = WV_ON;
                igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
                if (igprim < 0) {
                    SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
                } else {
                    SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                    cntxt->gPrims[igprim].lWidth = 1.0;
                }

                /* add plotdata to meta data (if there is room) */
                addToSgMetaData("\"%s\":[],", gpname);

            /* unknown type */
            } else {
                break;
            }

            /* break if we hit an end of file */
            if (feof(fp)) break;
        }

        fclose(fp);
    }

    /* BDF data to plot */
    if (STRLEN(BDFname) > 0) {
        fp = fopen(BDFname, "r");
        if (fp == NULL) {
            SPRINT1(0, "ERROR:: BDFname \"%s\" does not exist", BDFname);
            goto cleanup;
        } else {
            SPRINT1(1, "Opening \"%s\"", BDFname);
        }

        /* count the number of GRIDs in the file */
        rewind(fp);
        ngrid = 0;
        while (1) {
            if (fgets(text, 80, fp) == NULL) break;
            if (feof(fp)) break;

            if (strncmp(text, "GRID    ", 8) == 0) {
                ngrid++;
            }
        }
        SPRINT1(1, "   there are %d GRIDs", ngrid);

        if (ngrid <= 0) {
            fclose(fp);
            goto cleanup;
        }

        MALLOC(plotdata, float, 3*ngrid);

        nitems = 0;

        rewind(fp);
        i = 0;
        while (1) {
            if (fgets(text, 80, fp) == NULL) break;
            if (feof(fp)) break;

            if (strncmp(text, "GRID    ", 8) != 0) continue;

            sscanf(text, "%s %d %d %f %f %f", dum, &imax, &jmax, &plotdata[3*i  ],
                                                                 &plotdata[3*i+1],
                                                                 &plotdata[3*i+2]);
            i++;
        }

        /* name */
        snprintf(gpname, MAX_STRVAL_LEN-1, "PlotPoints: BDF_GRIDs");

        /* points in plotdata */
        status = wv_setData(WV_REAL32, ngrid, plotdata, WV_VERTICES, &(items[nitems]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(plotdata) -> status=%d", status);
        }

        wv_adjustVerts(&(items[nitems]), sgFocus);
        nitems++;

        /* point color */
        color[0] = 0;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }
        nitems++;

        /* make graphic primitive */
        attrs  = WV_ON;
        igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, nitems, items);
        if (igprim < 0) {
            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
        } else {
            SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

            cntxt->gPrims[igprim].pSize = 3.0;
        }

        /* add plotdata to meta data (if there is room) */
        addToSgMetaData("\"%s\":[],", gpname);

        /* count the number of CRODs in the file */
        rewind(fp);
        ncrod = 0;
        while (1) {
            if (fgets(text, 80, fp) == NULL) break;
            if (feof(fp)) break;

            if (strncmp(text, "CROD    ", 8) == 0) {
                ncrod++;
            }
        }
        SPRINT1(1, "   there are %d CRODs", ncrod);

        if (ncrod > 0) {
            MALLOC(segments, float, 6*ncrod);

            rewind(fp);
            i = 0;
            while (1) {
                if (fgets(text, 80, fp) == NULL) break;
                if (feof(fp)) break;

                if (strncmp(text, "CROD    ", 8) != 0) continue;

                sscanf(text, "%s %d %d %d %d", dum, &imax, &jmax, &ibeg, &iend);

                segments[6*i  ] = plotdata[3*ibeg-3];
                segments[6*i+1] = plotdata[3*ibeg-2];
                segments[6*i+2] = plotdata[3*ibeg-1];

                segments[6*i+3] = plotdata[3*iend-3];
                segments[6*i+4] = plotdata[3*iend-2];
                segments[6*i+5] = plotdata[3*iend-1];

                i++;
            }

            nitems = 0;

            /* name */
            snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: BDF_CRODs");

            /* segments */
            status = wv_setData(WV_REAL32, 2*i, (void*)segments, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(segments) -> status=%d", status);
            }

            FREE(segments);

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* line color */
            color[0] = 1.0;   color[1] = 0.5;   color[2] = 0.5;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive and set line width */
            attrs = WV_ON;
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* add plotdata to meta data (if there is room) */
            addToSgMetaData("\"%s\":[],", gpname);
        }

        /* count the number of CTRI3s in the file */
        rewind(fp);
        nctri3 = 0;
        while (1) {
            if (fgets(text, 80, fp) == NULL) break;
            if (feof(fp)) break;

            if (strncmp(text, "CTRI3   ", 8) == 0) {
                nctri3++;
            }
        }
        SPRINT1(1, "   there are %d CTRI3s", nctri3);

        if (nctri3 > 0) {
            MALLOC(segments, float, 18*nctri3);

            rewind(fp);
            i = 0;
            while (1) {
                if (fgets(text, 80, fp) == NULL) break;
                if (feof(fp)) break;

                if (strncmp(text, "CTRI3   ", 8) != 0) continue;

                sscanf(text, "%s %d %d %d %d %d", dum, &imax, &jmax, &isw, &ise, &ine);

                segments[18*i   ] = plotdata[3*isw-3];
                segments[18*i+ 1] = plotdata[3*isw-2];
                segments[18*i+ 2] = plotdata[3*isw-1];

                segments[18*i+ 3] = plotdata[3*ise-3];
                segments[18*i+ 4] = plotdata[3*ise-2];
                segments[18*i+ 5] = plotdata[3*ise-1];

                segments[18*i+ 6] = plotdata[3*ise-3];
                segments[18*i+ 7] = plotdata[3*ise-2];
                segments[18*i+ 8] = plotdata[3*ise-1];

                segments[18*i+ 9] = plotdata[3*ine-3];
                segments[18*i+10] = plotdata[3*ine-2];
                segments[18*i+11] = plotdata[3*ine-1];

                segments[18*i+12] = plotdata[3*ine-3];
                segments[18*i+13] = plotdata[3*ine-2];
                segments[18*i+14] = plotdata[3*ine-1];

                segments[18*i+15] = plotdata[3*isw-3];
                segments[18*i+16] = plotdata[3*isw-2];
                segments[18*i+17] = plotdata[3*isw-1];

                i++;
            }

            nitems = 0;

            /* name */
            snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: BDF_CTRI4s");

            /* segments */
            status = wv_setData(WV_REAL32, 6*i, (void*)segments, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(segments) -> status=%d", status);
            }

            FREE(segments);

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* line color */
            color[0] = 0.5;   color[1] = 1.0;   color[2] = 0.5;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive and set line width */
            attrs = WV_ON;
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* add plotdata to meta data (if there is room) */
            addToSgMetaData("\"%s\":[],", gpname);
        }

        /* count the number of CQUAD4s in the file */
        rewind(fp);
        ncquad4 = 0;
        while (1) {
            if (fgets(text, 80, fp) == NULL) break;
            if (feof(fp)) break;

            if (strncmp(text, "CQUAD4  ", 8) == 0) {
                ncquad4++;
            }
        }
        SPRINT1(1, "   there are %d CQUAD4s", ncquad4);

        if (ncquad4 > 0) {
            MALLOC(segments, float, 24*ncquad4);

            rewind(fp);
            i = 0;
            while (1) {
                if (fgets(text, 80, fp) == NULL) break;
                if (feof(fp)) break;

                if (strncmp(text, "CQUAD4  ", 8) != 0) continue;

                sscanf(text, "%s %d %d %d %d %d %d", dum, &imax, &jmax, &isw, &ise, &ine, &inw);

                segments[24*i   ] = plotdata[3*isw-3];
                segments[24*i+ 1] = plotdata[3*isw-2];
                segments[24*i+ 2] = plotdata[3*isw-1];

                segments[24*i+ 3] = plotdata[3*ise-3];
                segments[24*i+ 4] = plotdata[3*ise-2];
                segments[24*i+ 5] = plotdata[3*ise-1];

                segments[24*i+ 6] = plotdata[3*ise-3];
                segments[24*i+ 7] = plotdata[3*ise-2];
                segments[24*i+ 8] = plotdata[3*ise-1];

                segments[24*i+ 9] = plotdata[3*ine-3];
                segments[24*i+10] = plotdata[3*ine-2];
                segments[24*i+11] = plotdata[3*ine-1];

                segments[24*i+12] = plotdata[3*ine-3];
                segments[24*i+13] = plotdata[3*ine-2];
                segments[24*i+14] = plotdata[3*ine-1];

                segments[24*i+15] = plotdata[3*inw-3];
                segments[24*i+16] = plotdata[3*inw-2];
                segments[24*i+17] = plotdata[3*inw-1];

                segments[24*i+18] = plotdata[3*inw-3];
                segments[24*i+19] = plotdata[3*inw-2];
                segments[24*i+20] = plotdata[3*inw-1];

                segments[24*i+21] = plotdata[3*isw-3];
                segments[24*i+22] = plotdata[3*isw-2];
                segments[24*i+23] = plotdata[3*isw-1];

                i++;
            }

            nitems = 0;

            /* name */
            snprintf(gpname, MAX_STRVAL_LEN-1, "PlotLine: BDF_CQUAD4s");

            /* segments */
            status = wv_setData(WV_REAL32, 8*i, (void*)segments, WV_VERTICES, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(segments) -> status=%d", status);
            }

            FREE(segments);

            wv_adjustVerts(&(items[nitems]), sgFocus);
            nitems++;

            /* line color */
            color[0] = 0.5;   color[1] = 0.5;   color[2] = 1.0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[nitems]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }
            nitems++;

            /* make graphic primitive and set line width */
            attrs = WV_ON;
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, nitems, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            } else {
                SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

                cntxt->gPrims[igprim].lWidth = 1.0;
            }

            /* add plotdata to meta data (if there is room) */
            addToSgMetaData("\"%s\":[],", gpname);
        }

        FREE(plotdata);

        fclose(fp);
    }

    /* finish the scene graph meta data */
    sgMetaDataUsed--;
    addToSgMetaData("}");

cleanup:
    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    FREE(plotdata);
    FREE(pcolors );
    FREE(segments);
    FREE(segs    );
    FREE(ivrts   );
    FREE(tuft    );
    FREE(Tris    );

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   buildSceneGraphBody - build scene graph of one Body for wv        */
/*                                                                     */
/***********************************************************************/

static int
buildSceneGraphBody(int    ibody)       /* Body index (bias-1) */
{
    int       status = SUCCESS;         /* return status */

    int       iface, iedge, inode, atype, alen, attrs, head[3];
    int       npnt, ipnt, ntri, itri, igprim, nseg, k;
    int       *segs=NULL, *ivrts=NULL;
    CINT      *ptype, *pindx, *tris, *tric;
    float     color[18];
    double    xyz_dum[6], bbox[6], size, params[3];
    CDOUBLE   *xyz, *uv, *t;
    char      gpname[MAX_STRVAL_LEN];
    ego       etess;

    CINT      *tempIlist;
    CDOUBLE   *tempRlist;
    CCHAR     *tempClist;

    wvData    items[5];

    modl_T    *MODL = (modl_T*)modl;

    ROUTINE(buildSceneGraphBody);

    /* --------------------------------------------------------------- */

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    etess = MODL->body[ibody].etess;

    /* tessellate Body if not already tessellated */
    if (etess == NULL) {
        status = EG_attributeRet(MODL->body[ibody].ebody, "_tParams",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == EGADS_NOTFOUND) {
            status = EG_attributeRet(MODL->body[ibody].ebody, ".tParams",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == SUCCESS && alen == 3) {
                params[0] = tempRlist[0];
                params[1] = tempRlist[1];
                params[2] = tempRlist[2];

            } else {
                status = EG_getBoundingBox(MODL->body[ibody].ebody, bbox);
                CHECK_STATUS(EG_getBoundingBox);

                size = sqrt(SQR(bbox[3]-bbox[0]) + SQR(bbox[4]-bbox[1]) + SQR(bbox[5]-bbox[2]));

                params[0] = TESS_PARAM_0 * size;
                params[1] = TESS_PARAM_1 * size;
                params[2] = TESS_PARAM_2;
            }

            status = EG_attributeAdd(MODL->body[ibody].ebody, "_tParams",
                                     ATTRREAL, 3, NULL, params, NULL);
            CHECK_STATUS(EG_attributeAdd);
        }

        status = EG_attributeRet(MODL->body[ibody].ebody, "_tParams",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == SUCCESS && atype == ATTRREAL && alen == 3) {
            status = EG_makeTessBody(MODL->body[ibody].ebody, (double*)tempRlist,
                                     &(etess));
            if (etess != NULL) {
                MODL->body[ibody].etess = etess;
            }
        } else {
            SPRINT1(0, "ERROR:: cannot tessellate ibody %d", ibody);
        }
    }
    SPLINT_CHECK_FOR_NULL(etess);

    /* loop through the Faces within the current Body */
    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {

        /* name and attributes */
        snprintf(gpname, MAX_STRVAL_LEN-1, "Face %d", iface);
        attrs = WV_ON | WV_ORIENTATION;

        /* render the Triangles */
        status = EG_getTessFace(etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: EG_getTessFace -> status=%d", status);
        }

        /* skip if no Triangles either */
        if (ntri <= 0) continue;

        /* vertices */
        status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(xyz) -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), sgFocus);

        /* loop through the triangles and build up the segment table (bias-1) */
        nseg = 0;
        for (itri = 0; itri < ntri; itri++) {
            for (k = 0; k < 3; k++) {
                if (tric[3*itri+k] < itri+1) {
                    nseg++;
                }
            }
        }

        assert (nseg > 0);
        MALLOC(segs, int, 2*nseg);

        nseg = 0;
        for (itri = 0; itri < ntri; itri++) {
            for (k = 0; k < 3; k++) {
                if (tric[3*itri+k] < itri+1) {
                    segs[2*nseg  ] = tris[3*itri+(k+1)%3];
                    segs[2*nseg+1] = tris[3*itri+(k+2)%3];
                    nseg++;
                }
            }
        }

        /* triangles */
        status = wv_setData(WV_INT32, 3*ntri, (void*)tris, WV_INDICES, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(tris) -> status=%d", status);
        }

        color[0] = 1;   color[1] = 1;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[2]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }

        /* segment indices */
        status = wv_setData(WV_INT32, 2*nseg, (void*)segs, WV_LINDICES, &(items[3]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(segs) -> status=%d", status);
        }

        FREE(segs);

        /* segment colors */
        color[0] = 0;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[4]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }

        /* make graphic primitive */
        igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, 5, items);
        if (igprim < 0) {
            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
        } else {
            SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

            cntxt->gPrims[igprim].lWidth = 1.0;
        }
    }

    /* loop through the Edges within the current Body */
    for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
        status = EG_getTessEdge(etess, iedge, &npnt, &xyz, &t);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: EG_getTessEdge -> status=%d", status);
        }

        /* name and attributes */
        snprintf(gpname, MAX_STRVAL_LEN-1, "Edge %d", iedge);
        attrs = WV_ON;

        /* vertices */
        status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(xyz) -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), sgFocus);

        /* segments (bias-1) */
        MALLOC(ivrts, int, 2*(npnt-1));

        for (ipnt = 0; ipnt < npnt-1; ipnt++) {
            ivrts[2*ipnt  ] = ipnt + 1;
            ivrts[2*ipnt+1] = ipnt + 2;
        }

        status = wv_setData(WV_INT32, 2*(npnt-1), (void*)ivrts, WV_INDICES, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(ivrts) -> status=%d", status);
        }

        FREE(ivrts);

        /* line colors */
        color[0] = 0;   color[1] = 1;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[2]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }

        /* points */
        MALLOC(ivrts, int, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            ivrts[ipnt] = ipnt + 1;
        }

        status = wv_setData(WV_INT32, npnt, (void*)ivrts, WV_PINDICES, &(items[3]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(ivrts) -> status=%d", status);
        }

        FREE(ivrts);

        /* point colors */
        color[0] = 0;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[4]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }

        /* make graphic primitive */
        igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, 5, items);
        if (igprim < 0) {
            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
        } else {
            SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

            /* make line width 2 (does not work for ANGLE) */
            cntxt->gPrims[igprim].lWidth = 2.0;

            /* make point size 5 */
            cntxt->gPrims[igprim].pSize  = 5.0;

            /* add arrow heads (requires that WV_ORIENTATION be set above) */
            head[0] = npnt - 1;
            status = wv_addArrowHeads(cntxt, igprim, 0.10/sgFocus[3], 1, head);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_addArrowHeads -> status=%d", status);
            }
        }
    }

    /* if ibody is a NodeBody, draw the Node */
    if (MODL->body[ibody].botype == OCSM_NODE_BODY) {
        inode = 1;

        /* name and attributes */
        snprintf(gpname, MAX_STRVAL_LEN-1, "Node %d", inode);
        attrs = WV_ON;

        /* vertices */
        xyz_dum[0] = MODL->body[ibody].node[inode].x;
        xyz_dum[1] = MODL->body[ibody].node[inode].y;
        xyz_dum[2] = MODL->body[ibody].node[inode].z;
        xyz_dum[3] = MODL->body[ibody].node[inode].x;
        xyz_dum[4] = MODL->body[ibody].node[inode].y;
        xyz_dum[5] = MODL->body[ibody].node[inode].z;

        status = wv_setData(WV_REAL64, 2, (void*)xyz_dum, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(xyz) -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), sgFocus);

        /* point colors */
        color[0] = 0;   color[1] = 0;   color[2] = 0;
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
        }

        /* make graphic primitive */
        igprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, 2, items);
        if (igprim < 0) {
            SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
        } else {
            SPLINT_CHECK_FOR_NULL(cntxt->gPrims);

            /* make point size 5 */
            cntxt->gPrims[igprim].pSize  = 5.0;
        }
    }

    /* draw Edges for last SheetBody or SolidBody */
    ibody = MODL->nbody;
    while (ibody > 1) {
        if (MODL->body[ibody].botype == OCSM_SHEET_BODY ||
            MODL->body[ibody].botype == OCSM_SOLID_BODY   ) {
            break;
        } else {
            ibody--;
        }
    }

    if (ibody > 0) {
        etess = MODL->body[ibody].etess;

        /* loop through the Edges within the current Body */
        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(etess, iedge, &npnt, &xyz, &t);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: EG_getTessEdge -> status=%d", status);
            }

            /* name and attributes */
            snprintf(gpname, MAX_STRVAL_LEN-1, "Outline %d", iedge);
            attrs = WV_ON;

            /* vertices */
            status = wv_setData(WV_REAL64, npnt, (void*)xyz, WV_VERTICES, &(items[0]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(xyz) -> status=%d", status);
            }

            wv_adjustVerts(&(items[0]), sgFocus);

            /* segments (bias-1) */
            MALLOC(ivrts, int, 2*(npnt-1));

            for (ipnt = 0; ipnt < npnt-1; ipnt++) {
                ivrts[2*ipnt  ] = ipnt + 1;
                ivrts[2*ipnt+1] = ipnt + 2;
            }

            status = wv_setData(WV_INT32, 2*(npnt-1), (void*)ivrts, WV_INDICES, &(items[1]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(ivrts) -> status=%d", status);
            }

            FREE(ivrts);

            /* line colors */
            color[0] = 0.5;   color[1] = 0.5;   color[2] = 0.5;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[2]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }

            /* points */
            MALLOC(ivrts, int, npnt);

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                ivrts[ipnt] = ipnt + 1;
            }

            status = wv_setData(WV_INT32, npnt, (void*)ivrts, WV_PINDICES, &(items[3]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(ivrts) -> status=%d", status);
            }

            FREE(ivrts);

            /* point colors */
            color[0] = 0;   color[1] = 0;   color[2] = 0;
            status = wv_setData(WV_REAL32, 1, (void*)color, WV_PCOLOR, &(items[4]));
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: wv_setData(color) -> status=%d", status);
            }

            /* make graphic primitive */
            igprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, 5, items);
            if (igprim < 0) {
                SPRINT2(0, "ERROR:: wv_addGPrim(%s) -> igprim=%d", gpname, igprim);
            }
        }
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   cleanupMemory - clean up all memory in OpenCSM and EGADS          */
/*                                                                     */
/***********************************************************************/

static void
cleanupMemory(int    quiet)             /* (in)  =1 for no messages */
{
    int    status;
    ego    context;

    modl_T    *MODL = (modl_T*)modl;

    /* --------------------------------------------------------------- */

    if (MODL == NULL) return;

    /* remember the EGADS context (since it will be needed after
       MODL is removed) */
    context = MODL->context;

    if (quiet == 1) {
        outLevel = 0;
        (void) ocsmSetOutLevel(outLevel);
        (void) EG_setOutLevel(context, outLevel);
    }

    /* free up the modl */
    status = ocsmFree(MODL);
    SPRINT2(1, "--> ocsmFree() -> status=%d (%s)", status, ocsmGetText(status));

    /* clean up the udp storage */
    status = ocsmFree(NULL);
    SPRINT2(1, "--> ocsmFree(NULL) -> status=%d (%s)", status, ocsmGetText(status));

    /* remove the context */
    if (context != NULL) {
        status = EG_setOutLevel(context, 0);
        if (status < 0) {
            SPRINT1(0, "EG_setOutLevel -> status=%d", status);
        }

        status = EG_close(context);
        SPRINT1(1, "--> EG_close() -> status=%d", status);
    }
}


/***********************************************************************/
/*                                                                     */
/*   getToken - get a token from a string                              */
/*                                                                     */
/***********************************************************************/

static int
getToken(char   text[],                 /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   sep,                    /* (in)  separator character */
         char   token[])                /* (out) token */
{
    int    lentok, i, count, iskip;

    /* --------------------------------------------------------------- */

    token[0] = '\0';
    lentok   = 0;

    /* convert tabs to spaces */
    for (i = 0; i < STRLEN(text); i++) {
        if (text[i] == '\t') {
            text[i] = ' ';
        }
    }

    /* count the number of separators */
    count = 0;
    for (i = 0; i < STRLEN(text); i++) {
        if (text[i] == sep) {
            count++;
        }
    }

    if (count < nskip+1) return 0;

    /* skip over nskip tokens */
    i = 0;
    for (iskip = 0; iskip < nskip; iskip++) {
        while (text[i] != sep) {
            i++;
        }
        i++;
    }

    /* if token we are looking for is empty, return 0 */
    if (text[i] == sep) {
        token[0] = '0';
        token[1] = '\0';
    }

    /* extract the token we are looking for */
    while (text[i] != sep) {
        token[lentok++] = text[i++];
        token[lentok  ] = '\0';

        if (lentok >= MAX_EXPR_LEN-1) {
            SPRINT0(0, "ERROR:: token exceeds MAX_EXPR_LEN");
            break;
        }
    }

    return STRLEN(token);
}


/***********************************************************************/
/*                                                                     */
/*   maxDistance - compute maximum distance between two Bodys          */
/*                                                                     */
/***********************************************************************/

static int
maxDistance(modl_T  *MODL1,             /* (in)  first  MODL */
            modl_T  *MODL2,             /* (in)  second MODL */
            int     ibody,              /* (in)  Body index (bias-1) */
            double  *dist)              /* (out) maximum distance */
{
    int       status = SUCCESS;

    int       inode, jnode, iedge, jedge, iface, jface, ient;
    int       ipnt, jpnt, itype, jtype, oclass, mtype, nchild;
    int       npnt1, ntri1, npnt2, ntri2, *senses;
    int       *nMap=NULL, *eMap=NULL, *fMap=NULL;
    CINT      *ptype, *pindx, *tris, *tric, *tempIlist;
    double    dx, dy, dz, data1[18], data2[18];
    CDOUBLE   *xyz1, *t1, *uv1, *xyz2, *t2, *uv2, *tempRlist;
    CCHAR     *tempClist;
    ego       eref, *echilds;

    ROUTINE(maxDistance);

    /* --------------------------------------------------------------- */

    /* default return */
    *dist =  0;
    itype =  0;
    jtype = -1;
    ient  = -1;
    jpnt  = -1;

    /* check for valid MODLs and that Bodys match */
    if        (MODL1 == NULL || MODL2 == NULL) {
        status = OCSM_BODY_NOT_FOUND;
        goto cleanup;
    } else if (ibody < 1 || ibody > MODL1->nbody || ibody > MODL2->nbody) {
        status = OCSM_ILLEGAL_BODY_INDEX;
        goto cleanup;
    } else if (MODL1->body[ibody].nnode != MODL2->body[ibody].nnode) {
        SPRINT2(1, "MODL1->nnode=%d  MODL2->nnode=%d", MODL1->body[ibody].nnode, MODL2->body[ibody].nnode);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    } else if (MODL1->body[ibody].nedge != MODL2->body[ibody].nedge) {
        SPRINT2(1, "MODL1->nedge=%d  MODL2->nedge=%d", MODL1->body[ibody].nedge, MODL2->body[ibody].nedge);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    } else if (MODL1->body[ibody].nface != MODL2->body[ibody].nface) {
        SPRINT2(1, "MODL1->nface=%d  MODL2->nface=%d", MODL1->body[ibody].nface, MODL2->body[ibody].nface);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* get mapping info */
    status = EG_attributeRet(MODL1->body[ibody].ebody, ".nMap", &itype, &jtype,
                             &tempIlist, &tempRlist, &tempClist);
    if (status == SUCCESS) {
        nMap = (int*) tempIlist;
    }

    status = EG_attributeRet(MODL1->body[ibody].ebody, ".eMap", &itype, &jtype,
                             &tempIlist, &tempRlist, &tempClist);
    if (status == SUCCESS) {
        eMap = (int*) tempIlist;
    }

    status = EG_attributeRet(MODL1->body[ibody].ebody, ".fMap", &itype, &jtype,
                             &tempIlist, &tempRlist, &tempClist);
    if (status == SUCCESS) {
        fMap = (int*) tempIlist;
    }

    SPRINT3(1, "nMap=%llx  eMap=%llx  fMap=%llx",
            (long long)nMap, (long long)eMap, (long long)fMap);

    /* maximum distance between the Nodes */
    for (inode = 1; inode <= MODL1->body[ibody].nnode; inode++) {
        if (nMap == NULL) {
            jnode = inode;
        } else {
            jnode = nMap[inode];
        }

        status = EG_getTopology(MODL1->body[ibody].node[inode].enode,
                                &eref, &oclass, &mtype, data1, &nchild, &echilds, &senses);
        if (status < SUCCESS) goto cleanup;
        status = EG_getTopology(MODL2->body[ibody].node[jnode].enode,
                                &eref, &oclass, &mtype, data2, &nchild, &echilds, &senses);
        if (status < SUCCESS) goto cleanup;

        dx = data1[0] - data2[0];
        dy = data1[1] - data2[1];
        dz = data1[2] - data2[2];

        if (fabs(dx) > *dist) {
            jtype = 0;
            *dist = fabs(dx);
            itype = OCSM_NODE;
            ient  = inode;
        }
        if (fabs(dy) > *dist) {
            jtype = 1;
            *dist = fabs(dy);
            itype = OCSM_NODE;
            ient  = inode;
        }
        if (fabs(dz) > *dist) {
            jtype = 2;
            *dist = fabs(dz);
            itype = OCSM_NODE;
            ient  = inode;
        }
    }

    /* maximum distance between the Edges */
    for (iedge = 1; iedge <= MODL1->body[ibody].nedge; iedge++) {
        if (eMap == NULL) {
            jedge = iedge;
        } else {
            jedge = eMap[iedge];
        }

        status = EG_getTessEdge(MODL1->body[ibody].etess, iedge,
                                &npnt1, &xyz1, &t1);
        if (status < SUCCESS) goto cleanup;
        status = EG_getTessEdge(MODL2->body[ibody].etess, jedge,
                                &npnt2, &xyz2, &t2);
        if (status < SUCCESS) goto cleanup;

        if (npnt1 != npnt2) {
            SPRINT3(0, "ERROR:: iedge=%d: npnt1=%d, npnt2=%d", iedge, npnt1, npnt2);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }

        for (ipnt = 0; ipnt < npnt1; ipnt++) {
            dx = xyz1[3*ipnt  ] - xyz2[3*ipnt  ];
            dy = xyz1[3*ipnt+1] - xyz2[3*ipnt+1];
            dz = xyz1[3*ipnt+2] - xyz2[3*ipnt+2];

            if (fabs(dx) > *dist) {
                jtype = 0;
                *dist = fabs(dx);
                itype = OCSM_EDGE;
                ient  = iedge;
                jpnt  = ipnt;
            }
            if (fabs(dy) > *dist) {
                jtype = 1;
                *dist = fabs(dy);
                itype = OCSM_EDGE;
                ient  = iedge;
                jpnt  = ipnt;
            }
            if (fabs(dz) > *dist) {
                jtype = 2;
                *dist = fabs(dz);
                itype = OCSM_EDGE;
                ient  = iedge;
                jpnt  = ipnt;
            }
        }
    }

    /* maximum distance between the Faces */
    for (iface = 1; iface <= MODL1->body[ibody].nface; iface++) {
        if (fMap == NULL) {
            jface = iface;
        } else {
            jface = fMap[iface];
        }

        status = EG_getTessFace(MODL1->body[ibody].etess, iface,
                                &npnt1, &xyz1, &uv1, &ptype, &pindx,
                                &ntri1, &tris, &tric);
        if (status < SUCCESS) goto cleanup;
        status = EG_getTessFace(MODL2->body[ibody].etess, jface,
                                &npnt2, &xyz2, &uv2, &ptype, &pindx,
                                &ntri2, &tris, &tric);
        if (status < SUCCESS) goto cleanup;

        if (npnt1 != npnt2) {
            SPRINT3(0, "ERROR:: iface=%d: npnt1=%d, npnt2=%d", iface, npnt1, npnt2);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }

        for (ipnt = 0; ipnt < npnt1; ipnt++) {
            dx = xyz1[3*ipnt  ] - xyz2[3*ipnt  ];
            dy = xyz1[3*ipnt+1] - xyz2[3*ipnt+1];
            dz = xyz1[3*ipnt+2] - xyz2[3*ipnt+2];

            if (fabs(dx) > *dist) {
                jtype = 0;
                *dist = fabs(dx);
                itype = OCSM_FACE;
                ient  = iface;
                jpnt  = ipnt;
            }
            if (fabs(dy) > *dist) {
                jtype = 1;
                *dist = fabs(dy);
                itype = OCSM_FACE;
                ient  = iface;
                jpnt  = ipnt;
            }
            if (fabs(dz) > *dist) {
                jtype = 2;
                *dist = fabs(dz);
                itype = OCSM_FACE;
                ient  = iface;
                jpnt  = ipnt;
            }
        }
    }

    if        (jtype == 0) {
        SPRINT4(1, "maximum distance is dx=%e for %s %d (ipnt=%d)", *dist, ocsmGetText(itype), ient, jpnt);
    } else if (jtype == 1) {
        SPRINT4(1, "maximum distance is dy=%e for %s %d (ipnt=%d)", *dist, ocsmGetText(itype), ient, jpnt);
    } else if (jtype == 2) {
        SPRINT4(1, "maximum distance is dz=%e for %s %d (ipnt=%d)", *dist, ocsmGetText(itype), ient, jpnt);
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   mesgCallbackFromOpenCSM - post a message from OpenCSM             */
/*                                                                     */
/***********************************************************************/

void
mesgCallbackFromOpenCSM(char mesg[])    /* (in)  message */
{
    int  status=SUCCESS;

    ROUTINE(mesgCallbackFromOpenCSM);

    /* --------------------------------------------------------------- */

    /* make sure the messages buffer is big enough */
    if (messages_len+strlen(mesg) > max_mesg_len-2) {
        max_mesg_len += 4096;
        SPRINT1(2, "increasing max_mesg_len=%d", max_mesg_len);

        RALLOC(messages, char, max_mesg_len+1);
    }

    /* add the current message to the messages buffer */
    strcat(messages, mesg);
    strcat(messages, "\n");

    messages_len += strlen(mesg);

cleanup:
    if (status != SUCCESS) {
        SPRINT0(0, "ERROR:: max_mesg_len could not be increased");
    }

    return;
}


/***********************************************************************/
/*                                                                     */
/*   processBrowserToServer - process the message from client and create the response */
/*                                                                     */
/***********************************************************************/

static int
processBrowserToServer(char    text[])
{
    int       status = SUCCESS;

    int       i, ibrch, itype, nlist, builtTo, buildStatus, ichar, iundo;
    int       ipmtr, jpmtr, nrow, ncol, irow, icol, index, iattr, actv, itemp, linenum;
    int       itoken1, itoken2, itoken3, ibody, onstack, direction=1, nwarn;
    int       nclient;
    CINT      *tempIlist;
    double    scale, dihedral, value, dot;
    CDOUBLE   *tempRlist;
    char      *pEnd, bname[MAX_NAME_LEN+1], *bodyinfo=NULL;
    CCHAR     *tempClist;

    char      *name=NULL,  *type=NULL, *valu=NULL;
    char      *arg1=NULL,  *arg2=NULL, *arg3=NULL;
    char      *arg4=NULL,  *arg5=NULL, *arg6=NULL;
    char      *arg7=NULL,  *arg8=NULL, *arg9=NULL;
    char      *entry=NULL, *temp=NULL, *temp2=NULL;
    char      *matrix=NULL;

#define  MAX_TOKN_LEN  16384

    char      *begs=NULL, *vars=NULL, *cons=NULL, *segs=NULL, *vars_out=NULL;

    static FILE  *fp=NULL;

    modl_T    *MODL = (modl_T*)modl;
    modl_T    *saved_MODL = NULL;

    ROUTINE(processBrowserToServer);

    /* --------------------------------------------------------------- */

    MALLOC(name,     char, MAX_EXPR_LEN);
    MALLOC(type,     char, MAX_EXPR_LEN);
    MALLOC(valu,     char, MAX_EXPR_LEN);
    MALLOC(arg1,     char, MAX_EXPR_LEN);
    MALLOC(arg2,     char, MAX_EXPR_LEN);
    MALLOC(arg3,     char, MAX_EXPR_LEN);
    MALLOC(arg4,     char, MAX_EXPR_LEN);
    MALLOC(arg5,     char, MAX_EXPR_LEN);
    MALLOC(arg6,     char, MAX_EXPR_LEN);
    MALLOC(arg7,     char, MAX_EXPR_LEN);
    MALLOC(arg8,     char, MAX_EXPR_LEN);
    MALLOC(arg9,     char, MAX_EXPR_LEN);
    MALLOC(entry,    char, MAX_STR_LEN);
    MALLOC(temp,     char, MAX_EXPR_LEN);
    MALLOC(temp2,    char, MAX_EXPR_LEN);
    MALLOC(matrix,   char, MAX_EXPR_LEN);
    MALLOC(begs,     char, MAX_TOKN_LEN);
    MALLOC(vars,     char, MAX_TOKN_LEN);
    MALLOC(cons,     char, MAX_TOKN_LEN);
    MALLOC(segs,     char, MAX_TOKN_LEN);
    MALLOC(vars_out, char, MAX_TOKN_LEN);

    /* show message, except if one associated with syncing multiple clients */
    if (strncmp(text, "xform|",      6) != 0 &&
        strncmp(text, "lastPoint|", 10) != 0 &&
        strncmp(text, "toggle|",     7) != 0   ) {
        SPRINT1(1, "\n>>> browser2server(text=%s)", text);
    }

    /* initialize the response */
    response_len = 0;
    response[0] = '\0';

    /* NO-OP */
    if (STRLEN(text) == 0) {

    /* "identify|" */
    } else if (strncmp(text, "identify|", 9) == 0) {

        nclient = wv_nClientServer(serverNum);
        SPRINT0(1, "********************************************");
        SPRINT2(1, "server %d has %d clients", serverNum, nclient);
        SPRINT0(1, "********************************************");

        /* build the response */
        snprintf(response, max_resp_len, "identify|serveCSM|%d|", nclient);
        response_len = STRLEN(response);

    /* "userName|name|passTo|" */
    } else if (strncmp(text, "userName|", 9) == 0) {

        /* extract arguments */
        getToken(text, 1, '|', arg1);
        getToken(text, 2, '|', arg2);

        /* if second argument is *closed*, then remove user from usernames */
        if (strcmp(arg2, "*closed*") == 0) {
            strcpy(temp, usernames);
            strcpy(usernames, "|");
            for (i = 1; i < MAX_CLIENTS; i++) {
                getToken(temp, i, '|', arg3);

                /* end of users */
                if (STRLEN(arg3) <= 0) {
                    break;

                /* this user should be removed */
                } else if (strcmp(arg3, arg1) == 0) {

                    /* if this user has the ball, give it to first user in list */
                    if (i-1 == hasBall) {
                        hasBall = 0;

                    /* if the user with the ball is "above" the user to
                       be removed, decrement hasBall */
                    } else if (i <= hasBall) {
                        hasBall--;
                    }

                /* user still exists */
                } else {
                    strncat(usernames, arg3, 1023);
                    strncat(usernames, "|",  1023);
                }
            }

        /* otherwise add new user if not already in list */
        } else {
            snprintf(temp, MAX_EXPR_LEN,  "|%s|", arg1);
            if (strstr(usernames, temp) == 0) {
                strncat(usernames, arg1, 1023);
                strncat(usernames, "|",  1023);
            }
        }

        SPLINT_CHECK_FOR_NULL(usernames);

        /* update the browser with the ball if second argument is not blank */
        if (strlen(arg2) > 0) {
            for (i = 1; i < MAX_CLIENTS; i++) {
                getToken(usernames, i, '|', arg3);
                if (strcmp(arg2, arg3) == 0) {
                    hasBall = i-1;
                    break;
                }
            }
        }

        /* inform all browsers of changes */
        snprintf(response, max_resp_len, "userName|%d%s", hasBall, usernames);
        response_len = STRLEN(response);

    /* "xform|width|height|scale|matrix|" */
    } else if (strncmp(text, "xform|", 6) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "lastPoint|x|y|z|" or "lastPoint|off|" */
    } else if (strncmp(text, "lastPoint|", 10) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "toggle|inode|icol|state|" */
    } else if (strncmp(text, "toggle|", 7) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "resetMode|" */
    } else if (strncmp(text, "resetMode|", 10) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "sendState|" */
    } else if (strncmp(text, "sendState|", 10) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "message|text|" */
    } else if (strncmp(text, "message|", 8) == 0) {

        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    /* "nextStep|0|" */
    } else if (strncmp(text, "nextStep|0|", 11) == 0) {
        curStep = 0;
        buildSceneGraph();

        snprintf(response, max_resp_len, "nextStep|||");
        response_len = STRLEN(response);

    /* "nextStep|direction|" */
    } else if (strncmp(text, "nextStep|", 9) == 0) {

        if (getToken(text,  1, '|', arg1)) direction = strtol(arg1, &pEnd, 10);

        /* find next/previous SheetBody or SolidBody */
        if        (direction == +1 || direction == -1) {
            curStep  += direction;
        } else if (direction == +2) {
            curStep   = MODL->nbody;
            direction = -1;
        } else if (direction == -2) {
            curStep   = 1;
            direction = +1;
        } else {
            curStep   = 0;
        }

        while (curStep > 0 && curStep <= MODL->nbody) {
            if (MODL->body[curStep].botype == OCSM_NODE_BODY  ||
                MODL->body[curStep].botype == OCSM_WIRE_BODY  ||
                MODL->body[curStep].botype == OCSM_SHEET_BODY ||
                MODL->body[curStep].botype == OCSM_SOLID_BODY   ) {
                buildSceneGraphBody(curStep);

                ibrch = MODL->body[curStep].ibrch;

                status = EG_attributeRet(MODL->body[curStep].ebody, "_name", &itype, &nlist,
                                         &tempIlist, &tempRlist, &tempClist);
                if (status == SUCCESS && itype == ATTRSTRING) {
                    snprintf(bname, MAX_NAME_LEN-1, "%s", tempClist);
                } else {
                    snprintf(bname, MAX_NAME_LEN-1, "Body %d", curStep);
                }

                snprintf(response, max_resp_len, "nextStep|%d|%s|%s (%s)|",
                         ibrch, bname, MODL->brch[ibrch].name, ocsmGetText(MODL->brch[ibrch].type));
                response_len = STRLEN(response);
                break;
            } else {
                curStep += direction;
            }
        }

        /* if we did not find a WireBody, SheetBody, or SolidBody, we are
           done with StepThru mode */
        if (curStep < 1 || curStep > MODL->nbody) {
            curStep = 0;
            buildSceneGraph();

            snprintf(response, max_resp_len, "nextStep|||");
            response_len = STRLEN(response);
        }

    /* "getPmtrs|" */
    } else if (strncmp(text, "getPmtrs|", 9) == 0) {

        /* build the response in JSON format */
        snprintf(response, max_resp_len, "getPmtrs|[");
        response_len = STRLEN(response);

        /* constant Parameters first */
        if (MODL != NULL) {
            for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
                if (MODL->pmtr[ipmtr].type != OCSM_CONPMTR) continue;

                if (strlen(response) > 10) {
                    addToResponse(",");
                }

                snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                         MODL->pmtr[ipmtr].name,
                         MODL->pmtr[ipmtr].type,
                         MODL->pmtr[ipmtr].nrow,
                         MODL->pmtr[ipmtr].ncol);
                addToResponse(entry);

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].dot[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg]", MODL->pmtr[ipmtr].dot[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                addToResponse("}");
            }

            /* external and configuration Parameters second */
            for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
                if (MODL->pmtr[ipmtr].type != OCSM_DESPMTR &&
                    MODL->pmtr[ipmtr].type != OCSM_CFGPMTR   ) continue;

                if (strlen(response) > 10) {
                    addToResponse(",");
                }

                snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                         MODL->pmtr[ipmtr].name,
                         MODL->pmtr[ipmtr].type,
                         MODL->pmtr[ipmtr].nrow,
                         MODL->pmtr[ipmtr].ncol);
                addToResponse(entry);

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].dot[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg]", MODL->pmtr[ipmtr].dot[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                addToResponse("}");
            }

            /* internal Parameters last */
            for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
                if (MODL->pmtr[ipmtr].type != OCSM_LOCALVAR &&
                    MODL->pmtr[ipmtr].type != OCSM_OUTPMTR  ) continue;

                /* make sure the mass properties have been set up */
                if (MODL->pmtr[ipmtr].name[0] == '@') {
                    status = ocsmGetValu(MODL, ipmtr, 1, 1, &value, &dot);
                    if (status != SUCCESS) {
                        SPRINT2(0, "ERROR:: ocsmGetValu(%s) detected %s", MODL->pmtr[ipmtr].name, ocsmGetText(status));
                    }
                }

                /* skip if string-valued */
                if (MODL->pmtr[ipmtr].nrow == 0 || MODL->pmtr[ipmtr].ncol == 0) {
                    continue;
                }

                if (strlen(response) > 10) {
                    addToResponse(",");
                }

                snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                         MODL->pmtr[ipmtr].name,
                         MODL->pmtr[ipmtr].type,
                         MODL->pmtr[ipmtr].nrow,
                         MODL->pmtr[ipmtr].ncol);
                addToResponse(entry);

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                index = 0;
                for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                    for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                        if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                            snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].dot[index++]);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "%lg]", MODL->pmtr[ipmtr].dot[index++]);
                        }
                        addToResponse(entry);
                    }
                }

                addToResponse("}");
            }
            addToResponse("]");
        }

    /* "newPmtr|name|nrow|ncol|value1| ..." */
    } else if (strncmp(text, "newPmtr|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        nrow = 0;
        ncol = 0;

        if (getToken(text,  1, '|', name) == 0) name[0] = '\0';
        if (getToken(text,  2, '|', arg1)     ) nrow = strtol(arg1, &pEnd, 10);
        if (getToken(text,  3, '|', arg2)     ) ncol = strtol(arg2, &pEnd, 10);

        /* store an undo snapshot */
        status = storeUndo("newPmtr", name);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo(newPmtr) detected: %s", ocsmGetText(status));
        }

        /* build the response */
        status = ocsmNewPmtr(MODL, name, OCSM_DESPMTR, nrow, ncol);

        if (status == SUCCESS) {
            ipmtr = MODL->npmtr;

            i = 4;
            for (irow = 1; irow <= nrow; irow++) {
                for (icol = 1; icol <= ncol; icol++) {
                    if (getToken(text, i, '|', arg3)) {
                        (void) ocsmSetValu(MODL, ipmtr, irow, icol, arg3);
                    }

                    i++;
                }
            }

            snprintf(response, max_resp_len, "newPmtr|");
        } else {
            snprintf(response, max_resp_len, "ERROR:: newPmtr(%s,%s,%s) detected: %s",
                     name, arg1, arg2, ocsmGetText(status));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "setPmtr|pmtrname|irow|icol|value1| " */
    } else if (strncmp(text, "setPmtr|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        ipmtr = 0;
        irow  = 0;
        icol  = 0;

        getToken(text, 1, '|', arg1);

        for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
            if (strcmp(arg1, MODL->pmtr[jpmtr].name) == 0) {
                ipmtr = jpmtr;
                break;
            }
        }

        if (ipmtr > 0) {
            if (getToken(text, 2, '|', arg2)) irow  = strtol(arg2, &pEnd, 10);
            if (getToken(text, 3, '|', arg3)) icol  = strtol(arg3, &pEnd, 10);

            /* store an undo snapshot */
            status = storeUndo("setPmtr", MODL->pmtr[ipmtr].name);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: storeUndo(setPmtr) detected: %s", ocsmGetText(status));
            }

            if (getToken(text, 4, '|', arg4)) {
                status = ocsmSetValu(MODL, ipmtr, irow, icol, arg4);
                if (status != SUCCESS) {
                    SPRINT5(0, "ERROR:: ocsmSetValu(%d,%d,%d,%s) detected: %s",
                            ipmtr, irow, icol, arg4, ocsmGetText(status));
                }
            } else {
                status = -999;
            }

            /* build the response */
            if (status == SUCCESS) {
                snprintf(response, max_resp_len, "setPmtr|");
            } else {
                snprintf(response, max_resp_len, "setPmtr|ERROR:: %s",
                         ocsmGetText(status));
            }
        } else {
            snprintf(response, max_resp_len, "setPmtr|ERROR:: %s",
                     ocsmGetText(OCSM_NAME_NOT_FOUND));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "delPmtr|ipmtr|" */
    } else if (strncmp(text, "delPmtr|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        getToken(text, 1, '|', arg1);

        /* store an undo snapshot */
        status = storeUndo("delPmtr", arg1);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* delete the Parameter */
        ipmtr = 0;
        for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
            if (strcmp(MODL->pmtr[jpmtr].name, arg1) == 0) {
                ipmtr = jpmtr;
                break;
            }
        }

        if (ipmtr > 0) {
            status = ocsmDelPmtr(MODL, ipmtr);

            /* build the response */
            if (status == SUCCESS) {
                snprintf(response, max_resp_len, "delPmtr|");
            } else {
                snprintf(response, max_resp_len, "ERROR:: delPmtr(%s) detected: %s",
                         arg1, ocsmGetText(status));
            }
        } else {
            snprintf(response, max_resp_len, "ERROR:: delPmtr(%s detected: %s",
                     arg1, ocsmGetText(OCSM_NAME_NOT_FOUND));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "clrVels|mode|" */
    } else if (strncmp(text, "clrVels|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        getToken(text, 1, '|', arg1);

        if (strcmp(arg1, ".") == 0) {
            // do not change tessel flag
        } else if (strcmp(arg1, "tess") == 0) {
            tessel = 1;
        } else {                        // default is geometry sensitivities
            tessel = 0;
        }

        status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmSetVelD -> status=%d", status);
        }

        /* store an undo snapshot */
        status = storeUndo("clrVels", "");
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* build the response */
        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "clrVels|");
        } else {
            snprintf(response, max_resp_len, "ERROR:: clrVels() detected: %s",
                     ocsmGetText(status));
        }
        response_len = STRLEN(response);

    /* "setVel|pmtrname|irow|icol|vel|" */
    } else if (strncmp(text, "setVel|", 7) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        ipmtr = 0;
        irow  = 0;
        icol  = 0;

        getToken(text, 1, '|', arg1);

        for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
            if (strcmp(arg1, MODL->pmtr[jpmtr].name) == 0) {
                ipmtr = jpmtr;
                break;
            }
        }

        if (ipmtr > 0) {
            if (getToken(text, 2, '|', arg2)) irow  = strtol(arg2, &pEnd, 10);
            if (getToken(text, 3, '|', arg3)) icol  = strtol(arg3, &pEnd, 10);

            if (getToken(text, 4, '|', arg4)) {
                status = ocsmSetVel(MODL, ipmtr, irow, icol, arg4);
                if (status != SUCCESS) {
                    SPRINT1(0, "ERROR:: ocsmSetVel -> status=%d", status);
                }
            }

            /* store an undo snapshot */
            status = storeUndo("setVel", MODL->pmtr[ipmtr].name);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
            }

            /* build the response */
            if (status == SUCCESS) {
                snprintf(response, max_resp_len, "setVel|");
            } else {
                snprintf(response, max_resp_len, "ERROR:: setVel(%d,%d,%d) detected: %s",
                         ipmtr, irow, icol, ocsmGetText(status));
            }
        } else {
            snprintf(response, max_resp_len, "ERROR:: setVel(%s) detected: %s",
                     arg1, ocsmGetText(OCSM_NAME_NOT_FOUND));
        }
        response_len = STRLEN(response);

    /* "getBrchs|" */
    } else if (strncmp(text, "getBrchs|", 9) == 0) {

        /* build the response in JSON format */
        snprintf(response, max_resp_len, "getBrchs|[");
        response_len = STRLEN(response);

        if (MODL != NULL) {
            for (ibrch = 1; ibrch <= MODL->nbrch; ibrch++) {
                snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":\"%s\",\"actv\":%d,\"indent\":%d,\"level\":%d,\"attrs\":[",
                         MODL->brch[ibrch].name,
                         ocsmGetText(MODL->brch[ibrch].type),
                         MODL->brch[ibrch].actv,
                         MODL->brch[ibrch].indent,
                         MODL->brch[ibrch].level);
                addToResponse(entry);

                for (iattr = 0; iattr < MODL->brch[ibrch].nattr; iattr++) {
                    if (MODL->brch[ibrch].attr[iattr].type != ATTRCSYS) {
                        if (iattr < MODL->brch[ibrch].nattr-1) {
                            snprintf(entry, MAX_STR_LEN, "[\"%s\",\"(attr)\",\"%s\"],",
                                     MODL->brch[ibrch].attr[iattr].name,
                                     MODL->brch[ibrch].attr[iattr].defn);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "[\"%s\",\"(attr)\",\"%s\"]",
                                     MODL->brch[ibrch].attr[iattr].name,
                                     MODL->brch[ibrch].attr[iattr].defn);
                        }
                    } else {
                        if (iattr < MODL->brch[ibrch].nattr-1) {
                            snprintf(entry, MAX_STR_LEN, "[\"%s\",\"(csys)\",\"%s\"],",
                                     MODL->brch[ibrch].attr[iattr].name,
                                     MODL->brch[ibrch].attr[iattr].defn);
                        } else {
                            snprintf(entry, MAX_STR_LEN, "[\"%s\",\"(csys)\",\"%s\"]",
                                     MODL->brch[ibrch].attr[iattr].name,
                                     MODL->brch[ibrch].attr[iattr].defn);
                        }
                    }
                    addToResponse(entry);
                }

                snprintf(entry, MAX_STR_LEN, "],\"ileft\":%d,\"irite\":%d,\"ichld\":%d,\"args\":[",
                         MODL->brch[ibrch].ileft,
                         MODL->brch[ibrch].irite,
                         MODL->brch[ibrch].ichld);
                addToResponse(entry);

                if (MODL->brch[ibrch].narg >= 1) {
                    snprintf(entry, MAX_STR_LEN, "\"%s\"", MODL->brch[ibrch].arg1);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, "\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 2) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg2);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 3) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg3);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 4) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg4);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 5) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg5);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 6) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg6);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 7) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg7);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 8) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg8);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }
                if (MODL->brch[ibrch].narg >= 9) {
                    snprintf(entry, MAX_STR_LEN, ",\"%s\"", MODL->brch[ibrch].arg9);
                    addToResponse(entry);
                } else if (MODL->brch[ibrch].type == OCSM_UDPARG ||
                           MODL->brch[ibrch].type == OCSM_UDPRIM ||
                           MODL->brch[ibrch].type == OCSM_SELECT   ) {
                    snprintf(entry, MAX_STR_LEN, ",\"\"");
                    addToResponse(entry);
                }

                if (ibrch < MODL->nbrch) {
                    snprintf(entry, MAX_STR_LEN, "]},");
                } else {
                    snprintf(entry, MAX_STR_LEN, "]}]");
                }
                addToResponse(entry);
            }
        }

    /* "newBrch|ibrch|type|arg1|arg2|arg3|arg4|arg5|arg6|arg7|arg8|arg9|" */
    } else if (strncmp(text, "newBrch|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        ibrch = 0;
        itype = 0;

        if (getToken(text,  1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        if (getToken(text,  2, '|', type) >  0) itype = ocsmGetCode(type);
        if (getToken(text,  3, '|', arg1) == 0) arg1[0] = '\0';
        if (getToken(text,  4, '|', arg2) == 0) arg2[0] = '\0';
        if (getToken(text,  5, '|', arg3) == 0) arg3[0] = '\0';
        if (getToken(text,  6, '|', arg4) == 0) arg4[0] = '\0';
        if (getToken(text,  7, '|', arg5) == 0) arg5[0] = '\0';
        if (getToken(text,  8, '|', arg6) == 0) arg6[0] = '\0';
        if (getToken(text,  9, '|', arg7) == 0) arg7[0] = '\0';
        if (getToken(text, 10, '|', arg8) == 0) arg8[0] = '\0';
        if (getToken(text, 11, '|', arg9) == 0) arg9[0] = '\0';

        /* if this Branch is a udprim or udparg (which can have a variable
           number of arguments) change the zeros added by getToken to NULLs */
        if (itype == OCSM_UDPRIM || itype == OCSM_UDPARG) {
            if        (strcmp(arg2, "0") == 0) {
                arg2[0] = '\0';    arg3[0] = '\0';
                arg4[0] = '\0';    arg5[0] = '\0';
                arg6[0] = '\0';    arg7[0] = '\0';
                arg8[0] = '\0';    arg9[0] = '\0';
            } else if (strcmp(arg4, "0") == 0) {
                arg4[0] = '\0';    arg5[0] = '\0';
                arg6[0] = '\0';    arg7[0] = '\0';
                arg8[0] = '\0';    arg9[0] = '\0';
            } else if (strcmp(arg6, "0") == 0) {
                arg6[0] = '\0';    arg7[0] = '\0';
                arg8[0] = '\0';    arg9[0] = '\0';
            } else if (strcmp(arg8, "0") == 0) {
                arg8[0] = '\0';    arg9[0] = '\0';
            }
        }

        /* if this Branch is a select (which can have a variable number of
           arguments) change the trailing zeros added by getToken to NULLs */
        if (itype == OCSM_SELECT) {
            if (strcmp(arg9, "0") == 0) {
                arg9[0] = '\0';
                if (strcmp(arg8, "0") == 0) {
                    arg8[0] = '\0';
                    if (strcmp(arg7, "0") == 0) {
                        arg7[0] = '\0';
                        if (strcmp(arg6, "0") == 0) {
                            arg6[0] = '\0';
                            if (strcmp(arg5, "0") == 0) {
                                arg5[0] = '\0';
                                if (strcmp(arg4, "0") == 0) {
                                    arg4[0] = '\0';
                                    if (strcmp(arg3, "0") == 0) {
                                        arg3[0] = '\0';
                                        if (strcmp(arg2, "0") == 0) {
                                            arg2[0] = '\0';
                                            if (strcmp(arg1, "0") == 0) {
                                                arg1[0] = '\0';
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /* store an undo snapshot */
        status = storeUndo("newBrch", type);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* build the response */
        status = ocsmNewBrch(MODL, ibrch, itype, "", -1,
                             arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        if (status != SUCCESS) {
            snprintf(response, max_resp_len, "ERROR:: newBrch(%d,%d) detected: %s",
                     ibrch, itype, ocsmGetText(status));
            response_len = STRLEN(response);
            goto cleanup;
        }

        /* automatically add a SKEND if the Branch added was a SKBEG */
        if (itype == OCSM_SKBEG) {
            status = ocsmNewBrch(MODL, ibrch+1, OCSM_SKEND, "", -1,
                                 "0", "", "", "", "", "", "", "", "");
            if (status != SUCCESS) {
                snprintf(response, max_resp_len, "ERROR:: newBrch(%d,%d) detected: %s",
                         ibrch, OCSM_SKEND, ocsmGetText(status));
                response_len = STRLEN(response);
                goto cleanup;
            }
        }

        status = ocsmCheck(MODL);
        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "newBrch|");
        } else {
            snprintf(response, max_resp_len, "newBrch|WARNING:: %s",
                     ocsmGetText(status));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "setBrch|ibrch|name|actv|arg1|arg2|arg3|arg4|arg5|arg6|arg7|arg8|arg9| aname1|avalu1| ..." */
    } else if (strncmp(text, "setBrch|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        ibrch = 0;
        actv  = 0;

        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        /* store an undo snapshot */
        status = storeUndo("setBrch", MODL->brch[ibrch].name);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* build the response */
        if (ibrch >= 1 && ibrch <= MODL->nbrch) {
            if (getToken(text,  2, '|', name)) (void) ocsmSetName(MODL, ibrch, name);

            if (getToken(text,  3, '|', arg1)) {
                if (strcmp(arg1, "suppressed") == 0) {
                    (void) ocsmSetBrch(MODL, ibrch, OCSM_SUPPRESSED);
                    actv = 1;
                } else {
                    (void) ocsmSetBrch(MODL, ibrch, OCSM_ACTIVE);
                    actv = 1;
                }
            }

            if (getToken(text,  4, '|', arg1)) (void) ocsmSetArg(MODL, ibrch, 1, arg1);
            if (getToken(text,  5, '|', arg2)) (void) ocsmSetArg(MODL, ibrch, 2, arg2);
            if (getToken(text,  6, '|', arg3)) (void) ocsmSetArg(MODL, ibrch, 3, arg3);
            if (getToken(text,  7, '|', arg4)) (void) ocsmSetArg(MODL, ibrch, 4, arg4);
            if (getToken(text,  8, '|', arg5)) (void) ocsmSetArg(MODL, ibrch, 5, arg5);
            if (getToken(text,  9, '|', arg6)) (void) ocsmSetArg(MODL, ibrch, 6, arg6);
            if (getToken(text, 10, '|', arg7)) (void) ocsmSetArg(MODL, ibrch, 7, arg7);
            if (getToken(text, 11, '|', arg8)) (void) ocsmSetArg(MODL, ibrch, 8, arg8);
            if (getToken(text, 12, '|', arg9)) (void) ocsmSetArg(MODL, ibrch, 9, arg9);

            i = 13;
            while (1) {
                if (getToken(text, i++, '|', name) == 0) break;
                if (getToken(text, i++, '|', valu) == 0) break;

                if (strcmp(name, "0") == 0) break;

                (void) ocsmSetAttr(MODL, ibrch, name, valu);
            }

            if (actv > 0) {
                status = ocsmCheck(MODL);

                if (status >= SUCCESS) {
                    snprintf(response, max_resp_len, "setBrch|");
                } else {
                    snprintf(response, max_resp_len, "setBrch|WARNING:: %s",
                             ocsmGetText(status));
                }
            } else {
                snprintf(response, max_resp_len, "setBrch|");
            }
        } else {
            status = OCSM_ILLEGAL_BRCH_INDEX;
            snprintf(response, max_resp_len, "ERROR: setBrch(%d) detected: %s",
                     ibrch, ocsmGetText(status));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "delBrch|ibrch|" */
    } else if (strncmp(text, "delBrch|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        ibrch = 0;

        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        /* store an undo snapshot */
        status = storeUndo("delBrch", MODL->brch[ibrch].name);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* delete the Branch */
        status = ocsmDelBrch(MODL, ibrch);

        /* check that the Branches are properly ordered */
        if (status == SUCCESS) {
            status = ocsmCheck(MODL);

            /* build the response */
            if (status == SUCCESS) {
                snprintf(response, max_resp_len, "delBrch|");
            } else {
                snprintf(response, max_resp_len, "delBrch|WARNING:: %s",
                         ocsmGetText(status));
            }
        } else {
            snprintf(response, max_resp_len, "ERROR: delBrch(%d) detected: %s",
                     ibrch, ocsmGetText(status));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "setAttr|ibrch|aname|atype|avalue|" */
    } else if (strncmp(text, "setAttr|", 8) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        ibrch = 0;

        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);
        getToken(text, 2, '|', arg2);
        getToken(text, 3, '|', arg3);
        getToken(text, 4, '|', arg4);

        /* store an undo snapshot */
        status = storeUndo("setAttr", MODL->brch[ibrch].name);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: storeUndo -> status=%d", status);
        }

        /* special code to delete an Attribute */
        if (strcmp(arg4, "<DeLeTe>") == 0) {
            arg4[0] = '\0';
        }

        /* set the Attribute */
        if (strcmp(arg3, "2") == 0) {
            status = ocsmSetCsys(MODL, ibrch, arg2, arg4);
        } else {
            status = ocsmSetAttr(MODL, ibrch, arg2, arg4);
        }

        /* build the response */
        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "setAttr|");
        } else {
            snprintf(response, max_resp_len, "ERROR: setAttr(%d,%s,%s,%s) detected: %s",
                     ibrch, arg2, arg3, arg4, ocsmGetText(status));
        }
        response_len = STRLEN(response);

        /* write autosave file */
        status = ocsmSave(MODL, "autosave.csm");
        SPRINT1(2, "ocsmSave(autosave.csm) -> status=%d", status);

    /* "undo|" */
    } else if (strncmp(text, "undo|", 5) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* build the response */
        if (nundo <= 0) {
            snprintf(response, max_resp_len, "ERROR:: there is nothing to undo");

        } else {
            /* remove the current MODL */
            status = ocsmFree(modl);

            if (status < SUCCESS) {
                snprintf(response, max_resp_len, "ERROR:: undo() detected: %s",
                         ocsmGetText(status));
            } else {

                /* repoint MODL to the saved modl */
                modl = undo_modl[--nundo];
                snprintf(response, max_resp_len, "undo|%s|",
                         undo_text[nundo]);
            }
        }
        response_len = STRLEN(response);

    /* "new|" */
    } else if (strncmp(text, "new|", 4) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            rewind(jrnl_out);
        }

        /* remove previous undos (if any) */
        for (iundo = nundo-1; iundo >= 0; iundo--) {
            (void) ocsmFree(undo_modl[iundo]);
        }

        /* remove undo information */
        nundo = 0;

        /* free up the current MODL */
        status = ocsmFree(MODL);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmFree -> status=%d", status);
        }

        /* load an empty MODL */
        filename[0] = '\0';
        status = ocsmLoad(filename, &modl);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: osmLoad(NULL) -> status=%d", status);
        }

        MODL = (modl_T *)modl;
        if(filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        updatedFilelist = 1;

        status = ocsmLoadDict(modl, dictname);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmLoadDict -> status=%d", status);
        }

        status = ocsmRegMesgCB(modl, mesgCallbackFromOpenCSM);
        if (status < EGADS_SUCCESS) goto cleanup;

        status = ocsmRegSizeCB(modl, sizeCallbackFromOpenCSM);
        if (status < EGADS_SUCCESS) goto cleanup;

        if (strlen(despname) > 0) {
            status = ocsmUpdateDespmtrs(MODL, despname);
            if (status < EGADS_SUCCESS) goto cleanup;
        }

        status = buildBodys(0, &builtTo, &buildStatus, &nwarn);

        /* build the response */
        if (status == SUCCESS && buildStatus == SUCCESS) {
            snprintf(response, max_resp_len, "new|");
        } else {
            snprintf(response, max_resp_len, "ERROR:: new detected: %s",
                     ocsmGetText(status));
        }
        response_len = STRLEN(response);

    /* "open|filename|" */
    } else if (strncmp(text, "open|", 5) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            rewind(jrnl_out);
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* remove previous undos (if any) */
        for (iundo = nundo-1; iundo >= 0; iundo--) {
            (void) ocsmFree(undo_modl[iundo]);
        }

        /* remove undo information */
        nundo = 0;

        /* extract argument */
        getToken(text, 1, '|', filename);

        /* save the current MODL (to be deleted below) */
        saved_MODL = (modl_T *)modl;

        /* load the new MODL */
        status = ocsmLoad(filename, &modl);
        if (status != SUCCESS) {
            MODL = (modl_T *) modl;

            snprintf(response, max_resp_len, "%s||",
                     MODL->sigMesg);

            buildSceneGraph();
        } else {
            MODL = (modl_T *) modl;

            status = ocsmLoadDict(MODL, dictname);
            if (status != SUCCESS) {
                SPRINT2(0, "ERROR:: ocsmLoadDict(%s) detected %s",
                        dictname, ocsmGetText(status));
            }

            status = ocsmRegMesgCB(MODL, mesgCallbackFromOpenCSM);
            if (status < EGADS_SUCCESS) goto cleanup;

            status = ocsmRegSizeCB(MODL, sizeCallbackFromOpenCSM);
            if (status < EGADS_SUCCESS) goto cleanup;

            if (strlen(despname) > 0) {
                status = ocsmUpdateDespmtrs(MODL, despname);
                if (status < EGADS_SUCCESS) goto cleanup;
            }

            /* apply and free up a saved MODL (if it exists) */
            if (saved_MODL != NULL) {
                status = updateModl(saved_MODL, MODL);
                if (status < EGADS_SUCCESS) goto cleanup;

                status = ocsmFree(saved_MODL);
                if (status < EGADS_SUCCESS) goto cleanup;
            }

            status = buildBodys(0, &builtTo, &buildStatus, &nwarn);

            if (status != SUCCESS || buildStatus != SUCCESS) {
                snprintf(response, max_resp_len, "%s|%s|",
                         MODL->sigMesg, messages);
            } else {
                onstack = 0;
                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    onstack += MODL->body[ibody].onstack;
                }

                snprintf(response, max_resp_len, "build|%d|%d|%s|",
                         abs(builtTo), onstack, messages);
            }

            messages[0] = '\0';
            messages_len = 0;
        }

        if (filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        updatedFilelist = 1;

        response_len = STRLEN(response);

    /* "save|filename|" */
    } else if (strncmp(text, "save|", 5) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        getToken(text, 1, '|', filename);

        /* save the file */
        status = ocsmSave(MODL, filename);

        /* build the response */
        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "save|");
        } else {
            snprintf(response, max_resp_len, "ERROR:: save(%s) detected: %s",
                     filename, ocsmGetText(status));
        }
        response_len = STRLEN(response);

    /* "getFilenames|" */
    } else if (strncmp(text, "getFilenames|", 13) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        MODL = (modl_T *)modl;
        if (filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        SPLINT_CHECK_FOR_NULL(filelist);
        updatedFilelist = 0;

        /* build the response */
        snprintf(response, max_resp_len, "getFilenames|%s",
                 filelist);
        response_len = STRLEN(response);

    /* "getCsmFile|" */
    } else if (strncmp(text, "getCsmFile|", 11) == 0) {
        char subfilename[MAX_FILENAME_LEN];

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        getToken(text, 1, '|', subfilename);

        /* build the response */
        snprintf(response, max_resp_len, "getCsmFile|");
        response_len = STRLEN(response);

        /* add the input file, line by line */
        if (STRLEN(subfilename) > 0) {
            fp = fopen(subfilename, "r");
            if (fp != NULL) {
                while (1) {
                    if (fgets(entry, MAX_STR_LEN-1, fp) == NULL) break;
                    addToResponse(entry);
                    if (feof(fp) != 0) break;
                }

                fclose(fp);
                fp = NULL;
            }
        }

    /* "setCsmFileBeg|" */
    } else if (strncmp(text, "setCsmFileBeg|", 14) == 0) {
        char subfilename[MAX_FILENAME_LEN];

        /* extract argument */
        getToken(text, 1, '|', subfilename);

        /* over-writing the .csm file means that everything that was
           done previous cannot be re-done (since the original input
           file no longer exists).  therefore, start a new journal file */
        if (jrnl_out != NULL) {
            rewind(jrnl_out);
            response_len = 0;
            response[0]  = '\0';

            fprintf(jrnl_out, "open|%s|\n", subfilename);
        }

        /* open the casefile and overwrite */
        fp = fopen(subfilename, "w");
        SPLINT_CHECK_FOR_NULL(fp);

        ichar = 14;
        while (text[ichar++] != '|') {
        }
        while (text[ichar] != '\0') {
            fprintf(fp, "%c", text[ichar++]);
        }

        /* update filename if file was the .csm file */
        if (strstr(subfilename, ".csm") != NULL) {
            strcpy(filename, subfilename);
        }

    /* "setCsmFileMid|" */
    } else if (strncmp(text, "setCsmFileMid|", 14) == 0) {
        SPLINT_CHECK_FOR_NULL(fp);

        ichar = 14;
        while (text[ichar] != '\0') {
            fprintf(fp, "%c", text[ichar++]);
        }

    /* "setCsmFileEnd|" */
    } else if (strncmp(text, "setCsmFileEnd|", 14) == 0) {
        SPLINT_CHECK_FOR_NULL(fp);

        fclose(fp);
        fp = NULL;

        /* save the current MODL (to be deleted below) */
        saved_MODL = (modl_T *) modl;

        /* load the new MODL */
        status = ocsmLoad(filename, &modl);

        if (status != SUCCESS) {
            MODL = (modl_T *) modl;

            snprintf(response, max_resp_len, "%s||",
                     MODL->sigMesg);

            /* clear any previous builds from the scene graph */
            buildSceneGraph();
        } else {
            MODL = (modl_T *) modl;

            status = ocsmLoadDict(modl, dictname);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmLoadDict -> status=%d", status);
            }

            status = ocsmRegMesgCB(modl, mesgCallbackFromOpenCSM);
            if (status < EGADS_SUCCESS) goto cleanup;

            status = ocsmRegSizeCB(modl, sizeCallbackFromOpenCSM);
            if (status < EGADS_SUCCESS) goto cleanup;

            if (strlen(despname) > 0) {
                status = ocsmUpdateDespmtrs(MODL, despname);
                if (status < EGADS_SUCCESS) goto cleanup;
            }

            status = updateModl(saved_MODL, MODL);
            if (status < EGADS_SUCCESS) goto cleanup;

            snprintf(response, max_resp_len, "load|");
        }

        MODL = (modl_T *)modl;
        if (filelist != NULL) EG_free(filelist);
        status = ocsmGetFilelist(MODL, &filelist);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmGetFilelist -> status=%d", status);
        }
        updatedFilelist = 1;

        /* free up the saved MODL */
        status = ocsmFree(saved_MODL);
        if (status != SUCCESS) {
            SPRINT1(0, "ERROR:: ocsmFree -> status=%d", status);
        }

        /* disable -loadEgads */
        loadEgads = 0;

        response_len = STRLEN(response);

    /* "build|"  */
    } else if (strncmp(text, "build|", 6) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        ibrch = 0;
        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        /* if ibrch is negative, clear the velocities */
        if (ibrch < 0) {
            status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmSetVelD -> status=%d", status);
            }
        }

        /* build the response */
        status = buildBodys(ibrch, &builtTo, &buildStatus, &nwarn);

        if (status != SUCCESS || buildStatus != SUCCESS) {
            snprintf(response, max_resp_len, "%s|%s|",
                     MODL->sigMesg, messages);
        } else {
            onstack = 0;
            for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                onstack += MODL->body[ibody].onstack;
            }

            snprintf(response, max_resp_len, "build|%d|%d|%s|",
                     abs(builtTo), onstack, messages);
        }

        messages[0] = '\0';
        messages_len = 0;

        /* disable -loadEgads */
        loadEgads = 0;

        response_len = STRLEN(response);

    /* "recycle|"  */
    } else if (strncmp(text, "recycle|", 6) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        ibrch = 0;
        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        /* if ibrch is negative, clear the velocities */
        if (ibrch < 0) {
            status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
            if (status != SUCCESS) {
                SPRINT1(0, "ERROR:: ocsmSetVelD -> status=%d", status);
            }
        }

        /* build the response */
        status = buildBodys(ibrch, &builtTo, &buildStatus, &nwarn);

        goto cleanup;

    /* "getBodyDetails|filename|linenum||" */
    } else if (strncmp(text, "getBodyDetails|", 15) == 0) {

        /* extract arguments */
        linenum = 0;
        getToken(text, 1, '|', arg1);
        if (getToken(text, 2, '|', arg2)) linenum = strtol(arg2, &pEnd, 10);

        status = ocsmBodyDetails(MODL, arg1, linenum, &bodyinfo);

        SPLINT_CHECK_FOR_NULL(bodyinfo);

        if (status == SUCCESS) {
            itemp = 25 + STRLEN(bodyinfo);

            if (itemp > max_resp_len) {
                max_resp_len = itemp + 1;

                RALLOC(response, char, max_resp_len+1);
            }

            /* build the response */
            snprintf(response, max_resp_len, "getBodyDetails|%s|%d|%s|", arg1, linenum, bodyinfo);
            response_len = STRLEN(response);
        }

        FREE(bodyinfo);

    /* "loadSketch|" */
    } else if (strncmp(text, "loadSketch|", 11) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        ibrch = 0;
        if (getToken(text, 1, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        status = ocsmGetSketch(MODL, ibrch, MAX_TOKN_LEN, begs, vars, cons, segs);

        itemp = 20 + STRLEN(begs) + STRLEN(vars) + STRLEN(cons) + STRLEN(segs);

        if (itemp > max_resp_len) {
            max_resp_len = itemp + 1;

            RALLOC(response, char, max_resp_len+1);
        }

        if (status != SUCCESS) {
            snprintf(response, max_resp_len, "loadSketch|%s|",
                     MODL->sigMesg);
            response_len = STRLEN(response);
        } else {
            snprintf(response, max_resp_len, "loadSketch|%s|%s|%s|%s|",
                     begs, vars, cons, segs);
            response_len = STRLEN(response);
        }

    /* "solveSketch|" */
    } else if (strncmp(text, "solveSketch|", 12) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        itoken1 = -1;
        itoken2 = -1;
        itemp   = 12;

        while (itemp < STRLEN(text)) {
            if (text[itemp++] == '|') {
                itoken1 = itemp;
                break;
            }
        }
        if (itoken1 < 0) {
            snprintf(response, max_resp_len, "solveSketch|error extracting token1");
            response_len = STRLEN(response);
            goto cleanup;
        }

        while (itemp < STRLEN(text)) {
            if (text[itemp++] == '|') {
                itoken2 = itemp;
                break;
            }
        }
        if (itoken2 < 0) {
            snprintf(response, max_resp_len, "solveSketch|error extracting token2");
            response_len = STRLEN(response);
            goto cleanup;
        }

        text[itoken1-1] = '\0';
        text[itoken2-1] = '\0';

        /* solve the sketch (which always returns a response) */
        status = ocsmSolveSketch(MODL, &(text[12]), &(text[itoken1]), vars_out);

        if (status < SUCCESS) {
            snprintf(response, max_resp_len, "solveSketch|ERROR:: %s|", MODL->sigMesg);
        } else if (STRLEN(vars_out) == 0) {
            snprintf(response, max_resp_len, "solveSketch|%s|", MODL->sigMesg);
        } else {
            snprintf(response, max_resp_len, "solveSketch|%s|", vars_out);
        }
        response_len = STRLEN(response);

    /* "saveSketchBeg|" */
    } else if (strncmp(text, "saveSketchBeg|", 14) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* start saving into a buffer */
        strncpy(skbuff, &(text[14]), MAX_STR_LEN);

    /* "saveSketchMid|" */
    } else if (strncmp(text, "saveSketchMid|", 14) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        strncat(skbuff, &(text[14]), MAX_STR_LEN);

    /* "saveSketchEnd|" */
    } else if (strncmp(text, "saveSketchEnd|", 11) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        strncat(skbuff, &(text[14]), MAX_STR_LEN);

        /* extract arguments */
        ibrch = 0;
        if (getToken(skbuff, 0, '|', arg1)) ibrch = strtol(arg1, &pEnd, 10);

        itoken1 = -1;
        itoken2 = -1;
        itoken3 = -1;
        itemp   =  0;

        while (itemp < STRLEN(skbuff)) {
            if (skbuff[itemp++] == '|') {
                itoken1 = itemp;
                break;
            }
        }
        if (itoken1 < 0) {
            snprintf(response, max_resp_len, "saveSketch|error extracting token1");
            response_len = STRLEN(response);
            goto cleanup;
        }

        while (itemp < STRLEN(skbuff)) {
            if (skbuff[itemp++] == '|') {
                itoken2 = itemp;
                break;
            }
        }
        if (itoken2 < 0) {
            snprintf(response, max_resp_len, "saveSketch|error extracting token2");
            response_len = STRLEN(response);
            goto cleanup;
        }

        while (itemp < STRLEN(skbuff)) {
            if (skbuff[itemp++] == '|') {
                itoken3 = itemp;
                break;
            }
        }
        if (itoken3 < 0) {
            snprintf(response, max_resp_len, "saveSketch|error extracting token3");
            response_len = STRLEN(response);
            goto cleanup;
        }

        skbuff[itoken1-1] = '\0';
        skbuff[itoken2-1] = '\0';
        skbuff[itoken3-1] = '\0';

        /* save the sketch */
        status = ocsmSaveSketch(MODL, ibrch, &(skbuff[itoken1]), &(skbuff[itoken2]), &(skbuff[itoken3]));

        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "saveSketch|ok|");
            response_len = STRLEN(response);
        } else {
            snprintf(response, max_resp_len, "saveSketch|error|");
            response_len = STRLEN(response);
        }

    /* "makeEBody|ibody|dihedral|ents|" */
    } else if (strncmp(text, "makeEBody|", 10) == 0) {

        /* extract arguments */
        ibody    = 0;
        dihedral = 0;
        if (getToken(text, 1, '|', arg1)) ibody    = strtol(arg1, &pEnd, 10);
        if (getToken(text, 2, '|', arg2)) dihedral = strtod(arg2, &pEnd);
        getToken(text, 3, '|', arg3);

        /* make or delete the EBody */
        if (strcmp(arg3, ".") == 0) {
            status = ocsmMakeEBody(MODL, ibody, dihedral, NULL);
        } else {
            status = ocsmMakeEBody(MODL, ibody, dihedral, arg3);
        }

        /* update the scene graph before the response */
        if (status == SUCCESS && batch == 0) {
            buildSceneGraph();
        }

        /* build the response */
        if (status == SUCCESS) {
            snprintf(response, max_resp_len, "makeEBody|");
        } else {
            snprintf(response, max_resp_len, "makeEBody|ERROR:: unable to make EBody(s) status=%d", status);
        }
        response_len = STRLEN(response);

    /* "setLims|type|lo|hi|" */
    } else if (strncmp(text, "setLims|", 8) == 0) {

        /* extract arguments */
        if (getToken(text, 1, '|', arg1)) plotType = strtod(arg1, &pEnd);
        if (getToken(text, 2, '|', arg2)) lims[0]  = strtod(arg2, &pEnd);
        if (getToken(text, 3, '|', arg3)) lims[1]  = strtod(arg3, &pEnd);

        /* handle special case of Ereps */
        if (plotType < 7) {
            if (MODL->erepAtEnd == 1) {
                status = buildBodys(0, &builtTo, &buildStatus, &nwarn);
                if (status != SUCCESS) goto cleanup;
            }
        } else {
            if (MODL->erepAtEnd == 0) {
                status = buildBodys(0, &builtTo, &buildStatus, &nwarn);
                if (status != SUCCESS) goto cleanup;
            }
        }

        /* build the response */
        snprintf(response, max_resp_len, "setLims|");
        response_len = STRLEN(response);

        /* update the scene graph (after clearing the meta data) */
        if (batch == 0) {
            buildSceneGraph();
        }

    /* "saveView|viewfile|scale|array|" */
    } else if (strncmp(text, "saveView|", 9) == 0) {
        char viewfile[MAX_FILENAME_LEN];

        /* extract arguments */
        scale = 1;

        getToken(text, 1, '|', viewfile);
        if (getToken(text, 2, '|', arg2)) scale = strtod(arg2, &pEnd);
        getToken(text, 3, '|', matrix  );

        /* save the view in viewfile */
        fp = fopen(viewfile, "w");
        if (fp != NULL) {
            fprintf(fp, "%f\n", scale);
            fprintf(fp, "%s\n", matrix);
            fclose(fp);
        }

        /* build the response */
        snprintf(response, max_resp_len, "saveView|");
        response_len = STRLEN(response);

    /* "readView|viewfile|" */
    } else if (strncmp(text, "readView|", 9) == 0) {
        char viewfile[MAX_FILENAME_LEN];

        /* estract arguments */
        getToken(text, 1, '|', viewfile);

        /* read the view from viewfile */
        fp = fopen(viewfile, "r");
        if (fp != NULL) {
            fscanf(fp, "%lf", &scale);
            fscanf(fp, "%s", matrix);
            fclose(fp);

            /* build the response */
            snprintf(response, max_resp_len, "readView|%f|%s|",
                     scale, matrix);
            response_len = STRLEN(response);
        }

    /* "editor|...|" */
    } else if (strncmp(text, "editor|", 7) == 0) {
        snprintf(response, max_resp_len, "%s", text);
        response_len = STRLEN(response);

    }

    status = SUCCESS;

cleanup:
    FREE(vars_out);
    FREE(segs);
    FREE(cons);
    FREE(vars);
    FREE(begs);
    FREE(matrix);
    FREE(temp2);
    FREE(temp);
    FREE(entry);
    FREE(arg9);
    FREE(arg8);
    FREE(arg7);
    FREE(arg6);
    FREE(arg5);
    FREE(arg4);
    FREE(arg3);
    FREE(arg2);
    FREE(arg1);
    FREE(valu);
    FREE(type);
    FREE(name);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   sizeCallbackFromOpenCSM - change size of a DESPMTR                */
/*                                                                     */
/***********************************************************************/

void
sizeCallbackFromOpenCSM(void   *modl,   /* (in)  pointer to MODL */
                        int    ipmtr,   /* (in)  Parameter index (bias-1) */
                        int    nrow,    /* (in)  new number of rows */
                        int    ncol)    /* (in)  new number of columns */
{

    int    status=SUCCESS, irow, icol, index;

    double value, dot;
    char   *entry=NULL;

    modl_T *MODL = (modl_T*)modl;

    ROUTINE(sizeCallbackFromOpenCSM);

    /* --------------------------------------------------------------- */

    /* check that Modl is not NULL */
    if (MODL == NULL) {
        goto cleanup;
    }

    if (ipmtr >= 1 && ipmtr <= MODL->npmtr) {
        SPRINT3(2, "Size of %s changed to (%d,%d)", MODL->pmtr[ipmtr].name, nrow, ncol);
    }

    /* no further processor if webViewer is no up yet */
    if (cntxt == NULL) {
        SPRINT0(2, "WebViewer not up yet");
        goto cleanup;
    }

    MALLOC(entry, char, MAX_STR_LEN);

    /* build the response in JSON format */
    snprintf(response, max_resp_len, "getPmtrs|[");
    response_len = STRLEN(response);

    /* constant Parameters first */
    if (MODL != NULL) {
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (MODL->pmtr[ipmtr].type != OCSM_CONPMTR) continue;

            if (strlen(response) > 10) {
                addToResponse(",");
            }

            snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                     MODL->pmtr[ipmtr].name,
                     MODL->pmtr[ipmtr].type,
                     MODL->pmtr[ipmtr].nrow,
                     MODL->pmtr[ipmtr].ncol);
            addToResponse(entry);

            index = 0;
            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                    } else {
                        snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                    }
                    addToResponse(entry);
                    }
            }

            index = 0;
            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].dot[index++]);
                    } else {
                        snprintf(entry, MAX_STR_LEN, "%lg]", MODL->pmtr[ipmtr].dot[index++]);
                    }
                    addToResponse(entry);
                }
            }

            addToResponse("}");
        }

        /* external and configuration Parameters second */
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (MODL->pmtr[ipmtr].type != OCSM_DESPMTR &&
                MODL->pmtr[ipmtr].type != OCSM_CFGPMTR   ) continue;

            if (strlen(response) > 10) {
                addToResponse(",");
            }

            snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                     MODL->pmtr[ipmtr].name,
                     MODL->pmtr[ipmtr].type,
                     MODL->pmtr[ipmtr].nrow,
                     MODL->pmtr[ipmtr].ncol);
            addToResponse(entry);

            index = 0;
            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                    } else {
                        snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                    }
                    addToResponse(entry);
                }
            }

            index = 0;
            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].dot[index++]);
                    } else {
                        snprintf(entry, MAX_STR_LEN, "%lg]", MODL->pmtr[ipmtr].dot[index++]);
                    }
                    addToResponse(entry);
                }
            }

            addToResponse("}");
        }

        /* internal Parameters last */
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (MODL->pmtr[ipmtr].type != OCSM_LOCALVAR &&
                MODL->pmtr[ipmtr].type != OCSM_OUTPMTR  ) continue;

            /* make sure the mass properties have been set up */
            if (MODL->pmtr[ipmtr].name[0] == '@') {
                status = ocsmGetValu(MODL, ipmtr, 1, 1, &value, &dot);
                if (status != SUCCESS) {
                    SPRINT2(0, "ERROR:: ocsmGetValu(%s) detected %s", MODL->pmtr[ipmtr].name, ocsmGetText(status));
                }
            }

            /* skip if string-valued */
            if (MODL->pmtr[ipmtr].nrow == 0 || MODL->pmtr[ipmtr].ncol == 0) {
                continue;
            }

            if (strlen(response) > 10) {
                addToResponse(",");
            }

            snprintf(entry, MAX_STR_LEN, "{\"name\":\"%s\",\"type\":%d,\"nrow\":%d,\"ncol\":%d,\"value\":[",
                     MODL->pmtr[ipmtr].name,
                     MODL->pmtr[ipmtr].type,
                     MODL->pmtr[ipmtr].nrow,
                     MODL->pmtr[ipmtr].ncol);
            addToResponse(entry);

            index = 0;
            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "%lg,", MODL->pmtr[ipmtr].value[index++]);
                    } else {
                        snprintf(entry, MAX_STR_LEN, "%lg],\"dot\":[", MODL->pmtr[ipmtr].value[index++]);
                    }
                    addToResponse(entry);
                }
            }

            for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                    if (irow < MODL->pmtr[ipmtr].nrow || icol < MODL->pmtr[ipmtr].ncol) {
                        snprintf(entry, MAX_STR_LEN, "0,");
                    } else {
                        snprintf(entry, MAX_STR_LEN, "0]");
                    }
                    addToResponse(entry);
                }
            }

            addToResponse("}");
        }
        addToResponse("]");
    }

    /* send to browsers */
    SPRINT1(2, "\n<<< server2browser: %s", response);
    wv_broadcastText(response);

    response[0]  = '\0';
    response_len = 0;

cleanup:
    FREE(entry);

    if (status != SUCCESS) {
        SPRINT1(1, "sizeCallbackFromOpenCSM -> status=%d", status);
    }
}


/***********************************************************************/
/*                                                                     */
/*   spec_col - return color for a given scalar value                  */
/*                                                                     */
/***********************************************************************/

static void
spec_col(float  scalar,
         float  color[])
{
    int   indx;
    float frac;

    /* --------------------------------------------------------------- */

    if (lims[0] == lims[1]) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
    } else if (scalar <= lims[0]) {
        color[0] = color_map[0];
        color[1] = color_map[1];
        color[2] = color_map[2];
    } else if (scalar >= lims[1]) {
        color[0] = color_map[3*255  ];
        color[1] = color_map[3*255+1];
        color[2] = color_map[3*255+2];
    } else {
        frac  = 255.0 * (scalar - lims[0]) / (lims[1] - lims[0]);
        if (frac < 0  ) frac = 0;
        if (frac > 255) frac = 255;
        indx  = frac;
        frac -= indx;
        if (indx == 255) {
            indx--;
            frac += 1.0;
        }

        color[0] = frac * color_map[3*(indx+1)  ] + (1.0-frac) * color_map[3*indx  ];
        color[1] = frac * color_map[3*(indx+1)+1] + (1.0-frac) * color_map[3*indx+1];
        color[2] = frac * color_map[3*(indx+1)+2] + (1.0-frac) * color_map[3*indx+2];
    }
}


/***********************************************************************/
/*                                                                     */
/*   storeUndo - store an undo for the current command                 */
/*                                                                     */
/***********************************************************************/

static int
storeUndo(char   cmd[],                 /* (in)  current command */
          char   arg[])                 /* (in)  current argument */
{
    int       status = SUCCESS;         /* return status */

    int       iundo;
    char      *text=NULL;

    modl_T    *MODL = (modl_T*)modl;

    ROUTINE(storeUndo);

    /* --------------------------------------------------------------- */

    MALLOC(text, char, MAX_EXPR_LEN);

    /* if the undos are full, discard the most ancient one */
    if (nundo >= MAX_UNDOS) {
        status = ocsmFree(undo_modl[0]);
        CHECK_STATUS(ocsmFree);

        for (iundo = 0; iundo < nundo; iundo++) {
            undo_modl[iundo] = undo_modl[iundo+1];
            (void) STRNCPY(undo_text[iundo], undo_text[iundo+1], MAX_NAME_LEN-1);
        }

        nundo--;
    }

    /* store an undo snapshot */
    snprintf(text, MAX_EXPR_LEN, "%s %s", cmd, arg);
    STRNCPY(undo_text[nundo], text, 31);

    status = ocsmCopy(MODL, &undo_modl[nundo]);
    CHECK_STATUS(ocsmCopy);

    nundo++;

    SPRINT2(1, "~~> ocsmCopy() -> status=%d  (nundo=%d)",
            status, nundo);

cleanup:
    FREE(text);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   updateModl - update Bodys and mark Branches as dirty from prev MODL */
/*                                                                     */
/***********************************************************************/

static int
updateModl(modl_T *src_MODL,            /* (in)   source MODL */
           modl_T *tgt_MODL)            /* (both) target MODL */
{
    int       status = SUCCESS;         /* return status */

    int       ibrch, iattr;

    ROUTINE(updateModl);

    /* --------------------------------------------------------------- */

    /* move the Body info from the src_MODL into the
       new tgt_MODL so that recycling might happen */
    tgt_MODL->nbody = src_MODL->nbody;
    tgt_MODL->mbody = src_MODL->mbody;
    tgt_MODL->body  = src_MODL->body;

    src_MODL->nbody = 0;
    src_MODL->mbody = 0;
    src_MODL->body  = NULL;

    /* use src_MODL's context in tgt_MODL */
    if (tgt_MODL->context != NULL) {
        status = EG_close(tgt_MODL->context);
        CHECK_STATUS(EG_close);
    }

    tgt_MODL->context = src_MODL->context;

    /* starting at the beginning of the Branches, mark the Branch
       as  dirty if its type, name, arguements (except tmp_OpenCSM_),
       or attributes differ from the src_MODL */
    for (ibrch = 1; ibrch <= tgt_MODL->nbrch; ibrch++) {
        tgt_MODL->brch[ibrch].dirty = 0;

        if (ibrch > src_MODL->nbrch) {
            tgt_MODL->brch[ibrch].dirty = 1;
            goto recycle_message;
        } else if (tgt_MODL->brch[ibrch].type != src_MODL->brch[ibrch].type) {
            tgt_MODL->brch[ibrch].dirty = 1;
            goto recycle_message;
        } else if (strcmp(tgt_MODL->brch[ibrch].name, src_MODL->brch[ibrch].name) != 0) {
            tgt_MODL->brch[ibrch].dirty = 1;
            goto recycle_message;
        } else if (tgt_MODL->brch[ibrch].narg != src_MODL->brch[ibrch].narg) {
            tgt_MODL->brch[ibrch].dirty = 1;
            goto recycle_message;
        }

        if (tgt_MODL->brch[ibrch].narg >= 1) {
            if (strcmp(tgt_MODL->brch[ibrch].arg1, src_MODL->brch[ibrch].arg1) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 2) {
            if (strcmp(tgt_MODL->brch[ibrch].arg2, src_MODL->brch[ibrch].arg2) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 3) {
            if        (tgt_MODL->brch[ibrch].type == OCSM_UDPARG && strncmp(tgt_MODL->brch[ibrch].arg3, "tmp_OpenCSM_", 12) != 0) {
            } else if (tgt_MODL->brch[ibrch].type == OCSM_UDPRIM && strncmp(tgt_MODL->brch[ibrch].arg3, "tmp_OpenCSM_", 12) != 0) {
            } else if (strcmp(tgt_MODL->brch[ibrch].arg3, src_MODL->brch[ibrch].arg3) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 4) {
            if (strcmp(tgt_MODL->brch[ibrch].arg4, src_MODL->brch[ibrch].arg4) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 5) {
            if        (tgt_MODL->brch[ibrch].type == OCSM_UDPARG && strncmp(tgt_MODL->brch[ibrch].arg5, "tmp_OpenCSM_", 12) != 0) {
            } else if (tgt_MODL->brch[ibrch].type == OCSM_UDPRIM && strncmp(tgt_MODL->brch[ibrch].arg5, "tmp_OpenCSM_", 12) != 0) {
            } else if (strcmp(tgt_MODL->brch[ibrch].arg5, src_MODL->brch[ibrch].arg5) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 6) {
            if (strcmp(tgt_MODL->brch[ibrch].arg6, src_MODL->brch[ibrch].arg6) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 7) {
            if        (tgt_MODL->brch[ibrch].type == OCSM_UDPARG && strncmp(tgt_MODL->brch[ibrch].arg7, "tmp_OpenCSM_", 12) != 0) {
            } else if (tgt_MODL->brch[ibrch].type == OCSM_UDPRIM && strncmp(tgt_MODL->brch[ibrch].arg7, "tmp_OpenCSM_", 12) != 0) {
            } else if (strcmp(tgt_MODL->brch[ibrch].arg7, src_MODL->brch[ibrch].arg7) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 8) {
            if (strcmp(tgt_MODL->brch[ibrch].arg8, src_MODL->brch[ibrch].arg8) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }
        if (tgt_MODL->brch[ibrch].narg >= 9) {
            if        (tgt_MODL->brch[ibrch].type == OCSM_UDPARG && strncmp(tgt_MODL->brch[ibrch].arg9, "tmp_OpenCSM_", 12) != 0) {
            } else if (tgt_MODL->brch[ibrch].type == OCSM_UDPRIM && strncmp(tgt_MODL->brch[ibrch].arg9, "tmp_OpenCSM_", 12) != 0) {
            } else if (strcmp(tgt_MODL->brch[ibrch].arg9, src_MODL->brch[ibrch].arg9) != 0) {
                tgt_MODL->brch[ibrch].dirty = 1;
                goto recycle_message;
            }
        }

        if (tgt_MODL->brch[ibrch].nattr != src_MODL->brch[ibrch].nattr) {
            tgt_MODL->brch[ibrch].dirty = 1;
        } else {
            for (iattr = 0; iattr < tgt_MODL->brch[ibrch].nattr; iattr++) {
                if        (strcmp(tgt_MODL->brch[ibrch].attr[iattr].name, src_MODL->brch[ibrch].attr[iattr].name) != 0) {
                    tgt_MODL->brch[ibrch].dirty = 1;
                    goto recycle_message;
                } else if (strcmp(tgt_MODL->brch[ibrch].attr[iattr].defn, src_MODL->brch[ibrch].attr[iattr].defn) != 0) {
                    tgt_MODL->brch[ibrch].dirty = 1;
                    goto recycle_message;
                } else if (tgt_MODL->brch[ibrch].attr[iattr].type != src_MODL->brch[ibrch].attr[iattr].type) {
                    tgt_MODL->brch[ibrch].dirty = 1;
                    goto recycle_message;
                }
            }
        }

    recycle_message:
        if (tgt_MODL->brch[ibrch].dirty > 0) {
            SPRINT1(1, "    recycling disabled starting at Branch %d because of file differences", ibrch);
            break;
        }
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   addToHistogram - add entry to histogram                           */
/*                                                                     */
/***********************************************************************/

static int
addToHistogram(double entry,            /* (in)  entry to add */
               int    nhist,            /* (in)  number of histogram entries */
               double dhist[],          /* (in)  histogram entries */
               int    hist[])           /* (both)histogram counts */
{
    int    status = SUCCESS;            /* return status */

    int    ileft, imidl, irite;

    ROUTINE(addToHistogram);

    /* --------------------------------------------------------------- */

    /* binary search for correct histogram bin */
    ileft = 0;
    irite = nhist-1;

    while (irite-ileft > 1) {
        imidl = (ileft + irite) / 2;
        if (entry > dhist[imidl]) {
            ileft = imidl;
        } else {
            irite = imidl;
        }
    }

    hist[ileft]++;

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   printHistogram - print a histogram                                */
/*                                                                     */
/***********************************************************************/

static int
printHistogram(int    nhist,            /* (in)  number of histogram entries */
               double dhist[],          /* (in)  histogram entries */
               int    hist[])           /* (in)  histogram counts */
{
    int    status = SUCCESS;            /* return status */

    int    ihist, ntotal, ix;
    double percent;

    ROUTINE(printHistogram);

    /* --------------------------------------------------------------- */

    /* compute and print total entries in histogram */
    ntotal = 0;
    for (ihist = 0; ihist < nhist; ihist++) {
        ntotal += hist[ihist];
    }

    /* print basic histogram */
    percent = 100.0 * (double)(hist[0]) / (double)(ntotal);
    SPRINT3x(1, "    %9d (%5.1f%%)                    < %8.1e   |",
            hist[0], percent, dhist[1]);
    for (ix = 0; ix < 20; ix++) {
        if (5.0*ix >= percent) break;

        if (ix%5 == 4) {
            SPRINT0x(1, "+");
        } else {
            SPRINT0x(1, "-");
        }
    }
    SPRINT0(1, " ");

    for (ihist = 1; ihist < nhist-2; ihist++) {
        percent = 100.0 * (double)(hist[ihist]) / (double)(ntotal);
        SPRINT4x(1, "    %9d (%5.1f%%) between %8.1e and %8.1e   |",
                hist[ihist], percent, dhist[ihist], dhist[ihist+1]);
        for (ix = 0; ix < 20; ix++) {
            if (5.0*ix >= percent) break;

            if (ix%5 == 4) {
                SPRINT0x(1, "+");
            } else {
                SPRINT0x(1, "-");
            }
        }
        SPRINT0(1, " ");
    }

    percent = 100.0 * (double)(hist[nhist-2]) / (double)(ntotal);
    SPRINT3x(1, "    %9d (%5.1f%%)       > %8.1e                |",
            hist[nhist-2], percent, dhist[nhist-2]);
    for (ix = 0; ix < 20; ix++) {
        if (5.0*ix >= percent) break;

        if (ix%5 == 4) {
            SPRINT0x(1, "+");
        } else {
            SPRINT0x(1, "-");
        }
    }
    SPRINT0(1, " ");

    SPRINT1(1, "    %9d total", ntotal);

//cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   computeMassProps - compute mass properties discretely             */
/*                                                                     */
/***********************************************************************/

static int
computeMassProps(modl_T *MODL,          /* (in)  pointer to MODL */
                 int    ibody,          /* (in)  body index (bias-1) */
                 double props[])        /* (out) volume, area/length, xcg, ycg, zcg */
{
    int    status = SUCCESS;            /* return status */

    int      iedge, iface, npnt, ipnt, ntri, itri, ip0, ip1, ip2, nnode;
    CINT     *ptype, *pindx, *tris, *tric;
    double   len=0, area=0, vol=0, xcg=0, ycg=0, zcg=0;
    double   Ixx=0, Ixy=0, Ixz=0, Iyy=0, Iyz=0, Izz=0;
    double   data1[18], len1, area1, xa, ya, za, xb, yb, zb, xbar, ybar, zbar, areax, areay, areaz;
    CDOUBLE  *xyz, *uv;
    ego      *enodes;

    ROUTINE(computeMassProps);

    /* --------------------------------------------------------------- */

    /* for a NodeBody, just look up coordinates */
    if        (MODL->body[ibody].botype == OCSM_NODE_BODY) {
        status = EG_getBodyTopos(MODL->body[ibody].ebody, NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_evaluate(enodes[0], NULL, data1);
        CHECK_STATUS(EG_evaluate);

        EG_free(enodes);

        xcg  = data1[0];
        ycg  = data1[1];
        zcg  = data1[2];

    /* for a WireBody, accumulate statistics on each segment of each Edge */
    } else if (MODL->body[ibody].botype == OCSM_WIRE_BODY) {
        if (MODL->body[ibody].etess == NULL) {
            status = OCSM_NEED_TESSELLATION;
            goto cleanup;
        }

        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt, &xyz, &uv);

            for (ipnt = 1; ipnt < npnt; ipnt++) {
                ip0 = 3 * (ipnt - 1);
                ip1 = 3 * (ipnt    );

                len1 = sqrt((xyz[ip1  ]-xyz[ip0  ]) * (xyz[ip1  ]-xyz[ip0  ])
                            +(xyz[ip1+1]-xyz[ip0+1]) * (xyz[ip1+1]-xyz[ip0+1])
                            +(xyz[ip1+2]-xyz[ip0+2]) * (xyz[ip1+2]-xyz[ip0+2]));

                len += len1;
                xcg += (xyz[ip1  ] + xyz[ip0  ]) * len1 / 2;
                ycg += (xyz[ip1+1] + xyz[ip0+1]) * len1 / 2;
                zcg += (xyz[ip1+2] + xyz[ip0+2]) * len1 / 2;
            }
        }

        xcg /= len;
        ycg /= len;
        zcg /= len;

        area = len;     // store length in area for printout below

    /* for a SheetBody, accumulate statistics on each triangle of each face */
    } else if (MODL->body[ibody].botype == OCSM_SHEET_BODY) {
        if (MODL->body[ibody].etess == NULL) {
            status = OCSM_NEED_TESSELLATION;
            goto cleanup;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);

            for (itri = 0; itri < ntri; itri++) {
                ip0 = 3 * (tris[3*itri  ] - 1);
                ip1 = 3 * (tris[3*itri+1] - 1);
                ip2 = 3 * (tris[3*itri+2] - 1);

                xa = xyz[ip1  ] - xyz[ip0  ];
                ya = xyz[ip1+1] - xyz[ip0+1];
                za = xyz[ip1+2] - xyz[ip0+2];

                xb = xyz[ip2  ] - xyz[ip0  ];
                yb = xyz[ip2+1] - xyz[ip0+1];
                zb = xyz[ip2+2] - xyz[ip0+2];

                xbar = xyz[ip0  ] + xyz[ip1  ] + xyz[ip2  ];
                ybar = xyz[ip0+1] + xyz[ip1+1] + xyz[ip2+1];
                zbar = xyz[ip0+2] + xyz[ip1+2] + xyz[ip2+2];

                areax = ya * zb - za * yb;
                areay = za * xb - xa * zb;
                areaz = xa * yb - ya * xb;
                area1 = sqrt(areax*areax + areay*areay + areaz*areaz) / 2;

                area += area1;
                xcg  += (xbar * area1) / 3;
                ycg  += (ybar * area1) / 3;
                zcg  += (zbar * area1) / 3;
            }
        }

        xcg /= area;
        ycg /= area;
        zcg /= area;

    /* for a SolidBody, use the divergence theorem on each triangle
       of each Face to get the volume statistics */
    } else if (MODL->body[ibody].botype == OCSM_SOLID_BODY) {
        if (MODL->body[ibody].etess == NULL) {
            status = OCSM_NEED_TESSELLATION;
            goto cleanup;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);

            for (itri = 0; itri < ntri; itri++) {
                ip0 = 3 * (tris[3*itri  ] - 1);
                ip1 = 3 * (tris[3*itri+1] - 1);
                ip2 = 3 * (tris[3*itri+2] - 1);

                xa = xyz[ip1  ] - xyz[ip0  ];
                ya = xyz[ip1+1] - xyz[ip0+1];
                za = xyz[ip1+2] - xyz[ip0+2];

                xb = xyz[ip2  ] - xyz[ip0  ];
                yb = xyz[ip2+1] - xyz[ip0+1];
                zb = xyz[ip2+2] - xyz[ip0+2];

                xbar = xyz[ip0  ] + xyz[ip1  ] + xyz[ip2  ];
                ybar = xyz[ip0+1] + xyz[ip1+1] + xyz[ip2+1];
                zbar = xyz[ip0+2] + xyz[ip1+2] + xyz[ip2+2];

                areax = ya * zb - za * yb;
                areay = za * xb - xa * zb;
                areaz = xa * yb - ya * xb;

                area += sqrt(areax*areax + areay*areay + areaz*areaz) / 2;
                vol  += (                xbar * areax +                 ybar * areay +                 zbar * areaz) /  18;

                xcg  += (xbar/2        * xbar * areax + xbar          * ybar * areay + xbar          * zbar * areaz) /  54;
                ycg  += (ybar          * xbar * areax + ybar/2        * ybar * areay + ybar          * zbar * areaz) /  54;
                zcg  += (zbar          * xbar * areax + zbar          * ybar * areay + zbar/2        * zbar * areaz) /  54;

                Ixx  += (                               ybar   * ybar * ybar * areay + zbar   * zbar * zbar * areaz) / 162;
                Iyy  += (xbar   * xbar * xbar * areax                                + zbar   * zbar * zbar * areaz) / 162;
                Izz  += (xbar   * xbar * xbar * areax + ybar   * ybar * ybar * areay                               ) / 162;

                Ixy  -= (xbar/2 * ybar * xbar * areax + ybar/2 * xbar * ybar * areay + xbar   * ybar * zbar * areaz) / 162;
                Ixz  -= (xbar/2 * zbar * xbar * areax + xbar   * zbar * ybar * areay + zbar/2 * xbar * zbar * areaz) / 162;
                Iyz  -= (ybar   * zbar * xbar * areax + ybar/2 * zbar * ybar * areay + zbar/2 * ybar * zbar * areaz) / 162;
            }
        }

        xcg /= vol;
        ycg /= vol;
        zcg /= vol;

        /* parallel-axis theorem */
        Ixx -= vol * (            ycg * ycg + zcg * zcg);
        Iyy -= vol * (xcg * xcg             + zcg * zcg);
        Izz -= vol * (xcg * xcg + ycg * ycg            );

        Ixy += vol * (xcg * ycg      );
        Ixz += vol * (xcg *       zcg);
        Iyz += vol * (      ycg * zcg);
    }

    /* store the properties in the array */
    props[ 0] = vol;
    props[ 1] = area;
    props[ 2] = xcg;
    props[ 3] = ycg;
    props[ 4] = zcg;
    props[ 5] = Ixx;
    props[ 6] = Ixy;
    props[ 7] = Ixz;
    props[ 8] = Ixy;
    props[ 9] = Iyy;
    props[10] = Iyz;
    props[11] = Ixz;
    props[12] = Iyz;
    props[13] = Izz;

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   checkForGanged - check for possible ganged Boolean operations     */
/*                                                                     */
/***********************************************************************/

static int
checkForGanged(modl_T *MODL)            /* (in)  pointer to MODL */
{
    int    status = SUCCESS;            /* return status */

    int    ibody, jbody, nlist, i, *list=NULL;

    ROUTINE(checkForGanged);

    /* --------------------------------------------------------------- */

    SPRINT0(1, "\nChecking for opportunity for ganged Boolean operations...");

    MALLOC(list, int, MODL->nbody);

    /* look for possible ganged SUBTRACTs */
    for (ibody = MODL->nbody; ibody >= 1; ibody--) {

        nlist = 0;
        jbody     = ibody;
        while (MODL->body[jbody].brtype == OCSM_SUBTRACT) {
            list[nlist++] = jbody;
            jbody = MODL->body[jbody].ileft;
            if (jbody <= 0) break;
        }

        if (nlist > 1) {
            SPRINT0x(1, "   possible ganged SUBTRACTs that created Bodys:");
            for (i = nlist-1; i >= 0; i--) {
                SPRINT1x(1, " %d", list[i]);
            }
            SPRINT0(1, " ");

            /* move beyond this ganged operation */
            ibody = list[nlist-1];
        }
    }

    /* look for possible ganged UNIONs */
    for (ibody = MODL->nbody; ibody >= 1; ibody--) {

        nlist = 0;
        jbody     = ibody;
        while (MODL->body[jbody].brtype == OCSM_UNION) {
            list[nlist++] = jbody;
            jbody = MODL->body[jbody].ileft;
            if (jbody <= 0) break;
        }

        if (nlist > 1) {
            SPRINT0x(1, "   possible ganged UNIONs that created Bodys:");
            for (i = nlist-1; i >= 0; i--) {
                SPRINT1x(1, " %d", list[i]);
            }
            SPRINT0(1, " ");

            /* move beyond this ganged operation */
            ibody = list[nlist-1];
        }
    }

cleanup:
    FREE(list);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   checkParallel - compute parallalizability of case                 */
/*                                                                     */
/***********************************************************************/

/* maximum number of Bodys per processor */
#define  MBODY  100

typedef struct {
    int       iprnt;                    /* parent processor */
    int       ibody;                    /* index within bodys[] for parent */
    double    CPU;                      /* total CPU used by processor */
    int       nbody;                    /* number of Bodys */
    int       bodys[MBODY];             /* list of Bodys */
} proc_T;


static int
checkParallel(modl_T *MODL)             /* (in)  pointer to MODL */
{
    int    status = SUCCESS;            /* return status */

    int    nproc, iproc, jproc, ibody, jbody;
    int    itype, ichld, ileft, jleft, irite, jrite;
    double CPUttl, CPUmax;
    proc_T *proc=NULL;

    ROUTINE(checkParallel);

    /* --------------------------------------------------------------- */

    SPRINT0(1, "\nOpportunity for parallelism by executing feature tree in parallel");

    /* get an overly-large table of processors */
    MALLOC(proc, proc_T, MODL->nbody);

    /* initialize the processors */
    for (iproc = 0; iproc < MODL->nbody; iproc++) {
        proc[iproc].iprnt = -1;
        proc[iproc].ibody = -1;
        proc[iproc].CPU   =  0;
        proc[iproc].nbody =  0;
    }

    nproc = 0;

    /* assign each Body/task to a processor */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        itype = MODL->body[ibody].brtype;
        ichld = MODL->body[ibody].ichld;
        irite = MODL->body[ibody].irite;
        ileft = MODL->body[ibody].ileft;

        /* skip UDPARG statements */
        if (MODL->nbody > 1 && ichld == 0 && ileft == -1 && irite == -1) continue;

        /* if ileft is -1, this is top level, so add a new processor */
        if (ileft == -1) {
            proc[nproc].bodys[0] = ibody;
            proc[nproc].CPU      = MODL->body[ibody].CPU;
            proc[nproc].nbody++;
            nproc++;

       /* if this is a RESTORE, add a new processor */
        } else if (irite == -1 && itype == OCSM_RESTORE) {
            proc[nproc].bodys[0] = ibody;
            proc[nproc].iprnt    = -1;
            proc[nproc].ibody    = -1;
            proc[nproc].nbody++;

            for (jproc = 0; jproc < nproc; jproc++) {
                for (jbody = 0; jbody < proc[jproc].nbody; jbody++) {
                    if (proc[jproc].bodys[jbody] == ileft) {
                        proc[nproc].iprnt = jproc;
                        proc[nproc].ibody = jbody;
                        break;
                    }
                }
                if (proc[nproc].iprnt >= 0) break;
            }

            /* find time needed to get to the STORE */
            proc[nproc].CPU = MODL->body[ibody].CPU;
            jproc = proc[nproc].iprnt;
            jbody = proc[nproc].ibody;
            while (jproc >= 0) {
                while (jbody >= 0) {
                    proc[nproc].CPU += MODL->body[proc[jproc].bodys[jbody]].CPU;
                    jbody--;
                }
                jbody = proc[jproc].ibody;
                jproc = proc[jproc].iprnt;
            }

            nproc++;

        /* otherwise, find the entry in the tree that matches its one parent */
        } else if (irite == -1) {
            jleft = -1;
            for (jproc = 0; jproc < nproc; jproc++) {
                if (proc[jproc].bodys[proc[jproc].nbody-1] == ileft) {
                    jleft = jproc;
                }
            }

            assert(jleft >= 0);

            /* add to jproc */
            if (proc[jleft].nbody < MBODY-1) {
                proc[jleft].bodys[proc[jleft].nbody] = ibody;
                proc[jleft].CPU                     += MODL->body[ibody].CPU;
                proc[jleft].nbody++;

            /* create a new processor */
            } else {
                proc[nproc].bodys[0] = ibody;
                proc[nproc].CPU      = proc[jleft].CPU + MODL->body[ibody].CPU;
                proc[nproc].iprnt    = jleft;
                proc[nproc].ibody    = proc[jleft].nbody - 1;
                proc[nproc].nbody++;
                nproc++;
            }

        /* otherwise, find the two parents and add body to proc with larger CPU */
        } else {
            jleft = -1;
            jrite = -1;
            for (jproc = 0; jproc < nproc; jproc++) {
                if        (proc[jproc].bodys[proc[jproc].nbody-1] == ileft) {
                    jleft = jproc;
                } else if (proc[jproc].bodys[proc[jproc].nbody-1] == irite) {
                    jrite = jproc;
                }
            }

            assert(jleft >= 0);
            assert(jrite >= 0);

            if (proc[jleft].CPU >= proc[jrite].CPU) {
                proc[jleft].bodys[proc[jleft].nbody] = ibody;
                proc[jleft].CPU                     += MODL->body[ibody].CPU;
                proc[jleft].nbody++;
            } else {
                proc[jrite].bodys[proc[jrite].nbody] = ibody;
                proc[jrite].CPU                     += MODL->body[ibody].CPU;
                proc[jrite].nbody++;
            }
        }
    }

    /* print processor summary */
    for (iproc = 0; iproc < nproc; iproc++) {
        SPRINT4x(1, "proc %3d: prnt=%3d:%-3d CPU=%9.4f,", iproc, proc[iproc].iprnt, proc[iproc].ibody, proc[iproc].CPU);
        for (ibody = 0; ibody < proc[iproc].nbody; ibody++) {
            SPRINT1x(1, " %3d", proc[iproc].bodys[ibody]);
        }
        SPRINT0(1, " ");
    }

    /* print time summary */
    CPUmax = 0;
    for (iproc = 0; iproc < nproc; iproc++) {
        CPUmax = MAX(CPUmax, proc[iproc].CPU);
    }

    CPUttl = 0;
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        CPUttl += MODL->body[ibody].CPU;
    }

    SPRINT1(1, "total CPU=%9.4f", CPUttl);
    SPRINT1(1, "max   CPU=%9.4f", CPUmax);
    SPRINT1(1, "ttl/max  =%9.4f", CPUttl/CPUmax);

cleanup:
    FREE(proc);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plugsMain - main routine for PLUGS                                */
/*                                                                     */
/***********************************************************************/

static int
plugsMain(modl_T *MODL,                 /* (in)  pointer to MODL */
          int    npass,                 /* (in)  number of phase2 passes */
          int    ncloud,                /* (in)  number of points in cloud */
          double cloud[])               /* (in)  array  of points in cloud */
{
    int    status = SUCCESS;            /* return status */

    int     npmtr, ipmtr, ibody, jbody;
    int     *pmtrindx=NULL;
    double  rms;
    char    egadsname[81];
    ego     emodel;
    clock_t old_time, new_time;
    FILE    *fp_plugs=NULL;

    ROUTINE(plugsMain);

    /* --------------------------------------------------------------- */

    if (ncloud <= 0) {
        SPRINT0(0, "ERROR:: there needs to be at least one point in the cloud");
        status = -999;
        goto cleanup;
    }

    SPRINT0(0, "\n===================================");
    SPRINT1(0,   "PLUGS with %d points in the cloud", ncloud);
    SPRINT0(0,   "===================================\n");

    /* get array that will certainly be big enough */
    MALLOC(pmtrindx, int, MODL->npmtr);

    /* make a list of the DESPMTRs that will be varied */
    npmtr = 0;
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
            if (MODL->pmtr[ipmtr].nrow != 1 ||
                MODL->pmtr[ipmtr].ncol != 1   ) {
                SPRINT3(0, "ERROR:: DESPMTR %s is (%d*%d) and must be a scalar",
                        MODL->pmtr[ipmtr].name, MODL->pmtr[ipmtr].nrow, MODL->pmtr[ipmtr].ncol);
                status = -999;
                goto cleanup;
            }
            pmtrindx[npmtr] = ipmtr;

            SPRINT3(1, "initial DESPMTR %3d: %20s = %10.5f",
                    npmtr, MODL->pmtr[pmtrindx[npmtr]].name, MODL->pmtr[ipmtr].value[0]);

            npmtr++;
        }
    }

    /* open plugs history file and write header */
    fp_plugs = fopen("plugs.hist", "w");

    if (fp_plugs != NULL) {
        fprintf(fp_plugs, " pass iter");
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            fprintf(fp_plugs, " %11s", MODL->pmtr[pmtrindx[ipmtr]].name);
        }
        fprintf(fp_plugs, "        rms\n");
        fflush( fp_plugs);

        fprintf(fp_plugs, " %4d %4d", -1, 0);
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            fprintf(fp_plugs, " %11.6f", MODL->pmtr[pmtrindx[ipmtr]].value[0]);
        }
        fprintf(fp_plugs, "\n");
    }

    if (PLUGS_CREATE_EGADS_FILES > 0) {
        status = EG_makeTopology(MODL->context, NULL, MODEL, 0, NULL,
                                 1, &(MODL->body[MODL->nbody].ebody), NULL, &emodel);
        CHECK_STATUS(EG_makeTopology);

        snprintf(egadsname, 80, "plugs_%d_%d.egads", -1, 0);

        (void) remove(egadsname);

        status = EG_saveModel(emodel, egadsname);
        CHECK_STATUS(EG_saveModel);

        status = EG_deleteObject(emodel);
        CHECK_STATUS(EG_deleteObject);
    }

    /* ensure that there is only one Body on the stack */
    ibody = -1;
    for (jbody = 1; jbody <= MODL->nbody; jbody++) {
        if (MODL->body[jbody].onstack == 1) {
            if (ibody < 0) {
                ibody = jbody;
            } else {
                SPRINT0(0, "ERROR:: there can only be one Body on the stack");
                status = -999;
                goto cleanup;
            }
        }
    }

    old_time = clock();

    /* change design parameters to get correct bounding box */
    if (npass >= 0) {
        SPRINT0(0, "\nPLUGS phase1: match bounding boxes\n");
        status = plugsPhase1(MODL, ibody, npmtr, pmtrindx, ncloud, cloud, &rms, fp_plugs);
        CHECK_STATUS(plugsPhase1);
    } else {
        SPRINT0(0, "\nPLUGS phase1: skipped\n");
    }

    /* change design parameters to match points in cloud */
    if (npass != 0) {
        SPRINT0(0, "\nPLUGS phase2: match cloud points");
        status = plugsPhase2(MODL, abs(npass), ibody, npmtr, pmtrindx, ncloud, cloud, &rms, fp_plugs);
        CHECK_STATUS(plugsPhase2);
    }

    new_time = clock();

    /* close the history file */
    if (fp_plugs != NULL) {
        fclose(fp_plugs);
        fp_plugs = NULL;
    }

    /* print final DESPMTRs and write a file that can be read by ocsmUpdateDespmtrs */
    fp_plugs = fopen("plugs.despmtrs", "w");

    if (fp_plugs != NULL) {
        npmtr = 0;
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
                SPRINT3(1, "final  DESPMTR %3d: %20s = %10.5f",
                        npmtr, MODL->pmtr[ipmtr].name, MODL->pmtr[ipmtr].value[0]);

                fprintf(fp_plugs, "%-32s %15.8f\n", MODL->pmtr[ipmtr].name, MODL->pmtr[ipmtr].value[0]);

                npmtr++;
            }
        }

        fclose(fp_plugs);
    }

    SPRINT1(1, "final  rms distance to cloud               %10.5f", rms);

    SPRINT1(1, "\n==> PLUGS total CPUtime=%9.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    /* tell user how to get configuration with updated DESPMTRs */
    SPRINT0(1, "\nHit \"Up to date\" to show results of PLUGS\n");

    /* reset the plugs flag so that verification will be reanabled */
    plugs = -1;

cleanup:
    FREE(pmtrindx);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plugsPhase1 - phase 1 of PLUGS                                    */
/*                                                                     */
/***********************************************************************/

static int
plugsPhase1(modl_T *MODL,               /* (in)  pointer to MODL */
            int    ibody,               /* (in)  Body index (bias-1) */
            int    npmtr,               /* (in)  number  of DESPMTRs */
            int    pmtrindx[],          /* (in)  indices of DESPMTRs */
            int    ncloud,              /* (in)  number of points in cloud */
            double cloud[],             /* (in)  array  of points in cloud */
            double *rmsbest,            /* (out) best rms distance to cloud */
  /*@null@*/FILE   *fp_plugs)           /* (in)  file pointer to history file */
{
    int    status = SUCCESS;            /* return status */

    int     icloud, ipmtr, jpmtr, inode, builtTo, nbody, ierr, iter, niter=20, idone, oldOutLevel;
    double  bbox_cloud[6], bbox_modl[6], vel[3], qerr[6], rms, lambda, value, dot, lbound, ubound, dmax;
    double  *ajac=NULL, *ajtj=NULL, *ajtq=NULL, *delta=NULL, *pmtrbest=NULL, *W=NULL;
    char    egadsname[81];
    ego     emodel;
    clock_t old_time, new_time;

    ROUTINE(plugsPhase1);

    /* --------------------------------------------------------------- */

    *rmsbest = 0;

    assert(ncloud>0);

    /* if there are no DESPMTRs, simply return */
    if (npmtr == 0) {
        SPRINT0(1, "Phase1 will be skipped because npmtr=0");
        goto cleanup;
    }

    old_time = clock();

    /* in phase1, the DESPMTRs are changed (as little as possible) to match
       the bounding box of the cloud.  hence there are 6 errors that
       are being minimized (xmin, ymin, zmin, xnax, ymax, zmax) */

    /* get needed arrays */
    MALLOC(ajac,     double,     6*npmtr);    // jacobian
    MALLOC(ajtj,     double, npmtr*npmtr);    // Jt * J
    MALLOC(ajtq,     double,       npmtr);    // Jt * Q
    MALLOC(delta,    double,       npmtr);    // proposed change in DESPMTRs
    MALLOC(pmtrbest, double,       npmtr);    // DESPMTRs associated with best
    MALLOC(W,        double,       npmtr);    // singular values

    /* find the bounding box of the cloud */
    bbox_cloud[0] = cloud[0];
    bbox_cloud[1] = cloud[1];
    bbox_cloud[2] = cloud[2];
    bbox_cloud[3] = cloud[0];
    bbox_cloud[4] = cloud[1];
    bbox_cloud[5] = cloud[2];

    for (icloud = 0; icloud < ncloud; icloud++) {
        if (cloud[3*icloud  ] < bbox_cloud[0]) bbox_cloud[0] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] < bbox_cloud[1]) bbox_cloud[1] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] < bbox_cloud[2]) bbox_cloud[2] = cloud[3*icloud+2];
        if (cloud[3*icloud  ] > bbox_cloud[3]) bbox_cloud[3] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] > bbox_cloud[4]) bbox_cloud[4] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] > bbox_cloud[5]) bbox_cloud[5] = cloud[3*icloud+2];
    }

    SPRINT3(1, "bbox of cloud: %10.5f %10.5f %10.5f",   bbox_cloud[0], bbox_cloud[1], bbox_cloud[2]);
    SPRINT3(1, "               %10.5f %10.5f %10.5f\n", bbox_cloud[3], bbox_cloud[4], bbox_cloud[5]);

    /* compute initial objective function */
    bbox_modl[0] = +HUGEQ;
    bbox_modl[1] = +HUGEQ;
    bbox_modl[2] = +HUGEQ;
    bbox_modl[3] = -HUGEQ;
    bbox_modl[4] = -HUGEQ;
    bbox_modl[5] = -HUGEQ;

    for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
        if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
            bbox_modl[0] = MODL->body[ibody].node[inode].x;
        }
        if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
            bbox_modl[1] = MODL->body[ibody].node[inode].y;
        }
        if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
            bbox_modl[2] = MODL->body[ibody].node[inode].z;
        }
        if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
            bbox_modl[3] = MODL->body[ibody].node[inode].x;
        }
        if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
            bbox_modl[4] = MODL->body[ibody].node[inode].y;
        }
        if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
            bbox_modl[5] = MODL->body[ibody].node[inode].z;
        }
    }

    SPRINT3(1, "bbox of MODL:  %10.5f %10.5f %10.5f",   bbox_modl[0], bbox_modl[1], bbox_modl[2]);
    SPRINT3(1, "               %10.5f %10.5f %10.5f\n", bbox_modl[3], bbox_modl[4], bbox_modl[5]);

    /* compute the initial errors */
    rms = 0;
    for (ierr = 0; ierr < 6; ierr++) {
        qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
        rms       += qerr[ierr] * qerr[ierr];
    }
    rms = sqrt(rms / 6);

    /* print the initial rms and DESPMTRs */
    SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", -1, rms);
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValu);

        SPRINT1x(1, " %10.5f", pmtrbest[ipmtr]);
    }
    SPRINT0(1, " ");

    /* initialize for levenberg-marquardt */
    *rmsbest = rms;
    lambda   = 1;

    /* levenberg-marquardt iterations */
    for (iter = 0; iter < niter; iter++) {

        /* find the bounding box of the MODL and its velocities */
        bbox_modl[0] = +HUGEQ;
        bbox_modl[1] = +HUGEQ;
        bbox_modl[2] = +HUGEQ;
        bbox_modl[3] = -HUGEQ;
        bbox_modl[4] = -HUGEQ;
        bbox_modl[5] = -HUGEQ;

        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmSetVelD(MODL, 0,               0, 0, 0.0);
            CHECK_STATUS(ocsmSetVelD);

            status = ocsmSetVelD(MODL, pmtrindx[ipmtr], 1, 1, 1.0);
            CHECK_STATUS(ocsmSetVelD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);
            CHECK_STATUS(ocsmBuild);

            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
                    bbox_modl[0] = MODL->body[ibody].node[inode].x;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[        ipmtr] = vel[0];
                }
                if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
                    bbox_modl[1] = MODL->body[ibody].node[inode].y;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[  npmtr+ipmtr] = vel[1];
                }
                if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
                    bbox_modl[2] = MODL->body[ibody].node[inode].z;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[2*npmtr+ipmtr] = vel[2];
                }
                if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
                    bbox_modl[3] = MODL->body[ibody].node[inode].x;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[3*npmtr+ipmtr] = vel[0];
                }
                if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
                    bbox_modl[4] = MODL->body[ibody].node[inode].y;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[4*npmtr+ipmtr] = vel[1];
                }
                if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
                    bbox_modl[5] = MODL->body[ibody].node[inode].z;

                    status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vel);
                    CHECK_STATUS(ocsmGetVel);

                    ajac[5*npmtr+ipmtr] = vel[2];
                }
            }
        }

        /* compute the errors */
        for (ierr = 0; ierr < 6; ierr++) {
            qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
        }

        /* compute Jt * J and Jt * Q */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                ajtj[ipmtr*npmtr+jpmtr] = 0;
                for (ierr = 0; ierr < 6; ierr++) {
                    ajtj[ipmtr*npmtr+jpmtr] += ajac[ierr*npmtr+ipmtr] * ajac[ierr*npmtr+jpmtr];
                }
            }
            ajtj[ipmtr*npmtr+ipmtr] *= (1.0 + lambda);
            ajtq[ipmtr] = 0;
            for (ierr = 0; ierr < 6; ierr++) {
                ajtq[ipmtr] -= qerr[ierr] * ajac[ierr*npmtr+ipmtr];
            }
        }

        /* solve for potential change.  this is done with SVD because the
           matrix will be singular for any DESPMTR that does not (currently)
           affect the errors */
        status = solsvd(ajtj, ajtq, npmtr, npmtr, W, delta);
        CHECK_STATUS(solsvd);

        /* if all the delta's are small, no sense in any more iterations */
        dmax = 0;
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            if (fabs(delta[ipmtr]) > dmax) dmax = fabs(delta[ipmtr]);
        }

        if (dmax < EPS06) {
            SPRINT0(1, "maximum delta is small, so no more iterations");
            goto cleanup;
        }

        /* temporarily apply change */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            value = pmtrbest[ipmtr] + delta[ipmtr];

            status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
            CHECK_STATUS(ocsmGetBnds);

            if (value < lbound) value = lbound;
            if (value > ubound) value = ubound;

            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);
        }

        /* rebuild and evaluate the new objective function */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);
        if (status < SUCCESS) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                SPRINT3(0, "error  DESPMTR %3d: %20s = %10.5f",
                        ipmtr, MODL->pmtr[pmtrindx[ipmtr]].name, MODL->pmtr[pmtrindx[ipmtr]].value[0]);
            }
        }
        CHECK_STATUS(ocsmBuild);

        bbox_modl[0] = +HUGEQ;
        bbox_modl[1] = +HUGEQ;
        bbox_modl[2] = +HUGEQ;
        bbox_modl[3] = -HUGEQ;
        bbox_modl[4] = -HUGEQ;
        bbox_modl[5] = -HUGEQ;

        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            if (MODL->body[ibody].node[inode].x <= bbox_modl[0]) {
                bbox_modl[0] = MODL->body[ibody].node[inode].x;
            }
            if (MODL->body[ibody].node[inode].y <= bbox_modl[1]) {
                bbox_modl[1] = MODL->body[ibody].node[inode].y;
            }
            if (MODL->body[ibody].node[inode].z <= bbox_modl[2]) {
                bbox_modl[2] = MODL->body[ibody].node[inode].z;
            }
            if (MODL->body[ibody].node[inode].x >= bbox_modl[3]) {
                bbox_modl[3] = MODL->body[ibody].node[inode].x;
            }
            if (MODL->body[ibody].node[inode].y >= bbox_modl[4]) {
                bbox_modl[4] = MODL->body[ibody].node[inode].y;
            }
            if (MODL->body[ibody].node[inode].z >= bbox_modl[5]) {
                bbox_modl[5] = MODL->body[ibody].node[inode].z;
            }
        }

        rms = 0;
        for (ierr = 0; ierr < 6; ierr++) {
            qerr[ierr] = bbox_modl[ierr] - bbox_cloud[ierr];
            rms       += qerr[ierr] * qerr[ierr];
        }
        rms = sqrt(rms / 6);

        SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", iter, rms);
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
            CHECK_STATUS(ocsmGetValu);

            SPRINT1x(1, " %10.5f", value);
        }

        /* store history */
        if (fp_plugs != NULL) {
            fprintf(fp_plugs, " %4d %4d", 0, iter);
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
                CHECK_STATUS(ocsmGetPmtr);

                fprintf(fp_plugs, " %11.6f", value);
            }
            if (rms < *rmsbest) {
                fprintf(fp_plugs, "%11.5f ACCEPT\n", rms);
            } else {
                fprintf(fp_plugs, "%11.5f\n", rms);
            }
            fflush(fp_plugs);
        }

        if (PLUGS_CREATE_EGADS_FILES > 0) {
            status = EG_makeTopology(MODL->context, NULL, MODEL, 0, NULL,
                                     1, &(MODL->body[MODL->nbody].ebody), NULL, &emodel);
            CHECK_STATUS(EG_makeTopology);

            snprintf(egadsname, 80, "plugs_%d_%d.egads", 0, iter);

            (void) remove(egadsname);

            status = EG_saveModel(emodel, egadsname);
            CHECK_STATUS(EG_saveModel);

            status = EG_deleteObject(emodel);
            CHECK_STATUS(EG_deleteObject);
        }

        /* if this new solution is better, accept it and decrease lambda
           (making it more newton-like) */
        if (rms < *rmsbest) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
                CHECK_STATUS(ocsmGetPmtr);
            }
            *rmsbest = rms;
            lambda   = MAX(1.0e-10, lambda/2);
            SPRINT1(1, "  accepted: lambda=%10.3e", lambda);

            /* check for convergence */
            idone = 1;
            for (ierr = 0; ierr < 6; ierr++) {
                if (fabs(qerr[ierr]) > EPS06) {
                    idone = 0;
                    break;
                }
            }
            if (idone == 1) {
                SPRINT0(1, "Phase 1 converged");
                break;
            }

        /* otherwise reject the step and increase lambda
           (making it more steppest-descent-like) */
        } else {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                CHECK_STATUS(ocsmSetValuD);
            }
            lambda = MIN(1.0e+10, lambda*2);
            SPRINT1(1, "  rejected: lambda=%10.3e", lambda);
        }
    }

    new_time = clock();
    SPRINT1(1, "Phase 1 CPUtime=%9.3f sec", (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

cleanup:
    FREE(ajac);
    FREE(ajtj);
    FREE(ajtq);
    FREE(delta);
    FREE(pmtrbest);
    FREE(W);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plugsPhase2 - phase 2 of PLUGS                                    */
/*                                                                     */
/***********************************************************************/

static int
plugsPhase2(modl_T *MODL,               /* (in)  pointer to MODL */
            int    npass,               /* (in)  number of passes */
            int    ibody,               /* (in)  Body index (bias-1) */
            int    npmtr,               /* (in)  number  of DESPMTRs */
            int    pmtrindx[],          /* (in)  indices of DESPMTRs */
            int    ncloud,              /* (in)  number of points in cloud */
            double cloud[],             /* (in)  array  of points in cloud */
            double *rmsbest,            /* (out) best rms distance to cloud */
  /*@null@*/FILE   *fp_plugs)           /* (in)  FILE pointer to history file */
{
    int    status = SUCCESS;            /* return status */

    int     iface, icloud, ipmtr, jpmtr, npnt, ntri, itri, ip0, ip1, ip2;
    int     nerr, nvar, ivar, count, unclass, reclass, ireclass, ipass, restart, iter, niter=50;
    int     nbody, builtTo, oldOutLevel, periodic, ibest, scaleDiag, naccept, iperturb, limit;
    int     *face=NULL, *prevface=NULL;
    CINT    *ptype, *pindx, *tris, *tric;
    double  bbox_cloud[6], dtest, value, dot, lbound, ubound, rms, data[18], rmsperturb, valperturb;
    double  dmax, lambda, uvrange[4], dbest, scaleFact, uv_guess[2], xyz_guess[3];
    double  atotal, massprops[14];
    double  *dist=NULL, *uvface=NULL, *velface=NULL;
    double  *beta=NULL, *delta=NULL, *qerr=NULL, *qerrbest=NULL, *ajac=NULL;
    double  *atri=NULL, *btri=NULL, *ctri=NULL, *dtri=NULL, *xtri=NULL;
    double  *mat=NULL, *rhs=NULL, *xxx=NULL;
    double  *pmtrbest=NULL, *pmtrlast=NULL;
    CDOUBLE *xyz, *uv;
    char    egadsname[81];
    ego     emodel;
    clock_t old_time, new_time;
    FILE    *fp2;

    ROUTINE(plugsPhase2);

#define JAC(I,J)  ajac[(I)*npmtr+(J)]
#define JtJ(I,J)  ajtj[(I)*nvar +(J)]
#define JtQ(I)    ajtq[(I)]

    /* --------------------------------------------------------------- */

    nvar = 2 * ncloud + npmtr;
    nerr = 3 * ncloud;

    *rmsbest = 0;

    assert(ncloud>0);

    /* get necessary arrays */
    MALLOC(face,     int,      ncloud    );     // Face associated with each cloud point (bias-1)
    MALLOC(prevface, int,      ncloud    );     // Face associated on previous pass
    MALLOC(dist,     double,   ncloud    );     // min dist from each cloud to tessellation
    MALLOC(uvface,   double, 2*ncloud    );     // temp (u,v) for ocsmGetVel
    MALLOC(velface,  double, 3*ncloud    );     // temp vel   from ocsmGetVel

    MALLOC(beta,     double,   nvar      );     // optimization variables
    MALLOC(delta,    double,   nvar      );     // chg in optimization variables
    MALLOC(qerr,     double,   nerr      );     // error to be minimized
    MALLOC(qerrbest, double,   nerr      );     // error to be minimized for best so far
    MALLOC(ajac,     double,   nerr*npmtr);     // d(qerr)/d(DESPMTR)

    MALLOC(atri, double, 2*ncloud);
    MALLOC(btri, double, 2*ncloud);
    MALLOC(ctri, double, 2*ncloud);
    MALLOC(dtri, double, 2*ncloud);
    MALLOC(xtri, double, 2*ncloud);

    MALLOC(mat,  double, npmtr*npmtr);
    MALLOC(rhs,  double,       npmtr);
    MALLOC(xxx,  double,       npmtr);

    MALLOC(pmtrbest, double,  npmtr      );     // best DESPMTRs so far
    MALLOC(pmtrlast, double,  npmtr      );     // best DESPMTRS on last pass

    /* remember the best DESPMTRs */
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
        CHECK_STATUS(ocsmGetValuD);

        pmtrlast[ipmtr] = pmtrbest[ipmtr];
    }

    /* find the bounding box of the cloud (needed for classification) */
    bbox_cloud[0] = cloud[0];
    bbox_cloud[1] = cloud[1];
    bbox_cloud[2] = cloud[2];
    bbox_cloud[3] = cloud[0];
    bbox_cloud[4] = cloud[1];
    bbox_cloud[5] = cloud[2];

    for (icloud = 0; icloud < ncloud; icloud++) {
        if (cloud[3*icloud  ] < bbox_cloud[0]) bbox_cloud[0] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] < bbox_cloud[1]) bbox_cloud[1] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] < bbox_cloud[2]) bbox_cloud[2] = cloud[3*icloud+2];
        if (cloud[3*icloud  ] > bbox_cloud[3]) bbox_cloud[3] = cloud[3*icloud  ];
        if (cloud[3*icloud+1] > bbox_cloud[4]) bbox_cloud[4] = cloud[3*icloud+1];
        if (cloud[3*icloud+2] > bbox_cloud[5]) bbox_cloud[5] = cloud[3*icloud+2];
    }

    /* create .csm files to show progress at end of each pass */
    if (PLUGS_CREATE_CSM_FILES == 1) {
        char tempname[256];

        snprintf(tempname, MAX_FILENAME_LEN-1, "plugs_pass_%02d.csm", 0);
        (void) ocsmSave(MODL, tempname);
    }

    /* make multiple passes, since we can only really figure out which
       Face each cloud point is associated with when we get to a converged
       solution */
    restart = 0;
    for (ipass = 1; ipass <= npass; ipass++) {

        if (fp_plugs != NULL) {
            fprintf(fp_plugs, "\n");
        }

        /* note: odd printout is to always print pass number
           so as to void Jenkins timeouts */
        SPRINT0(1, "  ");
        SPRINT2(0, "Starting pass %d (of %d) of phase2", ipass, npass);
        SPRINT0(1, "  ");
        old_time = clock();

        /* rebuild and evaluate the new objective function */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);
        CHECK_STATUS(ocsmBuild);

        /* remember Face associations from the last time */
        for (icloud = 0; icloud < ncloud; icloud++) {
            prevface[icloud] = face[icloud];
        }

        /* only classify points that are within 0.25 of the bounding box size */
        dmax = 0.25 * MAX(MAX(bbox_cloud[3]-bbox_cloud[0],
                              bbox_cloud[4]-bbox_cloud[1]),
                              bbox_cloud[5]-bbox_cloud[2]);

        /* associate each cloud point with the closest tessellation point */
        for (icloud = 0; icloud < ncloud; icloud++) {
            face[  icloud  ] = 0;
            dist[  icloud  ] = dmax * dmax;
            beta[2*icloud  ] = 0;
            beta[2*icloud+1] = 0;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);
            CHECK_STATUS(EG_getFaceTess);

            for (icloud = 0; icloud < ncloud; icloud++) {

                /* compare with centroid of each triangle */
                for (itri = 0; itri < ntri; itri++) {
                    ip0 = tris[3*itri  ] - 1;
                    ip1 = tris[3*itri+1] - 1;
                    ip2 = tris[3*itri+2] - 1;

                    data[0] = (xyz[3*ip0  ] + xyz[3*ip1  ] + xyz[3*ip2  ]) / 3;
                    data[1] = (xyz[3*ip0+1] + xyz[3*ip1+1] + xyz[3*ip2+1]) / 3;
                    data[2] = (xyz[3*ip0+2] + xyz[3*ip1+2] + xyz[3*ip2+2]) / 3;

                    dtest = (cloud[3*icloud  ]-data[0]) * (cloud[3*icloud  ]-data[0])
                          + (cloud[3*icloud+1]-data[1]) * (cloud[3*icloud+1]-data[1])
                          + (cloud[3*icloud+2]-data[2]) * (cloud[3*icloud+2]-data[2]);
                    if (dtest < dist[icloud]) {
                        face[  icloud  ] = iface;
                        dist[  icloud  ] = dtest;
                        beta[2*icloud  ] = (uv[2*ip0  ] + uv[2*ip1  ] + uv[2*ip2  ]) / 3;
                        beta[2*icloud+1] = (uv[2*ip0+1] + uv[2*ip1+1] + uv[2*ip2+1]) / 3;
                    }
                }
            }
        }

        /* each Face must have its share of cloud points, defined as the larger of
           (arbitrarily) a tenth of the fractional area (as compared with the whole Body's area
           or (arbitrarily) 5, whichever is larger.  for Faces with fewer points,
           reclassify the cloud points that are closest to the center of the Face */
        status = EG_getMassProperties(MODL->body[ibody].ebody, massprops);
        CHECK_STATUS(EG_getMassProperties);

        atotal = massprops[1];

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            count = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] == iface) count++;
            }

            status = EG_getMassProperties(MODL->body[ibody].face[iface].eface, massprops);
            CHECK_STATUS(EG_getMassProperties);

            limit = 0.10 * ncloud * massprops[1] / atotal;
            if (limit < 5) limit = 5;

            for (ireclass = count; ireclass < limit; ireclass++) {
                status = EG_getRange(MODL->body[ibody].face[iface].eface, uvrange, &periodic);
                CHECK_STATUS(EG_getRange);

                uvrange[0] = (uvrange[0] + uvrange[1]) / 2;
                uvrange[1] = (uvrange[2] + uvrange[3]) / 2;

                status = EG_evaluate(MODL->body[ibody].face[iface].eface, uvrange, data);
                CHECK_STATUS(EG_evaluate);

                ibest = -1;
                dbest = HUGEQ;

                for (icloud = 0; icloud < ncloud; icloud++) {
                    if (face[icloud] == iface) continue;

                    dtest = (cloud[3*icloud  ]-data[0]) * (cloud[3*icloud  ]-data[0])
                          + (cloud[3*icloud+1]-data[1]) * (cloud[3*icloud+1]-data[1])
                          + (cloud[3*icloud+2]-data[2]) * (cloud[3*icloud+2]-data[2]);
                    if (dtest < dbest) {
                        ibest = icloud;
                        dbest = dtest;
                    }
                }

                face[  ibest  ] = iface;
                dist[  ibest  ] = dbest;
                beta[2*ibest  ] = uvrange[0];
                beta[2*ibest+1] = uvrange[1];

                SPRINT2(1, "WARNING:: reclassifying cloud point %5d to be associated with Face %d", ibest, iface);
            }
        }

        /* report the number of cloud points associated with each Face */
        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            count = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] == iface) count++;
            }
            SPRINT2(1, "Face %3d has %5d cloud points", iface, count);
        }

        unclass = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            if (face[icloud] <= 0) unclass++;
        }
        SPRINT1(1, "Unclassified %5d cloud points", unclass);

        for (icloud = 0; icloud < ncloud; icloud++) {
            dist[icloud] = sqrt(dist[icloud]);
        }

        /* if the faceID associated with all cloud points are the same as
           the previous pass, we are done */
        reclass = 0;
        if (ipass > 1) {
            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] != prevface[icloud]) {
                    if (reclass < 20) {
                        SPRINT6(1, "    cloud point %5d (%10.5f, %10.5f, %10.5f) was reclassified (%3d to %3d)",
                                icloud, cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2], prevface[icloud], face[icloud]);
                    }
                    reclass++;
                }
            }
            if (reclass >= 20) {
                SPRINT1(1, "    ... too many to list (%d total)", reclass);
            }
            if (restart == 0 && unclass == 0 && reclass == 0) {
                SPRINT0(0, "\n    Phase2 passes converged because points are classified same as previous pass\n");
                break;
            }
            restart = 0;
        } else {
            reclass = ncloud;
        }

        /* compute the errors and the rms */
        rms = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            iface = face[icloud];

            if (iface <= 0) {
                qerr[3*icloud  ] = 0;
                qerr[3*icloud+1] = 0;
                qerr[3*icloud+2] = 0;
                continue;
            }

            status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
            CHECK_STATUS(EG_evaluate);

            qerr[3*icloud  ] = cloud[3*icloud  ] - data[0];
            qerr[3*icloud+1] = cloud[3*icloud+1] - data[1];
            qerr[3*icloud+2] = cloud[3*icloud+2] - data[2];

            rms += qerr[3*icloud  ] * qerr[3*icloud  ]
                +  qerr[3*icloud+1] * qerr[3*icloud+1]
                +  qerr[3*icloud+2] * qerr[3*icloud+2];
        }

        rms = sqrt(rms / (3 * ncloud));

        /* print the rms and DESPMTRs */
        SPRINT2x(1, "\niter=%3d, rms=%10.3e, DESPMTRs=", -1, rms);
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
            CHECK_STATUS(ocsmGetValu);

            SPRINT1x(1, " %10.5f", pmtrbest[ipmtr]);
        }
        SPRINT0(1, " ");

        /* save initial DESPMTRs as best so far for this pass */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(beta[2*ncloud+ipmtr]), &dot);
            CHECK_STATUS(ocsmGetValu);
        }
        *rmsbest = rms;

        /* initialize for levenberg-marquardt */
        lambda    = 1;
        scaleDiag = 0;
        scaleFact = 1;
        naccept   = 0;

        /* levenberg-marquardt iterations */
        for (iter = 0; iter < niter; iter++) {
            if (scaleDiag == 0) {
                for (icloud = 0; icloud < ncloud; icloud++) {
                    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                        JAC(3*icloud  ,ipmtr) = 0;
                        JAC(3*icloud+1,ipmtr) = 0;
                        JAC(3*icloud+2,ipmtr) = 0;
                    }
                }

                /* find the sensitivities of the errors w.r.t. the DESPMTRs */
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    status = ocsmSetVelD(MODL, 0,               0, 0, 0.0);
                    CHECK_STATUS(ocsmSetVelD);

                    status = ocsmSetVelD(MODL, pmtrindx[ipmtr], 1, 1, 1.0);
                    CHECK_STATUS(ocsmSetVelD);

                    nbody = 0;
                    oldOutLevel = ocsmSetOutLevel(0);
                    status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
                    (void) ocsmSetOutLevel(oldOutLevel);
                    CHECK_STATUS(ocsmBuild);

                    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                        count = 0;
                        for (icloud = 0; icloud < ncloud; icloud++) {
                            if (face[icloud] == iface) {
                                uvface[2*count  ] = beta[2*icloud  ];
                                uvface[2*count+1] = beta[2*icloud+1];
                                count++;
                            }
                        }

                        status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, count, uvface, velface);
                        CHECK_STATUS(ocsmGetVel);

                        count = 0;
                        for (icloud = 0; icloud < ncloud; icloud++) {
                            if (face[icloud] == iface) {
                                JAC(3*icloud  ,ipmtr) = velface[3*count  ];
                                JAC(3*icloud+1,ipmtr) = velface[3*count+1];
                                JAC(3*icloud+2,ipmtr) = velface[3*count+2];
                                count++;
                            }
                        }
                    }
                }

                /* initialize matrices */
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                        mat[ipmtr*npmtr+jpmtr] = 0;
                    }
                    rhs[ipmtr] = 0;
                }

                /* fill in matricies */
                for (icloud = 0; icloud < ncloud; icloud++) {
                    iface = face[icloud];

                    /* tridiagonal matrix for (U,V)s */
                    if (iface <= 0) {
                        atri[2*icloud  ] = 0;
                        btri[2*icloud  ] = 1;
                        ctri[2*icloud  ] = 0;
                        dtri[2*icloud  ] = 0;

                        atri[2*icloud+1] = 0;
                        btri[2*icloud+1] = 1;
                        ctri[2*icloud+1] = 0;
                        dtri[2*icloud+1] = 0;
                    } else {
                        status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
                        CHECK_STATUS(EG_evaluate);

                        atri[2*icloud  ] =  0;
                        btri[2*icloud  ] = (data[3] * data[3] + data[4] * data[4] + data[5] * data[5]) * (1 + lambda);
                        ctri[2*icloud  ] =  data[3] * data[6] + data[4] * data[7] + data[5] * data[8];
                        dtri[2*icloud  ] =  data[3] * qerr[3*icloud  ]
                                         +  data[4] * qerr[3*icloud+1]
                                         +  data[5] * qerr[3*icloud+2];

                        atri[2*icloud+1] =  data[3] * data[6] + data[4] * data[7] + data[5] * data[8];
                        btri[2*icloud+1] = (data[6] * data[6] + data[7] * data[7] + data[8] * data[8]) * (1 + lambda);
                        ctri[2*icloud+1] =  0;
                        dtri[2*icloud+1] =  data[6] * qerr[3*icloud  ]
                                         +  data[7] * qerr[3*icloud+1]
                                         +  data[8] * qerr[3*icloud+2];
                    }

                    /* full matrix for DESPMTRs */
                    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                        for (jpmtr = 0; jpmtr < npmtr; jpmtr++) {
                            mat[ipmtr*npmtr+jpmtr] += JAC(3*icloud  ,ipmtr) * JAC(3*icloud  ,jpmtr)
                                                    + JAC(3*icloud+1,ipmtr) * JAC(3*icloud+1,jpmtr)
                                                    + JAC(3*icloud+2,ipmtr) * JAC(3*icloud+2,jpmtr);
                        }

                        rhs[ipmtr] += JAC(3*icloud  ,ipmtr) * qerr[3*icloud  ]
                                    + JAC(3*icloud+1,ipmtr) * qerr[3*icloud+1]
                                    + JAC(3*icloud+2,ipmtr) * qerr[3*icloud+2];
                    }
                }

                /* modify diagonal of full matrix for Levenberg-Marquardt */
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    mat[ipmtr*npmtr+ipmtr] *= (1 + lambda);
                }

            /* we have just REJECTed a change, so only the diagonal must be scaled
               because of a new lambda */
            } else {
                for (icloud = 0; icloud < ncloud; icloud++) {
                    btri[2*icloud  ] *= scaleFact;
                    btri[2*icloud+1] *= scaleFact;
                }
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    mat[ipmtr*npmtr+ipmtr] *= scaleFact;
                }
            }

            /* solve tridiagonal system to update the (U,V)s */
            status = tridiag(2*ncloud, atri, btri, ctri, dtri, xtri);
            CHECK_STATUS(tridiag);

            /* solve the full system to update the DESPMTRs */
            status = matsol(mat, rhs, npmtr, xxx);
            CHECK_STATUS(matsol);

            /* assemble the soltions into the delta array */
            for (icloud = 0; icloud < ncloud; icloud++) {
                delta[2*icloud  ] = xtri[2*icloud  ];
                delta[2*icloud+1] = xtri[2*icloud+1];
            }
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                delta[2*ncloud+ipmtr] = xxx[ipmtr];
            }

            /* check for no change in DESPMTRs (since it is unlikely that rms will
               ever go to zero) */
            dmax = 0;
            for (ivar = 0; ivar < 2*ncloud+npmtr; ivar++) {
                if (fabs(delta[ivar]) > dmax) {
                    dmax = fabs(delta[ivar]);
                }
            }
            if (dmax < EPS06 && lambda <= 1) {
                SPRINT5(0, "    Pass %2d converged, dmax=%10.3e,     rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                        ipass, dmax, *rmsbest, reclass, unclass);
                break;
            }

            /* temporarily apply change */
            for (icloud = 0; icloud < ncloud; icloud++) {
                if (face[icloud] > 0) {
                    beta[2*icloud  ] += delta[2*icloud  ];
                    beta[2*icloud+1] += delta[2*icloud+1];
                }
            }
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                value = pmtrbest[ipmtr] + delta[2*ncloud+ipmtr];

                status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
                CHECK_STATUS(ocsmGetBnds);

                if (value < lbound) value = lbound;
                if (value > ubound) value = ubound;
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
                CHECK_STATUS(ocsmSetValuD);
            }

            /* rebuild and evaluate the new objective function */
            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);

            /* compute the errors and the rms */
            if (status != SUCCESS) {
                rms = 1 + *rmsbest;
            } else {
                rms = 0;
                for (icloud = 0; icloud < ncloud; icloud++) {
                    qerrbest[3*icloud  ] = qerr[3*icloud  ];
                    qerrbest[3*icloud+1] = qerr[3*icloud+1];
                    qerrbest[3*icloud+2] = qerr[3*icloud+2];

                    iface = face[icloud];

                    if (iface <= 0) {
                        qerr[3*icloud  ] = 0;
                        qerr[3*icloud+1] = 0;
                        qerr[3*icloud+2] = 0;
                        continue;
                    }

                    status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
                    CHECK_STATUS(EG_evaluate);

                    qerr[3*icloud  ] = cloud[3*icloud  ] - data[0];
                    qerr[3*icloud+1] = cloud[3*icloud+1] - data[1];
                    qerr[3*icloud+2] = cloud[3*icloud+2] - data[2];

                    rms += qerr[3*icloud  ] * qerr[3*icloud  ]
                        +  qerr[3*icloud+1] * qerr[3*icloud+1]
                        +  qerr[3*icloud+2] * qerr[3*icloud+2];
                }

                rms = sqrt(rms / (3 * ncloud));
            }

            /* print the rms and DESPMTRs */
            SPRINT2x(1, "iter=%3d, rms=%10.3e, DESPMTRs=", iter, rms);
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
                CHECK_STATUS(ocsmGetValu);

                SPRINT1x(1, " %10.5f", value);
            }

            /* store history */
            if (fp_plugs != NULL) {
                fprintf(fp_plugs, " %4d %4d", ipass, iter);
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &value, &dot);
                    CHECK_STATUS(ocsmGetPmtr);

                    fprintf(fp_plugs, " %11.6f", value);
                }
                if (rms < *rmsbest) {
                    fprintf(fp_plugs, "%11.5f ACCEPT\n", rms);
                } else {
                    fprintf(fp_plugs, "%11.5f\n", rms);
                }
                fflush( fp_plugs);
            }

            if (PLUGS_CREATE_EGADS_FILES > 0) {
                status = EG_makeTopology(MODL->context, NULL, MODEL, 0, NULL,
                                         1, &(MODL->body[MODL->nbody].ebody), NULL, &emodel);
                CHECK_STATUS(EG_makeTopology);

                snprintf(egadsname, 80, "plugs_%d_%d.egads", ipass, iter);

                (void) remove(egadsname);

                status = EG_saveModel(emodel, egadsname);
                CHECK_STATUS(EG_saveModel);

                status = EG_deleteObject(emodel);
                CHECK_STATUS(EG_deleteObject);
            }

            /* if this new solution is better, accept it and decrease lambda
               (making it more newton-like) */
            if (rms < *rmsbest) {
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    status = ocsmGetValu(MODL, pmtrindx[ipmtr], 1, 1, &(pmtrbest[ipmtr]), &dot);
                    CHECK_STATUS(ocsmGetPmtr);
                }
                *rmsbest  = rms;
                scaleDiag = 0;
                naccept++;
                lambda    = MAX(1.0e-10, lambda/2);
                SPRINT1(1, "  accepted: lambda=%10.3e", lambda);

            /* otherwise reject the step and increase lambda
               (making it more steepest-descent-like) */
            } else {
                for (icloud = 0; icloud < ncloud; icloud++) {
                    if (face[icloud] > 0) {
                        beta[2*icloud  ] -= delta[2*icloud  ];
                        beta[2*icloud+1] -= delta[2*icloud+1];
                    }

                    qerr[3*icloud  ] = qerrbest[3*icloud  ];
                    qerr[3*icloud+1] = qerrbest[3*icloud+1];
                    qerr[3*icloud+2] = qerrbest[3*icloud+2];
                }
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                    CHECK_STATUS(ocsmSetValuD);
                }
                scaleDiag  = 1;
                scaleFact  = 1 / (1 + lambda);
                lambda     = MIN(1.0e+10, lambda*2);
                scaleFact *= (1 + lambda);
                SPRINT1(1, "  rejected: lambda=%10.3e", lambda);
            }

            /* if lambda gets very big, stop iterations */
            if (lambda > 100) {
                SPRINT5(0, "    Pass %2d has stalled, lambda=%10.3e, rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                        ipass, lambda, *rmsbest, reclass, unclass);
                break;
            }
        }

        if (iter >= niter-1) {
            SPRINT4(0, "    Pass %2d ran out of iterations,          rmsbest=%10.3e, reclass=%5d, unclass=%5d",
                    ipass, *rmsbest, reclass, unclass);
        }

        /* if the last build was a rejection, rebuild with the best */
        if (rms >= *rmsbest) {
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                CHECK_STATUS(ocsmSetValuD);
            }

            status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
            CHECK_STATUS(ocsmSetVelD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);
            if (status < SUCCESS) {
                for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                    SPRINT3(0, "error  DESPMTR %3d: %20s = %10.5f",
                            ipmtr, MODL->pmtr[pmtrindx[ipmtr]].name, MODL->pmtr[pmtrindx[ipmtr]].value[0]);
                }
            }
            CHECK_STATUS(ocsmBuild);
        }

        new_time = clock();
        SPRINT2(1, "Phase 2 (pass %2d) CPUtime=%9.3f sec", ipass,
                (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

        /* create .csm files to show progress at end of each pass */
        if (PLUGS_CREATE_CSM_FILES == 1) {
            char tempname[256];

            snprintf(tempname, MAX_FILENAME_LEN-1, "plugs_pass_%02d.csm", ipass);
            (void) ocsmSave(MODL, tempname);
        }

        /* create final.plot to show how close cloud is to final surface at
           the end of this pass */
        if (PLUGS_CREATE_FINAL_PLOT == 1 || naccept == 0) {
            fp2 = fopen("final.plot", "w");

            if (fp2 != NULL) {
                count = 0;
                for (icloud = 0; icloud < ncloud; icloud++) {
                    if (face[icloud] == 0) count++;
                }

                if (count > 0) {
                    fprintf(fp2, "%5d%5d Unclassified_cloud_points\n", count, 0);

                    for (icloud = 0; icloud < ncloud; icloud++) {
                        if (face[icloud] == 0) {
                            fprintf(fp2, " %9.5f %9.5f %9.5f\n", cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2]);
                        }
                    }
                }

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    count = 0;
                    for (icloud = 0; icloud < ncloud; icloud++) {
                        if (face[icloud] == iface) count++;
                    }

                    if (count > 0) {
                        fprintf(fp2, "%5d%5d Face_%d_cloud_points\n", count, 0, iface);

                        for (icloud = 0; icloud < ncloud; icloud++) {
                            if (face[icloud] == iface) {
                                fprintf(fp2, " %9.5f %9.5f %9.5f\n", cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2]);
                            }
                        }
                    }
                }

                for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                    count = 0;
                    for (icloud = 0; icloud < ncloud; icloud++) {
                        if (face[icloud] == iface) count++;
                    }

                    if (count > 0) {
                        fprintf(fp2, "%5d%5d Face_%d_distances\n", count, -1, iface);

                        for (icloud = 0; icloud < ncloud; icloud++) {
                            if (face[icloud] == iface) {
                                status = EG_evaluate(MODL->body[ibody].face[iface].eface, &(beta[2*icloud]), data);
                                CHECK_STATUS(EG_evaluate);

                                fprintf(fp2, " %9.5f %9.5f %9.5f\n", cloud[3*icloud], cloud[3*icloud+1], cloud[3*icloud+2]);
                                fprintf(fp2, " %9.5f %9.5f %9.5f\n", data[0],         data[1],           data[2]          );
                            }
                        }
                    }
                }

                fprintf(fp2, "%5d%5d end\n", 0, 0);
                fclose(fp2);
            }
        }

        /* if none of the LM itertions were accepted, check if perturbing a DESPMTR will allow us to restart */
        if (naccept == 0) {
            SPRINT0(0, "\n    ERROR:: no LM iterations were accepted, so checking if any perturbations will be better\n");

            /* get the baseline */
            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            CHECK_STATUS(ocsmBuild);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms     = 0;
            unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    unclass++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-unclass)));

            SPRINT1(0, "    baseline                     rms=%12.5e\n", rms);
            iperturb   = 0;
            rmsperturb = rms;
            valperturb = HUGEQ;

            /* now pertrub each of the DESPMTRs, one at a time, and compute the rms by inverse evaluations */
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {

                status = ocsmGetBnds(MODL, pmtrindx[ipmtr], 1, 1, &lbound, &ubound);
                CHECK_STATUS(ocsmGetBnds);

                /* decrease by 5% */
                if (pmtrbest[ipmtr] > 0) {
                    value = MAX(pmtrbest[ipmtr]/1.05, lbound);
                } else {
                    value = MIN(pmtrbest[ipmtr]/1.05, ubound);
                }
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
                CHECK_STATUS(ocsmSetValuD);

                nbody = 0;
                oldOutLevel = ocsmSetOutLevel(0);
                status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
                CHECK_STATUS(ocsmBuild);
                (void) ocsmSetOutLevel(oldOutLevel);

                rms     = 0;
                unclass = 0;
                for (icloud = 0; icloud < ncloud; icloud++) {
                    iface       = face[  icloud  ];
                    uv_guess[0] = beta[2*icloud  ];
                    uv_guess[1] = beta[2*icloud+1];

                    if (iface > 0) {
                        status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                        CHECK_STATUS(EG_invEvaluateGuess);

                        rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                             + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                             + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                    } else {
                        unclass++;
                    }
                }
                rms = sqrt(rms / (3*(ncloud-unclass)));

                SPRINT4(1, "    ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/best=%10.5f", ipmtr, value, rms, rms/rmsperturb);

                if (rms < rmsperturb) {
                    iperturb   = -ipmtr;
                    rmsperturb = rms;
                    valperturb = value;
                }

                /* increase by 5% */
                if (pmtrbest[ipmtr] > 0) {
                    value = MIN(pmtrbest[ipmtr]*1.05, ubound);
                } else {
                    value = MAX(pmtrbest[ipmtr]*1.05, lbound);
                }
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
                CHECK_STATUS(ocsmSetValuD);

                nbody = 0;
                oldOutLevel = ocsmSetOutLevel(0);
                status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
                CHECK_STATUS(ocsmBuild);
                (void) ocsmSetOutLevel(oldOutLevel);

                rms     = 0;
                unclass = 0;
                for (icloud = 0; icloud < ncloud; icloud++) {
                    iface       = face[  icloud  ];
                    uv_guess[0] = beta[2*icloud  ];
                    uv_guess[1] = beta[2*icloud+1];

                    if (iface > 0) {
                        status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                        CHECK_STATUS(EG_invEvaluateGuess);

                        rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                             + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                             + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                    } else {
                        unclass++;
                    }
                }
                rms = sqrt(rms / (3*(ncloud-unclass)));

                SPRINT4(1, "    ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/best=%10.5f\n", ipmtr, value, rms, rms/rmsperturb);

                if (rms < rmsperturb) {
                    iperturb   = +ipmtr;
                    rmsperturb = rms;
                    valperturb = value;
                }

                /* set back to nominal value */
                status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
                CHECK_STATUS(ocsmSetValuD);
            }

            if (iperturb > 0) {
                pmtrbest[+iperturb] = valperturb;

                status = ocsmSetValuD(MODL, pmtrindx[+iperturb], 1, 1, pmtrbest[+iperturb]);
                CHECK_STATUS(ocsmSetValuD);

                SPRINT2(1, "    restarting with ipmtr=%2d perturbed to %12.5f", +iperturb, pmtrbest[+iperturb]);

                restart = 1;
            } else if (iperturb < 0) {
                pmtrbest[-iperturb] = valperturb;

                SPRINT2(1, "    restarting with ipmtr=%2d perturbed to %12.5f", -iperturb, pmtrbest[-iperturb]);

                status = ocsmSetValuD(MODL, pmtrindx[-iperturb], 1, 1, pmtrbest[-iperturb]);
                CHECK_STATUS(ocsmSetValuD);

                restart = 1;
            } else {
                SPRINT0(1, "    no perturbation succeeded, so stopping");

                break;
            }

        /* if the DESPMTRs are essentially the same on last two passes and all
           cloud points are classified, we are done */
        } else if (ipass > 0) {
            dmax = 0;
            for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
                if (fabs(pmtrbest[ipmtr]-pmtrlast[ipmtr]) > dmax) {
                    dmax = fabs(pmtrbest[ipmtr]-pmtrlast[ipmtr]);
                }
            }

            if (unclass == 0 && dmax < 1.0e-5) {
                SPRINT3(0, "\n    Phase2 passes converged because maximum DESPMTR change is %10.3e, reclass=%5d, unclass=%5d\n",
                        dmax, reclass, unclass);
                break;
            }
        }

        /* remember DESPMTRs for next pass */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
            pmtrlast[ipmtr] = pmtrbest[ipmtr];
        }
    }

    /* set the DESPMTRs to the final converged values */
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
        CHECK_STATUS(ocsmSetValuD);
    }

    /* find the sensitivity of the objective function to perturbations in
       the various DESPMTRs */
    if (PLUGS_SHOW_SENSITIVITY == 1) {
        double rmsuv, baseline;

        rms = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            rms += qerr[3*icloud]*qerr[3*icloud] + qerr[3*icloud+1]*qerr[3*icloud+1] + qerr[3*icloud+2]*qerr[3*icloud+2];
        }
        rms = sqrt(rms / (3*ncloud));

        /* baseline rms by inverse evaluations */
        nbody = 0;
        oldOutLevel = ocsmSetOutLevel(0);
        status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
        (void) ocsmSetOutLevel(oldOutLevel);

        rms     = 0;
        rmsuv   = 0;
        unclass = 0;
        for (icloud = 0; icloud < ncloud; icloud++) {
            iface       = face[  icloud  ];
            uv_guess[0] = beta[2*icloud  ];
            uv_guess[1] = beta[2*icloud+1];

            if (iface > 0) {
                status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                CHECK_STATUS(EG_invEvaluateGuess);

                rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                     + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                     + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                rmsuv += (uv_guess[0]-beta[2*icloud  ]) * (uv_guess[0]-beta[2*icloud  ])
                       + (uv_guess[1]-beta[2*icloud+1]) * (uv_guess[1]-beta[2*icloud+1]);
            } else {
                unclass++;
            }
        }
        rms   = sqrt(rms   / (3*(ncloud-unclass)));
        rmsuv = sqrt(rmsuv / (2*(ncloud-unclass)));

        SPRINT2(1, "baseline                     rms=%12.5e, rmsuv=%12.5e\n", rms, rmsuv);
        baseline = rms;

        /* now pertrub each of the DESPMTRs, one at a time, and compute the rms by inverse evaluations */
        for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {

            /* decrease by 5% */
            if (pmtrbest[ipmtr] > 0) {
                value = MAX(pmtrbest[ipmtr]/1.05, lbound);
            } else {
                value = MIN(pmtrbest[ipmtr]/1.05, ubound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms     = 0;
            unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    unclass++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-unclass)));

            SPRINT4(1, "ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/baseline=%10.5f", ipmtr, value, rms, rms/baseline);

            /* decrease by 1% */
            if (pmtrbest[ipmtr] > 0) {
                value = MAX(pmtrbest[ipmtr]/1.01, lbound);
            } else {
                value = MIN(pmtrbest[ipmtr]/1.01, ubound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms     = 0;
            unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    unclass++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-unclass)));

            SPRINT4(1, "ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/baseline=%10.5f", ipmtr, value, rms, rms/baseline);

            SPRINT4(1, "ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/baseline=%10.5f", ipmtr, pmtrbest[ipmtr], baseline, 1.0);

            /* increase by 1% */
            if (pmtrbest[ipmtr] > 0) {
                value = MIN(pmtrbest[ipmtr]*1.01, ubound);
            } else {
                value = MAX(pmtrbest[ipmtr]*1.01, lbound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms     = 0;
            unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                         + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                         + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    unclass++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-unclass)));

            SPRINT4(1, "ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/baseline=%10.5f", ipmtr, value, rms, rms/baseline);

            /* increase by 5% */
            if (pmtrbest[ipmtr] > 0) {
                value = MIN(pmtrbest[ipmtr]*1.05, ubound);
            } else {
                value = MAX(pmtrbest[ipmtr]*1.05, lbound);
            }
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            nbody = 0;
            oldOutLevel = ocsmSetOutLevel(0);
            status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
            (void) ocsmSetOutLevel(oldOutLevel);

            rms     = 0;
            unclass = 0;
            for (icloud = 0; icloud < ncloud; icloud++) {
                iface       = face[  icloud  ];
                uv_guess[0] = beta[2*icloud  ];
                uv_guess[1] = beta[2*icloud+1];

                if (iface > 0) {
                    status = EG_invEvaluateGuess(MODL->body[MODL->nbody].face[iface].eface, &(cloud[3*icloud]), uv_guess, xyz_guess);
                    CHECK_STATUS(EG_invEvaluateGuess);

                    rms += (cloud[3*icloud  ]-xyz_guess[0]) * (cloud[3*icloud  ]-xyz_guess[0])
                        + (cloud[3*icloud+1]-xyz_guess[1]) * (cloud[3*icloud+1]-xyz_guess[1])
                        + (cloud[3*icloud+2]-xyz_guess[2]) * (cloud[3*icloud+2]-xyz_guess[2]);
                } else {
                    unclass++;
                }
            }
            rms = sqrt(rms / (3*(ncloud-unclass)));

            SPRINT4(1, "ipmtr=%2d, valu=%12.7f, rms=%12.5e, rms/baseline=%10.5f\n", ipmtr, value, rms, rms/baseline);

            /* set back to nominal value */
            status = ocsmSetValuD(MODL, pmtrindx[ipmtr], 1, 1, pmtrbest[ipmtr]);
            CHECK_STATUS(ocsmSetValuD);
        }
    }

cleanup:
    FREE(face    );
    FREE(prevface);
    FREE(dist    );
    FREE(uvface  );
    FREE(velface );

    FREE(beta    );
    FREE(delta   );
    FREE(qerr    );
    FREE(qerrbest);
    FREE(ajac    );

    FREE(atri);
    FREE(btri);
    FREE(ctri);
    FREE(dtri);
    FREE(xtri);

    FREE(mat);
    FREE(rhs);
    FREE(xxx);

    FREE(pmtrbest);
    FREE(pmtrlast);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   matsol - Gaussian elimination with partial pivoting                *
 *                                                                      *
 ************************************************************************
 */

static int
matsol(double    A[],                   /* (in)  matrix to be solved (stored rowwise) */
                                        /* (out) upper-triangular form of matrix */
       double    b[],                   /* (in)  right hand side */
                                        /* (out) right-hand side after swapping */
       int       n,                     /* (in)  size of matrix */
       double    x[])                   /* (out) solution of A*x=b */
{
    int       status = SUCCESS;         /* (out) return status */

    int       ir, jc, kc, imax;
    double    amax, swap, fact;

    ROUTINE(matsol);

    /* --------------------------------------------------------------- */

    /* reduce each column of A */
    for (kc = 0; kc < n; kc++) {

        /* find pivot element */
        imax = kc;
        amax = fabs(A[kc*n+kc]);

        for (ir = kc+1; ir < n; ir++) {
            if (fabs(A[ir*n+kc]) > amax) {
                imax = ir;
                amax = fabs(A[ir*n+kc]);
            }
        }

        /* check for possibly-singular matrix (ie, near-zero pivot) */
        if (amax < EPS12) {
            status = OCSM_SINGULAR_MATRIX;
            goto cleanup;
        }

        /* if diagonal is not pivot, swap rows in A and b */
        if (imax != kc) {
            for (jc = 0; jc < n; jc++) {
                swap         = A[kc  *n+jc];
                A[kc  *n+jc] = A[imax*n+jc];
                A[imax*n+jc] = swap;
            }

            swap    = b[kc  ];
            b[kc  ] = b[imax];
            b[imax] = swap;
        }

        /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
        for (ir = kc+1; ir < n; ir++) {
            fact = A[ir*n+kc] / A[kc*n+kc];

            for (jc = kc+1; jc < n; jc++) {
                A[ir*n+jc] -= fact * A[kc*n+jc];
            }

            b[ir] -= fact * b[kc];

            A[ir*n+kc] = 0;
        }
    }

    /* back-substitution pass */
    x[n-1] = b[n-1] / A[(n-1)*n+(n-1)];

    for (jc = n-2; jc >= 0; jc--) {
        x[jc] = b[jc];
        for (kc = jc+1; kc < n; kc++) {
            x[jc] -= A[jc*n+kc] * x[kc];
        }
        x[jc] /= A[jc*n+jc];
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   solsvd - solve  A * x = b  via singular value decomposition        *
 *                                                                      *
 ************************************************************************
 */

static int
solsvd(double A[],                      /* (in)  mrow*ncol matrix */
       double b[],                      /* (in)  mrow-vector (right-hand sides) */
       int    mrow,                     /* (in)  number of rows */
       int    ncol,                     /* (in)  number of columns */
       double W[],                      /* (out) ncol-vector of singular values */
       double x[])                      /* (out) ncol-vector (solution) */
{
    int    status = SUCCESS;            /* (out) return status */

    int    irow, jcol, i, j, k, flag, its, jj, ip1=0, nm=0;
    double *U=NULL, *V=NULL, *r=NULL, *t=NULL;
    double wmin, wmax, s, anorm, c, f, g, h, scale, xx, yy, zz;

    ROUTINE(solsvd);

    /* --------------------------------------------------------------- */

    /* this routine is an adaptation of svf.f found in the netlib
       repository.  it is a modification of a routine from the eispack
       collection, which in turn is a translation of the algol procedure svd,
       num. math. 14, 403-420(1970) by golub and reinsch.
       handbook for auto. comp., vol ii-linear algebra, 134-151(1971). */

    /* default return */
    for (jcol = 0; jcol < ncol; jcol++) {
        x[jcol] = 0;
    }

    /* check for legal size for A */
    if (ncol <= 0) {
        status = OCSM_ILLEGAL_VALUE;
        goto cleanup;
    } else if (mrow < ncol) {
        status = OCSM_ILLEGAL_VALUE;
        goto cleanup;
    }

    MALLOC(U, double, mrow*ncol);
    MALLOC(V, double, ncol*ncol);
    MALLOC(r, double, ncol     );
    MALLOC(t, double, ncol     );

    /* initializations needed to avoid clang warning */
    for (i = 0; i < ncol; i++) {
        W[i] = 0;
        r[i] = 0;
        t[i] = 0;
    }

    /* initialize U to the original A */
    for (irow = 0; irow < mrow; irow++) {
        for (jcol = 0; jcol < ncol; jcol++) {
            U[irow*ncol+jcol] = A[irow*ncol+jcol];
        }
    }

    /* decompose A into U*W*V' */
    g     = 0;
    scale = 0;
    anorm = 0;

    /* Householder reduction of U to bidiagonal form */
    for (i = 0; i < ncol; i++) {
        ip1   = i + 1;
        r[i]  = scale * g;
        g     = 0;
        s     = 0;
        scale = 0;
        if (i < mrow) {
            for (k = i; k < mrow; k++) {
                scale += fabs(U[k*ncol+i]);
            }
            if (scale != 0) {
                for (k = i; k < mrow; k++) {
                    U[k*ncol+i] /= scale;
                    s           += U[k*ncol+i] * U[k*ncol+i];
                }
                f           = U[i*ncol+i];
                g           = -FSIGN(sqrt(s), f);
                h           = f * g - s;
                U[i*ncol+i] = f - g;
                for (j = ip1; j < ncol; j++) {
                    s = 0;
                    for (k = i; k < mrow; k++) {
                        s += U[k*ncol+i] * U[k*ncol+j];
                    }
                    f = s / h;
                    for (k = i; k < mrow; k++) {
                        U[k*ncol+j] += f * U[k*ncol+i];
                    }
                }
                for (k = i; k < mrow; k++) {
                    U[k*ncol+i] *= scale;
                }
            }
        }
        W[i]  = scale  * g;
        g     = 0;
        s     = 0;
        scale = 0;
        if (i < mrow && i+1 != ncol) {
            for (k = ip1; k < ncol; k++) {
                scale += fabs(U[i*ncol+k]);
            }
            if (scale != 0) {
                for (k = ip1; k < ncol; k++) {
                    U[i*ncol+k] /= scale;
                    s           += U[i*ncol+k] * U[i*ncol+k];
                }
                f           = U[i*ncol+ip1];
                g           = -FSIGN(sqrt(s), f);
                h           = f * g - s;
                U[i*ncol+ip1] = f - g;
                for (k = ip1; k < ncol; k++) {
                    r[k] = U[i*ncol+k] / h;
                }
                for (j = ip1; j < mrow; j++) {
                    s = 0;
                    for (k = ip1; k < ncol; k++) {
                        s += U[j*ncol+k] * U[i*ncol+k];
                    }
                    for (k = ip1; k < ncol; k++) {
                        U[j*ncol+k] += s * r[k];
                    }
                }
                for (k = ip1; k < ncol; k++) {
                    U[i*ncol+k] *= scale;
                }
            }
        }
        anorm = MAX(anorm, (fabs(W[i]) + fabs(r[i])));
    }

    /* accumulation of right-hand transformations */
    for (i = ncol-1; i >= 0; i--) {
        if (i < ncol-1) {
            if (g != 0) {
                for (j = ip1; j < ncol; j++) {
                    V[j*ncol+i] = (U[i*ncol+j] / U[i*ncol+ip1]) / g; /* avoid possible underflow */
                }
                for (j = ip1; j < ncol; j++) {
                    s = 0;
                    for (k = ip1; k < ncol; k++) {
                        s += U[i*ncol+k] * V[k*ncol+j];
                    }
                    for (k = ip1; k < ncol; k++) {
                        V[k*ncol+j] += s * V[k*ncol+i];
                    }
                }
            }
            for (j = ip1; j < ncol; j++) {
                V[i*ncol+j] = 0;
                V[j*ncol+i] = 0;
            }
        }
        V[i*ncol+i] = 1;
        g           = r[i];
        ip1         = i;
    }

    /* accumulation of left-side transformations */
    for (i = MIN(mrow, ncol)-1; i >= 0; i--) {
        ip1 = i + 1;
        g = W[i];
        for (j = ip1; j < ncol; j++) {
            U[i*ncol+j] = 0;
        }
        if (g != 0) {
            g = 1 / g;
            for (j = ip1; j < ncol; j++) {
                s = 0;
                for (k = ip1; k < mrow; k++) {
                    s += U[k*ncol+i] * U[k*ncol+j];
                }
                f = (s / U[i*ncol+i]) * g;
                for (k = i; k < mrow; k++) {
                    U[k*ncol+j] += f * U[k*ncol+i];
                }
            }
            for (j = i; j < mrow; j++) {
                U[j*ncol+i] *= g;
            }
        } else {
            for (j = i; j < mrow; j++) {
                U[j*ncol+i] = 0;
            }
        }
        ++U[i*ncol+i];
    }

    /* diagonalization of the bidiagonal form */

    /* loop over singular values */
    for (k = ncol-1; k >= 0; k--) {

        /* loop over allowed iterations */
        for (its = 0; its < 30; its++) {

            /* test for splitting */
            flag = 1;
            for (ip1 = k; ip1 >= 0; ip1--) {
                nm = ip1 - 1;

                if ((double)(fabs(r[ip1]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }

                assert (nm >= 0);                 /* needed to avoid clang warning */
                assert (nm < ncol);               /* needed to avoid clang warning */

                if ((double)(fabs(W[nm]) + anorm) == anorm) break;
            }
            if (flag) {
                c = 0;
                s = 1;
                for (i = ip1; i < k+1; i++) {
                    f    = s * r[i];
                    r[i] = c * r[i];
                    if ((double)(fabs(f) + anorm) == anorm) break;
                    g    = W[i];
                    if (fabs(f) > fabs(g)) {
                        h = fabs(f) * sqrt(1 + (g/f) * (g/f));
                    } else if (fabs(g) == 0) {
                        h = 0;
                    } else {
                        h = fabs(g) * sqrt(1 + (f/g) * (f/g));
                    }
                    W[i] = h;
                    h    = 1 / h;
                    c    = g * h;
                    s    = -f * h;
                    for (j = 0; j < mrow; j++) {
                        yy           = U[j*ncol+nm];
                        zz           = U[j*ncol+i ];
                        U[j*ncol+nm] = yy * c + zz * s;
                        U[j*ncol+i ] = zz * c - yy * s;
                    }
                }
            }

            /* test for convergence */
            zz = W[k];
            if (ip1 == k) {

                /* make singular values non-negative */
                if (zz < 0) {
                    W[k] = -zz;
                    for (j = 0; j < ncol; j++) {
                        V[j*ncol+k] = -V[j*ncol+k];
                    }
                }
                break;
            }

            assert (ip1 >= 0);                    /* needed to avoid clang warning */
            assert (ip1 < ncol);                  /* needed to avoid clang warning */

            /* shift from bottom 2*2 minor */
            xx = W[ip1];
            nm = k - 1;
            yy = W[nm];
            g  = r[nm];
            h  = r[k];
            f  = ((yy - zz) * (yy + zz) + (g - h) * (g + h)) / (2 * h * yy);
            g  = sqrt(f * f + 1);
            f  = ((xx - zz) * (xx + zz) + h * ((yy / (f + FSIGN(g, f))) - h)) / xx;

            /* next QR transformation */
            c = 1;
            s = 1;
            for (j = ip1; j <= nm; j++) {
                i    = j + 1;
                g    = r[i];
                yy   = W[i];
                h    = s * g;
                g    = c * g;
                if (fabs(f) > fabs(h)) {
                    zz = fabs(f) * sqrt(1 + (h/f) * (h/f));
                } else if (fabs(h) == 0) {
                    zz = 0;
                } else {
                    zz = fabs(h) * sqrt(1 + (f/h) * (f/h));
                }
                r[j] = zz;
                c    = f / zz;
                s    = h / zz;
                f    = xx * c + g * s;
                g    = g * c - xx * s;
                h    = yy * s;
                yy  *= c;
                for (jj = 0; jj < ncol; jj++) {
                    xx           = V[jj*ncol+j];
                    zz           = V[jj*ncol+i];
                    V[jj*ncol+j] = xx * c + zz * s;
                    V[jj*ncol+i] = zz * c - xx * s;
                }
                if (fabs(f) > fabs(h)) {
                    zz = fabs(f) * sqrt(1 + (h/f) * (h/f));
                } else if (fabs(h) == 0) {
                    zz = 0;
                } else {
                    zz = fabs(h) * sqrt(1 + (f/h) * (f/h));
                }

                /* rotation can be arbitrary if zz=0 */
                W[j] = zz;
                if (zz != 0) {
                    zz = 1 / zz;
                    c  = f * zz;
                    s  = h * zz;
                }
                f  = c * g  + s * yy;
                xx = c * yy - s * g;
                for (jj = 0; jj < mrow; jj++) {
                    yy           = U[jj*ncol+j];
                    zz           = U[jj*ncol+i];
                    U[jj*ncol+j] = yy * c + zz * s;
                    U[jj*ncol+i] = zz * c - yy * s;
                }
            }
            r[ip1] = 0;
            r[k  ] = f;
            W[k  ] = xx;
        }
    }

    /* find the largest singular value (for scaling) */
    wmax = 0;
    for (jcol = 0; jcol < ncol; jcol++) {
        if (W[jcol] > wmax) {
            wmax = W[jcol];
        }
    }

    /* set all singular values less than wmin to zero */
    wmin = wmax * 1.0e-6;
    for (jcol = 0; jcol < ncol; jcol++) {
        if (W[jcol] < wmin) {
            W[jcol] = 0;
        }
    }

    /* perform the back-substitution */
    for (j = 0; j < ncol; j++) {
        s = 0;
        if (W[j] != 0) {
            for (i = 0; i < mrow; i++) {
                s += U[i*ncol+j] * b[i];
            }
            s /= W[j];
        }
        t[j] = s;
    }

    for (j = 0; j < ncol; j++) {
        s = 0;
        for (k = 0; k < ncol; k++) {
            s += V[j*ncol+k] * t[k];
        }
        x[j] = s;
    }

cleanup:
    FREE(t);
    FREE(r);
    FREE(V);
    FREE(U);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   tridiag - solve tridiaginal system                                 *
 *                                                                      *
 ************************************************************************
 */

static int
tridiag(int    n,                       /* (in)  size of system */
        double a[],                     /* (in)  sub-diagonals */
        double b[],                     /* (in)  diagonals */
        double c[],                     /* (in)  super-diagonals */
        double d[],                     /* (in)  right-hand side */
        double x[])                     /* (out) solution */
{
    int    status = SUCCESS;            /* (out) return status */

    int    i;
    double W, *p=NULL, *q=NULL;

    ROUTINE(tridiag);

    /* --------------------------------------------------------------- */

    MALLOC(p, double, n);
    MALLOC(q, double, n);

    /* forward elimination */
    p[0] = -c[0] / b[0];
    q[0] =  d[0] / b[0];

    for (i = 1; i < n; i++) {
        W     = b[i] + c[i] * p[i-1];
        p[i] = -a[i] / W;
        q[i] = (d[i] - c[i] * q[i-1]) / W;
    }

    /* final solution */
    x[n-1] = q[n-1];

    /* back substitution */
    for (i = n-2; i >= 0; i--) {
        x[i] = p[i] * x[i+1] + q[i];
    }

cleanup:
    FREE(p);
    FREE(q);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   writeSensFile - write ASCII .sens file                             *
 *                                                                      *
 ************************************************************************
 */

/*
 * NOTE: this is a duplicate of a routine in OpenCSM.c
 */

static int
writeSensFile(modl_T *MODL,             /* (in)  pointer to MODL */
              int    ibody,             /* (in)  Body index (1:nbody) */
              char   filename[])        /* (in)  filename */
{
    int    status = SUCCESS;            /* (out) return status */

    int     count, ipmtr, i, inode, iedge, iface, npnt, ipnt, ntri, itri;
    CINT    *pindx, *ptype, *tris, *tric;
    double  *vels=NULL;
    CDOUBLE *xyz, *uv, *Vels;
    FILE    *fp;

    ROUTINE(writeSensFile);

    /* --------------------------------------------------------------- */

    /* count the number od DESPMTRs */
    count = 0;
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
            count++;
        }
    }

    /* open the file */
    fp = fopen(filename, "w");
    if (fp == NULL) {
        status = OCSM_FILE_NOT_FOUND;
        goto cleanup;
    }

    /* write the DESPMTRs in the header of the file */
    fprintf(fp, "%8d\n", count);

    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
            fprintf(fp, "%8d %s\n",
                    MODL->pmtr[ipmtr].nrow*MODL->pmtr[ipmtr].ncol,
                    MODL->pmtr[ipmtr].name);

            for (i = 0; i < MODL->pmtr[ipmtr].nrow*MODL->pmtr[ipmtr].ncol; i++) {
                fprintf(fp, "     %22.15e %22.15e\n",
                        MODL->pmtr[ipmtr].value[i], MODL->pmtr[ipmtr].dot[i]);
            }
        }
    }

    fprintf(fp, "%8d %8d %8d\n",
            MODL->body[ibody].nnode,
            MODL->body[ibody].nedge,
            MODL->body[ibody].nface);

    /* write Nodes to the files */
    for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
        MALLOC(vels, double, 3);

        status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, vels);
        CHECK_STATUS(ocsmGetVel);

        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
                MODL->body[ibody].node[inode].x,
                MODL->body[ibody].node[inode].y,
                MODL->body[ibody].node[inode].z,
                vels[0], vels[1], vels[2]);

        FREE(vels);
    }

    /* write Edges to the file */
    for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
        status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        status = ocsmGetTessVel(MODL, ibody, OCSM_EDGE, iedge, &Vels);
        CHECK_STATUS(ocsmGetTessVel);

        fprintf(fp, "%8d\n", npnt);
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
                    xyz[ 3*ipnt], xyz[ 3*ipnt+1], xyz[ 3*ipnt+2],
                    Vels[3*ipnt], Vels[3*ipnt+1], Vels[3*ipnt+2],
                    uv[    ipnt]);
        }
    }

    /* write Faces to the file */
    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
        status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                &npnt, &xyz, &uv, &pindx, &ptype,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        status = ocsmGetTessVel(MODL, ibody, OCSM_FACE, iface, &Vels);
        CHECK_STATUS(ocsmGetTessVel);

        fprintf(fp, "%8d %8d\n", npnt, ntri);
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %8d %8d\n",
                    xyz[ 3*ipnt], xyz[ 3*ipnt+1], xyz[ 3*ipnt+2],
                    Vels[3*ipnt], Vels[3*ipnt+1], Vels[3*ipnt+2],
                    uv[  2*ipnt], uv[  2*ipnt+1],
                    ptype[ ipnt], pindx[ ipnt  ]);
        }
        for (itri = 0; itri < ntri; itri++) {
            fprintf(fp, "%8d %8d %8d %8d %8d %8d\n",
                    tris[3*itri], tris[3*itri+1], tris[3*itri+2],
                    tric[3*itri], tric[3*itri+1], tric[3*itri+2]);
        }
    }

    /* finalize the file */
    fclose(fp);

cleanup:
    return status;
}
