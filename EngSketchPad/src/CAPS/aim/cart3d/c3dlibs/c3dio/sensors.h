
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: sensors.h,v 1.6 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __SENSORINFO_H_
#define __SENSORINFO_H_

#include <string.h>

#include "c3d_global.h"
#include "basicTypes.h"
#include "stateVector.h"
#include "memory_util.h"

/* return error code */
#define C3D_SENSOR_NOT_FOUND -1

/* sensor types: equivalent area (EA) sensor is frequently used in supersonics
 * and is very similar to line sensor
 */
typedef enum{UNSET_SENSOR, POINT_SENSOR, LINE_SENSOR, EA_SENSOR} sensorType;
typedef enum{POST_PROC, CONV_HIS} sensorInfo;

typedef struct sensorDataStructure tsSensorData;
typedef tsSensorData               *p_tsSensorData;

struct sensorDataStructure {
  double  val;    /* sensor data value, eg (p-pinf)/pinf */
  double  dl;     /* line segment length in cell */
  tsState Up;     /* prims at line seg centroid */
  double  T;      /* temperature */
  double  ssp;    /* sound speed */
  dpoint3 loc;    /* line centroid location or point loc (x,y,z) */
};

typedef struct fieldSensorStructure tsFieldSensor;
typedef tsFieldSensor               *p_tsFieldSensor;

struct fieldSensorStructure {
  sensorType     type;
  sensorInfo     info;
  char           name[STRING_LEN]; /* user defined name */
  dpoint3        orig;   /* line origin */
  dpoint3        dest;   /* destination */
  int            maxR;   /* maximum cell refinement along sensor */ 
  int            nSegs;
  p_tsSensorData a_segs; /* sensor line segments */
  double         radius; /* distance from sensor to body for Ae functionals */
  double         eac;    /* constant for Ae functionals */
};

typedef struct SensorStructure tsSensor;
typedef tsSensor               *p_tsSensor;

struct SensorStructure {
  int             nSensors;
  bool            convHisMonitor;
  p_tsFieldSensor a_sensors;
};

/* prototypes */

#ifdef __cplusplus
extern "C" {
#endif
  p_tsSensor C3D_newSensors(const int nSensors);
  unsigned   C3D_freeSensors(p_tsSensor p_sensor);
  unsigned   C3D_line_cube_intersect(const dpoint3 orig, const dpoint3 dest);
  unsigned   C3D_getCode(const dpoint3 pt, const dpoint3 target,
                         const unsigned *const faceCode);
  unsigned   C3D_ck_pt_face(const dpoint3 O, const dpoint3 D, dpoint3 pt,
                            const double *const target, const unsigned int mask,
                            const double fraction, const unsigned *const faceCode);
#ifdef __cplusplus
}
#endif

#endif /* __SENSORINFO_H_ */
