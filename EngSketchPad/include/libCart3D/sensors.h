/*
 * $Id: sensors.h,v 1.4 2014/07/28 19:14:28 mnemec Exp $
 */

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

p_tsSensor C3D_newSensors(const int nSensors);
unsigned   C3D_freeSensors(p_tsSensor p_sensor);
unsigned   C3D_line_cube_intersect(const dpoint3 orig, const dpoint3 dest);
unsigned   C3D_getCode(const dpoint3 pt, const dpoint3 target,
                       const unsigned *const faceCode);
unsigned   C3D_ck_pt_face(const dpoint3 O, const dpoint3 D, dpoint3 pt,
                          const double *const target, const unsigned int mask,
                          const double fraction, const unsigned *const faceCode);

#endif /* __SENSORINFO_H_ */
