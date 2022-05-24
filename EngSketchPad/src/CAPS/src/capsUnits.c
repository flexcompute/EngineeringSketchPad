/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Units Utility Functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "capsBase.h"

#include "udunits.h"

/*@-incondefs@*/
extern void ut_free(/*@only@*/ ut_unit* const unit);
extern void cv_free(/*@only@*/ cv_converter* const converter);
/*@+incondefs@*/

#define UNIT_BUFFER_MAX 257

static ut_system *utsystem = NULL;


static void caps_free_utsystem(void)
{
  ut_free_system(utsystem);
  utsystem = NULL;
}


/*@null@*/
void *caps_initUnits()
{
  char *env;
#ifdef WIN32
  int  oldrive, drive;
#endif

  if (utsystem != NULL) {
    return utsystem;
  }

#ifdef WIN32
  oldrive = _getdrive();
  env     =  getenv("UDUNITS2_XML_PATH");
  if (env != NULL)
    if (env[1] == ':') {
      if (env[0] >= 97) {
        drive = env[0] - 96;
      } else {
        drive = env[0] - 64;
      }
      _chdrive(drive);
    }
#endif

  /* suppress error messages */
  (void) ut_set_error_message_handler(ut_ignore);

/*@-nullpass@*/
  utsystem = ut_read_xml(NULL);
/*@+nullpass@*/

  if (utsystem == NULL) {
    printf("Error: Failed to load UDUNITS XML file!\n"
           "       This might be resolved by setting the environment variable:\n"
           "       UDUNITS2_XML_PATH=$ESP_ROOT/../udunits/udunits2.xml\n");
    env = getenv("UDUNITS2_XML_PATH");
    if (env != NULL) {
      printf("       Currently: UDUNITS2_XML_PATH=%s\n",env);
    }
  }

#ifdef WIN32
  _chdrive(oldrive);
#endif
  atexit(caps_free_utsystem);

  return utsystem;
}


int
caps_convert(const int count, /*@null@*/ const char  *inUnit, double  *inVal,
                              /*@null@*/ const char *outUnit, double *outVal)
{
  int          i;
  ut_unit      *utunit1, *utunit2;
  cv_converter *converter;

  if (outVal  == NULL) return CAPS_NULLVALUE;

  /* check for simple equality */
  if ((inUnit == NULL) && (outUnit == NULL)) {
    for (i = 0; i < count; i++) outVal[i] = inVal[i];
    return CAPS_SUCCESS;
  }
  if (inUnit  == NULL) return CAPS_UNITERR;
  if (outUnit == NULL) return CAPS_UNITERR;
  if (strcmp(outUnit, inUnit) == 0) {
    for (i = 0; i < count; i++) outVal[i] = inVal[i];
    return CAPS_SUCCESS;
  }

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1   = ut_parse(utsystem, outUnit, UT_ASCII);
  utunit2   = ut_parse(utsystem, inUnit,  UT_ASCII);
  converter = ut_get_converter(utunit2, utunit1); /* From -> To */
  if (converter == NULL || ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit2);
    return CAPS_UNITERR;
  }

  (void)cv_convert_doubles(converter, inVal, count, outVal);
  cv_free(converter);
  ut_free(utunit2);
  ut_free(utunit1);

  if (ut_get_status() != UT_SUCCESS) {
    return CAPS_UNITERR;
  }

  return CAPS_SUCCESS;
}


int
caps_unitParse(/*@null@*/ const char *unit)
{
  ut_unit *utunit;

  if (unit == NULL) return CAPS_SUCCESS;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit = ut_parse(utsystem, unit, UT_ASCII);
  if (utunit == NULL || ut_get_status() != UT_SUCCESS) {
    return CAPS_UNITERR;
  }
  ut_free(utunit);

  return CAPS_SUCCESS;
}


int
caps_unitConvertible( const char *unit1, const char *unit2)
{
  int     status;
  ut_unit *utunit1, *utunit2;

  /* check for simple equality */
  if        ((unit1 == NULL) && (unit2 == NULL)) {
    return CAPS_SUCCESS;
  } else if ((unit1 == NULL) || (unit2 == NULL)) {
    return CAPS_UNITERR;
  } else if (strcmp(unit1, unit2) == 0) {
    return CAPS_SUCCESS;
  }

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, unit1, UT_ASCII);
  utunit2 = ut_parse(utsystem, unit2, UT_ASCII);
  status  = ut_are_convertible(utunit1, utunit2);
  ut_free(utunit1);
  ut_free(utunit2);

  if (status == UT_SUCCESS) { /* yeah... success means failure */
    return CAPS_UNITERR;
  }

  return CAPS_SUCCESS;
}


int
caps_unitCompare( const char *unit1, const char *unit2, int *compare)
{
  ut_unit *utunit1, *utunit2;

  if (compare  == NULL) return CAPS_NULLVALUE;
  *compare = 0;

  /* check for simple equality */
  if        ((unit1 == NULL) && (unit2 == NULL)) {
    return CAPS_SUCCESS;
  } else if ((unit1 == NULL) || (unit2 == NULL)) {
    return CAPS_UNITERR;
  } else if (strcmp(unit1, unit2) == 0) {
    return CAPS_SUCCESS;
  }

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1  = ut_parse(utsystem, unit1, UT_ASCII);
  utunit2  = ut_parse(utsystem, unit2, UT_ASCII);
  *compare = ut_compare(utunit1, utunit2);
  ut_free(utunit1);
  ut_free(utunit2);

  return CAPS_SUCCESS;
}


int
caps_unitMultiply(const char *inUnits1, const char *inUnits2,
                  char **outUnits)
{
  int         status;
  ut_unit     *utunit1, *utunit2, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnits1 == NULL) ||
      (inUnits2 == NULL) || (outUnits == NULL)) return CAPS_NULLVALUE;

  *outUnits = NULL;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, inUnits1, UT_ASCII);
  utunit2 = ut_parse(utsystem, inUnits2, UT_ASCII);
  utunit  = ut_multiply(utunit1, utunit2);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit2);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit2);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
caps_unitDivide(const char *inUnits1, const char *inUnits2,
                char **outUnits)
{
  int         status;
  ut_unit     *utunit1, *utunit2, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnits1 == NULL) ||
      (inUnits2 == NULL) || (outUnits == NULL)) return CAPS_NULLVALUE;

  *outUnits = NULL;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, inUnits1, UT_ASCII);
  utunit2 = ut_parse(utsystem, inUnits2, UT_ASCII);
  utunit  = ut_divide(utunit1, utunit2);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit2);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit2);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
caps_unitInvert(const char *inUnit, char **outUnits)
{
  int         status;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLVALUE;

  *outUnits = NULL;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, inUnit, UT_ASCII);
  utunit  = ut_invert(utunit1);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
caps_unitRaise(const char  *inUnit, const int power, char **outUnits)
{
  int         status;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLVALUE;

  *outUnits = NULL;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, inUnit, UT_ASCII);
  utunit  = ut_raise(utunit1, power);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}


int
caps_unitOffset(const char  *inUnit, const double offset, char **outUnits)
{
  int         status;
  ut_unit     *utunit1, *utunit;
  char        buffer[UNIT_BUFFER_MAX];

  if ((inUnit == NULL) || (outUnits == NULL)) return CAPS_NULLVALUE;

  *outUnits = NULL;

  if (utsystem == NULL) {
    (void) caps_initUnits();
    if (utsystem == NULL) return CAPS_UNITERR;
  }

  utunit1 = ut_parse(utsystem, inUnit, UT_ASCII);
  utunit  = ut_offset(utunit1, offset);
  if (ut_get_status() != UT_SUCCESS) {
    ut_free(utunit1);
    ut_free(utunit);
    return CAPS_UNITERR;
  }

  status = ut_format(utunit, buffer, UNIT_BUFFER_MAX, UT_ASCII);

  ut_free(utunit1);
  ut_free(utunit);

  if (ut_get_status() != UT_SUCCESS || status >= UNIT_BUFFER_MAX) {
    return CAPS_UNITERR;
  } else {
    *outUnits = EG_strdup(buffer);
    if (*outUnits == NULL) return EGADS_MALLOC;
    return CAPS_SUCCESS;
  }
}
