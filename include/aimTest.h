/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifndef AIM_TEST_H_
#define AIM_TEST_H_

#include "acutest.h"
#include "capsErrors.h"

#define TEST_CHECK_SUCCESS( statTest ) \
    TEST_CHECK_(CAPS_SUCCESS == statTest, "expected status CAPS_SUCCESS != %s", #statTest)

#define TEST_CHECK_STATUS( statTrue, statTest ) \
    TEST_CHECK_(statTrue == statTest,"expected status %d != %s", statTrue, #statTest)

#define TEST_ASSERT_EQUAL_INT( expected, actual ) \
    TEST_ASSERT_( expected == actual, "expected %d != %d", expected, actual )

#define TEST_CHECK_EQUAL_INT( expected, actual ) \
    TEST_CHECK_( expected == actual, "expected %d != %d", expected, actual )

#define TEST_CHECK_EQUAL_PTR( expected, actual ) \
    TEST_CHECK_( expected == actual, "pointers do not match" )

#define TEST_CHECK_EQUAL_DBL( expected, actual ) \
    TEST_CHECK_( expected == actual, "expected %f != %f", expected, actual )

#define TEST_CHECK_STRING( strTrue, strTest ) \
    if (strTrue == NULL) { TEST_CHECK_(0, "Null string!"); return; } \
    if (strTest == NULL) { TEST_CHECK_(0, "Null string!"); return; } \
    TEST_CHECK_( strcmp(strTrue, strTest) == 0, "expected %s != %s", strTrue, strTest )

#define TEST_CHECK_CLOSE( expected, actual, abserr ) \
    TEST_CHECK_( fabs(expected - actual) < abserr, "expected %e != %e within %e", expected, actual, abserr)

#endif // AIM_TEST_H_
