/*
 *   hex.h
 */

#ifndef __HEX_H_
#define __HEX_H_

/*
 * $Id: hex.h,v 1.3 2005/03/29 20:12:40 smurman Exp $
 */

#include "int64.h"
#include "basicTypes.h"

/* -----prototypes of public and inteface routines --- */

bool within(const INT64 p_a_name, const char *p_refA,
	    const INT64 p_b_name, const char *p_refB,
	    const int *p_Rmax, const int nBits);
  

/*    parentName()-----...returns parent int64 for this hex
                         (this is the name of the cell which this
                          hex coarsens into)------------------------------ */
INT64 parentName(const INT64 name, const char *refLevel,
		 const int nBits, const int *maxAllowRef,
		 const int *nLevels);

/* unpack the packed vertex into a dpoint3 */
void getVertex3(const INT64 vertex, const int nBits,
		const float min[DIM], const double dx[DIM], double geom[DIM]);

#endif /* __HEX_H_ */
