/*
 * general routines for handling bit-wise encoding 
 * and decoding 
 */

#ifndef __ENCODING_H_
#define __ENCODING_H_

/*
 * $Id: encoding.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 */

#include "c3d_global.h"
#include "int64.h"

/*       ---- pack 3 dimensions, x1, x2, x3 into a single "INT64" var ----

    packing scheme:
                        |<----num_bits----->|
     | --|-------------||-------------------||------------------||
     | 0 | ...x1...    ||    ...x2...       ||   ...x3...       ||
     | --|-------------||-------------------||------------------||
*/


INT64 pack(const int x1,const int x2,const int x3,const int num_bits );
/*       ---- unpack the i_dim th coord of a packed_var ----

                        |<----num_bits----->|
     | --|-------------||-------------------||------------------||
     | 0 | i_dim=0     ||    i_dim=1        ||    i_dim=2       ||
     | --|-------------||-------------------||------------------||
*/

int unpack(const int i_dim, const INT64 packed_var,const int num_bits);


void bitprint(const INT64 p, const INT64 nbits);

/*        ---- return the 8 vertices of a cartesian cell ---- */
/*
  
     A                                             Vertex Ordering:  
   Y |			 V6 O----------------O V7   ===============
     |   7X		   /:               /|        The vertex ordering scheme
     |  /  		  / :		   / |	      stems from the packing of
     | /   		 /  :		  /  |	      INT64 names, since Z occu-
     |/	       		/   :		 /   |	      pies the least sig. bits,
     o-------->Z       /    :           /    | 	      it varies the fastest.
  Origin       	   V2 O----------------O V3  |	       V0 = (0,0,0) = 000
		      |     :          |     |	       V1 = (0,0,1) = 001
		      |     :          |     |	       V2 = (0,1,0) = 010
		      |     :          |     |	       V3 = (0,1,1) = 011
		      |     O..........|.....O	       V4 = (1,0,0) = 100
		      |    . V4        |    / V5       V5 = (1,0,1) = 101
		      |   .            |   /	       V6 = (1,1,0) = 110
		      |  .             |  /	       V7 = (1,1,1) = 111
		      | .              | /
		      |.               |/
       	       	   V0 O----------------O V1
*/		   
void  _Get_verts(const INT64 v0       , INT64 v[8], const char R[DIM],
		 const int   Rmax[DIM],             const int nBits );

#endif /* __ENCODING_H_ */





