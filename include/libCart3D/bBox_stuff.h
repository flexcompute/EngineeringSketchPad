/*
 *
 */

#ifndef __BBOX_STUFF_H_
#define __BBOX_STUFF_H_

/*
 * $Id: bBox_stuff.h,v 1.5 2005/06/20 20:07:31 smurman Exp $
 */

#include "c3d_global.h"
#include "cartCells.h"
#include "geomTypes.h"
#include "GridInfo.h"

typedef struct BBoxStructure tsBBox;    /* --  a bounding box struct---------*/
typedef tsBBox  *p_tsBBox;

struct BBoxStructure {
  double x[DIM*2];
};


void init_BBox(p_tsBBox object_BBox, const int nKeys);

void copy_BBox(const p_tsBBox source_BBox, p_tsBBox dest_BBox);

void get_hex_BBox(p_tsBBox object_BBox, const p_tsTinyHex hex,
		  const p_tsGinfo p_Ginfo);

void get_tri_BBox(p_tsBBox p_box,  const p_tsVertex p_V0,
		  const p_tsVertex p_V1, const p_tsVertex p_V2); 

/* anisotropic scaling */
void scale_BBox(p_tsBBox p_box, const double scale_x, const double scale_y,
		const double scale_z);

/* true if they intersect (touching == true) */
bool intersect_BBox(const p_tsBBox p_box1, const p_tsBBox p_box2);

/* union of box1 and box2.  result is returned in box1 */
void union_BBox(p_tsBBox p_box1, const p_tsBBox p_box2);

/* true if it contains the point (touching == true) */
bool BBox_contains(const p_tsBBox p_box1, const double x, const double y,
                   const double z);

#endif /* __BBOX_STUFF_H_ */
