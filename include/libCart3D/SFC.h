/*
 * $Id: SFC.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 *
 *  Interface for the space-filling curves
 */
#ifndef __SFC_H_
#define __SFC_H_

#include "int64.h"

void  output_all_positions(int level);
INT64 cell_peano_order(int level, int x, int y, int z);

INT64 cell_morton_order(int level, int x, int y, int z);

#endif

