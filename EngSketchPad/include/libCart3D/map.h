/*
 *
 */

#ifndef __MAP_H_
#define __MAP_H_

/*
 * $Id: map.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 */

#include "gridStructures.h"

int* mk_CV2HexMap(const p_tsTinyHex a_Cells, const p_tsCutCell a_cCells,
		  const int nVolHexes, const int nCutHexes, const int nSplitCells);

#endif /* __MAP_H_ */
