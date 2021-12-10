/*
 *
 */

#ifndef __CONVINFO_H_
#define __CONVINFO_H_

/*
 * $Id: ConvInfo.h,v 1.5 2007/09/21 14:26:00 aftosmis Exp $
 */

/* #include "IOinfo.h"  <-- removed 07.09.21 no dependency here
                            and IOinfo breaks clients of ConvInfo.h 
                            Include IOinfo separately if you need it. M.A. */
#include "sensors.h"

typedef struct ConvergenceHistoryInfoStructure tsHinfo;
typedef tsHinfo *p_tsHinfo;

struct  ConvergenceHistoryInfoStructure {
  double     initResidual;
  double   targetResidual;
  double     timeArray[3];              /* ...for use with the dtime() timer */
  double   elapsedCPUtime;
  int               iHist;              /* How often to Update HistoryFile?  */
  int              iForce;              /* How often to update Force&MomFile */
  double          refArea;              /* the reference area for force calcs*/
  double        refLength;           /* the reference Length for moment calcs*/
  double            xStop; /* ...X-loc to stop mom. calc(if used) aka Xsting */
  double           xStart; /* ...X-loc to start mom. calc(if used) aka Xsting*/
  double      netForce[3];
  double     momentCtr[3];              /* Moment center for moment calc.    */
  double     netMoment[3];
  p_tsSensor    p_sensors;                    /* field point or line sensors */
};


#endif /* __CONVINFO_H_ */
