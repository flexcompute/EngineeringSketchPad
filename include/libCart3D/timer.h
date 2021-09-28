/*
 * $Id: timer.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 */

#ifndef __TIMER_H_
#define __TIMER_H_

/* dtime(p) outputs the elapsed time seconds in p[1] */
/* from a call of dtime(p) to the next call of       */
/* dtime(p). */
/* p should have at least 3 elements */
int dtime(double p[]);

#endif /* __TIMER_H_ */
