/*
 *
 * $Id: c3d_timeout.h,v 1.8.2.1 2015/09/03 17:49:16 maftosmi Exp $
 *
 */

#ifndef __C3D_TIMEOUT_H__
#define __C3D_TIMEOUT_H__

/*                         ....predefiend variables */
#define FOREVER        999999999
#define NEVER_TIMEOUT  -1
/*         ....set time periods you want to build in */
#define GRACE_PERIOD   60
#define WARN_PERIOD    30
#define ACTIVE_PERIOD  800



/**
 *  --c3d_timeout()---...public API -----------------------------------------
 *
 *           return value : (+): Integer value with num. of days remaining
 *                           0 : Date is expired
 *
 *           sample usage:
 *               Allow execution for a period of 370 days from code build
 *
 *           if (0 == c3d_timeout(370)){
 *              printf("Time limit expired. stopping execution.\n");
 *              fflush(stdout);
 *              fflush(stderr);
 *           exit(-1);
 *
 *  ENVIRONMENT: Setting the following environment variables to non-null
 *               values provides additional functionalitya
 *               o C3D_GRACEPERIOD --Permits runnign during grace period
 *                                   after official expiration
 *                 % setenv C3D_GRACEPERIOD  1
 *
 *               o C3D_TIMEREMAINING --Report number of days till
 *                                     expiration to stdout at runtime
 *                 % setenv C3D_TIMEREMAINING 1
 *
 *
 *  COMPILATION: The compile line should put the build date into the
 *               COMP_EPOCH preprocessor constant.
 *               Here is a sample compile line:
 *
 *               % gcc -c c3d_timeout.c -DCOMP_EPOCH=`date "+%s"`
 *
 *                            -Jun'09  -Michael Aftosmis & Marian Nemec
 *                                                 NASA Ames & STC
 * -------------------------------------------------------------------- */

int c3d_timeout(const int num_days_valid);


#endif /*  __C3D_TIMEOUT_H__ */
