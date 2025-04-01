
#ifndef _CAPS_AIM_CART3D_XDDM_UTIL_H_
#define _CAPS_AIM_CART3D_XDDM_UTIL_H_

#include "xddm.h"

int xddm_addAeroFunForce(p_tsXddmAFun p_afun,
                         const char *const name,
                         const int force,
                         const int frame,
                         const int J,
                         const int N,
                         const double target,
                         const double weight,
                         const int bnd,
                         const char *const comp);
int xddm_addAeroFunMoment_Point(p_tsXddmAFun p_afun,
                                const char *const name,
                                const int index,
                                const int moment,
                                const int frame,
                                const int J,
                                const int N,
                                const double target,
                                const double weight,
                                const int bnd,
                                const char *const comp);
int xddm_addAeroFunLoD(p_tsXddmAFun p_afun,
                       const char *const name,
                       const int frame,
                       const int J,
                       const int N,
                       const double A,
                       const double bias,
                       const double target,
                       const double weight,
                       const int bnd,
                       const char *const comp);


#endif //_CAPS_AIM_CART3D_XDDM_UTIL_H_
