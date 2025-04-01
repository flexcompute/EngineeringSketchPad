
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: memory_util.h,v 1.11 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __MEMORY_UTIL_H_
#define __MEMORY_UTIL_H_

#include <stdio.h>
#include <stdlib.h>

#include "basicTypes.h"

#define MEMORY_ERROR   -11

#ifdef __cplusplus
extern "C" {
#endif

#define NEW(p,type)     if ((p=(type *)                                \
                           malloc(sizeof(type))) == NULL) {            \
                           printf ("Malloc out of Memory!\n");         \
                           exit(MEMORY_ERROR);                         \
                        }                                              \
        else {                                 \
                 }


#define FREE(p,_nBytes)  if (p) {free ((char *) p);p=NULL;  }

#define NEW_ARRAY(p,type,_n){ ASSERT(0 != (_n));                         \
           if ((p=(type *) malloc                                        \
             (((size_t)_n)*sizeof(type))) == NULL) {                     \
             printf ("Array malloc out of Memory!\n");                   \
             printf ("While trying to allocate %.2gMb\n",                \
                    (double)(sizeof(type)*(double)(_n)/(double)ONE_MB)); \
             printf ("(Malloc tried to alloc %ld items for a total of %zu bytes)\n",\
                     ((unsigned long int)(_n)),(((size_t)_n)*sizeof(type)));        \
             if ( 0 > (_n)){ERR("Negative value passed to NEW_ARRAY\n");}\
             exit(MEMORY_ERROR);                                         \
           }                                                             \
        else {                                                           \
                  }                                                      \
           }

                                          /* ...cover function for realloc */
bool ResizeArray(void **ppv, size_t sizeNew);

#define RESIZE_ARRAY(pp,type,_n) { \
        ASSERT(0 != (_n)); \
        if (!ResizeArray((void*)(pp), sizeof(type)*(size_t)(_n))) {    \
             ERR("Resizing array failed, while trying to get %.2g Mb", \
                (double)(sizeof(type)*(double)(_n)/(double)ONE_MB));   \
             printf ("(Malloc tried to alloc %ld items for a total of %zu bytes)\n",\
                     ((unsigned long int)(_n)),(((size_t)_n)*sizeof(type)));        \
             if ( 0 > (_n)){ERR("Negative value passed to NEW_ARRAY\n");}\
             exit(MEMORY_ERROR);                                         \
        }                                                                \
     }

#ifndef NO_ALLOCA

#ifndef DARWIN   /* on DARWIN alloca is in stdlib.h */
#include <alloca.h>  /* <- _ALLOC_A_ non-STD, if your sys doesnt have this
            You'll need to use malloc where we use alloca() */
#endif

#ifdef SGI
/* alloca only takes up to an unsigned int on SGI */
typedef unsigned int alloca_t;
#else
typedef size_t alloca_t;
#endif

#define NEW_ALLOCA(p,type,_n) { \
    ASSERT(0 != (_n));                             \
    if ((p= (type *)alloca((alloca_t)(_n*sizeof(type)))) == NULL) {\
      printf ("Array alloca out of Memory!\n");    \
      printf ("While trying to allocate %4.2gMb\n",\
      (sizeof(type)*(double)(_n)/(double)ONE_MB));\
      exit(MEMORY_ERROR);                          \
    } else { /* debugging clause printf("allocated %4.2gMb\n", sizeof(type)*(double)(_n)/(double)ONE_MB);*/ } \
    }
#define FREE_ALLOCA(p,_nBytes) ASSERT(0 != _nBytes);

#else  /* no stack alloca use heap malloc instead */

#define NEW_ALLOCA(p,type,_n) NEW_ARRAY(p,type,_n)

#define FREE_ALLOCA(p,_nBytes) FREE(p, _nBytes)

#endif  /* end of alloca pre-processing */

#ifdef __cplusplus
}
#endif

#endif /* __MEMORY_UTIL_H_ */
