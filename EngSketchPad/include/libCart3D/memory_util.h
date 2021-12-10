/*
 *                    -----    MEMORY MANAGEMENT      -----
 *
 * $Id: memory_util.h,v 1.10 2015/01/23 15:33:22 maftosmi Exp $
 *
 */


#ifndef __MEMORY_UTIL_H_
#define __MEMORY_UTIL_H_

#include <stdio.h>
#include <stdlib.h>

#include "basicTypes.h"

#define MEMORY_ERROR   -11


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

/* malloc takes a size_t argument */
typedef size_t alloca_t;

#define NEW_ALLOCA(p,type,_n) NEW_ARRAY(p,type,_n)

#define FREE_ALLOCA(p,_nBytes) FREE(p, _nBytes)
 
#endif  /* end of alloca pre-processing */


#endif /* __MEMORY_UTIL_H_ */
