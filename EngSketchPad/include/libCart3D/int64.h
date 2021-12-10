/*
 * $Id: int64.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 *
 */
/*                                                     M. Aftosmis, Apr 96
  
           ...begin       int64.h    checks machine info required to define
			             the INT64 and INT64_FMT macros to the
				     proper types */

#ifndef __INT64_H_
#define __INT64_H_

#ifdef __cplusplus
extern "C" {
#endif
  
#include <limits.h>                           /*     info on sizes and limits */
#ifndef _WIN32
#include <sys/param.h>                        /* #  bits per byte(NBBY) etc.. */
#endif
                                /*    ---- synonoms for machine names ---- */
#if defined(SUN) | defined(sun) | defined( __sun__)
#  if !defined (SUN)
#     define SUN
#  endif
#  if !defined (_LONGLONG)
#     define _LONGLONG
#  endif
#endif
  
                    /*     make sure  _WORD_BIT and WORD_BIT are synonomous */
#if defined (_WORD_BIT) && !defined (WORD_BIT)
#   define WORD_BIT _WORD_BIT
#endif
#if defined (WORD_BIT) && !defined (_WORD_BIT)
#   define _WORD_BIT WORD_BIT
#endif
                   /*         make sure NBBY and BITSPERBYTE are synonomous */
#if defined (NBBY) && !defined (BITSPERBYTE)
#   define BITSPERBYTE NBBY
#endif

#if !defined (WORD_BIT) && defined (SUN)
#   define WORD_BIT 32
#endif  

  /* #if !defined (WORD_BIT) && defined (BITSPERBYTE)
#   define WORD_BIT BITSPERBYTE*sizeof(int)
            Shishir: you can't evaluate sizeof(int) during preprocessing step */
#if !defined (_WORD_BIT) && defined (BITSPERBYTE) && defined (SIZEOF_INT)
#   define WORD_BIT BITSPERBYTE*SIZEOF_INT
#else
# if !defined (WORD_BIT)
#   ifdef _WIN64
#     define   WORD_BIT 64
#   else
#     define   WORD_BIT 32         /* take a guess that were on a 32 bit mach */
#   endif 
#ifndef _WIN32
# pragma ident "==> guesing at bits/word = 32, correct me if Im wrong"
#endif
# endif
#endif

/* Shishir: need this to define INT64 */
#if defined (WORD_BIT) && !defined (_WORD_BIT)
#   define _WORD_BIT WORD_BIT
#endif

                                        /*       Define INT64 nad INT64_FMT  */

#                                       /*  (1) were on a 32 bit w/ long long */
#
#if defined(_LONGLONG)  && (_WORD_BIT == 32)
#  define INT64     unsigned long long int
#  if !defined(SUN)
#      define INT64_FMT  "%lld"
#  else
#      define INT64_FMT  "%ld"  /* sun uses plain old long fmt */
#  endif
#endif
                                        /*  (2)  were on a 64 bit machine */
#if (_WORD_BIT == 64)                 
#  if defined(DEC)
#    define INT64      unsigned long
#    define INT64_FMT  "%d"
#  else
#    define INT64      unsigned long     /* ...used to be 'unsigned int' */
#    define INT64_FMT  "%d"
#  endif
#endif

#if ((!defined(_NO_LONGLONG))  && (_WORD_BIT == 32)) && !defined(INT64)
#  if defined(_WIN32)
#    define INT64      unsigned long
#    define INT64_FMT  "%ld"
#  else
#    define INT64    unsigned long long int
#    define INT64_FMT  "%lld"
#  endif
#endif

#if (!defined(INT64))
#  pragma ident  "==> int64.h real problems no way to get 64 bit integers"
#endif
  

#ifdef __cplusplus
}
#endif                         /*                       ...end int64.h */

/*-----------------------------------------------------------------------------
  NOTES:
    o   I think that this is going to need revision for some machines,
        but for now its a good starting point.
    o   there is a chance that we could use 8 contiguous characters to
        get the 64 bits we need for this but it would all just be a re-def
	of INT64 anyway, so were safe in doing it like this for now..
    o   seems to work fine on (1)  SGI IRIX 5.3/6.0 32/64 bit
                              (2)  SUN SunOS 4.1.3
			      (3)  UNICOS 8.0.4
	                                               -Michael Aftosmis, Apr 96
------------------------------------------------------------------------------*/

#endif /* __INT64_H_ */
