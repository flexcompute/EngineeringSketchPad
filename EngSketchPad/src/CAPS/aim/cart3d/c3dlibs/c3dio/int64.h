
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
 * $Id: int64.h,v 1.2 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

/**
 * checks machine info required to define the INT64 and INT64_FMT
 * macros to the proper types
 */

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

#if defined (WORD_BIT) && !defined (_WORD_BIT)
#   define _WORD_BIT WORD_BIT
#endif

                                        /*       Define INT64 and INT64_FMT  */

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
#endif

/*-----------------------------------------------------------------------------
  NOTES:
  o   May need revision for some machines, but it is a good starting point
  o   there is a chance that we could use 8 contiguous characters to
      get the 64 bits we need for this but it would all just be a re-def
	    of INT64 anyway, so were safe in doing it like this for now
------------------------------------------------------------------------------*/

#endif /* __INT64_H_ */
