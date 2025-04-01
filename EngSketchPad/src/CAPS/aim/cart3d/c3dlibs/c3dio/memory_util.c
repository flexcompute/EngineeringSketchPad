
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
 * $Id: memory_util.c,v 1.3 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#include "memory_util.h"

bool ResizeArray(void **ppv, size_t sizeNew) {                                            
  byte **ppb = (byte **) ppv;
  byte *pbNew;

  ASSERT(NULL != ppb && 0 != sizeNew);                       /* idiot checks */
  
  pbNew = (byte *) realloc (*ppb, sizeNew);
 
  /* did the block move? */
  /*
    #ifdef DEBUG
    if (pbNew != *ppb){
    NOTE("ResizeMemory MOVED block of %.3g Mb\n",
    sizeNew/(float)ONE_MB);
    }
    #endif
  */
  
  if (pbNew != NULL) *ppv = (void*) pbNew;
    
  return((bool)(pbNew != NULL));
}
/* ---------------------------------------------------------------------------
   NOTES:
   o   this is set up so that if the realloc fails, it just leaves the
   original block undisturbed, rather than assign it to NULL and
   kill the program. Its a modified version of the wrapper that
   S.Maguire recommends in Writing Solid Code
   ---------------------------------------------------------------------------*/

void _Assert(const char *strFile,const unsigned uLine) {
  fflush(stdout);
  fprintf(stderr,"\n Assertion failed: %s, line %u, exiting with status = %d\n",
          strFile,uLine,ASSERT_ERROR);
  fflush(stderr);
  strFile = NULL;                /* just to prevent compiler from complaining */
  (void) exit(ASSERT_ERROR);
}
