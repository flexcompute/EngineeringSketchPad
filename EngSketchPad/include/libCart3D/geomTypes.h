/*
 * $Id: geomTypes.h,v 1.2 2011/07/26 17:30:57 mnemec Exp $
 *
 */

/* -----------
   geomtypes.h
   -----------
                ...typedefs and structure defs for basic geometric entities
*/   

#ifndef __GEOMTYPES_H_
#define __GEOMTYPES_H_

#include "basicTypes.h"

#define X  0                   /*      ...dimension names */
#define Y  1
#define Z  2

#define V0 0                   /* ...vertices in tsVertex */
#define V1 1
#define V2 2

typedef int      iquad[4];              /* --  define a indexed quad --------*/
typedef iquad   *p_iquad;

typedef struct VertexStructure tsVertex;/*--  a vertex ----------------------*/
typedef tsVertex * p_tsVertex;

struct VertexStructure {
  float x[DIM];
};

typedef struct DoublePrecisionVertexStructure tsDPVertex;
typedef tsDPVertex * p_tsDPVertex;
/* Double precision vertices - added so we can store and work with double
 * precision geometry, MN July 2011
 */
struct DoublePrecisionVertexStructure {
  double x[DIM];
};

typedef struct TriStructure tsTri;      /*--   an annotated triangle --------*/
typedef tsTri * p_tsTri;

struct TriStructure { 
 int      vtx[3];                             /* limits compnums to +-32767  */
 int      Comp:16;                            /* 2 because mtype is an enum  */
 mtype    mark:2, mark2:2;                    /* mark/unmark... etc          */
 unsigned :0;                                 /* pad to word boundary        */
};

/* ------ */

#endif /* __GEOMTYPES_H_ */
