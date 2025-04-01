
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
 * $Id: xddm.h,v 1.5 2022/11/07 23:01:39 mnemec Exp $
 */

/* open source */

/**
 * Data structures and public functions of XDDM library
 */

#ifndef __XDDM_H_
#define __XDDM_H_

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> /* for fsync */

#ifndef UNSET
#define UNSET -888888
#endif

#define MAX_STR_LEN 4096

#ifndef ERR
#define ERR printf(" ===> ERROR:  ");printf /* standardize IO msgs */
#endif

#ifndef WARN
#define WARN printf(" ===> WARNING:  ");printf /* standardize IO msgs */
#endif

#ifdef __cplusplus
extern "C" {
#endif
                                  /* bits for flags variable in xddm element */
typedef enum {
  XDDM_VERBOSE     = 1,  /* be verbose */
  XDDM_DEBUG       = 2,  /* debug on */
  XDDM_LIN         = 4,  /* linearization required */
  XDDM_NOLIN       = 8,  /* no linearization */
  XDDM_BOUND_UPPER = 16, /* constraint upper bound */
  XDDM_BOUND_LOWER = 32  /* constraint lower bound */
} xddmOptions;

typedef struct xddm tsXddm;
typedef tsXddm *p_tsXddm;

struct xddmAttribute {
  char   *p_name;
  char   *p_value;
};
typedef struct xddmAttribute tsXddmAttr;
typedef tsXddmAttr *p_tsXddmAttr;

struct xddmElement {
  char         *p_nn;      /* node name, e.g. Variable */
  char         *p_id;      /* local id */
  char         *p_xn;      /* xddm name */
  char         *p_comment;
  p_tsXddmAttr  p_attr;    /* external attributes */
  size_t        nAttr;
  unsigned      flags;     /* bit variable for configure options */
};
typedef struct xddmElement tsXddmElem;
typedef tsXddmElem *p_tsXddmElem;

struct xddmVariable {
  p_tsXddmElem p_e;
  double       val;
  double       typicalSize;
  double       minVal;
  double       maxVal;
  double       fdstep;
};
typedef struct xddmVariable tsXddmVar;
typedef tsXddmVar *p_tsXddmVar;

struct xddmAeroFun {
  p_tsXddmElem p_e;
  char *p_text;
};
typedef struct xddmAeroFun tsXddmAFun;
typedef tsXddmAFun *p_tsXddmAFun;

struct xddmGeometry {
  p_tsXddmElem p_e;
  char *p_parts;
  char *p_comp2tri;
  char *p_cutout;
  char *p_overlap;
  char *p_ps;
};
typedef struct xddmGeometry tsXddmGeom;
typedef tsXddmGeom *p_tsXddmGeom;

struct xddmAnalysis {   /* analysis parameters */
  p_tsXddmElem  p_e;
  char        **pa_dvs; /* design variable names */
  double       *a_lin;  /* sensitivities */
  p_tsXddmAFun  p_afun; /* AeroFun kid */
  double        val;
  double        derr;   /* discretization error */
  size_t        ndvs;   /* # of dv's and sensitivities */
};
typedef struct xddmAnalysis tsXddmAPar;
typedef tsXddmAPar *p_tsXddmAPar;

struct xddmFunction {
  p_tsXddmElem  p_e;
  char         *p_expr;
  double       *a_lin;  /* derivatives */
  char        **pa_dvs; /* design variable names */
  double        val;
  double        min;
  double        max;
  size_t        ndvs;   /* # of sensitivities */
};
typedef struct xddmFunction tsXddmFun;
typedef tsXddmFun *p_tsXddmFun;

struct xddmSum {
  p_tsXddmElem  p_e;
  char         *p_expr;
  double       *a_lin;  /* derivatives */
  char        **pa_dvs; /* design variable names */
  double        val;
  size_t        ndvs;   /* # of sensitivities */
};
typedef struct xddmSum tsXddmSum;
typedef tsXddmSum *p_tsXddmSum;

struct xddm {
  p_tsXddmElem p_e;     /* root data: node name and attributes */
  p_tsXddmVar  a_v;     /* array of variables  */
  p_tsXddmVar  a_c;     /* array of constants  */
  p_tsXddmAPar a_ap;    /* array of analysis params */
  p_tsXddmAFun a_afun;  /* array of AeroFun elements */
  p_tsXddmElem a_t;     /* array of tessellate params */
  p_tsXddmSum  a_s;     /* array of sums */
  p_tsXddmFun  a_f;     /* array of functions */
  p_tsXddmFun  a_j;     /* array of objectives */
  p_tsXddmFun  a_con;   /* array of constraints */
  p_tsXddmGeom a_geo;   /* array of intersect nodes */
  p_tsXddm     a_kids;  /* array of xddm kids */
  size_t       nv;      /* number of variables */
  size_t       nc;      /* number of constants */
  size_t       na;      /* number of analysis params */
  size_t       naf;     /* number of AeroFun nodes */
  size_t       nt;      /* number of tessellate nodes */
  size_t       ns;      /* number of sum functions */
  size_t       nf;      /* number of functions */
  size_t       nj;      /* number of objectives */
  size_t       ncon;    /* number of constraints */
  size_t       ng;      /* number of intersect nodes */
  size_t       nk;      /* number of kids */
};

/* public functions */

p_tsXddm xddm_new(const size_t n);
void     xddm_free(/*@null@*/ /*@only@*/ p_tsXddm p_xddm, unsigned flags);

p_tsXddmElem xddm_newElement(    const size_t n);
p_tsXddmAttr xddm_newAttribute(  const size_t n);
p_tsXddmVar  xddm_newVariable(   const size_t n);
p_tsXddmAPar xddm_newAnalysis(   const size_t n);
p_tsXddmAFun xddm_newAeroFun(    const size_t n);
p_tsXddmFun  xddm_newFunction(   const size_t n);
p_tsXddmSum  xddm_newSum(        const size_t n);
p_tsXddmGeom xddm_newIntersect(  const size_t n);

void xddm_addAttribute(const char * p_name, const char * p_value,
                       size_t *p_nAttr, p_tsXddmAttr *pp_attr);

void xddm_setID(const p_tsXddmElem p_e, const char *const p_nn,
                const char *const p_id, const char *const p_comment);
void xddm_setElementXN(const p_tsXddmElem p_e, const char *const p_xn);
void xddm_setElementFlags(const p_tsXddmElem p_e, const xddmOptions flag);

void xddm_setVariable(const p_tsXddmVar p_v, const double *const p_val,
                      const double *const p_typicalSize,
                      const double *const p_minVal,
                      const double *const p_maxVal,
                      const double *const p_fdstep);

void xddm_setFunction(p_tsXddmFun p_f, const char *const p_expr,
                      const double *const p_val, const double *const p_min,
                      const double *const p_max);

p_tsXddm xddm_readFile(const char *const p_fileName,
                       const char *const xpathExpr,
                       unsigned *p_options);

void xddm_echo(const p_tsXddm p_xddm, unsigned indent);

int xddm_writeFile(const char *const p_fileName, const
                   p_tsXddm p_xddm,
                   const unsigned options);

int xddm_updateAnalysisParams(const char *const p_fileName,
                              const p_tsXddm p_xddm,
                              const unsigned options);

#ifdef __cplusplus
}
#endif

#endif /* __XDDM_H_ */
