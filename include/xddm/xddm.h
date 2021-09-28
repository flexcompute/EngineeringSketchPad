/*
 * $Id: xddm.h,v 1.2 2014/05/07 23:55:32 mnemec Exp $
 */

/**
 * Data-structures and public functions of XDDM library
 *
 * Marian Nemec, STC, June 2013
 */

#ifndef __XDDM_H_
#define __XDDM_H_

#define XDDM_VERBOSE 1

#ifndef UNSET
#define UNSET       -888888
#endif

#define MAX_STR_LEN       4096

#ifndef ERR
#define ERR printf(" ===> ERROR:  ");printf /* standardize IO msgs */
#endif

#ifndef WARN
#define WARN printf(" ===> WARNING:  ");printf /* standardize IO msgs */
#endif

struct xddmVariable {
  char   *p_comment;
  char   *p_id;
  double  val;
  double  typicalSize;
  double  minVal;
  double  maxVal;
};
typedef struct xddmVariable tsXddmVar;
typedef tsXddmVar *p_tsXddmVar;

struct xddmAnalysis {  /* analysis parameters */
  char   *p_comment;
  char   *p_id;
  double  val;
  int     lin;
  int     ndvs;   /* # of dv's and sensitivities */
  double *a_lin;  /* sensitivities */
  char  **pa_dvs; /* design variable names */
};
typedef struct xddmAnalysis tsXddmAPar;
typedef tsXddmAPar *p_tsXddmAPar;

struct xddmAttribute {
  char   *p_name;
  char   *p_value;
};
typedef struct xddmAttribute tsXddmAttr;
typedef tsXddmAttr *p_tsXddmAttr;

struct xddmTessellate {
  char        *p_comment;
  char        *p_id;
  int          lin;
  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmTessellate tsXmTess;
typedef tsXmTess *p_tsXmTess;

struct xddmParent {
  char         *name;
  int           nAttr;
  p_tsXddmAttr  p_attr;
};
typedef struct xddmParent tsXmParent;
typedef tsXmParent *p_tsXmParent;

struct xddm {
  char        *fileName;
  char        *xpathExpr;
  p_tsXmParent p_parent; /** parent data: element name and attributes */
  int          nv;       /** number of variables */
  p_tsXddmVar  a_v;      /** array of variables  */
  int          nc;       /** number of constants */
  p_tsXddmVar  a_c;      /** array of constants  */
  int          na;       /** number of analysis params */
  p_tsXddmAPar a_ap;     /** array of analysis params */
  int          ntess;    /** number of tessellates */
  p_tsXmTess   a_tess;   /** array of tessellate params */
};
typedef struct xddm tsXddm;
typedef tsXddm *p_tsXddm;

/* public functions */

p_tsXddm xddm_alloc(void);
void xddm_free(p_tsXddm p_xddm);

p_tsXddm xddm_readFile(const char *const p_fileName,
                       const char *const xpathExpr,
                       const int options);

void xddm_echo(const p_tsXddm p_xddm);

int xddm_writeFile(const char *const p_fileName, const
                   p_tsXddm p_xddm,
                   const int options);

int xddm_updateAnalysisParams(const char *const p_fileName,
                              const p_tsXddm p_xddm,
                              const int options);

#endif
