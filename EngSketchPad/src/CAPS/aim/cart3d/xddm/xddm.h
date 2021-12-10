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

struct xddmAttribute {
  char   *p_name;
  char   *p_value;
};
typedef struct xddmAttribute tsXddmAttr;
typedef tsXddmAttr *p_tsXddmAttr;

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

struct xddmFunctional {
  char   *p_comment;
  char   *p_id;
  char   *p_expr;
  double  val;

  int     lin;
  int     ndvs;   /* # of dv's and sensitivities */
  double *a_lin;  /* sensitivities */
  char  **pa_dvs; /* design variable names */

  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmFunctional tsXddmFun;
typedef tsXddmFun *p_tsXddmFun;

struct xddmAeroFun {  /* AeroFun element */
  char   *p_comment;
  char   *p_id;
  char   *p_options;

  int          nt;
  char       **pa_text;
  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmAeroFun tsXddmAFun;
typedef tsXddmAFun *p_tsXddmAFun;

struct xddmAnalysis {  /* analysis parameters */
  char   *p_comment;
  char   *p_id;
  char   *p_sens;
  double  val;
  double  discrErr;
  int     lin;
  int     ndvs;   /* # of dv's and sensitivities */
  double *a_lin;  /* sensitivities */
  char  **pa_dvs; /* design variable names */

  /* AeroFun Element */
  p_tsXddmAFun p_afun;
};
typedef struct xddmAnalysis tsXddmAPar;
typedef tsXddmAPar *p_tsXddmAPar;

struct xddmDesignPoint {  /* design point parameters */
  char        *p_comment;
  char        *p_id;
  char        *p_geometry;
  int          nv;       /** number of variables */
  p_tsXddmVar  a_v;      /** array of variables  */
  int          nc;       /** number of constants */
  p_tsXddmVar  a_c;      /** array of constants  */
  int          na;       /** number of analysis params */
  p_tsXddmAPar a_ap;     /** array of analysis params */
  p_tsXddmFun  p_obj;    /** objective functional */
  int          ncr;      /** number of constraint functionals */
  p_tsXddmFun  a_cr;     /** array of constraint functionals */
  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmDesignPoint tsXddmDesP;
typedef tsXddmDesP *p_tsXddmDesP;

struct xddmComponent {
  char   *p_comment;
  char   *p_name;
  char   *p_parent;
  char   *p_type;
  char   *p_data;
  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmComponent tsXddmComp;
typedef tsXddmComp *p_tsXddmComp;

struct xddmElement {
  char        *p_comment;
  int          nAttr;
  p_tsXddmAttr p_attr;
};
typedef struct xddmElement tsXddmElem;
typedef tsXddmElem *p_tsXddmElem;

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
  p_tsXddmElem p_config; /** configure */
  p_tsXddmElem p_inter;  /** intersect */
  int          nv;       /** number of variables */
  p_tsXddmVar  a_v;      /** array of variables  */
  int          nc;       /** number of constants */
  p_tsXddmVar  a_c;      /** array of constants  */
  int          na;       /** number of analysis params */
  p_tsXddmAPar a_ap;     /** array of analysis params */
  int          nd;       /** number of design points */
  p_tsXddmDesP a_dp;     /** array of design points */
  int          ncmp;     /** number of components */
  p_tsXddmComp a_cmp;    /** array of components */
  int          nafun;    /** number of AeroFun */
  p_tsXddmAFun a_afun;   /** array of AeroFun */
  int          ntess;    /** number of tessellates */
  p_tsXmTess   a_tess;   /** array of tessellate params */
};
typedef struct xddm tsXddm;
typedef tsXddm *p_tsXddm;

/* public functions */

p_tsXddm xddm_alloc(void);
void xddm_free(/*@null@*/ /*@only@*/ p_tsXddm p_xddm);
void xddm_echo(const p_tsXddm p_xddm);

p_tsXddmAttr xddm_allocAttribute(const int n);
p_tsXddmElem xddm_allocElement(const int n);
p_tsXddmVar  xddm_allocVariable(const int n);
p_tsXddmFun  xddm_allocFunctinoal(const int n);
p_tsXddmAFun xddm_allocAeroFun(const int n);
p_tsXddmAPar xddm_allocAnalysis(const int n);
p_tsXddmDesP xddm_allocDesignPoint(const int n);
p_tsXmTess   xddm_allocTessellate(const int n);
p_tsXddmComp xddm_allocComponent(const int n);

int xddm_addAttribute(const char *const name, const char *const value, int *n, p_tsXddmAttr *p_a);
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

p_tsXddm xddm_readFile(const char *const p_fileName,
                       const char *const xpathExpr,
                       const int options);

int xddm_writeFile(const char *const p_fileName, const
                   p_tsXddm p_xddm,
                   const int options);

int xddm_updateAnalysisParams(const char *const p_fileName,
                              const p_tsXddm p_xddm,
                              const int options);

#endif
