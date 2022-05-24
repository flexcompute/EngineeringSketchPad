#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egads.h"
#include "capsTypes.h"
#include "aimUtil.h"

/* these are defined in STATE.INC */
#define NSTATI  31
#define NSTATS 203
#define NBX      1
#define NSCTX    4
#define NBITX  (79/30 + 1)


typedef struct {
  char   code[33];
  char   name[33];
  int    kalfa;
  int    kmach;
  int    kreyn;
  int    ldepma;
  int    ldepre;
  double alfa;
  double mach;
  double reyn;
  double dnrms;
  double drrms;
  double dvrms;
  double dnmax;
  double drmax;
  double dvmax;
  int    ii;
  int    nbl;
  int    nmod;
  int    npos;
  int    *ileb;
  int    *iteb;
  int    *iend;
  double *xleb;
  double *yleb;
  double *xteb;
  double *yteb;
  double *sblegn;
  double cl;
  double cm;
  double cdw;
  double cdv;
  double cdf;
  double al_alfa;
  double cl_alfa;
  double cm_alfa;
  double cdw_alfa;
  double cdv_alfa;
  double cdf_alfa;
  double al_mach;
  double cl_mach;
  double cm_mach;
  double cdw_mach;
  double cdv_mach;
  double cdf_mach;
  double al_reyn;
  double cl_reyn;
  double cm_reyn;
  double cdw_reyn;
  double cdv_reyn;
  double cdf_reyn;
  double **xbi;
  double **ybi;
  double **cp;
  double **hk;
  double **cp_alfa;
  double **hk_alfa;
  double **cp_mach;
  double **hk_mach;
  double **cp_reyn;
  double **hk_reyn;
  double *modn;
  double *al_mod;
  double *cl_mod;
  double *cm_mod;
  double *cdw_mod;
  double *cdv_mod;
  double *cdf_mod;
  double ***gn;
  double ***xbi_mod;
  double ***ybi_mod;
  double ***cp_mod;
  double ***hk_mod;
  int    *nposel;
  int    **nbpos;
  double **xbpos;
  double **ybpos;
  double **abpos;
  double *posn;
  double *al_pos;
  double *cl_pos;
  double *cm_pos;
  double *cdw_pos;
  double *cdv_pos;
  double *cdf_pos;
  double **xbi_pos;
  double **ybi_pos;
  double **cp_pos;
  double **hk_pos;
} msesSensx;

typedef struct {
  char   name[33];
  int    istate[NSTATI];
  double sstate[NSTATS];
  int    jbld[NBX];
  int    ninl[NBX];
  int    nout[NBX];
  int    nbld[NBX];
  int    iib[NBX];
  int    ible[NBX];
  double *mfract;
  double **x;
  double **y;
  double **r;
  double **h;
  double **xb;
  double **yb;
  double **xpb;
  double **ypb;
  double **sb;
  double **sginl;
  double **sgout;
  double **xw;
  double **yw;
  double **wgap;
  double **vcen;
  double **sg;
  double **disp;
  double **pspec;
  double **thet;
  double **dstr;
  double **uedg;
  double **ctau;
  double **entr;
  double **tauw;
  double **dint;
  double **tstr;
  int    *nfreq;
  double **freq;
  double **famp;
  double **alfr;
  int    *knor;
  double **snor;
  double **xnor;
  double **xsnor;
  double **ynor;
  double **ysnor;
  double blift[NBX];
  double bdrag[NBX];
  double bmomn[NBX];
  double bdragv[NBX];
  double bdragf[NBX];
  double sble[NBX];
  double sblold[NBX];
  double swak[NBX];
  double sbcmax[NBX];
  double sbnose[NBX];
  double xbnose[NBX];
  double ybnose[NBX];
  double xbtail[NBX];
  double ybtail[NBX];
  double pxx0[2*NBX];
  double pxx1[2*NBX];
  double xtr[2*NBX];
  double xitran[2*NBX];
  double *cl_mod;
  double *cm_mod;
  double *cdw_mod;
  double *cdv_mod;
  double *cdf_mod;
  double *modn;
  double *dmspn;
  double *cl_pos;
  double *cm_pos;
  double *cdw_pos;
  double *cdv_pos;
  double *cdf_pos;
  double *posn;
  double *dpspn;
  int    igfix[2*NBX];
  int    igcorn[2*NBX];
  int    itran[2*NBX];
  int    issuct[NSCTX];
  double cqsuct[NSCTX];
  double sgsuct[NSCTX][2];
  int    *isbits[NBITX];
} msesMdat;


extern void msesSensxFree(msesSensx **sensx);
extern int  msesSensxRead(void *aimInfo, const char *fname, msesSensx **sensx);

extern void msesMdatFree(msesMdat **mdat);
extern int  msesMdatRead(void *aimInfo, const char *fname, msesMdat **mdat);
