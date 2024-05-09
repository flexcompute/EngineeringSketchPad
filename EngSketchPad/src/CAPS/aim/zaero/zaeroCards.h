#ifndef __ZAERO_CARDS_H__
#define __ZAERO_CARDS_H__

#include "zaeroTypes.h"


int zaero_card_aeroz(void *aimInfo,
                     FILE *fp, int acsid, char *xzsym,
                     char *flip, /*@null@*/ char *fmmunit, char *fmlunit,
                     double refc, double refb, double refs,
                     double refx, double refy, double refz,
                     feaFileTypeEnum formatType);

int zaero_card_caero7(void *aimInfo,
                      FILE *fp, int wid, /*@null@*/ char *label, int acoord,
                      int nspan, int nchord, int lspan, int ztaic, int pafoil7,
                      double rl[3], double rch, int lrchd, int attchr, int acordr,
                      double tl[3], double tch, int ltchd, int attcht, int acordt,
                      feaFileTypeEnum formatType);

int zaero_card_fixmden(void *aimInfo,
                       FILE *fp, int setid, int idmk, double den,
                       char *ftmunit, char *ftlunit, double vref, int fluttf,
                       int print, int numVelocities, double *velocities,
                       feaFileTypeEnum formatType);

int zaero_card_flutter(void *aimInfo,
                       FILE *fp, int setid, char *sym, int fix, int nmode,
                       int tabdmp, int mlist, int conmlst, int nkstep,
                       feaFileTypeEnum formatType);

int zaero_card_mkaeroz(void *aimInfo,
                       FILE *fp, int idmk, double mach, int method, int idflt,
                       const char *save, const char *filenm, int print,
                       int numFreq, double freq[],
                       feaFileTypeEnum formatType);

int zaero_card_panlst2(void *aimInfo,
                       FILE *fp, int setid, int macroid, int boxBegin, int boxEnd,
                       feaFileTypeEnum formatType);

int zaero_card_pltaero(void *aimInfo,
                       FILE *fp, int setid, char *femgrid, int offset,
                       char *form, char *filenm, char *cell, char *vct,
                       feaFileTypeEnum formatType);

int zaero_card_pltflut(void *aimInfo,
                       FILE *fp, int setid, int idflut, int mode, int ntime,
                       double maxdisp, char *form, char *filenm, char *aeronm,
                       feaFileTypeEnum formatType);

int zaero_card_pltmode(void *aimInfo,
                       FILE *fp, int setid, char *sym, int mode,
                       double maxdisp, char *form, char *filenm, char *aeronm,
                       feaFileTypeEnum formatType);

int zaero_card_plttrim(void *aimInfo,
                       FILE *fp, int setid, int idtrim, char *flex, char *type,
                       char *form, char *filenm, double scale, char *aeronm,
                       feaFileTypeEnum formatType);

int zaero_card_pltvg(void *aimInfo,
                     FILE *fp, int setid, int idflut, int nmode,
                     char *xaxis, char *form, char *filenm, double refrho,
                     feaFileTypeEnum formatType);

int zaero_card_set1(void *aimInfo,
                    FILE *fp, int setid, int numIntegers, int integers[],
                    feaFileTypeEnum formatType);

int zaero_card_spline1(void *aimInfo,
                       FILE *fp, int eid, char *model, int cp,
                       int setk, int setg, double dz, double eps,
                       feaFileTypeEnum formatType);

int zaero_card_spline3(void *aimInfo,
                       FILE *fp, int eid, char *model,
                       int setk, int setg, double eps,
                       feaFileTypeEnum formatType);

int zaero_card_tabdmp1(void *aimInfo,
                       FILE *fp, int tid, char *type,
                       int numDamping,
                       double *dampingFreq,
                       double *dampingValues,
                       feaFileTypeEnum formatType);

int zaero_card_trim(void *aimInfo,
                    FILE *fp, int trimid, int idmk, double qinf,
                    int idobj, int idcons, double rho[3], double wtmass,
                    double weight, double I[6], char *trnacc,
                    char *nx, char *ny, char *nz, char *pdot, char *qdot,
                    char *rdot, int loadset, int numVars, int *varIDs,
                    char **vals,
                    feaFileTypeEnum formatType);

int zaero_card_trimvar(void *aimInfo,
                       FILE *fp, int idvar, /*@null@*/ char *label, double lower,
                       double upper, int trimlnk, /*@null@*/ char *dmi, char *sym,
                       int guessInitial, double initial, char *dc[6],
                       feaFileTypeEnum formatType);

#endif // __ZAERO_CARDS_H__
