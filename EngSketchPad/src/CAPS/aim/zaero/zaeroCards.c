#include <string.h>

#include "aimUtil.h"
#include "cardUtils.h"
#include "miscUtils.h"
#include "zaeroCards.h"
#include "zaeroUtils.h"


int zaero_card_aeroz(void *aimInfo,
                     FILE *fp, int acsid, char *xzsym,
                     char *flip, char *fmmunit, char *fmlunit,
                     double refc, double refb, double refs,
                     double refx, double refy, double refz,
                     feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "AEROZ", formatType);
  AIM_STATUS(aimInfo, status);

  // ACSID (Integer > 0 or Blank)
  status = card_addInteger(&card, acsid);
  AIM_STATUS(aimInfo, status);

  // XZSYM (Character), can be blank
  status = card_addString(&card, xzsym);
  AIM_STATUS(aimInfo, status);

  // FLIP (Character), can be blank
  status = card_addString(&card, flip);
  AIM_STATUS(aimInfo, status);

  // FMMUNIT(Character), can be blank
  status = card_addString(&card, fmmunit);
  AIM_STATUS(aimInfo, status);

  // FMLUNIT (Character), can be blank
  status = card_addString(&card, fmlunit);
  AIM_STATUS(aimInfo, status);

  // REFC (Real >= 0), can be blank
  status = card_addDouble(&card, refc);
  AIM_STATUS(aimInfo, status);

  // REFB (Real >= 0), can be blank
  status = card_addDouble(&card, refb);
  AIM_STATUS(aimInfo, status);

  // REFS (Real >= 0), can be blank
  status = card_addDouble(&card, refs);
  AIM_STATUS(aimInfo, status);

  // REFX, REFY, REFZ (Real)
  status = card_addDouble(&card, refx);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, refy);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, refz);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:
  card_destroy(&card);

  return status;
}

// write MKAEROZ card
int zaero_card_mkaeroz(void *aimInfo,
                       FILE *fp, int idmk, double mach, int method, int idflt,
                       const char *save, const char *filenm, int print,
                       int numFreq, double freq[],
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "MKAEROZ", formatType);
  AIM_STATUS(aimInfo, status);

  // IDMK (Integer > 0)
  status = card_addInteger(&card, idmk);
  AIM_STATUS(aimInfo, status);

  // MACH (Real >= 0)
  status = card_addDouble(&card, mach);
  AIM_STATUS(aimInfo, status);

  // METHOD (Integer >= 0), can be blank
  status = card_addInteger(&card, method);
  AIM_STATUS(aimInfo, status);

  // IDFLT (Integer >= 0)
  status = card_addInteger(&card, idflt);
  AIM_STATUS(aimInfo, status);

  // SAVE (Characters or blank)
  status = card_addString(&card, save);
  AIM_STATUS(aimInfo, status);

  // FILENM (Characters or blank)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // PRINT (Integer)
  status = card_addInteger(&card, print);
  AIM_STATUS(aimInfo, status);

  // FREQi (Real)
  status = card_addDoubleArray(&card, numFreq, freq);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write FLUTTER card
int zaero_card_flutter(void *aimInfo,
                       FILE *fp, int setid, char *sym, int fix, int nmode,
                       int tabdmp, int mlist, int conmlst, int nkstep,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "FLUTTER", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // SYM
  status = card_addString(&card, sym);
  AIM_STATUS(aimInfo, status);

  // FIX (Integer > 0)
  status = card_addInteger(&card, fix);
  AIM_STATUS(aimInfo, status);

  // NMODE (Integer >= 0)
  status = card_addInteger(&card, nmode);
  AIM_STATUS(aimInfo, status);

  // TABDMP (Integer >= 0)
  status = card_addInteger(&card, tabdmp);
  AIM_STATUS(aimInfo, status);

  // MLIST (Integer >= 0), can be blank
  status = card_addInteger(&card, mlist);
  AIM_STATUS(aimInfo, status);

  // CONMLST (Integer >= 0), can be blank
  status = card_addInteger(&card, conmlst);
  AIM_STATUS(aimInfo, status);

  // NKSTEP (Integer >= 0), can be blank
  status = card_addInteger(&card, nkstep);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write FIXMDEN card
int zaero_card_fixmden(void *aimInfo,
                       FILE *fp, int setid, int idmk, double den,
                       char *ftmunit, char *ftlunit, double vref, int fluttf,
                       int print, int numVelocities, double *velocities,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  card_initiate(&card, "FIXMDEN", formatType);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // IDMK (Integer > 0)
  status = card_addInteger(&card, idmk);
  AIM_STATUS(aimInfo, status);

  // DEN (Real > 0.0)
  status = card_addDouble(&card, den);
  AIM_STATUS(aimInfo, status);

  // FTMUNIT (Character)
  status = card_addString(&card, ftmunit);
  AIM_STATUS(aimInfo, status);

  // FTLUNIT (Character)
  status = card_addString(&card, ftlunit);
  AIM_STATUS(aimInfo, status);

  // VREF (Real)
  status = card_addDouble(&card,vref);
  AIM_STATUS(aimInfo, status);

  // FLUTTF (Integer >= 0), can be blank
  status = card_addInteger(&card, fluttf);
  AIM_STATUS(aimInfo, status);

  // PRINT (Integer)
  status = card_addInteger(&card, print);
  AIM_STATUS(aimInfo, status);

  // Vi (Real > 0.0) // TODO: verify ascending order
  status = card_addDoubleArray(&card, numVelocities, velocities);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write TABDMP1 card
int zaero_card_tabdmp1(void *aimInfo,
                       FILE *fp, int tid, char *type,
                       int numDamping,
                       double *dampingFreq,
                       double *dampingValues,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;
  int i;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "TABDMP1", formatType);
  AIM_STATUS(aimInfo, status);

  // TID
  status = card_addInteger(&card, tid);
  AIM_STATUS(aimInfo, status);

  // TYPE, can be blank
  status = card_addString(&card, type);
  AIM_STATUS(aimInfo, status);

  // 6 empty fields
  for (i = 0; i < 6; i++) {
    card_addBlank(&card);
  }

  // fi, gi
  for (i = 0; i < numDamping; i++) {
    status = card_addDouble(&card, dampingFreq[i]);
    AIM_STATUS(aimInfo, status);

    status = card_addDouble(&card, dampingValues[i]);
    AIM_STATUS(aimInfo, status);
  }

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write CAERO7 card
int zaero_card_caero7(void *aimInfo,
                      FILE *fp, int wid, char *label, int acoord,
                      int nspan, int nchord, int lspan, int ztaic, int pafoil7,
                      double rl[3], double rch, int lrchd, int attchr, int acordr,
                      double tl[3], double tch, int ltchd, int attcht, int acordt,
                      feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "CAERO7", formatType);
  AIM_STATUS(aimInfo, status);

  // WID (Integer > 0)
  status = card_addInteger(&card, wid);
  AIM_STATUS(aimInfo, status);

  // LABEL (Character)
  status = card_addString(&card, label);
  AIM_STATUS(aimInfo, status);

  // ACOORD (Integer >= 0 or Blank, default = 0)
  status = card_addInteger(&card, acoord);
  AIM_STATUS(aimInfo, status);

  // NSPAN (Integer >= 2)
  status = card_addInteger(&card, nspan);
  AIM_STATUS(aimInfo, status);

  // NCHORD (Integer >= 2)
  status = card_addInteger(&card, nchord);
  AIM_STATUS(aimInfo, status);

  // LSPAN (Integer >= 0)
  status = card_addInteger(&card, lspan);
  AIM_STATUS(aimInfo, status);

  // ZTAIC (Integer >= 0)
  status = card_addInteger(&card, ztaic);
  AIM_STATUS(aimInfo, status);

  // PAFOIL7 (Integer >= 0)
  status = card_addInteger(&card, pafoil7);
  AIM_STATUS(aimInfo, status);

  // XRL, YRL, ZRL (Real)
  status = card_addDouble(&card, rl[0]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, rl[1]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, rl[2]);
  AIM_STATUS(aimInfo, status);

  // RCH (Real)
  status = card_addDouble(&card, rch);
  AIM_STATUS(aimInfo, status);

  // LRCHD (Integer >= 0)
  status = card_addInteger(&card, lrchd);
  AIM_STATUS(aimInfo, status);

  // ATTCHR (Integer >= 0)
  status = card_addInteger(&card, attchr);
  AIM_STATUS(aimInfo, status);

  // ACORDR (Integer > 0 or Blank)
  status = card_addInteger(&card, acordr);
  AIM_STATUS(aimInfo, status);

  // <empty>
  status = card_addBlank(&card);
  AIM_STATUS(aimInfo, status);

  // XTL, YTL, ZTL (Real)
  status = card_addDouble(&card, tl[0]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, tl[1]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, tl[2]);
  AIM_STATUS(aimInfo, status);

  // TCH (Real)
  status = card_addDouble(&card, tch);
  AIM_STATUS(aimInfo, status);

  // LTCHD (Integer >= 0)
  status = card_addInteger(&card, ltchd);
  AIM_STATUS(aimInfo, status);

  // ATTCHT (Integer >= 0)
  status = card_addInteger(&card, attcht);
  AIM_STATUS(aimInfo, status);

  // ACORDT (Integer > 0 or Blank)
  status = card_addInteger(&card, acordt);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write SPLINE1 card
int zaero_card_spline1(void *aimInfo,
                       FILE *fp, int eid, char *model, int cp, int setk, int setg, double dz, double eps,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "SPLINE1", formatType);
  AIM_STATUS(aimInfo, status);

  // EID
  status = card_addInteger(&card, eid);
  AIM_STATUS(aimInfo, status);

  // MODEL (not used)
  status = card_addString(&card, model);
  AIM_STATUS(aimInfo, status);

  // CP
  status = card_addInteger(&card, cp);
  AIM_STATUS(aimInfo, status);

  // SETK
  status = card_addInteger(&card, setk);
  AIM_STATUS(aimInfo, status);

  // SETG
  status = card_addInteger(&card, setg);
  AIM_STATUS(aimInfo, status);

  // DZ
  status = card_addDouble(&card, dz);
  AIM_STATUS(aimInfo, status);

  // EPS
  status = card_addDouble(&card, eps);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write SPLINE3 card
int zaero_card_spline3(void *aimInfo,
                       FILE *fp, int eid, char *model,
                       int setk, int setg, double eps,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "SPLINE3", formatType);
  AIM_STATUS(aimInfo, status);

  // EID
  status = card_addInteger(&card, eid);
  AIM_STATUS(aimInfo, status);

  // MODEL (not used)
  status = card_addString(&card, model);
  AIM_STATUS(aimInfo, status);

  // CP (not used)
  status = card_addBlank(&card);
  AIM_STATUS(aimInfo, status);

  // SETK
  status = card_addInteger(&card, setk);
  AIM_STATUS(aimInfo, status);

  // SETG
  status = card_addInteger(&card, setg);
  AIM_STATUS(aimInfo, status);

  // DZ (not used)
  status = card_addBlank(&card);
  AIM_STATUS(aimInfo, status);

  // EPS
  status = card_addDouble(&card, eps);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PANLST2 card
int zaero_card_panlst2(void *aimInfo,
                       FILE *fp, int setid, int macroid, int boxBegin, int boxEnd,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PANLST2", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // MACROID
  status = card_addInteger(&card, macroid);
  AIM_STATUS(aimInfo, status);

  // BOX 1
  status = card_addInteger(&card, boxBegin);
  AIM_STATUS(aimInfo, status);

  // THRU
  status = card_addField(&card, "THRU", 1);
  AIM_STATUS(aimInfo, status);

  // BOX 3
  status = card_addInteger(&card, boxEnd);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write SET1 card
int zaero_card_set1(void *aimInfo,
                    FILE *fp, int sid, int numIntegers, int integers[],
                    feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "SET1", formatType);
  AIM_STATUS(aimInfo, status);

  // SID
  status = card_addInteger(&card, sid);
  AIM_STATUS(aimInfo, status);

  // Gi
  status = card_addIntegerArray(&card, numIntegers, integers);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PLTAERO card
int zaero_card_pltaero(void *aimInfo,
                       FILE *fp, int setid, char *femgrid, int offset,
                       char *form, char *filenm, char *cell, char *vct,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PLTAERO", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer != 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // FEMGRID (Character)
  status = card_addString(&card, femgrid);
  AIM_STATUS(aimInfo, status);

  // OFFSET (Integer >= 0, or Blank)
  status = card_addInteger(&card, offset);
  AIM_STATUS(aimInfo, status);

  // FORM (Character, Default = "TECPLOT")
  status = card_addString(&card, form);
  AIM_STATUS(aimInfo, status);

  // FILENM (Character)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // CELL (Character; Default = "NO")
  status = card_addString(&card, cell);
  AIM_STATUS(aimInfo, status);

  // VCT (Character, Default = "NO")
  status = card_addString(&card, vct);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PLTFLUT card
int zaero_card_pltflut(void *aimInfo,
                       FILE *fp, int setid, int idflut, int mode, int ntime,
                       double maxdisp, char *form, char *filenm, char *aeronm,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PLTFLUT", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // IDFLUT (Integer > 0)
  status = card_addInteger(&card, idflut);
  AIM_STATUS(aimInfo, status);

  // MODE (Integer > 0)
  status = card_addInteger(&card, mode);
  AIM_STATUS(aimInfo, status);

  // NTIME (Integer > 0, Default = 1)
  status = card_addInteger(&card, ntime);
  AIM_STATUS(aimInfo, status);

  // MAXDISP (Real > 0.0, Default = 1.0)
  status = card_addDouble(&card, maxdisp);
  AIM_STATUS(aimInfo, status);

  // FORM (Character, Default = "TECPLOT")
  status = card_addString(&card, form);
  AIM_STATUS(aimInfo, status);

  // FILENM (Character)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // AERONM (Character, default = "AEROGEOM.PAT")
  status = card_addString(&card, aeronm);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PLTMODE card
int zaero_card_pltmode(void *aimInfo,
                       FILE *fp, int setid, char *sym, int mode,
                       double maxdisp, char *form, char *filenm, char *aeronm,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PLTMODE", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // SYM (Character)
  status = card_addString(&card, sym);
  AIM_STATUS(aimInfo, status);

  // MODE (Integer > 0)
  status = card_addInteger(&card, mode);
  AIM_STATUS(aimInfo, status);

  // TYPE (NOT USED)
  status = card_addBlank(&card);
  AIM_STATUS(aimInfo, status);

  // MAXDISP (Real > 0.0, Default = 1.0)
  status = card_addDouble(&card, maxdisp);
  AIM_STATUS(aimInfo, status);

  // FORM (Character, Default = "TECPLOT")
  status = card_addString(&card, form);
  AIM_STATUS(aimInfo, status);

  // FILENM (Character)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // AERONM (Character, default = "AEROGEOM.PAT")
  status = card_addString(&card, aeronm);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PLTTRIM card
int zaero_card_plttrim(void *aimInfo,
                       FILE *fp, int setid, int idtrim, char *flex, char *type,
                       char *form, char *filenm, double scale, char *aeronm,
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  int numTypeOptions = 6;
  char *typeOptions[6] = {"FORCE", "AERO", "INERTIAL",  // TODO: with current file format, "INERTIAL" will break
      "CP", "DEFORM", "ELASTIC"};

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PLTTRIM", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // IDTRIM (Integer > 0)
  status = card_addInteger(&card, idtrim);
  AIM_STATUS(aimInfo, status);

  // FLEX (Character, Default = "FLEX")
  status = card_addString(&card, flex);
  AIM_STATUS(aimInfo, status);

  // TYPE (Character)
  if (!string_isInArray(type, numTypeOptions, typeOptions)) {
    AIM_ERROR(aimInfo, "Unknown PLTTRIM TYPE field value: %s\n", type);
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  status = card_addString(&card, type);
  AIM_STATUS(aimInfo, status);

  // FORM (Character, Default = "TECPLOT")
  status = card_addString(&card, form);
  AIM_STATUS(aimInfo, status);

  // FILENM (Character)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // SCALE (Real > 0.0, Default = 1.0)
  status = card_addDouble(&card, scale);
  AIM_STATUS(aimInfo, status);

  // AERONM (Character, default = "AEROGEOM.PAT")
  status = card_addString(&card, aeronm);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write PLTVG card
int zaero_card_pltvg(void *aimInfo,
                     FILE *fp, int setid, int idflut, int nmode,
                     char *xaxis, char *form, char *filenm, double refrho,
                     feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "PLTVG", formatType);
  AIM_STATUS(aimInfo, status);

  // SETID (Integer > 0)
  status = card_addInteger(&card, setid);
  AIM_STATUS(aimInfo, status);

  // IDFLUT (Integer > 0)
  status = card_addInteger(&card, idflut);
  AIM_STATUS(aimInfo, status);

  // NMODE (Integer, Default = all modes)
  status = card_addInteger(&card, nmode);
  AIM_STATUS(aimInfo, status);

  // XAXIS (Character)
  status = card_addString(&card, xaxis);
  AIM_STATUS(aimInfo, status);

  // FORM (Character, Default = "TABLE")
  status = card_addString(&card, form);
  AIM_STATUS(aimInfo, status);

  // FILENM (Character)
  status = card_addString(&card, filenm);
  AIM_STATUS(aimInfo, status);

  // REFRHO (Real > 0.0, Default = 0.0023769 slug/ft^3)
  status = card_addDouble(&card, refrho);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return status;
}

// write TRIM card
int zaero_card_trim(void *aimInfo,
                    FILE *fp, int trimid, int idmk, double qinf,
                    int idobj, int idcons, double rho[3], double wtmass,
                    double weight, double I[6], char *trnacc,
                    char *nx, char *ny, char *nz, char *pdot, char *qdot,
                    char *rdot, int loadset, int numVars, int *varIDs,
                    char **vals,
                    feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;
  int i;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "TRIM", formatType);
  AIM_STATUS(aimInfo, status);

  // TRIMID (Integer > 0)
  status = card_addInteger(&card, trimid);
  AIM_STATUS(aimInfo, status);

  // IDMK (Integer > 0)
  status = card_addInteger(&card, idmk);
  AIM_STATUS(aimInfo, status);

  // QINF (Real > 0.0)
  status = card_addDouble(&card, qinf);
  AIM_STATUS(aimInfo, status);

  // IDOBJ (Integer >= 0)
  status = card_addInteger(&card, idobj);
  AIM_STATUS(aimInfo, status);

  // IDCONS (Integer >= 0)
  status = card_addInteger(&card, idcons);
  AIM_STATUS(aimInfo, status);

  // RHOX, RHOY, RHOZ (Real)
  status = card_addDouble(&card, rho[0]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, rho[1]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, rho[2]);
  AIM_STATUS(aimInfo, status);

  // WTMASS (Real > 0.0)
  status = card_addDouble(&card, wtmass);
  AIM_STATUS(aimInfo, status);

  // WEIGHT (Real > 0.0)
  status = card_addDouble(&card, weight);
  AIM_STATUS(aimInfo, status);

  // IXX, IXY, IYY, IXZ, IYZ, IZZ (Real)
  status = card_addDouble(&card, I[I11]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, I[I21]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, I[I22]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, I[I31]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, I[I32]);
  AIM_STATUS(aimInfo, status);
  status = card_addDouble(&card, I[I33]);
  AIM_STATUS(aimInfo, status);

  // TRNACC (Character)
  status = card_addString(&card, trnacc);
  AIM_STATUS(aimInfo, status);

  // NX, NY, NZ (Character or Real)
  status = card_addString(&card, nx);
  AIM_STATUS(aimInfo, status);
  status = card_addString(&card, ny);
  AIM_STATUS(aimInfo, status);
  status = card_addString(&card, nz);
  AIM_STATUS(aimInfo, status);

  // PDOT, QDOT, RDOT (Character or Real)
  status = card_addString(&card, pdot);
  AIM_STATUS(aimInfo, status);
  status = card_addString(&card, qdot);
  AIM_STATUS(aimInfo, status);
  status = card_addString(&card, rdot);
  AIM_STATUS(aimInfo, status);

  // LOADSET (Integer >= 0)
  status = card_addInteger(&card, loadset);
  AIM_STATUS(aimInfo, status);

  // IDVARi (Integer > 0)
  // VALi (Character or Real)
  for (i = 0; i < numVars; i++) {
    status = card_addInteger(&card, varIDs[i]);
    AIM_STATUS(aimInfo, status);

    status = card_addString(&card, vals[i]);
    AIM_STATUS(aimInfo, status);
  }

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return CAPS_SUCCESS;
}

// write TRIMVAR card
int zaero_card_trimvar(void *aimInfo,
                       FILE *fp, int idvar, char *label, double lower,
                       double upper, int trimlnk, char *dmi, char *sym,
                       int guessInitial, double initial, char *dc[6],
                       feaFileTypeEnum formatType)
{
  int status = CAPS_SUCCESS;

  cardStruct card;

  if (fp == NULL) return CAPS_IOERR;

  // begin card
  status = card_initiate(&card, "TRIMVAR", formatType);
  AIM_STATUS(aimInfo, status);

  // IDVAR (Integer > 0)
  status = card_addInteger(&card, idvar);
  AIM_STATUS(aimInfo, status);

  // LABEL (Character, cannot be Blank)
  status = card_addString(&card, label);
  AIM_STATUS(aimInfo, status);

  // LOWER (Real)
  status = card_addDouble(&card, lower);
  AIM_STATUS(aimInfo, status);

  // UPPER (Real)
  status = card_addDouble(&card, upper);
  AIM_STATUS(aimInfo, status);

  // TRIMLNK (Integer > 0)
  status = card_addIntegerOrBlank(&card, trimlnk > 0 ? &trimlnk : NULL);
  AIM_STATUS(aimInfo, status);

  // DMI (Character)
  status = card_addString(&card, dmi);
  AIM_STATUS(aimInfo, status);

  // DMI (Character)
  status = card_addString(&card, sym);
  AIM_STATUS(aimInfo, status);

  // INITIAL (Real)
  status = card_addDoubleOrBlank(&card, guessInitial == (int)true ? &initial : NULL);
  AIM_STATUS(aimInfo, status);

  // DCD (Character)
  status = card_addString(&card, dc[0]);
  AIM_STATUS(aimInfo, status);

  // DCY (Character)
  status = card_addString(&card, dc[1]);
  AIM_STATUS(aimInfo, status);

  // DCL (Character)
  status = card_addString(&card, dc[2]);
  AIM_STATUS(aimInfo, status);

  // DCR (Character)
  status = card_addString(&card, dc[3]);
  AIM_STATUS(aimInfo, status);

  // DCM (Character)
  status = card_addString(&card, dc[4]);
  AIM_STATUS(aimInfo, status);

  // DCN (Character)
  status = card_addString(&card, dc[5]);
  AIM_STATUS(aimInfo, status);

  // write card to file
  card_write(&card, fp);

cleanup:

  card_destroy(&card);

  return CAPS_SUCCESS;
}
