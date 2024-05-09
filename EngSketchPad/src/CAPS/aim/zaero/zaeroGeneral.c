#include <string.h>
#include <math.h>

#include "aimUtil.h"
#include "vlmTypes.h"
#include "jsonUtils.h"
#include "zaeroCards.h"
#include "zaeroGeneral.h"

static int _getCaero7Parameters(void *aimInfo,
                                const feaAeroStruct *feaAero, int compIndex,
                                int *wid, char **label, int *acoord,
                                int *nspan, int *nchord, int *lspan,
                                int *ztaic, int *pafoil7,
                                double rl[3], double *rch, int *lrchd,
                                int *attchr, int *acordr,
                                double tl[3], double *tch, int *ltchd,
                                int *attcht, int *acordt)
{
  int status = CAPS_SUCCESS;

  // indexing
  int sectionIndex;

  // compIndex
  int compIndexLen;

  if (feaAero == NULL) return CAPS_NULLVALUE;
  if (nspan   == NULL) return CAPS_NULLVALUE;

  *wid = feaAero->surfaceID;
  *acoord = feaAero->coordSystemID;

  // append compIndex to feaAero.name to ensure unique label
  /*@-nullpass@*/
  compIndexLen  = snprintf(NULL, 0, "%d", compIndex);
  compIndexLen += strlen(feaAero->name) + 1;
  AIM_ALLOC(*label, compIndexLen, char, aimInfo, status);
  /*@+nullpass@*/
  // append to front to ensure does not get trimmed off in card
  snprintf(*label, compIndexLen, "%d%s", compIndex, feaAero->name);

  // get number of spanwise divisions
  if (feaAero->vlmSurface.Sspace == 0) {
    *nspan = feaAero->vlmSurface.NspanTotal+1;
  } else {
    *nspan = 0;
  }

  // get number of chordwise divisions
  if (feaAero->vlmSurface.Cspace == 0.0) {
    *nchord = feaAero->vlmSurface.Nchord+1;
  } else {
    *nchord = 0;
  }

  if (*nspan == 0) {
    printf("Not implemented yet!\n");
    return CAPS_NOTIMPLEMENT;
  }

  if (*nchord == 0) {
    printf("Not implemented yet!\n");
    return CAPS_NOTIMPLEMENT;
  }

  // TODO: we're making these 0/blank for now... NOT permanent
  *lspan = 0;
  *ztaic = 0;
  *pafoil7 = 0;

  // get root chord variables
  sectionIndex = feaAero->vlmSurface.vlmSection[0].sectionIndex;

  rl[0] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[0];
  rl[1] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[1];
  rl[2] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[2];

  *rch = feaAero->vlmSurface.vlmSection[sectionIndex].chord;

  // TODO: we're making these 0 or blank for now... NOT permanent
  *lrchd = 0;
  *attchr = 0;
  *acordr = -1; // blank

  // get tip chord variables
  sectionIndex = feaAero->vlmSurface.vlmSection[1].sectionIndex;

  tl[0] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[0];
  tl[1] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[1];
  tl[2] = feaAero->vlmSurface.vlmSection[sectionIndex].xyzLE[2];

  *tch = feaAero->vlmSurface.vlmSection[sectionIndex].chord;

  // TODO: we're making these 0 or blank for now... NOT permanent
  *ltchd = 0;
  *attcht = 0;
  *acordt = -1; // blank

cleanup:
  return status;
}

// create (write cards for) aerodynamic wing component from feaAeroStruct
int zaero_aerodynamicWingComponent(void *aimInfo,
                                   FILE *fp,
                                   const feaAeroStruct *feaAero,
                                   const int compIndex,
                                   const feaFileTypeEnum formatType)
{
  // status return
  int status;

  // caero7 parameters
  char *label = NULL;
  int wid, acoord, nspan, nchord, lspan, ztaic, pafoil7;
  double rl[3], rch;
  int lrchd, attchr, acordr;
  double tl[3], tch;
  int ltchd, attcht, acordt;

  status = _getCaero7Parameters(aimInfo,
                                feaAero, compIndex,
                                &wid, &label, &acoord, &nspan,
                                &nchord, &lspan, &ztaic, &pafoil7,
                                rl, &rch, &lrchd, &attchr, &acordr,
                                tl, &tch, &ltchd, &attcht, &acordt);
  AIM_STATUS(aimInfo, status);

  status = zaero_card_caero7(aimInfo,
                             fp, wid, label, acoord, nspan,
                             nchord, lspan, ztaic, pafoil7,
                             rl, rch, lrchd, attchr, acordr,
                             tl, tch, ltchd, attcht, acordt,
                             formatType);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(label);
  return status;
}

// create (write cards for) spline method from a feaAeroStruct
int zaero_splineMethod(void *aimInfo,
                       FILE *fp,
                       const feaAeroStruct *feaAero,
                       const zaeroSplineStruct *spline,
                       const feaFileTypeEnum formatType)
{
  // status return
  int status;

  // set ids
  int panlst2ID, set1ID;

  int numSpanWise;

  if (fp == NULL) return CAPS_IOERR;
  if (feaAero == NULL) return CAPS_NULLVALUE;

  // PANLST2 id
  panlst2ID = feaAero->surfaceID + 1;

  // SET1 id
  set1ID = feaAero->surfaceID + 2;

  fprintf(fp, "$\n$ %s\n$\n", feaAero->name);

  if (spline->method == 0) {
    AIM_ERROR(aimInfo, "Zero Displacement Spline Method (SPLINE0) is not implemented yet.\n");
    return CAPS_NOTIMPLEMENT;

  } else if (spline->method == 1) {
    status = zaero_card_spline1(aimInfo,
                                fp,
                                feaAero->surfaceID, // eid
                                " ", // model (not used)
                                -1, // cp (blank), TODO: hardcoded
                                panlst2ID, // setk
                                set1ID, // setg
                                spline->attachFlex, // dz
                                spline->eps, // eps
                                formatType);
    AIM_STATUS(aimInfo, status);

  } else if (spline->method == 2) {
    AIM_ERROR(aimInfo, "Beam Spline Method (SPLINE2) is not implemented yet.\n");
    return CAPS_NOTIMPLEMENT;

  } else if (spline->method == 3) {
    status = zaero_card_spline3(aimInfo,
                                fp,
                                feaAero->surfaceID, // eid
                                " ", // model (not used)
                                panlst2ID, // setk
                                set1ID, // setg
                                spline->eps, // eps
                                formatType);
    AIM_STATUS(aimInfo, status);

  } else {
    AIM_ERROR(aimInfo, "Unknown Spline Method (%d).\n", spline->method);
    return CAPS_BADVALUE;
  }


  // determine id of last aero box in spline
  if (feaAero->vlmSurface.NspanTotal > 0)
    numSpanWise = feaAero->vlmSurface.NspanTotal;
  else if (feaAero->vlmSurface.NspanSection > 0)
    numSpanWise = (feaAero->vlmSurface.numSection-1)*feaAero->vlmSurface.NspanSection;
  else {
    AIM_ERROR(  aimInfo, "Only one of numSpanTotal and numSpanPerSection can be non-zero!\n");
    AIM_ADDLINE(aimInfo, "       numSpanTotal      = %d\n", feaAero->vlmSurface.NspanTotal);
    AIM_ADDLINE(aimInfo, "       numSpanPerSection = %d\n", feaAero->vlmSurface.NspanSection);
    return CAPS_BADVALUE;
  }

  status = zaero_card_panlst2(aimInfo,
                              fp,
                              panlst2ID, // setid
                              feaAero->surfaceID, // macroid
                              feaAero->surfaceID, // box1
                              feaAero->surfaceID + numSpanWise*feaAero->vlmSurface.Nchord - 1, // thru box3
                              formatType);
  AIM_STATUS(aimInfo, status);

  status = zaero_card_set1(aimInfo,
                           fp,
                           set1ID, // sid
                           feaAero->numGridID,
                           feaAero->gridIDSet, // Gi
                           formatType);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

// create (write cards for) modeling physical data
int zaero_modelPhysicalData(void *aimInfo,
                            FILE *fp,
                            const zaeroHFGStruct *hfg,
                            const cfdUnitsStruct *units,
                            const feaFileTypeEnum formatType)
{
  int status;

  status = zaero_card_aeroz(aimInfo,
                            fp,
                            0, // acsid TODO: hardcoded
                            hfg->XZSymmetry, // xzsym
                            hfg->flip, // flip
                            NULL, // fmmunit (not used anyway)
                            units->length != NULL ? units->length : "m", // fmlunit
                            hfg->refChord, // refc
                            hfg->refSpan, // refb
                            hfg->refArea, // refs
                            hfg->refCenter[0], // refx
                            hfg->refCenter[1], // refy
                            hfg->refCenter[2], // refz
                            formatType);
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}

// create (write cards for) unsteady aerodynamics generator from zaeroUAICStruct
int zaero_unsteadyAerodynamicsGenerator(void *aimInfo,
                                        FILE *fp,
                                        const zaeroUAICStruct *uaic,
                                        const feaFileTypeEnum formatType)
{
  int status;

  char *save = NULL, *saveSAVE = "SAVE", *saveACQUIRE = "ACQUIRE";
  int mkaerozID, trimfltID;

  mkaerozID = uaic->id;
  trimfltID = 0; // TODO: hardcoded

  if (uaic->saveFlag == 0) {
    save = saveSAVE;
  }
  else if (uaic->saveFlag == 1) {
    save = saveACQUIRE;
  }
  else {
    AIM_ERROR(aimInfo, "Unknown UAIC saveFlag: %d\n", uaic->saveFlag);
    return CAPS_BADVALUE;
  }

  status = zaero_card_mkaeroz(aimInfo,
                              fp,
                              mkaerozID, // idmk
                              uaic->machNumber, // mach
                              uaic->methodFlag, // method
                              trimfltID, // idflt
                              save,
                              uaic->aicFilename, // filenm
                              uaic->printFlag, // print
                              uaic->numReducedFreq,
                              uaic->reducedFreq, // freqi
                              formatType);
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}

