// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>

#include "cardUtils.h" // Card writing utility
#include "astrosCards.h"


/*
 * Write AEFACT card
 */
int astrosCard_aefact(FILE *fp, int *sid, int numD, double *D,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AEFACT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID (Unique Integer > 0)
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Di (Real)
    status = card_addDoubleArray(&card, numD, D);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write AERO card
 */
int astrosCard_aero(FILE *fp, int *acsid, double *refc, double *rhoref,
                     feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AERO", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACSID (Integer >= 0 or Blank)
    if (acsid != NULL) {
        status = card_addInteger(&card, *acsid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // REFC (Real >= 0)
    status = card_addDouble(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RHOREF (Real >= 0)
    status = card_addDouble(&card, *rhoref);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write AEROS card
 */
int astrosCard_aeros(FILE *fp, int *acsid, int *rcsid, double *refc,
                     double *refb, double *refs, /*@null@*/ int *gref,
                     /*@null@*/ double *refd, /*@null@*/ double *refl,
                     feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AEROS", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACSID (Integer > 0)
    status = card_addInteger(&card, *acsid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RCSID (Integer > 0, or blank)
    if (rcsid != NULL) {
        status = card_addInteger(&card, *rcsid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // REFC (Real > 0.0)
    status = card_addDouble(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFB (Real > 0.0)
    status = card_addDouble(&card, *refb);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFS (Real > 0.0)
    status = card_addDouble(&card, *refs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GREF (Integer > 0)
    if (gref != NULL)
        status = card_addInteger(&card, *gref);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFD (Real > 0.0, or blank)
    if (refd != NULL)
        status = card_addDouble(&card, *refd);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFL (Real > 0.0, or blank)
    if (refl != NULL)
        status = card_addDouble(&card, *refl);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write AIRFOIL card
 */
int astrosCard_airfoil(FILE *fp, int *acid, char *cmpnt, /*@null@*/ int *cp,
                       int *chord, int *usothk, int *lso, /*@null@*/ int *cam,
                       /*@null@*/ double *radius, double x1y1z1[3], double *x12,
                       /*@null@*/ int *ipanel, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AIRFOIL", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACID
    status = card_addInteger(&card, *acid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CMPNT
    status = card_addString(&card, cmpnt);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (cp != NULL)
        status = card_addInteger(&card, *cp);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CHORD
    status = card_addInteger(&card, *chord);
    if (status != CAPS_SUCCESS) goto cleanup;

    // USO/THK
    if (usothk != NULL)
        status = card_addInteger(&card, *usothk);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LSO
    if (lso != NULL)
        status = card_addInteger(&card, *lso);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CAM
    if (cam != NULL)
        status = card_addInteger(&card, *cam);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RADIUS
    if (radius != NULL)
        status = card_addDouble(&card, *radius);
    else
        card_addBlank(&card);;
    if (status != CAPS_SUCCESS) goto cleanup;

    // X1, Y1, Z1
    status = card_addDoubleArray(&card, 3, x1y1z1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X12
    if (x12 != NULL)
        status = card_addDouble(&card, *x12);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IPANEL
    if (ipanel != NULL)
        status = card_addInteger(&card, *ipanel);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CAERO6 card
 */
int astrosCard_caero6(FILE *fp, int *acid, char *cmpnt, /*@null@*/ int *cp,
                      int *igrp, int *lchord, int *lspan,
                      feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CAERO6", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACID
    status = card_addInteger(&card, *acid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CMPNT
    status = card_addString(&card, cmpnt);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CP
    if (cp != NULL)
        status = card_addInteger(&card, *cp);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IGRP
    status = card_addInteger(&card, *igrp);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LCHORD
    if (lchord != NULL)
        status = card_addInteger(&card, *lchord);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LSPAN
    if (lspan != NULL)
        status = card_addInteger(&card, *lspan);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CELAS2 card
 */
int astrosCard_celas2(FILE *fp, int *eid, double *k, int G[2], int C[2],
                      double *ge, /*@unused@*/ double *s,
                      /*@null@*/ double *tmin, /*@null@*/ double *tmax,
                      feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CELAS2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // K
    status = card_addDouble(&card, *k);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi, Ci
    for (i = 0; i < 2; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
        status = card_addInteger(&card, C[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // GE
    status = card_addDouble(&card, *ge);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN, TMAX
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (tmax != NULL)
        status = card_addDouble(&card, *tmax);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

int astrosCard_cihex1(FILE *fp, int *eid, int *pid, int G[8],
                      feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CIHEX1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID (Integer > 0)
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID (Integer > 0), defaults to EID
    if (pid == NULL) {
        card_addBlank(&card);
    }
    else {
        status = card_addInteger(&card, *pid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Gi (Integer > 0)
    // TODO: check Gi values unique?
    for (i = 0; i < 8; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CMASS2 card
 */
int astrosCard_cmass2(FILE *fp, int *eid, double *m, int G[2], int C[2],
                      /*@null@*/ double *tmin, /*@null@*/ double *tmax,
                       feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CMASS2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    status = card_addDouble(&card, *m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi, Ci
    for (i = 0; i < 2; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
        status = card_addInteger(&card, C[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // TMIN, TMAX
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (tmax != NULL)
        status = card_addDouble(&card, *tmax);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CONM2 card
 */
int astrosCard_conm2(FILE *fp, int *eid, int *g, int *cid, double *m,
                     double X[3], double I[6], /*@null@*/ double *tmin,
                     /*@null@*/ double *tmax, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CONM2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addInteger(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    status = card_addDouble(&card, *m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Xi
    if (X != NULL)
        status = card_addDoubleArray(&card, 3, X);
    else
        card_addBlanks(&card, 3);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // Iij
    if (I != NULL)
        status = card_addDoubleArray(&card, 6, I);
    else
        card_addBlanks(&card, 6);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN, TMAX
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (tmax != NULL)
        status = card_addDouble(&card, *tmax);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CQDMEM1 card
 */
int astrosCard_cqdmem1(FILE *fp, int *eid, int *pid, int G[4], int tmType,
                       /*@null@*/ void *tm, /*@null@*/ double *tmax,
                       feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CQDMEM1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID (Integer > 0)
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID (Integer > 0), defaults to EID
    if (pid == NULL) {
        card_addBlank(&card);
    }
    else {
        status = card_addInteger(&card, *pid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Gi (Integer > 0)
    for (i = 0; i < 4; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // TM (Real or blank; or 0 <= Integer <= 1,000,000)
    if (tm == NULL) {
        card_addBlank(&card);
    }
    else if (tmType == Integer) {
        status = card_addInteger(&card, *((int *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else if (tmType == Double) {
        status = card_addDouble(&card, *((double *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        PRINT_ERROR("CQDMEM1 tmType must be Integer or Real");
        status = CAPS_BADVALUE;
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMAX (Real > 0.0), may be ignored
    if (tmax != NULL) {
        status = card_addDouble(&card, *tmax);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CQUAD4 card
 */
int astrosCard_cquad4(FILE *fp, int *eid, int *pid, int G[4], int tmType,
                      /*@null@*/ void *tm, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmax, /*@null@*/ double *T,
                      feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CQUAD4", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID (Integer > 0)
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID (Integer > 0), defaults to EID
    if (pid == NULL) {
        card_addBlank(&card);
    }
    else {
        status = card_addInteger(&card, *pid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Gi (Integer > 0)
    for (i = 0; i < 4; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // TM (Real or blank; or 0 <= Integer <= 1,000,000)
    if (tm == NULL) {
        card_addBlank(&card);
    }
    else if (tmType == Integer) {
        status = card_addInteger(&card, *((int *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else if (tmType == Double) {
        status = card_addDouble(&card, *((double *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        PRINT_ERROR("CQUAD4 tmType must be Integer or Real");
        status = CAPS_BADVALUE;
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFF (Real or Blank)
    if (zoff != NULL) {
        status = card_addDouble(&card, *zoff);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    card_addBlank(&card); // 9th field is blank

    // TMAX (Real > 0.0), may be ignored
    if (tmax != NULL) {
        status = card_addDouble(&card, *tmax);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // Ti (Real or blank) 
    if (T != NULL) {
        for (i = 0; i < 4; i++) {
            status = card_addDouble(&card, T[i]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

int astrosCard_crod(FILE *fp, int *eid, int *pid, int G[2],
                    /*@null@*/ double *tmax, feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CROD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID (Integer > 0)
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID (Integer > 0), defaults to EID
    if (pid == NULL) {
        card_addBlank(&card);
    }
    else {
        status = card_addInteger(&card, *pid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Gi (Integer > 0)
    for (i = 0; i < 2; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // TMAX (Real > 0.0 or blank)
    if (tmax != NULL) {
        status = card_addDouble(&card, *tmax);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write CSHEAR card
 */
int astrosCard_cshear(FILE *fp, int *eid, int *pid, int G[4],
                      /*@null@*/ double *tmax, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CSHEAR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 4, G);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMAX
    if (tmax != NULL)
        status = card_addDouble(&card, *tmax);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

int astrosCard_ctria3(FILE *fp, int *eid, int *pid, int G[3], int tmType,
                      /*@null@*/ void *tm, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmax, /*@null@*/ double *T,
                      feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CTRIA3", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID (Integer > 0)
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID (Integer > 0), defaults to EID
    if (pid == NULL) {
        card_addBlank(&card);
    }
    else {
        status = card_addInteger(&card, *pid);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Gi (Integer > 0)
    for (i = 0; i < 3; i++) {
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // TM (Real or blank; or 0 <= Integer <= 1,000,000)
    if (tm == NULL) {
        card_addBlank(&card);
    }
    else if (tmType == Integer) {
        // TODO: do we really have to check < 1,000,000?
        status = card_addInteger(&card, *((int *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else if (tmType == Double) {
        status = card_addDouble(&card, *((double *) tm));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        PRINT_ERROR("CTRIA3 tmType must be Integer or Real");
        status = CAPS_BADVALUE;
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFF (Real or Blank)
    if (zoff != NULL) {
        status = card_addDouble(&card, *zoff);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    card_addBlanks(&card, 2); // 8th and 9th field are blank

    // TMAX (Real > 0.0), may be ignored
    if (tmax != NULL) {
        status = card_addDouble(&card, *tmax);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // Ti (Real or blank) 
    if (T != NULL) {
        for (i = 0; i < 4; i++) {
            status = card_addDouble(&card, T[i]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write DCONFLT
 */
int astrosCard_dconflt(FILE *fp, int *sid, char *vtype, double *gfact,
                       int numVal, double *v, double *gam, 
                       feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONFLT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VTYPE
    status = card_addString(&card, vtype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GFACT
    status = card_addDoubleOrBlank(&card, gfact);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numVal; i++) {

        // Vi
        status = card_addDouble(&card, v[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        // GAMi
        status = card_addDouble(&card, gam[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    } 

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/* 
 * Write DCONFLT card
 */
int astrosCard_dconfrq(FILE *fp, int *sid, int *mode, char *ctype, 
                       double *frqall, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONFRQ", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MODE
    status = card_addInteger(&card, *mode);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CTYPE
    status = card_addString(&card, ctype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // FRQALL
    status = card_addDouble(&card, *frqall);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write DCONTWP card
 */
int astrosCard_dcontwp(FILE *fp, int *sid, double *xt, double *xc,
                        double *yt, double *yc, double *ss, double *f12, 
                        char *ptype, int *layrnum, int numPID, int *PID,
                        feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONTWP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // XT
    status = card_addDouble(&card, *xt);
    if (status != CAPS_SUCCESS) goto cleanup;

    // XC
    if (xc != NULL)
        status = card_addDouble(&card, *xc);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // YT
    status = card_addDouble(&card, *yt);
    if (status != CAPS_SUCCESS) goto cleanup;

    // YC
    if (yc != NULL)
        status = card_addDouble(&card, *yc);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SS
    status = card_addDouble(&card, *ss);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F12
    status = card_addDouble(&card, *f12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PTYPE
    status = card_addString(&card, ptype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LAYRNUM
    if (layrnum != NULL)
        status = card_addInteger(&card, *layrnum);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PIDi
    status = card_addIntegerArray(&card, numPID, PID);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write DCONVMP card
 */
int astrosCard_dconvmp(FILE *fp, int *sid, double *st, /*@null@*/ double *sc,
                       /*@null@*/ double *ss, char *ptype,
                       /*@null@*/ int *layrnum, int numPID, int *PID,
                       feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONVMP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ST
    if (st != NULL)
        status = card_addDouble(&card, *st);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SC
    if (sc != NULL)
        status = card_addDouble(&card, *sc);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SS
    if (ss != NULL)
        status = card_addDouble(&card, *ss);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PTYPE
    status = card_addString(&card, ptype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LAYRNUM
    if (layrnum != NULL)
        status = card_addInteger(&card, *layrnum);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PIDi
    status = card_addIntegerArray(&card, numPID, PID);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write DESVARP card
 */
int astrosCard_desvarp(FILE *fp, int *dvid, int *linkid, double *vmin,
                       double *vmax, double *vinit, /*@null@*/ int *layrnum,
                       /*@null@*/ int *layrlst, /*@null@*/ char *label,
                       feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DESVARP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DVID
    status = card_addInteger(&card, *dvid);

    // LINKID
    if (linkid != NULL)
        status = card_addInteger(&card, *linkid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VMIN
    if (vmin != NULL)
        status = card_addDouble(&card, *vmin);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VMAX
    if (vmax != NULL)
        status = card_addDouble(&card, *vmax);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VINIT
    if (vinit != NULL)
        status = card_addDouble(&card, *vinit);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LAYRNUM
    if (layrnum != NULL)
        status = card_addInteger(&card, *layrnum);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LAYRLST
    if (layrlst != NULL)
        status = card_addInteger(&card, *layrlst);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

/*@-nullpass@*/
    status = card_addString(&card, label);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

// /*
//  * Write DVGRID card
//  */
// int astrosCard_(FILE *fp,
//                        feaFileTypeEnum formatType) {
    
//     int status;

//     cardStruct card;

//     if (fp == NULL) return CAPS_IOERR;

//     // begin card
//     status = card_initiate(&card, "DVGRID", formatType);
//     if (status != CAPS_SUCCESS) goto cleanup;

//     // write card to file
//     card_write(&card, fp);

//     cleanup:

//         card_destroy(&card);

//         return status;
// }

/*
 * Write EIGR card
 */
int astrosCard_eigr(FILE *fp, int *sid, char *method, double *f1,
                    double *f2, int *ne, int *nd, /*@null@*/ double *e,
                    char *norm, int *gid, int *dof, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "EIGR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // METHOD
    status = card_addString(&card, method);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F1, F2
    status = card_addDouble(&card, *f1);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, *f2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NE
    // v20: Not Used
    if (ne != NULL)
        status = card_addInteger(&card, *ne);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ND
    // v20: NVEC
    if (nd != NULL)
        status = card_addInteger(&card, *nd);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // E
    if (e != NULL)
        status = card_addDouble(&card, *e);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NORM
    status = card_addString(&card, norm);
    if (status != CAPS_SUCCESS) goto cleanup;

    if ((strcmp(norm, "POINT") == 0) || ((strcmp(norm, "point") == 0))) {

        // GID
        status = card_addInteger(&card, *gid);
        if (status != CAPS_SUCCESS) goto cleanup;

        // DOF
        status = card_addInteger(&card, *dof);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write FLFACT card
 */
int astrosCard_flfact(FILE *fp, int *sid, int numF, double *F,
                      feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "FLFACT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Fi
    status = card_addDoubleArray(&card, numF, F);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write FLUTTER card
 */
int astrosCard_flutter(FILE *fp, int *sid, char *method, int *dens,
                       double *mach, int *vel, /*@null@*/ int *mlist,
                       /*@null@*/ int *klist, /*@null@*/ int *effid,
                       int *symxz, int *symxy, /*@null@*/ double *eps,
                       /*@null@*/ char *curfit, /*@null@*/ int *nroot,
                       /*@null@*/ char *vtype, /*@null@*/ double *gflut,
                       /*@null@*/ double *gfilter, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "FLUTTER", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // METHOD
    status = card_addString(&card, method);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DENS
    status = card_addInteger(&card, *dens);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MACH
    status = card_addDouble(&card, *mach);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VEL
    status = card_addInteger(&card, *vel);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MLIST
    if (mlist != NULL)
        status = card_addInteger(&card, *mlist);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // KLIST
    if (klist != NULL)
        status = card_addInteger(&card, *klist);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EFFID
    if (effid != NULL)
        status = card_addInteger(&card, *effid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXZ, SYMXY
    status = card_addInteger(&card, *symxz);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addInteger(&card, *symxy);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EPS
    if (eps != NULL)
        status = card_addDouble(&card, *eps);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CURFIT
/*@-nullpass@*/
    status = card_addString(&card, curfit);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // NROOT
    if (nroot != NULL)
        status = card_addInteger(&card, *nroot);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VTYPE
/*@-nullpass@*/
    status = card_addString(&card, vtype);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // GFLUT
    if (gflut != NULL)
        status = card_addDouble(&card, *gflut);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GFILTER
    if (gfilter != NULL)
        status = card_addDouble(&card, *gfilter);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write FORCE card
 */
int astrosCard_force(FILE *fp, int *sid, int *g, int *cid, double *f,
                     double N[3], feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "FORCE", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addInteger(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    if (cid != NULL)
        status = card_addInteger(&card, *cid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F
    status = card_addDouble(&card, *f);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ni
    status = card_addDoubleArray(&card, 3, N);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write GRAV card
 */
int astrosCard_grav(FILE *fp, int *sid, int *cid, double *g, double N[3],
                    feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "GRAV", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addDouble(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ni
    status = card_addDoubleArray(&card, 3, N);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

int astrosCard_grid(FILE *fp, int *id, /*@null@*/ int *cp,
                     double xyz[3], /*@null@*/ int *cd, /*@null@*/ int *ps,
                     feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "GRID", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID (Integer > 0)
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CP (Integer > 0 or Blank)
    if (cp != NULL) {
        status = card_addInteger(&card, *cp);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // Xi (Real)
    status = card_addDouble(&card, xyz[0]);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, xyz[1]);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, xyz[2]);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CD (Integer > 0 or Blank)
    if (cd != NULL) {
        status = card_addInteger(&card, *cd);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // PS (Integer > 0 or Blank)
    if (ps != NULL) {
        status = card_addInteger(&card, *ps);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlank(&card);
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write MKAERO1 card
 */
int astrosCard_mkaero1(FILE *fp, int *symxz, int *symxy, int numM, double *M,
                       int numK, double *K, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MKAERO1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXZ
    status = card_addInteger(&card, *symxz);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXY
    status = card_addInteger(&card, *symxy);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Mi
    if (numM > 6) // at most 6 Mach numbers
        numM = 6;
    status = card_addDoubleArray(&card, numM, M);
    if (status != CAPS_SUCCESS) goto cleanup;

    card_addBlanks(&card, 6 - numM);

    // Ki
    if (numM > 8) // at most 8 reduced frequencies
        numM = 8;
    status = card_addDoubleArray(&card, numK, K);
    if (status != CAPS_SUCCESS) goto cleanup;

    card_addBlanks(&card, 8 - numM);

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write MOMENT card
 */
int astrosCard_moment(FILE *fp, int *sid, int *g, int *cid, double *m,
                      double N[3], feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MOMENT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addInteger(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    status = card_addDouble(&card, *m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ni
    status = card_addDoubleArray(&card, 3, N);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PBAR card
 */
int astrosCard_pbar(FILE *fp, int *pid, int *mid, double *a, double *i1,
                    double *i2, double *j, double *nsm, /*@null@*/ double *tmin,
                    /*@null@*/ double *k1, /*@null@*/double *k2,
                    /*@null@*/ double *C, /*@null@*/ double *D,
                    /*@null@*/ double *E, /*@null@*/ double *F,
                    /*@null@*/ double *r12, /*@null@*/ double *r22,
                    /*@null@*/ double *alpha, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PBAR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // A
    status = card_addDouble(&card, *a);
    if (status != CAPS_SUCCESS) goto cleanup;

    // I1, I2
    status = card_addDouble(&card, *i1);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, *i2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // J
    status = card_addDouble(&card, *j);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ci, Di, Ei, Fi
    if (C != NULL)
        status = card_addDoubleArray(&card, 2, C);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (D != NULL)
        status = card_addDoubleArray(&card, 2, D);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (E != NULL)
        status = card_addDoubleArray(&card, 2, E);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (F != NULL)
        status = card_addDoubleArray(&card, 2, F);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // K1, K2
    if (k1 != NULL)
        status = card_addDoubleArray(&card, 2, k1);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (k2 != NULL)
        status = card_addDoubleArray(&card, 2, k2);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // R12, R22, ALPHA
    if (r12 != NULL)
        status = card_addDoubleArray(&card, 2, r12);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (r22 != NULL)
        status = card_addDoubleArray(&card, 2, r22);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (alpha != NULL)
        status = card_addDoubleArray(&card, 2, alpha);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PCOMP card
 */
int astrosCard_pcomp(FILE *fp, int *pid, /*@null@*/ double *z0,
                     /*@null@*/ double *nsm, double *sbond, /*@null@*/ char *ft,
                     /*@null@*/ double *tmin, /*@null@*/ char *lopt, int numLayers,
                     int *MID, double *T, double *TH, /*@null@*/ char **SOUT,
                     int symmetricLaminate, feaFileTypeEnum formatType)
{
    int i, j, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PCOMP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Z0
    if (z0 != NULL)
        status = card_addDouble(&card, *z0);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SBOND
    status = card_addDouble(&card, *sbond);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F.T.
/*@-nullpass@*/
    status = card_addString(&card, ft);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // LOPT
/*@-nullpass@*/
    status = card_addString(&card, lopt);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // for each layer
    for (i = 0; i < numLayers; i++) {

        // MIDi
        if (MID != NULL)
            status = card_addInteger(&card, MID[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Ti
        if (T != NULL)
            status = card_addDouble(&card, T[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // THi
        if (TH != NULL)
            status = card_addDouble(&card, TH[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // SOUTi
        if (SOUT != NULL)
            status = card_addString(&card, SOUT[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // if symmetric laminate layers, loop over layers in reverse
    if (symmetricLaminate) {

        if (numLayers % 2 == 0 ) { // even
            j = numLayers - 1;
        } 
        else { // odd - don't repeat the last layer
            j = numLayers - 2;
        }

        for (i = j; i >= 0; i--) {

            // MIDi
            if (MID != NULL)
                status = card_addInteger(&card, MID[i]);
            else 
                card_addBlank(&card);
            if (status != CAPS_SUCCESS) goto cleanup;

            // Ti
            if (T != NULL)
                status = card_addDouble(&card, T[i]);
            else 
                card_addBlank(&card);
            if (status != CAPS_SUCCESS) goto cleanup;

            // THi
            if (TH != NULL)
                status = card_addDouble(&card, TH[i]);
            else 
                card_addBlank(&card);
            if (status != CAPS_SUCCESS) goto cleanup;

            // SOUTi
            if (SOUT != NULL)
                status = card_addString(&card, SOUT[i]);
            else 
                card_addBlank(&card);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PIHEX card
 */
int astrosCard_pihex(FILE *fp, int *pid, int *mid, /*@null@*/ int *cid,
                     /*@null@*/ int *nip, /*@null@*/ double *ar,
                     /*@null@*/ double *alpha, /*@null@*/ double *beta,
                     feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PIHEX", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    if (cid != NULL)
        status = card_addInteger(&card, *cid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NIP
    if (nip != NULL)
        status = card_addInteger(&card, *nip);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // AR
    if (ar != NULL)
        status = card_addDouble(&card, *ar);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ALPHA
    if (alpha != NULL)
        status = card_addDouble(&card, *alpha);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // BETA
    if (beta != NULL)
        status = card_addDouble(&card, *beta);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PLIST card
 */
int astrosCard_plist(FILE *fp, int *linkid, /*@null@*/ char *ptype, int numPID,
                     int *PID, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLIST", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LINKID
    status = card_addInteger(&card, *linkid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PTYPE
/*@-nullpass@*/
    status = card_addString(&card, ptype);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID1, PID2, PID3, ...
    status = card_addIntegerArray(&card, numPID, PID);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PLOAD card
 */
int astrosCard_pload(FILE *fp, int *sid, double *p, int numG, int *G,
                     feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLOAD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // P
    status = card_addDouble(&card, *p);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addIntegerArray(&card, numG, G);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PLOAD4 card
 */
int astrosCard_pload4(FILE *fp, int *lid, int *eid, int numP, double *P,
                      /*@null@*/ int *cid, /*@null@*/ double *V,
                      feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLOAD4", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LID
    status = card_addInteger(&card, *lid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Pi
    status = card_addDoubleArray(&card, numP, P);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    if (cid != NULL)
        status = card_addInteger(&card, *cid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Vi
    if (V != NULL)
        status = card_addDoubleArray(&card, 3, V);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PLYLIST card
 */
int astrosCard_plylist(FILE *fp, int *sid, int numP, int *P,
                       feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLYLIST", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Pi
    status = card_addIntegerArray(&card, numP, P);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PROD card
 */
int astrosCard_prod(FILE *fp, int *pid, int *mid, double *a, double *j,
                    double *c, double *nsm, /*@null@*/ double *tmin,
                    feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PROD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // A
    if (a != NULL)
        status = card_addDouble(&card, *a);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // J
    if (j != NULL)
        status = card_addDouble(&card, *j);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C
    if (c != NULL)
        status = card_addDouble(&card, *c);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PSHEAR card
 */
int astrosCard_pshear(FILE *fp, int *pid, int *mid, double *t, 
                      /*@null@*/ double *nsm, /*@null@*/ double *tmin,
                      feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PSHEAR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;
    
    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // T
    status = card_addDouble(&card, *t);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PQDMEM1 card
 */
int astrosCard_pqdmem1(FILE *fp, int *pid, int *mid, double *t, 
                       /*@null@*/ double *nsm, /*@null@*/ double *tmin,
                       feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PQDMEM1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;
    
    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // T
    status = card_addDouble(&card, *t);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write PSHELL card
 */
int astrosCard_pshell(FILE *fp, int *pid, int *mid1, double *t,
                      /*@null@*/ int *mid2, /*@null@*/ double *i12t3,
                      /*@null@*/ int *mid3, /*@null@*/ double *tst,
                      /*@null@*/ double *nsm, /*@null@*/ double *z1,
                      /*@null@*/double *z2, /*@null@*/ int *mid4, int mcsidType,
                      /*@null@*/ void *mcsid, int scsidType,
                      /*@null@*/ void *scsid, /*@null@*/ double *zoff,
                      /*@null@*/ double *tmin, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PSHELL", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID1
    if (mid1 != NULL)
        status = card_addInteger(&card, *mid1);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // T
    if (t != NULL)
        status = card_addDouble(&card, *t);
    else 
        card_addBlank(&card);

    // MID2
    if (mid2 != NULL)
        status = card_addInteger(&card, *mid2);
    else 
        card_addBlank(&card);

    // 12I/T3
    if (i12t3 != NULL)
        status = card_addDouble(&card, *i12t3);
    else 
        card_addBlank(&card);

    // MID3
    if (mid3 != NULL)
        status = card_addInteger(&card, *mid3);
    else 
        card_addBlank(&card);

    // TS/T
    if (tst != NULL)
        status = card_addDouble(&card, *tst);
    else 
        card_addBlank(&card);

    // NSM
    if (nsm != NULL)
        status = card_addDouble(&card, *nsm);
    else 
        card_addBlank(&card);

    // Z1, Z2
    if (z1 != NULL)
        status = card_addDouble(&card, *z1);
    else 
        card_addBlank(&card);

    if (z2 != NULL)
        status = card_addDouble(&card, *z2);
    else 
        card_addBlank(&card);

    // MID4
    if (mid4 != NULL)
        status = card_addInteger(&card, *mid4);
    else 
        card_addBlank(&card);

    // MCSID
    if (mcsid == NULL) {
        card_addBlank(&card);
    }
    else if (mcsidType == Integer) {
        status = card_addInteger(&card, *((int *) mcsid));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else if (mcsidType == Double) {
        status = card_addInteger(&card, *((double *) mcsid));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        status = CAPS_BADVALUE;
        PRINT_ERROR("mcsidType must be 1 (Integer) or 2(Double)");
        goto cleanup;
    }

    // SCSID
    if (scsid == NULL) {
        card_addBlank(&card);
    }
    else if (scsidType == Integer) {
        status = card_addInteger(&card, *((int *) scsid));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else if (scsidType == Double) {
        status = card_addInteger(&card, *((double *) scsid));
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        status = CAPS_BADVALUE;
        PRINT_ERROR("scsidType must be 1 (Integer) or 2 (Double)");
        goto cleanup;
    }

    // ZOFF
    if (zoff != NULL)
        status = card_addDouble(&card, *zoff);
    else 
        card_addBlank(&card);

    // TMIN
    if (tmin != NULL)
        status = card_addDouble(&card, *tmin);
    else 
        card_addBlank(&card);

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write RBE2 card
 */
int astrosCard_rbe2(FILE *fp, int *setid, int *eid, int *gn, int *cm,
                    int numGM, int *GM, feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "RBE2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SETID
    status = card_addInteger(&card, *setid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GN
    status = card_addInteger(&card, *gn);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CM
    status = card_addInteger(&card, *cm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GMi
    status = card_addIntegerArray(&card, numGM, GM);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write RBE3 card
 */
int astrosCard_rbe3(FILE *fp, int *setid, int *eid, int *refg, int *refc,
                    int numG, double *wt, int *c, int *g,
                    int numM, /*@null@*/ int *gm, /*@null@*/ int *cm,
                    feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "RBE3", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SETID
    status = card_addInteger(&card, *setid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFG
    status = card_addInteger(&card, *refg);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFC
    status = card_addInteger(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // WTi, Ci, Gi,j
    // These values must be contained within data 
    // fields 3 through 9
    for (i = 0; i < numG; i++) {

        if ((card.size % 8) == 0)
            card_addBlank(&card);

        status = card_addDouble(&card, wt[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if ((card.size % 8) == 0)
            card_addBlank(&card);
        
        status = card_addInteger(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if ((card.size % 8) == 0)
            card_addBlank(&card);
        
        status = card_addInteger(&card, g[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // if GM and CM defined
    if ((gm != NULL) && (cm != NULL)) {

        // force continuation
        status = card_continue(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // "UM"
        status = card_addString(&card, "UM");
        if (status != CAPS_SUCCESS) goto cleanup;

        // GMi, CMi
        for (i = 0; i < numM; i++) {

            // UM section has weird formatting, six values per line
            if ((card.size % 8) == 7)
                card_addBlanks(&card, 2);

            status = card_addInteger(&card, gm[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = card_addInteger(&card, cm[i]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write SPC card
 */
int astrosCard_spc(FILE *fp, int *sid, int numSPC, int *G, int *C, double *D,
                   feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPC", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numSPC; i++) {

        // G
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        // C
        if (C != NULL)
            status = card_addInteger(&card, C[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // D
        status = card_addDouble(&card, D[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write SPC1 card
 */
int astrosCard_spc1(FILE *fp, int *sid, int *c, int numSPC, int *G,
                     feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPC1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C
    status = card_addInteger(&card, *c);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, numSPC, G);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write SPLINE1 card
 */
int astrosCard_spline1(FILE *fp, int *eid, /*@null@*/ int *cp, int *macroid,
                       int *box1, int *box2, int *setg, /*@null@*/ double *dz,
                       feaFileTypeEnum formatType)
{
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPLINE1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CP
    if (cp != NULL)
        status = card_addInteger(&card, *cp);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MACROID
    status = card_addInteger(&card, *macroid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // BOX1
    status = card_addInteger(&card, *box1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // BOX2
    status = card_addInteger(&card, *box2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SETG
    status = card_addInteger(&card, *setg);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DZ
    if (dz != NULL)
        status = card_addDouble(&card, *dz);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write SUPORT card
 */
int astrosCard_suport(FILE *fp, int *setid, int numSUP, int *ID, int *C,
                      feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SUPORT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SETID
    status = card_addInteger(&card, *setid);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numSUP; i++) {

        // ID
        status = card_addInteger(&card, ID[i]);

        // C
        if (C != NULL)
            status = card_addInteger(&card, C[i]);
        else
            card_addBlank(&card);
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write TEMP card
 */
int astrosCard_temp(FILE *fp, int *sid, int numT, int *G, double *T,
                    feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TEMP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numT; i++) {

        // G
        status = card_addInteger(&card, G[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
        
        // T
        status = card_addDouble(&card, T[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write TEMPD card
 */
int astrosCard_tempd(FILE *fp, int numT, int *SID, double *T,
                     feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TEMPD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numT; i++) {

        // SID
        status = card_addInteger(&card, SID[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        // T
        status = card_addDouble(&card, T[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}
/*
 * TRIM functions
 */
int astrosCard_initiateTrimParam(astrosCardTrimParamStruct *param)
{

    if (param == NULL) return CAPS_NULLVALUE;

    param->label = NULL;
    param->value = 0.0;
    param->isFree = (int) false;

    return CAPS_SUCCESS;
}

int astrosCard_destroyTrimParam(astrosCardTrimParamStruct *param)
{

    if (param == NULL) return CAPS_NULLVALUE;

    if (param->label != NULL) EG_free(param->label);
    param->label = NULL;
    param->value = 0.0;
    param->isFree = (int) false;

    return CAPS_SUCCESS;
}

int astrosCard_trim(FILE *fp, int *trimid, double *mach, double *qdp,
                    /*@null@*/ char *trmtyp, /*@null@*/ int *effid,
                    /*@null@*/ double *vo, int numParam,
                    astrosCardTrimParamStruct *param, feaFileTypeEnum formatType)
{
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TRIM", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TRIMID
    status = card_addInteger(&card, *trimid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MACH
    if (mach != NULL)
        status = card_addDouble(&card, *mach);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // QDP
    status = card_addDouble(&card, *qdp);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TRMTYP
/*@-nullpass@*/
    status = card_addString(&card, trmtyp);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    // EFFID
    if (effid != NULL)
        status = card_addInteger(&card, *effid);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VO
    if (vo != NULL)
        status = card_addDouble(&card, *vo);
    else 
        card_addBlank(&card);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <two blanks>
    card_addBlanks(&card, 2);

    for (i = 0; i < numParam; i++) {

        // LABELi
        status = card_addString(&card, param[i].label);
        if (status != CAPS_SUCCESS) goto cleanup;

        // VALi
        if (param[i].isFree)
            status = card_addString(&card, "FREE");
        else
            status = card_addDouble(&card, param[i].value);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}
