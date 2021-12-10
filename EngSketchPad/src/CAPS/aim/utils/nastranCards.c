// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>

#include "aimUtil.h"

#include "cardUtils.h" // Card writing utility
#include "nastranCards.h"


/*
 * Write AELIST card
 */
int nastranCard_aelist(FILE *fp, int *sid, int numE, int *e,
                       feaFileTypeEnum formatType){
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AELIST", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ei
    status = card_addIntegerArray(&card, numE, e);
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
int nastranCard_aero(FILE *fp, int *acsid, double *velocity, double *refc,
                     double *rhoref, int *symxz, int *symxy,
                     feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AERO", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACSID
    status = card_addInteger(&card, *acsid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // VELOCITY
    status = card_addDoubleOrBlank(&card, velocity);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFC
    status = card_addDouble(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RHOREF
    status = card_addDouble(&card, *rhoref);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXZ
    status = card_addIntegerOrBlank(&card, symxz);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXY
    status = card_addIntegerOrBlank(&card, symxy);
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
int nastranCard_aeros(FILE *fp, int *acsid, int *rcsid, double *refc,
                      double *refb, double *refs, int *symxz,
                      int *symxy, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AEROS", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ACSID
    status = card_addInteger(&card, *acsid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RCSID
    status = card_addInteger(&card, *rcsid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFC
    status = card_addDouble(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFB
    status = card_addDouble(&card, *refb);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFS
    status = card_addDouble(&card, *refs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXZ
    status = card_addIntegerOrBlank(&card, symxz);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SYMXY
    status = card_addIntegerOrBlank(&card, symxy);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write AESURF card
 */
int nastranCard_aesurf(FILE *fp, int *id, char *label, int *cid, int *alid,
                       double *eff, char *ldw, int *crefc, int *crefs, 
                       double *pllim, double *pulim, double *hmllim, double *hmulim,
                       int *tqllim, int *tqulim, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "AESURF", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addIntegerOrBlank(&card, id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABEL
    status = card_addString(&card, label);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID1
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ALID1
    status = card_addInteger(&card, *alid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID2 == blank
    card_addBlank(&card);

    // ALID2 == blank
    card_addBlank(&card);

    // EFF
    status = card_addDoubleOrBlank(&card, eff);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LDW
    status = card_addString(&card, ldw);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CREFC
    status = card_addIntegerOrBlank(&card, crefc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CREFS
    status = card_addIntegerOrBlank(&card, crefs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PLLIM, PULIM
    status = card_addDoubleOrBlank(&card, pllim);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDoubleOrBlank(&card, pulim);
    if (status != CAPS_SUCCESS) goto cleanup;

    // HMLLIM, HMULIM
    status = card_addDoubleOrBlank(&card, hmllim);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDoubleOrBlank(&card, hmulim);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TQLLIM, TQULIM
    status = card_addIntegerOrBlank(&card, tqllim);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addIntegerOrBlank(&card, tqulim);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CAERO1 card
 */
int nastranCard_caero1(FILE *fp, int *eid, int *pid, int *cp,
                       int *nspan, int *nchord, int *lspan, int *lchord,
                       int *igid, double xyz1[3], double xyz4[3],
                       double *x12, double *x43, 
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CAERO1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CP
    status = card_addIntegerOrBlank(&card, cp);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSPAN
    status = card_addIntegerOrBlank(&card, nspan);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NCHORD
    status = card_addIntegerOrBlank(&card, nchord);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LSPAN
    status = card_addIntegerOrBlank(&card, lspan);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LCHORD
    status = card_addIntegerOrBlank(&card, lchord);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IGID
    status = card_addInteger(&card, *igid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X1, Y1, Z1
    status = card_addDoubleArray(&card, 3, xyz1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X12
    status = card_addDouble(&card, *x12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X4, Y4, Z4
    status = card_addDoubleArray(&card, 3, xyz4);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X43
    status = card_addDouble(&card, *x43);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CBAR card
 */
int nastranCard_cbar(FILE *fp, int *eid, int *pid, int g[2], 
                     double x[3], int *g0, int *pa, int *pb,
                     double wa[3], double wb[3],
                     feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CBAR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GA, GB
    status = card_addIntegerArray(&card, 2, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // if Xi, then X1, X2, X3; else G0 <blank> <blank>
    if (x != NULL) {
        status = card_addDoubleArray(&card, 3, x);
    }
    else {
        status = card_addInteger(&card, *g0);
        card_addBlanks(&card, 2);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // PA
    status = card_addIntegerOrBlank(&card, pa);
    if (status != CAPS_SUCCESS) goto cleanup;


    // PB
    status = card_addIntegerOrBlank(&card, pb);
    if (status != CAPS_SUCCESS) goto cleanup;


    // W1A, W2A, W3A
    if (wa != NULL) {
        status = card_addDoubleArray(&card, 3, wa);
    }
    else {
        status = card_addBlanks(&card, 3);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // W1B, W2B, W3B
    if (wa != NULL) {
        status = card_addDoubleArray(&card, 3, wb);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CDAMP2 card
 */
int nastranCard_cdamp2(FILE *fp, int *eid, double *b, int *g1, int *g2,
                       int *c1, int *c2, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CDAMP2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // B
    status = card_addDouble(&card, *b);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G1
    status = card_addIntegerOrBlank(&card, g1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C1
    status = card_addIntegerOrBlank(&card, c1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G2
    status = card_addIntegerOrBlank(&card, g2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C2
    status = card_addIntegerOrBlank(&card, c2);
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
int nastranCard_celas2(FILE *fp, int *eid, double *k, int *g1, int *g2,
                       int *c1, int *c2, double *ge, double *s,
                       feaFileTypeEnum formatType) {
    
    int status;

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

    // G1
    status = card_addIntegerOrBlank(&card, g1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C1
    status = card_addIntegerOrBlank(&card, c1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G2
    status = card_addIntegerOrBlank(&card, g2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C2
    status = card_addIntegerOrBlank(&card, c2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GE
    status = card_addDoubleOrBlank(&card, ge);
    if (status != CAPS_SUCCESS) goto cleanup;

    // S
    status = card_addDouble(&card, *s);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CMASS2 card
 */
int nastranCard_cmass2(FILE *fp, int *eid, double *m, int *g1, int *g2,
                       int *c1, int *c2, feaFileTypeEnum formatType) {
    
    int status;

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

    // G1
    status = card_addIntegerOrBlank(&card, g1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C1
    status = card_addIntegerOrBlank(&card, c1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G2
    status = card_addIntegerOrBlank(&card, g2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C2
    status = card_addIntegerOrBlank(&card, c2);
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
int nastranCard_conm2(FILE *fp, int *eid, int *g, int *cid, double *m,
                      double x[3], double i[6],
                      feaFileTypeEnum formatType) {
    
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
    status = card_addIntegerOrBlank(&card, cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    status = card_addDouble(&card, *m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // X1, X2, X3
    if (x != NULL) {
        status = card_addDoubleArray(&card, 3, x);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    else {
        card_addBlanks(&card, 3);
    }

    // <blank>
    card_addBlank(&card);

    // I11, I21, I22, I31, I32, I33
    if (i != NULL) {
        status = card_addDoubleArray(&card, 6, i);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CORD2* card
 */
static int _cord2Card(const char *cardname, FILE *fp, int *cid, int *rid, 
                          double a[3], double b[3], double c[3],
                          feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, cardname, formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RID
    status = card_addIntegerOrBlank(&card, rid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ai
    status = card_addDoubleArray(&card, 3, a);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Bi
    status = card_addDoubleArray(&card, 3, b);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ci
    status = card_addDoubleArray(&card, 3, c);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CORD2C card
 */
int nastranCard_cord2c(FILE *fp, int *cid, int *rid, 
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType) {
    return _cord2Card("CORD2C", fp, cid, rid, a, b, c, formatType);
}

/*
 * Write CORD2R card
 */
int nastranCard_cord2r(FILE *fp, int *cid, int *rid, 
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType) {
    return _cord2Card("CORD2R", fp, cid, rid, a, b, c, formatType);
}

/*
 * Write CORD2S card
 */
int nastranCard_cord2s(FILE *fp, int *cid, int *rid, 
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType) {
    return _cord2Card("CORD2S", fp, cid, rid, a, b, c, formatType);
}

/*
 * Write CQUAD4 card
 */
int nastranCard_cquad4(FILE *fp, int *eid, int *pid, int g[4], 
                       double *theta, int *mcid, double *zoffs,
                       double t[4], feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CQUAD4", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 4, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // THETA or MCID
    if (theta != NULL) {
        status = card_addDouble(&card, *theta);
    }
    else if (mcid != NULL) {
        status = card_addInteger(&card, *mcid);
    }
    else {
        status = card_addBlank(&card);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFFS
    status = card_addDoubleOrBlank(&card, zoffs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ti
    if (t != NULL) {

        // <2 blanks> before T values
        card_addBlanks(&card, 2);

        status = card_addDoubleArray(&card, 4, t);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CQUAD8 card
 */
int nastranCard_cquad8(FILE *fp, int *eid, int *pid, int g[8], 
                       double *theta, int *mcid, double *zoffs,
                       double t[4], feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CQUAD8", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 8, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ti
    if (t != NULL) {
        status = card_addDoubleArray(&card, 4, t);
    }
    else {
        status = card_addBlanks(&card, 4);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // THETA or MCID
    if (theta != NULL) {
        status = card_addDouble(&card, *theta);
    }
    else if (mcid != NULL) {
        status = card_addInteger(&card, *mcid);
    }
    else {
        status = card_addBlank(&card);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFFS
    status = card_addDoubleOrBlank(&card, zoffs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CSHEAR card
 */
int nastranCard_cshear(FILE *fp, int *eid, int *pid, int g[4],
                       feaFileTypeEnum formatType) {
    
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
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 4, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CTRIA3 card
 */
int nastranCard_ctria3(FILE *fp, int *eid, int *pid, int g[3], 
                       double *theta, int *mcid, double *zoffs,
                       double t[3], feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CTRIA3", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 3, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // THETA or MCID
    if (theta != NULL) {
        status = card_addDouble(&card, *theta);
    }
    else if (mcid != NULL) {
        status = card_addInteger(&card, *mcid);
    }
    else {
        status = card_addBlank(&card);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFFS
    status = card_addDoubleOrBlank(&card, zoffs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // Ti
    if (t != NULL) {

        // <2 blanks> before T values
        card_addBlanks(&card, 2);

        status = card_addDoubleArray(&card, 3, t);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write CTRIA6 card
 */
int nastranCard_ctria6(FILE *fp, int *eid, int *pid, int g[6], 
                       double *theta, int *mcid, double *zoffs,
                       double t[3], feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "CTRIA6", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addIntegerOrBlank(&card, pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, 6, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // THETA or MCID
    if (theta != NULL) {
        status = card_addDouble(&card, *theta);
    }
    else if (mcid != NULL) {
        status = card_addInteger(&card, *mcid);
    }
    else {
        status = card_addBlank(&card);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ZOFFS
    status = card_addDoubleOrBlank(&card, zoffs);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Ti
    if (t != NULL) {
        status = card_addDoubleArray(&card, 3, t);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/* 
 * Write DCONADD card
 */
int nastranCard_dconadd(FILE *fp, int *dcid, int numDC, int *dc,
                        feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONADD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DCID
    status = card_addInteger(&card, *dcid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DCi
    status = card_addIntegerArray(&card, numDC, dc);

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DCONSTR card
 */
int nastranCard_dconstr(FILE *fp, int *dcid, int *rid,
                        double *lallow, double *uallow,
                        feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DCONSTR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DCID
    status = card_addInteger(&card, *dcid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RID
    status = card_addInteger(&card, *rid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LALLOW
    status = card_addDoubleOrBlank(&card, lallow);
    if (status != CAPS_SUCCESS) goto cleanup;

    // UALLOW
    status = card_addDoubleOrBlank(&card, uallow);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DDVAL card
 */
int nastranCard_ddval(FILE *fp, int *id, int numDV, double *dval, 
                      feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DDVAL", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DVALi
    status = card_addDoubleArray(&card, numDV, dval);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

// Split `string` into `numFields` 8-char fields. If shorter, remaining fields will be blank.
// TODO: add to cardUtils ?
/*@null@*/
static char ** _splitIntoFields(const char *string, const int numFields) {

    int i, j, numChars, fieldIndex, charIndex;
    char *field = NULL, **fields = NULL;

    numChars = strlen(string);

    fields = EG_alloc(sizeof(char *) * numFields);
    if (fields == NULL) {
        return NULL;
    }

    fieldIndex = 0;
    for (i = 0; i < numFields; i++) {

        field = EG_strdup("        "); // blank field
        if (field == NULL) {
            AIM_FREE(fields);
            return NULL; // TODO: we should clean up nicely in this event
        }

        // copy chars from string into field
        charIndex = i*8;
        for (j = 0; j < 8; j++) {
            if (charIndex + j < numChars) {
                field[j] = string[charIndex + j];
            }
        }

        fields[fieldIndex++] = field;
    }

    return fields;
}

static void _freeFields(int numFields, char **fields) {

    int i;

    if (fields != NULL) {
        for (i = 0; i < numFields; i++) {
            if (fields[i] != NULL) {
                EG_free(fields[i]);
            }
        }
    }
}

/*
 * Write DEQATN card
 */
int nastranCard_deqatn(FILE *fp, int *eqid, int numEquation,
                       char **equation) {

    
    int i, status;
    char **fields = NULL;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DEQATN", SmallField);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EQID
    status = card_addInteger(&card, *eqid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // break each equation string into 8-char fields
    if (equation != NULL && numEquation > 0) {

        // first line, max 56 chars
        fields = _splitIntoFields(equation[0], 7);
        status = card_addStringArray(&card, 7, fields);
        _freeFields(7, fields);
        AIM_FREE(fields);
        if (status != CAPS_SUCCESS) goto cleanup;

        // rest of lines, max 64 chars
        for (i = 1; i < numEquation; i++) {
            fields = _splitIntoFields(equation[i], 8);
            status = card_addStringArray(&card, 8, fields);
            _freeFields(8, fields);
            AIM_FREE(fields);
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
 * Write DESVAR card
 */
int nastranCard_desvar(FILE *fp, int *id, char *label, double *xinit, 
                       double *xlb, double *xub, double *delxv, int *ddval,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DESVAR", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABEL
    status = card_addString(&card, label);
    if (status != CAPS_SUCCESS) goto cleanup;

    // XINIT
    status = card_addDouble(&card, *xinit);
    if (status != CAPS_SUCCESS) goto cleanup;

    // XLB
    status = card_addDoubleOrBlank(&card, xlb);
    if (status != CAPS_SUCCESS) goto cleanup;

    // XUB
    status = card_addDoubleOrBlank(&card, xub);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DELXV
    status = card_addDoubleOrBlank(&card, delxv);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DDVAL
    status = card_addIntegerOrBlank(&card, ddval);
    if (status != CAPS_SUCCESS) goto cleanup;


    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DLINK card
 */
int nastranCard_dlink(FILE *fp, int *id, int *ddvid, double *c0, 
                      double *cmult, int numDV, int *idv, double *c,
                      feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DLINK", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DDVID
    status = card_addInteger(&card, *ddvid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C0
    status = card_addDoubleOrBlank(&card, c0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CMULT
    status = card_addDoubleOrBlank(&card, cmult);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IDV1, C1, IDV2, C2, -etc-
    for (i = 0; i < numDV; i++) {

        status = card_addInteger(&card, idv[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DMI card
 */
int nastranCard_dmi(FILE *fp, char *name, int *form, 
                    int *tin, int* tout, int m, int n, 
                    double *a, double *b,
                    feaFileTypeEnum formatType) {
    
    int i, j, k, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin header card
    status = card_initiate(&card, "DMI", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NAME
    status = card_addString(&card, name);
    if (status != CAPS_SUCCESS) goto cleanup;

    // J always "0" in header
    status = card_addInteger(&card, 0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // FORM
    status = card_addInteger(&card, *form);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TIN
    status = card_addInteger(&card, *tin);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TOUT
    status = card_addInteger(&card, *tout);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // M
    status = card_addInteger(&card, m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // N
    status = card_addInteger(&card, n);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write header card to file
    card_write(&card, fp);

    // columns
    for (i = 0; i < n; i++) {

        // begin column card
        status = card_initiate(&card, "DMI", formatType);
        if (status != CAPS_SUCCESS) goto cleanup;

        // NAME
        status = card_addString(&card, name);
        if (status != CAPS_SUCCESS) goto cleanup;

        // J
        status = card_addInteger(&card, i+1); // 1-bias
        if (status != CAPS_SUCCESS) goto cleanup;

        // I1
        status = card_addInteger(&card, 1);
        if (status != CAPS_SUCCESS) goto cleanup;

        // values
        for (j = 0; j < m; j++) {

            k = i*n + j;

            // A
            status = card_addDouble(&card, a[k]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // B
            if ((*tin > 2) && (b != NULL)) {
                status = card_addDouble(&card, b[k]);
                if (status != CAPS_SUCCESS) goto cleanup;
            }
        }


        // write column card to file
        card_write(&card, fp);
        card_destroy(&card);
    }

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DOPTPRM card
 */
int nastranCard_doptprm(FILE *fp, int numParam, char **param, 
                        int *paramType, void **val,
                        feaFileTypeEnum formatType){
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DOPTPRM", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PARAMi, VALi
    for (i = 0; i < numParam; i++) {

        status = card_addString(&card, param[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (paramType[i] == Double) {
            status = card_addDouble(&card, *((double *) val[i]));
            if (status != CAPS_SUCCESS) goto cleanup;
        }
        else if (paramType[i] == Integer) {
            status = card_addInteger(&card,  *((int *) val[i]));
            if (status != CAPS_SUCCESS) goto cleanup;
        }
        else {
            status = CAPS_BADVALUE;
            PRINT_ERROR("Unrecognized param type enum: %d", paramType[i]);
            goto cleanup;
        }
    }

    // write card to file
    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

/*
 * Write DRESP1 card
 */
int nastranCard_dresp1(FILE *fp, int *id, char *label, char *rtype,
                       char *ptype, int *region, int attaType, void *atta,    
                       int attbType, void *attb, 
                       int attsType, int numAtts, void *atts,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DRESP1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABEL
    status = card_addString(&card, label);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RTYPE
    status = card_addString(&card, rtype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PTYPE
    status = card_addString(&card, ptype);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REGION
    status = card_addIntegerOrBlank(&card, region);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ATTA
    if (atta == NULL) {
        card_addBlank(&card);
    }
    else if (attaType == Integer) {
        status = card_addIntegerOrBlank(&card, (int*) atta);
    }
    else if (attaType == Double) {
        status = card_addDoubleOrBlank(&card, (double*) atta);
    }
    else {
        status = CAPS_BADVALUE;
        PRINT_ERROR("DRESP1 attaType must be 1 (Integer) "
                    "or 2 (Double)");
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ATTB
    if (attb == NULL) {
        card_addBlank(&card);
    }
    else if (attbType == Integer) {
        status = card_addIntegerOrBlank(&card, (int*) attb);
    }
    else if (attbType == Double) {
        status = card_addDoubleOrBlank(&card, (double*) attb);
    }
    else {
        status = CAPS_BADVALUE;
        PRINT_ERROR("DRESP1 attbType must be 1 (Integer) "
                    "or 2 (Double)");
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // ATTi
    if (atts == NULL) {
        card_addBlank(&card);
    }
    else if (attsType == Integer) {
        status = card_addIntegerArray(&card, numAtts, (int*) atts);
    }
    else if (attsType == Double) {
        status = card_addDoubleArray(&card, numAtts, (double*) atts);
    }
    else {
        status = CAPS_BADVALUE;
        PRINT_ERROR("DRESP1 attsType must be 1 (Integer) "
                    "or 2 (Double)");
    }
    if (status != CAPS_SUCCESS) goto cleanup;


    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DRESP2 card
 */
int nastranCard_dresp2(FILE *fp, int *id, char *label, int *eqid,
                       int *region, int numDV, int *dvid, int numLabl,
                       char **labl, int numNR, int *nr, int numG, 
                       int *g, int *c, int numNRR, int *nrr,
                       feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DRESP2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABEL
    status = card_addString(&card, label);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EQID
    status = card_addInteger(&card, *eqid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REGION
    status = card_addIntegerOrBlank(&card, region);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <4 blanks>
    card_addBlanks(&card, 4);

    // DESVAR section
    if (dvid != NULL && numDV > 0) {

        status = card_addString(&card, "DESVAR");
        if (status != CAPS_SUCCESS) goto cleanup;

        // DVIDi
        for (i = 0; i < numDV; i++) {

            status = card_addInteger(&card, dvid[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // DESVAR section has weird formatting, seven values per line
            if ((card.size % 8) == 0)
                card_addBlank(&card);
        }
        // force continuation
        status = card_continue(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // DTABLE section
    if (labl != NULL && numLabl > 0) {

        status = card_addString(&card, "DTABLE");
        if (status != CAPS_SUCCESS) goto cleanup;

        // LABLj
        for (i = 0; i < numLabl; i++) {

            status = card_addString(&card, labl[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // DTABLE section has weird formatting, seven values per line
            if ((card.size % 8) == 0)
                card_addBlank(&card);
        }

        // force continuation
        status = card_continue(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // DRESP1 section
    if (nr != NULL && numNR > 0) {

        status = card_addString(&card, "DRESP1");
        if (status != CAPS_SUCCESS) goto cleanup;

        // NRk
        for (i = 0; i < numNR; i++) {

            status = card_addInteger(&card, nr[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // DRESP1 section has weird formatting, seven values per line
            if ((card.size % 8) == 0)
                card_addBlank(&card);
        }

        // force continuation
        status = card_continue(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // DNODE section
    if ((g != NULL) && (c != NULL) && numG > 0) {

        status = card_addString(&card, "DNODE");

        // Gi, Ci
        for (i = 0; i < numG; i++) {

            status = card_addInteger(&card, g[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = card_addInteger(&card, c[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // DNODE section has weird formatting, six values per line
            if ((card.size % 8) == 7)
                card_addBlanks(&card, 2);
        }

    }

    // DRESP2 section
    if (nrr != NULL && numNRR > 0) {

        status = card_addString(&card, "DRESP2");
        if (status != CAPS_SUCCESS) goto cleanup;

        // NRk
        for (i = 0; i < numNRR; i++) {

            status = card_addInteger(&card, nrr[i]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // DRESP2 section has weird formatting, seven values per line
            if ((card.size % 8) == 0)
                card_addBlank(&card);
        }

        // force continuation
        status = card_continue(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DTABLE card
 */
int nastranCard_dtable(FILE *fp, int numVal, char **labl, double *valu,
                       feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DTABLE", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABLi, VALUi
    for (i = 0; i < numVal; i++) {

        status = card_addString(&card, labl[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, valu[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DVCREL1 card
 */
int nastranCard_dvcrel1(FILE *fp, int *id, char *type, int *eid,
                        char *cpname, double *cpmin, double *cpmax,
                        double *c0, int numDV, int *dvid, double *coeff,
                        feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DVCREL1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TYPE
    status = card_addString(&card, type);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CPNAME
    status = card_addString(&card, cpname);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CPMIN
    status = card_addDoubleOrBlank(&card, cpmin);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CPMAX
    status = card_addDoubleOrBlank(&card, cpmax);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C0
    status = card_addDoubleOrBlank(&card, c0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // DVIDi, COEFFi
    for (i = 0; i < numDV; i++) {

        status = card_addInteger(&card, dvid[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, coeff[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DVMREL1 card
 */
int nastranCard_dvmrel1(FILE *fp, int *id, char *type, int *mid, 
                        char *mpname, double *mpmin, double *mpmax, 
                        double *c0, int numDV, int *dvid, double *coeff,
                        feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DVMREL1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TYPE
    status = card_addString(&card, type);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MPNAME
    status = card_addString(&card, mpname);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MPMIN
    status = card_addDoubleOrBlank(&card, mpmin);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MPMAX
    status = card_addDoubleOrBlank(&card, mpmax);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C0
    status = card_addDoubleOrBlank(&card, c0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // DVIDi, COEFFi
    for (i = 0; i < numDV; i++) {

        status = card_addInteger(&card, dvid[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, coeff[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write DVPREL1 card
 */
int nastranCard_dvprel1(FILE *fp, int *id, char *type, int *pid,
                        int *fid, char *pname, double *pmin, 
                        double *pmax, double *c0,
                        int numDV, int *dvid, double *coef,
                        feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "DVPREL1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TYPE
    status = card_addString(&card, type);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // FID/PNAME
    if (fid != NULL)
        status = card_addInteger(&card, *fid);
    else 
        status = card_addString(&card, pname);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PMIN
    status = card_addDoubleOrBlank(&card, pmin);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PMAX
    status = card_addDoubleOrBlank(&card, pmax);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C0
    status = card_addDoubleOrBlank(&card, c0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    for (i = 0; i < numDV; i++) {

        // DVIDi
        status = card_addInteger(&card, dvid[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        // COEFi
        status = card_addDouble(&card, coef[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write EIGR card
 */
int nastranCard_eigr(FILE *fp, int *sid, char *method, double *f1,
                     double *f2, int *ne, int *nd, char *norm,
                     int *g, int *c, feaFileTypeEnum formatType) {
    
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

    // F1
    status = card_addDoubleOrBlank(&card, f1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F2
    status = card_addDoubleOrBlank(&card, f2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NE
    status = card_addIntegerOrBlank(&card, ne);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ND
    status = card_addIntegerOrBlank(&card, nd);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <2 blanks>
    card_addBlanks(&card, 2);

    // NORM
    status = card_addString(&card, norm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addIntegerOrBlank(&card, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C
    status = card_addIntegerOrBlank(&card, c);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write EIGRL card
 */
int nastranCard_eigrl(FILE *fp, int *sid, double *v1, double *v2,
                      int *nd, int *msglvl, int *maxset,
                      double *shfscl, char *norm,
                      feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "EIGRL", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // V1
    status = card_addDouble(&card, *v1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // V2
    status = card_addDouble(&card, *v2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ND
    status = card_addIntegerOrBlank(&card, nd);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MSGLVL
    status = card_addIntegerOrBlank(&card, msglvl);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MAXSET
    status = card_addIntegerOrBlank(&card, maxset);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SHFSCL
    status = card_addDoubleOrBlank(&card, shfscl);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NORM
    status = card_addString(&card, norm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write FLFACT card
 */
int nastranCard_flfact(FILE *fp, int *sid, int numF, double *f,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "FLFACT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F
    status = card_addDoubleArray(&card, numF, f);
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
int nastranCard_flutter(FILE *fp, int *sid, char *method, int *dens,
                        int *mach, int *rfreq, char *imeth, 
                        int *nvalue, double *eps, 
                        feaFileTypeEnum formatType) {
    
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
    status = card_addInteger(&card, *mach);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RFREQ (or VEL)
    status = card_addInteger(&card, *rfreq);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IMETH
    status = card_addString(&card, imeth);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NVALUE
    status = card_addIntegerOrBlank(&card, nvalue);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EPS
    status = card_addDoubleOrBlank(&card, eps);
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
int nastranCard_force(FILE *fp, int *sid, int *g, int *cid,
                      double *f, double n[3], 
                      feaFileTypeEnum formatType) {
    
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
    status = card_addInteger(&card, *cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F
    status = card_addDouble(&card, *f);
    if (status != CAPS_SUCCESS) goto cleanup;

    // N
    status = card_addDoubleArray(&card, 3, n);
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
int nastranCard_grav(FILE *fp, int *sid, int *cid,
                     double *g, double n[3],
                     feaFileTypeEnum formatType) {
    
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
    status = card_addIntegerOrBlank(&card, cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addDouble(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // N
    status = card_addDoubleArray(&card, 3, n);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write LOAD card
 */
int nastranCard_load(FILE *fp, int *sid, double *s, int numL, 
                     double *ls, int *l, feaFileTypeEnum formatType)  {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "LOAD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // S
    status = card_addDouble(&card, *s);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (i = 0; i < numL; i++) {

        // Si
        status = card_addDouble(&card, ls[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Li
        status = card_addInteger(&card, l[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write MAT1 card
 */
int nastranCard_mat1(FILE *fp, int *mid, double *e, double* g,
                     double *nu, double *rho, double *a, double *tref,
                     double *ge, double *st, double *sc, double *ss,
                     int *mcsid, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MAT1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // E
    status = card_addDoubleOrBlank(&card, e);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addDoubleOrBlank(&card, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NU
    status = card_addDoubleOrBlank(&card, nu);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RHO
    status = card_addDouble(&card, *rho);
    if (status != CAPS_SUCCESS) goto cleanup;

    // A
    status = card_addDouble(&card, *a);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TREF
    status = card_addDoubleOrBlank(&card, tref);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GE
    status = card_addDouble(&card, *ge);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ST
    status = card_addDoubleOrBlank(&card, st);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SC
    status = card_addDoubleOrBlank(&card, sc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SS
    status = card_addDoubleOrBlank(&card, ss);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MCSID
    status = card_addIntegerOrBlank(&card, mcsid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write MAT8 card
 */
int nastranCard_mat8(FILE *fp, int *mid, double *e1, double *e2,
                     double *nu12, double *g12, double *g1z, 
                     double *g2z, double *rho, double *a1, double *a2,
                     double *tref, double *xt, double *xc, 
                     double *yt, double *yc, double *s, double *ge,
                     double *f12, double *strn,
                     feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MAT8", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // E1, E2
    status = card_addDouble(&card, *e1);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, *e2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NU12
    status = card_addDouble(&card, *nu12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G12
    status = card_addDouble(&card, *g12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G1Z
    status = card_addDoubleOrBlank(&card, g1z);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G2Z
    status = card_addDoubleOrBlank(&card, g2z);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RHO
    status = card_addDouble(&card, *rho);
    if (status != CAPS_SUCCESS) goto cleanup;

    // A1, A2
    status = card_addDouble(&card, *a1);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDouble(&card, *a2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TREF
    status = card_addDoubleOrBlank(&card, tref);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Xt, Xc
    status = card_addDoubleOrBlank(&card, xt);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDoubleOrBlank(&card, xc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Yt, Yc
    status = card_addDoubleOrBlank(&card, yt);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addDoubleOrBlank(&card, yc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // S
    status = card_addDoubleOrBlank(&card, s);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GE
    status = card_addDoubleOrBlank(&card, ge);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F12
    status = card_addDoubleOrBlank(&card, f12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // STRN
    status = card_addDoubleOrBlank(&card, strn);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write MKAERO1 card
 */
int nastranCard_mkaero1(FILE *fp, int numM, double *m, int numK,
                        double *k, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MKAERO1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    if (numM > 8) {
        status = CAPS_BADVALUE;
        PRINT_ERROR("Number of mach values must be less than 9");
        goto cleanup;
    }
    status = card_addDoubleArray(&card, numM, m);
    if (status != CAPS_SUCCESS) goto cleanup;
    status = card_addBlanks(&card, 8 - numM);
    if (status != CAPS_SUCCESS) goto cleanup;

    // K
    if (numK > 8) {
        status = CAPS_BADVALUE;
        PRINT_ERROR("Number of reduced freq values must be less than 9");
        goto cleanup;
    }
    status = card_addDoubleArray(&card, numK, k);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write MOMENT card
 */
int nastranCard_moment(FILE *fp, int *sid,
           /*@unused@*/int *g,
                       int *cid, double *m,
                       double n[3], feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "MOMENT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addIntegerOrBlank(&card, cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // M
    status = card_addDouble(&card, *m);
    if (status != CAPS_SUCCESS) goto cleanup;

    // N
    status = card_addDoubleArray(&card, 3, n);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write PAERO1 card
 */
int nastranCard_paero1(FILE *fp, int *pid, int numB, int *b,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PAERO1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Bi
    status = card_addIntegerArray(&card, numB, b);
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
int nastranCard_pbar(FILE *fp, int *pid, int *mid, double *a,
                     double *i1, double *i2, double *i12, double *j,
                     double *nsm, double c[2], double d[2],
                     double e[2], double f[2], double *k1, double *k2, 
                     feaFileTypeEnum formatType) {
    
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

    // I1
    status = card_addDoubleOrBlank(&card, i1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // I2
    status = card_addDoubleOrBlank(&card, i2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // J
    status = card_addDoubleOrBlank(&card, j);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    status = card_addDoubleOrBlank(&card, nsm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // C1, C2
    if (c != NULL) {
        status = card_addDoubleArray(&card, 2, c);
    }
    else {
        status = card_addBlanks(&card, 2);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // D1, D2
    if (d != NULL) {
        status = card_addDoubleArray(&card, 2, d);
    }
    else {
        status = card_addBlanks(&card, 2);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // E1, E2
    if (e != NULL) {
        status = card_addDoubleArray(&card, 2, e);
    }
    else {
        status = card_addBlanks(&card, 2);
    }
    if (status != CAPS_SUCCESS) goto cleanup;
    
    // F1, F2
    if (f != NULL) {
        status = card_addDoubleArray(&card, 2, f);
    }
    else {
        status = card_addBlanks(&card, 2);
    }
    if (status != CAPS_SUCCESS) goto cleanup;

    // K1
    status = card_addDoubleOrBlank(&card, k1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // K2
    status = card_addDoubleOrBlank(&card, k2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // I12
    status = card_addDoubleOrBlank(&card, i12);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write PBARL card
 */
int nastranCard_pbarl(FILE *fp, int *pid, int *mid, char *type, 
                      double *f0, int numDim, double *dim, double *nsm,
                      feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PBARL", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // TYPE
    status = card_addString(&card, type);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <3 blanks>
    card_addBlanks(&card, 3);

    // F0
    status = card_addDoubleOrBlank(&card, f0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // DIMi
    status = card_addDoubleArray(&card, numDim, dim);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    status = card_addDoubleOrBlank(&card, nsm);
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
int nastranCard_pcomp(FILE *fp, int *pid, double *z0, double* nsm,
                      double *sb, char *ft, double *tref, double *ge,
                      char *lam, int numLayers, int *mid, double *t,
                      double *theta, char **sout, 
                      feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PCOMP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Z0
    status = card_addDoubleOrBlank(&card, z0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    status = card_addDoubleOrBlank(&card, nsm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SB
    status = card_addDoubleOrBlank(&card, sb);
    if (status != CAPS_SUCCESS) goto cleanup;

    // FT
    status = card_addString(&card, ft);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TREF
    status = card_addDoubleOrBlank(&card, tref);
    if (status != CAPS_SUCCESS) goto cleanup;

    // GE
    status = card_addDoubleOrBlank(&card, ge);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LAM
    status = card_addString(&card, lam);
    if (status != CAPS_SUCCESS) goto cleanup;

    // for each layer
    for (i = 0; i < numLayers; i++) {

        // MIDi
        if (mid != NULL)
            status = card_addInteger(&card, mid[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Ti
        if (t != NULL)
            status = card_addDouble(&card, t[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // THi
        if (theta != NULL)
            status = card_addDouble(&card, theta[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;

        // SOUTi
        if (sout != NULL)
            status = card_addString(&card, sout[i]);
        else 
            card_addBlank(&card);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write PLOAD2 card
 */
int nastranCard_pload2(FILE *fp, int *sid, double *p, int numE,
                       int *eid, feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLOAD2", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // P
    status = card_addDouble(&card, *p);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EIDi
    status = card_addIntegerArray(&card, numE, eid);
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
int nastranCard_pload4(FILE *fp, int *sid, int *eid, double p[4],
                       int *g1, int *g3, int *cid, double n[3],
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PLOAD4", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // P
    status = card_addDoubleArray(&card, 4, p);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G1
    status = card_addIntegerOrBlank(&card, g1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G3 or G4
    status = card_addIntegerOrBlank(&card, g3);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addIntegerOrBlank(&card, cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // N1, N2, N3
    if (n != NULL) {
        status = card_addDoubleArray(&card, 3, n);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write PROD card
 */
int nastranCard_prod(FILE *fp, int *pid, int *mid, double *a, double *j,
                     double *c, double *nsm, 
                     feaFileTypeEnum formatType) {
    
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
    status = card_addDouble(&card, *a);
    if (status != CAPS_SUCCESS) goto cleanup;

    // J
    status = card_addDouble(&card, *j);
    if (status != CAPS_SUCCESS) goto cleanup;

    // C
    status = card_addDoubleOrBlank(&card, c);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    status = card_addDoubleOrBlank(&card, nsm);
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
int nastranCard_pshear(FILE *fp, int *pid, int *mid, double *t, 
                       double *nsm, double *f1, double *f2,
                       feaFileTypeEnum formatType) {
    
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
    status = card_addDoubleOrBlank(&card, nsm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F1
    status = card_addDoubleOrBlank(&card, f1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // F2
    status = card_addDoubleOrBlank(&card, f2);
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
int nastranCard_pshell(FILE *fp, int *pid, int *mid1, double *t,
                       int *mid2, double *i12t3, int *mid3, 
                       double *tst, double *nsm, double *z1, 
                       double *z2, int *mid4, 
                       feaFileTypeEnum formatType) {
    
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
    status = card_addIntegerOrBlank(&card, mid1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // T
    status = card_addDoubleOrBlank(&card, t);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID2
    status = card_addIntegerOrBlank(&card, mid2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // I12/T3
    status = card_addDoubleOrBlank(&card, i12t3);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID3
    status = card_addIntegerOrBlank(&card, mid3);
    if (status != CAPS_SUCCESS) goto cleanup;

    // TS/T
    status = card_addDoubleOrBlank(&card, tst);
    if (status != CAPS_SUCCESS) goto cleanup;

    // NSM
    status = card_addDoubleOrBlank(&card, nsm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Z1
    status = card_addDoubleOrBlank(&card, z1);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Z2
    status = card_addDoubleOrBlank(&card, z2);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID4
    status = card_addIntegerOrBlank(&card, mid4);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write PSOLID card
 */
int nastranCard_psolid(FILE *fp, int *pid, int *mid, int *cordm,
                       char *in, char *stress, char *isop, char *fctn, 
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "PSOLID", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // PID
    status = card_addInteger(&card, *pid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MID
    status = card_addInteger(&card, *mid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CORDM
    status = card_addIntegerOrBlank(&card, cordm);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IN
    status = card_addString(&card, in);
    if (status != CAPS_SUCCESS) goto cleanup;

    // STRESS
    status = card_addString(&card, stress);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ISOP
    status = card_addString(&card, isop);
    if (status != CAPS_SUCCESS) goto cleanup;

    // FCTN
    status = card_addString(&card, fctn);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write RBE2 card
 */
int nastranCard_rbe2(FILE *fp, int *eid, int *gn, int *cm, 
                     int numGM, int *gm, feaFileTypeEnum formatType) {   
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "RBE2", formatType);
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
    status = card_addIntegerArray(&card, numGM, gm);
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
int nastranCard_rbe3(FILE *fp, int *eid, int *refgrid, int *refc,
                     int numG, double *wt, int *c, int *g,
                     int numGM, int *gm, int *cm,
                     feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "RBE3", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // <blank>
    card_addBlank(&card);

    // REFGRID
    status = card_addInteger(&card, *refgrid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // REFC
    status = card_addInteger(&card, *refc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // WTi, Ci, Gi,j
    for (i = 0; i < numG; i++) {

        status = card_addDouble(&card, wt[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
        
        status = card_addInteger(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
        
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
        for (i = 0; i < numGM; i++) {

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
 * Write RFORCE card
 */
int nastranCard_rforce(FILE *fp, int *sid, int *g, int *cid, double *a,
                       double r[3], int *method, double *racc,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "RFORCE", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // G
    status = card_addInteger(&card, *g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CID
    status = card_addIntegerOrBlank(&card, cid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // A
    status = card_addDouble(&card, *a);
    if (status != CAPS_SUCCESS) goto cleanup;

    // R1, R2, R3
    status = card_addDoubleArray(&card, 3, r);
    if (status != CAPS_SUCCESS) goto cleanup;

    // METHOD
    status = card_addIntegerOrBlank(&card, method);
    if (status != CAPS_SUCCESS) goto cleanup;

    // RACC
    status = card_addDoubleOrBlank(&card, racc);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write SET1 card
 */
int nastranCard_set1(FILE *fp, int *sid, int numG, int *g,
                     feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SET1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi
    status = card_addIntegerArray(&card, numG, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write SPC card
 */
int nastranCard_spc(FILE *fp, int *sid, int numSPC, int *g, int *c,
                    double *d, feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPC", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Gi, Ci, Di
    for (i = 0; i < numSPC; i++) {

        status = card_addInteger(&card, g[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addInteger(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, d[i]);
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
int nastranCard_spc1(FILE *fp, int *sid, int *c, int numSPC, int *g,
                     feaFileTypeEnum formatType) {
    
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
    status = card_addIntegerArray(&card, numSPC, g);
    if (status != CAPS_SUCCESS) goto cleanup;

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write SPCADD card
 */
int nastranCard_spcadd(FILE *fp, int *sid, int numSPC, int *s,
                       feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPCADD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Si
    status = card_addIntegerArray(&card, numSPC, s);
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
int nastranCard_spline1(FILE *fp, int *eid, int *caero, int *box1, 
                        int *box2, int *setg, double *dz,
                        feaFileTypeEnum formatType) {
    
    int status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SPLINE1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // EID
    status = card_addInteger(&card, *eid);
    if (status != CAPS_SUCCESS) goto cleanup;

    // CAERO
    status = card_addInteger(&card, *caero);
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
    status = card_addDoubleOrBlank(&card, dz);
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
int nastranCard_suport(FILE *fp, int numID, int *id, int *c,
                       feaFileTypeEnum formatType){
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SUPORT", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // IDi, Ci
    for (i = 0; i < numID; i++) {

        status = card_addInteger(&card, id[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addInteger(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write SUPORT1 card
 */
int nastranCard_suport1(FILE *fp, int *sid, int numID, int *id, int *c,
                       feaFileTypeEnum formatType){
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "SUPORT1", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);

    // IDi, Ci
    for (i = 0; i < numID; i++) {

        status = card_addInteger(&card, id[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addInteger(&card, c[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
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
int nastranCard_temp(FILE *fp, int *sid, int numG, int *g, double *t,
                     feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TEMP", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SID
    status = card_addInteger(&card, *sid);

    // Gi, Ti
    for (i = 0; i < numG; i++) {

        status = card_addInteger(&card, g[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addInteger(&card, t[i]);
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
int nastranCard_tempd(FILE *fp, int numSID, int *sid, double *t,
                      feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TEMPD", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // SIDi, Ti
    for (i = 0; i < numSID; i++) {

        status = card_addInteger(&card, sid[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addInteger(&card, t[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

/*
 * Write TRIM card
 */
int nastranCard_trim(FILE *fp, int *id, double *mach, double *q,
                     int numVar, char **label, double *ux,
                     feaFileTypeEnum formatType) {
    
    int i, status;

    cardStruct card;

    if (fp == NULL) return CAPS_IOERR;

    // begin card
    status = card_initiate(&card, "TRIM", formatType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // ID
    status = card_addInteger(&card, *id);
    if (status != CAPS_SUCCESS) goto cleanup;

    // MACH
    status = card_addDouble(&card, *mach);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Q
    status = card_addDouble(&card, *q);
    if (status != CAPS_SUCCESS) goto cleanup;

    // LABELi, UXi
    for (i = 0; i < numVar; i++) {

        status = card_addString(&card, label[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = card_addDouble(&card, ux[i]);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // write card to file
    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}
