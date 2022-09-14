#include "msesUtils.h"


void msesSensxFree(msesSensx **sensx)
{
  int       i, j, k;
  msesSensx *sStruct = *sensx;
  
  if (sStruct == NULL) return;
  
  if (sStruct->ileb   != NULL) EG_free(sStruct->ileb);
  if (sStruct->iteb   != NULL) EG_free(sStruct->iteb);
  if (sStruct->iend   != NULL) EG_free(sStruct->iend);
  if (sStruct->xleb   != NULL) EG_free(sStruct->xleb);
  if (sStruct->yleb   != NULL) EG_free(sStruct->yleb);
  if (sStruct->xteb   != NULL) EG_free(sStruct->xteb);
  if (sStruct->yteb   != NULL) EG_free(sStruct->yteb);
  if (sStruct->sblegn != NULL) EG_free(sStruct->sblegn);
  
  if (sStruct->xbi != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->xbi[i] != NULL) EG_free(sStruct->xbi[i]);
    EG_free(sStruct->xbi);
  }
  if (sStruct->ybi != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->ybi[i] != NULL) EG_free(sStruct->ybi[i]);
    EG_free(sStruct->ybi);
  }
  if (sStruct->cp != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->cp[i] != NULL) EG_free(sStruct->cp[i]);
    EG_free(sStruct->cp);
  }
  if (sStruct->hk != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->hk[i] != NULL) EG_free(sStruct->hk[i]);
    EG_free(sStruct->hk);
  }
  if (sStruct->cp_alfa != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->cp_alfa[i] != NULL) EG_free(sStruct->cp_alfa[i]);
    EG_free(sStruct->cp_alfa);
  }
  if (sStruct->hk_alfa != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->hk_alfa[i] != NULL) EG_free(sStruct->hk_alfa[i]);
    EG_free(sStruct->hk_alfa);
  }
  if (sStruct->cp_mach != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->cp_mach[i] != NULL) EG_free(sStruct->cp_mach[i]);
    EG_free(sStruct->cp_mach);
  }
  if (sStruct->hk_mach != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->hk_mach[i] != NULL) EG_free(sStruct->hk_mach[i]);
    EG_free(sStruct->hk_mach);
  }
  if (sStruct->cp_reyn != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->cp_reyn[i] != NULL) EG_free(sStruct->cp_reyn[i]);
    EG_free(sStruct->cp_reyn);
  }
  if (sStruct->hk_reyn != NULL) {
    for (i = 0; i < 2*sStruct->nbl; i++)
      if (sStruct->hk_reyn[i] != NULL) EG_free(sStruct->hk_reyn[i]);
    EG_free(sStruct->hk_reyn);
  }
  
  if (sStruct->modn    != NULL) EG_free(sStruct->modn);
  if (sStruct->al_mod  != NULL) EG_free(sStruct->al_mod);
  if (sStruct->cl_mod  != NULL) EG_free(sStruct->cl_mod);
  if (sStruct->cm_mod  != NULL) EG_free(sStruct->cm_mod);
  if (sStruct->cdw_mod != NULL) EG_free(sStruct->cdw_mod);
  if (sStruct->cdv_mod != NULL) EG_free(sStruct->cdv_mod);
  if (sStruct->cdf_mod != NULL) EG_free(sStruct->cdf_mod);
  
  j = 2*sStruct->nbl;
  if (sStruct->gn != NULL) {
    for (k = 0; k < sStruct->nmod; k++) {
      if (sStruct->gn[k] == NULL) continue;
      for (i = 0; i < j; i++)
        EG_free(sStruct->gn[k][i]);
      EG_free(sStruct->gn[k]);
    }
    EG_free(sStruct->gn);
  }
  if (sStruct->xbi_mod != NULL) {
    for (k = 0; k < sStruct->nmod; k++) {
      if (sStruct->xbi_mod[k] == NULL) continue;
      for (i = 0; i < j; i++)
        EG_free(sStruct->xbi_mod[k][i]);
      EG_free(sStruct->xbi_mod[k]);
    }
    EG_free(sStruct->xbi_mod);
  }
  if (sStruct->ybi_mod != NULL) {
    for (k = 0; k < sStruct->nmod; k++) {
      if (sStruct->ybi_mod[k] == NULL) continue;
      for (i = 0; i < j; i++)
        EG_free(sStruct->ybi_mod[k][i]);
      EG_free(sStruct->ybi_mod[k]);
    }
    EG_free(sStruct->ybi_mod);
  }
  if (sStruct->cp_mod != NULL) {
    for (k = 0; k < sStruct->nmod; k++) {
      if (sStruct->cp_mod[k] == NULL) continue;
      for (i = 0; i < j; i++)
        EG_free(sStruct->cp_mod[k][i]);
      EG_free(sStruct->cp_mod[k]);
    }
    EG_free(sStruct->cp_mod);
  }
  if (sStruct->hk_mod != NULL) {
    for (k = 0; k < sStruct->nmod; k++) {
      if (sStruct->hk_mod[k] == NULL) continue;
      for (i = 0; i < j; i++)
        EG_free(sStruct->hk_mod[k][i]);
      EG_free(sStruct->hk_mod[k]);
    }
    EG_free(sStruct->hk_mod);
  }

  if (sStruct->nposel  != NULL) EG_free(sStruct->nposel);
  if (sStruct->nbpos   != NULL) {
    for (i = 0; i < sStruct->npos; i++)
      if (sStruct->nbpos[i] != NULL) EG_free(sStruct->nbpos[i]);
    EG_free(sStruct->nbpos);
  }
  if (sStruct->xbpos    != NULL) {
    for (i = 0; i < sStruct->npos; i++)
      if (sStruct->xbpos[i] != NULL) EG_free(sStruct->xbpos[i]);
    EG_free(sStruct->xbpos);
  }
  if (sStruct->ybpos    != NULL) {
    for (i = 0; i < sStruct->npos; i++)
      if (sStruct->ybpos[i] != NULL) EG_free(sStruct->ybpos[i]);
    EG_free(sStruct->ybpos);
  }
  if (sStruct->abpos   != NULL) {
    for (i = 0; i < sStruct->npos; i++)
      if (sStruct->abpos[i] != NULL) EG_free(sStruct->abpos[i]);
    EG_free(sStruct->abpos);
  }
  if (sStruct->posn    != NULL) EG_free(sStruct->posn);
  if (sStruct->al_pos  != NULL) EG_free(sStruct->al_pos);
  if (sStruct->cl_pos  != NULL) EG_free(sStruct->cl_pos);
  if (sStruct->cm_pos  != NULL) EG_free(sStruct->cm_pos);
  if (sStruct->cdw_pos != NULL) EG_free(sStruct->cdw_pos);
  if (sStruct->cdv_pos != NULL) EG_free(sStruct->cdv_pos);
  if (sStruct->cdf_pos != NULL) EG_free(sStruct->cdf_pos);
  j = 2*sStruct->nbl*sStruct->npos;
  if (sStruct->xbi_pos != NULL) {
    for (i = 0; i < j; i++)
      if (sStruct->xbi_pos[i] != NULL) EG_free(sStruct->xbi_pos[i]);
    EG_free(sStruct->xbi_pos);
  }
  if (sStruct->ybi_pos != NULL) {
    for (i = 0; i < j; i++)
      if (sStruct->ybi_pos[i] != NULL) EG_free(sStruct->ybi_pos[i]);
    EG_free(sStruct->ybi_pos);
  }
  if (sStruct->cp_pos  != NULL) {
    for (i = 0; i < j; i++)
      if (sStruct->cp_pos[i] != NULL) EG_free(sStruct->cp_pos[i]);
    EG_free(sStruct->cp_pos);
  }
  if (sStruct->hk_pos  != NULL) {
    for (i = 0; i < j; i++)
      if (sStruct->hk_pos[i] != NULL) EG_free(sStruct->hk_pos[i]);
    EG_free(sStruct->hk_pos);
  }
  
  EG_free(sStruct);
  *sensx = NULL;
}


int msesSensxRead(void *aimInfo, const char *filename, msesSensx **sensx)
{
  int       status, i, j, k, m, n, begMarker, endMarker, ib, ioff, ix;
  FILE      *fp;
  msesSensx *sStruct = NULL;
  
  *sensx = NULL;
  fp = aim_fopen(aimInfo, filename, "rb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }
  AIM_ALLOC(sStruct, 1, msesSensx, aimInfo, status);
  sStruct->ileb    = NULL;
  sStruct->iteb    = NULL;
  sStruct->iend    = NULL;
  sStruct->xleb    = NULL;
  sStruct->yleb    = NULL;
  sStruct->xteb    = NULL;
  sStruct->yteb    = NULL;
  sStruct->sblegn  = NULL;
  sStruct->xbi     = NULL;
  sStruct->ybi     = NULL;
  sStruct->cp      = NULL;
  sStruct->hk      = NULL;
  sStruct->cp_alfa = NULL;
  sStruct->hk_alfa = NULL;
  sStruct->cp_mach = NULL;
  sStruct->hk_mach = NULL;
  sStruct->cp_reyn = NULL;
  sStruct->hk_reyn = NULL;
  sStruct->modn    = NULL;
  sStruct->al_mod  = NULL;
  sStruct->cl_mod  = NULL;
  sStruct->cm_mod  = NULL;
  sStruct->cdw_mod = NULL;
  sStruct->cdv_mod = NULL;
  sStruct->cdf_mod = NULL;
  sStruct->gn      = NULL;
  sStruct->xbi_mod = NULL;
  sStruct->ybi_mod = NULL;
  sStruct->cp_mod  = NULL;
  sStruct->hk_mod  = NULL;
  sStruct->nposel  = NULL;
  sStruct->nbpos   = NULL;
  sStruct->xbpos   = NULL;
  sStruct->ybpos   = NULL;
  sStruct->abpos   = NULL;
  sStruct->posn    = NULL;
  sStruct->al_pos  = NULL;
  sStruct->cl_pos  = NULL;
  sStruct->cm_pos  = NULL;
  sStruct->cdw_pos = NULL;
  sStruct->cdv_pos = NULL;
  sStruct->cdf_pos = NULL;
  sStruct->xbi_pos = NULL;
  sStruct->ybi_pos = NULL;
  sStruct->cp_pos  = NULL;
  sStruct->hk_pos  = NULL;
  
  status = CAPS_IOERR;
  n = fread(&begMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  n = fread(sStruct->code, sizeof(char), 32, fp);
  if (n != 32) goto cleanup;
  n = fread(&endMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  for (i = 31; i >= 0; i--)
    if (sStruct->code[i] != ' ') break;
  sStruct->code[i+1] = 0;
  
  n = fread(&begMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  n = fread(sStruct->name, sizeof(char), 32, fp);
  if (n != 32) goto cleanup;
  n = fread(&endMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  for (i = 31; i >= 0; i--)
    if (sStruct->name[i] != ' ') break;
  sStruct->name[i+1] = 0;
  
  n = fread(&begMarker,      sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->kalfa, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->kmach, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->kreyn, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,      sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,       sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->ldepma, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->ldepre, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,       sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,     sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->alfa, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->mach, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->reyn, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,     sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,     sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->dnrms, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->drrms, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->dvrms, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->dnmax, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->drmax, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->dvmax, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,     sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,     sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->ii,   sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->nbl,  sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->nmod, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->npos, sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,     sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  AIM_ALLOC(sStruct->ileb, sStruct->nbl, int, aimInfo, status);
  AIM_ALLOC(sStruct->iteb, sStruct->nbl, int, aimInfo, status);
  AIM_ALLOC(sStruct->iend, sStruct->nbl, int, aimInfo, status);

  n = fread(&begMarker,          sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  for (i = 0; i < sStruct->nbl; i++) {
    n = fread(&sStruct->ileb[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->iteb[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,          sizeof(int), 1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  AIM_ALLOC(sStruct->xleb, sStruct->nbl, double, aimInfo, status);
  AIM_ALLOC(sStruct->yleb, sStruct->nbl, double, aimInfo, status);
  AIM_ALLOC(sStruct->xteb, sStruct->nbl, double, aimInfo, status);
  AIM_ALLOC(sStruct->yteb, sStruct->nbl, double, aimInfo, status);
  AIM_ALLOC(sStruct->sblegn, sStruct->nbl, double, aimInfo, status);

  n = fread(&begMarker,            sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  for (i = 0; i < sStruct->nbl; i++) {
    n = fread(&sStruct->xleb[i],   sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->yleb[i],   sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->xteb[i],   sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->yteb[i],   sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->sblegn[i], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  if (sStruct->kalfa < 1) sStruct->kalfa = 1;
  if (sStruct->kmach < 1) sStruct->kmach = 1;
  if (sStruct->kreyn < 1) sStruct->kreyn = 1;
  
  /* shift indices so i=1 at LE  (prevents wasted points upstream of LE) */
  for (i = 0; i < sStruct->nbl; i++) {
    ioff = sStruct->ileb[i] - 1;
    sStruct->ileb[i] -= ioff;
    sStruct->iteb[i] -= ioff;
    sStruct->iend[i]  = sStruct->ii - ioff;
  }
  
  n = fread(&begMarker,         sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cl,       sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cm,       sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdw,      sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdv,      sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdf,      sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->al_alfa,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cl_alfa,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cm_alfa,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdw_alfa, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdv_alfa, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdf_alfa, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->al_mach,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cl_mach,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cm_mach,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdw_mach, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdv_mach, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdf_mach, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->al_reyn,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cl_reyn,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cm_reyn,  sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdw_reyn, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdv_reyn, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&sStruct->cdf_reyn, sizeof(double), 1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker,         sizeof(int),    1, fp);
  if (n != 1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  AIM_ALLOC(sStruct->xbi, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->xbi[i] = NULL;
  AIM_ALLOC(sStruct->ybi, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->ybi[i] = NULL;
  AIM_ALLOC(sStruct->cp, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->cp[i] = NULL;
  AIM_ALLOC(sStruct->hk, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->hk[i] = NULL;
  AIM_ALLOC(sStruct->cp_alfa, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->cp_alfa[i] = NULL;
  AIM_ALLOC(sStruct->hk_alfa, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->hk_alfa[i] = NULL;
  AIM_ALLOC(sStruct->cp_mach, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->cp_mach[i] = NULL;
  AIM_ALLOC(sStruct->hk_mach, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->hk_mach[i] = NULL;
  AIM_ALLOC(sStruct->cp_reyn, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->cp_reyn[i] = NULL;
  AIM_ALLOC(sStruct->hk_reyn, 2*sStruct->nbl, double*, aimInfo, status);
  for (i = 0; i < 2*sStruct->nbl; i++) sStruct->hk_reyn[i] = NULL;
  
  for (i = 0; i < sStruct->nbl; i++) {
    j = sStruct->iend[i] - sStruct->ileb[i] + 1;
    AIM_ALLOC(sStruct->xbi[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->xbi[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->ybi[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->ybi[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->cp[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->cp[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->hk[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->hk[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->cp_alfa[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->cp_alfa[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->hk_alfa[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->hk_alfa[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->cp_mach[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->cp_mach[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->hk_mach[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->hk_mach[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->cp_reyn[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->cp_reyn[2*i+1], j, double, aimInfo, status);

    AIM_ALLOC(sStruct->hk_reyn[2*i  ], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->hk_reyn[2*i+1], j, double, aimInfo, status);
  }

  /* 105 DO loop */
  for (i = 0; i < sStruct->nbl; i++) {
    j = sStruct->iend[i] - sStruct->ileb[i];
    for (k = 0; k < 2; k++) {
      n = fread(&begMarker,                    sizeof(int),    1, fp);
      if (n != 1) goto cleanup;
      for (m = 0; m < j; m++) {
        n = fread(&sStruct->xbi[2*i+k][m],     sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->ybi[2*i+k][m],     sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->cp[2*i+k][m],      sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->hk[2*i+k][m],      sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->cp_alfa[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->hk_alfa[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->cp_mach[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->hk_mach[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->cp_reyn[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
        n = fread(&sStruct->hk_reyn[2*i+k][m], sizeof(double), 1, fp);
        if (n != 1) goto cleanup;
      }
      n = fread(&endMarker,                    sizeof(int),    1, fp);
      if (n != 1) goto cleanup;
      if (begMarker != endMarker) goto cleanup;
    }
  }
  
  if (sStruct->nmod == 0) goto NoMOD;
  AIM_ALLOC(sStruct->modn   , sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->al_mod , sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->cl_mod , sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->cm_mod , sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdw_mod, sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdv_mod, sStruct->nmod, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdf_mod, sStruct->nmod, double, aimInfo, status);
  
  j = sStruct->nmod;
  AIM_ALLOC(sStruct->gn, j, double**, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->gn[i] = NULL;
  AIM_ALLOC(sStruct->xbi_mod, j, double**, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->xbi_mod[i] = NULL;
  AIM_ALLOC(sStruct->ybi_mod, j, double**, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->ybi_mod[i] = NULL;
  AIM_ALLOC(sStruct->cp_mod, j, double**, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->cp_mod[i] = NULL;
  AIM_ALLOC(sStruct->hk_mod, j, double**, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->hk_mod[i] = NULL;
  
  j = 2*sStruct->nbl;
  for (k = 0; k < sStruct->nmod; k++) {
    AIM_ALLOC(sStruct->gn[k], j, double*, aimInfo, status);
    for (i = 0; i < j; i++) sStruct->gn[k][i] = NULL;
    AIM_ALLOC(sStruct->xbi_mod[k], j, double*, aimInfo, status);
    for (i = 0; i < j; i++) sStruct->xbi_mod[k][i] = NULL;
    AIM_ALLOC(sStruct->ybi_mod[k], j, double*, aimInfo, status);
    for (i = 0; i < j; i++) sStruct->ybi_mod[k][i] = NULL;
    AIM_ALLOC(sStruct->cp_mod[k], j, double*, aimInfo, status);
    for (i = 0; i < j; i++) sStruct->cp_mod[k][i] = NULL;
    AIM_ALLOC(sStruct->hk_mod[k], j, double*, aimInfo, status);
    for (i = 0; i < j; i++) sStruct->hk_mod[k][i] = NULL;
  }

  for (k = 0; k < sStruct->nmod; k++) {
    for (ib = 0; ib < sStruct->nbl; ib++) {
      j = sStruct->iteb[ib] - sStruct->ileb[ib] + 1;
      AIM_ALLOC(sStruct->gn[k][2*ib  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->gn[k][2*ib+1], j, double, aimInfo, status);

      j = sStruct->iend[ib] - sStruct->ileb[ib];
      AIM_ALLOC(sStruct->xbi_mod[k][2*ib  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->xbi_mod[k][2*ib+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->ybi_mod[k][2*ib  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->ybi_mod[k][2*ib+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->cp_mod[k][2*ib  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->cp_mod[k][2*ib+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->hk_mod[k][2*ib  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->hk_mod[k][2*ib+1], j, double, aimInfo, status);
    }
  }
    
  /* 20 DO loop */
  for (k = 0; k < sStruct->nmod; k++) {
    n = fread(&begMarker,           sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->modn[k],    sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->al_mod[k],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cl_mod[k],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cm_mod[k],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdw_mod[k], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdv_mod[k], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdf_mod[k], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&endMarker,           sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    for (ib = 0; ib < sStruct->nbl; ib++) {
      for (ix = 0; ix < 2; ix++) {
        j = sStruct->iteb[ib] - sStruct->ileb[ib] + 1;
        n = fread(&begMarker,              sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        n = fread(sStruct->gn[k][2*ib+ix], sizeof(double), j, fp);
        if (n != j) goto cleanup;
        n = fread(&endMarker,              sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        if (begMarker != endMarker) goto cleanup;
      
        j = sStruct->iend[ib] - sStruct->ileb[ib];
        n = fread(&begMarker,                     sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        for (m = 0; m < j; m++) {
          n = fread(&sStruct->xbi_mod[k][2*ib+ix][m], sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->ybi_mod[k][2*ib+ix][m], sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->cp_mod[k][2*ib+ix][m],  sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->hk_mod[k][2*ib+ix][m],  sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
        }
        n = fread(&endMarker,                     sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        if (begMarker != endMarker) goto cleanup;
      }
    }
  }
  
NoMOD:
  if (sStruct->npos) goto NoPOS;
  AIM_ALLOC(sStruct->nposel, sStruct->npos, int, aimInfo, status);
  AIM_ALLOC(sStruct->nbpos, sStruct->npos, int*, aimInfo, status);
  for (i = 0; i < sStruct->npos; i++) sStruct->nbpos[i] = NULL;
  AIM_ALLOC(sStruct->xbpos, sStruct->npos, double*, aimInfo, status);
  for (i = 0; i < sStruct->npos; i++) sStruct->xbpos[i] = NULL;
  AIM_ALLOC(sStruct->ybpos, sStruct->npos, double*, aimInfo, status);
  for (i = 0; i < sStruct->npos; i++) sStruct->ybpos[i] = NULL;
  AIM_ALLOC(sStruct->abpos, sStruct->npos, double*, aimInfo, status);
  for (i = 0; i < sStruct->npos; i++) sStruct->abpos[i] = NULL;
  
  AIM_ALLOC(sStruct->posn   , sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->al_pos , sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->cl_pos , sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->cm_pos , sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdw_pos, sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdv_pos, sStruct->npos, double, aimInfo, status);
  AIM_ALLOC(sStruct->cdf_pos, sStruct->npos, double, aimInfo, status);
  
  j = 2*sStruct->nbl*sStruct->npos;
  AIM_ALLOC(sStruct->xbi_pos, j, double*, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->xbi_pos[i] = NULL;
  AIM_ALLOC(sStruct->ybi_pos, j, double*, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->ybi_pos[i] = NULL;
  AIM_ALLOC(sStruct->cp_pos, j, double*, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->cp_pos[i] = NULL;
  AIM_ALLOC(sStruct->hk_pos, j, double*, aimInfo, status);
  for (i = 0; i < j; i++) sStruct->hk_pos[i] = NULL;
  for (i = 0; i < sStruct->npos; i++) {
    for (ib = 0; ib < sStruct->nbl; ib++) {
      ix = i*sStruct->nbl + ib;
      j  = sStruct->iend[ib] - sStruct->ileb[ib] + 1;
      AIM_ALLOC(sStruct->xbi_pos[2*ix  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->xbi_pos[2*ix+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->ybi_pos[2*ix  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->ybi_pos[2*ix+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->cp_pos[2*ix  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->cp_pos[2*ix+1], j, double, aimInfo, status);

      AIM_ALLOC(sStruct->hk_pos[2*ix  ], j, double, aimInfo, status);
      AIM_ALLOC(sStruct->hk_pos[2*ix+1], j, double, aimInfo, status);
    }
  }
  
  /* 30 DO loop -- untested! */
  for (i = 0; i < sStruct->npos; i++) {
    n = fread(&begMarker,          sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->nposel[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&endMarker,          sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;

    j = sStruct->nposel[i];
    AIM_ALLOC(sStruct->nbpos[i], j, int, aimInfo, status);

    n = fread(&begMarker,        sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(sStruct->nbpos[i], sizeof(int), j, fp);
    if (n != j) goto cleanup;
    n = fread(&endMarker,        sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    AIM_ALLOC(sStruct->xbpos[i], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->ybpos[i], j, double, aimInfo, status);
    AIM_ALLOC(sStruct->abpos[i], j, double, aimInfo, status);

    n = fread(&begMarker,               sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    for (ib = 0; ib < j; ib++) {
      n = fread(&sStruct->xbpos[i][ib], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&sStruct->ybpos[i][ib], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&sStruct->abpos[i][ib], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,               sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,           sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->posn[i],    sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->al_pos[i],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cl_pos[i],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cm_pos[i],  sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdw_pos[i], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdv_pos[i], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&sStruct->cdf_pos[i], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&endMarker,           sizeof(int),    1, fp);
    if (n != 1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    for (ib = 0; ib < sStruct->nbl; ib++) {
      ix = i*sStruct->nbl + ib;
      for (k = 0; k < 2; k++) {
        j = sStruct->iend[ib] - sStruct->ileb[ib];
        n = fread(&begMarker,                     sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        for (m = 0; m < j; m++) {
          n = fread(&sStruct->xbi_pos[2*ix+k][m], sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->ybi_pos[2*ix+k][m], sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->cp_pos[2*ix+k][m],  sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
          n = fread(&sStruct->hk_pos[2*ix+k][m],  sizeof(double), 1, fp);
          if (n != 1) goto cleanup;
        }
        n = fread(&endMarker,                     sizeof(int),    1, fp);
        if (n != 1) goto cleanup;
        if (begMarker != endMarker) goto cleanup;
      }
    }
  }
  
NoPOS:
  fclose(fp);
  fp     = NULL;
  *sensx = sStruct;
  status = CAPS_SUCCESS;
  
cleanup:
  if (fp != NULL) fclose(fp);
  if (status != CAPS_SUCCESS) msesSensxFree(&sStruct);
  return status;
}


void msesMdatFree(msesMdat **mdat)
{
  int      j, jj, k, nbl, ns;
  msesMdat *mStruct = *mdat;
  
  if (mStruct == NULL) return;

  jj  = mStruct->istate[1];
  nbl = mStruct->istate[2];
  ns  = 2*nbl;
  
  if (mStruct->mfract != NULL) EG_free(mStruct->mfract);
  if (mStruct->x != NULL) {
    for (j = 0; j < jj; j++)
      if (mStruct->x[j] != NULL) EG_free(mStruct->x[j]);
    EG_free(mStruct->x);
  }
  if (mStruct->y != NULL) {
    for (j = 0; j < jj; j++)
      if (mStruct->y[j] != NULL) EG_free(mStruct->y[j]);
    EG_free(mStruct->y);
  }
  if (mStruct->r != NULL) {
    for (j = 0; j < jj; j++)
      if (mStruct->r[j] != NULL) EG_free(mStruct->r[j]);
    EG_free(mStruct->r);
  }
  if (mStruct->h != NULL) {
    for (j = 0; j < jj; j++)
      if (mStruct->h[j] != NULL) EG_free(mStruct->h[j]);
    EG_free(mStruct->h);
  }
  
  if (mStruct->xb != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->xb[j] != NULL) EG_free(mStruct->xb[j]);
    EG_free(mStruct->xb);
  }
  if (mStruct->yb != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->yb[j] != NULL) EG_free(mStruct->yb[j]);
    EG_free(mStruct->yb);
  }
  if (mStruct->xpb != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->xpb[j] != NULL) EG_free(mStruct->xpb[j]);
    EG_free(mStruct->xpb);
  }
  if (mStruct->ypb != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->ypb[j] != NULL) EG_free(mStruct->ypb[j]);
    EG_free(mStruct->ypb);
  }
  if (mStruct->sb != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->sb[j] != NULL) EG_free(mStruct->sb[j]);
    EG_free(mStruct->sb);
  }
  if (mStruct->sginl != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->sginl[j] != NULL) EG_free(mStruct->sginl[j]);
    EG_free(mStruct->sginl);
  }
  if (mStruct->sgout != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->sgout[j] != NULL) EG_free(mStruct->sgout[j]);
    EG_free(mStruct->sgout);
  }
  if (mStruct->xw != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->xw[j] != NULL) EG_free(mStruct->xw[j]);
    EG_free(mStruct->xw);
  }
  if (mStruct->yw != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->yw[j] != NULL) EG_free(mStruct->yw[j]);
    EG_free(mStruct->yw);
  }
  if (mStruct->wgap != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->wgap[j] != NULL) EG_free(mStruct->wgap[j]);
    EG_free(mStruct->wgap);
  }
  if (mStruct->vcen != NULL) {
    for (j = 0; j < nbl; j++)
      if (mStruct->vcen[j] != NULL) EG_free(mStruct->vcen[j]);
    EG_free(mStruct->vcen);
  }
  
  if (mStruct->sg != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->sg[j] != NULL) EG_free(mStruct->sg[j]);
    EG_free(mStruct->sg);
  }
  if (mStruct->disp != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->disp[j] != NULL) EG_free(mStruct->disp[j]);
    EG_free(mStruct->disp);
  }
  if (mStruct->pspec != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->pspec[j] != NULL) EG_free(mStruct->pspec[j]);
    EG_free(mStruct->pspec);
  }
  if (mStruct->thet != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->thet[j] != NULL) EG_free(mStruct->thet[j]);
    EG_free(mStruct->thet);
  }
  if (mStruct->dstr != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->dstr[j] != NULL) EG_free(mStruct->dstr[j]);
    EG_free(mStruct->dstr);
  }
  if (mStruct->uedg != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->uedg[j] != NULL) EG_free(mStruct->uedg[j]);
    EG_free(mStruct->uedg);
  }
  if (mStruct->ctau != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->ctau[j] != NULL) EG_free(mStruct->ctau[j]);
    EG_free(mStruct->ctau);
  }
  if (mStruct->entr != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->entr[j] != NULL) EG_free(mStruct->entr[j]);
    EG_free(mStruct->entr);
  }
  if (mStruct->tauw != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->tauw[j] != NULL) EG_free(mStruct->tauw[j]);
    EG_free(mStruct->tauw);
  }
  if (mStruct->dint != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->dint[j] != NULL) EG_free(mStruct->dint[j]);
    EG_free(mStruct->dint);
  }
  if (mStruct->tstr != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->tstr[j] != NULL) EG_free(mStruct->tstr[j]);
    EG_free(mStruct->tstr);
  }
  
  if (mStruct->freq  != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->freq[j] != NULL) EG_free(mStruct->freq[j]);
    EG_free(mStruct->freq);
  }
  if (mStruct->famp != NULL) {
    for (k = j = 0; j < ns; j++) k += (mStruct->nfreq[j]+1)*ns;
    for (j = 0; j < k; j++)
      if (mStruct->famp[j] != NULL) EG_free(mStruct->famp[j]);
    EG_free(mStruct->famp);
  }
  if (mStruct->alfr != NULL) {
    for (k = j = 0; j < ns; j++) k += (mStruct->nfreq[j]+1)*ns;
    for (j = 0; j < k; j++)
      if (mStruct->alfr[j] != NULL) EG_free(mStruct->alfr[j]);
    EG_free(mStruct->alfr);
  }
  if (mStruct->nfreq != NULL) EG_free(mStruct->nfreq);
  
  if (mStruct->snor  != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->snor[j] != NULL) EG_free(mStruct->snor[j]);
    EG_free(mStruct->snor);
  }
  if (mStruct->xnor  != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->xnor[j] != NULL) EG_free(mStruct->xnor[j]);
    EG_free(mStruct->xnor);
  }
  if (mStruct->xsnor != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->xsnor[j] != NULL) EG_free(mStruct->xsnor[j]);
    EG_free(mStruct->xsnor);
  }
  if (mStruct->ynor  != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->ynor[j] != NULL) EG_free(mStruct->ynor[j]);
    EG_free(mStruct->ynor);
  }
  if (mStruct->ysnor != NULL) {
    for (j = 0; j < ns; j++)
      if (mStruct->ysnor[j] != NULL) EG_free(mStruct->ysnor[j]);
    EG_free(mStruct->ysnor);
  }
  if (mStruct->knor != NULL) EG_free(mStruct->knor);
  
  if (mStruct->cl_mod  != NULL) EG_free(mStruct->cl_mod);
  if (mStruct->cm_mod  != NULL) EG_free(mStruct->cm_mod);
  if (mStruct->cdw_mod != NULL) EG_free(mStruct->cdw_mod);
  if (mStruct->cdv_mod != NULL) EG_free(mStruct->cdv_mod);
  if (mStruct->cdf_mod != NULL) EG_free(mStruct->cdf_mod);
  if (mStruct->modn    != NULL) EG_free(mStruct->modn);
  if (mStruct->dmspn   != NULL) EG_free(mStruct->dmspn);
  
  if (mStruct->cl_pos  != NULL) EG_free(mStruct->cl_pos);
  if (mStruct->cm_pos  != NULL) EG_free(mStruct->cm_pos);
  if (mStruct->cdw_pos != NULL) EG_free(mStruct->cdw_pos);
  if (mStruct->cdv_pos != NULL) EG_free(mStruct->cdv_pos);
  if (mStruct->cdf_pos != NULL) EG_free(mStruct->cdf_pos);
  if (mStruct->posn    != NULL) EG_free(mStruct->posn);
  if (mStruct->dpspn   != NULL) EG_free(mStruct->dpspn);
  
  for (j = 0; j < NBITX; j++)
    if (mStruct->isbits[j] != NULL) EG_free(mStruct->isbits[j]);
  
  EG_free(mStruct);
  *mdat = NULL;
}


int msesMdatRead(void *aimInfo, const char *filename, msesMdat **mdat)
{
  int      status, i, j, k, n, ii, jj, nbl, ns, nmodn, nposn, nsuct;
  int      begMarker, endMarker;
  FILE     *fp;
  msesMdat *mStruct = NULL;
  
  *mdat = NULL;
  fp    = aim_fopen(aimInfo, filename, "rb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }
  AIM_ALLOC(mStruct, 1, msesMdat, aimInfo, status);
  mStruct->mfract  = NULL;
  mStruct->x       = NULL;
  mStruct->y       = NULL;
  mStruct->r       = NULL;
  mStruct->h       = NULL;
  mStruct->xb      = NULL;
  mStruct->yb      = NULL;
  mStruct->xpb     = NULL;
  mStruct->ypb     = NULL;
  mStruct->sb      = NULL;
  mStruct->sginl   = NULL;
  mStruct->sgout   = NULL;
  mStruct->xw      = NULL;
  mStruct->yw      = NULL;
  mStruct->wgap    = NULL;
  mStruct->vcen    = NULL;
  mStruct->sg      = NULL;
  mStruct->disp    = NULL;
  mStruct->pspec   = NULL;
  mStruct->thet    = NULL;
  mStruct->dstr    = NULL;
  mStruct->uedg    = NULL;
  mStruct->ctau    = NULL;
  mStruct->entr    = NULL;
  mStruct->tauw    = NULL;
  mStruct->dint    = NULL;
  mStruct->tstr    = NULL;
  mStruct->nfreq   = NULL;
  mStruct->freq    = NULL;
  mStruct->famp    = NULL;
  mStruct->alfr    = NULL;
  mStruct->knor    = NULL;
  mStruct->snor    = NULL;
  mStruct->xnor    = NULL;
  mStruct->xsnor   = NULL;
  mStruct->ynor    = NULL;
  mStruct->ysnor   = NULL;
  mStruct->cl_mod  = NULL;
  mStruct->cm_mod  = NULL;
  mStruct->cdw_mod = NULL;
  mStruct->cdv_mod = NULL;
  mStruct->cdf_mod = NULL;
  mStruct->modn    = NULL;
  mStruct->dmspn   = NULL;
  mStruct->cl_pos  = NULL;
  mStruct->cm_pos  = NULL;
  mStruct->cdw_pos = NULL;
  mStruct->cdv_pos = NULL;
  mStruct->cdf_pos = NULL;
  mStruct->posn    = NULL;
  mStruct->dpspn   = NULL;
  for (j = 0; j < NBITX; j++) mStruct->isbits[j] = NULL;
  
  status = CAPS_IOERR;
  n = fread(&begMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  n = fread(mStruct->name, sizeof(char), 32, fp);
  if (n != 32) goto cleanup;
  n = fread(&endMarker,    sizeof(int),   1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  for (i = 31; i >= 0; i--)
    if (mStruct->name[i] != ' ') break;
  mStruct->name[i+1] = 0;
  
  n = fread(&begMarker,      sizeof(int),      1, fp);
  if (n !=  1) goto cleanup;
  n = fread(mStruct->istate, sizeof(int), NSTATI, fp);
  if (n != NSTATI) goto cleanup;
  n = fread(&endMarker,      sizeof(int),      1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,      sizeof(int),         1, fp);
  if (n !=  1) goto cleanup;
  n = fread(mStruct->sstate, sizeof(double), NSTATS, fp);
  if (n != NSTATS) goto cleanup;
  n = fread(&endMarker,      sizeof(int),         1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  ii    = mStruct->istate[ 0];
  jj    = mStruct->istate[ 1];
  nbl   = mStruct->istate[ 2];
  nmodn = mStruct->istate[ 9];
  nposn = mStruct->istate[10];
  if (nbl > NBX) {
    AIM_ERROR(aimInfo, "Increase NBX to at least %d\n", nbl);
    status = EGADS_INDEXERR;
    goto cleanup;
  }
  ns = 2*nbl;
  
  n = fread(&begMarker,          sizeof(int), 1, fp);
  if (n !=  1) goto cleanup;
  for (i = 0; i < nbl; i++) {
    n = fread(&mStruct->jbld[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->ninl[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->nout[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->nbld[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->iib[i],  sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->ible[i], sizeof(int), 1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,          sizeof(int), 1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  AIM_ALLOC(mStruct->mfract, jj, double, aimInfo, status);

  AIM_ALLOC(mStruct->x, jj, double*, aimInfo, status);
  for (j = 0; j < jj; j++) mStruct->x[j] = NULL;
  AIM_ALLOC(mStruct->y, jj, double*, aimInfo, status);
  for (j = 0; j < jj; j++) mStruct->y[j] = NULL;
  AIM_ALLOC(mStruct->r, jj, double*, aimInfo, status);
  for (j = 0; j < jj; j++) mStruct->r[j] = NULL;
  AIM_ALLOC(mStruct->h, jj, double*, aimInfo, status);
  for (j = 0; j < jj; j++) mStruct->h[j] = NULL;
  for (j = 0; j < jj; j++) {
    AIM_ALLOC(mStruct->x[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->y[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->r[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->h[j], ii, double, aimInfo, status);
  }
  
  /* 10 DO loop */
  for (j = 0; j < jj; j++) {
    n = fread(&begMarker,          sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    n = fread(&mStruct->mfract[j], sizeof(double), 1, fp);
    if (n != 1) goto cleanup;
    n = fread(&endMarker,          sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,          sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < ii; i++) {
      n = fread(&mStruct->x[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->y[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->r[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->h[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,          sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }
  
  AIM_ALLOC(mStruct->xb, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->xb[j] = NULL;
  AIM_ALLOC(mStruct->yb, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->yb[j] = NULL;
  AIM_ALLOC(mStruct->xpb, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->xpb[j] = NULL;
  AIM_ALLOC(mStruct->ypb, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->ypb[j] = NULL;
  AIM_ALLOC(mStruct->sb, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->sb[j] = NULL;
  AIM_ALLOC(mStruct->sginl, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->sginl[j] = NULL;
  AIM_ALLOC(mStruct->sgout, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->sgout[j] = NULL;
  AIM_ALLOC(mStruct->xw, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->xw[j] = NULL;
  AIM_ALLOC(mStruct->yw, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->yw[j] = NULL;
  AIM_ALLOC(mStruct->wgap, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->wgap[j] = NULL;
  AIM_ALLOC(mStruct->vcen, nbl, double*, aimInfo, status);
  for (j = 0; j < nbl; j++) mStruct->vcen[j] = NULL;

  for (j = 0; j < nbl; j++) {
    AIM_ALLOC(mStruct->xb[j] , mStruct->iib[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->yb[j] , mStruct->iib[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->xpb[j], mStruct->iib[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->ypb[j], mStruct->iib[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->sb[j] , mStruct->iib[j], double, aimInfo, status);

    AIM_ALLOC(mStruct->sginl[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->sgout[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->xw[j]   , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->yw[j]   , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->wgap[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->vcen[j] , ii, double, aimInfo, status);
  }

  /* 20 DO loop */
  for (j = 0; j < nbl; j++) {
    n = fread(&begMarker,            sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < mStruct->iib[j]; i++) {
      n = fread(&mStruct->xb[j][i],  sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->yb[j][i],  sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->xpb[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->ypb[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->sb[j][i],  sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,            sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,              sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < ii; i++) {
      n = fread(&mStruct->sginl[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->sgout[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->xw[j][i],    sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->yw[j][i],    sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->wgap[j][i],  sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,              sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,       sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    n = fread(mStruct->vcen[j], sizeof(double), ii, fp);
    if (n != ii) goto cleanup;
    n = fread(&endMarker,       sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }

  AIM_ALLOC(mStruct->sg, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->sg[j] = NULL;
  AIM_ALLOC(mStruct->disp, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->disp[j] = NULL;
  AIM_ALLOC(mStruct->pspec, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->pspec[j] = NULL;
  AIM_ALLOC(mStruct->thet, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->thet[j] = NULL;
  AIM_ALLOC(mStruct->dstr, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->dstr[j] = NULL;
  AIM_ALLOC(mStruct->uedg, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->uedg[j] = NULL;
  AIM_ALLOC(mStruct->ctau, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->ctau[j] = NULL;
  AIM_ALLOC(mStruct->entr, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->entr[j] = NULL;
  AIM_ALLOC(mStruct->tauw, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->tauw[j] = NULL;
  AIM_ALLOC(mStruct->dint, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->dint[j] = NULL;
  AIM_ALLOC(mStruct->tstr, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->tstr[j] = NULL;
  for (j = 0; j < ns; j++) {
    AIM_ALLOC(mStruct->sg[j]   , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->disp[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->pspec[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->thet[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->dstr[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->uedg[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->ctau[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->entr[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->tauw[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->dint[j] , ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->tstr[j] , ii, double, aimInfo, status);
  }
  
  /* 30 DO loop */
  for (j = 0; j < ns; j++) {
    n = fread(&begMarker,              sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < ii; i++) {
      n = fread(&mStruct->sg[j][i],    sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->disp[j][i],  sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->pspec[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,              sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,             sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < ii; i++) {
      n = fread(&mStruct->thet[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->dstr[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->uedg[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->ctau[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->entr[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->tauw[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->dint[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->tstr[j][i], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,             sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }

  AIM_ALLOC(mStruct->nfreq, ns, int, aimInfo, status);

  n = fread(&begMarker,     sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  n = fread(mStruct->nfreq, sizeof(int), ns, fp);
  if (n != ns) goto cleanup;
  n = fread(&endMarker,     sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;

  for (k = j = 0; j < ns; j++) k += (mStruct->nfreq[j]+1)*ns;
  
  AIM_ALLOC(mStruct->freq, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->freq[j] = NULL;
  AIM_ALLOC(mStruct->famp, k, double*, aimInfo, status);
  for (j = 0; j < k; j++) mStruct->famp[j] = NULL;
  AIM_ALLOC(mStruct->alfr, k, double*, aimInfo, status);
  for (j = 0; j < k; j++) mStruct->alfr[j] = NULL;
  for (j = 0; j < ns; j++) {
    AIM_ALLOC(mStruct->freq[j], (mStruct->nfreq[j]+1), double, aimInfo, status);
  }
  for (j = 0; j < k; j++) {
    AIM_ALLOC(mStruct->famp[j], ii, double, aimInfo, status);
    AIM_ALLOC(mStruct->alfr[j], ii, double, aimInfo, status);
  }
  
  /* 35 DO loop */
  for (k = j = 0; j < ns; j++) {
    for (i = 0; i < mStruct->nfreq[j]+1; i++, k++) {
      n = fread(&begMarker,           sizeof(int),     1, fp);
      if (n !=  1) goto cleanup;
      n = fread(&mStruct->freq[j][i], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(mStruct->famp[k],     sizeof(double), ii, fp);
      if (n != ii) goto cleanup;
      n = fread(mStruct->alfr[k],     sizeof(double), ii, fp);
      if (n != ii) goto cleanup;
      n = fread(&endMarker,           sizeof(int),     1, fp);
      if (n !=  1) goto cleanup;
      if (begMarker != endMarker) goto cleanup;
    }
  }

  AIM_ALLOC(mStruct->knor, ns, int, aimInfo, status);

  n = fread(&begMarker,    sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  n = fread(mStruct->knor, sizeof(int), ns, fp);
  if (n != ns) goto cleanup;
  n = fread(&endMarker,    sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  AIM_ALLOC(mStruct->snor, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->snor[j] = NULL;
  AIM_ALLOC(mStruct->xnor, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->xnor[j] = NULL;
  AIM_ALLOC(mStruct->xsnor, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->xsnor[j] = NULL;
  AIM_ALLOC(mStruct->ynor, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->ynor[j] = NULL;
  AIM_ALLOC(mStruct->ysnor, ns, double*, aimInfo, status);
  for (j = 0; j < ns; j++) mStruct->ysnor[j] = NULL;
  for (j = 0; j < ns; j++) {
    AIM_ALLOC(mStruct->snor[j] , mStruct->knor[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->xnor[j] , mStruct->knor[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->xsnor[j], mStruct->knor[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->ynor[j] , mStruct->knor[j], double, aimInfo, status);
    AIM_ALLOC(mStruct->ysnor[j], mStruct->knor[j], double, aimInfo, status);
  }
  
  /* 40 DO loop */
  for (j = 0; j < ns; j++) {
    n = fread(&begMarker,              sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    for (i = 0; i < mStruct->knor[j]; i++) {
      n = fread(&mStruct->snor[j][i],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->xnor[j][i],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->xsnor[j][i], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->ynor[j][i],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->ysnor[j][i], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,              sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }
  
  n = fread(&begMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  for (j = 0; j < nbl; j++) {
    n = fread(&mStruct->blift[j],  sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->bdrag[j],  sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->bmomn[j],  sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->bdragv[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->bdragf[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  for (j = 0; j < nbl; j++) {
    n = fread(&mStruct->sble[j],   sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->sblold[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->swak[j],   sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->sbcmax[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->sbnose[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  for (j = 0; j < nbl; j++) {
    n = fread(&mStruct->xbnose[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->ybnose[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->xbtail[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->ybtail[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  for (j = 0; j < ns; j++) {
    n = fread(&mStruct->pxx0[j],   sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->pxx1[j],   sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->xtr[j],    sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->xitran[j], sizeof(double),  1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),     1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  if (nmodn > 0) {
    mStruct->cl_mod = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->cl_mod == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CL_MODs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cm_mod = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->cm_mod == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CM_MODs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdw_mod = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->cdw_mod == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDW_MODs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdv_mod = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->cdv_mod == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDV_MODs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdf_mod = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->cdf_mod == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDF_MODs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->modn = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->modn == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d MODNs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->dmspn = (double *) EG_alloc(nmodn*sizeof(double));
    if (mStruct->dmspn == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d DMSPNs\n", nmodn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    
    n = fread(&begMarker,             sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    for (j = 0; j < nmodn; j++) {
      n = fread(&mStruct->cl_mod[j],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cm_mod[j],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdw_mod[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdv_mod[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdf_mod[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->modn[j],    sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->dmspn[j],   sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,             sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }
  
  if (nposn > 0) {
    mStruct->cl_pos = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->cl_pos == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CL_POSs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cm_pos = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->cm_pos == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CM_POSs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdw_pos = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->cdw_pos == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDW_POSs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdv_pos = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->cdv_pos == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDV_POSs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->cdf_pos = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->cdf_pos == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d CDF_POSs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->posn = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->posn == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d POSNs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    mStruct->dpspn = (double *) EG_alloc(nposn*sizeof(double));
    if (mStruct->dpspn == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d DPSPNs\n", nposn);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    
    n = fread(&begMarker,             sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    for (j = 0; j < nposn; j++) {
      n = fread(&mStruct->cl_pos[j],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cm_pos[j],  sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdw_pos[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdv_pos[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->cdf_pos[j], sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->posn[j],    sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->dpspn[j],   sizeof(double),  1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,             sizeof(int),     1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }
  
  n = fread(&begMarker,            sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  for (j = 0; j < ns; j++) {
    n = fread(&mStruct->igfix[j],  sizeof(int),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->igcorn[j], sizeof(int),  1, fp);
    if (n != 1) goto cleanup;
    n = fread(&mStruct->itran[j],  sizeof(int),  1, fp);
    if (n != 1) goto cleanup;
  }
  n = fread(&endMarker,            sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  n = fread(&begMarker, sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  n = fread(&nsuct,     sizeof(int),  1, fp);
  if (n != 1) goto cleanup;
  n = fread(&endMarker, sizeof(int),  1, fp);
  if (n !=  1) goto cleanup;
  if (begMarker != endMarker) goto cleanup;
  
  if (nsuct > 0) {
    n = fread(&begMarker,      sizeof(int),      1, fp);
    if (n !=  1) goto cleanup;
    n = fread(mStruct->issuct, sizeof(int),  nsuct, fp);
    if (n != 1) goto cleanup;
    n = fread(&endMarker,      sizeof(int),      1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
    
    n = fread(&begMarker,               sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    for (j = 0; j < nsuct; j++) {
      n = fread(&mStruct->cqsuct[j],    sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->sgsuct[j][0], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
      n = fread(&mStruct->sgsuct[j][1], sizeof(double), 1, fp);
      if (n != 1) goto cleanup;
    }
    n = fread(&endMarker,               sizeof(int),    1, fp);
    if (n !=  1) goto cleanup;
    if (begMarker != endMarker) goto cleanup;
  }
  
  for (j = 0; j < NBITX; j++) {
    mStruct->isbits[j] = (int *) EG_alloc(ii*sizeof(int));
    if (mStruct->isbits[j] == NULL) {
      AIM_ERROR(aimInfo, "Cannot allocate %d %d ISBITSs\n", j, ii);
      status = EGADS_MALLOC;
      goto cleanup;
    }
  }
  
  /* 50 Do loop */
  for (j = 0; j < NBITX; j++) {
    n = fread(&begMarker,         sizeof(int),  1, fp);
    if (n !=  1) goto ignore;
    n = fread(mStruct->isbits[j], sizeof(int), ii, fp);
    if (n != ii) goto ignore;
    n = fread(&endMarker,         sizeof(int),  1, fp);
    if (n !=  1) goto ignore;
  }

ignore:
  fclose(fp);
  fp     = NULL;
  *mdat  = mStruct;
  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);
  if (status != CAPS_SUCCESS) msesMdatFree(&mStruct);
  return status;
}
