/* HSM unit tester */

#include "egads.h"
#include <math.h>
#include <string.h>

#ifndef WIN32
#define HSMSOL hsmsol_
#define HSMDEP hsmdep_
#define HSMOUT hsmout_
#endif

/* FORTRAN logicals */
static int ffalse = 0, ftrue = 1;

/* these need to be consistent with "index.inc" */
#define IVTOT  7
#define IRTOT  6
#define LVTOT 42
#define LGTOT 24
#define LBTOT 22
#define LPTOT 16
#define JVTOT 25

#define CROSS(a,b,c)  a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                      a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                      a[2] = (b[0]*c[1]) - (b[1]*c[0])

#define PI            3.1415926535897931159979635


extern void HSMSOL(int *lvinit, int *lprint,
                   int *itmax, double *dref,
                   double *dlim, double *dtol, double *ddel,
                   double *alim, double *atol, double *adel,
                   double *parg,
                   int *nnode, double *par, double *var, double *dep,
                   int *nelem, int *kelem,
                   int *nbcedge, int *kbcedge, double *parb,
                   int *nbcnode, int *kbcnode, double *parp, int *lbcnode,
                   int *njoint, int *kjoint,
                   int *kdim, int *ldim, int *nedim, int *nddim, int *nmdim,
                   double *bf, double *bf_dj, double *bm, double *bm_dj,
                   int *ibr1, int *ibr2, int *ibr3, int *iba1, int *iba2,
                   double *resc, double *resc_var,
                   double *resp, double *resp_var, double *resp_dvp,
                   int *kdvp, int *ndvp, double *ares,
                   int *ifrst, int *ilast, int *mfrst,
                   /*@null@*/ double *amat, /*@null@*/ int *ipp, double *dvar);
extern void HSMDEP(int *leinit, int *lprint,
                   int *itmax,
                   double *elim, double *etol, double *edel,
                   int *nnode, double *par, double *var, double *dep,
                   int *nelem, int *kelem,
                   int *kdim, int *ldim, int *nedim, int *nddim, int *nmdim,
                   double *acn, double *resn, double *rese, double *rese_de,
                   double *rest, double *rest_t,
                   int *kdt, int *ndt,
                   int *ifrstt, int *ilastt, int *mfrstt, double *amatt,
                   double *resv, double *resv_v, int *kdv, int *ndv,
                   int *ifrstv, int *ilastv, int *mfrstv, double *amatv);
extern void HSMOUT(int *nelem, int *kelem, double *var, double *dep,
                   double *par, double *parg, int *kdim, int *ldim, int *nedim,
                   int *nddim, int *nmdim);

extern int  triangle_rcm(int node_num, int triangle_num, int *triangle_node,
                         /*@null@*/ int *perm, /*@null@*/ int *perm_inv,
                         /*@null@*/ int *segs);
extern int  maxValence(int nvert, int ntri, int *tris);


static void
ortmat(double  e1, double  e2, double g12, double v12, double c16, double c26,
       double g13, double g23, double tsh, double zrf, double srfac,
       double  *a, double  *b, double  *d, double *s)
{
/* -------------------------------------------------------------------------
       Sets stiffness matrices A,B,D  and shear-compliance matrix S
       for an orthotropic shell, augmented with shear/extension coupling
 
       Inputs:
         e1    modulus in 1 direction
         e2    modulus in 2 direction
         g12   shear modulus
         v12   Poisson's ratio
 
         c16   12-shear / 1-extension coupling modulus
         c26   12-shear / 2-extension coupling modulus
 
         g13   1-direction transverse-shear modulus
         g23   2-direction transverse-shear modulus
 
         tsh   shell thickness
         zrf   reference surface location parameter -1 .. +1
         srfac transverse-shear strain energy reduction factor
                ( = 5/6 for isotropic shell)
 
      Outputs:
         a(.)   stiffness  tensor components  A11, A22, A12, A16, A26, A66
         b(.)   stiffness  tensor components  B11, B22, B12, B16, B26, B66
         d(.)   stiffness  tensor components  D11, D22, D12, D16, D26, D66
         s(.)   compliance tensor components  S55, S44
    ------------------------------------------------------------------------- */
  int    i;
  double den, tfac1, tfac2, tfac3, econ[6], scon[2];
  
  /* ---- in-plane stiffnesses  */
  den     = 1.0 - v12*v12 * e2/e1;
  econ[0] = e1/den;                /* c11 */
  econ[1] = e2/den;                /* c22 */
  econ[2] = e2/den * v12;          /* c12 */
  econ[3] = c16;                   /* c16 */
  econ[4] = c26;                   /* c26 */
  econ[5] = 2.0*g12;               /* c66 */
  
  /* ---- transverse shear complicances */
  scon[0] = 1.0/g13;               /* s55 */
  scon[1] = 1.0/g23;               /* s44 */
  
  /*---- elements of in-plane stiffness matrices A,B,D for homogeneous shell */
  tfac1 =  tsh;
  tfac2 = -tsh*tsh     * zrf / 2.0;
  tfac3 =  tsh*tsh*tsh * (1.0 + 3.0*zrf*zrf) / 12.0;
  for (i = 0; i < 6; i++) {
    a[i] = econ[i] * tfac1;
    b[i] = econ[i] * tfac2;
    d[i] = econ[i] * tfac3;
  }
  
  s[0] = scon[0] / (srfac*tsh);
  s[1] = scon[1] / (srfac*tsh);
}


int main(int argc, char *argv[])
{
  int          i, j, k, status, oclass, mtype, nbody, nvert, nface, ntri, nseg;
  int          plen, tlen, pty1, pin1, pty2, pin2, atype, alen, nddim;
  int          *vtable, *trin, *perm, *perm_inv, *segs, *sens;
  float        arg;
  double       params[3], box[6], result[18], norm[3], size;
  ego          context, model, geom, *bodies, tess, *dum, *faces = NULL;
  const char   *OCCrev, *string;
  const double *points, *uv, *real;
  const int    *tris, *tric, *ptype, *pindex, *integer;

  if ((argc != 2) && (argc != 5)) {
    printf(" Usage: testHSM Model [angle relSide relSag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  status = EG_open(&context);
  if (status != EGADS_SUCCESS) {
    printf(" EG_open = %d!\n\n", status);
    return 1;
  }
  status = EG_loadModel(context, 0, argv[1], &model);
  if (status != EGADS_SUCCESS) {
    printf(" EG_loadModel = %d\n\n", status);
    return 1;
  }
  status = EG_getBoundingBox(model, box);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getBoundingBox = %d\n\n", status);
    return 1;
  }
  size = sqrt((box[0]-box[3])*(box[0]-box[3]) + (box[1]-box[4])*(box[1]-box[4]) +
              (box[2]-box[5])*(box[2]-box[5]));

  /* get all bodies */
  status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                          &bodies, &sens);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n\n", status);
    return 1;
  }

  params[0] =  0.025*size;
  params[1] =  0.001*size;
  params[2] = 15.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
           params[2], params[0], params[1]);
    params[0] *= size;
    params[1] *= size;
  }

  printf(" Number of Bodies = %d\n\n", nbody);
  
  for (i = 0; i < nbody; i++) {
    
    status = EG_getTopology(bodies[i], &geom, &oclass, &mtype, NULL, &j, &dum,
                            &sens);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getTopology Body %d = %d\n", i+1, status);
      continue;
    }
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", i+1);
      continue;
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", i+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", i+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", i+1);
    }
    
    status = EG_makeTessBody(bodies[i], params, &tess);
    if (status != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", i, status);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos Face %d = %d\n", i+1, status);
      EG_deleteObject(tess);
      continue;
    }
    
    for (nvert = ntri = j = 0; j < nface; j++) {
      status = EG_getTessFace(tess, j+1, &plen, &points, &uv, &ptype, &pindex,
                              &tlen, &tris, &tric);
      if (status != EGADS_SUCCESS) {
        printf(" %d EG_getTessFace %d = %d\n", i+1, j, status);
        continue;
      }
      ntri  += tlen;
      nvert += plen;
    }
    printf(" nvert = %d, ntris = %d\n", nvert, ntri);
    
    vtable   = (int *) EG_alloc(4*nvert*sizeof(int));
    trin     = (int *) EG_alloc(3*ntri*sizeof(int));
    perm     = (int *) EG_alloc( nvert*sizeof(int));
    perm_inv = (int *) EG_alloc( nvert*sizeof(int));
    if ((vtable   == NULL) || (trin == NULL) || (perm == NULL) ||
        (perm_inv == NULL)) {
      if (perm_inv != NULL) EG_free(perm_inv);
      if (perm     != NULL) EG_free(perm);
      if (trin     != NULL) EG_free(trin);
      if (vtable   != NULL) EG_free(vtable);
      EG_deleteObject(tess);
      EG_free(faces);
      printf(" Error:: cannot allocate memory!\n");
      continue;
    }
    
    for (nvert = ntri = j = 0; j < nface; j++) {
      status = EG_getTessFace(tess, j+1, &plen, &points, &uv, &ptype, &pindex,
                              &tlen, &tris, &tric);
      if (status != EGADS_SUCCESS) {
        printf(" %d EG_getTessFace %d = %d\n", i+1, j, status);
        continue;
      }
      for (k = 0; k < tlen; k++, ntri++) {
        trin[3*ntri  ] = nvert + tris[3*k  ];
        trin[3*ntri+1] = nvert + tris[3*k+1];
        trin[3*ntri+2] = nvert + tris[3*k+2];
      }
      for (k = 0; k < plen; k++, nvert++)
        if (ptype[k] < 0) {
          vtable[4*nvert  ] = -(k+1);
          vtable[4*nvert+1] =   j+1;
          vtable[4*nvert+2] = vtable[4*nvert+3] = 0;
        } else {
          vtable[4*nvert  ] = ptype[k];
          vtable[4*nvert+1] = pindex[k];
          vtable[4*nvert+2] = vtable[4*nvert+3] = 0;
        }
    }
    
    /* mark up shared vertices */
    for (j = 0; j < nvert; j++)
      for (k = 0; k < j; k++) {
        if ((vtable[4*j  ] == vtable[4*k  ]) &&
            (vtable[4*j+1] == vtable[4*k+1])) {
          vtable[4*j+2] = k+1;
          if (vtable[4*k+2] == 0) vtable[4*k+2] = -1;
          break;
        }
      }
    
    /* get the maximum valence */
    nddim = maxValence(nvert, ntri, trin);
    if (nddim <= 0) {
      EG_free(perm_inv);
      EG_free(perm);
      EG_free(trin);
      EG_free(vtable);
      EG_free(faces);
      EG_deleteObject(tess);
      printf(" Error:: maxValence = %d!\n", nddim);
      continue;
    }
    /* count ourself */
    nddim += 10;
    printf(" max Valence = %d\n", nddim);
    
    /* get the number of open segments */
    nseg = triangle_rcm(-nvert, -ntri, trin, NULL, NULL, NULL);
    printf(" nOpen Segments = %d\n\n", nseg);
    segs = NULL;
    if (nseg > 0) {
      segs = (int *) EG_alloc(3*nseg*sizeof(int));
      if (segs == NULL) {
        EG_free(perm_inv);
        EG_free(perm);
        EG_free(trin);
        EG_free(vtable);
        EG_free(faces);
        EG_deleteObject(tess);
        printf(" Error:: cannot allocate segments!\n");
        continue;
      }
    }
    
    /* get the numbering! */
    status = triangle_rcm(-nvert, ntri, trin, perm, perm_inv, segs);
    if (status < 0)
      printf(" Error: triangle_rcm = %d\n", status);
    
    /* check for self-consistency & uniqueness */
    if (status >= 0) {
      for (j = 0; j < nvert; j++)
        if (j+1 != perm_inv[perm[j]-1])
          printf("  BAD Index %d: %d %d\n", j+1, perm[j], perm_inv[perm[j]-1]);
      
      for (j = 0; j < nvert; j++) perm_inv[j] = 0;
      for (j = 0; j < nvert; j++) perm_inv[perm[j]-1]++;
      for (j = 0; j < nvert; j++)
        if (perm_inv[j] != 1)
          printf("  BAD Index %d: #hits = %d\n", j+1, perm_inv[j]);
      printf("\n");

      /* output open segments */
      if ((nseg > 0) && (segs != NULL)) {
        status = EG_getBodyTopos(bodies[i], NULL, EDGE, &k, &dum);
        if (status == EGADS_SUCCESS) {
          for (j = 0; j < nseg; j++) {
            k     = segs[3*j  ] - 1;
            pty1  = vtable[4*k  ];
            pin1  = vtable[4*k+1];
            atype = vtable[4*k+2];
            k     = segs[3*j+1] - 1;
            pty2  = vtable[4*k  ];
            pin2  = vtable[4*k+1];
            if ((atype != 0) && (vtable[4*k+2] != 0)) continue;
            printf(" Open side = %4d: %4d (%d,%d)",
                   segs[3*j+2], segs[3*j  ], pty1, pin1);

            printf(" %4d (%d,%d)", segs[3*j+1], pty2, pin2);
            if ((pty1 == 0) && (pty2 == 0)) {
              printf("       No Edge Indicator!\n");
            } else {
              k = pin1;
              if (pty1 == 0) k = pin2;
              status = EG_attributeRet(dum[k-1], "HSMbc", &atype, &alen,
                                       &integer, &real, &string);
              if (status != EGADS_SUCCESS) {
                printf("      attribute status = %d\n", status);
              } else {
                if (atype == ATTRSTRING) {
                  printf("      Edge attribute = %s\n", string);
                  /* fix the root "nodes" */
                  if (strcmp(string,"root") == 0) {
                    k = segs[3*j  ] - 1;
                    vtable[4*k+3] = 3 + 20;
                    k = segs[3*j+1] - 1;
                    vtable[4*k+3] = 3 + 20;
                  }
                } else {
                  printf("      attribute type = %d\n", atype);
                }
              }
            }
          }
          EG_free(dum);
        }
      }
    }
  
    /* mimic case 8 in "hsmrun.f" */
    
    if ((status >= 0) && (faces != NULL)) {
      int    itmaxv, itmaxe, ldim, nmdim, jj, kk;
      int    nbcedge = 0, nbcnode = 0, njoint = 0;
      int    *kbcedge = NULL, *kbcnode = NULL, *lbcnode = NULL, *kjoint = NULL;
      int    *kelem = NULL, *ibx = NULL, *kdvp = NULL, *ndvp = NULL;
      int    *frst = NULL, *ipp = NULL, *idt = NULL, *frstt = NULL;
      int    *kdv = NULL, *ndv = NULL, *frstv = NULL;
      double *var = NULL, *dep = NULL, *par = NULL, *parb = NULL;
      double *parp = NULL;
      double *bf = NULL, *bf_dj = NULL, *bm = NULL, *bm_dj = NULL;
      double *resc = NULL, *resc_var = NULL, *ares = NULL;
      double *resp = NULL, *resp_var = NULL, *resp_dvp = NULL;
      double *amat = NULL, *dvar = NULL;
      double *res = NULL, *rest = NULL, *rest_t = NULL;
      double *amatt = NULL, *resv = NULL, *resv_v = NULL, *amatv = NULL;
      double atol, dtol, dref, dlim, alim, adel, ddel, elim, etol, edel, cp;
      double mate1, mate2, matg12, matv12, matc16, matc26, matg13, matg23;
      double tshell, zetref, srfac;
      double parg[LGTOT];
      
      /* global parameters */
      for (j = 0; j < LGTOT; j++) parg[j] = 0.0;
  
      /* Newton convergence tolerances */
      dtol = 1.0e-11;  /* relative displacements, |d|/dref */
      atol = 1.0e-11;  /* angles (unit vector changes) */
      
      /* reference length for displacement limiting, convergence checks
         (should be comparable to size of the geometry) */
      dref = 1.0;
      /* max Newton changes (dimensionless) */
      dlim = 1.0;
      alim = 1.0;
      
      mate1  = 4.0e5;
      mate2  = 4.0e5;
      matv12 = 0.3;
      matg12 = mate1*0.5/(1.0+matv12);
      matg13 = matg12;
      matg23 = matg12;
      matc16 = 0.0;
      matc26 = 0.0;
      tshell = 0.075;
      zetref = 0.0;
      srfac  = 5.0/6.0;

      itmaxv = 35;
      kelem  = (int *) EG_alloc(4*ntri*sizeof(int));
      if (kelem == NULL) goto out;
      for (j = 0; j < ntri; j++) {
        kelem[4*j  ] = perm[trin[3*j  ]-1];
        kelem[4*j+1] = perm[trin[3*j+1]-1];
        kelem[4*j+2] = perm[trin[3*j+2]-1];
        kelem[4*j+3] = 0;
      }
      for (j = 0; j < nvert; j++) {
        if (vtable[4*j+2] >  0) njoint++;
        if (vtable[4*j+3] != 0) nbcnode++;
      }
      printf(" nJoint = %d  nBCnode = %d\n", njoint, nbcnode);
      
                          ldim = nbcedge;
      if (ldim < nbcnode) ldim = nbcnode;
      if (ldim < njoint)  ldim = njoint;
      if (ldim == 0)      ldim = 1;
      
      var      = (double *) EG_alloc(IVTOT*nvert*sizeof(double));
      if (var      == NULL) goto out;
      dep      = (double *) EG_alloc(JVTOT*nvert*sizeof(double));
      if (dep      == NULL) goto out;
      par      = (double *) EG_alloc(LVTOT*nvert*sizeof(double));
      if (par      == NULL) goto out;
      parb     = (double *) EG_alloc(LBTOT*ldim *sizeof(double));
      if (parb     == NULL) goto out;
      parp     = (double *) EG_alloc(LPTOT*ldim *sizeof(double));
      if (parp     == NULL) goto out;
      kbcedge  = (int *)    EG_alloc(    2*ldim *sizeof(int));
      if (kbcedge  == NULL) goto out;
      kbcnode  = (int *)    EG_alloc(      ldim *sizeof(int));
      if (kbcnode  == NULL) goto out;
      lbcnode  = (int *)    EG_alloc(      ldim *sizeof(int));
      if (lbcnode  == NULL) goto out;
      kjoint   = (int *)    EG_alloc(    2*ldim *sizeof(int));
      if (kjoint   == NULL) goto out;
      
      /* fill in the joints */
      for (njoint = j = 0; j < nvert; j++) {
        if (vtable[4*j+2] <= 0) continue;
        kjoint[2*njoint  ] = perm[j];
        kjoint[2*njoint+1] = perm[vtable[4*j+2]-1];
        njoint++; // Durscher change.
      }
      
      /* temporary storage */
      bf       = (double *) EG_alloc(    9*nvert*sizeof(double));
      if (bf       == NULL) goto out;
      bf_dj    = (double *) EG_alloc(   27*nvert*sizeof(double));
      if (bf_dj    == NULL) goto out;
      bm       = (double *) EG_alloc(    6*nvert*sizeof(double));
      if (bm       == NULL) goto out;
      bm_dj    = (double *) EG_alloc(   18*nvert*sizeof(double));
      if (bm_dj    == NULL) goto out;
      resc     = (double *) EG_alloc(IVTOT*nvert*sizeof(double));
      if (resc     == NULL) goto out;
      resc_var = (double *) EG_alloc(IVTOT*IVTOT*nvert*nddim*sizeof(double));
      if (resc_var == NULL) goto out;
      resp     = (double *) EG_alloc(IRTOT*nvert*sizeof(double));
      if (resp     == NULL) goto out;
      resp_var = (double *) EG_alloc(IRTOT*IVTOT*nddim*sizeof(double));
      if (resp_var == NULL) goto out;
      resp_dvp = (double *) EG_alloc(IRTOT*IRTOT*nvert*nddim*sizeof(double));
      if (resp_dvp == NULL) goto out;
      ares     = (double *) EG_alloc(nvert*sizeof(double));
      if (ares     == NULL) goto out;
      dvar     = (double *) EG_alloc(IVTOT*nvert*sizeof(double));
      if (dvar     == NULL) goto out;
      res      = (double *) EG_alloc(6*nvert*sizeof(double));
      if (res      == NULL) goto out;
      rest     = (double *) EG_alloc(3*4*nvert*sizeof(double));
      if (rest     == NULL) goto out;
      rest_t   = (double *) EG_alloc(3*3*nddim*nvert*sizeof(double));
      if (rest_t   == NULL) goto out;
      resv     = (double *) EG_alloc(2*2*nvert*sizeof(double));
      if (resv     == NULL) goto out;
      resv_v   = (double *) EG_alloc(2*2*nddim*nvert*sizeof(double));
      if (resv_v   == NULL) goto out;
      
      ibx   = (int *) EG_alloc(5*ldim*sizeof(int));
      if (ibx   == NULL) goto out;
      kdvp  = (int *) EG_alloc(nddim*nvert*sizeof(int));
      if (kdvp  == NULL) goto out;
      ndvp  = (int *) EG_alloc(nvert*sizeof(int));
      if (ndvp  == NULL) goto out;
      frst  = (int *) EG_alloc((3*nvert+1)*sizeof(int));
      if (frst  == NULL) goto out;
      idt   = (int *) EG_alloc((nddim+1)*nvert*sizeof(int));
      if (idt   == NULL) goto out;
      frstt = (int *) EG_alloc((3*nvert+1)*sizeof(int));
      if (frstt == NULL) goto out;
      kdv   = (int *) EG_alloc(nddim*nvert*sizeof(int));
      if (kdv   == NULL) goto out;
      ndv   = (int *) EG_alloc(nvert*sizeof(int));
      if (ndv   == NULL) goto out;
      frstv = (int *) EG_alloc((3*nvert+1)*sizeof(int));
      if (frstv == NULL) goto out;
      
      for (j = 0; j < nvert; j++) {
        for (k = 0; k < IVTOT; k++) var[j*IVTOT+k] = 0.0;
        for (k = 0; k < JVTOT; k++) dep[j*JVTOT+k] = 0.0;
        for (k = 0; k < LVTOT; k++) par[j*LVTOT+k] = 0.0;
      }
      for (j = 0; j < LPTOT*ldim; j++) parp[j] = 0.0;
      for (j = 0; j < LBTOT*ldim; j++) parb[j] = 0.0;
      
      /* fill in the BCs */
      
      for (j = k = 0; k < nface; k++) {
        status = EG_getTessFace(tess, k+1, &plen, &points, &uv, &ptype, &pindex,
                                &tlen, &tris, &tric);
        if (status != EGADS_SUCCESS) {
          printf(" %d EG_getTessFace %d = %d\n", i+1, k, status);
          continue;
        }
        for (jj = 0; jj < plen; jj++, j++) {
          double *e1, *e2, normalize;
        
          e1     = &result[3];
          e2     = &result[6];
          kk     = perm[j] - 1;
          status = EG_evaluate(faces[k], &uv[2*jj], result);
          if (status != EGADS_SUCCESS) {
            printf(" %d EG_evaluate %d %d = %d\n", i+1, k+1, jj+1, status);
            continue;
          }
          par[kk*LVTOT   ] = result[0];
          par[kk*LVTOT+ 1] = result[1];
          par[kk*LVTOT+ 2] = result[2];
          normalize        = sqrt(result[3]*result[3] + result[4]*result[4] +
                                  result[5]*result[5]);
          if (normalize == 0.0) {
            printf(" %d  Face %d Vert %d e0_1 is degenerate!\n", i+1,k+1,jj+1);
            normalize = 1.0;
          }
          par[kk*LVTOT+ 3] = result[3]/normalize;
          par[kk*LVTOT+ 4] = result[4]/normalize;
          par[kk*LVTOT+ 5] = result[5]/normalize;
          normalize        = sqrt(result[6]*result[6] + result[7]*result[7] +
                                  result[8]*result[8]);
          if (normalize == 0.0) {
            printf(" %d  Face %d Vert %d e0_2 is degenerate!\n", i+1,k+1,jj+1);
            normalize = 1.0;
          }
          par[kk*LVTOT+ 6] = result[6]/normalize;
          par[kk*LVTOT+ 7] = result[7]/normalize;
          par[kk*LVTOT+ 8] = result[8]/normalize;
          CROSS(norm, e1, e2);
          if (faces[k]->mtype == SREVERSE) {
            norm[0] = -norm[0];
            norm[1] = -norm[1];
            norm[2] = -norm[2];
          }
          normalize        = sqrt(norm[0]*norm[0] + norm[1]*norm[1] +
                                  norm[2]*norm[2]);
          if (normalize == 0.0) {
            printf(" %d  Face %d Vert %d norm is degenerate!\n", i+1,k+1,jj+1);
            normalize = 1.0;
          }
          par[kk*LVTOT+ 9] = norm[0]/normalize;
          par[kk*LVTOT+10] = norm[1]/normalize;
          par[kk*LVTOT+11] = norm[2]/normalize;
          
          /* set material properties */
          ortmat(mate1, mate2, matg12, matv12, matc16, matc26, matg13, matg23,
                 tshell, zetref, srfac, &par[kk*LVTOT+12], &par[kk*LVTOT+18],
                                        &par[kk*LVTOT+24], &par[kk*LVTOT+30]);
          
          /* normal-loading Cp distribution */
/*        cp = sin(2.0*pi*u) - sin(pi*u) + 2.0*sin(pi*u)**32  */
          cp = 1.0;
          par[kk*LVTOT+32] = -cp;
        }
      }
      for (k = j = 0; j < nvert; j++) {
        if (vtable[4*j+3] == 0) continue;
        kbcnode[k] = perm[j];
        lbcnode[k] = vtable[4*j+3];
        for (kk = 0; kk < 3; kk++)
          parp[k*LPTOT+kk]   = par[(perm[j]-1)*LVTOT+kk];
        for (kk = 3; kk < 9; kk++)
          parp[k*LPTOT+kk+6] = par[(perm[j]-1)*LVTOT+kk];
        k++;
      }
      
      //for (i = 0; i < LVTOT*nvert; i++) printf("i = %d, %f\n", i, par[i]);
      //for (i = 0; i < LPTOT*ldim; i++) printf("i = %d, %f\n", i, parp[i]);
      //for (i = 0; i < ldim; i++) printf("i = %d, %d\n", i, kbcnode[i]);
      //for (i = 0; i < ldim; i++) printf("i = %d, %d\n", i, lbcnode[i]);
      //for (i = 0; i < LBTOT*ldim; i++) printf("i = %d, %f\n", i, parb[i]);

      for (i = 0; i < 2*ldim; i++) printf("i = %d, %d\n", i, kjoint[i]);

      /* get the size of the matrix */
      k     = -2;
      nmdim =  1;
      printf("%d %d %d %d %d\n", nvert, ldim, ntri, nddim, nmdim);

      HSMSOL(&ffalse, &ftrue, &k, &dref, &dlim, &dtol, &ddel, &alim, &atol,
             &adel, parg, &nvert, par, var, dep, &ntri, kelem, &nbcedge, kbcedge,
             parb, &nbcnode, kbcnode, parp, lbcnode, &njoint, kjoint,
             &nvert, &ldim, &ntri, &nddim, &nmdim, bf, bf_dj, bm, bm_dj,
             &ibx[0], &ibx[ldim], &ibx[2*ldim], &ibx[3*ldim], &ibx[4*ldim],
             resc, resc_var, resp, resp_var, resp_dvp, kdvp, ndvp, ares,
             &frst[0], &frst[nvert], &frst[2*nvert], amat, ipp, dvar);
      
      /* allocate the larger matrix storage after probing for the size */
      printf(" Matrix Size = %d\n", nmdim);
      amat  = (double *) EG_alloc(IRTOT*IRTOT*nmdim*sizeof(double));
      if (amat  == NULL) goto out;
      amatt = (double *) EG_alloc(3*3*nmdim*sizeof(double));
      if (amatt == NULL) goto out;
      amatv = (double *) EG_alloc(2*2*nmdim*sizeof(double));
      if (amatv == NULL) goto out;
      ipp   = (int *)    EG_alloc(IRTOT*nmdim*sizeof(int));
      if (ipp   == NULL) goto out;
      
      HSMSOL(&ffalse, &ftrue, &itmaxv, &dref, &dlim, &dtol, &ddel, &alim, &atol,
             &adel, parg, &nvert, par, var, dep, &ntri, kelem, &nbcedge, kbcedge,
             parb, &nbcnode, kbcnode, parp, lbcnode, &njoint, kjoint,
             &nvert, &ldim, &ntri, &nddim, &nmdim, bf, bf_dj, bm, bm_dj,
             &ibx[0], &ibx[ldim], &ibx[2*ldim], &ibx[3*ldim], &ibx[4*ldim],
             resc, resc_var, resp, resp_var, resp_dvp, kdvp, ndvp, ares,
             &frst[0], &frst[nvert], &frst[2*nvert], amat, ipp, dvar);
      if (itmaxv >= 0) {
        itmaxe = 20;
        elim   = 1.0;
        etol   = atol;
        HSMDEP(&ffalse, &ffalse, &itmaxe, &elim, &etol, &edel, &nvert, par, var,
               dep, &ntri, kelem, &nvert, &ldim, &ntri, &nddim, &nmdim,
               &res[0], &res[nvert], &res[4*nvert], &res[5*nvert], rest, rest_t,
               &idt[nvert], &idt[0], &frstt[0], &frstt[nvert], &frstt[2*nvert],
               amatt, resv, resv_v, kdv, ndv, &frstv[0], &frstv[nvert],
               &frstv[2*nvert], amatv);
        HSMOUT(&ntri, kelem, var, dep, par, parg,
               &nvert, &ldim, &ntri, &nddim, &nmdim);
      }
      
      /* use perm_inv to get hsm indices back to the tessellation object */
      
      
    out:
      if (frstv    != NULL) EG_free(frstv);
      if (ndv      != NULL) EG_free(ndv);
      if (kdv      != NULL) EG_free(kdv);
      if (frstt    != NULL) EG_free(frstt);
      if (idt      != NULL) EG_free(idt);
      if (ipp      != NULL) EG_free(ipp);
      if (frst     != NULL) EG_free(frst);
      if (ndvp     != NULL) EG_free(ndvp);
      if (kdvp     != NULL) EG_free(kdvp);
      if (ibx      != NULL) EG_free(ibx);
      
      if (amatv    != NULL) EG_free(amatv);
      if (resv_v   != NULL) EG_free(resv_v);
      if (resv     != NULL) EG_free(resv);
      if (amatt    != NULL) EG_free(amatt);
      if (rest_t   != NULL) EG_free(rest_t);
      if (rest     != NULL) EG_free(rest);
      if (res      != NULL) EG_free(res);
      if (dvar     != NULL) EG_free(dvar);
      if (amat     != NULL) EG_free(amat);
      if (ares     != NULL) EG_free(ares);
      if (resp_dvp != NULL) EG_free(resp_dvp);
      if (resp_var != NULL) EG_free(resp_var);
      if (resp     != NULL) EG_free(resp);
      if (resc_var != NULL) EG_free(resc_var);
      if (resc     != NULL) EG_free(resc);
      if (bm_dj    != NULL) EG_free(bm_dj);
      if (bm       != NULL) EG_free(bm);
      if (bf_dj    != NULL) EG_free(bf_dj);
      if (bf       != NULL) EG_free(bf);
      
      if (kjoint   != NULL) EG_free(kjoint);
      if (lbcnode  != NULL) EG_free(lbcnode);
      if (kbcnode  != NULL) EG_free(kbcnode);
      if (kbcedge  != NULL) EG_free(kbcedge);
      if (parp     != NULL) EG_free(parp);
      if (parb     != NULL) EG_free(parb);
      if (par      != NULL) EG_free(par);
      if (dep      != NULL) EG_free(dep);
      if (var      != NULL) EG_free(var);
      if (kelem    != NULL) EG_free(kelem);
    }

    if (segs != NULL) EG_free(segs);
    EG_free(perm_inv);
    EG_free(perm);
    EG_free(trin);
    EG_free(vtable);
    
    EG_deleteObject(tess);
    EG_free(faces);
  }
  printf("\n");

  EG_deleteObject(model);
  EG_close(context);
  
  return 0;
}
