
       subroutine hsmrfm(neln, lrcurv, ldrill,
     &                   par1    , par2    , par3    , par4    , 
     &                   var1    , var2    , var3    , var4    , 
     &                   dep1    , dep2    , dep3    , dep4    , 
     &             reseps, 
     &             reskap, 
     &             resfor, 
     &             resmom, rest_t,
     &             resfn ,
     &             resgn , resv_v, ac )
c---------------------------------------------------------------------------
c     Sets up residuals and Jacobians for f,m tensor definitions 
c     in node e1,e2,n bases, for a tri or quad element.
c
c Inputs:
c   neln        number of nodes (3 for tri, 4 for quad)
c   par#(.)   node parameters
c   var#(.)   node primary variables
c   dep#(.)   node dependent variables
c
c   Note: "4" quantities will not influence results for tri element
c
c
c Outputs:
c   reseps(.1)  eps11 residuals         int int [ Re  ] Wi/J^2 dA
c   reseps(.2)  eps22 residuals           
c   reseps(.3)  eps12 residuals           
c
c   reskap(.1)  kap11 residuals         int int [ Rk  ] Wi/J^2 dA
c   reskap(.2)  kap22 residuals           
c   reskap(.3)  kap12 residuals           
c
c   resfor(.1)  f11 residuals           int int [ Rf  ] Wi/J^2 dA
c   resfor(.2)  f22 residuals           
c   resfor(.3)  f12 residuals           
c
c   resmom(.1)  m11 residuals           int int [ Rm  ] Wi/J^2 dA
c   resmom(.2)  m22 residuals    
c   resmom(.3)  m12 residuals    

c   rest_t(....) = d(reseps)/d(e)    tensor Jacobian  (node mass matrix)
c                = d(reskap)/d(k)
c                = d(resfor)/d(f)
c                = d(resmom)/d(m)
c
c   resfn(.1)  f1n shear residuals    int int [ Rfn ] Wi/J^2 dA
c   resfn(.2)  f2n shear residuals    
c
c   resgn(.1)  gam1 tilt residuals    int int [ Rgn ] Wi/J^2 dA
c   resgn(.2)  gam2 tilt residuals    
c
c   resv_v(....) = d(resfn)/d(fn)    vector Jacobian  (node mass matrix)
c                = d(resga)/d(ga)
c
c   ac(.)    area integrals of weight functions,  int int [ 1 ] Wi/J^2 dA
c
c
c   Note: The Jacobians depend only on the undeformed geometry in par(..),
c         and NOT on the solution in var(..) or dep(..)
c
c---------------------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c---- gauss-integration locations and weights
      include 'gausst.inc'
      include 'gaussq.inc'

c----------------------------------------------
c---- passed variables
      logical lrcurv, ldrill

      real par1(lvtot),
     &     par2(lvtot),
     &     par3(lvtot),
     &     par4(lvtot)

      real var1(ivtot),
     &     var2(ivtot),
     &     var3(ivtot),
     &     var4(ivtot)

      real dep1(jvtot),
     &     dep2(jvtot),
     &     dep3(jvtot),
     &     dep4(jvtot)

      real reseps(3,4),
     &     reskap(3,4),
     &     resfor(3,4),
     &     resmom(3,4)
      real rest_t(3,4,3,4)

      real resfn(2,4),
     &     resgn(2,4)
      real resv_v(2,4,2,4)

      real ac(4)

c----------------------------------------------
c---- local work arrays
      real par(lvtot,4),
     &     var(ivtot,4),
     &     dep(jvtot,4)

      real e1a(2), e2a(2)

      real aaj(2,2,2,2,4), bbj(2,2,2,2,4),  ddj(2,2,2,2,4)
      real ssj(2,2,4)

      real rj(3,4), r0j(3,4)
      real dj(3,4), d0j(3,4)
      real pj(4), p0j(4)

      real ej(3,2,4), e0j(3,2,4)

      real epsj(2,2,4),
     &     kapj(2,2,4)
      real fj(2,2,4),
     &     mj(2,2,4)
      real fnj(2,4),
     &     gnj(2,4)
 
c---- element cartesian basis work arrays
      real d1r0(3), d2r0(3)

      real crs(3)
      real c1(3)
      real c2(3)
      real c3(3)

      real cai(2,2)

c- - - - - - - - - - - - - - - - - - - - - - - - - - 
c---- Gauss point quantities for calling HSMGEO
      parameter (ngdim = ngaussq)

      real xsig(ngdim), etag(ngdim), fwtg(ngdim)

      real r0(3,ngdim)
      real r(3,ngdim), 
     &     r_rj(3,3,4,ngdim), r_dj(3,3,4,ngdim), r_pj(3,4,ngdim)

      real d0(3,ngdim)
      real d(3,ngdim), 
     &     d_rj(3,3,4,ngdim), d_dj(3,3,4,ngdim), d_pj(3,4,ngdim)

      real a01(3,ngdim)
      real a1(3,ngdim), 
     &     a1_rj(3,3,4,ngdim), a1_dj(3,3,4,ngdim), a1_pj(3,4,ngdim)

      real a02(3,ngdim)
      real a2(3,ngdim),
     &     a2_rj(3,3,4,ngdim), a2_dj(3,3,4,ngdim), a2_pj(3,4,ngdim)

      real a01i(3,ngdim)
      real a1i(3,ngdim),
     &     a1i_rj(3,3,4,ngdim), a1i_dj(3,3,4,ngdim), a1i_pj(3,4,ngdim)

      real a02i(3,ngdim)
      real a2i(3,ngdim),
     &     a2i_rj(3,3,4,ngdim), a2i_dj(3,3,4,ngdim), a2i_pj(3,4,ngdim)

      real g0(2,2,ngdim)
      real g(2,2,ngdim),
     &     g_rj(2,2,3,4,ngdim), g_dj(2,2,3,4,ngdim), g_pj(2,2,4,ngdim)

      real g0det(ngdim)
      real gdet(ngdim),
     &     gdet_rj(3,4,ngdim), gdet_dj(3,4,ngdim), gdet_pj(4,ngdim)

      real g0inv(2,2,ngdim)
      real ginv(2,2,ngdim),
     &     ginv_rj(2,2,3,4,ngdim), ginv_dj(2,2,3,4,ngdim),
     &     ginv_pj(2,2,4,ngdim)

      real h0(2,2,ngdim)
      real h(2,2,ngdim),
     &     h_rj(2,2,3,4,ngdim), h_dj(2,2,3,4,ngdim), h_pj(2,2,4,ngdim)

      real l0(2,ngdim)
      real l(2,ngdim),
     &     l_rj(2,3,4,ngdim), l_dj(2,3,4,ngdim), l_pj(2,4,ngdim)

      real rxr0(3,ngdim)
      real rxr(3,ngdim),
     &     rxr_rj(3,3,4,ngdim)
      real arxr_rj(3,4)

c- - - - - - - - - - - - - - - - - - - - - - - - - - 
c---- quantities for evaluation of residual integrands
      real eps(2,2), eps_rj(2,2,3,4), eps_dj(2,2,3,4), eps_pj(2,2,4)
      real kap(2,2), kap_rj(2,2,3,4), kap_dj(2,2,3,4), kap_pj(2,2,4)
      
      real aac(2,2,2,2), bbc(2,2,2,2),  ddc(2,2,2,2)
      real ssc(2,2)

      real aa(2,2,2,2), bb(2,2,2,2),  dd(2,2,2,2)
      real ss(2,2)

      real fab(2,2), fab_rj(2,2,3,4), fab_dj(2,2,3,4), fab_pj(2,2,4)
      real mab(2,2), mab_rj(2,2,3,4), mab_dj(2,2,3,4), mab_pj(2,2,4)

      real fan(2), fan_rj(2,3,4), fan_dj(2,3,4), fan_pj(2,4)
      real gan(2), gan_rj(2,3,4), gan_dj(2,3,4), gan_pj(2,4)

      real nn(4)
      real d1nn(4), d2nn(4)

      real ajac0
      real ajac, ajac_rj(3,4), ajac_dj(3,4), ajac_pj(4)

      real eepse(2,2),
     &     ekape(2,2),
     &     efore(2,2),
     &     emome(2,2)
      real ete_tj(2,2,2,2,4)

      real efn(2),
     &     egn(2)
      real ev_vj(2,2,4)

      real eepser(2,2),
     &     ekaper(2,2),
     &     eforer(2,2),
     &     emomer(2,2) 

      real efnr(2),
     &     egnr(2)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------

      if(neln .eq. 3) then
c----- use tri element quadrature
       ngauss = ngausst
       do ig = 1, ngauss
         xsig(ig) = xsit(ig)
         etag(ig) = etat(ig)
         fwtg(ig) = fwtt(ig)
       enddo

      else
c----- use quad element quadrature
       ngauss = ngaussq
       do ig = 1, ngauss
         xsig(ig) = xsiq(ig)
         etag(ig) = etaq(ig)
         fwtg(ig) = fwtq(ig)
       enddo

      endif

c---- put passed node quantities into more convenient node-indexed arrays
      do lv = 1, lvtot
        par(lv,1) = par1(lv)
        par(lv,2) = par2(lv)
        par(lv,3) = par3(lv)
        par(lv,4) = par4(lv)
      enddo

      do iv = 1, ivtot
        var(iv,1) = var1(iv)
        var(iv,2) = var2(iv)
        var(iv,3) = var3(iv)
        var(iv,4) = var4(iv)
      enddo

      do jv = 1, jvtot
        dep(jv,1) = dep1(jv)
        dep(jv,2) = dep2(jv)
        dep(jv,3) = dep3(jv)
        dep(jv,4) = dep4(jv)
      enddo

c---- further unpack vectors
      do j = 1, 4
        r0j(1,j) = par(lvr0x,j)
        r0j(2,j) = par(lvr0y,j)
        r0j(3,j) = par(lvr0z,j)

        d0j(1,j) = par(lvn0x,j)
        d0j(2,j) = par(lvn0y,j)
        d0j(3,j) = par(lvn0z,j)

        p0j(j) = 0.

        e0j(1,1,j) = par(lve01x,j)
        e0j(2,1,j) = par(lve01y,j)
        e0j(3,1,j) = par(lve01z,j)

        e0j(1,2,j) = par(lve02x,j)
        e0j(2,2,j) = par(lve02y,j)
        e0j(3,2,j) = par(lve02z,j)

        rj(1,j) = var(ivrx,j)
        rj(2,j) = var(ivry,j)
        rj(3,j) = var(ivrz,j)

        dj(1,j) = var(ivdx,j)
        dj(2,j) = var(ivdy,j)
        dj(3,j) = var(ivdz,j)

        pj(j) = var(ivps,j)

        epsj(1,1,j) = dep(jve11,j)
        epsj(2,2,j) = dep(jve22,j)
        epsj(1,2,j) = dep(jve12,j)
        epsj(2,1,j) = dep(jve12,j)

        kapj(1,1,j) = dep(jvk11,j)
        kapj(2,2,j) = dep(jvk22,j)
        kapj(1,2,j) = dep(jvk12,j)
        kapj(2,1,j) = dep(jvk12,j)

        fj(1,1,j) = dep(jvf11,j)
        fj(2,2,j) = dep(jvf22,j)
        fj(1,2,j) = dep(jvf12,j)
        fj(2,1,j) = dep(jvf12,j)

        mj(1,1,j) = dep(jvm11,j)
        mj(2,2,j) = dep(jvm22,j)
        mj(1,2,j) = dep(jvm12,j)
        mj(2,1,j) = dep(jvm12,j)

        fnj(1,j) = dep(jvf1n,j)
        fnj(2,j) = dep(jvf2n,j)

        gnj(1,j) = dep(jvga1,j)
        gnj(2,j) = dep(jvga2,j)
      enddo

c---- element cartesian basis vectors
      if(neln .eq. 3) then
       xsi = 0.33333333333333333333333333333
       eta = 0.33333333333333333333333333333

       nn(1) = 1.0 - xsi - eta
       nn(2) = xsi
       nn(3) = eta
       nn(4) = 0.

       d1nn(1) = -1.0
       d1nn(2) =  1.0
       d1nn(3) =  0.
       d1nn(4) = 0.

       d2nn(1) = -1.0
       d2nn(2) =  0.
       d2nn(3) =  1.0
       d2nn(4) = 0.

      else
       xsi = 0.0
       eta = 0.0

       nn(1) = 0.25*(1.0 - xsi)*(1.0 - eta)
       nn(2) = 0.25*(1.0 + xsi)*(1.0 - eta)
       nn(3) = 0.25*(1.0 + xsi)*(1.0 + eta)
       nn(4) = 0.25*(1.0 - xsi)*(1.0 + eta)

       d1nn(1) = -0.25*(1.0 - eta)
       d1nn(2) =  0.25*(1.0 - eta)
       d1nn(3) =  0.25*(1.0 + eta)
       d1nn(4) = -0.25*(1.0 + eta)

       d2nn(1) = -0.25*(1.0 - xsi)
       d2nn(2) = -0.25*(1.0 + xsi)
       d2nn(3) =  0.25*(1.0 + xsi)
       d2nn(4) =  0.25*(1.0 - xsi)

      endif

c---- centroid dr0/dxsi, dr0/deta
      do k = 1, 3
        d1r0(k) = r0j(k,1)*d1nn(1)
     &          + r0j(k,2)*d1nn(2)
     &          + r0j(k,3)*d1nn(3)
     &          + r0j(k,4)*d1nn(4)
        d2r0(k) = r0j(k,1)*d2nn(1)
     &          + r0j(k,2)*d2nn(2)
     &          + r0j(k,3)*d2nn(3)
     &          + r0j(k,4)*d2nn(4)
      enddo

c---- element basis vectors
      d1ra = sqrt(d1r0(1)**2 + d1r0(2)**2 + d1r0(3)**2)
      c1(1) = d1r0(1)/d1ra
      c1(2) = d1r0(2)/d1ra
      c1(3) = d1r0(3)/d1ra

      crs(1) = d1r0(2)*d2r0(3) - d1r0(3)*d2r0(2)
      crs(2) = d1r0(3)*d2r0(1) - d1r0(1)*d2r0(3)
      crs(3) = d1r0(1)*d2r0(2) - d1r0(2)*d2r0(1)
      crsa = sqrt(crs(1)**2 + crs(2)**2 + crs(3)**2)
      c3(1) = crs(1)/crsa
      c3(2) = crs(2)/crsa
      c3(3) = crs(3)/crsa

      c2(1) = c3(2)*c1(3) - c3(3)*c1(2)
      c2(2) = c3(3)*c1(1) - c3(1)*c1(3)
      c2(3) = c3(1)*c1(2) - c3(2)*c1(1)

c---- put node stiffness matrices into element basis at each node
      call hsmabd(neln, par1    , par2    , par3    , par4    , 
     &            aaj, bbj, ddj, ssj )

c-------------------------------------------------------------------------
c---- set element geometry quantities at Gauss point
      call hsmgeo(neln, ngauss, xsig, etag, fwtg,
     &            lrcurv, ldrill,
     &            r0j(1,1), r0j(1,2), r0j(1,3), r0j(1,4),
     &            d0j(1,1), d0j(1,2), d0j(1,3), d0j(1,4),
     &            p0j(1)  , p0j(2)  , p0j(3)  , p0j(4)  ,
     &             r0,  r_rj,  r_dj,  r_pj,
     &             d0,  d_rj,  d_dj,  d_pj,
     &            a01, a1_rj, a1_dj, a1_pj,
     &            a02, a2_rj, a2_dj, a2_pj,
     &            a01i, a1i_rj, a1i_dj, a1i_pj,
     &            a02i, a2i_rj, a2i_dj, a2i_pj,
     &            g0  ,   g_rj,  g_dj, g_pj,
     &            h0  ,   h_rj,  h_dj, h_pj,
     &            l0  ,   l_rj,  l_dj, l_pj,
     &            g0inv, ginv_rj, ginv_dj, ginv_pj,
     &            g0det, gdet_rj, gdet_dj, gdet_pj,
     &            rxr0 ,  rxr_rj )

      call hsmgeo(neln, ngauss, xsig, etag, fwtg,
     &            lrcurv, ldrill,
     &            rj(1,1), rj(1,2), rj(1,3), rj(1,4),
     &            dj(1,1), dj(1,2), dj(1,3), dj(1,4),
     &            pj(1)  , pj(2)  , pj(3)  , pj(4)  ,
     &             r,  r_rj,  r_dj,  r_pj,
     &             d,  d_rj,  d_dj,  d_pj,
     &            a1, a1_rj, a1_dj, a1_pj,
     &            a2, a2_rj, a2_dj, a2_pj,
     &            a1i, a1i_rj, a1i_dj, a1i_pj,
     &            a2i, a2i_rj, a2i_dj, a2i_pj,
     &            g  ,   g_rj,  g_dj, g_pj,
     &            h  ,   h_rj,  h_dj, h_pj,
     &            l  ,   l_rj,  l_dj, l_pj,
     &            ginv, ginv_rj, ginv_dj, ginv_pj,
     &            gdet, gdet_rj, gdet_dj, gdet_pj,
     &            rxr ,  rxr_rj )

c=======================================================================
c---- clear residuals and Jacobians for summing over gauss points
      do ires = 1, 4
        do keq = 1, 3
          reseps(keq,ires) = 0.
          reskap(keq,ires) = 0.
          resfor(keq,ires) = 0.
          resmom(keq,ires) = 0.
          do jj = 1, 4
            do kv = 1, 3
              rest_t(keq,ires,kv,jj) = 0.
            enddo
          enddo
        enddo

        do keq = 1, 2
          resfn(keq,ires) = 0.
          resgn(keq,ires) = 0.
          do jj = 1, 4
            do kv = 1, 2
              resv_v(keq,ires,kv,jj) = 0.
            enddo
          enddo
        enddo

        ac(ires) = 0.
      enddo

c-------------------------------------------------------------------
c---- sum over gauss points
      do 500 ig = 1, ngauss

c---- bilinear interpolation weights "N" and their derivatives at gauss point
      if(neln .eq. 3) then
       xsi = xsit(ig)
       eta = etat(ig)
       fwt = fwtt(ig)

       nn(1) = 1.0 - xsi - eta
       nn(2) = xsi
       nn(3) = eta
       nn(4) = 0.

       d1nn(1) = -1.0
       d1nn(2) =  1.0
       d1nn(3) =  0.
       d1nn(4) = 0.

       d2nn(1) = -1.0
       d2nn(2) =  0.
       d2nn(3) =  1.0
       d2nn(4) = 0.

      else
       xsi = xsiq(ig)
       eta = etaq(ig)
       fwt = fwtq(ig)

       nn(1) = 0.25*(1.0 - xsi)*(1.0 - eta)
       nn(2) = 0.25*(1.0 + xsi)*(1.0 - eta)
       nn(3) = 0.25*(1.0 + xsi)*(1.0 + eta)
       nn(4) = 0.25*(1.0 - xsi)*(1.0 + eta)

       d1nn(1) = -0.25*(1.0 - eta)
       d1nn(2) =  0.25*(1.0 - eta)
       d1nn(3) =  0.25*(1.0 + eta)
       d1nn(4) = -0.25*(1.0 + eta)

       d2nn(1) = -0.25*(1.0 - xsi)
       d2nn(2) = -0.25*(1.0 + xsi)
       d2nn(3) =  0.25*(1.0 + xsi)
       d2nn(4) =  0.25*(1.0 - xsi)

      endif


      arxr0 = sqrt(  rxr0(1,ig)**2
     &             + rxr0(2,ig)**2
     &             + rxr0(3,ig)**2 )

      arxr  = sqrt(  rxr(1,ig)**2
     &             + rxr(2,ig)**2
     &             + rxr(3,ig)**2 )
      if(arxr .eq. 0.0) then
       arxri = 1.0
      else
       arxri = 1.0/arxr
      endif

      do jj = 1, neln
        do kk = 1, 3
          arxr_rj(kk,jj) = 
     &      (  rxr_rj(1,kk,jj,ig)*rxr(1,ig)
     &       + rxr_rj(2,kk,jj,ig)*rxr(2,ig)
     &       + rxr_rj(3,kk,jj,ig)*rxr(3,ig) ) * arxri
        enddo
      enddo

c---- covariant Green strain tensor
      eps(1,1) = (g(1,1,ig) - g0(1,1,ig))*0.5
      eps(1,2) = (g(1,2,ig) - g0(1,2,ig))*0.5
      eps(2,1) = (g(2,1,ig) - g0(2,1,ig))*0.5
      eps(2,2) = (g(2,2,ig) - g0(2,2,ig))*0.5

c---- covariant curvature-change tensor
      kap(1,1) =  h(1,1,ig) - h0(1,1,ig)
      kap(1,2) =  h(1,2,ig) - h0(1,2,ig)
      kap(2,1) =  h(2,1,ig) - h0(2,1,ig)
      kap(2,2) =  h(2,2,ig) - h0(2,2,ig)

c---- covariant transverse shear strain vector
      gan(1) = l(1,ig) - l0(1,ig)
      gan(2) = l(2,ig) - l0(2,ig)

      do jj = 1, neln
        do kk = 1, 3
          eps_rj(1,1,kk,jj) = g_rj(1,1,kk,jj,ig)*0.5
          eps_rj(1,2,kk,jj) = g_rj(1,2,kk,jj,ig)*0.5
          eps_rj(2,1,kk,jj) = g_rj(2,1,kk,jj,ig)*0.5
          eps_rj(2,2,kk,jj) = g_rj(2,2,kk,jj,ig)*0.5

          eps_dj(1,1,kk,jj) = g_dj(1,1,kk,jj,ig)*0.5
          eps_dj(1,2,kk,jj) = g_dj(1,2,kk,jj,ig)*0.5
          eps_dj(2,1,kk,jj) = g_dj(2,1,kk,jj,ig)*0.5
          eps_dj(2,2,kk,jj) = g_dj(2,2,kk,jj,ig)*0.5

          kap_rj(1,1,kk,jj) = h_rj(1,1,kk,jj,ig)
          kap_rj(1,2,kk,jj) = h_rj(1,2,kk,jj,ig)
          kap_rj(2,1,kk,jj) = h_rj(2,1,kk,jj,ig)
          kap_rj(2,2,kk,jj) = h_rj(2,2,kk,jj,ig)

          kap_dj(1,1,kk,jj) = h_dj(1,1,kk,jj,ig)
          kap_dj(1,2,kk,jj) = h_dj(1,2,kk,jj,ig)
          kap_dj(2,1,kk,jj) = h_dj(2,1,kk,jj,ig)
          kap_dj(2,2,kk,jj) = h_dj(2,2,kk,jj,ig)

          gan_rj(1,kk,jj) = l_rj(1,kk,jj,ig)
          gan_rj(2,kk,jj) = l_rj(2,kk,jj,ig)

          gan_dj(1,kk,jj) = l_dj(1,kk,jj,ig)
          gan_dj(2,kk,jj) = l_dj(2,kk,jj,ig)
        enddo
        eps_pj(1,1,jj) = g_pj(1,1,jj,ig)*0.5
        eps_pj(1,2,jj) = g_pj(1,2,jj,ig)*0.5
        eps_pj(2,1,jj) = g_pj(2,1,jj,ig)*0.5
        eps_pj(2,2,jj) = g_pj(2,2,jj,ig)*0.5

        kap_pj(1,1,jj) = h_pj(1,1,jj,ig)
        kap_pj(1,2,jj) = h_pj(1,2,jj,ig)
        kap_pj(2,1,jj) = h_pj(2,1,jj,ig)
        kap_pj(2,2,jj) = h_pj(2,2,jj,ig)

        gan_pj(1,jj) = l_pj(1,jj,ig)
        gan_pj(2,jj) = l_pj(2,jj,ig)
      enddo

c---- interpolate stiffness matrices to Gauss point
      do id = 1, 2
        do ic = 1, 2
          do ib = 1, 2
            do ia = 1, 2
              asum = 0.
              bsum = 0.
              dsum = 0.
              do j = 1, neln
                asum = asum + aaj(ia,ib,ic,id,j)*nn(j)
                bsum = bsum + bbj(ia,ib,ic,id,j)*nn(j)
                dsum = dsum + ddj(ia,ib,ic,id,j)*nn(j)
              enddo
              aac(ia,ib,ic,id) = asum
              bbc(ia,ib,ic,id) = bsum
              ddc(ia,ib,ic,id) = dsum
            enddo
          enddo
        enddo
      enddo

      do ib = 1, 2
        do ia = 1, 2
          ssum = 0.
          do j = 1, neln
            ssum = ssum + ssj(ia,ib,j)*nn(j)
          enddo
          ssc(ia,ib) = ssum
        enddo
      enddo

c---- basis vector dot products
      cai(1,1) = c1(1)*a01i(1,ig)
     &         + c1(2)*a01i(2,ig)
     &         + c1(3)*a01i(3,ig)
      cai(1,2) = c1(1)*a02i(1,ig)
     &         + c1(2)*a02i(2,ig)
     &         + c1(3)*a02i(3,ig)
      cai(2,1) = c2(1)*a01i(1,ig)
     &         + c2(2)*a01i(2,ig)
     &         + c2(3)*a01i(3,ig)
      cai(2,2) = c2(1)*a02i(1,ig)
     &         + c2(2)*a02i(2,ig)
     &         + c2(3)*a02i(3,ig)

c---- put stiffness matrices into element basis
      do id = 1, 2
        do ic = 1, 2
          do ib = 1, 2
            do ia = 1, 2
              aa(ia,ib,ic,id) =
     &              aac(1,1,1,1)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + aac(1,1,1,2)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + aac(1,1,2,1)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + aac(1,1,2,2)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + aac(1,2,1,1)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + aac(1,2,1,2)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + aac(1,2,2,1)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + aac(1,2,2,2)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
     &            + aac(2,1,1,1)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + aac(2,1,1,2)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + aac(2,1,2,1)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + aac(2,1,2,2)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + aac(2,2,1,1)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + aac(2,2,1,2)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + aac(2,2,2,1)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + aac(2,2,2,2)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
              bb(ia,ib,ic,id) =                                     
     &              bbc(1,1,1,1)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + bbc(1,1,1,2)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + bbc(1,1,2,1)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + bbc(1,1,2,2)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + bbc(1,2,1,1)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + bbc(1,2,1,2)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + bbc(1,2,2,1)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + bbc(1,2,2,2)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
     &            + bbc(2,1,1,1)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + bbc(2,1,1,2)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + bbc(2,1,2,1)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + bbc(2,1,2,2)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + bbc(2,2,1,1)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + bbc(2,2,1,2)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + bbc(2,2,2,1)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + bbc(2,2,2,2)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
              dd(ia,ib,ic,id) =                                     
     &              ddc(1,1,1,1)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + ddc(1,1,1,2)*cai(1,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + ddc(1,1,2,1)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + ddc(1,1,2,2)*cai(1,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + ddc(1,2,1,1)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + ddc(1,2,1,2)*cai(1,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + ddc(1,2,2,1)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + ddc(1,2,2,2)*cai(1,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
     &            + ddc(2,1,1,1)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(1,id)
     &            + ddc(2,1,1,2)*cai(2,ia)*cai(1,ib)*cai(1,ic)*cai(2,id)
     &            + ddc(2,1,2,1)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(1,id)
     &            + ddc(2,1,2,2)*cai(2,ia)*cai(1,ib)*cai(2,ic)*cai(2,id)
     &            + ddc(2,2,1,1)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(1,id)
     &            + ddc(2,2,1,2)*cai(2,ia)*cai(2,ib)*cai(1,ic)*cai(2,id)
     &            + ddc(2,2,2,1)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(1,id)
     &            + ddc(2,2,2,2)*cai(2,ia)*cai(2,ib)*cai(2,ic)*cai(2,id)
            enddo
          enddo
        enddo
      enddo

      do ib = 1, 2
        do ia = 1, 2
          ss(ia,ib) =
     &        + ssc(1,1)*cai(1,ia)*cai(1,ib)
     &        + ssc(1,2)*cai(1,ia)*cai(2,ib)
     &        + ssc(2,1)*cai(2,ia)*cai(1,ib)
     &        + ssc(2,2)*cai(2,ia)*cai(2,ib)
        enddo
      enddo

c---- in-surface stress and stress-moment resultants in element basis
      do ib = 1, 2
        do ia = 1, 2
          fab(ia,ib) = 
     &                 aa(ia,ib,1,1)*eps(1,1)
     &               + aa(ia,ib,1,2)*eps(1,2)
     &               + aa(ia,ib,2,1)*eps(2,1)
     &               + aa(ia,ib,2,2)*eps(2,2)
     &               + bb(ia,ib,1,1)*kap(1,1)
     &               + bb(ia,ib,1,2)*kap(1,2)
     &               + bb(ia,ib,2,1)*kap(2,1)
     &               + bb(ia,ib,2,2)*kap(2,2)
          mab(ia,ib) = 
     &                 bb(ia,ib,1,1)*eps(1,1)
     &               + bb(ia,ib,1,2)*eps(1,2)
     &               + bb(ia,ib,2,1)*eps(2,1)
     &               + bb(ia,ib,2,2)*eps(2,2)
     &               + dd(ia,ib,1,1)*kap(1,1)
     &               + dd(ia,ib,1,2)*kap(1,2)
     &               + dd(ia,ib,2,1)*kap(2,1)
     &               + dd(ia,ib,2,2)*kap(2,2)
          do jj = 1, neln
            do kk = 1, 3
              fab_rj(ia,ib,kk,jj) = 
     &                 aa(ia,ib,1,1)*eps_rj(1,1,kk,jj)
     &               + aa(ia,ib,1,2)*eps_rj(1,2,kk,jj)
     &               + aa(ia,ib,2,1)*eps_rj(2,1,kk,jj)
     &               + aa(ia,ib,2,2)*eps_rj(2,2,kk,jj)
     &               + bb(ia,ib,1,1)*kap_rj(1,1,kk,jj)
     &               + bb(ia,ib,1,2)*kap_rj(1,2,kk,jj)
     &               + bb(ia,ib,2,1)*kap_rj(2,1,kk,jj)
     &               + bb(ia,ib,2,2)*kap_rj(2,2,kk,jj)
              mab_rj(ia,ib,kk,jj) = 
     &                 bb(ia,ib,1,1)*eps_rj(1,1,kk,jj)
     &               + bb(ia,ib,1,2)*eps_rj(1,2,kk,jj)
     &               + bb(ia,ib,2,1)*eps_rj(2,1,kk,jj)
     &               + bb(ia,ib,2,2)*eps_rj(2,2,kk,jj)
     &               + dd(ia,ib,1,1)*kap_rj(1,1,kk,jj)
     &               + dd(ia,ib,1,2)*kap_rj(1,2,kk,jj)
     &               + dd(ia,ib,2,1)*kap_rj(2,1,kk,jj)
     &               + dd(ia,ib,2,2)*kap_rj(2,2,kk,jj)

              fab_dj(ia,ib,kk,jj) = 
     &                 aa(ia,ib,1,1)*eps_dj(1,1,kk,jj)
     &               + aa(ia,ib,1,2)*eps_dj(1,2,kk,jj)
     &               + aa(ia,ib,2,1)*eps_dj(2,1,kk,jj)
     &               + aa(ia,ib,2,2)*eps_dj(2,2,kk,jj)
     &               + bb(ia,ib,1,1)*kap_dj(1,1,kk,jj)
     &               + bb(ia,ib,1,2)*kap_dj(1,2,kk,jj)
     &               + bb(ia,ib,2,1)*kap_dj(2,1,kk,jj)
     &               + bb(ia,ib,2,2)*kap_dj(2,2,kk,jj)
              mab_dj(ia,ib,kk,jj) = 
     &                 bb(ia,ib,1,1)*eps_dj(1,1,kk,jj)
     &               + bb(ia,ib,1,2)*eps_dj(1,2,kk,jj)
     &               + bb(ia,ib,2,1)*eps_dj(2,1,kk,jj)
     &               + bb(ia,ib,2,2)*eps_dj(2,2,kk,jj)
     &               + dd(ia,ib,1,1)*kap_dj(1,1,kk,jj)
     &               + dd(ia,ib,1,2)*kap_dj(1,2,kk,jj)
     &               + dd(ia,ib,2,1)*kap_dj(2,1,kk,jj)
     &               + dd(ia,ib,2,2)*kap_dj(2,2,kk,jj)
            enddo
            fab_pj(ia,ib,jj) = 
     &               aa(ia,ib,1,1)*eps_pj(1,1,jj)
     &             + aa(ia,ib,1,2)*eps_pj(1,2,jj)
     &             + aa(ia,ib,2,1)*eps_pj(2,1,jj)
     &             + aa(ia,ib,2,2)*eps_pj(2,2,jj)
     &             + bb(ia,ib,1,1)*kap_pj(1,1,jj)
     &             + bb(ia,ib,1,2)*kap_pj(1,2,jj)
     &             + bb(ia,ib,2,1)*kap_pj(2,1,jj)
     &             + bb(ia,ib,2,2)*kap_pj(2,2,jj)
            mab_pj(ia,ib,jj) = 
     &               bb(ia,ib,1,1)*eps_pj(1,1,jj)
     &             + bb(ia,ib,1,2)*eps_pj(1,2,jj)
     &             + bb(ia,ib,2,1)*eps_pj(2,1,jj)
     &             + bb(ia,ib,2,2)*eps_pj(2,2,jj)
     &             + dd(ia,ib,1,1)*kap_pj(1,1,jj)
     &             + dd(ia,ib,1,2)*kap_pj(1,2,jj)
     &             + dd(ia,ib,2,1)*kap_pj(2,1,jj)
     &             + dd(ia,ib,2,2)*kap_pj(2,2,jj)
          enddo
        enddo
      enddo

c---- transverse shear stress in element basis, fn^a
      fan(1) = ss(1,1)*gan(1)
     &       + ss(1,2)*gan(2)
      fan(2) = ss(2,1)*gan(1)
     &       + ss(2,2)*gan(2)
      do jj = 1, neln
        do kk = 1, 3
          fan_rj(1,kk,jj) = ss(1,1)*gan_rj(1,kk,jj)
     &                    + ss(1,2)*gan_rj(2,kk,jj)
          fan_rj(2,kk,jj) = ss(2,1)*gan_rj(1,kk,jj)
     &                    + ss(2,2)*gan_rj(2,kk,jj)

          fan_dj(1,kk,jj) = ss(1,1)*gan_dj(1,kk,jj)
     &                    + ss(1,2)*gan_dj(2,kk,jj)
          fan_dj(2,kk,jj) = ss(2,1)*gan_dj(1,kk,jj)
     &                    + ss(2,2)*gan_dj(2,kk,jj)
        enddo

        fan_pj(1,jj) = ss(1,1)*gan_pj(1,jj)
     &               + ss(1,2)*gan_pj(2,jj)
        fan_pj(2,jj) = ss(2,1)*gan_pj(1,jj)
     &               + ss(2,2)*gan_pj(2,jj)
      enddo

c---- coordinate Jacobian J
      ajac0 = arxr0
      ajac = arxr
      do jj = 1, neln
        do kk = 1, 3
          ajac_rj(kk,jj) = arxr_rj(kk,jj)
          ajac_dj(kk,jj) = 0.
        enddo
        ajac_pj(jj) = 0.
      enddo

c---- residual area-weight for this gauss point
      awt = fwt * ajac0

c--------------------------------------------------
c---- go over the four residual weight functions
      do 100 ires = 1, neln

c---- weight function  Wi = Ni
      ww = nn(ires)
cc    ww = 1.0 + 0.125*nn(ires)

c---- tensor dot products with e1i,e2i at this Gauss point
      do ia = 1, 2
        do ib = 1, 2

c-------- tensors interpolated from nodal values
          eepse(ia,ib) = 0.
          ekape(ia,ib) = 0.
          efore(ia,ib) = 0.
          emome(ia,ib) = 0.

          do jj = 1, neln
            eke1 = e0j(1,ia,ires)*e0j(1,1,jj)
     &           + e0j(2,ia,ires)*e0j(2,1,jj)
     &           + e0j(3,ia,ires)*e0j(3,1,jj)
            eke2 = e0j(1,ia,ires)*e0j(1,2,jj)
     &           + e0j(2,ia,ires)*e0j(2,2,jj)
     &           + e0j(3,ia,ires)*e0j(3,2,jj)
            ele1 = e0j(1,ib,ires)*e0j(1,1,jj)
     &           + e0j(2,ib,ires)*e0j(2,1,jj)
     &           + e0j(3,ib,ires)*e0j(3,1,jj)
            ele2 = e0j(1,ib,ires)*e0j(1,2,jj)
     &           + e0j(2,ib,ires)*e0j(2,2,jj)
     &           + e0j(3,ib,ires)*e0j(3,2,jj)

            eepse(ia,ib) = eepse(ia,ib)
     &              + ( epsj(1,1,jj)*eke1*ele1
     &                + epsj(1,2,jj)*eke1*ele2
     &                + epsj(2,1,jj)*eke2*ele1
     &                + epsj(2,2,jj)*eke2*ele2 ) * nn(jj)

            ekape(ia,ib) = ekape(ia,ib)
     &              + ( kapj(1,1,jj)*eke1*ele1
     &                + kapj(1,2,jj)*eke1*ele2
     &                + kapj(2,1,jj)*eke2*ele1
     &                + kapj(2,2,jj)*eke2*ele2 ) * nn(jj)

            efore(ia,ib) = efore(ia,ib)
     &              + ( fj(1,1,jj)*eke1*ele1
     &                + fj(1,2,jj)*eke1*ele2
     &                + fj(2,1,jj)*eke2*ele1
     &                + fj(2,2,jj)*eke2*ele2 ) * nn(jj)

            emome(ia,ib) = emome(ia,ib)
     &              + ( mj(1,1,jj)*eke1*ele1
     &                + mj(1,2,jj)*eke1*ele2
     &                + mj(2,1,jj)*eke2*ele1
     &                + mj(2,2,jj)*eke2*ele2 ) * nn(jj)

            ete_tj(ia,ib,1,1,jj) = eke1*ele1 * nn(jj)
            ete_tj(ia,ib,1,2,jj) = eke1*ele2 * nn(jj)
            ete_tj(ia,ib,2,1,jj) = eke2*ele1 * nn(jj)
            ete_tj(ia,ib,2,2,jj) = eke2*ele2 * nn(jj)
          enddo

c-------- membrane and bending strains from geometry
          eka1i = e0j(1,ia,ires)*a01i(1,ig)
     &          + e0j(2,ia,ires)*a01i(2,ig)
     &          + e0j(3,ia,ires)*a01i(3,ig)
          eka2i = e0j(1,ia,ires)*a02i(1,ig)
     &          + e0j(2,ia,ires)*a02i(2,ig)
     &          + e0j(3,ia,ires)*a02i(3,ig)

          ela1i = e0j(1,ib,ires)*a01i(1,ig)
     &          + e0j(2,ib,ires)*a01i(2,ig)
     &          + e0j(3,ib,ires)*a01i(3,ig)
          ela2i = e0j(1,ib,ires)*a02i(1,ig)
     &          + e0j(2,ib,ires)*a02i(2,ig)
     &          + e0j(3,ib,ires)*a02i(3,ig)

          eepser(ia,ib) = eps(1,1)*eka1i*ela1i
     &                  + eps(1,2)*eka1i*ela2i
     &                  + eps(2,1)*eka2i*ela1i
     &                  + eps(2,2)*eka2i*ela2i

          ekaper(ia,ib) = kap(1,1)*eka1i*ela1i
     &                  + kap(1,2)*eka1i*ela2i
     &                  + kap(2,1)*eka2i*ela1i
     &                  + kap(2,2)*eka2i*ela2i


c-------- force and moment from geometry strains and stiffnesses
          eka1 = e0j(1,ia,ires)*a01(1,ig)
     &         + e0j(2,ia,ires)*a01(2,ig)
     &         + e0j(3,ia,ires)*a01(3,ig)
          eka2 = e0j(1,ia,ires)*a02(1,ig)
     &         + e0j(2,ia,ires)*a02(2,ig)
     &         + e0j(3,ia,ires)*a02(3,ig)

          ela1 = e0j(1,ib,ires)*a01(1,ig)
     &         + e0j(2,ib,ires)*a01(2,ig)
     &         + e0j(3,ib,ires)*a01(3,ig)
          ela2 = e0j(1,ib,ires)*a02(1,ig)
     &         + e0j(2,ib,ires)*a02(2,ig)
     &         + e0j(3,ib,ires)*a02(3,ig)

      
          eforer(ia,ib) = fab(1,1)*eka1*ela1
     &                  + fab(1,2)*eka1*ela2
     &                  + fab(2,1)*eka2*ela1
     &                  + fab(2,2)*eka2*ela2

          emomer(ia,ib) = mab(1,1)*eka1*ela1
     &                  + mab(1,2)*eka1*ela2
     &                  + mab(2,1)*eka2*ela1
     &                  + mab(2,2)*eka2*ela2
        enddo
      enddo

      do ia = 1, 2
c------ transverse shear and strain interpolated from nodal values
        efn(ia) = 0.
        egn(ia) = 0.

        do jj = 1, neln
          eke1 = e0j(1,ia,ires)*e0j(1,1,jj)
     &         + e0j(2,ia,ires)*e0j(2,1,jj)
     &         + e0j(3,ia,ires)*e0j(3,1,jj)
          eke2 = e0j(1,ia,ires)*e0j(1,2,jj)
     &         + e0j(2,ia,ires)*e0j(2,2,jj)
     &         + e0j(3,ia,ires)*e0j(3,2,jj)

          efn(ia) = efn(ia)
     &            + ( fnj(1,jj)*eke1
     &              + fnj(2,jj)*eke2 ) * nn(jj)
          egn(ia) = egn(ia)
     &            + ( gnj(1,jj)*eke1
     &              + gnj(2,jj)*eke2 ) * nn(jj)

          ev_vj(ia,1,jj) = eke1 * nn(jj)
          ev_vj(ia,2,jj) = eke2 * nn(jj)
        enddo

c------ transverse shear and strain calculated from stiffnesses and geometry
        eka1 = e0j(1,ia,ires)*a01(1,ig)
     &       + e0j(2,ia,ires)*a01(2,ig)
     &       + e0j(3,ia,ires)*a01(3,ig)
        eka2 = e0j(1,ia,ires)*a02(1,ig)
     &       + e0j(2,ia,ires)*a02(2,ig)
     &       + e0j(3,ia,ires)*a02(3,ig)

        eka1i = e0j(1,ia,ires)*a01i(1,ig)
     &        + e0j(2,ia,ires)*a01i(2,ig)
     &        + e0j(3,ia,ires)*a01i(3,ig)
        eka2i = e0j(1,ia,ires)*a02i(1,ig)
     &        + e0j(2,ia,ires)*a02i(2,ig)
     &        + e0j(3,ia,ires)*a02i(3,ig)

        efnr(ia) = fan(1)*eka1
     &           + fan(2)*eka2
        egnr(ia) = gan(1)*eka1i
     &           + gan(2)*eka2i
      enddo

c----------------------------------------------------------------------
c---- T11 residuals
      keq = 1
      ia = 1
      ib = 1
       reseps(keq,ires) = 
     & reseps(keq,ires) + (eepse(ia,ib) - eepser(ia,ib))*ww * awt
       reskap(keq,ires) = 
     & reskap(keq,ires) + (ekape(ia,ib) - ekaper(ia,ib))*ww * awt
       resfor(keq,ires) = 
     & resfor(keq,ires) + (efore(ia,ib) - eforer(ia,ib))*ww * awt
       resmom(keq,ires) = 
     & resmom(keq,ires) + (emome(ia,ib) - emomer(ia,ib))*ww * awt

      do jj = 1, neln
        rest_t(keq,ires,1,jj) = 
     &  rest_t(keq,ires,1,jj) + ete_tj(ia,ib,1,1,jj)*ww * awt
        rest_t(keq,ires,2,jj) = 
     &  rest_t(keq,ires,2,jj) + ete_tj(ia,ib,2,2,jj)*ww * awt
        rest_t(keq,ires,3,jj) = 
     &  rest_t(keq,ires,3,jj) + ete_tj(ia,ib,1,2,jj)*ww * awt
     &                        + ete_tj(ia,ib,2,1,jj)*ww * awt
      enddo

c----------------------------------------------------------------------
c---- T22 residuals
      keq = 2
      ia = 2
      ib = 2
       reseps(keq,ires) =
     & reseps(keq,ires) + (eepse(ia,ib) - eepser(ia,ib))*ww * awt
       reskap(keq,ires) =
     & reskap(keq,ires) + (ekape(ia,ib) - ekaper(ia,ib))*ww * awt
       resfor(keq,ires) =
     & resfor(keq,ires) + (efore(ia,ib) - eforer(ia,ib))*ww * awt
       resmom(keq,ires) =
     & resmom(keq,ires) + (emome(ia,ib) - emomer(ia,ib))*ww * awt

      do jj = 1, neln
        rest_t(keq,ires,1,jj) = 
     &  rest_t(keq,ires,1,jj) + ete_tj(ia,ib,1,1,jj)*ww * awt
        rest_t(keq,ires,2,jj) = 
     &  rest_t(keq,ires,2,jj) + ete_tj(ia,ib,2,2,jj)*ww * awt
        rest_t(keq,ires,3,jj) = 
     &  rest_t(keq,ires,3,jj) + ete_tj(ia,ib,1,2,jj)*ww * awt
     &                        + ete_tj(ia,ib,2,1,jj)*ww * awt
      enddo

c----------------------------------------------------------------------
c---- T12 residuals
      keq = 3
      ia = 2
      ib = 1
       reseps(keq,ires) =
     & reseps(keq,ires) + (eepse(ia,ib) - eepser(ia,ib))*ww * awt
       reskap(keq,ires) =
     & reskap(keq,ires) + (ekape(ia,ib) - ekaper(ia,ib))*ww * awt
       resfor(keq,ires) =
     & resfor(keq,ires) + (efore(ia,ib) - eforer(ia,ib))*ww * awt
       resmom(keq,ires) =
     & resmom(keq,ires) + (emome(ia,ib) - emomer(ia,ib))*ww * awt

      do jj = 1, neln
        rest_t(keq,ires,1,jj) = 
     &  rest_t(keq,ires,1,jj) + ete_tj(ia,ib,1,1,jj)*ww * awt
        rest_t(keq,ires,2,jj) = 
     &  rest_t(keq,ires,2,jj) + ete_tj(ia,ib,2,2,jj)*ww * awt
        rest_t(keq,ires,3,jj) = 
     &  rest_t(keq,ires,3,jj) + ete_tj(ia,ib,1,2,jj)*ww * awt
     &                        + ete_tj(ia,ib,2,1,jj)*ww * awt
      enddo

c----------------------------------------------------------------------
c---- V1 residuals
      keq = 1
      ia = 1
       resfn(keq,ires) =
     & resfn(keq,ires) + (efn(ia) - efnr(ia))*ww * awt
       resgn(keq,ires) =
     & resgn(keq,ires) + (egn(ia) - egnr(ia))*ww * awt
      do jj = 1, neln
        resv_v(keq,ires,1,jj) = 
     &  resv_v(keq,ires,1,jj) + ev_vj(ia,1,jj)*ww * awt
        resv_v(keq,ires,2,jj) = 
     &  resv_v(keq,ires,2,jj) + ev_vj(ia,2,jj)*ww * awt
      enddo

c----------------------------------------------------------------------
c---- V2 residuals
      keq = 2
      ia = 2
       resfn(keq,ires) =
     & resfn(keq,ires) + (efn(ia) - efnr(ia))*ww * awt
       resgn(keq,ires) =
     & resgn(keq,ires) + (egn(ia) - egnr(ia))*ww * awt
      do jj = 1, neln
        resv_v(keq,ires,1,jj) = 
     &  resv_v(keq,ires,1,jj) + ev_vj(ia,1,jj)*ww * awt
        resv_v(keq,ires,2,jj) = 
     &  resv_v(keq,ires,2,jj) + ev_vj(ia,2,jj)*ww * awt
      enddo

c----------------------------------------------------------------------
c---- integrate cell area
      ac(ires) = ac(ires)  + awt 
c
 100  continue ! next weighting function Wi
c--------------------------------------------------
c
 500  continue ! next Gauss point

      return
      end ! hsmrfm
