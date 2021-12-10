
       subroutine hsmeqn(neln, lrcurv, ldrill,
     &                   parg,
     &                   par1   , par2    , par3    , par4    , 
     &                   var1   , var2    , var3    , var4    , 
     &             res, res_var1, res_var2, res_var3, res_var4,
     &             ares, dres )
c---------------------------------------------------------------------------
c     Sets up area-integral parts of shell equation residuals and Jacobians
c     for a four-node shell element.  Edge integral terms are not included.
c
c     All vectors are assumed to be in common xyz cartesian axes
c
c Inputs:
c    neln      number of element nodes (3 = tri, 4 = quad)
c    parg(.)   global parameters
c    par#(.)   node parameters
c    var#(.)   node primary variables
c
c   Note: If neln=3, then par4(.) and var4(.) will not affect results
c
c Outputs:
c   res( ..)   residuals        int int [ R ] Wi dA
c      ( 1.)   x-force equilibrium
c      ( 2.)   y-force equilibrium
c      ( 3.)   z-force equilibrium
c      ( 4.)   x-moment equilibrium
c      ( 5.)   y-moment equilibrium
c      ( 6.)   z-moment equilibrium
c
c      res_var#(..)     d(res)/d(var#)
c
c      ares  area integral of weight function,  int int [ 1 ] Wi dA
c---------------------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c---- gauss-integration locations and weights
      include 'gausst.inc'
      include 'gaussq.inc'

c----------------------------------------------
c---- passed arguments
      logical lrcurv, ldrill
      real parg(lgtot)
      real par1(lvtot),
     &     par2(lvtot),
     &     par3(lvtot),
     &     par4(lvtot)

      real var1(ivtot),
     &     var2(ivtot),
     &     var3(ivtot),
     &     var4(ivtot)

      real res(irtot,4),
     &     res_var1(irtot,ivtot,4),
     &     res_var2(irtot,ivtot,4),
     &     res_var3(irtot,ivtot,4),
     &     res_var4(irtot,ivtot,4)

      real ares(4), dres(4)

c----------------------------------------------
c---- local work arrays
      parameter (ngdim = ngaussq)

      real gee(3), vel(3), rot(3), vac(3), rac(3)

      real pars(lvtot,4),
     &     vars(ivtot,4)

      real aaj(2,2,2,2,4), bbj(2,2,2,2,4),  ddj(2,2,2,2,4)
      real ssj(2,2,4)

      real r0j(3,4), rj(3,4)
      real d0j(3,4), dj(3,4) 
      real p0j(4), pj(4) 

c---- element cartesian basis work arrays
      real d1r0(3), d2r0(3)

      real crs(3)
      real c1(3)
      real c2(3)
      real c3(3)

      real cai(2,2)

c- - - - - - - - - - - - - - - - - - - - - - - - - - 
c---- Gauss point quantities for calling HSMGEO
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

c- - - - - - - - - - - - - - - - - - - - - - - - - - 
c---- quantities for evaluation of residual integrands
      real arxr_rj(3,4)

      real wxr(3), wxr_rj(3,3,4), wxr_dj(3,3,4), wxr_pj(3,4)

      real eps(2,2), eps_rj(2,2,3,4), eps_dj(2,2,3,4), eps_pj(2,2,4)
      real kap(2,2), kap_rj(2,2,3,4), kap_dj(2,2,3,4), kap_pj(2,2,4)

      real mu
      real qxyz(3)
      real txyz(3)

      real acc(3), acc_rj(3,3,4), acc_dj(3,3,4), acc_pj(3,4)
      real qld(3), qld_rj(3,3,4), qld_dj(3,3,4), qld_pj(3,4)
      real tld(3)
      
      real aac(2,2,2,2), bbc(2,2,2,2),  ddc(2,2,2,2)
      real ssc(2,2)

      real aa(2,2,2,2), bb(2,2,2,2),  dd(2,2,2,2)
      real ss(2,2)

      real fab(2,2), fab_rj(2,2,3,4), fab_dj(2,2,3,4), fab_pj(2,2,4)
      real mab(2,2), mab_rj(2,2,3,4), mab_dj(2,2,3,4), mab_pj(2,2,4)

      real fan(2), fan_rj(2,3,4), fan_dj(2,3,4), fan_pj(2,4)

      real gan(2), gan_rj(2,3,4), gan_dj(2,3,4), gan_pj(2,4)
      real gam(3), gam_rj(3,3,4), gam_dj(3,3,4), gam_pj(3,4)

      real nn(4)
      real d1nn(4), d2nn(4)


      real n(3), n_rj(3,3,4), n_dj(3,3,4), n_pj(3,4)

      real dxa1(3), dxa1_rj(3,3,4), dxa1_dj(3,3,4), dxa1_pj(3,4)
      real dxa2(3), dxa2_rj(3,3,4), dxa2_dj(3,3,4), dxa2_pj(3,4)

      real rxf, rxf_rj(3,4), rxf_dj(3,4), rxf_pj(4)

      real ajac0
      real ajac, ajac_rj(3,4), ajac_dj(3,4), ajac_pj(4)
      real awt, awt_rj(3,4), awt_dj(3,4), awt_pj(4)

      real g0w(3)
      real gw(3), gw_rj(3,3,4), gw_dj(3,3,4)

      real resi(irtot),
     &     resi_rj(irtot,3,4),
     &     resi_dj(irtot,3,4),
     &     resi_pj(irtot,4)

c      real
c     & lam(3), 
c     & d1lam(3),
c     & d2lam(3),
c     & lsq, 
c     & d1hll, 
c     & d2hll, 
c     & w, 
c     & d1w, 
c     & d2w, 
c     & v(3), 
c     & d1v(3),
c     & d2v(3),
c     & vxd(3), 
c     & d1vxd(3), 
c     & d2vxd(3)  
c
c      common / rlam / 
c     & lam, 
c     & d1lam,
c     & d2lam,
c     & lsq, 
c     & d1hll, 
c     & d2hll, 
c     & w, 
c     & d1w, 
c     & d2w, 
c     & v, 
c     & d1v,
c     & d2v,
c     & vxd, 
c     & d1vxd, 
c     & d2vxd  

c----------------------------------------------

      common / rcom /
     &  rg(3,ngaussq),
     &  gang(2,ngaussq),
     &  fang(2,ngaussq),
     &  gncg(3,ngaussq),
     &  fncg(3,ngaussq)


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

c---- unpack global parameters
      do k = 1, 3
        gee(k) = parg(lggeex-1+k)
        vel(k) = parg(lgvelx-1+k)
        rot(k) = parg(lgrotx-1+k)
        vac(k) = parg(lgvacx-1+k)
        rac(k) = parg(lgracx-1+k)
      enddo

c---- put passed node quantities into more convenient node-indexed arrays
      do lv = 1, lvtot
        pars(lv,1) = par1(lv)
        pars(lv,2) = par2(lv)
        pars(lv,3) = par3(lv)
        pars(lv,4) = par4(lv)
      enddo

      do iv = 1, ivtot
        vars(iv,1) = var1(iv)
        vars(iv,2) = var2(iv)
        vars(iv,3) = var3(iv)
        vars(iv,4) = var4(iv)
      enddo

c---- further unpack vectors
      do j = 1, 4
        do k = 1, 3
          r0j(k,j) = pars(lvr0x-1+k,j)
          d0j(k,j) = pars(lvn0x-1+k,j)
          rj(k,j) = vars(ivrx-1+k,j)
          dj(k,j) = vars(ivdx-1+k,j)
        enddo
        p0j(j) = 0.
        pj(j) = vars(ivps,j)
      enddo

c---- element cartesian basis vectors at element centroid
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

c---- element cartesian basis vectors
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

c---- element cartesian basis vectors
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

c---- put node stiffness and compliance matrices into element cartesian basis
      call hsmabd(neln,
     &            par1, par2, par3, par4, 
     &            aaj, bbj, ddj, ssj )

c-------------------------------------------------------------------------
c---- see element geometry quantities at Gauss point
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
        do keq = 1, irtot
          res(keq,ires) = 0.
          do iv = 1, ivtot
            res_var1(keq,iv,ires) = 0.
            res_var2(keq,iv,ires) = 0.
            res_var3(keq,iv,ires) = 0.
            res_var4(keq,iv,ires) = 0.
          enddo
        enddo
        ares(ires) = 0.
        dres(ires) = 0.
      enddo

c-------------------------------------------------------------------
c---- sum over gauss points
      do 100 ig = 1, ngauss

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

c      xc = 0.25*(rj(1,1)+rj(1,2)+rj(1,3)+rj(1,4))
c      yc = 0.25*(rj(2,1)+rj(2,2)+rj(2,3)+rj(2,4))
c      if(xc .gt. 0.49 .and. xc .lt. 0.51 .and.
c     &   yc .gt. 0.49 .and. yc .lt. 0.51       ) then
c         write(*,*) 'writing fort.90'
c         lu = 90
c      do jt = -10, 10
c      do it = -10, 10
c        ett = float(jt)/10.0
c        xit = float(it)/10.0
c      call hsmgeo(neln, xit, ett,
c     &            rj(1,1), rj(1,2), rj(1,3), rj(1,4),
c     &            dj(1,1), dj(1,2), dj(1,3), dj(1,4),
c     &             d,  d_rj,     d_dj,
c     &            a1, a1_rj,
c     &            a2, a2_rj,
c     &            a1i, a1i_rj,
c     &            a2i, a2i_rj,
c     &            g  ,   g_rj,
c     &            h  ,   h_rj,   h_dj,
c     &            l  ,   l_rj,   l_dj,
c     &            ginv, ginv_rj,
c     &            gdet, gdet_rj,
c     &            rxr ,  rxr_rj )
c      write(lu,33) 
c     &  xit,ett, d(1), d(2), lam(1), lam(2)
c 33   format(1x,20g12.4)
c      enddo
c      write(lu,*)
c      enddo
c      stop



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

c---- stress and stress-moment resultant tensor components, f^ab, m^ab
      do ib = 1, 2
        do ia = 1, 2
          fab(ia,ib) = aa(ia,ib,1,1)*eps(1,1)
     &               + aa(ia,ib,1,2)*eps(1,2)
     &               + aa(ia,ib,2,1)*eps(2,1)
     &               + aa(ia,ib,2,2)*eps(2,2)
     &               + bb(ia,ib,1,1)*kap(1,1)
     &               + bb(ia,ib,1,2)*kap(1,2)
     &               + bb(ia,ib,2,1)*kap(2,1)
     &               + bb(ia,ib,2,2)*kap(2,2)
          mab(ia,ib) = bb(ia,ib,1,1)*eps(1,1)
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
     &                 aa(ia,ib,1,1)*eps_pj(1,1,jj)
     &               + aa(ia,ib,1,2)*eps_pj(1,2,jj)
     &               + aa(ia,ib,2,1)*eps_pj(2,1,jj)
     &               + aa(ia,ib,2,2)*eps_pj(2,2,jj)
     &               + bb(ia,ib,1,1)*kap_pj(1,1,jj)
     &               + bb(ia,ib,1,2)*kap_pj(1,2,jj)
     &               + bb(ia,ib,2,1)*kap_pj(2,1,jj)
     &               + bb(ia,ib,2,2)*kap_pj(2,2,jj)
            mab_pj(ia,ib,jj) = 
     &                 bb(ia,ib,1,1)*eps_pj(1,1,jj)
     &               + bb(ia,ib,1,2)*eps_pj(1,2,jj)
     &               + bb(ia,ib,2,1)*eps_pj(2,1,jj)
     &               + bb(ia,ib,2,2)*eps_pj(2,2,jj)
     &               + dd(ia,ib,1,1)*kap_pj(1,1,jj)
     &               + dd(ia,ib,1,2)*kap_pj(1,2,jj)
     &               + dd(ia,ib,2,1)*kap_pj(2,1,jj)
     &               + dd(ia,ib,2,2)*kap_pj(2,2,jj)
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

      
c        xc = 0.5*(r0j(1,2) + r0j(1,1))
c        if(abs(xc-0.51) .lt. 0.001) then
c         write(*,*) 'ax g1 S11 fn1', a01(1), gan(1), ss(1,1), fan(1)
c        endif


c---- shear tilt vector in cartesian basis
      do k = 1, 3
        gam(k) = a1i(k,ig)*gan(1) + a2i(k,ig)*gan(2)
        do jj = 1, neln
          do kk = 1, 3
            gam_rj(k,kk,jj) = a1i_rj(k,kk,jj,ig)*gan(1)
     &                      + a2i_rj(k,kk,jj,ig)*gan(2)
     &                      + a1i(k,ig)*gan_rj(1,kk,jj)
     &                      + a2i(k,ig)*gan_rj(2,kk,jj)

            gam_dj(k,kk,jj) = a1i_dj(k,kk,jj,ig)*gan(1)
     &                      + a2i_dj(k,kk,jj,ig)*gan(2)
     &                      + a1i(k,ig)*gan_dj(1,kk,jj)
     &                      + a2i(k,ig)*gan_dj(2,kk,jj)
          enddo
          gam_pj(k,jj) = a1i_pj(k,jj,ig)*gan(1)
     &                 + a2i_pj(k,jj,ig)*gan(2)
     &                 + a1i(k,ig)*gan_pj(1,jj)
     &                 + a2i(k,ig)*gan_pj(2,jj)
        enddo
      enddo

c---- normal vector (director with shear tilt vector removed)
      do k = 1, 3
        n(k) = d(k,ig) - gam(k)
        do jj = 1, neln
          do kk = 1, 3
            n_rj(k,kk,jj) = d_rj(k,kk,jj,ig) - gam_rj(k,kk,jj)
            n_dj(k,kk,jj) = d_dj(k,kk,jj,ig) - gam_dj(k,kk,jj)
          enddo
          n_pj(k,jj) = d_pj(k,jj,ig) - gam_pj(k,jj)
        enddo
      enddo

c---- d x a1 ,  d x a2  vectors
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        dxa1(k) = d(ic,ig)*a1(jc,ig)
     &          - d(jc,ig)*a1(ic,ig)
        dxa2(k) = d(ic,ig)*a2(jc,ig)
     &          - d(jc,ig)*a2(ic,ig)
        do jj = 1, neln
          do kk = 1, 3
            dxa1_rj(k,kk,jj) = 
     &        d_rj(ic,kk,jj,ig)*a1(jc,ig) + d(ic,ig)*a1_rj(jc,kk,jj,ig)
     &      - d_rj(jc,kk,jj,ig)*a1(ic,ig) - d(jc,ig)*a1_rj(ic,kk,jj,ig)
            dxa2_rj(k,kk,jj) = 
     &        d_rj(ic,kk,jj,ig)*a2(jc,ig) + d(ic,ig)*a2_rj(jc,kk,jj,ig)
     &      - d_rj(jc,kk,jj,ig)*a2(ic,ig) - d(jc,ig)*a2_rj(ic,kk,jj,ig)

            dxa1_dj(k,kk,jj) = 
     &        d_dj(ic,kk,jj,ig)*a1(jc,ig) + d(ic,ig)*a1_dj(jc,kk,jj,ig)
     &      - d_dj(jc,kk,jj,ig)*a1(ic,ig) - d(jc,ig)*a1_dj(ic,kk,jj,ig)
            dxa2_dj(k,kk,jj) = 
     &        d_dj(ic,kk,jj,ig)*a2(jc,ig) + d(ic,ig)*a2_dj(jc,kk,jj,ig)
     &      - d_dj(jc,kk,jj,ig)*a2(ic,ig) - d(jc,ig)*a2_dj(ic,kk,jj,ig)
          enddo

          dxa1_pj(k,jj) = 
     &        d_pj(ic,jj,ig)*a1(jc,ig) + d(ic,ig)*a1_pj(jc,jj,ig)
     &      - d_pj(jc,jj,ig)*a1(ic,ig) - d(jc,ig)*a1_pj(ic,jj,ig)
          dxa2_pj(k,jj) = 
     &        d_pj(ic,jj,ig)*a2(jc,ig) + d(ic,ig)*a2_pj(jc,jj,ig)
     &      - d_pj(jc,jj,ig)*a2(ic,ig) - d(jc,ig)*a2_pj(ic,jj,ig)
        enddo
      enddo

c---- mass/area
      mu = pars(lvmu,1)*nn(1)
     &   + pars(lvmu,2)*nn(2)
     &   + pars(lvmu,3)*nn(3)
     &   + pars(lvmu,4)*nn(4)

c---- normal force/area
      qn = pars(lvqn,1)*nn(1)
     &   + pars(lvqn,2)*nn(2)
     &   + pars(lvqn,3)*nn(3)
     &   + pars(lvqn,4)*nn(4)

c---- xyz force/area
      do k = 1, 3
        lv = lvqx - 1 + k
        qxyz(k) = pars(lv,1)*nn(1)
     &          + pars(lv,2)*nn(2)
     &          + pars(lv,3)*nn(3)
     &          + pars(lv,4)*nn(4)
      enddo

c---- xyz moment/area
      do k = 1, 3
        lv = lvtx - 1 + k
        txyz(k) = pars(lv,1)*nn(1)
     &          + pars(lv,2)*nn(2)
     &          + pars(lv,3)*nn(3)
     &          + pars(lv,4)*nn(4)
      enddo

c---- total load/area vectors
      do k = 1, 3
        qld(k) = qxyz(k) + qn*n(k)
        tld(k) = txyz(k)
        do jj = 1, neln
          do kk = 1, 3
            qld_rj(k,kk,jj) = qn*n_rj(k,kk,jj)
            qld_dj(k,kk,jj) = qn*n_dj(k,kk,jj)
          enddo
          qld_pj(k,jj) = qn*n_pj(k,jj)
        enddo
      enddo

c      write(*,*) 
c      write(*,*) 'ig', ig
c      write(*,*) 'a1 ', a01
c      write(*,*) 'a1i', a01i
c      write(*,*) 'A',
c     &    aa(1,1,1,1),
c     &    aa(2,2,2,2),
c     &    aa(1,1,2,2),
c     &    aa(1,1,1,2),
c     &    aa(2,2,1,2),
c     &    aa(1,2,1,2)
c      write(*,*) 'eps11 22 12', eps(1,1), eps(2,2), eps(1,2)
c      write(*,*) 'fab11 22 12', fab(1,1), fab(2,2), fab(1,2)
cc     if(ig.eq.4) pause

c      do k = 1, 3
c        rg(k,ig) = r(k)
c        gncg(k,ig) = gam(k)
c        fncg(k,ig) = fan(1)*a1(k) + fan(2)*a2(k)
c      enddo
c      do ia = 1, 2
c        gang(ia,ig) = gan(ia)
c        fang(ia,ig) = fan(ia)
c      enddo

c---- W x r
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        wxr(k) = rot(ic)*r(jc,ig) - rot(jc)*r(ic,ig)
        do jj = 1, neln
          do kk = 1, 3
            wxr_rj(k,kk,jj) = rot(ic)*r_rj(jc,kk,jj,ig)
     &                      - rot(jc)*r_rj(ic,kk,jj,ig)
            wxr_dj(k,kk,jj) = rot(ic)*r_dj(jc,kk,jj,ig)
     &                      - rot(jc)*r_dj(ic,kk,jj,ig)
          enddo
          wxr_pj(k,jj) = rot(ic)*r_pj(jc,jj,ig)
     &                 - rot(jc)*r_pj(ic,jj,ig)
        enddo
      enddo

c---- absolute acceleration
c-    a  =  Udot  +  W x U  +  Wdot x r  +  W x (W x r)
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        acc(k) = vac(k)
     &         + rot(ic)*vel(jc)  - rot(jc)*vel(ic)
     &         + rac(ic)*r(jc,ig) - rac(jc)*r(ic,ig)
     &         + rot(ic)*wxr(jc)  - rot(jc)*wxr(ic)
        do jj = 1, neln
          do kk = 1, 3
            acc_rj(k,kk,jj) = 
     &           rac(ic)*r_rj(jc,kk,jj,ig) - rac(jc)*r_rj(ic,kk,jj,ig)
     &         + rot(ic)*wxr_rj(jc,kk,jj)  - rot(jc)*wxr_rj(ic,kk,jj)
            acc_dj(k,kk,jj) = 
     &           rac(ic)*r_dj(jc,kk,jj,ig) - rac(jc)*r_dj(ic,kk,jj,ig)
     &         + rot(ic)*wxr_dj(jc,kk,jj)  - rot(jc)*wxr_dj(ic,kk,jj)
          enddo
          acc_pj(k,jj) = 
     &           rac(ic)*r_pj(jc,jj,ig) - rac(jc)*r_pj(ic,jj,ig)
     &         + rot(ic)*wxr_pj(jc,jj)  - rot(jc)*wxr_pj(ic,jj)
        enddo
      enddo

c---- coordinate Jacobian J
c      ajac0 = sqrt(abs(g0det(ig)))
c      ajac = sqrt(abs(gdet(ig)))
c      ajac_gdet = 0.5/ajac
c      do jj = 1, neln
c        do kk = 1, 3
c          ajac_rj(kk,jj) = ajac_gdet*gdet_rj(kk,jj,ig)
c          ajac_dj(kk,jj) = ajac_gdet*gdet_dj(kk,jj,ig)
c        enddo
c        ajac_pj(jj) = ajac_gdet*gdet_pj(jj,ig)
c      enddo

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
      do jj = 1, 4
        do kk = 1, 3
          awt_rj(kk,jj) = 0.
          awt_dj(kk,jj) = 0.
        enddo
        awt_pj(jj) = 0.
      enddo

c--------------------------------------------------
c---- go over the four residual weight functions of this element
      do 10 ires = 1, neln

c---- clear working arrays for this weight function
      do keq = 1, irtot
        resi(keq) = 0.
        do jj = 1, neln
          do kk = 1, 3
            resi_rj(keq,kk,jj) = 0.
            resi_dj(keq,kk,jj) = 0.
          enddo
          resi_pj(keq,jj) = 0.
        enddo
      enddo

c---- weight function  Wi = Ni  and  d(Wi) = d(Ni)  at this gauss point
      ww = nn(ires)
      d1ww = d1nn(ires)
      d2ww = d2nn(ires)

c---- grad0(W) and grad(W)
      do k = 1, 3
        g0w(k) = a01i(k,ig)*d1ww
     &         + a02i(k,ig)*d2ww
        gw(k) = a1i(k,ig)*d1ww
     &        + a2i(k,ig)*d2ww
        do jj = 1, neln
          do kk = 1, 3
            gw_rj(k,kk,jj) = a1i_rj(k,kk,jj,ig)*d1ww
     &                     + a2i_rj(k,kk,jj,ig)*d2ww
            gw_dj(k,kk,jj) = a1i_dj(k,kk,jj,ig)*d1ww
     &                     + a2i_dj(k,kk,jj,ig)*d2ww
          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c---- force equilibrium
      do k = 1, 3
        keq = k
        resi(keq) = 
     &    - (fab(1,1)*d1ww + fab(1,2)*d2ww)*a1(k,ig)
     &    - (fab(2,1)*d1ww + fab(2,2)*d2ww)*a2(k,ig)
     &    - (fan(1)  *d1ww + fan(2)  *d2ww)* n(k)
     &    +  qld(k)*ww
     &    +  mu*(gee(k)-acc(k))*ww
        do jj = 1, neln
          do kk = 1, 3
            resi_rj(keq,kk,jj) = 
     &       - ( fab_rj(1,1,kk,jj)*d1ww
     &         + fab_rj(1,2,kk,jj)*d2ww )*a1(k,ig)
     &       - ( fab_rj(2,1,kk,jj)*d1ww
     &         + fab_rj(2,2,kk,jj)*d2ww )*a2(k,ig)
     &       - ( fan_rj(1  ,kk,jj)*d1ww
     &         + fan_rj(2  ,kk,jj)*d2ww )* n(k)
     &       - ( fab(1,1)*d1ww
     &         + fab(1,2)*d2ww )*a1_rj(k,kk,jj,ig)
     &       - ( fab(2,1)*d1ww
     &         + fab(2,2)*d2ww )*a2_rj(k,kk,jj,ig)
     &       - ( fan(1)  *d1ww
     &         + fan(2)  *d2ww )* n_rj(k,kk,jj)
     &       +  qld_rj(k,kk,jj)*ww
     &       +  mu*(-acc_rj(k,kk,jj))*ww

            resi_dj(keq,kk,jj) = 
     &       - ( fab_dj(1,1,kk,jj)*d1ww
     &         + fab_dj(1,2,kk,jj)*d2ww )*a1(k,ig)
     &       - ( fab_dj(2,1,kk,jj)*d1ww
     &         + fab_dj(2,2,kk,jj)*d2ww )*a2(k,ig)
     &       - ( fan_dj(1  ,kk,jj)*d1ww
     &         + fan_dj(2  ,kk,jj)*d2ww )* n(k)
     &       - ( fab(1,1)*d1ww
     &         + fab(1,2)*d2ww )*a1_dj(k,kk,jj,ig)
     &       - ( fab(2,1)*d1ww
     &         + fab(2,2)*d2ww )*a2_dj(k,kk,jj,ig)
     &       - ( fan(1)  *d1ww
     &         + fan(2)  *d2ww )* n_dj(k,kk,jj)
     &       +  qld_dj(k,kk,jj)*ww
     &       +  mu*(-acc_dj(k,kk,jj))*ww
          enddo
          resi_pj(keq,jj) = 
     &       - ( fab_pj(1,1,jj)*d1ww
     &         + fab_pj(1,2,jj)*d2ww )*a1(k,ig)
     &       - ( fab_pj(2,1,jj)*d1ww
     &         + fab_pj(2,2,jj)*d2ww )*a2(k,ig)
     &       - ( fan_pj(1  ,jj)*d1ww
     &         + fan_pj(2  ,jj)*d2ww )* n(k)
     &       - ( fab(1,1)*d1ww
     &         + fab(1,2)*d2ww )*a1_pj(k,jj,ig)
     &       - ( fab(2,1)*d1ww
     &         + fab(2,2)*d2ww )*a2_pj(k,jj,ig)
     &       - ( fan(1)  *d1ww
     &         + fan(2)  *d2ww )* n_pj(k,jj)
     &       +  qld_pj(k,jj)*ww
     &       +  mu*(-acc_pj(k,jj))*ww
        enddo
      enddo

c----------------------------------------------------------------------
c---- moment equilibrium
c      ww = 0.
      do k = 1, 3
        keq = 3 + k

c------ set all r x [ ] terms together, re-using force residual resi(1:3)
        ic = icrs(k)
        jc = jcrs(k)
        rxf = (r(ic,ig)-rj(ic,ires))*resi(jc)
     &      - (r(jc,ig)-rj(jc,ires))*resi(ic)
        do jj = 1, neln
          do kk = 1, 3
            rxf_rj(kk,jj) = 
     &         r_rj(ic,kk,jj,ig)*resi(jc)
     &      -  r_rj(jc,kk,jj,ig)*resi(ic)
     &      + (r(ic,ig)-rj(ic,ires))*resi_rj(jc,kk,jj)
     &      - (r(jc,ig)-rj(jc,ires))*resi_rj(ic,kk,jj)
            rxf_dj(kk,jj) = 
     &         r_dj(ic,kk,jj,ig)*resi(jc)
     &      -  r_dj(jc,kk,jj,ig)*resi(ic)
     &      + (r(ic,ig)-rj(ic,ires))*resi_dj(jc,kk,jj)
     &      - (r(jc,ig)-rj(jc,ires))*resi_dj(ic,kk,jj)
          enddo
          rxf_pj(jj) = 
     &         r_pj(ic,jj,ig)*resi(jc)
     &      -  r_pj(jc,jj,ig)*resi(ic)
     &      + (r(ic,ig)-rj(ic,ires))*resi_pj(jc,jj)
     &      - (r(jc,ig)-rj(jc,ires))*resi_pj(ic,jj)
        enddo
        rxf_rj(ic,ires) = rxf_rj(ic,ires) - resi(jc)
        rxf_rj(jc,ires) = rxf_rj(jc,ires) + resi(ic)

c------ overall moment residual
        resi(keq) = -(mab(1,1)*d1ww + mab(1,2)*d2ww)*dxa1(k)
     &              -(mab(2,1)*d1ww + mab(2,2)*d2ww)*dxa2(k)
     &            +  tld(k)*ww
     &            +  rxf
        r_dxa1    = -(mab(1,1)*d1ww + mab(1,2)*d2ww)
        r_dxa2    = -(mab(2,1)*d1ww + mab(2,2)*d2ww)
        do jj = 1, neln
          do kk = 1, 3
            resi_rj(keq,kk,jj) = 
     &      - ( mab_rj(1,1,kk,jj)*d1ww
     &        + mab_rj(1,2,kk,jj)*d2ww )*dxa1(k)
     &      - ( mab_rj(2,1,kk,jj)*d1ww
     &        + mab_rj(2,2,kk,jj)*d2ww )*dxa2(k)
     &      + r_dxa1*dxa1_rj(k,kk,jj)
     &      + r_dxa2*dxa2_rj(k,kk,jj)
     &      + rxf_rj(kk,jj)

            resi_dj(keq,kk,jj) = 
     &      - ( mab_dj(1,1,kk,jj)*d1ww
     &        + mab_dj(1,2,kk,jj)*d2ww )*dxa1(k)
     &      - ( mab_dj(2,1,kk,jj)*d1ww
     &        + mab_dj(2,2,kk,jj)*d2ww )*dxa2(k)
     &      + r_dxa1*dxa1_dj(k,kk,jj)
     &      + r_dxa2*dxa2_dj(k,kk,jj)
     &      + rxf_dj(kk,jj)
          enddo
          resi_pj(keq,jj) = 
     &      - ( mab_pj(1,1,jj)*d1ww
     &        + mab_pj(1,2,jj)*d2ww )*dxa1(k)
     &      - ( mab_pj(2,1,jj)*d1ww
     &        + mab_pj(2,2,jj)*d2ww )*dxa2(k)
     &      + r_dxa1*dxa1_pj(k,jj)
     &      + r_dxa2*dxa2_pj(k,jj)
     &      + rxf_pj(jj)
        enddo

      enddo

c---- keq slot used for checking linearization of any quantiity via testeqn
c      keq = 6
c      k = 1
c      j = 1
c      resi(keq) = a1i(k)
ccc     &    - (fab(1,1)*d1ww + fab(1,2)*d2ww)*a1(k)
ccc     &    - (fab(2,1)*d1ww + fab(2,2)*d2ww)*a2(k)
ccc     &    - (fan(1)  *d1ww + fan(2)  *d2ww)* n(k)
ccc     &    +  qld(k)*ww
ccc     &    +  mu*(gee(k)-acc(k))*ww
c      do jj = 1, neln
c        do kk = 1, 3
c          resi_rj(keq,kk,jj) = a1i_rj(k,kk,jj)
c          resi_dj(keq,kk,jj) = a1i_dj(k,kk,jj)
ccc          resi_rj(keq,kk,jj) = fab_rj(1,1,kk,jj)
ccc          resi_dj(keq,kk,jj) = fab_dj(1,1,kk,jj)
c        enddo
c        resi_pj(keq,jj) = a1i_pj(k,jj)
cccc       resi_pj(keq,jj) = fab_pj(1,1,jj)
c      enddo


c----------------------------------------------------------------------
c---- accumulate area-integrated residuals and Jacobians
      do keq = 1, irtot
        res(keq,ires) = 
     &  res(keq,ires)  +  resi(keq) * awt
        do kk = 1, 3
          ivr = ivrx-1 + kk
          ivd = ivdx-1 + kk

          res_var1(keq,ivr,ires) = 
     &    res_var1(keq,ivr,ires) + resi_rj(keq,kk,1) * awt
     &                           + resi(keq) * awt_rj(kk,1)
          res_var2(keq,ivr,ires) = 
     &    res_var2(keq,ivr,ires) + resi_rj(keq,kk,2) * awt
     &                           + resi(keq) * awt_rj(kk,2)
          res_var3(keq,ivr,ires) = 
     &    res_var3(keq,ivr,ires) + resi_rj(keq,kk,3) * awt
     &                           + resi(keq) * awt_rj(kk,3)
          res_var4(keq,ivr,ires) = 
     &    res_var4(keq,ivr,ires) + resi_rj(keq,kk,4) * awt
     &                           + resi(keq) * awt_rj(kk,4)

          res_var1(keq,ivd,ires) = 
     &    res_var1(keq,ivd,ires) + resi_dj(keq,kk,1) * awt
     &                           + resi(keq) * awt_dj(kk,1)
          res_var2(keq,ivd,ires) = 
     &    res_var2(keq,ivd,ires) + resi_dj(keq,kk,2) * awt
     &                           + resi(keq) * awt_dj(kk,2)
          res_var3(keq,ivd,ires) = 
     &    res_var3(keq,ivd,ires) + resi_dj(keq,kk,3) * awt
     &                           + resi(keq) * awt_dj(kk,3)
          res_var4(keq,ivd,ires) = 
     &    res_var4(keq,ivd,ires) + resi_dj(keq,kk,4) * awt
     &                           + resi(keq) * awt_dj(kk,4)
        enddo

        res_var1(keq,ivps,ires) = 
     &  res_var1(keq,ivps,ires) + resi_pj(keq,1) * awt
     &                          + resi(keq) * awt_pj(1)
        res_var2(keq,ivps,ires) = 
     &  res_var2(keq,ivps,ires) + resi_pj(keq,2) * awt
     &                          + resi(keq) * awt_pj(2)
        res_var3(keq,ivps,ires) = 
     &  res_var3(keq,ivps,ires) + resi_pj(keq,3) * awt
     &                          + resi(keq) * awt_pj(3)
        res_var4(keq,ivps,ires) = 
     &  res_var4(keq,ivps,ires) + resi_pj(keq,4) * awt
     &                          + resi(keq) * awt_pj(4)
      enddo ! next equation keq

c---- integrate element area
      ares(ires) = ares(ires) + ww * awt

  10  continue ! next weighting function Wi

 100  continue ! next Gauss point

      return
      end ! hsmeqn
