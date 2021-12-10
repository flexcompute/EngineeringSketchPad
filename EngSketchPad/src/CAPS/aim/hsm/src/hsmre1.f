
       subroutine hsmre1(neln, lrcurv, ldrill,
     &                   par1    , par2    , par3    , par4    , 
     &                   var1    , var2    , var3    , var4    ,
     &                   dep1    , dep2    , dep3    , dep4    ,
     &                   e1j1    , e1j2    , e1j3    , e1j4    ,
     &              rese, rese_e1j,
     &              ac )
c------------------------------------------------------------
c     Sets up residuals and Jacobians for e1 and en vector
c     definitions for a four-node shell element.
c
c Inputs:
c   par#(.)   node parameters
c   var#(.)   node primary variables
c   dep#(.)   node dependent (secondary) variables
c   e1j#(.)    node tangent unit vectors e1
c
c Outputs:
c   rese(.)   e residuals,                        int int [ Re ] Wi/J^2 dA
c   rese_e1j#(..)  d(rese)/d(e1j#)
c   ac(.)    area integrals of weight functions,  int int [ 1 ] Wi/J^2 dA
c------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

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

      real e1j1(3),
     &     e1j2(3),
     &     e1j3(3),
     &     e1j4(3)

      real rese(4), rese_e1j(4,3)

      real ac(4)

c----------------------------------------------
c---- local arrays
      real par(lvtot,4),
     &     var(ivtot,4),
     &     dep(jvtot,4)

      real e1j(3,4),
     &     enj(3,4)

      real e1ai(2), e2ai(2)
      real e1a(2), e2a(2)

      real d1r0(3),
     &     d2r0(3)

      real d1r(3),
     &     d2r(3)

      real rj(3,4), r0j(3,4)
      real dj(3,4)

      real e01j(3,4), e02j(3,4), e0nj(3,4)

      real r(3)
      real d(3)

      real a01(3), a1(3), a1_rj(3,3,4)
      real a02(3), a2(3), a2_rj(3,3,4)

      real a01i(3), a1i(3), a1i_rj(3,3,4)
      real a02i(3), a2i(3), a2i_rj(3,3,4)

      real g0(2,2)
      real g0inv(2,2)

      real g(2,2), g_rj(2,2,3,4)
      real gdet, gdet_rj(3,4)
      real ginv(2,2), ginv_rj(2,2,3,4)

      real eps(2,2), eps_rj(2,2,3,4)

      real nn(4)
      real d1nn(4), d2nn(4)


      real s0(3)
      real s(3)

      real sde, sde_ei(3)
      
      real sxe0(3)
      real sxe(3), sxe_ei(3,3)
      real sxen, sxen_ei(3)

      real th_ei(3)

      real rxr0(3)
      real arxr0

      real rxr(3), rxr_rj(3,3,4)
      real arxr, arxr_rj(3,4)

      real ajac0
      real ajac

      real resei, resei_ei(3)

c----------------------------------------------
c---- gauss-integration locations and weights
      include 'gausst.inc'
      include 'gaussq.inc'

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------

      if(neln .eq. 3) then
c----- use tri element quadrature
       ngauss = ngausst
      else
c----- use quad element quadrature
       ngauss = ngaussq
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

      do k = 1, 3
        e1j(k,1) = e1j1(k)
        e1j(k,2) = e1j2(k)
        e1j(k,3) = e1j3(k)
        e1j(k,4) = e1j4(k)
      enddo
c
c---- further unpack vectors
      do j = 1, 4
        r0j(1,j) = par(lvr0x,j)
        r0j(2,j) = par(lvr0y,j)
        r0j(3,j) = par(lvr0z,j)

        e01j(1,j) = par(lve01x,j)
        e01j(2,j) = par(lve01y,j)
        e01j(3,j) = par(lve01z,j)

        e02j(1,j) = par(lve02x,j)
        e02j(2,j) = par(lve02y,j)
        e02j(3,j) = par(lve02z,j)

        e0nj(1,j) = par(lvn0x,j)
        e0nj(2,j) = par(lvn0y,j)
        e0nj(3,j) = par(lvn0z,j)


        rj(1,j) = var(ivrx,j)
        rj(2,j) = var(ivry,j)
        rj(3,j) = var(ivrz,j)

        enj(1,j) = dep(jvnx,j)
        enj(2,j) = dep(jvny,j)
        enj(3,j) = dep(jvnz,j)
      enddo

c=======================================================================
c---- clear residuals and Jacobians for summing over gauss points
        do ires = 1, 4
          rese(ires) = 0.
          do kk = 1, 3
            rese_e1j(ires,kk) = 0.
          enddo
          ac(ires) = 0.
        enddo

c-------------------------------------------------------------------
c---- sum over gauss points
      do 500 igauss = 1, ngauss

c---- bilinear interpolation weights "N" and their derivatives at gauss point
      if(neln .eq. 3) then
       xsi = xsit(igauss)
       eta = etat(igauss)
       fwt = fwtt(igauss)

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
       xsi = xsiq(igauss)
       eta = etaq(igauss)
       fwt = fwtq(igauss)

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

c-------------------------------------------------------------------------
c---- r derivatives wrt element coordinates, dr/dxsi, dr/deta
c-    (covariant basis vectors)
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

      do k = 1, 3
        d1r(k) = rj(k,1)*d1nn(1)
     &         + rj(k,2)*d1nn(2)
     &         + rj(k,3)*d1nn(3)
     &         + rj(k,4)*d1nn(4)
        d2r(k) = rj(k,1)*d2nn(1)
     &         + rj(k,2)*d2nn(2)
     &         + rj(k,3)*d2nn(3)
     &         + rj(k,4)*d2nn(4)
      enddo

c---- element basis vectors
      do k = 1, 3
        a01(k) = d1r0(k)
        a02(k) = d2r0(k)

        a1(k) = d1r(k)
        a2(k) = d2r(k)
      enddo

c---- covariant metric tensors g_ij, and Jacobians wrt r
      g0(1,1) = a01(1)*a01(1)
     &        + a01(2)*a01(2)
     &        + a01(3)*a01(3)
      g0(2,2) = a02(1)*a02(1)
     &        + a02(2)*a02(2)
     &        + a02(3)*a02(3)
      g0(1,2) = a01(1)*a02(1)
     &        + a01(2)*a02(2)
     &        + a01(3)*a02(3)
      g0(2,1) = g0(1,2)
      g0det = g0(1,1)*g0(2,2) - g0(1,2)*g0(2,1)

      g(1,1) = a1(1)*a1(1)
     &       + a1(2)*a1(2)
     &       + a1(3)*a1(3)
      g(2,2) = a2(1)*a2(1)
     &       + a2(2)*a2(2)
     &       + a2(3)*a2(3)
      g(1,2) = a1(1)*a2(1)
     &       + a1(2)*a2(2)
     &       + a1(3)*a2(3)
      g(2,1) = g(1,2)

      gdet  =  g(1,1)*g(2,2) - g(1,2)*g(2,1)

c---- contravariant metric tensors  g^ij, and Jacobians wrt x,y
      g0inv(1,1) =  g0(2,2) / g0det
      g0inv(2,2) =  g0(1,1) / g0det
      g0inv(1,2) = -g0(1,2) / g0det
      g0inv(2,1) = -g0(2,1) / g0det
      
      ginv(1,1) =  g(2,2) / gdet
      ginv(2,2) =  g(1,1) / gdet
      ginv(1,2) = -g(1,2) / gdet
      ginv(2,1) = -g(2,1) / gdet

c---- contravariant basis vectors
      do k = 1, 3
        a01i(k) = g0inv(1,1)*a01(k) + g0inv(1,2)*a02(k)
        a02i(k) = g0inv(2,1)*a01(k) + g0inv(2,2)*a02(k)
      enddo

      do k = 1, 3
        a1i(k) = ginv(1,1)*a1(k) + ginv(1,2)*a2(k)
        a2i(k) = ginv(2,1)*a1(k) + ginv(2,2)*a2(k)
      enddo

c---- covariant Green strain tensor
      eps(1,1) = (g(1,1) - g0(1,1))*0.5
      eps(1,2) = (g(1,2) - g0(1,2))*0.5
      eps(2,1) = (g(2,1) - g0(2,1))*0.5
      eps(2,2) = (g(2,2) - g0(2,2))*0.5

c---- material-line angles
      do k = 1, 3
        s0(k) = d1r0(k)
        s(k) = d1r(k) - a1i(k)*eps(1,1)
     &                - a2i(k)*eps(1,2)
      enddo

c---- coordinate Jacobian J
      ajac0 = sqrt(abs(g0det))

      ajac = sqrt(abs(gdet))
      ajac_gdet = 0.5/ajac

c---- residual area-weight for this gauss point
cc    awt = fwt * ajac0
      awt = fwt / ajac0

c--------------------------------------------------
c---- go over the node residual weight functions
      do 100 ires = 1, neln

c---- clear working arrays for this weight function
      resei = 0.
      resei_ei(1) = 0.
      resei_ei(2) = 0.
      resei_ei(3) = 0.

c---- weight function  Wi = Ni at this gauss point
      ww = nn(ires)

      sde0 = s0(1)*e01j(1,ires)
     &     + s0(2)*e01j(2,ires)
     &     + s0(3)*e01j(3,ires)

      sde = s(1)*e1j(1,ires)
     &    + s(2)*e1j(2,ires)
     &    + s(3)*e1j(3,ires)
      do kk = 1, 3
        sde_ei(kk) = s(kk)
      enddo

      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        sxe0(k) = s0(ic)*e01j(jc,ires) - s0(jc)*e01j(ic,ires)

        sxe(k) = s(ic)*e1j(jc,ires) - s(jc)*e1j(ic,ires)
        sxe_ei(k,jc) =  s(ic)
        sxe_ei(k,ic) = -s(jc)
        sxe_ei(k,k)  = 0.
      enddo

      sxen0 = sxe0(1)*e0nj(1,ires)
     &      + sxe0(2)*e0nj(2,ires)
     &      + sxe0(3)*e0nj(3,ires)

      sxen = sxe(1)*enj(1,ires)
     &     + sxe(2)*enj(2,ires)
     &     + sxe(3)*enj(3,ires)
      do kk = 1, 3
        sxen_ei(kk) = sxe_ei(1,kk)*enj(1,ires)
     &              + sxe_ei(2,kk)*enj(2,ires)
     &              + sxe_ei(3,kk)*enj(3,ires)
      enddo

      th0 = atan2( sxen0 , sde0 )

      th = atanc( sxen , sde , th0 )
      th_sxen =   sde / (sxen**2 + sde**2)
      th_sde  = -sxen / (sxen**2 + sde**2)

      th_ei(1) = th_sxen*sxen_ei(1) + th_sde*sde_ei(1)
      th_ei(2) = th_sxen*sxen_ei(2) + th_sde*sde_ei(2)
      th_ei(3) = th_sxen*sxen_ei(3) + th_sde*sde_ei(3)

      resei       = (th - th0) * ww
      resei_ei(1) =  th_ei(1)  * ww
      resei_ei(2) =  th_ei(2)  * ww
      resei_ei(3) =  th_ei(3)  * ww

      aci = ww

c---------------------------------------------------
c---- accumulate area-integrated residuals and Jacobians
        rese(ires) = 
     &  rese(ires)  +  resei * awt

        do kk = 1, 3
          rese_e1j(ires,kk) = 
     &    rese_e1j(ires,kk) + resei_ei(kk) * awt
        enddo

c---- accumulate weighted cell area
      ac(ires) = ac(ires)  + aci * awt 
c
 100  continue ! next weighting function Wi
c--------------------------------------------------
c
 500  continue ! next Gauss point

      return
      end ! hsmre1
