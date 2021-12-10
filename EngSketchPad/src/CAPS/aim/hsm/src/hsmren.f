
       subroutine hsmren(neln, lrcurv, ldrill,
     &                   par1    , par2    , par3    , par4    , 
     &                   var1    , var2    , var3    , var4    ,
     &                   enj1    , enj2    , enj3    , enj4    ,
     &              resn, resn_enj , 
     &              ac )
c------------------------------------------------------------
c     Sets up residuals and Jacobians for n vector
c     definitions for a four-node shell element.
c
c Inputs:
c   nj        number of element nodes (3 = tri, 4 = quad)
c   par#(.)   node parameters
c   var#(.)   node primary variables
c   enj#(.)   node normal  unit vectors n
c
c   Note: If neln=3, then par4(.), var4(.), enj4(.) will not affect results
c
c Outputs:
c   resn(..)   n residuals,                       int int [ Rn ] Wi/J^2 dA
c   resn_enj#(....)  d(resn)/d(enj#)
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

      real enj1(3),
     &     enj2(3),
     &     enj3(3),
     &     enj4(3)

      real resn(3,4), resn_enj(3,4,3,4)

      real ac(4)

c----------------------------------------------
c---- local arrays
      real par(lvtot,4),
     &     var(ivtot,4)
      real enj(3,4)

      real d1r0(3),
     &     d2r0(3)

      real d1r(3),
     &     d2r(3)

      real rj(3,4), r0j(3,4)
      real dj(3,4), e0nj(3,4), e01j(3,4), e02j(3,4)

      real r(3), r0(3)
      real d(3), d0(3)
      real n(3), n0(3)
      real l(2), l0(2)

      real de(3), de0(3)
      real le(4), le0(4)

      real eps(2,2), eps_rj(2,2,3,4)

      real nn(4)
      real d1nn(4), d2nn(4)

      real nnc(4)

      real rxr0(3)
      real arxr0

      real rxr(3), rxr_rj(3,3,4)
      real arxr, arxr_rj(3,4)

      real ajac0

      real resni(3),
     &     resni_rj(3,3,4),
     &     resni_enj(3,3,4)

c----------------------------------------------
c---- gauss-integration locations and weights
      include 'gausst.inc'
      include 'gaussq.inc'

      logical lslerp
      logical lacurv
      logical lnproj

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------
c---- sqrt(2)
      data sq2 / 1.41421356237309514547462185873882845 /

c----------------------------------------------
c---- approximate machine zero
c      data epsm / 1.0e-9 /
      data epsm / 1.0e-18 /
c      data epsm / 100.0 /   ! use this to check small-lambda case Jacobians

c----------------------------------------------
      lslerp = .true.   ! use SLERP interpolation for d vectors
c     lslerp = .false.  ! use bilinear interpolation for d vectors

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

      do k = 1, 3
        enj(k,1) = enj1(k)
        enj(k,2) = enj2(k)
        enj(k,3) = enj3(k)
        enj(k,4) = enj4(k)
      enddo
c
c---- further unpack vectors
      do j = 1, 4
        r0j(1,j) = par(lvr0x,j)
        r0j(2,j) = par(lvr0y,j)
        r0j(3,j) = par(lvr0z,j)

        e0nj(1,j) = par(lvn0x,j)
        e0nj(2,j) = par(lvn0y,j)
        e0nj(3,j) = par(lvn0z,j)

        rj(1,j) = var(ivrx,j)
        rj(2,j) = var(ivry,j)
        rj(3,j) = var(ivrz,j)

        dj(1,j) = var(ivdx,j)
        dj(2,j) = var(ivdy,j)
        dj(3,j) = var(ivdz,j)
      enddo

c=======================================================================
c---- clear residuals and Jacobians for summing over gauss points
      do ires = 1, 4
        do k = 1, 3
          resn(k,ires) = 0.
          do jj = 1, 4
            do kk = 1, 3
              resn_enj(k,ires,kk,jj) = 0.
            enddo
          enddo
        enddo
        ac(ires) = 0.
      enddo

c---- transverse shear interpolation

c---- go over each edge
      do j = 1, neln
c------ nodes at edge endpoints
        jo = j
        jp = mod(j,neln) + 1

c------ set d using linear interpolation
        do k = 1, 3
          de(k)  = 0.5*(dj(k,jp) + dj(k,jo))
          de0(k) = 0.5*(e0nj(k,jp) + e0nj(k,jo))
        enddo

c------ d.a1 (j=1), d.a2 (j=2), -d.a1 (j=3), -d.a2 (j=4)
        le(j) = de(1)*(rj(1,jp)-rj(1,jo))
     &        + de(2)*(rj(2,jp)-rj(2,jo))
     &        + de(3)*(rj(3,jp)-rj(3,jo))
        le0(j) = de0(1)*(r0j(1,jp)-r0j(1,jo))
     &         + de0(2)*(r0j(2,jp)-r0j(2,jo))
     &         + de0(3)*(r0j(3,jp)-r0j(3,jo))
      enddo ! next edge

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
        r0(k) = r0j(k,1)*nn(1)
     &        + r0j(k,2)*nn(2)
     &        + r0j(k,3)*nn(3)
     &        + r0j(k,4)*nn(4)
        d1r0(k) = r0j(k,1)*d1nn(1)
     &          + r0j(k,2)*d1nn(2)
     &          + r0j(k,3)*d1nn(3)
     &          + r0j(k,4)*d1nn(4)
        d2r0(k) = r0j(k,1)*d2nn(1)
     &          + r0j(k,2)*d2nn(2)
     &          + r0j(k,3)*d2nn(3)
     &          + r0j(k,4)*d2nn(4)

        r(k) = rj(k,1)*nn(1)
     &       + rj(k,2)*nn(2)
     &       + rj(k,3)*nn(3)
     &       + rj(k,4)*nn(4)
        d1r(k) = rj(k,1)*d1nn(1)
     &         + rj(k,2)*d1nn(2)
     &         + rj(k,3)*d1nn(3)
     &         + rj(k,4)*d1nn(4)
        d2r(k) = rj(k,1)*d2nn(1)
     &         + rj(k,2)*d2nn(2)
     &         + rj(k,3)*d2nn(3)
     &         + rj(k,4)*d2nn(4)

        d0(k) = e0nj(k,1)*nn(1)
     &        + e0nj(k,2)*nn(2)
     &        + e0nj(k,3)*nn(3)
     &        + e0nj(k,4)*nn(4)
        d(k) = dj(k,1)*nn(1)
     &       + dj(k,2)*nn(2)
     &       + dj(k,3)*nn(3)
     &       + dj(k,4)*nn(4)

      enddo

c---- r_xi x r_eta
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        rxr0(k) = d1r0(ic)*d2r0(jc) - d1r0(jc)*d2r0(ic)

        rxr(k) = d1r(ic)*d2r(jc) - d1r(jc)*d2r(ic)
        do jj = 1, neln
          rxr_rj(k,k ,jj) = 0.
          rxr_rj(k,ic,jj) = 
     &           d1nn(jj)*d2r(jc) - d1r(jc)*d2nn(jj)
          rxr_rj(k,jc,jj) = 
     &           d1r(ic)*d2nn(jj) - d1nn(jj)*d2r(ic)
        enddo
      enddo

      arxr0 = sqrt(rxr0(1)**2 + rxr0(2)**2 + rxr0(3)**2)
      if(arxr0 .eq. 0.0) then
       arxr0i = 1.0
      else
       arxr0i = 1.0/arxr0
      endif

      arxr  = sqrt(rxr(1)**2 + rxr(2)**2 + rxr(3)**2)
      if(arxr .eq. 0.0) then
       arxri = 1.0
      else
       arxri = 1.0/arxr
      endif

      do jj = 1, neln
        do kk = 1, 3
          arxr_rj(kk,jj) = 
     &      (  rxr_rj(1,kk,jj)*rxr(1)
     &       + rxr_rj(2,kk,jj)*rxr(2)
     &       + rxr_rj(3,kk,jj)*rxr(3) ) * arxri
        enddo
      enddo

c---- coordinate Jacobian J
      ajac0 = arxr0

c---- residual area-weight for this gauss point
cc    awt = fwt * ajac0
      awt = fwt / ajac0

c---- interpolate d.a edge values to xsi,eta
      if(neln .eq. 3) then
c----- MITC3 interpolation
       l(1) =  le(1)*(1.0-eta) + (-le(3) - sq2*le(2))*eta
       l(2) = -le(3)*(1.0-xsi) + ( le(1) + sq2*le(2))*xsi
       l0(1) =  le0(1)*(1.0-eta) + (-le0(3) - sq2*le0(2))*eta
       l0(2) = -le0(3)*(1.0-xsi) + ( le0(1) + sq2*le0(2))*xsi

      else
c----- MITC4 interpolation
       nnc(1) =  0.25*(1.0-eta)
       nnc(3) = -0.25*(1.0+eta)

       nnc(2) =  0.25*(1.0+xsi)
       nnc(4) = -0.25*(1.0-xsi)

       l(1) = le(1)*nnc(1) + le(3)*nnc(3)
       l(2) = le(2)*nnc(2) + le(4)*nnc(4)
       l0(1) = le0(1)*nnc(1) + le0(3)*nnc(3)
       l0(2) = le0(2)*nnc(2) + le0(4)*nnc(4)
      endif

c---- metric for bilinear surface
      g11 = d1r(1)*d1r(1) + d1r(2)*d1r(2) + d1r(3)*d1r(3)
      g22 = d2r(1)*d2r(1) + d2r(2)*d2r(2) + d2r(3)*d2r(3)
      g12 = d1r(1)*d2r(1) + d1r(2)*d2r(2) + d1r(3)*d2r(3)
      gd = g11*g22 - g12**2

      g011 = d1r0(1)*d1r0(1) + d1r0(2)*d1r0(2) + d1r0(3)*d1r0(3)
      g022 = d2r0(1)*d2r0(1) + d2r0(2)*d2r0(2) + d2r0(3)*d2r0(3)
      g012 = d1r0(1)*d2r0(1) + d1r0(2)*d2r0(2) + d1r0(3)*d2r0(3)
      gd0 = g011*g022 - g012**2

c---- approximate n vector using approximate gamma
      do k = 1, 3
        n(k) = d(k)
     &       - l(1)*(g22*d1r(k) - g12*d2r(k))/gd
     &       - l(2)*(g11*d2r(k) - g12*d1r(k))/gd
        n0(k) = d0(k)
     &        - l0(1)*(g022*d1r0(k) - g012*d2r0(k))/gd0
     &        - l0(2)*(g011*d2r0(k) - g012*d1r0(k))/gd0
      enddo
      an  = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
      an0 = sqrt(n0(1)**2 + n0(2)**2 + n0(3)**2)

c@@@@@@@@@@
c      if(r(1) .lt. 0.12) then
c        write(*,*) r(1), r(2), l(1), l(2), n(1), n(2), n0(1), n0(2)
c      endif

c--------------------------------------------------
c---- go over the node residual weight functions
      do 100 ires = 1, neln

c---- weight function  Wi = Ni at this gauss point
      ww = nn(ires)

      do k = 1, 3
c       resni(k) = ( enj(k,ires) -  rxr(k)*arxri ) * ww
c    &           - (e0nj(k,ires) - rxr0(k)*arxr0i) * ww

        resni(k) = ( enj(k,ires) -  n(k)/an ) * ww
     &           - (e0nj(k,ires) - n0(k)/an0) * ww

        do jj = 1, neln
          do kk = 1, 3
            resni_enj(k,kk,jj) = 0.
c            resni_rj(k,kk,jj) =
c     &        ( - rxr_rj(k,kk,jj)
c     &          + rxr(k)*arxri*arxr_rj(kk,jj) )*arxri * ww
          enddo
        enddo

        jj = ires
        resni_enj(k,k,jj) = ww

      enddo

      aci = ww

c---------------------------------------------------
c---- accumulate area-integrated residuals and Jacobians
      do k = 1, 3
        resn(k,ires) = 
     &  resn(k,ires)  +  resni(k) * awt
        do jj = 1, neln
          do kk = 1, 3
            resn_enj(k,ires,kk,jj) = 
     &      resn_enj(k,ires,kk,jj) + resni_enj(k,kk,jj) * awt
          enddo
        enddo
      enddo

c---- accumulate weighted cell area
      ac(ires) = ac(ires)  + aci * awt 

c@@@@@@@@@@
c      if(r(1) .lt. 0.12) then
c       write(*,*) ires, aci, awt, ac(ires)
c      endif

c
 100  continue ! next weighting function Wi
c--------------------------------------------------
c
 500  continue ! next Gauss point

      return
      end ! hsmren
