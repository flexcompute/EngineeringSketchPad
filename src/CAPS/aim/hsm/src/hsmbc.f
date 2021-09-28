
      subroutine hsmfbc(lrcurv, ldrill,
     &                  pare    ,
     &                  par1    , par2    ,
     &                  var1    , var2    ,
     &             res, res_var1, res_var2,
     &             dres )
c------------------------------------------------------------
c     Sets up shell-edge force,moment loading BC residual 
c     terms for a shell element edge with two end nodes.
c
c Inputs:
c   pare(.)   boundary edge parameters
c   par#(.)   node parameters
c   var#(.)   node primary variables
c   
c
c Outputs:
c   res(..)   residuals    int [ R ] Wi dl  ,  i = 1, 2
c      ( 1)   x-force 
c      ( 2)   y-force 
c      ( 3)   z-force 
c      ( 4)   x-moment
c      ( 5)   y-moment
c      ( 6)   z-moment
c
c      res_vars(...)     d(res)/d(var)
c
c      dres(.)  length integrals of weight function,  int [ 1 ] Wi dl
c------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c---- gauss-integration locations and weights
      include 'gaussb.inc'

c----------------------------------------------
c---- passed variables
      logical lrcurv, ldrill
      real pare(letot)
      real par1(lvtot),
     &     par2(lvtot)
      real var1(ivtot),
     &     var2(ivtot)

      real res(irtot,2),
     &     res_var1(irtot,ivtot,2),
     &     res_var2(irtot,ivtot,2)
      real dres(2)

c----------------------------------------------
c---- local arrays
      parameter (ngdim = ngaussb)

      real pars(lvtot,2),
     &     vars(ivtot,2)

      real fxyzj(3,2),
     &     ftlnj(3,2),
     &     mxyzj(3,2),
     &     mtlnj(3,2)

      real r0j(3,2), rj(3,2)
      real d0j(3,2), dj(3,2)
      real p0j(2), pj(2)


      real r0(3,ngdim)
      real r(3,ngdim),
     &     r_rj(3,3,2,ngdim), r_dj(3,3,2,ngdim), r_pj(3,2,ngdim)

      real d0(3,ngdim)
      real d(3,ngdim),
     &     d_rj(3,3,2,ngdim), d_dj(3,3,2,ngdim), d_pj(3,2,ngdim)

      real a01(3,ngdim)
      real a1(3,ngdim),
     &     a1_rj(3,3,2,ngdim), a1_dj(3,3,2,ngdim), a1_pj(3,2,ngdim)

      real g0(ngdim)
      real g(ngdim),
     &     g_rj(3,2,ngdim), g_dj(3,2,ngdim), g_pj(2,ngdim)

      real h0(ngdim)
      real h(ngdim),
     &     h_rj(3,2,ngdim), h_dj(3,2,ngdim), h_pj(2,ngdim)

      real l0(ngdim)
      real l(ngdim),
     &     l_rj(3,2,ngdim), l_dj(3,2,ngdim), l_pj(2,ngdim)

      real th0(3,ngdim)
      real th(3,ngdim),
     &     th_rj(3,3,2,ngdim), th_dj(3,3,2,ngdim), th_pj(3,2,ngdim)

      real lh0(3,ngdim)
      real lh(3,ngdim),
     &     lh_rj(3,3,2,ngdim), lh_dj(3,3,2,ngdim), lh_pj(3,2,ngdim)


      real nn(2)
      real d1nn(2)

      real ajac0
      real ajac, ajac_rj(3,2), ajac_dj(3,2), ajac_pj(2)
      real awt, awt_rj(3,2), awt_dj(3,2), awt_pj(2)

      real fxyz(3),
     &     ftln(3),
     &     mxyz(3),
     &     mtln(3)

      real fbc(3), fbc_rj(3,3,2), fbc_dj(3,3,2), fbc_pj(3,2)
      real mbc(3), mbc_rj(3,3,2), mbc_dj(3,3,2), mbc_pj(3,2)


      real resi(irtot),
     &     resi_rj(irtot,3,2),
     &     resi_dj(irtot,3,2),
     &     resi_pj(irtot,2)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------
c---- machine zero (approximate)
      data epsm / 1.0e-9 /

c---- put passed node quantities into more convenient node-indexed arrays
      do lv = 1, lvtot
        pars(lv,1) = par1(lv)
        pars(lv,2) = par2(lv)
      enddo

      do iv = 1, ivtot
        vars(iv,1) = var1(iv)
        vars(iv,2) = var2(iv)
      enddo

c---- further unpack vectors
      do j = 1, 2
        do k = 1, 3
          r0j(k,j) = pars(lvr0x-1+k,j)
          d0j(k,j) = pars(lvn0x-1+k,j)
          rj(k,j) = vars(ivrx-1+k,j)
          dj(k,j) = vars(ivdx-1+k,j)
        enddo
        p0j(j) = 0.
        pj(j) = vars(ivps,j)
      enddo

      do k = 1, 3
        fxyzj(k,1) = pare(lef1x-1+k)
        fxyzj(k,2) = pare(lef2x-1+k)

        ftlnj(k,1) = pare(lef1t-1+k)
        ftlnj(k,2) = pare(lef2t-1+k)

        mxyzj(k,1) = pare(lem1x-1+k)
        mxyzj(k,2) = pare(lem2x-1+k)

        mtlnj(k,1) = pare(lem1t-1+k)
        mtlnj(k,2) = pare(lem2t-1+k)
      enddo

      call hsmgeo2(ngaussb, xsib, fwtb,
     &             lrcurv, ldrill,
     &            r0j(1,1), r0j(1,2),
     &            d0j(1,1), d0j(1,2),
     &            p0j(1)  , p0j(2)  ,
     &             r0,  r_rj,  r_dj,  r_pj,
     &             d0,  d_rj,  d_dj,  d_pj,
     &            a01, a1_rj, a1_dj, a1_pj,
     &            g0 ,  g_rj,  g_dj,  g_pj,
     &            h0 ,  h_rj,  h_dj,  h_pj,
     &            l0 ,  l_rj,  l_dj,  l_pj,
     &            th0, th_rj, th_dj, th_pj,
     &            lh0, lh_rj, lh_dj, lh_pj )

      call hsmgeo2(ngaussb, xsib, fwtb,
     &            lrcurv, ldrill,
     &            rj(1,1), rj(1,2),
     &            dj(1,1), dj(1,2),
     &            pj(1)  , pj(2)  ,
     &             r,  r_rj,  r_dj,  r_pj,
     &             d,  d_rj,  d_dj,  d_pj,
     &            a1, a1_rj, a1_dj, a1_pj,
     &            g ,  g_rj,  g_dj,  g_pj,
     &            h ,  h_rj,  h_dj,  h_pj,
     &            l ,  l_rj,  l_dj,  l_pj,
     &            th, th_rj, th_dj, th_pj,
     &            lh, lh_rj, lh_dj, lh_pj )


c=======================================================================
c---- clear residuals and Jacobians for summing over gauss points
      do ires = 1, 2
        do keq = 1, irtot
          res(keq,ires) = 0.
          do iv = 1, ivtot
            res_var1(keq,iv,ires) = 0.
            res_var2(keq,iv,ires) = 0.
          enddo
        enddo
        dres(ires) = 0.
      enddo

c-------------------------------------------------------------------
c---- sum over gauss points
      do 500 ig = 1, ngaussb
c
      xsi = xsib(ig)
      fwt = fwtb(ig)

c---- bilinear interpolation weights "N" at gauss point
      nn(1) = 0.5*(1.0 - xsi)
      nn(2) = 0.5*(1.0 + xsi)

c---- derivatives of interpolation weights
      d1nn(1) = -0.5
      d1nn(2) =  0.5

c-------------------------------------------------------------------------
c---- r derivatives wrt element coordinates, dr/dxsi, dr/deta

c---- interpolate edge force and moment from nodes
      do k = 1, 3
        fxyz(k) = fxyzj(k,1)*nn(1)
     &          + fxyzj(k,2)*nn(2)
        ftln(k) = ftlnj(k,1)*nn(1)
     &          + ftlnj(k,2)*nn(2)
        mxyz(k) = mxyzj(k,1)*nn(1)
     &          + mxyzj(k,2)*nn(2)
        mtln(k) = mtlnj(k,1)*nn(1)
     &          + mtlnj(k,2)*nn(2)
      enddo

c---- total edge force
      do k = 1, 3
        fbc(k) = fxyz(k) 
     &         + ftln(1)*th(k,ig)
     &         + ftln(2)*lh(k,ig)
     &         + ftln(3)* d(k,ig)
        do jj = 1, 2
          do kk = 1, 3
            fbc_rj(k,kk,jj) =
     &           ftln(1)*th_rj(k,kk,jj,ig)
     &         + ftln(2)*lh_rj(k,kk,jj,ig)
     &         + ftln(3)* d_rj(k,kk,jj,ig)
            fbc_dj(k,kk,jj) =
     &           ftln(1)*th_dj(k,kk,jj,ig)
     &         + ftln(2)*lh_dj(k,kk,jj,ig)
     &         + ftln(3)* d_dj(k,kk,jj,ig)
          enddo
          fbc_pj(k,jj) =
     &           ftln(1)*th_pj(k,jj,ig)
     &         + ftln(2)*lh_pj(k,jj,ig)
     &         + ftln(3)* d_pj(k,jj,ig)
        enddo
      enddo

c---- total edge moment
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        mbc(k) = mxyz(k)
     &         + mtln(1)*th(k,ig)
     &         + mtln(2)*lh(k,ig)
     &         + mtln(3)* d(k,ig)
        do jj = 1, 2
          do kk = 1, 3
            mbc_rj(k,kk,jj) =
     &           mtln(1)*th_rj(k,kk,jj,ig)
     &         + mtln(2)*lh_rj(k,kk,jj,ig)
     &         + mtln(3)* d_rj(k,kk,jj,ig)
            mbc_dj(k,kk,jj) =
     &           mtln(1)*th_dj(k,kk,jj,ig)
     &         + mtln(2)*lh_dj(k,kk,jj,ig)
     &         + mtln(3)* d_dj(k,kk,jj,ig)
          enddo
          mbc_pj(k,jj) =
     &           mtln(1)*th_pj(k,jj,ig)
     &         + mtln(2)*lh_pj(k,jj,ig)
     &         + mtln(3)* d_pj(k,jj,ig)
        enddo
      enddo

c---- coordinate Jacobian J
      ajac0 = sqrt(g0(ig))

c      ajac = sqrt(g(ig))
c      ajac_g = 0.5/sqrt(g(ig))
c      do jj = 1, 2
c        do kk = 1, 3
c          ajac_rj(kk,jj) = ajac_g*g_rj(kk,jj,ig)
c          ajac_dj(kk,jj) = ajac_g*g_dj(kk,jj,ig)
c        enddo
c        ajac_pj(jj) = ajac_g*g_pj(jj,ig)
c      enddo

c---- residual length-weight for this gauss point
      awt = fwt * ajac0
      do jj = 1, 2
        do kk = 1, 3
          awt_rj(kk,jj) = 0.
          awt_dj(kk,jj) = 0.
        enddo
        awt_pj(jj) = 0.
      enddo

c      awt = fwt * ajac
c      do jj = 1, 2
c        do kk = 1, 3
c          awt_rj(kk,jj) = fwt * ajac_rj(kk,jj)
c          awt_dj(kk,jj) = fwt * ajac_dj(kk,jj)
c        enddo
c        awt_pj(jj) = fwt * ajac_pj(jj)
c      enddo

c--------------------------------------------------
c---- go over the two residual weight functions
      do 100 ires = 1, 2

c---- clear working arrays for this weight function
      do keq = 1, irtot
        resi(keq) = 0.
        do jj = 1, 2
          do kk = 1, 3
            resi_rj(keq,kk,jj) = 0.
            resi_dj(keq,kk,jj) = 0.
          enddo
          resi_pj(keq,jj) = 0.
        enddo
      enddo

c---- weight function  Wi = Ni  at this gauss point
      ww = nn(ires)

c-----------------------------------------------------------------------
c---- edge force 
      do k = 1, 3
        keq = k
        resi(keq) = fbc(k) * ww
        do jj = 1, 2
          do kk = 1, 3
            resi_rj(keq,kk,jj) = fbc_rj(k,kk,jj) * ww
            resi_dj(keq,kk,jj) = fbc_dj(k,kk,jj) * ww
          enddo
          resi_pj(keq,jj) = fbc_pj(k,jj) * ww
        enddo
      enddo

c-----------------------------------------------------------------------
c---- edge moment
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        keq = 3 + k
        resi(keq) = ( mbc(k)
     &              + (r(ic,ig)-rj(ic,ires))*fbc(jc)
     &              - (r(jc,ig)-rj(jc,ires))*fbc(ic) ) * ww
        do jj = 1, 2
          do kk = 1, 3
            resi_rj(keq,kk,jj) = ( mbc_rj(k,kk,jj)
     &              + (r(ic,ig)-rj(ic,ires))*fbc_rj(jc,kk,jj)
     &              - (r(jc,ig)-rj(jc,ires))*fbc_rj(ic,kk,jj)
     &              + r_rj(ic,kk,jj,ig)*fbc(jc)
     &              - r_rj(jc,kk,jj,ig)*fbc(ic) ) * ww
            resi_dj(keq,kk,jj) = ( mbc_dj(k,kk,jj)
     &              + (r(ic,ig)-rj(ic,ires))*fbc_dj(jc,kk,jj)
     &              - (r(jc,ig)-rj(jc,ires))*fbc_dj(ic,kk,jj)
     &              + r_dj(ic,kk,jj,ig)*fbc(jc)
     &              - r_dj(jc,kk,jj,ig)*fbc(ic) ) * ww
          enddo
          resi_pj(keq,jj) = ( mbc_pj(k,jj)
     &              + (r(ic,ig)-rj(ic,ires))*fbc_pj(jc,jj)
     &              - (r(jc,ig)-rj(jc,ires))*fbc_pj(ic,jj)
     &              + r_pj(ic,jj,ig)*fbc(jc)
     &              - r_pj(jc,jj,ig)*fbc(ic) ) * ww
        enddo
        resi_rj(keq,ic,ires) = resi_rj(keq,ic,ires) - fbc(jc) * ww
        resi_rj(keq,jc,ires) = resi_rj(keq,jc,ires) + fbc(ic) * ww
      enddo

c---------------------------------------------------
c---- accumulate length-integrated residuals and Jacobians
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

          res_var1(keq,ivd,ires) = 
     &    res_var1(keq,ivd,ires) + resi_dj(keq,kk,1) * awt
     &                           + resi(keq) * awt_dj(kk,1)

          res_var2(keq,ivd,ires) = 
     &    res_var2(keq,ivd,ires) + resi_dj(keq,kk,2) * awt
     &                           + resi(keq) * awt_dj(kk,2)

        enddo

        res_var1(keq,ivps,ires) = 
     &  res_var1(keq,ivps,ires) + resi_pj(keq,1) * awt
     &                          + resi(keq) * awt_pj(1)

        res_var2(keq,ivps,ires) = 
     &  res_var2(keq,ivps,ires) + resi_pj(keq,2) * awt
     &                          + resi(keq) * awt_pj(2)

      enddo ! next equation keq


c---- integrate edge length
      dres(ires) = dres(ires) + ww * awt

 100  continue ! next weighting function Wi
c--------------------------------------------------

 500  continue ! next Gauss point

      return
      end ! hsmfbc



      subroutine hsmpbc(parp,
     &                  pars,
     &                  vars,
     &                  res, res_vars )
c------------------------------------------------------------
c     Sets up point loads on a node
c
c Inputs:
c   parp(.)    node BC or loading parameters
c   pars(.)    shell node parameters
c   vars(.)    shell node primary variables
c   
c
c Outputs:
c   res(.)   residuals    
c      (1)   x-force
c      (2)   y-force
c      (3)   z-force
c      (4)   x-moment
c      (5)   y-moment
c      (6)   z-moment
c
c      res_vars(..)     d(res)/d(var)
c
c------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c----------------------------------------------
c---- passed variables
      real parp(lptot)

      real pars(lvtot)
      real vars(ivtot)

      real res(irtot),
     &     res_vars(irtot,ivtot)

c----------------------------------------------
c---- local arrays
      real f(3),
     &     m(3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------
c---- machine zero (approximate)
      data epsm / 1.0e-9 /

c---- unpack specified force and moment load vectors
      do k = 1, 3
        f(k) = parp(lpfx+k-1)
        m(k) = parp(lpmx+k-1)
      enddo

c=======================================================================
c---- add force to force residuals
      do k = 1, 3
        keq = k
        res(keq) = res(keq) + f(k)
      enddo

c---- add moment to moment residuals
      do k = 1, 3
        keq = k + 3
        res(keq) = res(keq) + m(k)
      enddo

      return
      end ! hsmpbc


      subroutine hsmrbc(parp,
     &                  pars,
     &                  vars,
     &                  res, res_vars )
c------------------------------------------------------------
c     Sets up position,angle BC residuals for a shell node.
c
c Inputs:
c   parp(.)   node BC or loading parameters
c   pars(.)   shell node parameters
c   vars(.)   shell node primary variables
c   
c
c Outputs:
c   res(.)   residuals    
c      (1)   x-position
c      (2)   y-position
c      (3)   z-position
c      (4)   x-normal
c      (5)   y-normal
c      (6)   z-normal
c      (7)   psi
c
c      res_vars(..)     d(res)/d(var)
c
c------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c----------------------------------------------
c---- passed variables
      real parp(lptot)
      real pars(lvtot)
      real vars(ivtot)
      real res(ivtot),
     &     res_vars(ivtot,ivtot)

c----------------------------------------------
c---- local arrays
      real r(3),
     &     d(3)

      real rbc(3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------
c---- machine zero (approximate)
      data epsm / 1.0e-9 /

c---- unpack vectors
      do k = 1, 3
        r(k) = vars(ivrx-1+k)
        d(k) = vars(ivdx-1+k)
      enddo
      ps = vars(ivps)

      do k = 1, 3
        rbc(k) = parp(lprx-1+k)
      enddo

c=======================================================================
c---- clear residuals and Jacobians
      do keq = 1, ivtot
        res(keq) = 0.
        do iv = 1, ivtot
          res_vars(keq,iv) = 0.
        enddo
      enddo

c-----------------------------------------------------------------------
c---- node position
      do k = 1, 3
        keq = k
        iv  = ivrx-1 + k
        res(keq) = r(k) - rbc(k)
        res_vars(keq,iv) = 1.0
      enddo

c-----------------------------------------------------------------------
c---- node director vector
      do k = 1, 3
        keq = 3 + k
        ivd  = ivdx-1 + k
        res(keq) = d(k)
        res_vars(keq,ivd) = 1.0
      enddo

c---- drilling angle
      keq = 7
      res(keq) = ps
      res_vars(keq,ivps) = 1.0

      return
      end ! hsmrbc


      subroutine hsmjbc(par1 , par2 ,
     &                  var1 , var2 ,
     &                  res, res_var1, res_var2 )
c------------------------------------------------------------
c     Sets up shell-edge geometry Dirichlet residuals
c     for a joint point
c
c Inputs:
c   pars(.)    node parameters
c   vars(.)    node primary variables
c   
c
c Outputs:
c   res(.)   residuals    
c      (1)   x-position
c      (2)   y-position
c      (3)   z-position
c      (4)   x-rotation
c      (5)   y-rotation
c      (6)   z-rotation
c
c      res_vars(..)     d(res)/d(var)
c
c------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c----------------------------------------------
c---- passed variables
      real par1(lvtot), par2(lvtot)
      real var1(ivtot), var2(ivtot)

      real res(irtot),
     &     res_var1(irtot,ivtot), 
     &     res_var2(irtot,ivtot)

c----------------------------------------------
c---- local arrays
      real r1(3), r2(3),
     &     d1(3), d2(3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

c----------------------------------------------
c---- machine zero (approximate)
      data epsm / 1.0e-9 /

c---- unpack vectors
      do k = 1, 3
        r1(k) = var1(ivrx-1+k)
        r2(k) = var2(ivrx-1+k)
        d1(k) = var1(ivdx-1+k)
        d2(k) = var2(ivdx-1+k)
      enddo
      ps1 = var1(ivps)
      ps2 = var2(ivps)

c=======================================================================
c---- clear residuals and Jacobians
      do keq = 1, irtot
        res(keq) = 0.
        do iv = 1, ivtot
          res_var1(keq,iv) = 0.
          res_var2(keq,iv) = 0.
        enddo
      enddo

c-----------------------------------------------------------------------
c---- node position
      do k = 1, 3
        keq = k
        iv  = ivrx-1 + k
        res(keq) = r1(k) - r2(k)
        res_var1(keq,iv) =  1.0
        res_var2(keq,iv) = -1.0
      enddo

c-----------------------------------------------------------------------
c---- node virtual rotation increment
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        keq = 3 + k

c       res(keq) = d1(ic)*dd1(jc) - d1(jc)*dd1(ic) + d1(k)*dpsi1
c    &           - d2(ic)*dd2(jc) + d2(jc)*dd2(ic) - d2(k)*dpsi2

        res_var1(keq,ivdx-1+jc) =  d1(ic)
        res_var1(keq,ivdx-1+ic) = -d1(jc)
        res_var1(keq,ivps)      =  d1(k)

        res_var2(keq,ivdx-1+jc) = -d2(ic)
        res_var2(keq,ivdx-1+ic) =  d2(jc)
        res_var2(keq,ivps)      = -d2(k)
      enddo

      return
      end ! hsmjbc
