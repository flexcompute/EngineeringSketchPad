      subroutine hsmdep(leinit, lprint,
     &                  lrcurv, ldrill,
     &                  itmax, 
     &                  elim, etol, edel,
     &                  nnode,par,var,dep,
     &                  nelem,kelem,
     & kdim,ldim,nedim,nddim,nmdim,
     & acn,
     & resn, 
     & rese, rese_de,
     & rest, rest_t,
     & kdt, ndt,
     & ifrstt, ilastt, mfrstt,
     & amatt,
     & resv, resv_v,
     & kdv, ndv,
     & ifrstv, ilastv, mfrstv,
     & amatv )

c--------------------------------------------------------------------
c     Computes HSM secondary node variables in array dep.
c     All data arrays are referenced via pointers defined in index.inc
c
c  Inputs:
c      leinit = T  e vectors in dep(.k) are assumed to be initialized
c             = F  e vectors in dep(.k) are assumed to be not initialized
c
c      lprint = T  print out Newton convergence info
c             = F  work quietly
c
c      itmax     max number of Newton iterations for e1 vectors
c      dtol      convergence tolerance on e1 vectors
c
c      nnode     number of nodes k = 1..nnode
c      par(.k)   node parameters
c      var(.k)   node primary variables
c
c      nelem      number of elements n = 1..nelem
c      kelem(.n)  nodes k of each element
c
c  Outputs:
c      dep(.k)   node primary variables
c
c--------------------------------------------------------------------
      include 'index.inc'

c--------------------------------------------------------------------
c---- passed variables and arrays
c-   (last dimensions currently set to physical sizes to enable OOB checking)
      logical leinit, lprint
      logical lrcurv, ldrill

      real par(lvtot,kdim),
     &     var(ivtot,kdim),
     &     dep(jvtot,kdim)

      integer kelem(4,nedim)

c--------------------------------------------------------------------
c---- passed-in work arrays
      parameter (nbdimt = 3,
     &           nbdimv = 2 )

c---- residuals and Jacobians for e1 equations
      real resn(3,kdim), acn(kdim)
      real rese(kdim)  , rese_de(kdim)

      real rest(nbdimt,kdim,4)
      real rest_t(nbdimt,nbdimt,nddim,kdim)
      integer kdt(nddim,kdim)
      integer ndt(kdim)

c---- tensor system matrix arrays
      integer ifrstt(kdim),
     &        ilastt(kdim),
     &        mfrstt(kdim+1)
      real amatt(nbdimt,nbdimt,nmdim)

c---- vector system matrix arrays
      real resv(nbdimv,kdim,2)
      real resv_v(nbdimv,nbdimv,nddim,kdim)
      integer kdv(nddim,kdim)
      integer ndv(kdim)

c---- Newton system matrix arrays
      integer ifrstv(kdim),
     &        ilastv(kdim),
     &        mfrstv(kdim+1)
      real amatv(nbdimv,nbdimv,nmdim)


c--------------------------------------------------------------------
c---- small local work arrays
      real e1(3), e2(3), en(3)

c---- arrays for calling HSMREN
      real resni(3,4), resni_eni(3,4,3,4)
      real aci(4)

c---- arrays for calling HSMRE1
      real resei(4), resei_e1i(4,3)

c---- arrays for calling HSMRFM
      real reseps(3,4),
     &     reskap(3,4),
     &     resfor(3,4),
     &     resmom(3,4)
      real rest_ti(3,4,3,4),
     &    resti_ti(3,  3,4)

      real resfn(2,4),
     &     resga(2,4)
      real resv_vi(2,4,2,4),
     &    resvi_vi(2,  2,4)

      real acf(4)

      real de1(3,kdim)

      logical lconv


      if(.not. leinit) then
c----- initialize e1, e2 vectors, ensuring orthonormality with n
       do k = 1, nnode
         en(1) = par(lvn0x,k)
         en(2) = par(lvn0y,k)
         en(3) = par(lvn0z,k)

         e1(1) = par(lve01x,k)
         e1(2) = par(lve01y,k)
         e1(3) = par(lve01z,k)

         edn = e1(1)*en(1) + e1(2)*en(2) + e1(3)*en(3)

         if(abs(edn) .gt. 0.99) then
c-------- initial e1 almost parallel to n... 
c-        so instead, set initial e1 along smallest n component
          if    (abs(en(1)) .lt. max( abs(en(2)) ,  abs(en(3)) )) then
           e1(1) = 1.0
           e1(2) = 0.
           e1(3) = 0.
          elseif(abs(en(2)) .lt. max( abs(en(1)) ,  abs(en(3)) )) then
           e1(1) = 0.
           e1(2) = 1.0
           e1(3) = 0.
          else
           e1(1) = 0.
           e1(2) = 0.
           e1(3) = 1.0
          endif
          edn = e1(1)*en(1) + e1(2)*en(2) + e1(3)*en(3)
         endif

         e1(1) = e1(1) - en(1)*edn
         e1(2) = e1(2) - en(2)*edn
         e1(3) = e1(3) - en(3)*edn

         e1a = sqrt(e1(1)**2 + e1(2)**2 + e1(3)**2)
         e1(1) = e1(1)/e1a
         e1(2) = e1(2)/e1a
         e1(3) = e1(3)/e1a

         call cross(en,e1,e2)

         dep(jvnx,k) = en(1)
         dep(jvny,k) = en(2)
         dep(jvnz,k) = en(3)

         dep(jve1x,k) = e1(1)
         dep(jve1y,k) = e1(2)
         dep(jve1z,k) = e1(3)

         dep(jve2x,k) = e2(1)
         dep(jve2y,k) = e2(2)
         dep(jve2z,k) = e2(3)
       enddo
      endif


c=============================================================
c---- Newton iteration loop for n vectors
      do 100 iter = 1, itmax

c---- clear residuals and Jacobians for all nodes
      do k = 1, nnode
        resn(1,k) = 0.
        resn(2,k) = 0.
        resn(3,k) = 0.
        acn(k) = 0.
      enddo

c=============================================================
c---- go over all elements, accumulating residuals 
c-     and Jacobians onto each node in each element.
c
c-    The n equation system is diagonal, so no system solve is needed
c
      do 10 n = 1, nelem
c------ nodes of this element
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

        if(k4 .eq. 0) then
c------- tri element
         neln = 3
         k4 = k3   ! use node 3 as dummy node 4
        else
c------- quad element
         neln = 4
        endif

c------ calculate this element's contribution to equation residuals and Jacobian
        call hsmren(neln, lrcurv, ldrill,
     &        par(1,k1),   par(1,k2),   par(1,k3),   par(1,k4), 
     &        var(1,k1),   var(1,k2),   var(1,k3),   var(1,k4), 
     &     dep(jvnx,k1),dep(jvnx,k2),dep(jvnx,k3),dep(jvnx,k4),
     &    resni, resni_eni, aci )

c------ accumulate residuals, and integrated cell areas (Jacobians)
        do ires = 1, neln
          k = kelem(ires,n)
          do kk = 1, 3
            resn(kk,k) = resn(kk,k) + resni(kk,ires)
          enddo
          acn(k) = acn(k) + aci(ires)
c         if(k.eq.1) write(*,*) k, ires, aci(ires), acn(k)
c         if(k.eq.2) write(*,*) k, ires, aci(ires), acn(k)
        enddo
 10   continue    ! next element

      rlx = 1.0

c@@@@@@@@@@@@@@
c     k=1
c     write(*,*) 'k A R123', k, acn(k), resn(1,k), resn(2,k), resn(3,k)
c     k=2 
c     write(*,*) 'k A R123', k, acn(k), resn(1,k), resn(2,k), resn(3,k)
c     pause

c---- perform Newton update 
c-   (this is a linear problem, so should converge in one step)
      edel = 0.
      kmax = 0
      do k = 1, nnode
        dn1 = -resn(1,k) / acn(k)
        dn2 = -resn(2,k) / acn(k)
        dn3 = -resn(3,k) / acn(k)

        dep(jvnx,k) = dep(jvnx,k) + dn1
        dep(jvny,k) = dep(jvny,k) + dn2
        dep(jvnz,k) = dep(jvnz,k) + dn3

        dn = sqrt(dn1**2 + dn2**2 + dn3**2)
        if(abs(dn) .gt. edel) then
         edel = dn
         kmax = k
        endif
      enddo

      if(lprint) then
c----- print out Newton iteration changes
       if(iter.eq.1) then
        write(*,*)
        write(*,*) 'Converging n vectors ...'
       endif
       write(*,1100) iter, edel, rlx, kmax
 1100  format(1x, 'it =',i3, 
     &        '    dn =',e11.3,
     &        '   rlx =', f6.3, '  (', i4,' )')
      endif

c---- set convergence tolerance flag
      lconv = edel .lt. etol

c---- jump out of Newton loop if converged on last iteration or limit reached
      if(lconv .or. iter.ge.itmax) go to 110

c      call qpause(ians)
c      if(ians.eq.-1) return

 100  continue  ! next Newton iteration

c---- pick up here after the Newton loop
 110  continue


c=============================================================
c---- Newton iteration loop for e1,e2 vectors
      do 200 iter = 1, itmax

c---- clear residuals and Jacobians for all nodes
      do k = 1, nnode
        rese(k) = 0.
        rese_de(k) = 0.
      enddo

c=============================================================
c---- go over all elements, accumulating residuals 
c-     and Jacobians onto each node in each element.
c
c-    The e1 equation system is diagonal, so no system solve is needed
c
      do 20 n = 1, nelem
c------ nodes of this element
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

        if(k4 .eq. 0) then
c------- tri element
         neln = 3
         k4 = k3   ! use node 3 as dummy node 4
        else
c------- quad element
         neln = 4
        endif

c------ calculate this element's contribution to equation residuals and Jacobian
        call hsmre1(neln, lrcurv, ldrill,
     &        par(1,k1),    par(1,k2),    par(1,k3),    par(1,k4), 
     &        var(1,k1),    var(1,k2),    var(1,k3),    var(1,k4), 
     &        dep(1,k1),    dep(1,k2),    dep(1,k3),    dep(1,k4), 
     &    dep(jve1x,k1),dep(jve1x,k2),dep(jve1x,k3),dep(jve1x,k4),
     &    resei, resei_e1i, aci )

c------ accumulate residuals and Jacobians
        do ires = 1, neln
          k = kelem(ires,n)
          rese(k)    = rese(k)
     &               + resei(ires)
c-------- accumulate projected dRi/de1i Jacobians onto e2 of node i
          rese_de(k) = rese_de(k)
     &               + resei_e1i(ires,1)*dep(jve2x,k)
     &               + resei_e1i(ires,2)*dep(jve2y,k)
     &               + resei_e1i(ires,3)*dep(jve2z,k)
        enddo
 20   continue    ! next element

      edel = 0.
      kmax = 0
      do k = 1, nnode
c------ solve for Newton delta (e1 change along e2)
        de = -rese(k) / rese_de(k)

c------ undo projections to get cartesian e1 Newton deltas
        de1(1,k) = de*dep(jve2x,k)
        de1(2,k) = de*dep(jve2y,k)
        de1(3,k) = de*dep(jve2z,k)

        if(abs(de) .gt. edel) then
         edel = abs(de)
         kmax = k
        endif
      enddo

      rlx = 1.0
      if(edel*rlx .gt. elim) rlx = elim/edel

c---- perform Newton update
      do k = 1, nnode
        en(1) = dep(jvnx,k)
        en(2) = dep(jvny,k)
        en(3) = dep(jvnz,k)

        e1(1) = dep(jve1x,k) + rlx*de1(1,k)
        e1(2) = dep(jve1y,k) + rlx*de1(2,k)
        e1(3) = dep(jve1z,k) + rlx*de1(3,k)

c------ project out any component along n (might be due to bit-creep)
        edn = e1(1)*en(1) + e1(2)*en(2) + e1(3)*en(3)
        e1(1) = e1(1) - en(1)*edn
        e1(2) = e1(2) - en(2)*edn
        e1(3) = e1(3) - en(3)*edn

c------ normalize
        sqee = sqrt(e1(1)**2 + e1(2)**2 + e1(3)**2)
        e1(1) = e1(1)/sqee
        e1(2) = e1(2)/sqee
        e1(3) = e1(3)/sqee

c------ store updated e1
        dep(jve1x,k) = e1(1)
        dep(jve1y,k) = e1(2)
        dep(jve1z,k) = e1(3)

c------ store corresponding e2 = n x e1
        dep(jve2x,k) = en(2)*e1(3) - en(3)*e1(2)
        dep(jve2y,k) = en(3)*e1(1) - en(1)*e1(3)
        dep(jve2z,k) = en(1)*e1(2) - en(2)*e1(1)
      enddo

      if(lprint) then
c----- print out Newton iteration changes
       if(iter.eq.1) then
        write(*,*)
        write(*,*) 'Converging e1,e2 vectors ...'
       endif
       write(*,2100) iter, edel, rlx, kmax
 2100  format(1x, 'it =',i3, 
     &        '    de =',e11.3,
     &        '   rlx =', f6.3, '  (', i4,' )')
      endif

c---- set convergence tolerance flag
      lconv = edel .lt. etol

c---- jump out of Newton loop if converged on last iteration or limit reached
      if(lconv .or. iter.ge.itmax) go to 210

c      call qpause(ians)
c      if(ians.eq.-1) return

 200  continue  ! next Newton iteration

c---- pick up here after the Newton loop
 210  continue

c==========================================================================
c---- solve for all node tensor and vector values
c
      neqt = 3   ! 3 tensor equations for t11, t22, t12
      neqv = 2   ! 2 vector equations for v1, v2

      if(lprint) then
       write(*,*)
       write(*,*) 'Calculating  strains, stress resultants ...'
      endif

c---- clear tensor residuals rest, vector residuals resv for all nodes
      do k = 1, nnode
        do lt = 1, 4
          do ieq = 1, neqt
            rest(ieq,k,lt) = 0.0
          enddo
        enddo

        do lv = 1, 2
          do ieq = 1, neqv
            resv(ieq,k,lv) = 0.0
          enddo
        enddo

c------ also clear Jacobian node indices for both types of residual
        ndt(k) = 0
        ndv(k) = 0
        do id = 1, nddim
          kdv(id,k) = 0
          kdt(id,k) = 0
        enddo
      enddo

c---- clear all tensors and vectors to be determined
c-    (their equations are linear, so initial guess is arbitrary)
      do k = 1, nnode
        dep(jvf11,k) = 0.
        dep(jvf22,k) = 0.
        dep(jvf12,k) = 0.

        dep(jvm11,k) = 0.
        dep(jvm22,k) = 0.
        dep(jvm12,k) = 0.

        dep(jve11,k) = 0.
        dep(jve22,k) = 0.
        dep(jve12,k) = 0.

        dep(jvk11,k) = 0.
        dep(jvk22,k) = 0.
        dep(jvk12,k) = 0.

        dep(jvf1n,k) = 0.
        dep(jvf2n,k) = 0.

        dep(jvga1,k) = 0.
        dep(jvga2,k) = 0.
      enddo

c=============================================================
c---- go over all elements, accumulating residuals 
c-     and Jacobians onto each node of each element.
      do 300 n = 1, nelem

c------ nodes of this element
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

        if(k4 .eq. 0) then
c------- tri element
         neln = 3
         k4 = k3
        else
c------- quad element
         neln = 4
        endif

c------ calculate this element's contribution to equation residuals and Jacobian
        call hsmrfm(neln, lrcurv, ldrill,
     &        par(1,k1), par(1,k2), par(1,k3), par(1,k4), 
     &        var(1,k1), var(1,k2), var(1,k3), var(1,k4), 
     &        dep(1,k1), dep(1,k2), dep(1,k3), dep(1,k4), 
     &  reseps,
     &  reskap,
     &  resfor,
     &  resmom, rest_ti,
     &  resfn,
     &  resga, resv_vi, acf )

c------ project and accumulate residuals and Jacobians for this element,
c-      for each residulal weighting function  Wi 
        do 30 ires = 1, neln

        kres = kelem(ires,n)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ accumulate to overall tensor system
        do ieq = 1, neqt
          rest(ieq,kres,1) = rest(ieq,kres,1) + reseps(ieq,ires)
          rest(ieq,kres,2) = rest(ieq,kres,2) + reskap(ieq,ires)
          rest(ieq,kres,3) = rest(ieq,kres,3) + resfor(ieq,ires)
          rest(ieq,kres,4) = rest(ieq,kres,4) + resmom(ieq,ires)
        enddo

c------ set Jacobians for only this ires, for passing to JADD
        do j = 1, 4
          do ieq = 1, neqt
            resti_ti(ieq,  1,j) = rest_ti(ieq,ires,1,j)
            resti_ti(ieq,  2,j) = rest_ti(ieq,ires,2,j)
            resti_ti(ieq,  3,j) = rest_ti(ieq,ires,3,j)
          enddo
        enddo


        call jadd(1,neqt, neqt, nbdimt, nddim,
     &         resti_ti(1,1,1), k1,
     &           rest_t(1,1,1,kres), kdt(1,kres), ndt(kres))
        call jadd(1,neqt, neqt, nbdimt, nddim,
     &         resti_ti(1,1,2), k2,
     &           rest_t(1,1,1,kres), kdt(1,kres), ndt(kres))
        call jadd(1,neqt, neqt, nbdimt, nddim,
     &         resti_ti(1,1,3), k3,
     &           rest_t(1,1,1,kres), kdt(1,kres), ndt(kres))
        if(neln.eq.4) then
        call jadd(1,neqt, neqt, nbdimt, nddim,
     &         resti_ti(1,1,4), k4,
     &           rest_t(1,1,1,kres), kdt(1,kres), ndt(kres))
        endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ accumulate to overall vector system
        do ieq = 1, neqv
          resv(ieq,kres,1) = resv(ieq,kres,1) + resfn(ieq,ires)
          resv(ieq,kres,2) = resv(ieq,kres,2) + resga(ieq,ires)
        enddo

c------ set Jacobians for this ires, for passing to JADD
        do j = 1, neln
          do ieq = 1, neqv
            resvi_vi(ieq,  1,j) = resv_vi(ieq,ires,1,j)
            resvi_vi(ieq,  2,j) = resv_vi(ieq,ires,2,j)
          enddo
        enddo

        call jadd(1,neqv, neqv, nbdimv, nddim,
     &         resvi_vi(1,1,1), k1,
     &           resv_v(1,1,1,kres), kdv(1,kres), ndv(kres))
        call jadd(1,neqv, neqv, nbdimv, nddim,
     &         resvi_vi(1,1,2), k2,
     &           resv_v(1,1,1,kres), kdv(1,kres), ndv(kres))
        call jadd(1,neqv, neqv, nbdimv, nddim,
     &         resvi_vi(1,1,3), k3,
     &           resv_v(1,1,1,kres), kdv(1,kres), ndv(kres))
        if(neln.eq.4) then
        call jadd(1,neqv, neqv, nbdimv, nddim,
     &         resvi_vi(1,1,4), k4,
     &           resv_v(1,1,1,kres), kdv(1,kres), ndv(kres))
        endif

 30     continue ! next residual in this element

 300  continue  ! next element

c---- create tensor matrix in pointered block-column array amatt
      call colset(.true. ,
     &            kdim,nddim,nmdim,
     &            nbdimt,neqt, nnode, rest_t,kdt,ndt,
     &            namatt,ifrstt,ilastt,mfrstt,amatt)

c---- total number of blocks in array amatt
      ntmmat = mfrstt(namatt+1)

c---- create tensor matrix in pointered block-column array amatt
      call colset(.true. ,
     &            kdim,nddim,nmdim,
     &            nbdimv,neqv, nnode, resv_v,kdv,ndv,
     &            namatv,ifrstv,ilastv,mfrstv,amatv)

c---- total number of blocks in array amatt
      nvmmat = mfrstv(namatv+1)


      if(.false.) then
       write(*,*) 'Writing tensor block matrix to fort.70 ...'
       lu = 70
       call bmdump(lu, nnode, neqt,
     &             nbdimt,kdim,nmdim,
     &             ntmmat,ifrstt,ilastt,mfrstt,amatt)
       write(*,*) 'Writing vector block matrix to fort.60 ...'
       lu = 60
       call bmdump(lu, nnode, neqv,
     &             nbdimv,kdim,nmdim,
     &             nvmmat,ifrstv,ilastv,mfrstv,amatv)
      endif

c---------------------------------------------------------------------------
c---- LU-factor systems
      call sblud(nbdimt,neqt, namatt,ifrstt,ilastt,mfrstt, amatt)
      call sblud(nbdimv,neqv, namatv,ifrstv,ilastv,mfrstv, amatv)

c---- fwd,back substitution of residual vectors
      do lt = 1, 4
        call sbbks(nbdimt,neqt, namatt,ifrstt,ilastt,mfrstt,
     &             amatt,rest(1,1,lt))
      enddo

      do lv = 1, 2
        call sbbks(nbdimv,neqv, namatv,ifrstv,ilastv,mfrstv,
     &             amatv,resv(1,1,lv))
      enddo

c---- store solution variables in dep array
      do k = 1, nnode
        do kv = 1, neqt
          jv = jve11-1 + kv
          dep(jv,k) = dep(jv,k) - rest(kv,k,1)

          jv = jvk11-1 + kv
          dep(jv,k) = dep(jv,k) - rest(kv,k,2)

          jv = jvf11-1 + kv
          dep(jv,k) = dep(jv,k) - rest(kv,k,3)

          jv = jvm11-1 + kv
          dep(jv,k) = dep(jv,k) - rest(kv,k,4)
        enddo
      enddo

      do k = 1, nnode
        do kv = 1, neqv
          jv = jvf1n-1 + kv
          dep(jv,k) = dep(jv,k) - resv(kv,k,1)

          jv = jvga1-1 + kv
          dep(jv,k) = dep(jv,k) - resv(kv,k,2)
        enddo
      enddo


      return

      end ! hsmdep
