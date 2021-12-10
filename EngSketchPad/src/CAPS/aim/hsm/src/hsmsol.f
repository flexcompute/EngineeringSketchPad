
      subroutine hsmsol(lvinit, lprint,
     &                  lrcurv, ldrill,
     &                  itmax, rref,
     &                  rlim, rtol, rdel,
     &                  alim, atol, adel,
     &                  damem,rtolm,
     &                  parg,
     &                  nnode,pars,vars,
     &                  nvarg, varg,
     &                  nelem,kelem,
     &                  nbcedge,kbcedge,pare,
     &                  nbcnode,kbcnode,parp,lbcnode,
     &                  njoint,kjoint,
     & kdim,ldim,nedim,nddim,nmdim,
     & bf, bf_dj,
     & bm, bm_dj,
     & ibr1,ibr2,ibr3,ibd1,ibd2,ibd3,
     & resc, resc_vars,
     & resp, resp_vars, resp_dvp,
     & kdvp, ndvp,
     & ares,
     & ifrst,ilast, mfrst,
     & amat, ipp, dvars )

c--------------------------------------------------------------------
c     Version 1.00                                      21 Jul 2019
c
c     Solves HSM analysis problem.
c     All data arrays are referenced via pointers defined in index.inc
c
c  Inputs:
c  -------
c      lvinit = T  vars(.k) are assumed to be initialized
c             = F  vars(.k) are assumed to be not initialized
c
c      lprint = T  print out Newton convergence info
c             = F  work quietly
c
c      lrcurv = T  define vectors on curved surface   (4th order accurate)
c             = F  define vectors on bilinear surface (2nd order accurate)
c
c      ldrill = T  include drilling DOF p
c             = F  exclude drilling DOF p
c
c      itmax  max allowed number of Newton iterations
c             =  0  do equation setup and solve, but no Newton update
c             = -1  do equation setup only
c             = -2  return only required matrix size nmdim
c
c      rref    reference length for r position Newton changes
c
c      rlim    max allowed position Newton change  (fraction of rref)
c      rtol    convergence tolerance on position   (fraction of rref)
c
c      alim    max allowed angle Newton change (radians)
c      atol    convergence tolerance on angles
c
c      damem   angle change threshold which triggers membrane-only
c               sub-iterations
c      rtolm   convergence tolerance on position 
c               for membrane-only sub-iterations   (fraction of rref)
c
c      parg(.)   global parameters (can be all zero)
c
c      nnode     number of nodes k = 1..nnode
c      pars(.k)   node parameters
c      vars(.k)   node primary variables  (inputs only if lvinit = T)
c
c      nvarg      number of global variables  (not used currently)
c      varg(.)    global variables            (not used currently)
c
c      nelem      number of elements n = 1..nelem
c      kelem(.n)  node k indices of element n
c
c      nbcedge     number of BC edges  i = 1..nbcedge
c      kbcedge(.i) node k indices of BC edge i
c      pare(.i)    edge-BC parameters
c
c      nbcnode     number of BC nodes  i = 1..nbcnode
c      kbcnode(i)  node k of BC node i
c      parp(.i)    node-BC parameters
c      lbcnode(i)  node-BC type indicator
c
c      njoint      number of joint node pairs i = 1..njoint
c      kjoint(.i)  node k indices of the joint pair
c
c    kdim    |  work arrays and dimensions
c     .      |
c     .      |
c    dvar    |
c
c
c  Outputs:
c  --------
c      vars(.k)   node primary variables
c
c      rdel    last max position Newton change  (fraction of rref)
c      adel    last max angle Newton change
c
c               ( converged if  rdel < rtol  and  adel < atol )
c
c--------------------------------------------------------------------
c---- define mnemonic pointers for parg,pars,vars arrays
      include 'index.inc'

c--------------------------------------------------------------------
c---- passed variables and arrays
c-   (last dimensions currently set to physical sizes to enable OOB checking)
      logical lvinit, lprint
      logical lrcurv, ldrill

      integer itmax
      real rref
      real rlim, rtol, rdel
      real alim, atol, adel

      real parg(lgtot)

      integer nnode
      real pars(lvtot,kdim),
     &     vars(ivtot,kdim)

      integer nvarg
      real varg(*)

      integer nelem
      integer kelem(4,nedim)

      integer nbcedge
      integer kbcedge(2,ldim)

      integer nbcnode
      integer kbcnode(ldim)
      integer lbcnode(ldim)

      real pare(letot,ldim)
      real parp(lptot,ldim)

      integer njoint
      integer kjoint(2,ldim)

c--------------------------------------------------------------------
c---- passed-in work arrays

c---- residual and Newton delta projection vectors
      real bf(3,3,kdim),
     &     bm(3,3,kdim)
      real bf_dj(3,3,3,kdim),
     &     bm_dj(3,3,3,kdim)

c---- BC-type flags
      integer ibr1(ldim), ibr2(ldim), ibr3(ldim),
     &        ibd1(ldim), ibd2(ldim), ibd3(ldim)

c---- residuals and Jacobians in cartesian basis 
      real resc(irtot,kdim)
      real resc_vars(irtot,ivtot,nddim,kdim)

c---- residuals and Jacobians in b1,b2,bn projection basis
      real resp0(irtot,kdim)
      real dvp0(irtot,kdim)

      real resp(irtot,kdim)
      real resp_dvp(irtot,irtot,nddim,kdim)

c---- Jacobian element pointers
      integer kdvp(nddim,kdim)
      integer ndvp(kdim)

c---- work array for intermediate projection operations
      real resp_vars(irtot,ivtot,nddim)


c---- residual-weight areas
      real ares(kdim)

c---- Newton system matrix arrays
      integer ifrst(kdim), 
     &        ilast(kdim), 
     &        mfrst(kdim+1)
      real amat(irtot,irtot,nmdim)
      integer ipp(irtot,nmdim)

c---- primary variable changes
      real dvars(ivtot,kdim)

c--------------------------------------------------------------------
c---- small local work arrays
      real dxt(3)

      real dl(3), dla, dla_dl(3)

      real t0j(3), t0j_n01(3,3), t0j_n02(3,3)
      real l0j(3), l0j_n01(3,3), l0j_n02(3,3)
      real n0j(3), n0j_n01(3,3), n0j_n02(3,3)

      real tj(3), tj_d1(3,3), tj_d2(3,3)
      real lj(3), lj_d1(3,3), lj_d2(3,3)
      real nj(3), nj_d1(3,3), nj_d2(3,3)

      real da0, da0_n01(3), da0_n02(3)
      real da, da_d1(3), da_d2(3)

c---- arrays for calling HSMEQN
      real resi(irtot,4),
     &     resi_vars(irtot,ivtot,4,4)
      real aresi(4), dresi(4)

c---- arrays for calling HSMFBC
      real resf(irtot,2),
     &     resf_vars(irtot,ivtot,2,2)

c---- arrays for calling HSMFPR, HSMVPR
      real resp_bf(irtot,3,3),
     &     resp_bm(irtot,3,3),
     &     resp_dum(irtot,3)
      real resp_dj(3)

c---- arrays for calling HSMRBC
      real resr(ivtot),
     &     resr_vars(ivtot,ivtot)

c---- work arrays for joint setup
      real resd_var1(irtot,ivtot),
     &     resd_var2(irtot,ivtot)
      real resp_var1(irtot,ivtot),
     &     resp_var2(irtot,ivtot)
      real resp_dvp1(irtot,irtot),
     &     resp_dvp2(irtot,irtot)

c---- arrays for projecting resr
      real resd(irtot),
     &     resd_vars(irtot,ivtot),
     &     resd_dvp(irtot,irtot),
     &     resd_dum(irtot,3)

      real resd_bf(irtot,3,3),
     &     resd_bm(irtot,3,3)
      real resd_dj(3)

c---- Newton variable changes
      real dvp(irtot)

c---- Newton convergence info
      real dvmax(ivtot)
      integer kvmax(ivtot)

      real dr1old(3), dr2old(3), dxdold(3)
      real dr1new(3), dr2new(3), dxdnew(3)

      logical ldfix
      logical lconv
      logical lmset

      character*1 cfix

c------------------------------------------------------
cc---- global-variable solution stuff (tentative)
c      parameter(ngdim=500)
c      real gres_vars(ivtot,kdim)
c      real gres_varg(ngdim)
c      real rgl(ngdim), agl(ngdim,ngdim)
c      real fsum(ngdim), asum(ngdim)
c      real belem(3,ngdim)
c      real aelem(3,ngdim)
c      real felem(3,ngdim)
c      real vel(3), vrel(3)
c      real ftot(3)

c------------------------------------------------------
c      include 'gaussq.inc'
c      common / rcom /
c     &  rg(3,ngaussq),
c     &  gang(2,ngaussq),
c     &  fang(2,ngaussq),
c     &  gncg(3,ngaussq),
c     &  fncg(3,ngaussq)

c----------------------------------------------
c---- variable labels
      character*2 vlabel(8)
      data vlabel 
     & / 'rx', 
     &   'ry', 
     &   'rz', 
     &   'dx',
     &   'dy',
     &   'dz',
     &   'ps',
     &   'Ac'  /

c---- minimum allowable element area ratio Anew/Aold over a Newton update
      armin = 0.8

c---- number of equations to be solved
      if(ldrill) then 
       neq = 6
      else
       neq = 5
      endif
      
      if(.not. lvinit) then
c----- initialize solution to undeformed values
       do k = 1, nnode
         vars(ivrx,k) = pars(lvr0x,k)
         vars(ivry,k) = pars(lvr0y,k)
         vars(ivrz,k) = pars(lvr0z,k)

         vars(ivdx,k) = pars(lvn0x,k)
         vars(ivdy,k) = pars(lvn0y,k)
         vars(ivdz,k) = pars(lvn0z,k)

         vars(ivps,k) = 0.
       enddo

      endif

c---- zero out "old" residuals
      do k = 1, nnode
        do ir = irm1, irm2
          resp0(ir,k) = 0.
        enddo
      enddo


c========================================================================
c---- Newton iteration loop  (start at least one pass)
      ldfix = .false.
      do 900 iter = 1, max(itmax,1)

c---- pickup here for another Newton iteration with same iter counter
   5  continue

c---- clear residuals, Jacobians, and Jacobian node indices
      do k = 1, nnode
c------ cartesian Residuals and Jacobians
        do keq = 1, irtot
          resc(keq,k) = 0.0
          do id = 1, nddim
            do iv = 1, ivtot
              resc_vars(keq,iv,id,k) = 0.
            enddo
          enddo
        enddo

        ndvp(k) = 0
        do id = 1, nddim
          kdvp(id,k) = 0
        enddo

        ares(k) = 0.
      enddo

c---- clear geometric BC flags
      do ibc = 1, nbcnode
        ibr1(ibc) = 0
        ibr2(ibc) = 0
        ibr3(ibc) = 0
        ibd1(ibc) = 0
        ibd2(ibc) = 0
        ibd3(ibc) = 0
      enddo

c=============================================================
c---- set default projection vectors b1,b2,bn
      do k = 1, nnode
c------ set bn = d
        ibn = 3
        bf(1,ibn,k) = vars(ivdx,k)
        bf(2,ibn,k) = vars(ivdy,k)
        bf(3,ibn,k) = vars(ivdz,k)

        bm(1,ibn,k) = vars(ivdx,k)
        bm(2,ibn,k) = vars(ivdy,k)
        bm(3,ibn,k) = vars(ivdz,k)
        do kk = 1, 3
          bf_dj( 1,kk,ibn,k) = 0.
          bf_dj( 2,kk,ibn,k) = 0.
          bf_dj( 3,kk,ibn,k) = 0.
          bf_dj(kk,kk,ibn,k) = 1.0

          bm_dj( 1,kk,ibn,k) = 0.
          bm_dj( 2,kk,ibn,k) = 0.
          bm_dj( 3,kk,ibn,k) = 0.
          bm_dj(kk,kk,ibn,k) = 1.0
        enddo

c------ set b1,b2 orthogonal to d
        ib1 = 1
        ib2 = 2
        call hsmbb(vars(ivdx,k),
     &             bf(1,ib1,k), bf_dj(1,1,ib1,k),
     &             bf(1,ib2,k), bf_dj(1,1,ib2,k) )
        call hsmbb(vars(ivdx,k),
     &             bm(1,ib1,k), bm_dj(1,1,ib1,k),
     &             bm(1,ib2,k), bm_dj(1,1,ib2,k) )
      enddo

c---------------------------------------------------------------
c---- reset b1,b2,bn vectors for boundary nodes, to enable geometry BCs
      do ibc = 1, nbcnode
        k = kbcnode(ibc)

c------ find whether nBC or tBC vectors are needed at this boundary node
        lbc = lbcnode(ibc)
        lbcm = ( lbc           /1000) * 1000
        lbcf = ((lbc-lbcm     )/100 ) * 100
        lbcd = ((lbc-lbcm-lbcf)/10  ) * 10
        lbcr =   lbc-lbcm-lbcf-lbcd

c------------------------------------------
c------ check to make sure required BC vectors have been set and are suitable
        if(lbcr .eq. 1 .or. lbcr .eq. 2) then
          an1sq = parp(lpn1x,ibc)**2
     &          + parp(lpn1y,ibc)**2
     &          + parp(lpn1z,ibc)**2
          if(an1sq .le. 0.0) then
           write(*,*) 
     &      'HSMSOL: Zero nBC1 vector for position-BC node ibc =', ibc
           stop
          endif
        endif

        if(lbcr .eq. 2) then
          an2sq = parp(lpn2x,ibc)**2
     &          + parp(lpn2y,ibc)**2
     &          + parp(lpn2z,ibc)**2
          if(an2sq .eq. 0.0) then
           write(*,*) 
     &      'HSMSOL: Zero nBC2 vector for position-BC node ibc =', ibc
           stop
          endif
        endif

c- - - - - - - - - - - - - - - - - - - - - - - -
        if(lbcd .eq. 20 .or. lbcd .eq. 30) then
          at1sq = parp(lpt1x,ibc)**2
     &          + parp(lpt1y,ibc)**2
     &          + parp(lpt1z,ibc)**2
          if(at1sq .le. 0.0) then
           write(*,*) 
     &      'HSMSOL: Zero tBC1 vector for angle-BC node ibc =', ibc
           stop
          endif

          dt1 = ( vars(ivdx,k)*parp(lpt1x,ibc)
     &          + vars(ivdy,k)*parp(lpt1y,ibc)
     &          + vars(ivdz,k)*parp(lpt1z,ibc) ) / sqrt(at1sq)
          if(abs(dt1) .gt. 0.9) then
           write(*,*)
     &      'HSMSOL Warning: tBC1 vector has large component',
     &        ' along d for angle-BC node ibc =', ibc
          endif
        endif

        if(lbcd .eq. 30) then
          at2sq = parp(lpt2x,ibc)**2
     &          + parp(lpt2y,ibc)**2
     &          + parp(lpt2z,ibc)**2
          if(at2sq .le. 0.0) then
           write(*,*) 
     &      'HSMSOL: Zero tBC2 vector for angle-BC node ibc =', ibc
           stop
          endif

          dt2 = ( vars(ivdx,k)*parp(lpt2x,ibc)
     &          + vars(ivdy,k)*parp(lpt2y,ibc)
     &          + vars(ivdz,k)*parp(lpt2z,ibc) ) / sqrt(at2sq)
          if(abs(dt2) .gt. 0.9) then
           write(*,*)
     &      'HSMSOL Warning: tBC2 vector has large component',
     &        ' along d for strong-BC node ibc =', ibc
          endif
        endif

c------------------------------------------------
        if(lbcr .eq. 1) then
c-------- reset bf1,bf2,bf3 which is closest to nBC1 to be equal to nBC1,
c-         orthogonalize the rest
          call bort3(parp(lpn1x,ibc),
     &               bf(1,1,k),bf(1,2,k),bf(1,3,k),ibr1(ibc))

c- - - - - - - - - - - - - - - - - - - - - - - -
        elseif(lbcr .eq. 2) then
c-------- reset bf1,bf2,bf3 which is closest to nBC1 to be equal to nBC1,
          call bset3(parp(lpn1x,ibc),
     &               bf(1,1,k),bf(1,2,k),bf(1,3,k),ibr1(ibc))

c-------- reset onr of remaining bf`s closest to nBC2 to be equal to nBC2,
c-         and then set remaining 3rd bf to be normal to the first two.
          if    (ibr1(ibc) .eq. 1) then
c--------- bf1 was reset above, so reset remaining bf2 or bfn
           call bset2(parp(lpn2x,ibc),
     &                bf(1,2,k),bf(1,3,k),ibtmp)
           if(ibtmp.eq.1) then
c---------- bf2 was reset, so reset bfn = bf1 x bf2
            ibr2(ibc) = 2
            call cross(bf(1,1,k),bf(1,2,k),bf(1,3,k))
           else
c---------- bfn was reset, so reset bf2 = bfn x bf1
            ibr2(ibc) = 3
            call cross(bf(1,3,k),bf(1,1,k),bf(1,2,k))
           endif
              
          elseif(ibr1(ibc) .eq. 2) then
c--------- bf2 was reset above, so reset remaining bf1 or bfn
           call bort2(parp(lpn2x,ibc),
     &                bf(1,1,k),bf(1,3,k),ibtmp)
           if(ibtmp.eq.1) then
c---------- bf1 was reset, so reset bfn = bf1 x bf2
            ibr2(ibc) = 1
            call cross(bf(1,1,k),bf(1,2,k),bf(1,3,k))
           else
c---------- bfn was reset, so reset bf1 = bf2 x bfn
            ibr2(ibc) = 3
            call cross(bf(1,2,k),bf(1,3,k),bf(1,1,k))
           endif

          else ! ibr1(ibc) .eq. 3
c--------- bfn was reset above, so reset remaining bf1 or bf2
           call bort2(parp(lpn2x,ibc),
     &                bf(1,1,k),bf(1,2,k),ibtmp)
           if(ibtmp.eq.1) then
c---------- bf1 was reset, so reset bf2 = bfn x bf1
            ibr2(ibc) = 1
            call cross(bf(1,3,k),bf(1,1,k),bf(1,2,k))
           else
c---------- bf2 was reset, so reset bf1 = bf2 x bf3
            ibr2(ibc) = 2
            call cross(bf(1,2,k),bf(1,3,k),bf(1,1,k))
           endif
          endif

c- - - - - - - - - - - - - - - - - - - - - - - -
        elseif(lbcr .eq. 3) then
c------- all position components fixed... bf vectors will be used unmodified
         ibr1(ibc) = 1
         ibr2(ibc) = 2
         ibr3(ibc) = 3

        endif
c- - - - - - - - - - - - - - - - - - - - - - - -

        if(lbcr.eq.1 .or. lbcr.eq.2) then
c------- bf vectors no longer depend on d, so zero out their Jacobians
         do kk = 1, 3
           do ib = 1, 3
             bf_dj(1,kk,ib,k) = 0.
             bf_dj(2,kk,ib,k) = 0.
             bf_dj(3,kk,ib,k) = 0.
           enddo
         enddo
        endif

c------------------------------------------------
        if(lbcd .eq. 10) then
          ibd3(ibc) = 3

c------------------------------------------------
        elseif(lbcd .eq. 20) then
c-------- reset bm1,bm2 which is closest to t1 to be equal to t1
          call bset2(parp(lpt1x,ibc),
     &               bm(1,1,k),bm(1,2,k),ibd1(ibc))
          do ib = 1, 2
            do kk = 1, 3
              bm_dj(1,kk,ib,k) = 0.
              bm_dj(2,kk,ib,k) = 0.
              bm_dj(3,kk,ib,k) = 0.
            enddo
          enddo

c          if(k .eq. 10) then
c            write(*,*)
c            write(*,*) 'ib', ibd1(ibc)
c            write(*,*) 'dxt', dxt
c            write(*,*) 'b1 ', bm(1,1,k),bm(2,1,k),bm(3,1,k)
c            write(*,*) 'b2 ', bm(1,2,k),bm(2,2,k),bm(3,2,k)
c            write(*,*) 'b3 ', bm(1,3,k),bm(2,3,k),bm(3,3,k)
c            write(*,*) 'd  ', vars(ivdx,k),vars(ivdy,k),vars(ivdz,k)
c          endif

c-------- set other bm orthogonal to d and already-set bm1 or bm2
          ibset = ibd1(ibc)
          if(ibset.eq.1) then
c--------- bm1 was already set, so now set bm2
           ib = 2
          else
c--------- bm2 was already set, so now set bm1
           ib = 1
          endif

          bm(1,ib,k) = vars(ivdy,k)*bm(3,ibset,k)
     &               - vars(ivdz,k)*bm(2,ibset,k)
          bm(2,ib,k) = vars(ivdz,k)*bm(1,ibset,k)
     &               - vars(ivdx,k)*bm(3,ibset,k)
          bm(3,ib,k) = vars(ivdx,k)*bm(2,ibset,k)
     &               - vars(ivdy,k)*bm(1,ibset,k)

          bm_dj(1,2,ib,k) =  bm(3,ibset,k)
          bm_dj(1,3,ib,k) = -bm(2,ibset,k)
          bm_dj(2,3,ib,k) =  bm(1,ibset,k)
          bm_dj(2,1,ib,k) = -bm(3,ibset,k)
          bm_dj(3,1,ib,k) =  bm(2,ibset,k)
          bm_dj(3,2,ib,k) = -bm(1,ibset,k)

c          if(k .eq. 10) then
c            write(*,*)
c            write(*,*) 'ib ibset', ib, ibset
c            write(*,*) 'b1 ', bm(1,1,k),bm(2,1,k),bm(3,1,k)
c            write(*,*) 'b2 ', bm(1,2,k),bm(2,2,k),bm(3,2,k)
c            write(*,*)
c          endif

          ibd3(ibc) = 3

        elseif(lbcd .eq. 30) then
c-------- reset bm1,bm2 parallel to tBC1,tBC2
          ibd1(ibc) = 1
          t1a = sqrt( parp(lpt1x,ibc)**2
     &              + parp(lpt1y,ibc)**2
     &              + parp(lpt1z,ibc)**2 )
          bm(1,1,k) = parp(lpt1x,ibc) / t1a
          bm(2,1,k) = parp(lpt1y,ibc) / t1a
          bm(3,1,k) = parp(lpt1z,ibc) / t1a

          ibd2(ibc) = 2
          t2a = sqrt( parp(lpt2x,ibc)**2
     &              + parp(lpt2y,ibc)**2
     &              + parp(lpt2z,ibc)**2 )
          bm(1,2,k) = parp(lpt2x,ibc) / t2a
          bm(2,2,k) = parp(lpt2y,ibc) / t2a
          bm(3,2,k) = parp(lpt2z,ibc) / t2a

c-------- reset bm3 = bm1 x bm2
          ibd3(ibc) = 3
          call cross(bm(1,1,k),bm(1,2,k),bm(1,3,k))

c-------- check to make sure bm3 has significant magnitude 
          bm3sq = bm(1,3,k)**2 + bm(2,3,k)**2 + bm(3,3,k)**2
          if(bm3sq .lt. 0.01) then
           write(*,*)
     &    'HSMSOL: tBC1 and tBC2 are close to parallel, at node k=', k
          endif

          do ib = 1, 2
            do kk = 1, 3
              bm_dj(1,kk,ib,k) = 0.
              bm_dj(2,kk,ib,k) = 0.
              bm_dj(3,kk,ib,k) = 0.
            enddo
          enddo

        endif

      enddo ! next ibc


c      do k = 1, nnode
c        do ib = 1, 3
c          do kk = 1, 3
c            bf(kk,ib,k) = 0.
c            bm(kk,ib,k) = 0.
c            bf_dj(kk,ib,1,k) = 0.
c            bf_dj(kk,ib,2,k) = 0.
c            bf_dj(kk,ib,3,k) = 0.
c            bm_dj(kk,ib,1,k) = 0.
c            bm_dj(kk,ib,2,k) = 0.
c            bm_dj(kk,ib,3,k) = 0.
c          enddo
c          bf(ib,ib,k) = 1.0
c          bm(ib,ib,k) = 1.0
c        enddo
c      enddo


c      if(iter.eq.1) then
      if(iter.eq.-1) then
      efac = 0.04*rref
      write(*,*) 'Writing to fort.21:   bf1   vectors'
      write(*,*) 'Writing to fort.22:   bf2   vectors'
      write(*,*) 'Writing to fort.23:   bfn   vectors'
      write(*,*) 'Writing to fort.24:   bm1   vectors'
      write(*,*) 'Writing to fort.25:   bm2   vectors'
      write(*,*) 'Writing to fort.26:   bmn   vectors'
      call hsmoutv(21,efac, nnode,vars, 9,1,bf)
      call hsmoutv(22,efac, nnode,vars, 9,4,bf)
      call hsmoutv(23,efac, nnode,vars, 9,7,bf)
      call hsmoutv(24,efac, nnode,vars, 9,1,bm)
      call hsmoutv(25,efac, nnode,vars, 9,4,bm)
      call hsmoutv(26,efac, nnode,vars, 9,7,bm)
      stop
      endif

c=============================================================
c---- go over all elements, accumulating residuals 
c-     and Jacobians onto each node in each element.
      do 100 n = 1, nelem
c------ nodes of this element
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

c------ calculate this element's contribution to equation residuals and Jacobian
        if(k4.eq.0) then
c------- tri element
         nres = 3
         k4 = k3  ! avoids any array OOB problems with k4=0
        else
c------- quad element
         nres = 4
        endif

        call hsmeqn(nres, lrcurv, ldrill,
     &        parg,
     &        pars(1,k1), pars(1,k2), pars(1,k3), pars(1,k4), 
     &        vars(1,k1), vars(1,k2), vars(1,k3), vars(1,k4), 
     &   resi, 
     &   resi_vars(1,1,1,1),
     &   resi_vars(1,1,1,2),
     &   resi_vars(1,1,1,3),
     &   resi_vars(1,1,1,4),
     &   aresi, dresi )


c        if(iter.eq.8) then
c          lu = 31
c          do ig = 1, ngaussq
c          write(lu,'(1x,19g13.5)')
c     &      (rg(k,ig), k=1,3),
c     &      (gncg(k,ig), k=1,3),
c     &      (fncg(k,ig), k=1,3),
c     &      (gang(ia,ig), ia=1,2),
c     &      (fang(ia,ig), ia=1,2)
c          enddo
c          write(lu,*)
c          write(lu,*)
c        endif

c------ accumulate residuals and Jacobians for this element,
c-      for each residual weighting function  Wi 
        do 10 ires = 1, nres
          kres = kelem(ires,n)

c-------- accumulate residuals and Jacobians for nodes of this element
          do keq = 1, irtot
            resc(keq,kres) = resc(keq,kres) + resi(keq,ires)
          enddo

          call jadd(1,irtot, ivtot, irtot, nddim,
     &           resi_vars(1,1,ires,1), k1,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))

          call jadd(1,irtot, ivtot, irtot, nddim,
     &           resi_vars(1,1,ires,2), k2,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))

          call jadd(1,irtot, ivtot, irtot, nddim,
     &           resi_vars(1,1,ires,3), k3,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))

          if(nres.eq.4) then
           call jadd(1,irtot, ivtot, irtot, nddim,
     &           resi_vars(1,1,ires,4), k4,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))
          endif

c-------- accumulate residual areas
          ares(kres) = ares(kres) + aresi(ires)

 10     continue ! next ires residual for this element
 100  continue ! next element

c     do k = 1, nnode
c       if(vars(ivrx,k) .gt. 0.29 .and. 
c    &     vars(ivrx,k) .lt. 0.31        ) then
c        keq = 1
c        write(80,*) vars(ivry,k), resc(keq,k)
c       endif
c     enddo

c      if(iter.eq.itmax) then
c       do k = 1, nnode
c         write(*,44) k, (resc(in,k), in=1, irtot)
c 44      format(1x,i4, 6e12.4)
c       enddo
c      endif


c==================================================================
c---- go over nodes to set point loads
      do 200 ibc = 1, nbcnode
        k = kbcnode(ibc)

c------ find whether a force or moment is applied at this node
        lbc = lbcnode(ibc)
        lbcm = ( lbc           /1000) * 1000
        lbcf = ((lbc-lbcm     )/100 ) * 100
ccc     lbcd = ((lbc-lbcm-lbcf)/10  ) * 10
ccc     lbcr =   lbc-lbcm-lbcf-lbcd

        if(lbcf.eq.100 .or. lbcm.eq.1000) then
c-------- add force and/or moment to node residuals
          call hsmpbc(parp(1,ibc),
     &                pars(1,k) ,
     &                vars(1,k) ,
     &                resc(1,k),
     &                resc_vars(1,1,1,k))
        endif
 200  continue   ! next node with load

c=============================================================
c---- go over edges to set Neumann (loading) BCs
      do 300 ibc = 1, nbcedge
c------ endpoint nodes of this edge
        k1 = kbcedge(1,ibc)
        k2 = kbcedge(2,ibc)

        call hsmfbc(lrcurv, ldrill,
     &              pare(1,ibc),
     &              pars(1,k1) , pars(1,k2) ,
     &              vars(1,k1) , vars(1,k2) ,
     &              resf, 
     &              resf_vars(1,1,1,1),
     &              resf_vars(1,1,1,2),
     &              dresi )

c------ accumulate residuals and Jacobians for nodes of this edge
        do 30 ires = 1, 2
          kres = kbcedge(ires,ibc)

          do keq = 1, irtot
            resc(keq,kres) = resc(keq,kres) + resf(keq,ires)
          enddo

          call jadd(1,irtot, ivtot, irtot, nddim,
     &           resf_vars(1,1,ires,1), k1,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))

          call jadd(1,irtot, ivtot, irtot, nddim,
     &           resf_vars(1,1,ires,2), k2,
     &           resc_vars(1,1,1,kres), kdvp(1,kres), ndvp(kres))

 30     continue ! next ires residual in this edge
 300  continue   ! next boundary edge with loading

c------------------------------------------------------------
c---- add force and moment residuals for joint node pairs
      do 400 j = 1, njoint
c------ node pair at this joint point
        k1 = kjoint(1,j)
        k2 = kjoint(2,j)

c------ sum node 1,2 force and moment equations, put into node 2 equations slots
        keq1 = 1
        keq2 = 6
        do keq = keq1, keq2
          resc(keq,k2) = resc(keq,k2) + resc(keq,k1)
        enddo

c------ also sum Jacobians, making new matrix block elements as needed
        do id = 1, ndvp(k1)
          kd1 = kdvp(id,k1)
          call jadd(keq1,keq2, ivtot, irtot, nddim,
     &              resc_vars(1,1,id,k1), kd1,
     &              resc_vars(1,1, 1,k2), kdvp(1,k2), ndvp(k2))
        enddo
 400  continue ! next joint node


c      do k = 1, nnode
c        write(*,'(1x,3f10.6,6e13.5)')
c     &      (vars(iv,k),iv=1,3), (resc(iv,k), iv=1, ivtot)
c        if(mod(k,6) .eq. 0) write(*,*)
c      enddo

cc----- display stencil node indices
c      do kres = 1, nnode
c        write(*,'(1x,i5,a,18i5)') 
c     &        kres, ':', (kdvp(id,kres), id=1,ndvp(kres))
c      enddo
c      pause

c=============================================================
c---- go over all residuals for projection
      do 500 kres = 1, nnode
c------ project cartesian-basis residuals resc onto b vectors,
c-        giving b-vector basis residuals resp
        call hsmfpr(bf(1,1,kres),
     &              bf(1,2,kres),
     &              bf(1,3,kres),
     &              bm(1,1,kres),
     &              bm(1,2,kres),
     &              bm(1,3,kres),
     &              resc(1,kres),
     &              resp(1,kres),
     &              resp_bf(1,1,1),
     &              resp_bf(1,1,2),
     &              resp_bf(1,1,3),
     &              resp_bm(1,1,1),
     &              resp_bm(1,1,2),
     &              resp_bm(1,1,3) )

c------ project the Jacobians the same way, put into local work array  resp_vars
        do id = 1, ndvp(kres)
          do iv = 1, ivtot
            call hsmfpr(bf(1,1,kres),
     &                  bf(1,2,kres),
     &                  bf(1,3,kres),
     &                  bm(1,1,kres),
     &                  bm(1,2,kres),
     &                  bm(1,3,kres),
     &                  resc_vars(1,iv,id,kres),
     &                  resp_vars(1,iv,id),
     &                  resp_dum,
     &                  resp_dum,
     &                  resp_dum,
     &                  resp_dum,
     &                  resp_dum,
     &                  resp_dum )
          enddo
        enddo

c------ locate id index of kres (diagonal) Jacobian element
        idres = 0
        do id = 1, ndvp(kres)
          if(kdvp(id,kres) .eq. kres) then
           idres = id
          endif
        enddo

        if(idres .eq. 0) then
         write(*,*) '? HSMSOL: No (k,k) Jacobian element for k =', kres

        else
c------- add contribution from dependence of b vectors on primary variable d
         do keq = 1, neq
           do kk = 1, 3
             resp_dj(kk) = resp_bf(keq,1,1)*bf_dj(1,kk,1,kres)
     &                   + resp_bf(keq,2,1)*bf_dj(2,kk,1,kres)
     &                   + resp_bf(keq,3,1)*bf_dj(3,kk,1,kres)
     &                   + resp_bf(keq,1,2)*bf_dj(1,kk,2,kres)
     &                   + resp_bf(keq,2,2)*bf_dj(2,kk,2,kres) 
     &                   + resp_bf(keq,3,2)*bf_dj(3,kk,2,kres) 
     &                   + resp_bf(keq,1,3)*bf_dj(1,kk,3,kres)
     &                   + resp_bf(keq,2,3)*bf_dj(2,kk,3,kres)
     &                   + resp_bf(keq,3,3)*bf_dj(3,kk,3,kres)
     &                   + resp_bm(keq,1,1)*bm_dj(1,kk,1,kres)
     &                   + resp_bm(keq,2,1)*bm_dj(2,kk,1,kres)
     &                   + resp_bm(keq,3,1)*bm_dj(3,kk,1,kres)
     &                   + resp_bm(keq,1,2)*bm_dj(1,kk,2,kres)
     &                   + resp_bm(keq,2,2)*bm_dj(2,kk,2,kres) 
     &                   + resp_bm(keq,3,2)*bm_dj(3,kk,2,kres) 
           enddo
           resp_vars(keq,ivdx,idres) = 
     &     resp_vars(keq,ivdx,idres) + resp_dj(1)
           resp_vars(keq,ivdy,idres) = 
     &     resp_vars(keq,ivdy,idres) + resp_dj(2)
           resp_vars(keq,ivdz,idres) = 
     &     resp_vars(keq,ivdz,idres) + resp_dj(3)
         enddo

        endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ project xyz variables at each Jacobian node onto b1,b2,bn of that node,
c-       install final result in projected-variable Jacobian resp_dvp array
        do id = 1, ndvp(kres)
          k = kdvp(id,kres)
          call hsmvpr(bf(1,1,k),
     &                bf(1,2,k),
     &                bf(1,3,k),
     &                bm(1,1,k),
     &                bm(1,2,k),
     &                resp_vars(1,1,id),
     &                resp_dvp(1,1,id,kres) )
        enddo

 500  continue  ! next residual node

c==================================================================
c---- set joint position and angle continuity
      do 600 ij = 1, njoint
        k1 = kjoint(1,ij)
        k2 = kjoint(2,ij)

c------ put kinematic equations into node k1
c-     (k1 and k2 force and moment residuals were already summed into node k2)
        kres = k1

c------ set up joint Dirichlet residuals and Jacobians
        call hsmjbc(pars(1,k1) , pars(1,k2) ,
     &              vars(1,k1) , vars(1,k2) ,
     &              resd, resd_var1, resd_var2 )

c        write(*,*)
c        write(*,*) 'resd_var'
c        do k = 1, 3
c        write(*,5522) (resd_var1(k,l), l=1, 3),
c     &                (resd_var2(k,l), l=1, 3)
c        enddo
 5522   format(1x,3f8.3,3x,3f8.3)

c------ set joint edge, tangent, and bisector vectors
        call hsmjtl(pars(lvn0x,k1), pars(lvn0x,k2),
     &              t0j, t0j_n01, t0j_n02,
     &              l0j, l0j_n01, l0j_n02,
     &              n0j, n0j_n01, n0j_n02 )

        call hsmjtl(vars(ivdx,k1), vars(ivdx,k2),
     &              tj, tj_d1, tj_d2,
     &              lj, lj_d1, lj_d2,
     &              nj, nj_d1, nj_d2 )

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ project xyz residuals onto tj,lj,nj , 
c-      put into node kres equation slots
        call hsmfpr(tj,
     &              lj,
     &              nj,
     &              tj,
     &              lj,
     &              nj,
     &              resd,
     &              resp(1,kres),
     &              resp_bf(1,1,1),
     &              resp_bf(1,1,2),
     &              resp_bf(1,1,3),
     &              resp_bm(1,1,1),
     &              resp_bm(1,1,2),
     &              resp_bm(1,1,3) )
c----- note:  
c-      resp_bf(..1) is really resp_tj
c-      resp_bf(..2) is really resp_lj
c-      resp_bf(..3) is really resp_nj
c-      resp_bm(..1) is really resp_tj
c-      resp_bm(..2) is really resp_lj
c-      resp_bm(..3) is really resp_nj

c------ project the Jacobians the same way
        do iv = 1, ivtot
          call hsmfpr(tj,
     &                lj,
     &                nj,
     &                tj,
     &                lj,
     &                nj,
     &                resd_var1(1,iv),
     &                resp_var1(1,iv),
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum )
          call hsmfpr(tj,
     &                lj,
     &                nj,
     &                tj,
     &                lj,
     &                nj,
     &                resd_var2(1,iv),
     &                resp_var2(1,iv),
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum,
     &                resp_dum )
        enddo

c        write(*,*)
c        write(*,*) 'resp_var'
c        do k = 1, 3
c        write(*,5522) (resp_var1(k,l), l=1, 3),
c     &                (resp_var2(k,l), l=1, 3)
c        enddo

c------ add contributions from dependence of tj,lj,nj vectors on d1,d2
        do keq = 1, neq
          do kk = 1, 3
            iv = ivdx-1 + kk
            resp_var1(keq,iv) = 
     &      resp_var1(keq,iv) + resp_bf(keq,1,1)*tj_d1(1,kk)
     &                        + resp_bf(keq,2,1)*tj_d1(2,kk)
     &                        + resp_bf(keq,3,1)*tj_d1(3,kk)
     &                        + resp_bf(keq,1,2)*lj_d1(1,kk)
     &                        + resp_bf(keq,2,2)*lj_d1(2,kk) 
     &                        + resp_bf(keq,3,2)*lj_d1(3,kk) 
     &                        + resp_bf(keq,1,3)*nj_d1(1,kk)
     &                        + resp_bf(keq,2,3)*nj_d1(2,kk)
     &                        + resp_bf(keq,3,3)*nj_d1(3,kk)
     &                        + resp_bm(keq,1,1)*tj_d1(1,kk)
     &                        + resp_bm(keq,2,1)*tj_d1(2,kk)
     &                        + resp_bm(keq,3,1)*tj_d1(3,kk)
     &                        + resp_bm(keq,1,2)*lj_d1(1,kk)
     &                        + resp_bm(keq,2,2)*lj_d1(2,kk) 
     &                        + resp_bm(keq,3,2)*lj_d1(3,kk) 
     &                        + resp_bm(keq,1,3)*lj_d1(1,kk)
     &                        + resp_bm(keq,2,3)*lj_d1(2,kk) 
     &                        + resp_bm(keq,3,3)*lj_d1(3,kk) 
            resp_var2(keq,iv) = 
     &      resp_var2(keq,iv) + resp_bf(keq,1,1)*tj_d2(1,kk)
     &                        + resp_bf(keq,2,1)*tj_d2(2,kk)
     &                        + resp_bf(keq,3,1)*tj_d2(3,kk)
     &                        + resp_bf(keq,1,2)*lj_d2(1,kk)
     &                        + resp_bf(keq,2,2)*lj_d2(2,kk) 
     &                        + resp_bf(keq,3,2)*lj_d2(3,kk) 
     &                        + resp_bf(keq,1,3)*nj_d2(1,kk)
     &                        + resp_bf(keq,2,3)*nj_d2(2,kk)
     &                        + resp_bf(keq,3,3)*nj_d2(3,kk)
     &                        + resp_bm(keq,1,1)*tj_d2(1,kk)
     &                        + resp_bm(keq,2,1)*tj_d2(2,kk)
     &                        + resp_bm(keq,3,1)*tj_d2(3,kk)
     &                        + resp_bm(keq,1,2)*lj_d2(1,kk)
     &                        + resp_bm(keq,2,2)*lj_d2(2,kk) 
     &                        + resp_bm(keq,3,2)*lj_d2(3,kk) 
     &                        + resp_bm(keq,1,3)*lj_d2(1,kk)
     &                        + resp_bm(keq,2,3)*lj_d2(2,kk) 
     &                        + resp_bm(keq,3,3)*lj_d2(3,kk) 
          enddo
        enddo

c------ overwrite lj constraint with nonlinear constraint on angle between d1,d2
        call hsmdda(pars(lvn0x,k1),pars(lvn0x,k2),
     &              l0j, l0j_n01, l0j_n02,
     &              da0, da0_n01, da0_n02)
        call hsmdda(vars(ivdx ,k1),vars(ivdx ,k2),
     &              lj, lj_d1, lj_d2,
     &              da, da_d1, da_d2)
        keq = 4
        resp(keq,kres) = da - da0
        do k = 1, 3
          iv = ivrx-1 + k
          resp_var1(keq,iv) = 0.
          resp_var2(keq,iv) = 0.

          iv = ivdx-1 + k
          resp_var1(keq,iv) = da_d1(k)
          resp_var2(keq,iv) = da_d2(k)
        enddo
        resp_var1(keq,ivps) = 0.
        resp_var2(keq,ivps) = 0.

c        write(*,*)
c        write(*,*) 'resp_var amod'
c        do k = 1, 3
c        write(*,5522) (resp_var1(k,l), l=1, 3),
c     &                (resp_var2(k,l), l=1, 3)
c        enddo

c------ project variables, thus modifying Jacobians
        call hsmvpr(bf(1,1,k1),
     &              bf(1,2,k1),
     &              bf(1,3,k1),
     &              bm(1,1,k1),
     &              bm(1,2,k1),
     &              resp_var1,
     &              resp_dvp1 )

        call hsmvpr(bf(1,1,k2),
     &              bf(1,2,k2),
     &              bf(1,3,k2),
     &              bm(1,1,k2),
     &              bm(1,2,k2),
     &              resp_var2,
     &              resp_dvp2 )

cc     real resp_dvp(irtot,irtot,nddim,kdim)
cc     real resp_dvp1(irtot,irtot),

       if(itmax.ne.-2) then
c------ eliminate existing Jacobian blocks of residual node kres,
c-       since new ones will be set for geometric-continuity equation
c-       (skip this for matrix-size call, to get right size for HSMDEP)
        ndvp(kres) = 0
       endif

c----- install new residuals and Jacobians
       call jadd(1,neq,irtot,irtot,nddim,
     &           resp_dvp1, k1,
     &           resp_dvp(1,1,1,kres), kdvp(1,kres),ndvp(kres))
       call jadd(1,neq,irtot,irtot,nddim,
     &           resp_dvp2, k2,
     &           resp_dvp(1,1,1,kres), kdvp(1,kres),ndvp(kres))

 600  continue ! next joint node pair

c==================================================================
c---- go over nodes to set Dirichlet BCs
      do 700 ibc = 1, nbcnode
        kres = kbcnode(ibc)

c------ set all usable Dirichlet residuals and Jacobians
        call hsmrbc(parp(1,ibc),
     &              pars(1,kres) ,
     &              vars(1,kres) ,
     &              resr, 
     &              resr_vars )

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ project xyz residuals at kres node onto bf,bm of kres node
        call hsmdpr(bf(1,1,kres),
     &              bf(1,2,kres),
     &              bf(1,3,kres),
     &              bm(1,1,kres),
     &              bm(1,2,kres),
     &              resr,
     &              resd,
     &              resd_bf(1,1,1),
     &              resd_bf(1,1,2),
     &              resd_bf(1,1,3),
     &              resd_bm(1,1,1),
     &              resd_bm(1,1,2) )

c------ must project the Jacobians the same way
        do iv = 1, ivtot
          call hsmdpr(bf(1,1,kres),
     &                bf(1,2,kres),
     &                bf(1,3,kres),
     &                bm(1,1,kres),
     &                bm(1,2,kres),
     &                resr_vars(1,iv),
     &                resd_vars(1,iv),
     &                resd_dum,
     &                resd_dum,
     &                resd_dum,
     &                resd_dum,
     &                resd_dum )

c         if(iv.eq.4) then
c         write(*,*) 
c         write(*,*) '  resr_vars', (resr_vars(kv,iv), kv=4, 6)
c         write(*,*) '  bm1      ', (bm(k,1,kres), k=1, 3)
c         write(*,*) '1 resd_vars', (resd_vars(kv,iv), kv=4, 6)
c         pause
c         endif
        enddo

c        iv = 4
c         write(*,*) 
c         write(*,*) '2 resd_var', (resd_vars(kv,iv), kv=4, 6)
c         pause


        if(.false.) then
c------ contribution from dependence of b vectors on primary variable d
        do keq = 1, neq
          do kk = 1, 3
            resd_dj(kk) = resd_bf(keq,1,1)*bf_dj(1,kk,1,kres)
     &                  + resd_bf(keq,2,1)*bf_dj(2,kk,1,kres)
     &                  + resd_bf(keq,3,1)*bf_dj(3,kk,1,kres)
     &                  + resd_bf(keq,1,2)*bf_dj(1,kk,2,kres)
     &                  + resd_bf(keq,2,2)*bf_dj(2,kk,2,kres) 
     &                  + resd_bf(keq,3,2)*bf_dj(3,kk,2,kres) 
     &                  + resd_bf(keq,1,3)*bf_dj(1,kk,3,kres)
     &                  + resd_bf(keq,2,3)*bf_dj(2,kk,3,kres)
     &                  + resd_bf(keq,3,3)*bf_dj(3,kk,3,kres)
     &                  + resd_bm(keq,1,1)*bm_dj(1,kk,1,kres)
     &                  + resd_bm(keq,2,1)*bm_dj(2,kk,1,kres)
     &                  + resd_bm(keq,3,1)*bm_dj(3,kk,1,kres)
     &                  + resd_bm(keq,1,2)*bm_dj(1,kk,2,kres)
     &                  + resd_bm(keq,2,2)*bm_dj(2,kk,2,kres) 
     &                  + resd_bm(keq,3,2)*bm_dj(3,kk,2,kres) 
          enddo
c          write(*,*) 'resd_dj', keq, resd_dj 
          resd_vars(keq,ivdx) = 
     &    resd_vars(keq,ivdx) + resd_dj(1)
          resd_vars(keq,ivdy) = 
     &    resd_vars(keq,ivdy) + resd_dj(2)
          resd_vars(keq,ivdz) = 
     &    resd_vars(keq,ivdz) + resd_dj(3)
        enddo
        endif

c        iv = 4
c         write(*,*) 
c         write(*,*) '3 resd_var', (resd_vars(kv,iv), kv=4, 6)
c         pause

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ project xyz variables at each node onto b1,b2,bn of kres node
        call hsmvpr(bf(1,1,kres),
     &              bf(1,2,kres),
     &              bf(1,3,kres),
     &              bm(1,1,kres),
     &              bm(1,2,kres),
     &              resd_vars,
     &              resd_dvp )

c        iv = 4
c         write(*,*) 
c         write(*,*) '4 resd_var', (resd_vars(kv,iv), kv=4, 6)
c         pause

c--------------------------------------------------------------------
c------ overwrite existing equations for this node with selected 
c-      Dirichlet BC equations, as dictated by lbcnode(ibc) index
        lbc = lbcnode(ibc)

c------ extract  moment,force,angle,position digits
        lbcm = ( lbc           /1000) * 1000
        lbcf = ((lbc-lbcm     )/100 ) * 100
        lbcd = ((lbc-lbcm-lbcf)/10  ) * 10
        lbcr =   lbc-lbcm-lbcf-lbcd

        if(.false.) then
 77     format(1x,i4, 8f8.3)
        write(*,*)
        write(*,*) 'ibc lbcr lbcd', ibc, lbcr, lbcd
        write(*,*) 'ibr', ibr1(ibc), ibr2(ibc), ibr3(ibc)
        write(*,*) 'ibd', ibd1(ibc), ibd2(ibc), ibd3(ibc)
        write(*,*)
        write(*,*) 'resr_vars'
        do keq = 1, 7
          write(*,77) keq, (resr_vars(keq,k), k=1, 7)
        enddo
        write(*,*)
        write(*,*) 'bm1', (bm(k,1,kres),  k=1, 3)
        write(*,*) 'bm2', (bm(k,2,kres),  k=1, 3)
        write(*,*) 'bm3', (bm(k,3,kres),  k=1, 3)
        write(*,*)
        write(*,*) 'resd_var'
        do keq = 1, 6
          write(*,77) keq, (resd_vars(keq,k), k=1, 7)
        enddo
        write(*,*)
        write(*,*) 'resd_dvp'
        do keq = 1, 6
          write(*,77) keq, (resd_dvp(keq,k), k=1, 6)
        enddo
        endif

c------------------------------------------------------------
c------ set position BCs
        if(ibr1(ibc) .eq. 1 .or.
     &     ibr2(ibc) .eq. 1 .or.
     &     ibr3(ibc) .eq. 1     ) then
c-------- replace Rf1 with r1 constraint
          call dbcset(neq, irtot, nddim,
     &                irf1, resd, resd_dvp, kres,
     &                irf1, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

        if(ibr1(ibc) .eq. 2 .or.
     &     ibr2(ibc) .eq. 2 .or.
     &     ibr3(ibc) .eq. 2     ) then
c-------- replace Rf2 with r2 constraint
          call dbcset(neq, irtot, nddim,
     &                irf2, resd, resd_dvp, kres,
     &                irf2, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

        if(ibr1(ibc) .eq. 3 .or.
     &     ibr2(ibc) .eq. 3 .or.
     &     ibr3(ibc) .eq. 3     ) then
c-------- replace Rfn with rn constraint
          call dbcset(neq, irtot, nddim,
     &                irfn, resd, resd_dvp, kres,
     &                irfn, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

c        k = kbcnode(ibc)
c        if(k .eq. 10) then
c          write(*,*) 'ibd123', ibd1(ibc), ibd2(ibc), ibd3(ibc)
c          write(*,*) 'bm1', (bm(kk,1,k), kk=1, 3)
c          write(*,*) 'bm2', (bm(kk,2,k), kk=1, 3)
c          pause
c        endif
c------------------------------------------------------------
c------ set angle BCs
        if(ibd1(ibc) .eq. 1 .or.
     &     ibd2(ibc) .eq. 1 .or.
     &     ibd3(ibc) .eq. 1      ) then
c-------- replace Rm1 with d1 constraint
          call dbcset(neq, irtot, nddim,
     &                irm1, resd, resd_dvp, kres,
     &                irm1, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

        if(ibd1(ibc) .eq. 2 .or.
     &     ibd2(ibc) .eq. 2 .or.
     &     ibd3(ibc) .eq. 2      ) then
c-------- replace Rm2 with d2 constraint
          call dbcset(neq, irtot, nddim,
     &                irm2, resd, resd_dvp, kres,
     &                irm2, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

        if(ibd1(ibc) .eq. 3 .or.
     &     ibd2(ibc) .eq. 3 .or.
     &     ibd3(ibc) .eq. 3     ) then
c-------- replace Rmn with psi constraint
          call dbcset(neq, irtot, nddim,
     &                irmn, resd, resd_dvp, kres,
     &                irmn, resp, resp_dvp, 
     &                kdvp, ndvp)
        endif

 700  continue   ! next boundary node


cc----- display stencil node indices
c      do kres = 1, nnode
c        write(*,'(1x,i5,a,18i5)') 
c     &        kres, ':', (kdvp(id,kres), id=1,ndvp(kres))
c      enddo
c      pause

c      do k = 1, nnode
c        write(*,'(1x,3f10.6,6e13.5)')
c     &      (vars(iv,k),iv=1,3), (resp(in,k), in=1, irtot)
c        if(mod(k,6) .eq. 0) write(*,*)
c      enddo

c---- save residuals for checking convergence in next iteration 
      do k = 1, nnode
        do ir = 1, irtot
          resp0(ir,k) = resp(ir,k)
        enddo
      enddo

      if(ldfix) then
c      if(mod(iter,3) .ne. 1) then
c      if(.false.) then

c------ fix all d and psi variables
        do k = 1, nnode
c-------- find Jacobian slot id which corresponds to this node k
          do id = 1, ndvp(k)
            if(kdvp(id,k) .eq. k) go to 771
          enddo
          write(*,*) 'HSMSOL:  diagonal element not found for node', k
 771      iddiag = id

          do ir = irm1, irmn
c---------- clear residual and Jacobians
            resp(ir,k) = 0.
            do id = 1, ndvp(k)
              do iv = 1, irtot
                 resp_dvp(ir,iv,id,k) = 0.
              enddo
            enddo
c---------- put 1 on diagonal
            resp_dvp(ir,ir,iddiag,k) = 1.0
          enddo
        enddo
      endif


c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c      k = 4
c      do id = 1, ndvp(k)
c        if(kdvp(id,k) .eq. k) go to 772
c      enddo
c      write(*,*) 'HSMSOL:  diagonal element not found for node', k
c 772  iddiag = id
c      do keq = 1, neq
c        write(*,3377) 
c     &   (resp_dvp(keq,iv,iddiag,k), iv=1, neq)
c 3377   format(1x,6f14.2)
c      enddo
c      pause
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


c------------------------------------------------------------
c---- create unstructured Newton matrix
c
c-    if lmset == T,  create pointers and fill matrix array amat
c-    if lmset == F,  create pointers only
      lmset = itmax .ne. -2

      if(ldfix) then
       neqsol = 3
      else
       neqsol = neq
      endif

      call colset(lmset,
     &            kdim,nddim,nmdim,
     &            irtot,neqsol, nnode, resp_dvp,kdvp,ndvp,
     &            namat,ifrst,ilast,mfrst,amat)

c---- total number of blocks in array amat
      nmmat = mfrst(namat+1) - 1

c---- return with required matrix size
      if(itmax .eq. -2) then
       nmdim = nmmat
       return
      endif

      if(iter.eq.1 .and. itmax.le.0 .and. nnode.le.500) then
c     if(iter.eq.1) then
       write(*,*) 'Writing block matrix to fort.70 ...'
       lu = 70
       call bmdump(lu,
     &             irtot,kdim,nmdim,
     &             neq,nnode,nmmat, ifrst,ilast,mfrst,amat)
       stop
      endif

c---- skip solve if itmax = -1
      if(itmax .eq. -1) return

c---------------------------------------------------------------------------
c---- solve Newton system

c----- Newton matrix LU decomposition
c      call sbludi(irtot,neq, namat,ifrst,ilast,mfrst,amat,ipp)

c----- fwd,back substitution of residual vector
c      call sbbksi(irtot,neq, namat,ifrst,ilast,mfrst,amat,ipp,resp)

c      real resp(irtot,kdim)
c      real resp_dvp(irtot,irtot,nddim,kdim)
c      integer kdvp(nddim,kdim)
c      integer ndvp(kdim)

c     if(iter.eq.1) then
      if(.false.) then
       do k = 1, nnode
         idm = 0
         ido = 0
         idp = 0

         km = k-1
         ko = k
         kp = k+1

         do id=1, ndvp(k)
           if(kdvp(id,k) .eq. km) then
            idm = id
           endif
           if(kdvp(id,k) .eq. ko) then
            ido = id
           endif
           if(kdvp(id,k) .eq. kp) then
            idp = id
           endif
          enddo

          write(*,*)
          do id = 1, ndvp(k)
            kj = kdvp(id,k)
            write(*,44) k, kj, 
     &      vars(ivrx,kj), vars(ivry,kj), 0.0,
     &      (resp_dvp(keq,keq,id,k), keq=1, 3),
     &      (resp_dvp(keq,keq,id,k), keq=4, 6),
     &      (bm(kk,1,k), kk=1, 3),
     &      (bm(kk,2,k), kk=1, 3),
     &      (bm(kk,3,k), kk=1, 3)
 44         format(1x, 2i3, 3f10.4,4x, 2(3f12.4,2x), 3(3f7.2,2x))
          enddo
        
c        write(*,55) k, kdvp(1,k), 
c    &         (resp(ir,k), ir=1, irtot), 
c    &         resp_dvp(5,5,idm,k),
c    &         resp_dvp(5,5,ido,k) 
c    &         resp_dvp(5,5,idp,k)
c55      format(1x,2i4, 6e12.4, 4x, 6f10.5)
       enddo
      endif

c---- combined SBLUD,SBBKS operation  (pivots within a block)
      call sbludbi(irtot,neqsol, namat,ifrst,ilast,mfrst,amat,ipp,resp)

c---- combined SBLUD,SBBKS operation  (does NOT pivot within a block)
c      call sbludb(irtot,neq, namat,ifrst,ilast,mfrst,amat,resp)

      if(iter.eq.1 .and. itmax.eq.0 .and. nnode.le.500) then
       write(*,*) 'Writing filled block matrix to fort.71 ...'
       lu = 71
       call bmdump(lu,
     &             irtot,kdim,nmdim,
     &             neqsol,nnode,nmmat, ifrst,ilast,mfrst,amat)
      endif

c---- skip update if itmax = 0
      if(itmax .eq. 0) return

      do k = 1, nnode
c------ projected-variable Newton deltas for this node
        dvp(irmn) = 0.   ! no drilling DOF case
        do ir = 1, neqsol
          dvp(ir) = -resp(ir,k)
          dvp0(ir,k) = -resp(ir,k)
        enddo
        do ir = neqsol+1, neq
          dvp(ir) = 0.
          dvp0(ir,k) = 0.
        enddo

c------ undo projections to get primary-variable Newton deltas
        dvars(ivrx,k) = bf(1,1,k)*dvp(irf1)
     &                + bf(1,2,k)*dvp(irf2)
     &                + bf(1,3,k)*dvp(irfn)
        dvars(ivry,k) = bf(2,1,k)*dvp(irf1)
     &                + bf(2,2,k)*dvp(irf2)
     &                + bf(2,3,k)*dvp(irfn)
        dvars(ivrz,k) = bf(3,1,k)*dvp(irf1)
     &                + bf(3,2,k)*dvp(irf2)
     &                + bf(3,3,k)*dvp(irfn)

        dvars(ivdx,k) = bm(1,1,k)*dvp(irm1)
     &                + bm(1,2,k)*dvp(irm2)
        dvars(ivdy,k) = bm(2,1,k)*dvp(irm1)
     &                + bm(2,2,k)*dvp(irm2)
        dvars(ivdz,k) = bm(3,1,k)*dvp(irm1)
     &                + bm(3,2,k)*dvp(irm2)

        dvars(ivps,k) = dvp(irmn)
      enddo

c---- set max Newton changes
      do iv = 1, ivtot
        dvmax(iv) = 0.
        kvmax(iv) = 0
      enddo

      do k = 1, nnode
        do iv = 1, ivtot
          if(abs(dvars(iv,k)) .ge. dvmax(iv)) then
           dvmax(iv) = abs(dvars(iv,k))
           kvmax(iv) = k
          endif
        enddo
      enddo ! next node k

c---- find max change among all variables, note location and variable
      dmax = 0.0
      kmax = 0
      ivmax = 0
      do iv = 1, ivtot
        if(abs(dvmax(iv)) .ge. abs(dmax)) then
         dmax = dvmax(iv)
         kmax = kvmax(iv)
         ivmax = iv
        endif
      enddo

c---- set Newton change limiter rlx to keep max changes within limits
c-    ivrlx index indicates which variable triggered the limiting, if any
      rlx = 1.0
      ivrlx = 0

      do k = 1, 3
        iv = ivrx-1 + k
        if(rlx*dvmax(iv) .gt.  rlim*rref) then
         rlx =  rlim*rref/dvmax(iv)
         ivrlx = iv
        endif

        iv = ivdx-1 + k
        if(rlx*dvmax(iv) .gt.  alim) then
         rlx =  alim/dvmax(iv)
         ivrlx = iv
        endif
      enddo ! next vector component k

      iv = ivps
      if(rlx*dvmax(iv) .gt.  alim) then
       rlx =  alim/dvmax(iv)
       ivrlx = iv
      endif

      if(ivrlx .ne. 0) then
c----- redefine max-change variable to be one which set limiter
       ivmax = ivrlx
       kmax = kvmax(ivrlx)
      endif

c---- also check for underrelaxaion required by element area ratio changes

c---- rlx-halving loop
      do ipass = 1, 5

      rlxbeg = rlx

      do n = 1, nelem
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

        if(k4 .eq. 0) then
c------- tri element
         dr1old(1) = vars(ivrx,k2) - vars(ivrx,k1)
         dr1old(2) = vars(ivry,k2) - vars(ivry,k1)
         dr1old(3) = vars(ivrz,k2) - vars(ivrz,k1)

         dr2old(1) = vars(ivrx,k3) - vars(ivrx,k2)
         dr2old(2) = vars(ivry,k3) - vars(ivry,k2)
         dr2old(3) = vars(ivrz,k3) - vars(ivrz,k2)

         dr1new(1) = rlx*(dvars(ivrx,k2) - dvars(ivrx,k1)) + dr1old(1)
         dr1new(2) = rlx*(dvars(ivry,k2) - dvars(ivry,k1)) + dr1old(2)
         dr1new(3) = rlx*(dvars(ivrz,k2) - dvars(ivrz,k1)) + dr1old(3)
                                                             
         dr2new(1) = rlx*(dvars(ivrx,k3) - dvars(ivrx,k2)) + dr2old(1)
         dr2new(2) = rlx*(dvars(ivry,k3) - dvars(ivry,k2)) + dr2old(2)
         dr2new(3) = rlx*(dvars(ivrz,k3) - dvars(ivrz,k2)) + dr2old(3)

        else
c------- quad element
         dr1old(1) = vars(ivrx,k3) - vars(ivrx,k1)
         dr1old(2) = vars(ivry,k3) - vars(ivry,k1)
         dr1old(3) = vars(ivrz,k3) - vars(ivrz,k1)

         dr2old(1) = vars(ivrx,k4) - vars(ivrx,k2)
         dr2old(2) = vars(ivry,k4) - vars(ivry,k2)
         dr2old(3) = vars(ivrz,k4) - vars(ivrz,k2)

         dr1new(1) = rlx*(dvars(ivrx,k3) - dvars(ivrx,k1)) + dr1old(1)
         dr1new(2) = rlx*(dvars(ivry,k3) - dvars(ivry,k1)) + dr1old(2)
         dr1new(3) = rlx*(dvars(ivrz,k3) - dvars(ivrz,k1)) + dr1old(3)
                                                             
         dr2new(1) = rlx*(dvars(ivrx,k4) - dvars(ivrx,k2)) + dr2old(1)
         dr2new(2) = rlx*(dvars(ivry,k4) - dvars(ivry,k2)) + dr2old(2)
         dr2new(3) = rlx*(dvars(ivrz,k4) - dvars(ivrz,k2)) + dr2old(3)
        endif

        call cross(dr1old,dr2old,dxdold)
        call cross(dr1new,dr2new,dxdnew)

        dotp = dxdold(1)*dxdnew(1)
     &       + dxdold(2)*dxdnew(2)
     &       + dxdold(3)*dxdnew(3)

        if(dotp .lt. 0.0) then
c------- old and new area vectors point in opposite directions... 
c-        the element likely inverted
          rlx = 0.5*rlx
          ivmax = ivtot+1
          kmax = k1
          go to 50
        endif

        asqold = dxdold(1)**2 + dxdold(2)**2 + dxdold(3)**2
        asqnew = dxdnew(1)**2 + dxdnew(2)**2 + dxdnew(3)**2 

        if(asqnew/asqold .lt. armin**2) then
c-------- area changed by less than a factor of armin
          rlx = 0.5*rlx
          ivmax = ivtot+1
          kmax = k1
        endif

 50     continue
      enddo

      if(rlx .eq. rlxbeg) then
c----- no new change in rlx was made... exit rlx-halving loop
       go to 60
      endif

      enddo ! next ipass

 60   continue


      if(lprint) then
c----- print out Newton iteration info
       if(iter.eq.1 .and. .not.ldfix) then
        write(*,*)
        write(*,*) 'Converging HSM equation system ...'
       endif

       if(ldfix) then
        cfix = '*'
       else
        cfix = ' '
       endif

       if(mod(iter,10).eq.1 .and. .not.ldfix) then
        write(*,*) 
     & ' iter     dr         dd         dp         rlx   max'
c        000*  2345678901 2345678901 2345678901    1.000
       endif
 2000  format(1x,i3,a1, 1x, 3e11.3, 3x, f6.3,
     &      '  d', a, ' @ (',3f7.2,')')

       write(*,2000) iter, cfix, 
     &         max( dvmax(ivrx) , dvmax(ivry) , dvmax(ivrz) ) ,
     &         max( dvmax(ivdx) , dvmax(ivdy) , dvmax(ivdz) ) ,
     &              dvmax(ivps),
     &              rlx, vlabel(ivmax), 
     &           pars(lvr0x,kmax),pars(lvr0y,kmax),pars(lvr0z,kmax)
 2001  format(1x, 'it =',i3, a1,
     &     '  rlx =',f6.3,2x,
     &     '   dr =',e11.3,
     &     '   dd =',e11.3,
     &     '   dp =',e11.3,
     &      2x, a, '@(',3f7.2,')')
      endif

c-----------------------------------------------------------
c-    Newton update
      do k = 1, nnode
        do iv = 1, ivtot
          vars(iv,k) = vars(iv,k) + rlx*dvars(iv,k)
        enddo

c------ renormalize d to unit length
        ad = sqrt( vars(ivdx,k)**2
     &           + vars(ivdy,k)**2
     &           + vars(ivdz,k)**2 )
        vars(ivdx,k) = vars(ivdx,k) / ad
        vars(ivdy,k) = vars(ivdy,k) / ad
        vars(ivdz,k) = vars(ivdz,k) / ad
      enddo

c---- max relative displacement change
      rdel = max( dvmax(ivrx) , 
     &            dvmax(ivry) , 
     &            dvmax(ivrz)   ) / rref

c---- max director (angle) change
      adel = max( dvmax(ivdx) , 
     &            dvmax(ivdy) , 
     &            dvmax(ivdz)  )

c---- max drill angle change
      pdel = dvmax(ivps)

ccc      write(*,*) ldfix, rdel, rtolm, adel, damem

      if(ldfix) then
c----- d is currently fixed...
       if(rdel .lt. rtolm) then
c------ ... and r is weakly converged --> make d free again in next iteration
        ldfix = .false.

c------ skip convergence check, since moment residuals may not be converged,
        if(iter.lt.itmax) go to 900
       endif

      else
c----- d is currently free...
       if(adel .gt. damem .or.
     &    pdel .gt. damem      ) then
c------ big change in d or psi...  freeze them during next iteration
        ldfix = .true.

       endif
      endif

      if(ldfix) then
c----- go do another Newton iteration loop without incrementing iter
       go to 5
      endif

c---- set convergence tolerance flag
      lconv = rdel .lt. rtol .and.
     &        adel .lt. atol .and.
     &        pdel .lt. atol 

c---- jump out of Newton loop if converged on last iteration or limit reached
      if(lconv .or. iter.ge.itmax) go to 901

c      call qpause(ians)
c      if(ians.eq.-1) return

 900  continue  ! next Newton iteration
 901  continue


c      if(nvarg.eq.0) return
c
c      do igres = 1, nvarg
c        n = igres
c        call hsmglr(igres, parg,
c     &              nnode, pars, vars,
c     &              nvarg, varg,
c     &              nelem,kelem,
c     &              kdim,nedim,
c     &              gres,
c     &              nkgres, kgres, 
c     &              gres_vars,
c     &              gres_varg,
c     &              aelem(1,n),
c     &              belem(1,n) )
c        rgl(igres) = gres
c        do lvarg = 1, nvarg
c          agl(igres,lvarg) = gres_varg(lvarg)
c        enddo
cc        write(*,4455) (agl(igres,lvarg), lvarg=1, nvarg), rgl(igres)
cc 4455   format(1x,20f10.5)
c      enddo
c
c
c      call ludcmp(ngdim,nvarg,agl)
c      call baksub(ngdim,nvarg,agl,rgl)
c
c      dvgmax = 0.
c      do ivarg = 1, nvarg
c        dvgmax = max( dvgmax , abs(rgl(ivarg)) )
c      enddo
c
c      do ivarg = 1, nvarg
cc         write(*,*) ivarg, varg(ivarg),  - rgl(ivarg)
c        varg(ivarg) = varg(ivarg) - rgl(ivarg)
c      enddo
c
c      write(*,*) 'dvgmax', dvgmax
c      if(dvgmax .lt. 1.0e-7) return
c
c      do k = 1, nnode
c        fsum(k) = 0.
c        asum(k) = 0.
c      enddo
c
c      vel(1) = parg(lgvelx)
c      vel(2) = parg(lgvely)
c      vel(3) = parg(lgvelz)
c
c      do n = 1, nelem
c        vrel(1) = -vel(1)
c        vrel(2) = -vel(2)
c        vrel(3) = -vel(3)
c
c        ivarg = n
c        gamma = varg(ivarg)
c
c        felem(1,n) = gamma*( vrel(2)*belem(3,n)
c     &                     - vrel(3)*belem(2,n))
c        felem(2,n) = gamma*( vrel(3)*belem(1,n)
c     &                     - vrel(1)*belem(3,n))
c        felem(3,n) = gamma*( vrel(1)*belem(2,n)
c     &                     - vrel(2)*belem(1,n))
c
c        amag = sqrt(aelem(1,n)**2 + aelem(2,n)**2 + aelem(3,n)**2)
c        fnorm = ( felem(1,n)*aelem(1,n)
c     &          + felem(2,n)*aelem(2,n)
c     &          + felem(3,n)*aelem(3,n) ) / amag
c
cc        write(*,*) 'Fn', n, fnorm
c
c        k1 = kelem(1,n)
c        k2 = kelem(2,n)
c        k3 = kelem(3,n)
c        k4 = kelem(4,n)
c        fsum(k1) = fsum(k1) + 0.25*fnorm
c        fsum(k2) = fsum(k2) + 0.25*fnorm
c        fsum(k3) = fsum(k3) + 0.25*fnorm
c        fsum(k4) = fsum(k4) + 0.25*fnorm
c        asum(k1) = asum(k1) + 0.25*amag
c        asum(k2) = asum(k2) + 0.25*amag
c        asum(k3) = asum(k3) + 0.25*amag
c        asum(k4) = asum(k4) + 0.25*amag
c      enddo
c
c      dqnmax = 0
c      do k = 1, nnode
c        if(asum(k) .eq. 0.0) then
c         write(*,*) 'A=0, k=', k
c         stop
c        endif
c
c        qold = pars(lvqn,k)
c        qnew = fsum(k) / asum(k)
c
c        dqn = qnew - qold
c
c        dqnmax = max(dqnmax , abs(dqn))
c
c        pars(lvqn,k) = qnew
cc        write(*,*) 'qnew', qnew
c      enddo
c 
cc      return
cc      pause

      return
      end ! hsmsol




      subroutine jadd(keq1, keq2, niv, ibdim, nddim,
     &             resi_vars, ki,
     &             resc_vars, kdvp, ndvp)
c----------------------------------------------------------------------
c     Sets or adds Jacobian block element to Jacobian block row.
c     If element does not exist, it is created with input element values.
c     If element exists, the input values are added to it.
c
c  Inputs:
c     keq1    store only rows keq1..keq2 in block element
c     keq2
c     niv     number of columns in block element
c     ibdim   dimension of block elements
c     nddim   dimension of number of elements per block row
c     resi_vars(..)  input block element to be stored
c     ki      block index of block element
c
c  Outputs:
c     resc_vars(...)  block row of Jacobian matrix array 
c                      with created or modified block element 
c     kdvp(.)        column numbers of block elements
c     ndvp           number of block elements in this row
c
c----------------------------------------------------------------------

      real resi_vars(ibdim,niv)
      real resc_vars(ibdim,niv,nddim)
      integer kdvp(nddim)

c---- check if node ki Jacobian block entry already exists
c      write(*,*) 'ndvp =', ndvp
      do id = 1, ndvp
c        write(*,*) 'id kdvp', id, kdvp(id)
        if(kdvp(id) .eq. ki) go to 10
      enddo

      if(ndvp .ge. nddim) then
       write(*,*) 'JADD: Too many nodes influencing residual'
       write(*,*) '      Increase parameter nddim in index.inc'
       stop
      endif

c------------------------------------------------------
c---- Jacobian entry does not exist... so create it
      id = ndvp + 1
      ndvp = id
      kdvp(id) = ki

c      write(*,*) 'new id k =', id, ki

c---- set new Jacobian block for specified equations
      do keq = keq1, keq2
        do iv = 1, niv
          resc_vars(keq,iv,id) = resi_vars(keq,iv)
        enddo
      enddo

      return
 
c------------------------------------------------------
c---- pick up here if Jacobian entry already exists
 10   continue

c---- add on Jacobian elements for specified equations
      do keq = keq1, keq2
        do iv = 1, niv
          resc_vars(keq,iv,id) = 
     &    resc_vars(keq,iv,id) + resi_vars(keq,iv)
        enddo
      enddo

      return
      end ! jadd



      subroutine dbcset(neq, ibdim, nddim,
     &            keq1, res1, res1_dvp, kres,
     &            keqp, resp, resp_dvp, kdvp, ndvp)
c-------------------------------------------------------------------------------
c     Uses new Dirichlet system row res1(keq1), res1_dvp(keq1..) to
c     overwrite current  system row resp(keqp), resp_dvp(keqp..) for node kres
c-------------------------------------------------------------------------------
      real res1(ibdim)  , res1_dvp(ibdim,ibdim)
      real resp(ibdim,*), resp_dvp(ibdim,ibdim,nddim,*)
      integer kdvp(nddim,*), ndvp(*)

c---- overwrite residual
      resp(keqp,kres) = res1(keq1)

c---- clear current Jacobians
      do id = 1, ndvp(kres)
        do iv = 1, neq
          resp_dvp(keqp,iv,id,kres) = 0.
        enddo
      enddo

c---- find Jacobian slot id which corresponds to residual node kres
      do id = 1, ndvp(kres)
        if(kdvp(id,kres) .eq. kres) go to 10
      enddo
      write(*,*) '? DBCSET: No diagonal Jacobian element for k =', kres

      if(ndvp(kres) .ge. nddim) then
       write(*,*) 'DBCSET: Too many nodes influencing residual'
       write(*,*) '        Increase parameter nddim'
       stop
      endif

c---- Jacobian entry does not exsts... so create it
      id = ndvp(kres) + 1
      ndvp(kres) = id
      kdvp(id,kres) = kres

c----------------------------------------------------------------
c---- if we arrived here, we have an available Jacobian slot id
 10   continue

c---- install new Jacobian
      do iv = 1, neq
        resp_dvp(keqp,iv,id,kres) = res1_dvp(keq1,iv)
      enddo

      return
      end ! dbcset


      subroutine colset(lmset,
     &                  kdim,nddim,nmdim,
     &                  nbdim,nb, nnode, resp_dvp,kdvp,ndvp,
     &                  namat,ifrst,ilast,mfrst,amat)
c---------------------------------------------------------------------
c     Creates banded matrix of blocks in columns
c
c  Inputs:
c     lmset     T   set pointers and fill matrix array amat
c               F   set pointers only
c     nbdim           block dimension
c     nb              block size
c     nnode           number of nodes (number of block rows and columns)
c     resp_dvp(..dk)  matrix block elements of row k , d = 1 .. ndvp(k)
c     kdvp(dk)        column indices of blocks in row k
c     ndvp(k)         number of blocks in row k
c
c  Outputs:
c     namat       matrix size  ( = nnode )
c     ifrst(j)    i  row  index of topmost element in column j
c     ilast(j)    i  row  index of bottommost element in column j
c     mfrst(j)    m array index of topmost element in column j
c     amat(..m)   block matrix array, sequential, ordered by columns
c
c---------------------------------------------------------------------

c------------------------------------
c---- passed arrays
      logical lmset
      real resp_dvp(nbdim,nbdim,nddim,kdim)
      integer kdvp(nddim,kdim)
      integer ndvp(kdim)

      integer ifrst(kdim),
     &        ilast(kdim),
     &        mfrst(kdim+1)
      real amat(nbdim,nbdim,nmdim)
c------------------------------------

      namat = nnode

c---- initialize row indices for each column as most-optimistic diagonal case
      do j = 1, nnode
        ifrst(j) = j
        ilast(j) = j
      enddo

c---- go over each equation (matrix row)
      do i = 1, nnode
c------ go over Jacobian elements in this equation
        k = i

        jmin = namat + 1
        jmax = 0
        do id = 1, ndvp(k)
c-------- column index of this element
          j = kdvp(id,k)

c-------- adjust min and max indices with this equation's row index
          ifrst(j) = min( i , ifrst(j) )
          ilast(j) = max( i , ilast(j) )

          jmin = min( j , jmin )
          jmax = max( j , jmax )
        enddo

c------ allocate any additional blocks which will be needed for fill-in
        do j = jmin, jmax
          ilast(j) = max( i , ilast(j) )
        enddo
      enddo

c---- set element pointers for linear matrix array
      mfrst(1) = 1
      do j = 1, nnode
        inum = ilast(j) - ifrst(j)
        mfrst(j+1) = mfrst(j) + inum + 1
      enddo

c---- skip matrix fill if lmset == F
      if(.not.lmset) return

      nmmat = mfrst(nnode+1) - 1
      if(nmmat .gt. nmdim) then
       write(*,*)
       write(*,*) 'COLSET: Overflow on matrix-storage array'
       write(*,*) '        Increase  nmdim  to at least', nmmat
       stop
      endif

c---- clear matrix
      do m = 1, nmmat
        do jj = 1, nb
          do ii = 1, nb
            amat(ii,jj,m) = 0.
          enddo
        enddo
      enddo
          
c---- fill matrix

c---- go over each equation (matrix row)
      do i = 1, nnode
        k = i

c------ go over Jacobian elements in this equation
        do id = 1, ndvp(k)
c-------- column index of this element
          j = kdvp(id,k)

c-------- linear array index for this element
          m = mfrst(j) - ifrst(j) + i

          do jj = 1, nb
            do ii = 1, nb
              amat(ii,jj,m) = resp_dvp(ii,jj,id,k)
            enddo
          enddo

c        write(*,*) 
c        write(*,*) id, k
c        do ii = 1, nb
ccc        write(*,'(1x,8g12.4)') (amat(ii,jj,m), jj=1, nb)
c          write(*,'(1x,8g12.4)') (resp_dvp(ii,jj,id,k), jj=1, nb)
c        enddo
c         pause

        enddo
      enddo

      return
      end ! colset


      subroutine bort3(v,b1,b2,b3,ib)
c---------------------------------------------------------------------
c     Resets b1,b2,b3 vectors such that
c     * The b vector closest to v is set to v, noted by ib
c     * From the two remaining b vectors...
c        the b which is nearest normal to v is made normal to v
c     * The one remaining b is set to be normal to the other two b`s
c---------------------------------------------------------------------
      real v(3)
      real b1(3), b2(3), b3(3)

      real vh(3)

      va = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      vh(1) = v(1)/va
      vh(2) = v(2)/va
      vh(3) = v(3)/va

      vb1 = vh(1)*b1(1) + vh(2)*b1(2) + vh(3)*b1(3)
      vb2 = vh(1)*b2(1) + vh(2)*b2(2) + vh(3)*b2(3)
      vb3 = vh(1)*b3(1) + vh(2)*b3(2) + vh(3)*b3(3)

      avb1 = abs(vb1)
      avb2 = abs(vb2)
      avb3 = abs(vb3)

      if(avb1 .ge. max(avb2,avb3)) then
c----- v is closest to b1
       ib = 1
       b1(1) = vh(1)
       b1(2) = vh(2)
       b1(3) = vh(3)

       if(avb2 .lt. avb3) then
c------ set b2 by projection, set b3 = b1 x b2
        b2(1) = b2(1) - vh(1)*vb2
        b2(2) = b2(2) - vh(2)*vb2
        b2(3) = b2(3) - vh(3)*vb2
        call cross(b1,b2,b3)
       else
c------ set b3 by projection, set b2 = b3 x b1
        b3(1) = b3(1) - vh(1)*vb3
        b3(2) = b3(2) - vh(2)*vb3
        b3(3) = b3(3) - vh(3)*vb3
        call cross(b3,b1,b2)
       endif

      elseif(avb2 .ge. max(avb1,avb3)) then
c----- v is closest to b2
       ib = 2
       b2(1) = vh(1)
       b2(2) = vh(2)
       b2(3) = vh(3)

       if(avb1 .lt. avb3) then
c------ set b1 by projection, set b3 = b1 x b2
        b1(1) = b1(1) - vh(1)*vb1
        b1(2) = b1(2) - vh(2)*vb1
        b1(3) = b1(3) - vh(3)*vb1
        call cross(b1,b2,b3)
       else
c------ set b3 by projection, set b1 = b2 x b3
        b3(1) = b3(1) - vh(1)*vb3
        b3(2) = b3(2) - vh(2)*vb3
        b3(3) = b3(3) - vh(3)*vb3
        call cross(b2,b3,b1)
       endif

      else
c----- v is closest to b3
       ib = 3
       b3(1) = vh(1)
       b3(2) = vh(2)
       b3(3) = vh(3)

       if(avb1 .lt. avb2) then
c------ set b1 by projection, set b2 = b3 x b1
        b1(1) = b1(1) - vh(1)*vb1
        b1(2) = b1(2) - vh(2)*vb1
        b1(3) = b1(3) - vh(3)*vb1
        call cross(b3,b1,b2)
       else
c------ set b2 by projection, set b1 = b2 x b3
        b2(1) = b2(1) - vh(1)*vb2
        b2(2) = b2(2) - vh(2)*vb2
        b2(3) = b2(3) - vh(3)*vb2
        call cross(b2,b3,b1)
       endif

      endif

      return
      end ! bort3


      subroutine bset3(v,b1,b2,b3,ib)
c-----------------------------------------------------------------
c     Like bort3, but doesn`t orthogonalize remaining vectors
c-----------------------------------------------------------------
      real v(3)
      real b1(3), b2(3), b3(3)

      real vh(3)

      va = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      vh(1) = v(1)/va
      vh(2) = v(2)/va
      vh(3) = v(3)/va

      vb1 = vh(1)*b1(1) + vh(2)*b1(2) + vh(3)*b1(3)
      vb2 = vh(1)*b2(1) + vh(2)*b2(2) + vh(3)*b2(3)
      vb3 = vh(1)*b3(1) + vh(2)*b3(2) + vh(3)*b3(3)

      avb1 = abs(vb1)
      avb2 = abs(vb2)
      avb3 = abs(vb3)

      if    (avb1 .ge. max(avb2,avb3)) then
c----- v is closest to b1
       ib = 1
       b1(1) = vh(1)
       b1(2) = vh(2)
       b1(3) = vh(3)

      elseif(avb2 .ge. max(avb1,avb3)) then
c----- v is closest to b2
       ib = 2
       b2(1) = vh(1)
       b2(2) = vh(2)
       b2(3) = vh(3)

      else
c----- v is closest to b3
       ib = 3
       b3(1) = vh(1)
       b3(2) = vh(2)
       b3(3) = vh(3)

      endif

      return
      end ! bset3


      subroutine bort2(v,b1,b2,ib)
c--------------------------------------------------------------
c     Resets b1,b2 vectors such that
c     * the b vector closest to v is set to v, noted by ib
c     * the other b vector is orthogonalized against v
c--------------------------------------------------------------
      real v(3)
      real b1(3), b2(3)

      real vh(3)

      va = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      vh(1) = v(1)/va
      vh(2) = v(2)/va
      vh(3) = v(3)/va

      vb1 = vh(1)*b1(1) + vh(2)*b1(2) + vh(3)*b1(3)
      vb2 = vh(1)*b2(1) + vh(2)*b2(2) + vh(3)*b2(3)

      avb1 = abs(vb1)
      avb2 = abs(vb2)

      if(avb1 .ge. avb2) then
c----- v is closest to b1
       ib = 1
       b1(1) = vh(1)
       b1(2) = vh(2)
       b1(3) = vh(3)

c----- set b2 by projection
       b2(1) = b2(1) - vh(1)*vb2
       b2(2) = b2(2) - vh(2)*vb2
       b2(3) = b2(3) - vh(3)*vb2

      else
c----- v is closest to b2
       ib = 2
       b2(1) = vh(1)
       b2(2) = vh(2)
       b2(3) = vh(3)

c----- set b1 by projection, set b2 = b3 x b1
       b1(1) = b1(1) - vh(1)*vb1
       b1(2) = b1(2) - vh(2)*vb1
       b1(3) = b1(3) - vh(3)*vb1
      endif

      return
      end ! bort2



      subroutine bset2(v,b1,b2,ib)
c-----------------------------------------------------------------
c     Like bort2, but doesn`t orthogonalize remaining vector
c-----------------------------------------------------------------
      real v(3)
      real b1(3), b2(3)

      real vh(3)

      va = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      vh(1) = v(1)/va
      vh(2) = v(2)/va
      vh(3) = v(3)/va

      vb1 = vh(1)*b1(1) + vh(2)*b1(2) + vh(3)*b1(3)
      vb2 = vh(1)*b2(1) + vh(2)*b2(2) + vh(3)*b2(3)

      avb1 = abs(vb1)
      avb2 = abs(vb2)

      if(avb1 .ge. avb2) then
c----- v is closest to b1
       ib = 1
       b1(1) = vh(1)
       b1(2) = vh(2)
       b1(3) = vh(3)

      else
c----- v is closest to b2
       ib = 2
       b2(1) = vh(1)
       b2(2) = vh(2)
       b2(3) = vh(3)

      endif

      return
      end ! bset2



