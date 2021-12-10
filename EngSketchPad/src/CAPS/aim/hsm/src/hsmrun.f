
      program hsmrun
c-----------------------------------------------------------------------------
c     Example driver program for exercising the various HSM routines 
c     to solve a shell elasticity problem using the HSM formulation.
c
c  Primary routines called here
c   HSMSOL   primary-variable solution by invoking routines below
c   HSMDEP   dependent variables from primary variable solution
c   HSMOUT   outputs solution on all elements to various plain-text files
c
c  Helper routines used here for case setup
c   ISOMAT   calculated stiffness matrices for isotropic shell
c   ORTMAT   calculated stiffness matrices for orthotropic shell, 
c              with shear/extension couplings
c
c  Utility routines used by primary routines
c   HSMGEO   manifold basis vectors and metric tensors of one element
c   HSMGEO2  manifold basis vectors and metric tensors of one edge
c   HSMABD   element node stiffness matrices A,B,D in element basis
c   HSMEQN   equation residual + Jacobians for one element
c   HSMFBC   BC residual + Jacobian for one edge segment
c   HSMRBC   BC residual + Jacobian for one edge node
c   HSMREN   residuals and Jacobians for n  node basis vectors
c   HSMRE1   residuals and Jacobians for e1 node basis vectors
c   HSMRFM   residuals and Jacobians for eps,kap,f,m node values
c
c
c     All routines above assume an arbitrary unstructured grid.
c     However, this driver program a logically-rectangular i,j surfaces 
c     to construct the unstructured grid arrays.  Sign of the icase 
c     argument specifies a quad or triangle grid.
c
c     Usage:
c         
c       % hsmrun  icase  imf  ni-1 nj-1  itmax  fload nload  lrcurv ldrill
c
c     Arguments ( "-" or missing argument will take on default value ):
c
c       icase = case index  (see below for list of cases implemented)
c               (if icase < 0, use triangular mesh, then ignore sign)
c       imf   = manufactured solution index (ignored if not implemented)
c       ni,nj = grid size 
c       itmax = max number of Newton iterations 
c             = -1  calculate only residuals and post-process initial guess
c       fload = load scaling factor  (default=1.0)
c       nload = number of loading steps (default=1)
c
c-----------------------------------------------------------------------------
      implicit real (a-h,m,o-z)

      include 'index.inc'

c----------------------------------------
c---- Primary array dimensions
      parameter (kdim = 10501)     ! max number of nodes
      parameter (ldim = 750)       ! max number of BC nodes or BC edges
      parameter (nedim = 12000)    ! max number of elements
      parameter (nddim = 18)       ! max number of nodes influencing a residual
      parameter (nmdim = 2100000)  ! max number of matrix block elements

c----------------------------------------

      logical lvinit, leinit
      logical lprint
      logical lrcurv, ldrill

c---- global parameters
      real parg(lgtot)

c---- global variables (in development)
      real varg(nedim)

c---- shell node primary variables
      real vars(ivtot,kdim)

c---- shell node parameters
      real pars(lvtot,kdim)
      real par1(lvtot,kdim)

c---- shell node secondary (dependent) variables
      real deps(jvtot,kdim)

c---- edge boundary condition parameters
      real pare(letot,ldim)
      real pare1(letot,ldim)

c---- node boundary condition parameters
      real parp(lptot,ldim)
      real parp1(lptot,ldim)

c---- element -> node  pointers
      integer kelem(4,nedim)

cc---- perimeter-edge -> node  pointers
c      integer nperim
c      integer kperim(2,ldim)

c---- BC-edge -> node pointers
      integer nbcedge
      integer kbcedge(2,ldim)

c---- BC-node -> node pointers , BC type flags
      integer nbcnode
      integer kbcnode(ldim)
      integer lbcnode(ldim)

c---- joint-node pointers
      integer njoint
      integer kjoint(2,ldim)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc---- common blocks for debugging
c      real resp0(irtot,kdim)
c      real dres(kdim)
c      common / rcom / resc, resp0, resp, ares, dres
c      common / icom / k0, k1, k2, k3, k4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c--------------------------------------------------------------------
c---- work arrays for HSMSOL

c---- residual projection vectors
      real bf(3,3,kdim),
     &     bm(3,3,kdim)

      real bf_dj(3,3,3,kdim),
     &     bm_dj(3,3,3,kdim)

      integer ibr1(ldim), ibr2(ldim), ibr3(ldim),
     &        ibd1(ldim), ibd2(ldim), ibd3(ldim)

c---- residuals and Jacobians in cartesian basis 
      real resc(irtot,kdim)
      real resc_vars(irtot,ivtot,nddim,kdim)

c---- residuals and Jacobians in b1,b2,bn projection basis
      real resp(irtot,kdim)
      real resp_vars(irtot,ivtot,nddim)
      real resp_dvp(irtot,irtot,nddim,kdim)

c---- Jacobian element pointers
      integer kdvp(nddim,kdim)
      integer ndvp(kdim)

      real ares(kdim)

c---- Newton system matrix arrays
      integer ifrst(kdim),
     &        ilast(kdim),
     &        mfrst(kdim+1)
      real amat(irtot,irtot,nmdim)
      integer ipp(irtot,nmdim)

c---- shell primary variable changes
      real dvars(ivtot,kdim)

c--------------------------------------------------------------------
c---- work arrays for HSMDEP
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
     &        jtfrst(kdim),
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
c---- convenient node i,j --> pointers for case setup
      parameter (idim = 401, 
     &           jdim =  81 )
      integer kij(idim,jdim)

c---- pointers for multi-surface case setup
      parameter (nsdim = 5)  ! max number of i,j surfaces
      integer kijs(idim,jdim,nsdim)

      logical lquad

      real ru(3), rv(3), rxr(3)

      real r0j(3,4), rj(3,4), r0c(3), rc(3)
c------------------------------------------------

      logical lconv, lsout

      character*80 arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
      character*1 cdum

      call getarg(1,arg1)
      call getarg(2,arg2)
      call getarg(3,arg3)
      call getarg(4,arg4)
      call getarg(5,arg5)
      call getarg(6,arg6)
      call getarg(7,arg7)
      call getarg(8,arg8)
      call getarg(9,arg9)

      if(arg1 .eq. ' ') then
       write(*,*)
       write(*,*) 'Usage:'
       write(*,*)
       write(*,*) 
     &    '  % hsmrun  icase  imf  ni-1  nj-1  itmax  fload  nload ',
     &    ' lrcurv ldrill'
       write(*,*)
       stop
      endif

      pi = 4.0*atan(1.0)

c----------------------------------------------
      lrcurv = .true.   ! quadratic element surface
c     lrcurv = .false.  ! bilinear element surface
c- - - - - - -
      ldrill = .true.   ! include drilling DOF
c     ldrill = .false.  ! exclude drilling DOF
c----------------------------------------------

c---- default grid type (overridden with sign of icase parameter)
      lquad = .true.    ! quads
c     lquad = .false.   ! triangles

c---- Newton convergence printout flag
c     lprint = .false.
      lprint = .true.

c---- number of element-interior output points per element side (to fort.51)
      nsout = 0
c     nsout = 8

c---- Newton iteration limit for primary-variable solution
      itmaxv = 35

      if(arg5 .ne. '-' .and.
     &   arg5 .ne. ' '       ) then
       read(arg5,*) itmaxv
      endif

c---- Newton iteration limit for dependent-variable solution
      itmaxd = 15

c---- Newton convergence tolerances
c      rtol1 = 1.0e-8  ! relative displacements, |d|/rref
c      atol1 = 1.0e-8  ! angles (unit vector changes)
      rtol1 = 1.0e-11  ! relative displacements, |d|/rref
      atol1 = 1.0e-11  ! angles (unit vector changes)

c---- load ramp-up tolerances
      rtolf = 1.0e-11
      atolf = 1.0e-11

      rtol = rtol1
      atol = atol1

c--------------------------------------------------------------------
c---- d,psi change threshold to trigger membrane-only sub-iterations
      damem = 0.025

c---- displacement convergence tolerance for membrane-only sub-iterations
      rtolm = 1.0e-4

c==============================================================================
c---- default setup for all cases

c---- reference length for displacement limiting, convergence checks
c-     (should be comparable to size of the geometry)
      rref = 1.0

c---- max Newton changes for position (dimensionless), angle
      rlim = 1.0
      alim = 0.2

c---- vector-length/rref ratio for output of unit vectors e1,e2,n
      vfac = 0.2

c---- shell thickness
      tshell = 0.05

c---- relative location of reference surface within shell thickness
c-    +1.0 = top
c-     0.0 = center
c-    -1.0 = bottom
      zetref = 0.0

c---- coefficient for random noise in initial solution guess
c     ranfac = 1.0e-5
      ranfac = 0.0

c---- mass density (only significant in presence of gravity or acceleration)
      rho = 1.0

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c---- Young's modulus and Poisson's ratio for isotropic material ...
c     mate = 1.0e10
      mate = 1.0e7
c     mate = 1.0e6
c     mate = 1.0e5
c     mate = 1.0e4
c     mate = 1.0e3
c     matv = 0.5
c     matv = 0.3
      matv = 0.

c---- ... and corresponding isotropic shear modulus
      matg = mate*0.5/(1.0+matv)

c---- orthotropic material constants which will actually be used
c-    (here we use the isotropic special case for convenience)
      mate1  = mate
      mate2  = mate
      matv12 = matv
      matg12 = matg

      matg13 = matg
      matg23 = matg

c---- in-surface extension-shear coupling stiffnesses 
c-     (modification to orthotropic case)
c     matc16 = 0.1*mate
c     matc26 = 0.2*mate
      matc16 = 0.
      matc26 = 0.

c---- transverse-shear strain energy reduction factor
c-    (corrects for assumed constant transverse shear stress across thickness)
      srfac = 5.0/6.0

c---- default zero global parameters
      do lg = 1, lgtot
        parg(lg) = 0.
      enddo

c---- load scaling factor
      fload = 1.0
      if(arg6 .ne. '-' .and.
     &   arg6 .ne. ' '       ) then
       read(arg6,*) fload
      endif

c---- number of loading steps
      nload = 1
      if(arg7 .ne. '-' .and.
     &   arg7 .ne. ' '       ) then
       read(arg7,*) nload
      endif

      if(arg8 .ne. '-' .and.
     &   arg8 .ne. ' '       ) then
       read(arg8,*) lrcurv
c      lrcurv = ircurv .ne. 0
      endif

      if(arg9 .ne. '-' .and.
     &   arg9 .ne. ' '       ) then
       read(arg9,*) ldrill
c      ldrill = idrill .ne. 0
      endif

c---- number of global variables
      nvarg = 0

c==============================================================================
c---- case setup

c---- default case index for setup below
c     icase = 1  ! x-load stretched plate
c     icase = 2  ! y-load sheared plate (Cook's membrane problem)
c     icase = 3  ! cantilevered plate, left edge clamped
c     icase = 4  ! z-load bent plate, all edges pinned , manufactured solution
c     icase = 5  ! z-load bent plate, all edges clamped, manufactured solution
c     icase = 6  ! same as 5, but with symmetry planes on East and North sides
c     icase = 7  ! sheared plate
c     icase = 8  ! wing
c     icase = 9  ! round plate
c     icase = 10 ! hemisphere
c     icase = 11 ! tube beam
c     icase = 12 ! plate wing
      icase = 13 ! beam with joint

      imf = 0

c---- grid size
c      ni = 61
c      nj = 61

c      ni = 41
c      nj = 41

c      ni = 31
c      nj = 31

       ni = 21
       nj = 21

c      ni = 11
c      nj = 11

      ni = 6
      nj = 6

c      ni = 3
c      nj = 3

      if(arg1 .ne. '-' .and. arg1 .ne. ' ') then
       read(arg1,*) icase
      endif

c---- sign of icase indicates quad or tri grid
      lquad = icase .gt. 0
      icase = iabs(icase)

      if(arg2 .ne. '-' .and. arg2 .ne. ' ') then
       read(arg2,*) imf
      endif

      if(arg3 .ne. '-' .and. arg3 .ne. ' ') then
       read(arg3,*) ni1
       ni = ni1 + 1
      endif

      if(arg4 .ne. '-' .and. arg4 .ne. ' ') then
       read(arg4,*) nj1
       nj = nj1 + 1
      endif

c---- set  i,j  <-->  k  node  pointers for case setup convenience
      k = 0
c---- set node index k ordering for smallest matrix bandwidth
c     do j = 1, nj        !  best if ni < nj
c       do i = 1, ni      !
      do i = 1, ni        !  best if ni > nj
        do j = 1, nj      !
          k = k + 1
          kij(i,j) = k
        enddo
      enddo

      if(icase.eq.12) then
       k = 0
       do j = 1, nj        !  best if ni < nj
         do i = 1, ni      !
           k = k + 1
           kij(i,j) = k
         enddo
       enddo
      endif


c---- total number of nodes
      nnode = k


      if(icase.eq.9) then
c----- circular disk: replace i=ni line with i=1 line
       i = ni
       do j = 1, nj
         kij(i,j) = j
       enddo
c
c----- i=ni nodes no longer exist
       nnode = nnode - nj
      endif


      if(icase.eq.11) then
c------ tube beam: replace j=nj line with j=1 line
        k = 0
        do i = 1, ni 
          do j = 1, nj-1 
            k = k + 1
            kij(i,j) = k
          enddo
          kij(i,nj) = kij(i,1)
        enddo

        nnode = k
      endif


      if(nnode .gt. kdim) then
       write(*,*)
       write(*,*) 'Too many nodes.',
     &            '  Increase  kdim  to at least', nnode
       stop
      endif

c------------------------------------------------------------
c---- set element corner node indices
c-     (element order has no effect on matrix structure)
      nelem = 0
c      do j = 1, nj-1  
c        do i = 1, ni-1
      do i = 1, ni-1
        do j = 1, nj-1
          if(lquad) then
c--------- quads
           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i  ,j  )
           kelem(2,n) = kij(i+1,j  )
           kelem(3,n) = kij(i+1,j+1)
           kelem(4,n) = kij(i  ,j+1)

          else
c--------- triangles
           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i  ,j  )
           kelem(2,n) = kij(i+1,j  )
           kelem(3,n) = kij(i+1,j+1)
           kelem(4,n) = 0

           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i+1,j+1)
           kelem(2,n) = kij(i  ,j+1)
           kelem(3,n) = kij(i  ,j  )
           kelem(4,n) = 0
          endif
        enddo
      enddo

      if(icase.eq.12) then
       nelem = 0
       do j = 1, nj-1  
         do i = 1, ni-1
c       do i = 1, ni-1
c         do j = 1, nj-1
           if(lquad) then
c---------- quads
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = kij(i  ,j+1)

           else
c---------- triangles
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = 0

            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i+1,j+1)
            kelem(2,n) = kij(i  ,j+1)
            kelem(3,n) = kij(i  ,j  )
            kelem(4,n) = 0
           endif
         enddo
       enddo
      endif


      if(icase.eq.13) then
c----- beam with joint between lines  i = ijoint  and  i = ijoint+1
       nelem = 0
       ijoint = ni/2

       do j = 1, nj-1  
         do i = 1, ijoint-1
           if(lquad) then
c---------- quads
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = kij(i  ,j+1)

           else
c---------- triangles
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = 0

            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i+1,j+1)
            kelem(2,n) = kij(i  ,j+1)
            kelem(3,n) = kij(i  ,j  )
            kelem(4,n) = 0
           endif
         enddo

         do i = ijoint+1, ni-1
           if(lquad) then
c---------- quads
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = kij(i  ,j+1)

           else
c---------- triangles
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = 0

            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i+1,j+1)
            kelem(2,n) = kij(i  ,j+1)
            kelem(3,n) = kij(i  ,j  )
            kelem(4,n) = 0
           endif
         enddo
       enddo
      endif

      if(nelem .gt. nedim) then
       write(*,*)
       write(*,*) 'Too many elements.',
     &            '  Increase  nedim  to at least', nelem
       stop
      endif

c---- default zero variables, dependent variables, parameters
      do k = 1, nnode
        do iv = 1, ivtot
          vars(iv,k) = 0.
        enddo
        do jv = 1, jvtot
          deps(jv,k) = 0.
        enddo
        do lv = 1, lvtot
          pars(lv,k) = 0.
        enddo
      enddo

c---- default zero boundary-load (Neumann) BCs
      nbcedge = 0
      do ibc = 1, ldim
        do le = 1, letot
          pare(le,ibc) = 0.
        enddo
      enddo

c---- default no boundary point-load or position (Dirichlet) BCs
      nbcnode = 0
      do ibc = 1, ldim
        lbcnode(ibc) = 0
        do lp = 1, lptot
          parp(lp,ibc) = 0.
        enddo
      enddo

c---- default no joint data
      njoint = 0
      do ij = 1, ldim
        kjoint(1,ij) = 0
        kjoint(2,ij) = 0
      enddo

c---- primary variables are not initialized
      lvinit = .false.

c------------------------------------------------------------
c---- set default local basis vectors, stiffness matrices
      do k = 1, nnode
c------ node basis vectors e01, n0,  e02 = n0 x e01
        pars(lve01x,k) = 1.0
        pars(lve01y,k) = 0.
        pars(lve01z,k) = 0.
c       pars(lve01x,k) = 0.8
c       pars(lve01y,k) = 0.6
c       pars(lve01z,k) = 0.

        pars(lvn0x,k) = 0.
        pars(lvn0y,k) = 0.
        pars(lvn0z,k) = 1.0

        call cross(pars(lvn0x,k), pars(lve01x,k), pars(lve02x,k))

c------ shell mass/area, thickness, reference surface location
        pars(lvmu,k)  = tshell*rho
        pars(lvtsh,k) = tshell
        pars(lvzrf,k) = zetref

c------ set A,B,D,S matrices, store in pars(..) array
        call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lva11,k),pars(lva55,k))
      enddo

c==============================================================================
c---- setup of hard-wired cases
      if    (icase.eq.1) then
c----- end-loaded plate
       xlen = 1.0
       ylen = 0.1

       rlim = 0.1
       alim = 0.1

       mate = 10000.0
       matv = 0.333
       matv = 0.
       matg = mate*0.5/(1.0+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matc16 = 0.
       matc26 = 0.
       matv12 = matv

       matg13 = matg
       matg23 = matg

       tshell = 0.1 * 1.2**(1.0/3.0)

       den = 1.0 - matv**2
       a11 = mate     /den
       a22 = mate     /den
       a12 = mate*matv/den
       a66 = 2.0*matg

       dpar = tshell**3 * mate / (12.0 * (1.0-matv**2))
       write(*,*) 't D', tshell, dpar

       fxtip = 0.
       fytip = 0.
       fztip = 0.

       if(imf.eq.1) then
c------ manufacture area loading qx,qy from assumed deflection
        fxtip = 0.
        dfac = 0.05 * fload

       elseif(imf .eq. 0) then
c------ end load only
        fxtip = 0.0
        fztip = 1.0
        dfac = 0.

       else
        write(*,*) 'No manufactured solution  imf =', imf

        fxtip = 1.0 * fload
        dfac = 0.

       endif

c-------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
           fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           fx = fi - 1.40*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.40*fj*(0.5-fj)*(1.0-fj)

           fx = fi 
           fy = fj

           k = kij(i,j)


           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           x = xlen * fx
           y = ylen * fy

           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = 0.

           if(imf .eq. 1) then
c---------- set assumed displacements u,v
            xt = x/xlen
            yt = y/ylen

            u = xlen * xt
            v = 0.

            ux = xt
            uxx = 0.
            uxy = 0.

            uy = 0.
            uyx = 0.
            uyy = 0.
            
            vx = 0.
            vxx = 0.
            vxy = 0.

            vy = 0.
            vyx = 0.
            vyy = 0.


            exx = ux 
            eyy = vy
            exy = 0.5*(uy + vx)

            exx_x = uxx
            exx_y = uxy

            eyy_x = vyx
            eyy_y = vyy

            exy_x = 0.5*(uyx + vxx)
            exy_y = 0.5*(uyy + vxy)

            fxx   = a11*exx   + a12*eyy
            fxx_x = a11*exx_x + a12*eyy_x

            fyy   = a22*eyy   + a12*exx  
            fyy_y = a22*eyy_y + a12*exx_y

            fxy   = a66*exy
            fxy_x = a66*exy_x
            fxy_y = a66*exy_y

c---------- area loads required for x,y force equilibrium
c-          qx = -(dfxx/dx + dfxy/dy)
c-          qy = -(dfxy/dx + dfyy/dy)
            qx = -( fxx_x + fxy_y )
            qy = -( fxy_x + fyy_y )

            vars(ivrx,k) = x + u*dfac
            vars(ivry,k) = y + v*dfac
            vars(ivrz,k) = 0.

            vars(ivdx,k) = 0.
            vars(ivdy,k) = 0.
            vars(ivdz,k) = 1.0

            pars(lvqx,k) = qx * dfac
            pars(lvqy,k) = qy * dfac

            deps(jve11,k) = exx * dfac
            deps(jve22,k) = eyy * dfac
            deps(jve12,k) = exy * dfac

            deps(jvf11,k) = fxx * dfac
            deps(jvf22,k) = fyy * dfac
            deps(jvf12,k) = fxy * dfac

            lvinit = .true.

           else
c---------- displacements are not initialized
            lvinit = .false.

           endif

         enddo
       enddo

c----------------------------------------------------------------------
       if(imf.eq.0) then
c------ set load on right edge
        i = ni
        do j = 1, nj-1
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          pare(lef1x,ibc) = fload*fxtip / ylen
          pare(lef2x,ibc) = fload*fxtip / ylen

          pare(lef1y,ibc) = fload*fytip / ylen
          pare(lef2y,ibc) = fload*fytip / ylen

          pare(lef1z,ibc) = fload*fztip / ylen
          pare(lef2z,ibc) = fload*fztip / ylen
        enddo

c----------------------------------------------------------------------
       else ! if(.false.) then
c------ set loads on all edges to match stresses of manufactured solution

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on right edge
        i = ni
        do j = 1, nj-1
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on left edge (this will be overwritten with Dirichlet BCs)
        i = 1
        do j = 1, nj-1
          j1 = j+1
          j2 = j
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx = -1.0
          ty =  0.0

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on bottom edge
        j = 1
        do i = 1, ni-1
          i1 = i
          i2 = i+1
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx =  ylen / sqrt(xlen**2 + ylen**2)
          ty = -xlen / sqrt(xlen**2 + ylen**2)

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on top edge
        j = nj
        do i = 1, ni-1
          i1 = i+1
          i2 = i
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx = -ylen / sqrt(xlen**2 + ylen**2)
          ty =  xlen / sqrt(xlen**2 + ylen**2)

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

       endif

c-------------------------------------------------------------
c----- Dirichlet BCs on left edge
       i = 1
       do j = 1, nj
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3 + 30

c------- specified position, r-rBC = (0,0,0)
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- specified tangential directions, d.(tBC1,tBC2) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo


c==============================================================================
      elseif(icase.eq.2) then
c----- plate sheared with end y-load (Cook's membrane problem)
       xlen = 48.0
       ylen0  = 44.0
       ylen1  = 16.0

       rlim = 10.0


       mate = 1.0
       matv = 0.33
c       matv = 0.
       matg = mate*0.5/(1.0+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matc16 = 0.
       matc26 = 0.
       matv12 = matv

       matg13 = matg
       matg23 = matg

       tshell = 1.0

       den = 1.0 - matv**2
       a11 = mate     /den
       a22 = mate     /den
       a12 = mate*matv/den
       a66 = 2.0*matg

       if(imf.eq.0) then
c------ standard Cook problem... end load only
        fytip = 1.0 * fload
        dfac = 0.

       else
c------ manufacture area loading qx,qy from assumed deflection
        fytip = 0.
        dfac = 0.05 * fload

       endif

c-------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
           fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           fx = fi - 1.40*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.80*fj*(0.5-fj)*(1.0-fj)

           fx = fi 
           fy = fj

           k = kij(i,j)


           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           x =  xlen *fx
           y = (ylen0*fx)*(1.0-fy) + (ylen1*fx + ylen0)*fy

           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           if(imf .eq. 0) then
c---------- displacements are not initialized
            lvinit = .false.

           else
c---------- set assumed displacements u,v
            xt = x/xlen
            yt = y/xlen

            u = xlen*( -xt*yt/4.0 )
            v = xlen*( xt**2/2.0 - xt**3/6.0 )

            ux = -yt/4.0
            uxx = 0.
            uxy = -1.0/4.0 / xlen

            uy = -xt/4.0
            uyx = -1.0/4.0 / xlen
            uyy = 0.
            
            vx = xt - xt**2/2.0
            vxx = (1.0 - xt) / xlen
            vxy = 0.

            vy = 0.
            vyx = 0.
            vyy = 0.

            exx = ux 
            eyy = vy
            exy = 0.5*(uy + vx)

            exx_x = uxx
            exx_y = uxy

            eyy_x = vyx
            eyy_y = vyy

            exy_x = 0.5*(uyx + vxx)
            exy_y = 0.5*(uyy + vxy)

            fxx   = a11*exx   + a12*eyy
            fxx_x = a11*exx_x + a12*eyy_x

            fyy   = a22*eyy   + a12*exx  
            fyy_y = a22*eyy_y + a12*exx_y

            fxy   = a66*exy
            fxy_x = a66*exy_x
            fxy_y = a66*exy_y

c---------- area loads required for x,y force equilibrium
c-          qx = -(dfxx/dx + dfxy/dy)
c-          qy = -(dfxy/dx + dfyy/dy)
            qx = -( fxx_x + fxy_y )
            qy = -( fxy_x + fyy_y )

            vars(ivrx,k) = x + u*dfac
            vars(ivry,k) = y + v*dfac
            vars(ivrz,k) = 0.

            vars(ivdx,k) = 0.
            vars(ivdy,k) = 0.
            vars(ivdz,k) = 1.0

            pars(lvqx,k) = qx * dfac
            pars(lvqy,k) = qy * dfac

            deps(jve11,k) = exx * dfac
            deps(jve22,k) = eyy * dfac
            deps(jve12,k) = exy * dfac

            deps(jvf11,k) = fxx * dfac
            deps(jvf22,k) = fyy * dfac
            deps(jvf12,k) = fxy * dfac

            lvinit = .true.
           endif

         enddo
       enddo

c----------------------------------------------------------------------
       if(imf.eq.0) then
c------ set y-load on right edge (Cook's problem)
        i = ni
        do j = 1, nj-1
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          pare(lef1y,ibc) = fytip / ylen1
          pare(lef2y,ibc) = fytip / ylen1
        enddo

c----------------------------------------------------------------------
       else ! if(.false.) then
c------ set loads on all edges to match stresses of manufactured solution

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on right edge
        i = ni
        do j = 1, nj-1
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on left edge (this will be overwritten with Dirichlet BCs)
        i = 1
        do j = 1, nj-1
          j1 = j+1
          j2 = j
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx = -1.0
          ty =  0.0

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on bottom edge
        j = 1
        do i = 1, ni-1
          i1 = i
          i2 = i+1
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx =  ylen0 / sqrt(xlen**2 + ylen0**2)
          ty = -xlen  / sqrt(xlen**2 + ylen0**2)

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on top edge
        j = nj
        do i = 1, ni-1
          i1 = i+1
          i2 = i
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          tx = -ylen1 / sqrt(xlen**2 + ylen1**2)
          ty =  xlen  / sqrt(xlen**2 + ylen1**2)

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          dx = pars(lvr0x,k2) - pars(lvr0x,k1)
          dy = pars(lvr0y,k2) - pars(lvr0y,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty
        enddo

       endif

c-------------------------------------------------------------
c----- Dirichlet BCs on left edge
       i = 1
       do j = 1, nj
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3 + 30

c------- specified position, r-rBC = (0,0,0)
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- specified tangential directions, d.(tBC1,tBC2) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo

c==============================================================================
      elseif(icase.eq.3) then
c----- cantilevered plate
       xlen = 1.0
       ylen = 0.2
       twist = 0.

       rlim = 0.05
       alim = 0.05


       mate = 10000.0
       matv = 0.333
       matv = 0.
       matg = mate*0.5/(1.0+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matc16 = 0.
       matc26 = 0.
       matv12 = matv

       matg13 = matg
       matg23 = matg

       tshell = 0.1 * 1.2**(1.0/3.0)

c--    D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell**3 * mate1 / (12.0 * (1.0-matv12**2))

c----- clear all loads, for different loading setups below
       qzload = 0.
       qnload = 0.

       mxload = 0.
       myload = 0.
       mzload = 0.

       mtload = 0.
       mlload = 0.

       fxload = 0.
       fyload = 0.
       fzload = 0.

       fnload = 0.

       if(imf.eq.0) then
c------ distributed load
        qzload = 1.0
c       qnload = 1.0
        rlim = 0.5
        alim = 0.05

       elseif(imf.eq. 1) then
c------ tip load, buckling at fxload = -pi^2/4 = -2.467
        fxload = -pi**2 / 4.0
        myload = 0.001 * fxload*xlen
        rlim = 0.05
        alim = 0.05

        ranfac = 0.9
        ranfac = 0.

       elseif(imf.eq. -1) then
c------ tip load on arc beam
c        fxload = 1.0
        fzload = 1.0
        myload = 0.
        rlim = 0.1
        alim = 0.05

       elseif(imf.eq. 2 .or.
     &        imf.eq.12      ) then

c------ end bending moment to bend beam into a circle segment
c       myload = -2.0*pi*dpar/xlen * 0.3
        myload = -2.0*pi*dpar/xlen
        ranfac = 0.75
        ranfac = 0.
        rlim = 0.05 
        alim = 0.025

        rlim = 0.5
        alim = 0.5

        write(*,*) 'myload =', myload

       elseif(imf.eq.-2) then
c------ end bending moment to bend circular beam straight
        myload = 2.0*pi*dpar/xlen
        rlim = 0.01
        alim = 0.01

       elseif(imf.eq.3) then
c------ end torsion moment to twist beam into a spiral
        ylen = 0.2
        mxload = dpar/xlen * 0.01

        mate1  = mate1 * 0.001
        mate2  = mate2 * 0.001
        tshell = 0.2

       elseif(imf.eq.4 .or.
     &        imf.eq.5 .or.
     &        imf.eq.6     ) then
c------- twisted beam test case
        xlen = 12.0
        ylen = 1.1
        twist = 0.5*pi

        rref = 10.0
        rlim = 0.1
        alim = 0.1

        tshell = 0.0032
        pload = 1.0

        mate = 29.0e6
        matv = 0.22
c        matv = 0.
        matg = mate*0.5/(1.0+matv)

        mate1  = mate
        mate2  = mate
        matv12 = matv
        matg12 = matg

        matg13 = matg
        matg23 = matg

        if(imf.eq.4) then
c------- y load
         fyload = pload/ylen

        elseif(imf.eq.5) then
         fzload = pload/ylen

        elseif(imf.eq.6) then
         mxload = pload/ylen

        endif

       elseif(imf.eq.7) then

c------ end shear x or y or z load on twisted beam
        xlen = 2.0
        ylen = 0.5
        twist = 0.1 * fload

        rref = 2.0
        rlim = 0.5
        alim = 0.5

       else
        write(*,*) 'Unrecognized imf =', imf

       endif

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)

         if(i.eq.1 .or. i.eq.ni) then
          dfx = 0.
          dfy = 0.
         else
          dfx = ranfac/float(ni-1) * (rand(0) - 0.5)
          dfy = 0.
c           write(*,*) dfx
         endif

         do j = 1, nj
           fj = float(j-1) / float(nj-1)
c- - - - - - - - - - - - - - - - - - - -
c--------- spacing fractions
c          fx = fi - 0.3*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 1.8*fj*(0.5-fj)*(1.0-fj)
cc
c          fx = fi - 0.1*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)
cc
           fx = fi - 0.0*(fi + fi**6 - 2.0*fi**5)
           fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

c           dfx = 0.
c           dfy = 0.
           if(i.ne.1 .and. i.ne.ni) then
            dfx = ranfac/float(ni-1) * (rand(0) - 0.5)
           endif
c           if(j.ne.1 .and. j.ne.nj) then
c            dfy = ranfac/float(nj-1) * (rand(0) - 0.5)
c           endif

c           write(*,*) i, dfx, dfy

           fx = fx + dfx
           fy = fy + dfy

c          dtshell = -0.2*tshell
           dtshell = -0.1*tshell
c          dtshell = -0.08*tshell
c          dtshell = -0.05*tshell
           dtshell = 0.

c          tc0 = (tshell-0.5*dtshell)**3
c          tc1 = (tshell+0.5*dtshell)**3

           tc0 =  tshell-0.5*dtshell
           tc1 =  tshell+0.5*dtshell

           tci = tc0*(1.0-fx) + tc1*fx
           tshelli = tci**(1.0/3.0)

           tshelli = tci

c           rlim = 0.02 
c           alim = 0.02 

           pars(lvmu,k)  = tshelli*rho
           pars(lvtsh,k) = tshelli

c          write(*,*) fx, tshelli, rlim, alim

           if(i.eq.1 .and. j.eq.1) then
            write(*,*) 'E v t D =', mate1, matv12, tshelli, dpar
           endif

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshelli, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- default plate geometry
           x = xlen*fx
           y = ylen*(fy-0.5)
           z = 0.

c--------- default r0
           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = z

c--------- undeformed position r0 for individual cases
           if(imf .eq. 0) then
c---------- uniform loading (no manufactured solution set up)
            pars(lvqz,k) = qzload * fload ! * (fy-0.5)
            pars(lvqn,k) = qnload * fload ! * (fy-0.5)

            lvinit = .false.

           elseif(imf .eq. 99) then
c           elseif(imf .eq. 2) then
c---------- beam with end bending moment rolled into circle
            s = x

            curv = -myload * fload / dpar
            rad = 1.0/curv
            theta = curv * s

            vars(ivrx,k) = rad * sin(theta)
            vars(ivry,k) = pars(lvr0y,k)
            vars(ivrz,k) = rad * (1.0 - cos(theta))

            vars(ivdx,k) = -sin(theta)
            vars(ivdy,k) = 0.0
            vars(ivdz,k) =  cos(theta)

c            deps(jvk11,k) = -curv
c            deps(jvm11,k) = -dpar*curv
c            deps(jvm22,k) = -dpar*curv*matv

            lvinit = .true.
c           lvinit = .false.

           elseif(imf .eq. -1 .or.
     &            imf .eq. -2     ) then
            s = x
            if(imf.eq.-1) then
c----------- arc beam
             curv = pi

            elseif(imf.eq.-2) then
c----------- circular beam straightened with end bending moment
             curv = myload * fload / dpar
            endif
            rad = 1.0/curv
            theta = curv * s

            pars(lvr0x,k) = rad * sin(theta)
            pars(lvr0y,k) = y
            pars(lvr0z,k) = rad * (1.0 - cos(theta))

            pars(lvn0x,k) = -sin(theta)
            pars(lvn0y,k) = 0.0
            pars(lvn0z,k) =  cos(theta)

            pars(lve01x,k) = cos(theta)
            pars(lve01y,k) = 0.0
            pars(lve01z,k) = sin(theta)

            pars(lve02x,k) = 0.
            pars(lve02y,k) = 1.0
            pars(lve02z,k) = 0.
            call cross(pars(lvn0x,k), pars(lve01x,k), pars(lve02x,k))

            if(imf.eq.-2) then
c----------- straighten circular beam for initial guess
             vars(ivrx,k) = x
             vars(ivry,k) = y
             vars(ivrz,k) = 0.
             vars(ivdx,k) = 0.
             vars(ivdy,k) = 0.
             vars(ivdz,k) = 1.0
             lvinit = .true.
            endif


           elseif(imf.eq.3) then
c---------- twisted beam with tip end torsion
            lvinit = .false.

           elseif(imf.eq.4 .or.
     &            imf.eq.5 .or.
     &            imf.eq.6      ) then

c            fx = fi - 0.0*(fi + fi**6 - 2.0*fi**5)
c            fy = fj - 1.5*fj*(0.5-fj)*(1.0-fj)

            ft = fx
            ft_fx = 1.0

c           ft = 3.0*fx**2 - 2.0*fx**3
c           ft_fx = 6.0*fx - 6.0*fx**2

            tw = twist*ft
            sint = sin(tw)
            cost = cos(tw)

            x = xlen*fx
            y = ylen*(fy-0.5)*cost
            z = ylen*(fy-0.5)*sint

            ru(1) = xlen
            rv(1) = 0.

            ru(2) = ylen*(fy-0.5)*(-sint*twist*ft_fx)
            rv(2) = ylen         *cost

            ru(3) = ylen*(fy-0.5)*( cost*twist*ft_fx)
            rv(3) = ylen         *sint

            rua = sqrt(ru(1)**2 + ru(2)**2 + ru(3)**2)
            rva = sqrt(rv(1)**2 + rv(2)**2 + rv(3)**2)

c---------- r0
            pars(lvr0x,k) = x
            pars(lvr0y,k) = y
            pars(lvr0z,k) = z

c---------- e01
            pars(lve01x,k) = ru(1)/rua
            pars(lve01y,k) = ru(2)/rua
            pars(lve01z,k) = ru(3)/rua

c---------- e02
            pars(lve02x,k) = rv(1)/rva
            pars(lve02y,k) = rv(2)/rva
            pars(lve02z,k) = rv(3)/rva

c---------- n0 = e01 x e02
            call cross(pars(lve01x,k),pars(lve02x,k),pars(lvn0x,k))

           elseif(imf.eq.7) then
c            fx = fi - 0.0*(fi + fi**6 - 2.0*fi**5)
c            fy = fj - 1.5*fj*(0.5-fj)*(1.0-fj)

            sint = sin(twist*fx)
            cost = cos(twist*fx)

            x = xlen*fx
            y = ylen*(fy-0.5)*cost
            w = ylen*(fy-0.5)*sint

            wx  = ylen*(fy-0.5)*cost*twist/xlen
            wy  = ylen         *sint      /ylen
            wxy = ylen         *cost*twist/(xlen*ylen)

            wxx = 0.
            wyy = 0.
            wxxx = 0.
            wxyy = 0.
            wxxy = 0.
            wyyy = 0.

c---------- set loading so the KL plate equation is satisified
c           q = dpar*(wxxxx + 2.0*wxxyy + wyyyy)
            q = 0.
            pars(lvqz,k) = q
cc          pars(lvqn,k) = q

c---------- set all other variables for this solution
            phi  = -dpar*(wxx  + wyy )
            phix = -dpar*(wxxx + wxyy)
            phiy = -dpar*(wxxy + wyyy)

            gamx = phix / (srfac*matg*tshell)
            gamy = phiy / (srfac*matg*tshell)

            vars(ivrx,k) =  x
            vars(ivry,k) =  y
            vars(ivrz,k) =  w
            vars(ivdx,k) = -wx + gamx
            vars(ivdy,k) = -wy + gamy
            vars(ivdz,k) = sqrt(1.0 - wx**2 - wy**2)

            deps(jvk11,k) = -wxx
            deps(jvk22,k) = -wyy
            deps(jvk12,k) = -wxy

            deps(jvm11,k) = -dpar*(wxx + matv*wyy)
            deps(jvm22,k) = -dpar*(wyy + matv*wxx)
            deps(jvm12,k) = -dpar*(1.0-matv)*wxy

            deps(jvf1n,k) = phix
            deps(jvf2n,k) = phiy

            deps(jvga1,k) = gamx
            deps(jvga2,k) = gamy

            lvinit = .true.
           endif

         enddo
       enddo

c-------------------------------------------------------------
c----- edge BCs

c------ set applied force and moment on right edge
        do j = 1, nj-1
          i = ni
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

c-------- set linear ramp function wf = fy-0.5 at nodes 1,2
          sint = sin(twist)
          cost = cos(twist)
          r1 = pars(lvr0y,k1)*cost + pars(lvr0z,k1)*sint
          r2 = pars(lvr0y,k2)*cost + pars(lvr0z,k2)*sint
          wf1 = r1 / (0.5*ylen)
          wf2 = r2 / (0.5*ylen)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          pare(lem1x,ibc) = mxload * fload
          pare(lem2x,ibc) = mxload * fload
          pare(lem1y,ibc) = myload * fload
          pare(lem2y,ibc) = myload * fload
          pare(lem1z,ibc) = mzload * fload
          pare(lem2z,ibc) = mzload * fload

          pare(lem1t,ibc) = mtload * fload
          pare(lem2t,ibc) = mtload * fload
          pare(lem1l,ibc) = mlload * fload
          pare(lem2l,ibc) = mlload * fload

          pare(lef1x,ibc) = fxload * fload
          pare(lef2x,ibc) = fxload * fload
          pare(lef1y,ibc) = fyload * fload
          pare(lef2y,ibc) = fyload * fload
          pare(lef1z,ibc) = fzload * fload
          pare(lef2z,ibc) = fzload * fload

c         pare(lef1n,ibc) = fnload * fload * wf1
c         pare(lef2n,ibc) = fnload * fload * wf2
c         pare(lef1z,ibc) = fzload * fload * wf1
c         pare(lef2z,ibc) = fzload * fload * wf2

c          pare(lem1l,ibc) = -fnload * fload * 0.5*wf1**2 * 0.5*ylen
c          pare(lem2l,ibc) = -fnload * fload * 0.5*wf2**2 * 0.5*ylen

c         pare(lef1z,ibc) = fnload * fload * wf1
c         pare(lef2z,ibc) = fnload * fload * wf2
        enddo

        if(imf.eq.7) then
c------- manufactured solution BCs on all edges

         do j = 1, nj-1
c--------- left edge
           i = 1
           j1 = j+1
           j2 = j
 
           k1 = kij(i,j1)
           k2 = kij(i,j2)

           nbcedge = nbcedge+1
           ibc = min( nbcedge , ldim )
           kbcedge(1,ibc) = k1
           kbcedge(2,ibc) = k2

           pare(lem1y,ibc) = -deps(jvm12,k1)
           pare(lem2y,ibc) = -deps(jvm12,k2)

c--------- right edge
           i = ni
           j1 = j
           j2 = j+1
 
           k1 = kij(i,j1)
           k2 = kij(i,j2)

           nbcedge = nbcedge+1
           ibc = min( nbcedge , ldim )
           kbcedge(1,ibc) = k1
           kbcedge(2,ibc) = k2

           pare(lem1y,ibc) = deps(jvm12,k1)
           pare(lem2y,ibc) = deps(jvm12,k2)
         enddo

         do i = 1, ni-1
c--------- bottom edge
           j = 1
           i1 = i
           i2 = i+1

           k1 = kij(i1,j)
           k2 = kij(i2,j)

           nbcedge = nbcedge+1
           ibc = min( nbcedge , ldim )
           kbcedge(1,ibc) = k1
           kbcedge(2,ibc) = k2

           pare(lem1x,ibc) = -deps(jvm12,k1)
           pare(lem2x,ibc) = -deps(jvm12,k2)

c--------- top edge
           j = nj
           i1 = i+1
           i2 = i

           k1 = kij(i1,j)
           k2 = kij(i2,j)

           nbcedge = nbcedge+1
           ibc = min( nbcedge , ldim )
           kbcedge(1,ibc) = k1
           kbcedge(2,ibc) = k2

           pare(lem1x,ibc) = deps(jvm12,k1)
           pare(lem2x,ibc) = deps(jvm12,k2)
         enddo
        endif

c-------------------------------------------------------------

c- - - - - - - - - - - - - - - - - 
c----- SW corner
       i = 1
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c       if(imf.eq.6 .or.
c     &    imf.eq.7       ) then
cc------ pinned end
c        lbcnode(ibc) = 3 + 0
c       endif

c----- specified x,y,z position, r-rBC = (0,0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

c----- shell normal doesn`t move, d.(e01,e02) = (0,0)
       parp(lpt1x,ibc) = pars(lve01x,k)
       parp(lpt1y,ibc) = pars(lve01y,k)
       parp(lpt1z,ibc) = pars(lve01z,k)

       parp(lpt2x,ibc) = pars(lve02x,k)
       parp(lpt2y,ibc) = pars(lve02y,k)
       parp(lpt2z,ibc) = pars(lve02z,k)

c- - - - - - - - - - - - - - - - - 
c----- NW corner
       i = 1
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c       if(imf.eq.6 .or.
c     &    imf.eq.7       ) then
cc------ pinned end
c        lbcnode(ibc) = 2 + 0
c       endif


c----- specified x,z positions (r-r0).(x,z) = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = -1.0
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 0.

       parp(lpn2x,ibc) = 0.
       parp(lpn2y,ibc) = 0.
       parp(lpn2z,ibc) = 1.0

c----- shell normal doesn`t move, d.(e01,e02) = (0,0)
       parp(lpt1x,ibc) = pars(lve01x,k)
       parp(lpt1y,ibc) = pars(lve01y,k)
       parp(lpt1z,ibc) = pars(lve01z,k)

       parp(lpt2x,ibc) = pars(lve02x,k)
       parp(lpt2y,ibc) = pars(lve02y,k)
       parp(lpt2z,ibc) = pars(lve02z,k)


c- - - - - - - - - - - - - - - - - 
c       if(imf.eq.5) then
       if(.false.) then
c----- right edge specified z displacement
       i = ni
       do j = 1, nj
         k = kij(i,j)

         dtwist = fload
         yk = pars(lvr0y,k)
         wk = yk*dtwist
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1

c------- specified z position (r-r0).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k) + wk

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0
       enddo

       elseif(imf.eq.6 .or.
     &        imf.eq.7       ) then
c----- right point restraint
       i = ni
       j = (nj+1)/2
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 1

c----- specified z position (r-r0).z = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 0.
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 1.0
       
       endif

c-------------------------------------------------------------
c----- node BCs on left edge, except for corners
       i = 1
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
c        lbcnode(ibc) = 1 + 30
         lbcnode(ibc) = 3 + 30

c       if(imf.eq.5 .or.
c     &    imf.eq.7      ) then
       if(imf.eq.7) then
c------ pinned end
        lbcnode(ibc) = 1 + 0
       endif

c------- specified z-position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal doesn`t move, d.(e01,e02) = (0,0)
         parp(lpt1x,ibc) = pars(lve01x,k)
         parp(lpt1y,ibc) = pars(lve01y,k)
         parp(lpt1z,ibc) = pars(lve01z,k)

         parp(lpt2x,ibc) = pars(lve02x,k)
         parp(lpt2y,ibc) = pars(lve02y,k)
         parp(lpt2z,ibc) = pars(lve02z,k)
        enddo

c        if(imf.eq.-1) then
c        if(imf.eq.2) then
         if(.false.) then
c------- node BCs on bottom and top edges, except for corners
         do i = 2, ni
           j = 1
           k = kij(i,j)
           nbcnode = nbcnode+1
           ibc = min( nbcnode , ldim )
           kbcnode(ibc) = k
           lbcnode(ibc) = 3 + 20
           parp(lpt1x,ibc) = 0.
           parp(lpt1y,ibc) = 1.0
           parp(lpt1z,ibc) = 0.
           parp(lprx,k) = vars(ivrx,k)
           parp(lpry,k) = vars(ivry,k)
           parp(lprz,k) = vars(ivrz,k)

           j = nj
           k = kij(i,j)
           nbcnode = nbcnode+1
           ibc = min( nbcnode , ldim )
           kbcnode(ibc) = k
           lbcnode(ibc) = 3 + 20
           parp(lpt1x,ibc) = 0.
           parp(lpt1y,ibc) = 1.0
           parp(lpt1z,ibc) = 0.
           parp(lprx,k) = vars(ivrx,k)
           parp(lpry,k) = vars(ivry,k)
           parp(lprz,k) = vars(ivrz,k)

        enddo
       endif
c==============================================================================
      elseif(icase.eq.4) then
c----- z-loaded bent plate, all edges pinned, with imposed edge moment
       xlen = 1.0
       ylen = 1.0
       ainc = 0.   ! inclination angle
       rlim = 0.02

       tshell = 0.022894285 * (1.0-matv**2)**(1.0/3.0) ! * 0.2**(1.0/3.0)

c----- reduce shell thickness to decrease shear effects,
c-     (but membrane effects will inrease)
c      tfac = 0.1
c      tfac = 0.5
       tfac = 1.0
c      tfac = 2.0
c      tfac = 4.0
       tshell = tshell * tfac
       mate = mate / tfac**3
       matg = matg / tfac**3
       mate1  = mate
       mate2  = mate
       matg12 = matg
       matg13 = matg
       matg23 = matg

       qzload = 1.0
       qnload = 0.

       if(imf .lt.0 .or. imf .gt.3) then
        write(*,*) 'Unrecognized imf =', imf
        stop
       endif

c-     D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell**3 * mate / (12.0 * (1.0-matv**2))
       apar = tshell    * mate /         (1.0-matv**2)
       write(*,*) 't D =', tshell, dpar

c----- edge moment loadings (must be zero for manufactured solutions)
       mlload1 = 0.  ! left  edge
       mlload2 = 0.  ! right edge
       mlload3 = 0.  ! bot   edge
       mlload4 = 0.  ! top   edge

c-------------------------------------------------
       if(imf.eq.3) then
c------- Morely rhombic plate
         xlen = 100.0
         ylen = 100.0
         ainc = 60.0 * pi/180.0
         tshell = 1.0
         rref = 100.0

         mate = 1.0e7
         matv = 0.3
         matg = mate*0.5/(1.0+matv)

         mate1  = mate
         mate2  = mate
         matv12 = matv
         matg12 = matg
         matg13 = matg
         matg23 = matg

         qzload = 0.
         qnload = 1.0
       endif

c------------------------------------------------------------
c----- geometry setup
       sina = sin(ainc)
       cosa = cos(ainc)

       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
           if(imf.eq.3) then
c---------- bunch points near edges for Morely plate
            fx = fi - 1.95*fi*(0.5-fi)*(1.0-fi)
            fy = fj - 1.95*fj*(0.5-fj)*(1.0-fj)
c            fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c            fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)
           else
c---------- uniform spacing
            fx = fi
            fy = fj

c           fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c           fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)
           endif

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           x = ylen*fy*sina  + xlen*fx
           y = ylen*fy*cosa

           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           if(imf .eq. 0) then
c---------- uniform loading (no manufactured solution)
            pars(lvqz,k) = fload*qzload 
            pars(lvqn,k) = fload*qnload

            lvinit = .false.

           elseif(imf.eq.1 .or. imf.eq.2) then
c---------- manufactured solution, scaled to get q = fload
c-           at at plate center (fx,fy) = (0.5,0.5)

            if(imf .eq. 1) then
c----------- polynomial manufactured solution
              wfac = (fload/dpar) / 24.0
     &            / (1.0/xlen**4 + 0.75/(xlen**2*ylen**2) + 1.0/ylen**4)

             u     = ( fx -  2.0*fx**3 +      fx**4 )
             ux    = (1.0 -  6.0*fx**2 +  4.0*fx**3 )/xlen
             uxx   = (    - 12.0*fx    + 12.0*fx**2 )/xlen**2
             uxxx  = (    - 12.0       + 24.0*fx    )/xlen**3
             uxxxx =                     24.0        /xlen**4
             uxxxxx = 0.

             v     = ( fy -  2.0*fy**3 +      fy**4 )
             vy    = (1.0 -  6.0*fy**2 +  4.0*fy**3 )/ylen
             vyy   = (    - 12.0*fy    + 12.0*fy**2 )/ylen**2
             vyyy  = (    - 12.0       + 24.0*fy    )/ylen**3
             vyyyy =                     24.0        /ylen**4
             vyyyyy = 0.

            elseif(imf .eq. 2) then
c----------- sinusoidal manufactured solution
             wfac = (fload/dpar)/pi**4
     &            / (1.0/xlen**4 + 2.0/(xlen**2*ylen**2) + 1.0/ylen**4)

             sinx = sin(pi*fx)
             siny = sin(pi*fy)
             cosx = cos(pi*fx)
             cosy = cos(pi*fy)

             u      =  sinx
             ux     =  cosx * pi    / xlen
             uxx    = -sinx * pi**2 / xlen**2
             uxxx   = -cosx * pi**3 / xlen**3
             uxxxx  =  sinx * pi**4 / xlen**4
             uxxxxx =  cosx * pi**5 / xlen**5

             v      =  siny
             vy     =  cosy * pi    / ylen
             vyy    = -siny * pi**2 / ylen**2
             vyyy   = -cosy * pi**3 / ylen**3
             vyyyy  =  siny * pi**4 / ylen**4
             vyyyyy =  cosy * pi**5 / ylen**5

            endif

c---------- set manufactured z deflection and its derivatives
            w      = wfac*u     *v

            wx     = wfac*ux    *v
            wxx    = wfac*uxx   *v
            wxxx   = wfac*uxxx  *v
            wxxxx  = wfac*uxxxx *v
            wxxxxx = wfac*uxxxxx*v

            wy     = wfac*u     *vy
            wyy    = wfac*u     *vyy
            wyyy   = wfac*u     *vyyy
            wyyyy  = wfac*u     *vyyyy
            wyyyyy = wfac*u     *vyyyyy

            wxy    = wfac*ux    *vy
            wxxy   = wfac*uxx   *vy
            wxyy   = wfac*ux    *vyy
            wxxyy  = wfac*uxx   *vyy

            wxxxy  = wfac*uxxx  *vy
            wxyyy  = wfac*ux    *vyyy

            wxxxxy = wfac*uxxxx *vy
            wxyyyy = wfac*ux    *vyyyy
            wxxxyy = wfac*uxxx  *vyy
            wxxyyy = wfac*uxx   *vyyy


            e11 = 0.5*wx**2
            e22 = 0.5*wy**2
            e12 = 0.5*wx*wy

            e11x = wx*wxx
            e11y = wx*wxy
            e22x = wy*wxy
            e22y = wy*wyy
            e12x = 0.5*(wxx*wy + wx*wxy)
            e12y = 0.5*(wxy*wy + wx*wyy)

            f11 = apar*(e11 + matv*e22)
            f22 = apar*(e22 + matv*e11)
            f12 = apar*(1.0-matv)*e12

            f11x = apar*(e11x + matv*e22x)
            f11y = apar*(e11y + matv*e22y)
            f22x = apar*(e22x + matv*e11x)
            f22y = apar*(e22y + matv*e11y)
            f12x = apar*(1.0-matv)*e12x
            f12y = apar*(1.0-matv)*e12y

c---------- set in-plane loading so membrane equation is satisfied
            q1 = -(f11x + f12y)
            q2 = -(f12x + f22y)
            pars(lvqx,k) = q1
            pars(lvqy,k) = q2

c---------- set normal loading so the KL plate equation is satisfied
            qn = dpar*(wxxxx + 2.0*wxxyy + wyyyy)
            pars(lvqz,k) = qn
cc          pars(lvqn,k) = qn

c---------- set all other variables for this solution
            phi  = -dpar*(wxx  + wyy )
            phix = -dpar*(wxxx + wxyy)
            phiy = -dpar*(wxxy + wyyy)

            phixx = -dpar*(wxxxx + wxxyy)
            phiyy = -dpar*(wxxyy + wyyyy)
            phixy = -dpar*(wxxxy + wxyyy)

            phixxx = -dpar*(wxxxxx + wxxxyy)
            phixxy = -dpar*(wxxxxy + wxxyyy)
            phixyy = -dpar*(wxxxyy + wxyyyy)
            phiyyy = -dpar*(wxxyyy + wyyyyy)
            phixxy = -dpar*(wxxxxy + wxxyyy)
            phixyy = -dpar*(wxxxyy + wxyyyy)

            gam1 = phix / (srfac*matg*tshell)
            gam2 = phiy / (srfac*matg*tshell)

            gam1x = phixx / (srfac*matg*tshell)
            gam1y = phixy / (srfac*matg*tshell)
            gam2x = phixy / (srfac*matg*tshell)
            gam2y = phiyy / (srfac*matg*tshell)

            gam1xx = phixxx / (srfac*matg*tshell)
            gam1yy = phixyy / (srfac*matg*tshell)
            gam1xy = phixxy / (srfac*matg*tshell)
            gam2xx = phixxy / (srfac*matg*tshell)
            gam2yy = phiyyy / (srfac*matg*tshell)
            gam2xy = phixyy / (srfac*matg*tshell)

            d1 = -wx + gam1
            d2 = -wy + gam2

            d1x = -wxx + gam1x
            d1y = -wxy + gam1y
            d2x = -wxy + gam2x
            d2y = -wyy + gam2y

            d1xx = -wxxx + gam1xx
            d1yy = -wxyy + gam1yy
            d1xy = -wxxy + gam1xy
            d2xx = -wxxy + gam2xx
            d2yy = -wyyy + gam2yy
            d2xy = -wxyy + gam2xy

            m11 = dpar*(d1x + matv*d2y)
            m22 = dpar*(d2y + matv*d1x)
            m12 = dpar*(1.0-matv)*0.5*(d1y+d2x)

            m11x = dpar*(d1xx + matv*d2y)
            m11y = dpar*(d1xy + matv*d2y)
            m22x = dpar*(d2xy + matv*d1x)
            m22y = dpar*(d2yy + matv*d1x)
            m12x = dpar*(1.0-matv)*0.5*(d1xy+d2xx)
            m12y = dpar*(1.0-matv)*0.5*(d1yy+d2xy)

            t1 = -(m11x + m12y - phix)
            t2 = -(m12x + m22y - phiy)
            pars(lvtx,k) = -t2
            pars(lvty,k) =  t1
            pars(lvtz,k) = 0.

            dzsq = 1.0 - d1**2 - d2**2
            if(dzsq .lt. 0.0) then
             write(*,*) '? dzsq < 0  ,  -wx -wy gamx gamy =', 
     &                   -wx, -wy, gamx, gamy
             dzsq = 0.01
            endif


            vars(ivrx,k) =  x
            vars(ivry,k) =  y
            vars(ivrz,k) =  w
            vars(ivdx,k) = d1
            vars(ivdy,k) = d2
            vars(ivdz,k) = sqrt(dzsq)

            deps(jve11,k) = e11
            deps(jve22,k) = e22
            deps(jve12,k) = e12

            deps(jvf11,k) = f11
            deps(jvf22,k) = f22
            deps(jvf12,k) = f12

            deps(jvk11,k) = -wxx
            deps(jvk22,k) = -wyy
            deps(jvk12,k) = -wxy

            deps(jvm11,k) = -dpar*(wxx + matv*wyy)
            deps(jvm22,k) = -dpar*(wyy + matv*wxx)
            deps(jvm12,k) = -dpar*(1.0-matv)*wxy

            deps(jvf1n,k) = phix
            deps(jvf2n,k) = phiy

            deps(jvga1,k) = gam1
            deps(jvga2,k) = gam2

            lvinit = .true.

           elseif(imf .eq. 3) then
c---------- uniform loading for Morely plate
            pars(lvqz,k) = fload*qzload 
            pars(lvqn,k) = fload*qnload

            lvinit = .false.
           endif

         enddo
       enddo

c-------------------------------------------------------------
c----- loading BC setup

c----- set moment on left edge
       i = 1
       do j = 1, nj-1
         j1 = j+1
         j2 = j

         k1 = kij(i,j1)
         k2 = kij(i,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

c        pare(lem1l,ibc) = mlload1 * fload
c        pare(lem2l,ibc) = mlload1 * fload
         pare(lem1l,ibc) =  deps(jvm11,k1)
         pare(lem2l,ibc) =  deps(jvm11,k2)
         pare(lem1t,ibc) = -deps(jvm12,k1)
         pare(lem2t,ibc) = -deps(jvm12,k2)

         pare(lef1t,ibc) =  deps(jvf11,k1)
         pare(lef2t,ibc) =  deps(jvf11,k2)
         pare(lef1l,ibc) =  deps(jvf12,k1)
         pare(lef2l,ibc) =  deps(jvf12,k2)
         pare(lef1n,ibc) = -deps(jvf1n,k1)
         pare(lef2n,ibc) = -deps(jvf1n,k2)
       enddo

c----- set moment on right edge
       i = ni
       do j = 1, nj-1
         j1 = j
         j2 = j+1

         k1 = kij(i,j1)
         k2 = kij(i,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

c        pare(lem1t,ibc) = mtload2 * fload
c        pare(lem2t,ibc) = mtload2 * fload
         pare(lem1l,ibc) =  deps(jvm11,k1)
         pare(lem2l,ibc) =  deps(jvm11,k2)
         pare(lem1t,ibc) = -deps(jvm12,k1)
         pare(lem2t,ibc) = -deps(jvm12,k2)

         pare(lef1t,ibc) =  deps(jvf11,k1)
         pare(lef2t,ibc) =  deps(jvf11,k2)
         pare(lef1l,ibc) =  deps(jvf12,k1)
         pare(lef2l,ibc) =  deps(jvf12,k2)
         pare(lef1n,ibc) =  deps(jvf1n,k1)
         pare(lef2n,ibc) =  deps(jvf1n,k2)
       enddo

c----- set moment on bottom edge
       j = 1
       do i = 1, ni-1
         i1 = i
         i2 = i+1

         k1 = kij(i1,j)
         k2 = kij(i2,j)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

c        pare(lem1t,ibc) = mtload3 * fload
c        pare(lem2t,ibc) = mtload3 * fload
         pare(lem1l,ibc) =  deps(jvm22,k1)
         pare(lem2l,ibc) =  deps(jvm22,k2)
         pare(lem1t,ibc) =  deps(jvm12,k1)
         pare(lem2t,ibc) =  deps(jvm12,k2)

         pare(lef1l,ibc) =  deps(jvf12,k1)
         pare(lef2l,ibc) =  deps(jvf12,k2)
         pare(lef1t,ibc) =  deps(jvf22,k1)
         pare(lef2t,ibc) =  deps(jvf22,k2)
         pare(lef1n,ibc) = -deps(jvf2n,k1)
         pare(lef2n,ibc) = -deps(jvf2n,k2)
       enddo

c----- set moment on top edge
       j = nj
       do i = 1, ni-1
         i1 = i+1
         i2 = i

         k1 = kij(i1,j)
         k2 = kij(i2,j)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

c        pare(lem1t,ibc) = mtload4 * fload
c        pare(lem2t,ibc) = mtload4 * fload
         pare(lem1l,ibc) =  deps(jvm22,k1)
         pare(lem2l,ibc) =  deps(jvm22,k2)
         pare(lem1t,ibc) =  deps(jvm12,k1)
         pare(lem2t,ibc) =  deps(jvm12,k2)

         pare(lef1l,ibc) =  deps(jvf12,k1)
         pare(lef2l,ibc) =  deps(jvf12,k2)
         pare(lef1t,ibc) =  deps(jvf22,k1)
         pare(lef2t,ibc) =  deps(jvf22,k2)
         pare(lef1n,ibc) =  deps(jvf2n,k1)
         pare(lef2n,ibc) =  deps(jvf2n,k2)
       enddo

c-----------------------------------------------
c----- node BC setup

c- - - - - - - - - - - - - - - - - 
c----- SW corner
       i = 1
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c----- specified x,y,z position (r-r0) = (0,0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

c----- shell normal is in z direction,  d.(x,y) = 0
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- SE corner
       i = ni
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 2 + 30

c----- specified y,z position (r-r0).(y,z) = (0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 0.
       parp(lpn1y,ibc) = -1.0
       parp(lpn1z,ibc) = 0.

       parp(lpn2x,ibc) = 0.
       parp(lpn2y,ibc) = 0.
       parp(lpn2z,ibc) = 1.0

c----- shell normal is in z direction,  d.(x,y) = 0
       parp(lpt1x,ibc) = 1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.


c- - - - - - - - - - - - - - - - - 
c----- NE corner
       i = ni
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 1 + 30

c----- specified z position (r-r0).z = (0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 0.
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 1.0

c----- shell normal is in z direction,  d.(x,y) = 0
       parp(lpt1x,ibc) = 1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = 1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- NW corner
       i = 1
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 1 + 30

c----- specified z position (r-r0).z = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 0.
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 1.0

c----- shell normal is in z direction,  d.(x,y) = 0
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = 1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 

c----- left edge nodes, excluding corners
       i = 1
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 20

c------- specified z position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal is in xz plane,  d.y = 0
         parp(lpt1x,ibc) = -sina
         parp(lpt1y,ibc) = -cosa
         parp(lpt1z,ibc) = 0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- right edge nodes, excluding corners
       i = ni
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 20

c------- specified z position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal is in xz plane,  d.y = 0
         parp(lpt1x,ibc) = sina
         parp(lpt1y,ibc) = cosa
         parp(lpt1z,ibc) = 0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- bottom edge nodes, excluding corners
       j = 1
       do i = 2, ni-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 20

c------- specified z position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal is in yz plane,  d.x = 0
         parp(lpt1x,ibc) = 1.0
         parp(lpt1y,ibc) = 0.
         parp(lpt1z,ibc) = 0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- top edge nodes, excluding corners
       j = nj
       do i = 2, ni-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 20

c------- specified z position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal is in yz plane,  d.x = 0
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) = 0.
         parp(lpt1z,ibc) = 0.
       enddo

c==============================================================================
      elseif(icase.eq.5 .or. 
     &       icase.eq.6      ) then
c----- z-loaded bent plate
c-     icase=5:  all edges clamped 
c-     icase=6:  W,S edges clamped, E,N edges sliding
       xlen = 1.0
       ylen = 1.0

       rlim = 0.5

       tshell = 0.022894285 * (1.0-matv**2)**(1.0/3.0) *0.01**(1.0/3.0)
c       tshell = 0.022894285 * (1.0-matv**2)**(1.0/3.0) *0.1**(1.0/3.0)
c       tshell = 0.022894285 * (1.0-matv**2)**(1.0/3.0) *1.0**(1.0/3.0)

c----- scale shell thickness to decrease or increase shear effects
c      tfac = 0.5
       tfac = 1.0
c      tfac = 2.0
       tshell = tshell * tfac
       mate = mate / tfac**3
       matg = matg / tfac**3
       mate1  = mate
       mate2  = mate
       matg12 = matg
       matg13 = matg
       matg23 = matg

       if(imf .lt.0 .or. imf .gt.2) then
        write(*,*) 'Unrecognized imf =', imf
       endif

c-        D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell**3 * mate / (12.0 * (1.0-matv**2))
       write(*,*) 't D', tshell, dpar

       lvinit = .false.

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         do j = 1, nj
           fi = float(i-1) / float(ni-1)
           fj = float(j-1) / float(nj-1)

           if(icase.eq.6) then
            fi = 0.5*fi
            fj = 0.5*fj
           endif

c--------- spacing fractions
c          fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c          fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           fx = fi - 0.00*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 0.00*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))
           
c--------- reduce in-plane stiffness to make bending dominant
           do lv = lva11, lva66
             pars(lv,k) = pars(lv,k) * 0.1
           enddo

c--------- undeformed position r0
           x = xlen*fx
           y = ylen*fy

           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           if(imf .eq. 0) then
c---------- uniform loading (no manufactured solution)
            pars(lvqz,k) = fload
cc          pars(lvqn,k) = fload

            lvinit = .false.

           elseif(imf.eq.1 .or. imf.eq.2) then
c---------- manufactured solution, scaled to get q = fload
c-           at at plate center (fx,fy) = (0.5,0.5)

            if(imf .eq. 1) then
c----------- polynomial manufactured solution
             wfac = (fload/dpar) 
     &           / (24.0/xlen**4 + 2.0/(xlen**2*ylen**2) + 24.0/ylen**4)

             u     = (    fx**2 -  2.0*fx**3 +      fx**4)
             ux    = (2.0*fx    -  6.0*fx**2 +  4.0*fx**3)/xlen
             uxx   = (2.0       - 12.0*fx    + 12.0*fx**2)/xlen**2
             uxxx  = (          - 12.0       + 24.0*fx   )/xlen**3
             uxxxx = (                         24.0      )/xlen**4

             v     = (    fy**2 -  2.0*fy**3 +      fy**4)
             vy    = (2.0*fy    -  6.0*fy**2 +  4.0*fy**3)/ylen
             vyy   = (2.0       - 12.0*fy    + 12.0*fy**2)/ylen**2
             vyyy  = (          - 12.0       + 24.0*fy   )/ylen**3
             vyyyy = (                         24.0      )/ylen**4

            elseif(imf .eq. 2) then
c----------- sinusoidal manufactured solution
             wfac = (fload/dpar)/pi**4
     &           / (8.0/xlen**4 + 8.0/(xlen**2*ylen**2) + 8.0/ylen**4)

             sinx = sin(pi*fx)
             siny = sin(pi*fy)
             cosx = cos(pi*fx)
             cosy = cos(pi*fy)

             u     = sinx**2
             ux    =  2.0*sinx*cosx           * pi    / xlen
             uxx   =  2.0*(cosx**2 - sinx**2) * pi**2 / xlen**2
             uxxx  = -8.0*cosx*sinx           * pi**3 / xlen**3
             uxxxx = -8.0*(cosx**2 - sinx**2) * pi**4 / xlen**4

             v     = siny**2
             vy    =  2.0*siny*cosy           * pi    / ylen
             vyy   =  2.0*(cosy**2 - siny**2) * pi**2 / ylen**2
             vyyy  = -8.0*cosy*siny           * pi**3 / ylen**3
             vyyyy = -8.0*(cosy**2 - siny**2) * pi**4 / ylen**4

            endif

c---------- set manufactured z deflection and its derivatives
            w     = wfac*u    *v

            wx    = wfac*ux   *v
            wxx   = wfac*uxx  *v
            wxxx  = wfac*uxxx *v
            wxxxx = wfac*uxxxx*v

            wy    = wfac*u    *vy
            wyy   = wfac*u    *vyy
            wyyy  = wfac*u    *vyyy
            wyyyy = wfac*u    *vyyyy

            wxy   = wfac*ux   *vy
            wxxy  = wfac*uxx  *vy
            wxyy  = wfac*ux   *vyy
            wxxyy = wfac*uxx  *vyy

c---------- set loading so the KL plate equation is satisified
            q = dpar*(wxxxx + 2.0*wxxyy + wyyyy)
            pars(lvqz,k) = q
cc          pars(lvqn,k) = q

c---------- set all other variables for this solution
            phi  = -dpar*(wxx  + wyy )
            phix = -dpar*(wxxx + wxyy)
            phiy = -dpar*(wxxy + wyyy)

            pars(lvr0x,k) = x
            pars(lvr0y,k) = y
            pars(lvr0z,k) = 0.

            pars(lvqz,k) = q
cc          pars(lvqn,k) = q

            gamx = phix / (srfac*matg*tshell)
            gamy = phiy / (srfac*matg*tshell)

            vars(ivrx,k) =  x
            vars(ivry,k) =  y
            vars(ivrz,k) =  w
            vars(ivdx,k) = -wx + gamx
            vars(ivdy,k) = -wy + gamy
            vars(ivdz,k) = sqrt(1.0 - wx**2 - wy**2)

            deps(jvk11,k) = -wxx
            deps(jvk22,k) = -wyy
            deps(jvk12,k) = -wxy

            deps(jvm11,k) = -dpar*(wxx + matv*wyy)
            deps(jvm22,k) = -dpar*(wyy + matv*wxx)
            deps(jvm12,k) = -dpar*(1.0-matv)*wxy

            deps(jvf1n,k) = phix
            deps(jvf2n,k) = phiy

            deps(jvga1,k) = gamx
            deps(jvga2,k) = gamy

            lvinit = .true.
           endif

         enddo
       enddo

       i0 = ni/2 + 1
       j0 = nj/2 + 1
       k0 = kij(i0,j0)

c-------------------------------------------------------------
c----- edge BC setup

       if(.false.) then

c----- bottom edge
       j = 1
       do i = 1, ni-1
         i1 = i
         i2 = i+1

         j1 = j
         j2 = j

         k1 = kij(i1,j1)
         k2 = kij(i2,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lem1l,ibc) = deps(jvm22,k1)
         pare(lem2l,ibc) = deps(jvm22,k2)
       enddo

c----- right edge
       i = ni
       do j = 1, nj-1
         i1 = i
         i2 = i

         j1 = j
         j2 = j+1

         k1 = kij(i1,j1)
         k2 = kij(i2,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lem1l,ibc) = deps(jvm11,k1)
         pare(lem2l,ibc) = deps(jvm11,k2)
       enddo

c----- top edge
       j = nj
       do i = ni, 2, -1
         i1 = i
         i2 = i-1

         j1 = j
         j2 = j

         k1 = kij(i1,j1)
         k2 = kij(i2,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lem1l,ibc) = deps(jvm22,k1)
         pare(lem2l,ibc) = deps(jvm22,k2)
       enddo

c----- left edge
       i = 1
       do j = nj, 2, -1
         i1 = i
         i2 = i

         j1 = j
         j2 = j-1

         k1 = kij(i1,j1)
         k2 = kij(i2,j2)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lem1l,ibc) = deps(jvm11,k1)
         pare(lem2l,ibc) = deps(jvm11,k2)
       enddo

       endif

c-------------------------------------------------------------
c----- node Dirichlet BC setup

c- - - - - - - - - - - - - - - - - 
c----- SW corner
       i = 1
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k

c      lbcnode(ibc) = 1 + 30
       lbcnode(ibc) = 3 + 30

c----- specified z position (r-r0).z = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 0.
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 1.0

c----- shell normal is in z direction, d.(x,y) = (0,0)
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- SE corner
       i = ni
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k

       if(icase.eq.5) then
c------ full plate
c       lbcnode(ibc) = 2 + 30
        lbcnode(ibc) = 3 + 30

c------ specified x,z position (r-r0).(x,z) = 0
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 0.
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 1.0

        parp(lpn2x,ibc) = 1.0
        parp(lpn2y,ibc) = 0.
        parp(lpn2z,ibc) = 0.

       else
c------ quarter plate -- on yz  symmetry plane
        lbcnode(ibc) = 2 + 30

c------ specified x,z position (r-r0).(x,z) = (0,0)
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 0.
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 1.0

        parp(lpn2x,ibc) = 1.0
        parp(lpn2y,ibc) = 0.
        parp(lpn2z,ibc) = 0.

       endif

c----- shell normal is in z direction, d.(x,y) = (0,0)
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- NE corner
       i = ni
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k

       if(icase.eq.5) then
c------ full plate -- fixed corner
c       lbcnode(ibc) = 1 + 30
        lbcnode(ibc) = 3 + 30

c------ specified z position (r-r0).z = 0
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 0.
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 1.0

        parp(lpn2x,ibc) = 0.
        parp(lpn2y,ibc) = 1.0
        parp(lpn2z,ibc) = 0.

       else
c------ quarter plate -- double symmetry plane (clamped sliding corner)
        lbcnode(ibc) = 2 + 30

c------ specified x,y position (r-r0).(x,y) = (0,0)
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 1.0
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 0.

        parp(lpn2x,ibc) = 0.
        parp(lpn2y,ibc) = 1.0
        parp(lpn2z,ibc) = 0.

       endif

c----- shell normal is in z direction, d.(x,y) = (0,0)
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- NW corner
       i = 1
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k

       if(icase.eq.5) then
c------ full plate
c       lbcnode(ibc) = 2 + 30
        lbcnode(ibc) = 3 + 30

c------ specified y,z position (r-r0).(y,z) = 0
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 0.
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 1.0

        parp(lpn2x,ibc) = 0.
        parp(lpn2y,ibc) = 1.0
        parp(lpn2z,ibc) = 0.

       else
c------ quarter plate -- on xz  symmetry plane
        lbcnode(ibc) = 2 + 30

c------ specified y,z position (r-r0).(y,z) = (0,0)
        parp(lprx,ibc) = pars(lvr0x,k)
        parp(lpry,ibc) = pars(lvr0y,k)
        parp(lprz,ibc) = pars(lvr0z,k)

        parp(lpn1x,ibc) = 0.
        parp(lpn1y,ibc) = 0.
        parp(lpn1z,ibc) = 1.0

        parp(lpn2x,ibc) = 0.
        parp(lpn2y,ibc) = 1.0
        parp(lpn2z,ibc) = 0.

       endif

c----- shell normal is in z direction, d.(x,y) = (0,0)
       parp(lpt1x,ibc) = -1.0
       parp(lpt1y,ibc) = 0.
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 0.
       parp(lpt2y,ibc) = -1.0
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- left edge nodes, excluding corners
       i = 1
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 30

c------- specified z position (r-r0).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell-normal is in y-z plane,  d.x = 0
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- bottom edge nodes, excluding corners
       j = 1
       do i = 2, ni-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 1 + 30

c------- specified z position (r-r0).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal is in x-z plane,  d.y = 0
         parp(lpt1x,ibc) = 0.
         parp(lpt1y,ibc) = -1.0
         parp(lpt1z,ibc) = 0.

         parp(lpt2x,ibc) = 1.0
         parp(lpt2y,ibc) =  0.
         parp(lpt2z,ibc) =  0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- right edge nodes, excluding corners
       i = ni
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k

         if(icase.eq.5) then
c-------- full plate
          lbcnode(ibc) = 1 + 30

c-------- specified z position (r-r0).z = 0
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

          parp(lpn1x,ibc) = 0.
          parp(lpn1y,ibc) = 0.
          parp(lpn1z,ibc) = 1.0

c-------- shell normal is in y,z plane,  d.x = 0
          parp(lpt1x,ibc) = 1.0
          parp(lpt1y,ibc) = 0.
          parp(lpt1z,ibc) = 0.

          parp(lpt2x,ibc) =  0.
          parp(lpt2y,ibc) = -1.0
          parp(lpt2z,ibc) =  0.

         else
c-------- quarter plate -- on yz  symmetry plane
          lbcnode(ibc) = 1 + 20

c-------- specified x position (r-r0).x = 0
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

          parp(lpn1x,ibc) = 1.0
          parp(lpn1y,ibc) = 0.
          parp(lpn1z,ibc) = 0.

c-------- shell normal is in yz plane,  d.x = 0
          parp(lpt1x,ibc) = 1.0
          parp(lpt1y,ibc) = 0.
          parp(lpt1z,ibc) = 0.

         endif
       enddo

c- - - - - - - - - - - - - - - - - 
c----- top edge nodes, excluding corners
       j = nj
       do i = 2, ni-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k

         if(icase.eq.5) then
c-------- full plate
          lbcnode(ibc) = 1 + 30

c-------- specified z position (r-r0).z = 0
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

          parp(lpn1x,ibc) = 0.
          parp(lpn1y,ibc) = 0.
          parp(lpn1z,ibc) = 1.0

c-------- shell normal is in xz plane,  d.y = 0
          parp(lpt1x,ibc) = 0.
          parp(lpt1y,ibc) = 1.0
          parp(lpt1z,ibc) = 0.

          parp(lpt2x,ibc) = -1.0
          parp(lpt2y,ibc) =  0.
          parp(lpt2z,ibc) =  0.

         else
c-------- quarter plate -- on xz  symmetry plane
          lbcnode(ibc) = 1 + 20

c-------- specified y position (r-r0).y = 0
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

          parp(lpn1x,ibc) = 0.
          parp(lpn1y,ibc) = 1.0
          parp(lpn1z,ibc) = 0.

c-------- shell normal is in xz plane,  d.y = 0
          parp(lpt1x,ibc) = 0.
          parp(lpt1y,ibc) = 1.0
          parp(lpt1z,ibc) = 0.

         endif

       enddo

       if(imf.eq.3) then
c------ drilling moment at center node
        i = (ni+1)/2
        j = (nj+1)/2
        k = kij(i,j)

        nbcnode = nbcnode+1
        ibc = min( nbcnode , ldim )

        kbcnode(ibc) = k


        lbcnode(ibc) = 1000

        parp(lpmx,ibc) = 0.
        parp(lpmy,ibc) = 0.
        parp(lpmz,ibc) = 1.0 * fload

       endif

c       call solget(nnode1,vars,pars)
c       lvinit = .true.


c       if(icase.eq.6) then
       if(.false.) then
        do i = 1, ni
          do j = 1, nj
            k = kij(i,j)
            kold = (i-1)*((ni-1)*2+1) + j
            do iv = 1, ivtot
              vars(iv,k) = vars(iv,kold)
            enddo
            do lv = 1, lvtot
              pars(lv,k) = pars(lv,kold)
            enddo
          enddo
        enddo
       endif

c========================================================================
      elseif(icase.eq.7) then
c----- sheared plate
       xlen = 10.0
       ylen = 1.0

       rlim = 0.5


       mate = 10.0
       matv = 0.333
       matv = 0.
       matg = mate*0.5/(1.0+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matc16 = 0.
       matc26 = 0.
       matv12 = matv

       matg13 = matg
       matg23 = matg

       tshell = 1.0

       den = 1.0 - matv**2
       a11 = mate     /den
       a22 = mate     /den
       a12 = mate*matv/den
       a66 = 2.0*matg

       if(imf.eq.0) then
c------ top load only
        fxtop = 1.0 * fload
        dfac = 0.

       elseif(imf.eq.1) then
c------ manufacture area loading qx,qy from assumed deflection
        fxtop = 0.
        dfac = 0.5 * fload

       else
        write(*,*) 'No manufactured solution  imf =', imf

       endif

c-------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
           fx = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           fx = fi - 1.40*fi*(0.5-fi)*(1.0-fi)
           fy = fj - 1.40*fj*(0.5-fj)*(1.0-fj)

           fx = fi 
           fy = fj

           k = kij(i,j)


           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           x = xlen * fx ! - dfac*ylen*fy
           y = ylen * fy

           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           if(imf .eq. 0) then
c---------- displacements are not initialized
            lvinit = .false.

           else
c---------- set assumed displacements u,v
            xt = x/xlen
            yt = y/ylen

            u = ylen * yt
            v = 0.

            ux = 0.
            uxx = 0.
            uxy = 0.

            uy =  1.0
            uyx = 0.
            uyy = 0.
            
            vx = 0.
            vxx = 0.
            vxy = 0.

            vy = 0.
            vyx = 0.
            vyy = 0.


            exx = ux 
            eyy = vy
            exy = 0.5*(uy + vx)

            exx_x = uxx
            exx_y = uxy

            eyy_x = vyx
            eyy_y = vyy

            exy_x = 0.5*(uyx + vxx)
            exy_y = 0.5*(uyy + vxy)

            fxx   = a11*exx   + a12*eyy
            fxx_x = a11*exx_x + a12*eyy_x

            fyy   = a22*eyy   + a12*exx  
            fyy_y = a22*eyy_y + a12*exx_y

            fxy   = a66*exy
            fxy_x = a66*exy_x
            fxy_y = a66*exy_y

c---------- area loads required for x,y force equilibrium
c-          qx = -(dfxx/dx + dfxy/dy)
c-          qy = -(dfxy/dx + dfyy/dy)
            qx = -( fxx_x + fxy_y )
            qy = -( fxy_x + fyy_y )

            vars(ivrx,k) = x + u*dfac
            vars(ivry,k) = y + v*dfac
            vars(ivrz,k) = 0.

            vars(ivdx,k) = 0.
            vars(ivdy,k) = 0.
            vars(ivdz,k) = 1.0

            pars(lvqx,k) = qx * dfac
            pars(lvqy,k) = qy * dfac

            deps(jve11,k) = exx * dfac
            deps(jve22,k) = eyy * dfac
            deps(jve12,k) = exy * dfac

            deps(jvf11,k) = fxx * dfac
            deps(jvf22,k) = fyy * dfac
            deps(jvf12,k) = fxy * dfac

            lvinit = .true.
           endif

         enddo
       enddo

c----------------------------------------------------------------------
       if(imf.eq.0) then
c------ set x-load on top edge
        j = nj
        do i = 1, ni-1
          i1 = i+1
          i2 = i
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          pare(lef1x,ibc) = fxtop
          pare(lef2x,ibc) = fxtop
        enddo

c----------------------------------------------------------------------
       else ! if(.false.) then
c------ set loads on all edges to match stresses of manufactured solution

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on right edge
        i = ni
        do j = 1, nj-1
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty

c          write(*,*) 'R', pare(lef1x,ibc), pare(lef1y,ibc)
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on left edge
        i = 1
        do j = 1, nj-1
          j1 = j+1
          j2 = j
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty

c          write(*,*) 'L', pare(lef1x,ibc), pare(lef1y,ibc)
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on bottom edge  (this will be overwritten with Dirichlet BCs)
        j = 1
        do i = 1, ni-1
          i1 = i
          i2 = i+1
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty

c          write(*,*) 'B', pare(lef1x,ibc), pare(lef1y,ibc)
        enddo

c- - - - - - - - - - - - - - - - - - - - 
c------ set load on top edge
        j = nj
        do i = 1, ni-1
          i1 = i+1
          i2 = i
 
          k1 = kij(i1,j)
          k2 = kij(i2,j)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          dx = vars(ivrx,k2) - vars(ivrx,k1)
          dy = vars(ivry,k2) - vars(ivry,k1)
          tx =  dy/sqrt(dx**2 + dy**2)
          ty = -dx/sqrt(dx**2 + dy**2)

          fxx = deps(jvf11,k1)
          fyy = deps(jvf22,k1)
          fxy = deps(jvf12,k1)
          pare(lef1x,ibc) = fxx*tx + fxy*ty
          pare(lef1y,ibc) = fxy*tx + fyy*ty

          fxx = deps(jvf11,k2)
          fyy = deps(jvf22,k2)
          fxy = deps(jvf12,k2)
          pare(lef2x,ibc) = fxx*tx + fxy*ty
          pare(lef2y,ibc) = fxy*tx + fyy*ty

c          write(*,*) 'T', pare(lef1x,ibc), pare(lef1y,ibc)
        enddo

       endif

c-------------------------------------------------------------
c----- Dirichlet BCs on SW corner
       j = 1
       i = 1

         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3 + 30

c------- specified position, r-rBC = (0,0,0)
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- shell normal is vertical,  d.(x,y) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.

c-------------------------------------------------------------
c----- Dirichlet BCs on bottom edge
       j = 1
       do i = 2, ni
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k

         if(imf.eq.0) then
          lbcnode(ibc) = 3 + 30

c-------- specified position, r-rBC = (0,0,0)
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

         else
          lbcnode(ibc) = 2 + 30

c-------- specified y,z position, (r-rBC)(y,z) = 0,0
          parp(lprx,ibc) = pars(lvr0x,k)
          parp(lpry,ibc) = pars(lvr0y,k)
          parp(lprz,ibc) = pars(lvr0z,k)

          parp(lpt1x,ibc) =  0.
          parp(lpt1y,ibc) = -1.0
          parp(lpt1z,ibc) =  0.

          parp(lpt2x,ibc) =  0.
          parp(lpt2y,ibc) =  0.
          parp(lpt2z,ibc) = -1.0

         endif

c------- shell normal is vertical,  d.(x,y) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo

c-------------------------------------------------------------
c----- Dirichlet BCs on bottom edge
       if(imf.eq.0) then

       j = nj
       do i = 1, ni
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 2 + 30

c------- specified y,z position, (r-rBC)(y,z) = 0,0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpt1x,ibc) =  0.
         parp(lpt1y,ibc) = -1.0
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) =  0.
         parp(lpt2z,ibc) = -1.0

c------- shell normal is vertical,  d.(x,y) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo

       endif

c==============================================================================
      elseif(icase.eq.8) then
c----- wing

c----- i-major ordering (each chord grouped) 
c-      likely to give smaller bandwidth for this case
       k = 0
       do j = 1, nj   
         do i = 1, ni 
           k = k + 1
           kij(i,j) = k
         enddo
       enddo
       nnode = k

       nelem = 0
       do j = 1, nj-1  
         do i = 1, ni-1
           if(lquad) then
c---------- quads
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = kij(i  ,j+1)

           else
c---------- triangles
            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i  ,j  )
            kelem(2,n) = kij(i+1,j  )
            kelem(3,n) = kij(i+1,j+1)
            kelem(4,n) = 0

            nelem = nelem + 1
            n = min( nelem , nedim )
            kelem(1,n) = kij(i+1,j+1)
            kelem(2,n) = kij(i  ,j+1)
            kelem(3,n) = kij(i  ,j  )
            kelem(4,n) = 0
           endif
         enddo
       enddo


       span  = 5.0 
       croot = 1.0
       ctip  = 0.6
       toc   = 0.30
       dxtip = 0.2

       rlim = 0.2

c      tshell = 0.01
       tshell = 0.025
c      tshell = 0.05

       mate = 4.0e5
       matv = 0.3
c      matv = 0.0

       matg = mate*0.5/(1.0+matv)

c----- orthotropic material constants which will actually be used
c-     (here we use the isotropic special case for convenience)
       mate1  = mate
       mate2  = mate
       matv12 = matv
       matg12 = matg

       matg13 = matg
       matg23 = matg

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)
c- - - - - - - - - - - - - - - - - - - -
c--------- spacing fractions
           u = fi + 0.3*fi*(0.5-fi)*(1.0-fi)
           v = fj - 0.5*(fj + fj**6 - 2.0*fj**5)

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

           dx   = dxtip*v
           dx_v = dxtip
           cy   = (ctip-croot)*v + croot
           cy_v =  ctip-croot

           x = (1.0-sin(pi*u))*cy + dx
           y =  v * span
           z = -sin(2.0*pi*u) *cy * 0.5*toc

           ru(1) = -cos(pi*u)*cy*pi
           ru(2) =  0.
           ru(3) = -cos(2.0*pi*u)*cy*0.5*toc * 2.0*pi

           rv(1) = (1.0-sin(pi*u))*cy_v + dx_v
           rv(2) = span
           rv(3) = -sin(2.0*pi*u) *cy_v * 0.5*toc

           rua = sqrt(ru(1)**2 + ru(2)**2 + ru(3)**2)
           rva = sqrt(rv(1)**2 + rv(2)**2 + rv(3)**2)

c--------- r0
           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = z

c--------- e01
           pars(lve01x,k) = ru(1)/rua
           pars(lve01y,k) = ru(2)/rua
           pars(lve01z,k) = ru(3)/rua

c--------- e02
           pars(lve02x,k) = rv(1)/rva
           pars(lve02y,k) = rv(2)/rva
           pars(lve02z,k) = rv(3)/rva

c--------- n0 = e01 x e02
           call cross(pars(lve01x,k),pars(lve02x,k),pars(lvn0x,k))

c--------- normal-loading Cp distribution
           cp = sin(2.0*pi*u) - sin(pi*u) + 2.0*sin(pi*u)**32
           pars(lvqn,k) = -cp * fload

         enddo
       enddo

c-------------------------------------------------------------
c----- set joints all along TE
       njoint = 0
       do j = 2, nj
         k1 = kij( 1,j)
         k2 = kij(ni,j)

         njoint = njoint + 1
         ij = min( njoint , ldim )
         kjoint(1,ij) = k1
         kjoint(2,ij) = k2
       enddo

       if(njoint .gt. ldim) then
        write(*,*)
        write(*,*) 'Too many joints.',
     &             '  Increase  ldim  to at least', njoint
        stop
       endif

c-------------------------------------------------------------
c----- node BCs on bottom edge (wing root)
       j = 1
       do i = 1, ni
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3 + 30
c        lbcnode(ibc) = 3 

c------- specified position, r-rBC = (0,0,0)
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- shell normal is fixed, d.(e1,e2) = 0
         parp(lpt1x,ibc) = pars(lve01x,k)
         parp(lpt1y,ibc) = pars(lve01y,k)
         parp(lpt1z,ibc) = pars(lve01z,k)

         parp(lpt2x,ibc) = pars(lve02x,k)
         parp(lpt2y,ibc) = pars(lve02y,k)
         parp(lpt2z,ibc) = pars(lve02z,k)
        enddo

c==============================================================================
      elseif(icase.eq.9) then
c----- z-loaded round plate, all edges pinned, with imposed edge moment
       rout = 1.0
       rinn = 0.02
       rlim = 0.1

       tshell = 0.022894285 * (1.0-matv**2)**(1.0/3.0) ! * 0.2**(1.0/3.0)

c----- scale shell thickness to decrease or increase shear effects,
c-     (but membrane effects will go the other way)
       tfac = 0.2
c      tfac = 0.5
c      tfac = 1.0
c      tfac = 2.0

       tshell = tshell * tfac
       mate = mate / tfac**3
       matg = matg / tfac**3
       mate1  = mate
       mate2  = mate
       matg12 = matg
       matg13 = matg
       matg23 = matg

       qzload = 0.
       qnload = 0.1

c      qsload = 1.0
       qsload = 0.

c-     D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell**3 * mate / (12.0 * (1.0-matv**2))
c      write(*,*) 't D =', tshell, dpar

c----- edge moment loadings (must be zero for manufactured solutions)
       mlload = 0.
       mtload = 0.

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
c          ft = fi - 1.95*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.95*fj*(0.5-fj)*(1.0-fj)

c          ft = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           ft = fi - 0.00*fi*(0.5-fi)*(1.0-fi)
           fr = fj - 0.00*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           r = rout*(1.0-fr) + rinn*fr
           t = 2.0*pi*ft

           pars(lvr0x,k) = r*cos(t)
           pars(lvr0y,k) = r*sin(t)
           pars(lvr0z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           pars(lve01x,k) = -sin(t)
           pars(lve01y,k) =  cos(t)
           pars(lve01z,k) = 0.

           pars(lve02x,k) = -cos(t)
           pars(lve02y,k) = -sin(t)
           pars(lve02z,k) = 0.

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = 1.0

           pars(lvqz,k) = qzload * fload
           pars(lvqn,k) = qnload * fload

           if(r .le. 0.5) then
            pars(lvqx,k) = qsload * fload * (-r*sin(t))
            pars(lvqy,k) = qsload * fload * ( r*cos(t))
           endif

           lvinit = .false.

         enddo
       enddo

c-------------------------------------------------------------
c----- loading BC setup

c----- set moment on outer edge
       j = 1
       do i = 1, ni-1
         i1 = i
         i2 = i+1

         k1 = kij(i1,j)
         k2 = kij(i2,j)

         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lem1l,ibc) = mlload * fload
         pare(lem2l,ibc) = mlload * fload
         pare(lem1t,ibc) = mtload * fload
         pare(lem2t,ibc) = mtload * fload
       enddo

c-------------------------------------------------------------
c----- set joints all along seam
c       njoint = 0
c       do j = 2, nj
c         k1 = kij( 1,j)
c         k2 = kij(ni,j)
cc
c         njoint = njoint + 1
c         ij = min( njoint , ldim )
c         kjoint(1,ij) = k1
c         kjoint(2,ij) = k2
c       enddo
cc
c       if(njoint .gt. ldim) then
c        write(*,*)
c        write(*,*) 'Too many joints.',
c     &             '  Increase  ldim  to at least', njoint
c        stop
c       endif

c-----------------------------------------------
c----- node BC setup

c----- bottom edge (outer rim) nodes
       j = 1
       do i = 1, ni
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3

c------- specified z position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)
       enddo

c==============================================================================
      elseif(icase.eq.10) then
c----- hemispherical shell (MacNeal, Harder)
       rref = 10.0

       rout = 10.0
       thole = 18.0*pi/180.0

       rlim = 0.1
       alim = 0.1

       tshell = 0.04
       rho = 1.0


       mate = 6.825e7
       matv = 0.3
       matg = mate*0.5/(1+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matg13 = matg
       matg23 = matg

       frload = 1.0 * fload

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
c          ft = fi - 1.95*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.95*fj*(0.5-fj)*(1.0-fj)

c          ft = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

           fraci = fi - 0.00*fi*(0.5-fi)*(1.0-fi)
           fracj = fj - 0.00*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           psi = 0.5*pi*fraci           ! longitude
           phi = (0.5*pi-thole)*fracj   ! latitude

           pars(lvr0x,k) = rout*cos(phi)*cos(psi)
           pars(lvr0y,k) = rout*cos(phi)*sin(psi)
           pars(lvr0z,k) = rout*sin(phi)

           pars(lvn0x,k) = cos(phi)*cos(psi)
           pars(lvn0y,k) = cos(phi)*sin(psi)
           pars(lvn0z,k) = sin(phi)

           pars(lve01x,k) = -sin(psi)
           pars(lve01y,k) =  cos(psi)
           pars(lve01z,k) = 0.

           call cross(pars(lvn0x,k), pars(lve01x,k), pars(lve02x,k))

           lvinit = .false.

         enddo
       enddo

c-------------------------------------------------------------
c----- loading setup
       i = 1
       j = 1
       k = kij(i,j)
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 100

       parp(lpfx,ibc) = frload
       parp(lpfy,ibc) = 0.
       parp(lpfz,ibc) = 0.


       i = ni
       j = 1
       k = kij(i,j)
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 100

       parp(lpfx,ibc) = 0.
       parp(lpfy,ibc) = -frload
       parp(lpfz,ibc) = 0.

c-------------------------------------------------------------
c----- node Dirichlet BC setup

c- - - - - - - - - - - - - - - - - 
c----- left edge nodes
       i = 1
       do j = 1, nj
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k

c------- specified y position (r-r0).y = 0
         lbcnode(ibc) = 1 + 20
         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = -1.0
         parp(lpn1z,ibc) = 0.

         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- shell-normal is in x-z plane direction,  d.y = 0
         parp(lpt1x,ibc) =  0.
         parp(lpt1y,ibc) = -1.0
         parp(lpt1z,ibc) =  0.
       enddo

c- - - - - - - - - - - - - - - - - 
c----- right edge nodes
       i = ni
       do j = 1, nj
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k

         if(j.eq.1) then
c-------- specified x,z positions (r-r0).(x,z) = 0
          lbcnode(ibc) = 2 + 20
          parp(lpn2x,ibc) = 0.
          parp(lpn2y,ibc) = 0.
          parp(lpn2z,ibc) = 1.0
         else
c-------- specified x position (r-r0).x = 0
          lbcnode(ibc) = 1 + 20
         endif
          
c------- specified x position (r-r0).x = 0
         parp(lpn1x,ibc) = -1.0
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 0.

         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- shell-normal is in y-z plane direction,  d.x = 0
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0. 
         parp(lpt1z,ibc) =  0.
       enddo

c==============================================================================
      elseif(icase.eq.11) then
c----- tube beam
       rout = 0.5
       xlen = 10.0
       rlim = 0.5
       alim = 0.2

       tshell = 0.05

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matg13 = matg
       matg23 = matg

c----- incipient wrinkling with qzload = 1.2  (61 x 51 grid)
       qzload = fload * 0.
       qnload = fload * 0.

       fxload = fload * 0.
       fyload = fload * 0.
       fzload = fload * 1.0/pi

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         fi = float(i-1) / float(ni-1)
         do j = 1, nj
           fj = float(j-1) / float(nj-1)

c--------- spacing fractions
c          ft = fi - 1.95*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.95*fj*(0.5-fj)*(1.0-fj)

c          ft = fi - 1.85*fi*(0.5-fi)*(1.0-fi)
c          fr = fj - 1.85*fj*(0.5-fj)*(1.0-fj)

c          fx = fi - 0.98*fi*(1.0-fi)
           fx = fi - 0.70*fi*(1.0-fi)
           ft = fj - 0.00*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- undeformed position r0
           x = xlen*fx
           t = 2.0*pi*ft
           r = rout

           pars(lvr0x,k) = x
           pars(lvr0y,k) = r*cos(t)
           pars(lvr0z,k) = r*sin(t)

           pars(lvn0x,k) = 0.
           pars(lvn0y,k) = cos(t)
           pars(lvn0z,k) = sin(t)

           pars(lve01x,k) = 1.0
           pars(lve01y,k) = 0.
           pars(lve01z,k) = 0.

           pars(lve02x,k) = 0.
           pars(lve02y,k) = -sin(t)
           pars(lve02z,k) =  cos(t)


           pars(lvqz,k) = qzload * fload
           pars(lvqn,k) = qnload * fload


           lvinit = .false.

         enddo
       enddo

c-----------------------------------------------
c----- edge loading BC setup

c----- set load on right edge
       i = ni
       do j = 1, nj-1
         j1 = j
         j2 = j+1

         k1 = kij(i,j1)
         k2 = kij(i,j2)
 
         nbcedge = nbcedge+1
         ibc = min( nbcedge , ldim )
         kbcedge(1,ibc) = k1
         kbcedge(2,ibc) = k2

         pare(lef1x,ibc) = fxload
         pare(lef1y,ibc) = fyload
         pare(lef1z,ibc) = fzload
       enddo

c-----------------------------------------------
c----- node BC setup

c----- left edge nodes
       i = 1
       do j = 1, nj
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
         lbcnode(ibc) = 3 
c        lbcnode(ibc) = 3 + 30

c------- specified position, r-rBC = (0,0,0)
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

c------- specified tangential directions, d.(tBC1,tBC2) = (0,0)
         parp(lpt1x,ibc) = -1.0
         parp(lpt1y,ibc) =  0.
         parp(lpt1z,ibc) =  0.

         parp(lpt2x,ibc) =  0.
         parp(lpt2y,ibc) = -1.0
         parp(lpt2z,ibc) =  0.
       enddo

c==============================================================================
      elseif(icase.eq.12) then
c----- plate wing
       span = 10.0
       chord = 1.0
       sweep = 10.0
       twist = 0.
       taper = 0.7

       alpha = fload

       sina = sin(alpha*pi/180.0)
       cosa = cos(alpha*pi/180.0)

       parg(lgvelx) = -cosa
       parg(lgvely) = 0.
       parg(lgvelz) = -sina

       nvarg = nelem
       qnload = 0.
       qzload = 0.

c       nvarg = 0.
c       qnload = 0.
c       qzload = 0.1*fload

       do ivarg = 1, nvarg
          varg(ivarg) = 0.
       enddo

       rlim = 0.02
       alim = 0.02

       tshell0 = 0.2

c----- set to get D = 1
c       tshell0 = 2.2894285 * ((1.0-matv**2)/mate)**(1.0/3.0)

c-     D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell0**3 * mate / (12.0 * (1.0-matv**2))
       write(*,*) 't D =', tshell, dpar, mate

c------------------------------------------------------------
c----- geometry setup
       do j = 1, nj
         fj = float(j-1) / float(nj-1)

       do i = 1, ni
         fi = float(i-1) / float(ni-1)

c- - - - - - - - - - - - - - - - - - - -
c--------- spacing fractions
c          fx = fi - 0.3*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 1.8*fj*(0.5-fj)*(1.0-fj)
cc
c          fx = fi - 0.1*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)
cc
           fx = fi - 0.0*(fi + fi**6 - 2.0*fi**5)
           fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

           tshell = tshell0

           pars(lvmu,k)  = tshell*rho
           pars(lvtsh,k) = tshell

c          write(*,*) fx, tshelli, rlim, alim

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshell, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- default geometry
           cloc = chord*(1.0-taper*fj)
           y = span*fy
           x = cloc*fx + y*tan(sweep*pi/180.0)
           z = 0.

c--------- default r0
           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = z

           pars(lvqn,k) = qnload
           pars(lvqz,k) = qzload

         enddo
       enddo


       lvinit = .false.

c-------------------------------------------------------------
c- - - - - - - - - - - - - - - - - 
c----- SW corner
       i = 1
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c       if(imf.eq.6 .or.
c     &    imf.eq.7       ) then
cc------ pinned end
c        lbcnode(ibc) = 3 + 0
c       endif

c----- specified x,y,z position, r-rBC = (0,0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 1.0
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 0.

       parp(lpn2x,ibc) = 0.
       parp(lpn2y,ibc) = 1.0
       parp(lpn2z,ibc) = 0.

       parp(lpt1x,ibc) = 0.
       parp(lpt1y,ibc) = 1.0
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 1.0
       parp(lpt2y,ibc) = 0.
       parp(lpt2z,ibc) = 0.

c- - - - - - - - - - - - - - - - - 
c----- SE corner
       i = ni
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30


c----- specified x,z positions (r-r0).(x,z) = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = 1.0
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 0.

       parp(lpn2x,ibc) = 0.
       parp(lpn2y,ibc) = 1.0
       parp(lpn2z,ibc) = 0.

       parp(lpt1x,ibc) = 0.
       parp(lpt1y,ibc) = 1.0
       parp(lpt1z,ibc) = 0.

       parp(lpt2x,ibc) = 1.0
       parp(lpt2y,ibc) = 0.
       parp(lpt2z,ibc) = 0.

c-------------------------------------------------------------
c----- node BCs on bottom edge, except for corners
       j = 1
       do i = 2, ni-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
c        lbcnode(ibc) = 1 + 30
         lbcnode(ibc) = 3 + 30

c------- specified z-position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 1.0
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 0.

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 1.0
         parp(lpn1z,ibc) = 0.
  
         parp(lpt1x,ibc) = 0.
         parp(lpt1y,ibc) = 1.0
         parp(lpt1z,ibc) = 0.

         parp(lpt2x,ibc) = 1.0
         parp(lpt2y,ibc) = 0.
         parp(lpt2z,ibc) = 0.
        enddo

c==============================================================================
      elseif(icase.eq.13) then
c----- cantilevered plate with joint
       xlen1 = 0.5
       xlen2 = 0.5
       ylen = 0.2
       angle = 30.0 * pi/180
       angle = 1. * pi/180
       angle = 0. * pi/180

       rlim = 0.05
       alim = 0.05


       mate = 10000.0
       matv = 0.333
       matv = 0.
       matg = mate*0.5/(1.0+matv)

       mate1  = mate
       mate2  = mate
       matg12 = matg
       matc16 = 0.
       matc26 = 0.
       matv12 = matv

       matg13 = matg
       matg23 = matg

       tshell = 0.1 * 1.2**(1.0/3.0)

c--    D = t^3 E / [ 12 (1-v^2) ]
       dpar = tshell**3 * mate1 / (12.0 * (1.0-matv12**2))

c----- clear all loads, for different loading setups below
       qzload = 0.
       qnload = 0.

       mxload = 0.
       myload = 0.
       mzload = 0.

       mtload = 0.
       mlload = 0.

       fxload = 0.
       fyload = 0.
       fzload = 0.

       fnload = 0.

       if(imf.eq.0) then
c------ distributed load
        qzload = 1.0
c       qnload = 1.0
        rlim = 0.5
        alim = 0.05

       elseif(imf.eq. 1) then
c------ tip load, buckling at fxload = -pi^2/4 = -2.467
        fxload = -pi**2 / 4.0
        myload = 0.001 * fxload*(xlen1+xlen2)
        rlim = 0.05
        alim = 0.05

       elseif(imf.eq. 2) then
c------ end bending moment to bend beam into a circle segment
        myload = -2.0*pi*dpar/(xlen1+xlen2)

        rlim = 0.5
        alim = 0.5

        write(*,*) 'myload =', myload

       elseif(imf.eq.3) then
c------ end torsion moment to twist beam into a spiral
        mxload = dpar/(xlen1+xlen2) * 0.01

       else
        write(*,*) 'Unrecognized imf =', imf

       endif

c------------------------------------------------------------
c----- geometry setup
       do i = 1, ni
         if(i .le. ijoint) then
          fi = float(i-1) / float(ijoint-1)
         else
          fi = float(i-ijoint-1) / float(ni-ijoint-1)
         endif

         do j = 1, nj
           fj = float(j-1) / float(nj-1)
c- - - - - - - - - - - - - - - - - - - -
c--------- spacing fractions
c          fx = fi - 0.3*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 1.8*fj*(0.5-fj)*(1.0-fj)
cc
c          fx = fi - 0.1*(fi + fi**6 - 2.0*fi**5)
c          fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)
cc
           fx = fi - 0.0*(fi + fi**6 - 2.0*fi**5)
           fy = fj - 0.0*fj*(0.5-fj)*(1.0-fj)

           k = kij(i,j)

c          dtshell = -0.2*tshell
           dtshell = -0.1*tshell
c          dtshell = -0.08*tshell
c          dtshell = -0.05*tshell
           dtshell = 0.

c          tc0 = (tshell-0.5*dtshell)**3
c          tc1 = (tshell+0.5*dtshell)**3

           tc0 =  tshell-0.5*dtshell
           tc1 =  tshell+0.5*dtshell

           tci = tc0*(1.0-fx) + tc1*fx
           tshelli = tci**(1.0/3.0)

           tshelli = tci

c           rlim = 0.02 
c           alim = 0.02 

           pars(lvmu,k)  = tshelli*rho
           pars(lvtsh,k) = tshelli

c          write(*,*) fx, tshelli, rlim, alim

           if(i.eq.1 .and. j.eq.1) then
            write(*,*) 'E v t D =', mate1, matv12, tshelli, dpar
           endif

           call ortmat(mate1, mate2, matg12, matv12, 
     &              matc16, matc26,
     &              matg13, matg23,
     &              tshelli, zetref, srfac,
     &              pars(lva11,k),pars(lvb11,k),
     &              pars(lvd11,k),pars(lva55,k))

c--------- default plate geometry
           if(i .le. ijoint) then
            x = xlen1*fx
            y = ylen*(fy-0.5)
            z = 0.
            anx = 0.
            anz = 1.0
           else
            x = xlen1 + xlen2*fx * cos(angle)
            y = ylen*(fy-0.5)
            z =         xlen2*fx * sin(angle)
            anx = -sin(angle)
            anz =  cos(angle)
           endif

c--------- r0
           pars(lvr0x,k) = x
           pars(lvr0y,k) = y
           pars(lvr0z,k) = z

           pars(lvn0x,k) = anx
           pars(lvn0y,k) = 0.
           pars(lvn0z,k) = anz

c--------- undeformed position r0 for individual cases
           if(imf .eq. 0) then
c---------- uniform loading 
            pars(lvqz,k) = qzload * fload ! * (fy-0.5)
            pars(lvqn,k) = qnload * fload ! * (fy-0.5)

            lvinit = .false.
           endif

         enddo
       enddo

c-------------------------------------------------------------
c----- edge BCs

c------ set applied force and moment on right edge
        do j = 1, nj-1
          i = ni
          j1 = j
          j2 = j+1
 
          k1 = kij(i,j1)
          k2 = kij(i,j2)

c-------- set linear ramp function wf = fy-0.5 at nodes 1,2
          sint = sin(twist)
          cost = cos(twist)
          r1 = pars(lvr0y,k1)*cost + pars(lvr0z,k1)*sint
          r2 = pars(lvr0y,k2)*cost + pars(lvr0z,k2)*sint
          wf1 = r1 / (0.5*ylen)
          wf2 = r2 / (0.5*ylen)

          nbcedge = nbcedge+1
          ibc = min( nbcedge , ldim )
          kbcedge(1,ibc) = k1
          kbcedge(2,ibc) = k2

          pare(lem1x,ibc) = mxload * fload
          pare(lem2x,ibc) = mxload * fload
          pare(lem1y,ibc) = myload * fload
          pare(lem2y,ibc) = myload * fload
          pare(lem1z,ibc) = mzload * fload
          pare(lem2z,ibc) = mzload * fload

          pare(lem1t,ibc) = mtload * fload
          pare(lem2t,ibc) = mtload * fload
          pare(lem1l,ibc) = mlload * fload
          pare(lem2l,ibc) = mlload * fload

          pare(lef1x,ibc) = fxload * fload
          pare(lef2x,ibc) = fxload * fload
          pare(lef1y,ibc) = fyload * fload
          pare(lef2y,ibc) = fyload * fload
          pare(lef1z,ibc) = fzload * fload
          pare(lef2z,ibc) = fzload * fload

c         pare(lef1n,ibc) = fnload * fload * wf1
c         pare(lef2n,ibc) = fnload * fload * wf2
c         pare(lef1z,ibc) = fzload * fload * wf1
c         pare(lef2z,ibc) = fzload * fload * wf2

c          pare(lem1l,ibc) = -fnload * fload * 0.5*wf1**2 * 0.5*ylen
c          pare(lem2l,ibc) = -fnload * fload * 0.5*wf2**2 * 0.5*ylen

c         pare(lef1z,ibc) = fnload * fload * wf1
c         pare(lef2z,ibc) = fnload * fload * wf2
        enddo
c-------------------------------------------------------------

c- - - - - - - - - - - - - - - - - 
c----- SW corner
       i = 1
       j = 1
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c       if(imf.eq.6 .or.
c     &    imf.eq.7       ) then
cc------ pinned end
c        lbcnode(ibc) = 3 + 0
c       endif

c----- specified x,y,z position, r-rBC = (0,0,0)
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

c----- shell normal doesn`t move, d.(e01,e02) = (0,0)
       parp(lpt1x,ibc) = pars(lve01x,k)
       parp(lpt1y,ibc) = pars(lve01y,k)
       parp(lpt1z,ibc) = pars(lve01z,k)

       parp(lpt2x,ibc) = pars(lve02x,k)
       parp(lpt2y,ibc) = pars(lve02y,k)
       parp(lpt2z,ibc) = pars(lve02z,k)

c- - - - - - - - - - - - - - - - - 
c----- NW corner
       i = 1
       j = nj
       k = kij(i,j)
         
       nbcnode = nbcnode+1
       ibc = min( nbcnode , ldim )

       kbcnode(ibc) = k
       lbcnode(ibc) = 3 + 30

c       if(imf.eq.6 .or.
c     &    imf.eq.7       ) then
cc------ pinned end
c        lbcnode(ibc) = 2 + 0
c       endif


C----- Specified x,z positions (r-r0).(x,z) = 0
       parp(lprx,ibc) = pars(lvr0x,k)
       parp(lpry,ibc) = pars(lvr0y,k)
       parp(lprz,ibc) = pars(lvr0z,k)

       parp(lpn1x,ibc) = -1.0
       parp(lpn1y,ibc) = 0.
       parp(lpn1z,ibc) = 0.

       parp(lpn2x,ibc) = 0.
       parp(lpn2y,ibc) = 0.
       parp(lpn2z,ibc) = 1.0

c----- shell normal doesn`t move, d.(e01,e02) = (0,0)
       parp(lpt1x,ibc) = pars(lve01x,k)
       parp(lpt1y,ibc) = pars(lve01y,k)
       parp(lpt1z,ibc) = pars(lve01z,k)

       parp(lpt2x,ibc) = pars(lve02x,k)
       parp(lpt2y,ibc) = pars(lve02y,k)
       parp(lpt2z,ibc) = pars(lve02z,k)

c-------------------------------------------------------------
c----- node BCs on left edge, except for corners
       i = 1
       do j = 2, nj-1
         k = kij(i,j)
         
         nbcnode = nbcnode+1
         ibc = min( nbcnode , ldim )

         kbcnode(ibc) = k
c        lbcnode(ibc) = 1 + 30
         lbcnode(ibc) = 3 + 30

c------- specified z-position, (r-rBC).z = 0
         parp(lprx,ibc) = pars(lvr0x,k)
         parp(lpry,ibc) = pars(lvr0y,k)
         parp(lprz,ibc) = pars(lvr0z,k)

         parp(lpn1x,ibc) = 0.
         parp(lpn1y,ibc) = 0.
         parp(lpn1z,ibc) = 1.0

c------- shell normal doesn`t move, d.(e01,e02) = (0,0)
         parp(lpt1x,ibc) = pars(lve01x,k)
         parp(lpt1y,ibc) = pars(lve01y,k)
         parp(lpt1z,ibc) = pars(lve01z,k)

         parp(lpt2x,ibc) = pars(lve02x,k)
         parp(lpt2y,ibc) = pars(lve02y,k)
         parp(lpt2z,ibc) = pars(lve02z,k)
       enddo

c-------------------------------------------------------------
c----- set joints between i = ijoint and i = ijoint+1 lines
       njoint = 0
       do j = 1, nj
         k1 = kij(ijoint  ,j)
         k2 = kij(ijoint+1,j)

         njoint = njoint + 1
         ij = min( njoint , ldim )
         kjoint(1,ij) = k1
         kjoint(2,ij) = k2
       enddo

       if(njoint .gt. ldim) then
        write(*,*)
        write(*,*) 'Too many joints.',
     &             '  Increase  ldim  to at least', njoint
        stop
       endif


c==============================================================================
      else
       write(*,*) 'Undefined icase index', icase
       stop

      endif
c==============================================================================

      if(nbcedge .gt. ldim) then
       write(*,*) 'Array overflow on nbcedge.',
     &            ' Increase ldim to at least', nbcedge
       stop
      endif

      if(nbcnode .gt. ldim) then
       write(*,*) 'Array overflow on nbcnode.',
     &            ' Increase ldim to at least', nbcnode
       stop
      endif


c- - - - - - - - - - - - - - - 
      if(itmaxv .eq. -2) then
c---- "dummy" run to determine required size of  amat,ipp  work arrays
c-     (cannot use parameter nmdim in call list here -- use nmdim1 instead)
       nmdim1 = 1
       call hsmsol(lvinit, lprint, 
     &             lrcurv, ldrill,
     &            -2    , rref,
     &            rlim, rtol, rdel, 
     &            alim, atol, adel,
     &            damem,rtolm,
     &            parg,
     &            nnode,pars,vars,
     &            nvarg,varg,
     &            nelem,kelem,
     &            nbcedge,kbcedge,pare,
     &            nbcnode,kbcnode,parp,lbcnode,
     &            njoint,kjoint,
     & kdim,ldim,nedim,nddim, nmdim1,
     & bf, bf_dj,
     & bm, bm_dj,
     & ibr1,ibr2,ibr3,ibd1,ibd2,ibd3,
     & resc, resc_vars,
     & resp, resp_vars, resp_dvp,
     & kdvp, ndvp,
     & ares,
     & ifrst,ilast,mfrst,
     & amat, ipp, dvars )

       write(*,*)
       write(*,*) 'Minimum required  nmdim =', nmdim1
       write(*,*)
       stop
      endif
c- - - - - - - - - - - - - - - 

      do k = 1, nnode
        do lv = 1, lvtot
          par1(lv,k) = pars(lv,k)
        enddo
      enddo
      
c      write(*,*) 'nbcedge', nbcedge, ldim

      do ibc = 1, nbcedge
        do le = 1, letot
          pare1(le,ibc) = pare(le,ibc)
        enddo
      enddo

      do ibc = 1, nbcnode
        do lp = 1, lptot
          parp1(lp,ibc) = parp(lp,ibc)
        enddo
      enddo

      write(*,*) 'lrcurv =', lrcurv
      write(*,*) 'ldrill =', ldrill

c---- ramp up loading in nload steps
      do 800 iload = 1, nload
        frac = float(iload)/float(nload)

        write(*,*) 
        write(*,*) 'fload =', fload*frac

        do k = 1, nnode
          do lv = lvq1, lvtz
            pars(lv,k) = par1(lv,k) * frac
          enddo
        enddo

        do ibc = 1, nbcedge
          do le = 1, letot
            pare(le,ibc) = pare1(le,ibc) * frac
          enddo
        enddo

        do ibc = 1, nbcnode
          do lp = lpfx, lpmz
            parp(lp,ibc) = parp1(lp,ibc) * frac
          enddo
        enddo
       
        if(iload .eq. nload) then
         rtol = rtol1
         atol = atol1
        else
         rtol = rtolf
         atol = atolf
        endif

        write(*,*) 'rlim alim', rlim, alim


c---- solve HSM system for primary variables vars(..)
      call hsmsol(lvinit, lprint, 
     &            lrcurv, ldrill,
     &            itmaxv, rref,
     &            rlim, rtol, rdel,
     &            alim, atol, adel,
     &            damem,rtolm,
     &            parg,
     &            nnode,pars,vars,
     &            nvarg,varg,
     &            nelem,kelem,
     &            nbcedge,kbcedge,pare,
     &            nbcnode,kbcnode,parp,lbcnode,
     &            njoint,kjoint,
     & kdim,ldim,nedim,nddim,nmdim,
     & bf, bf_dj,
     & bm, bm_dj,
     & ibr1,ibr2,ibr3,ibd1,ibd2,ibd3,
     & resc, resc_vars,
     & resp, resp_vars, resp_dvp,
     & kdvp, ndvp,
     & ares,
     & ifrst,ilast,mfrst,
     & amat, ipp, dvars )

       lvinit = .true.
 800  continue

      
      if(itmaxv .lt. 0) then
       write(*,*) 'Skipping HSMDEP'
       
      else
c----- calculate dependent variables deps(..)
       leinit = .false.
       itmaxe = 20
       elim = 1.0
       etol = atol
       call hsmdep(leinit, lprint, 
     &             lrcurv, ldrill,
     &            itmaxe, 
     &            elim,etol,edel,
     &            nnode,pars,vars,deps,
     &            nelem,kelem,
     & kdim,ldim,nedim,nddim,nmdim,
     & acn,
     & resn, 
     & rese, rese_de,
     & rest, rest_t,
     & kdt, ndt,
     & ifrstt,ilastt,mfrstt,
     & amatt,
     & resv, resv_v,
     & kdv, ndv,
     & ifrstv,ilastv,mfrstv,
     & amatv )

      endif

      if(lprint) write(*,*)

 3140 format(1x,a,f12.4)
 3150 format(1x,a,f12.5)
 3160 format(1x,a,f12.6)
 3170 format(1x,a,f12.7)
 3180 format(1x,a,f13.8)
 3190 format(1x,a,f14.9)
 3200 format(1x,a,f15.10)

      if(icase.eq.2) then
c----- report displacement at middle of tip edge for Cook's problem
       i = ni
       j = (nj+1)/2
       k = kij(i,j)
       write(*,3200) 'tip dy =', (vars(ivry,k) - pars(lvr0y,k))/fload
      endif

      if(icase.eq.3) then
c----- report displacement at middle of tip edge for MacNeal's problem
       i = ni
       j = (nj+1)/2
       k = kij(i,j)
       ang0 = atan2( pars(lve02z,k) , pars(lve02y,k) )
       ang  = atan2( deps(jve2z ,k) , deps(jve2y ,k) )
       write(*,3200) 'tip  x =', vars(ivrx,k)
       write(*,3200) 'tip dx =', vars(ivrx,k) - pars(lvr0x,k)
       write(*,3200) 'tip dy =', vars(ivry,k) - pars(lvr0y,k)
       write(*,3200) 'tip dz =', vars(ivrz,k) - pars(lvr0z,k)
       write(*,3180) 'tip ax =', (ang-ang0)*180.0/pi
      endif

      if(icase.eq.4 .or. icase.eq.5) then
c----- report center displacement
       i = (ni+1)/2
       j = (nj+1)/2
       k = kij(i,j)
       write(*,3200) 'center dz =', vars(ivrz,k) - pars(lvr0z,k)
      endif

      if(icase.eq.6) then
c----- report corner displacement
       i = ni
       j = nj
       k = kij(i,j)
       write(*,3200) 'corner dz =', vars(ivrz,k) - pars(lvr0z,k)
      endif

      if(icase.eq.1 .or.
     &   icase.eq.2 ) then
c----- report displacement at tip
       i = ni
       j = (nj+1)/2
       k = kij(i,j)
       write(*,3200) 'tip dx =', vars(ivrx,k) - pars(lvr0x,k)
       write(*,3200) 'tip dz =', vars(ivrz,k) - pars(lvr0z,k)
      endif

      if(icase.eq.10) then
c----- report displacement at middle of tip edge for MacNeal's problem
       i = 1
       j = 1
       k = kij(i,j)
       write(*,3200) 'corner dx =', vars(ivrx,k) - pars(lvr0x,k)
       write(*,3200) 'corner dy =', vars(ivry,k) - pars(lvr0y,k)
       write(*,3200) 'corner dz =', vars(ivrz,k) - pars(lvr0z,k)
      endif

      if(icase.eq.12) then
c----- report displacement at tip
       i = (ni+1)/2
       j = nj
       k = kij(i,j)
       write(*,3200) 'tip dz =', vars(ivrz,k) - pars(lvr0z,k)
      endif

      if(lprint) write(*,*)

c---- dump solution to text files
      call hsmout(nelem,kelem,vars,deps,pars,parg,
     &            kdim,ldim,nedim,nddim,nmdim,
     &            idim,jdim,ni,nj,kij)

      if(nsout .gt. 0) then

      if(.false.) then
c     if(.true.) then
c----- change quads to triangles, and vice versa
      lquad = .not. lquad
      lquad = .false.    ! triangles
      nelem = 0
      do i = 1, ni-1
        do j = 1, nj-1
          if(lquad) then
c--------- quads
           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i  ,j  )
           kelem(2,n) = kij(i+1,j  )
           kelem(3,n) = kij(i+1,j+1)
           kelem(4,n) = kij(i  ,j+1)

          else
c--------- triangles
           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i  ,j  )
           kelem(2,n) = kij(i+1,j  )
           kelem(3,n) = kij(i+1,j+1)
           kelem(4,n) = 0

           nelem = nelem + 1
           n = min( nelem , nedim )
           kelem(1,n) = kij(i+1,j+1)
           kelem(2,n) = kij(i  ,j+1)
           kelem(3,n) = kij(i  ,j  )
           kelem(4,n) = 0
          endif
        enddo
      enddo
      endif

      lrcurv = .true.
      call hsmouts(lrcurv, ldrill,
     &            nelem,kelem,vars,deps,pars,parg,
     &            kdim,ldim,nedim,nddim,nmdim,
     &            idim,jdim,ni,nj,kij, nsout)
      endif

c---- vectors to text files
      vfacd = vfac * rref
      write(*,*) 'Writing to fort.11:   e1 vectors'
      write(*,*) 'Writing to fort.12:   e2 vectors'
      write(*,*) 'Writing to fort.13:   n  vectors'
      write(*,*) 'Writing to fort.14:   d  vectors'
      call hsmoutv(11,vfacd, nnode,vars, jvtot,jve1x,deps)
      call hsmoutv(12,vfacd, nnode,vars, jvtot,jve2x,deps)
      call hsmoutv(13,vfacd, nnode,vars, jvtot,jvnx, deps)
      call hsmoutv(14,vfacd, nnode,vars, ivtot,ivdx, vars)

      if(nvarg .ne. 0) then

      write(*,*) 'Writing to fort.55:   _r0_  _r_  Gamma'
      lu = 55
      do ivarg = 1, nvarg
        n = ivarg
        k1 = kelem(1,n)
        k2 = kelem(2,n)
        k3 = kelem(3,n)
        k4 = kelem(4,n)

        do k = 1, 3
          lv = lvr0x-1 + k
          iv = ivrx-1 + k
          r0j(k,1) = pars(lv,k1)
          r0j(k,2) = pars(lv,k2)
          r0j(k,3) = pars(lv,k3)
          r0j(k,4) = pars(lv,k4)
          rj(k,1) = vars(iv,k1)
          rj(k,2) = vars(iv,k2)
          rj(k,3) = vars(iv,k3)
          rj(k,4) = vars(iv,k4)

          r0c(k) = 0.25*(r0j(k,1)+r0j(k,2)+r0j(k,3)+r0j(k,4))
          rc(k)  = 0.25*( rj(k,1)+ rj(k,2)+ rj(k,3)+ rj(k,4))
        enddo

        gamma = varg(ivarg)
        write(lu,3344) r0c, rc, gamma
 3344   format(1x,9g13.5)

        if(mod(n,ni-1) .eq. 0) then
         write(lu,*) 
        endif
      enddo

      endif

c      call solsav(nnode,vars,pars)

      stop
      end ! hsmrun



      subroutine isomat(e, v, 
     &                  tsh, zrf, srfac,
     &                  a, b, d, s)
c---------------------------------------------------------------------------
c     Sets stiffness matrices A,B,D,S for an isotropic shell.
c
c   Inputs:
c     e     modulus
c     v     Poisson's ratio
c     tsh   shell thickness
c     zrf   reference surface location parameter -1 .. +1
c     srfac transverse-shear strain energy reduction factor
c            ( = 5/6 for isotropic shell)
c
c  Outputs:
c     a(.)   stiffness tensor components  A11, A22, A12, A16, A26, A66
c     b(.)   stiffness tensor components  B11, B22, B12, B16, B26, B66
c     d(.)   stiffness tensor components  D11, D22, D12, D16, D26, D66
c     s(.)   stiffness tensor components  S55, S44
c
c---------------------------------------------------------------------------
      real a(6), b(6), d(6), s(2)

      real econ(6), scon(2)

      g = e * 0.5/(1.0+v)

c---- in-plane stiffnesses
      den = 1.0 - v**2
      econ(1) = e/den        ! c11
      econ(2) = e/den        ! c22
      econ(3) = e/den * v    ! c12
      econ(4) = 0.           ! c16
      econ(5) = 0.           ! c26
      econ(6) = 2.0*g        ! c66

c---- transverse shear stiffnesses
      scon(1) = g      ! c55
      scon(2) = g      ! c44

c---- elements of in-plane stiffness matrices A,B,D for homogeneous shell
      tfac1 =  tsh
      tfac2 = -tsh**2 * zrf / 2.0
      tfac3 =  tsh**3 * (1.0 + 3.0*zrf**2) / 12.0
      do l = 1, 6
        a(l) = econ(l) * tfac1
        b(l) = econ(l) * tfac2
        d(l) = econ(l) * tfac3
      enddo

      sfac = tsh*srfac
      s(1) = scon(1) * sfac
      s(2) = scon(2) * sfac

      return
      end ! isomat



      subroutine ortmat(e1, e2, g12, v12, 
     &                  c16, c26,
     &                  g13, g23, 
     &                  tsh, zrf, srfac,
     &                  a, b, d, s)
c---------------------------------------------------------------------------
c     Sets stiffness matrices A,B,D,S for an orthotropic shell, 
c     augmented with shear/extension coupling
c
c   Inputs:
c     e1    modulus in 1 direction
c     e2    modulus in 2 direction
c     g12   shear modulus
c     v12   Poisson's ratio
c
c     c16   12-shear / 1-extension coupling modulus
c     c26   12-shear / 2-extension coupling modulus
c
c     g13   1-direction transverse-shear modulus
c     g23   2-direction transverse-shear modulus
c
c     tsh   shell thickness
c     zrf   reference surface location parameter -1 .. +1
c     srfac transverse-shear strain energy reduction factor
c            ( = 5/6 for isotropic shell)
c
c  Outputs:
c     a(.)   stiffness tensor components  A11, A22, A12, A16, A26, A66
c     b(.)   stiffness tensor components  B11, B22, B12, B16, B26, B66
c     d(.)   stiffness tensor components  D11, D22, D12, D16, D26, D66
c     s(.)   stiffness tensor components  A55, A44
c
c---------------------------------------------------------------------------
      real a(6), b(6), d(6), s(2)

      real econ(6), scon(2)

c---- in-plane stiffnesses
      den = 1.0 - v12**2 * e2/e1
      econ(1) = e1/den        ! c11
      econ(2) = e2/den        ! c22
      econ(3) = e2/den * v12  ! c12
      econ(4) = c16           ! c16
      econ(5) = c26           ! c26
      econ(6) = 2.0*g12       ! c66

c---- transverse shear stiffnesses
      scon(1) = g13    ! c55
      scon(2) = g23    ! c44

c---- elements of in-plane stiffness matrices A,B,D for homogeneous shell
      tfac1 =  tsh
      tfac2 = -tsh**2 * zrf / 2.0
      tfac3 =  tsh**3 * (1.0 + 3.0*zrf**2) / 12.0
      do l = 1, 6
        a(l) = econ(l) * tfac1
        b(l) = econ(l) * tfac2
        d(l) = econ(l) * tfac3
      enddo

      sfac = tsh*srfac
      s(1) = scon(1) * sfac
      s(2) = scon(2) * sfac

      return
      end ! ortmat




      real function dot(a,b)
c-------------------------------------
c     Returns dot product a.b
c-------------------------------------
      real a(3), b(3)

      dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      return
      end ! dot




      subroutine qpause(ians)
      character*1 ans
      write(*,*) 'Continue, Stop, Finish?   C'
      read(*,'(a)') ans
      if    (index('c ',ans) .ne. 0) then
       ians = 1
       return
      elseif(index('f',ans) .ne. 0) then
       ians = -1
       return
      else
       stop
      endif
      end ! qpause


      function randl(x1,x2)
      data im, ia, ic / 259200 , 7141 , 54773 /
      data jran / 1 /
      save jran
c
      jran = mod(jran*ia+ic,im)
      randl = x1 + (x2-x1)*float(jran)/float(im)
      return
      end


      subroutine solsav(nnode, vars, pars)
      include 'index.inc'
      real vars(ivtot,*)
      real pars(lvtot,*)

      lu = 77
      open(lu,file='sol.dat',status='new',form='unformatted')
      write(lu) nnode
      do k = 1, nnode
        write(lu) (vars(iv,k), iv=1, ivtot)
        write(lu) (pars(lv,k), lv=1, lvtot)
      enddo
      close(lu)

      return
      end


      subroutine solget(nnode, vars, pars)
      include 'index.inc'
      real vars(ivtot,*)
      real pars(lvtot,*)

      lu = 77
      open(lu,file='sol5.dat',status='old',form='unformatted')
      read(lu) nnode
      do k = 1, nnode
        read(lu) (vars(iv,k), iv=1, ivtot)
        read(lu) (pars(lv,k), lv=1, lvtot)
      enddo
      close(lu)

      return
      end

