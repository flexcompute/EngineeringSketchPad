 
      subroutine hsmout(nelem,kelem,vars,deps,pars,parg,
     &                  kdim,ldim,nedim,nddim,nmdim,
     &                  idim,jdim,ni,nj,kij)
c-----------------------------------------------------
c     Writes solution to files
c-----------------------------------------------------
      include 'index.inc'

c---- node primary variables
      real vars(ivtot,kdim)

c---- node parameters
      real pars(lvtot,kdim)

c---- node secondary (dependent) variables
      real deps(jvtot,kdim)

c---- global parameters
      real parg(lgtot)

c---- element -> node  pointers
      integer kelem(4,nedim)


      integer kij(idim,jdim)


      integer kel(5)

      real n0x, n0y, n0z

      real mp1, mp2

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc---- common blocks for debugging
c      real resc(ivtot,kdim)
c      real resp(irtot,kdim)
c      real resp0(irtot,kdim)
c      real ares(kdim), dres(kdim)
c      common / rcom / resc, resp0, resp, ares, dres
cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      rtod = 45.0 / atan(1.0)

 1100 format('#', a,a,a,a,a,a,a)
 1200 format('#', 20(6x,i2,7x))
 2000 format(1x,30g15.7)

      write(*,*) 'Writing to fort.1 :',
     &       '   _r0_   _r_  _d_ _qxyz_ qn  psi'
      lu = 1
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      dx             dy             dz      ',
     & '      qx             qy             qz      ',
     & '      qn             psi '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 14)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       vars(ivdx ,k)  ,  !  7
     &       vars(ivdy ,k)  ,  !  8 
     &       vars(ivdz ,k)  ,  !  9 
     &       pars(lvqx ,k)  ,  ! 10
     &       pars(lvqy ,k)  ,  ! 11
     &       pars(lvqz ,k)  ,  ! 12
     &       pars(lvqn ,k)  ,  ! 13
     &       vars(ivps ,k)     ! 14
        enddo

        write(lu,*)
        write(lu,*)
      enddo

      write(*,*) 'Writing to fort.2 :',
     &    ' r0+/-zd   r+/-zd   (top/bot surfaces)'
      lu = 2

      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 6)
      do ibt = 1, 2

      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          tsh    = pars(lvtsh,k)
          zetref = pars(lvzrf,k)

c-------- fractional distance of bottom or top from reference surface
          if(ibt .eq. 1) then
c--------- bottom surface
           zeta = -tsh*0.5*(1.0+zetref)
          else
c--------- top surface
           zeta =  tsh*0.5*(1.0-zetref)
          endif

c-------- undeformed material normal vector (same as n0)
          n0x = pars(lvn0x,k)
          n0y = pars(lvn0y,k)
          n0z = pars(lvn0z,k)

c-------- deformed material normal vector n~, tilted from n
          dtx = vars(ivdx,k)
          dty = vars(ivdy,k)
          dtz = vars(ivdz,k)

          write(lu,2000) 
     &       pars(lvr0x,k) + zeta*n0x   ,  !  1
     &       pars(lvr0y,k) + zeta*n0y   ,  !  2 
     &       pars(lvr0z,k) + zeta*n0z   ,  !  3 
     &       vars(ivrx ,k) + zeta*dtx   ,  !  4
     &       vars(ivry ,k) + zeta*dty   ,  !  5 
     &       vars(ivrz ,k) + zeta*dtz      !  6 
        enddo
        write(lu,*)
        write(lu,*)
      enddo ! next element

      enddo ! next bot/top side


      write(*,*) 'Writing to fort.3 :',
     &           '   _r0_   _r_   eps   kap   ga_'
      lu = 3
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '     eps11          eps22          eps12     ',
     & '     kap11          kap22          kap12     ',
     & '      gam1           gam2    ',
     & '       n1             n2     '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 14)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          d1 = vars(ivdx,k)*deps(jve1x,k)
     &       + vars(ivdy,k)*deps(jve1y,k)
     &       + vars(ivdz,k)*deps(jve1z,k)
          d2 = vars(ivdx,k)*deps(jve2x,k)
     &       + vars(ivdy,k)*deps(jve2y,k)
     &       + vars(ivdz,k)*deps(jve2z,k)

          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1   x0
     &       pars(lvr0y,k)  ,  !  2   y0
     &       pars(lvr0z,k)  ,  !  3   z0
     &       vars(ivrx ,k)  ,  !  4   x
     &       vars(ivry ,k)  ,  !  5   y
     &       vars(ivrz ,k)  ,  !  6   z
     &       deps(jve11,k)  ,  !  7   e11
     &       deps(jve22,k)  ,  !  8   e22
     &       deps(jve12,k)  ,  !  9   e12
     &       deps(jvk11,k)  ,  ! 10   k11
     &       deps(jvk22,k)  ,  ! 11   k22
     &       deps(jvk12,k)  ,  ! 12   k12
     &       deps(jvga1,k)  ,  ! 13   gam1
     &       deps(jvga2,k)     ! 14   gam2
        enddo
        write(lu,*)
        write(lu,*)
      enddo


      write(*,*) 'Writing to fort.4 :',
     &           '   _r0_   _r_   _f_   _m_   fn_'
      lu = 4
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      f11            f22            f12      ',
     & '      m11            m22            m12      ',
     & '      f1n            f2n     '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 14)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       deps(jvf11,k)  ,  !  7 
     &       deps(jvf22,k)  ,  !  8 
     &       deps(jvf12,k)  ,  !  9
     &       deps(jvm11,k)  ,  ! 10 
     &       deps(jvm22,k)  ,  ! 11 
     &       deps(jvm12,k)  ,  ! 12
     &       deps(jvf1n,k)  ,  ! 13
     &       deps(jvf2n,k)     ! 14 
        enddo
        write(lu,*)
        write(lu,*)
      enddo



      write(*,*) 'Writing to fort.7 :',
     &           '   _r0_   _r_   fp1   fp2   af1'
      lu = 7
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      fp1            fp2            af1      '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 9)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          call pval(deps(jvf11,k), deps(jvf22,k), deps(jvf12,k),
     &              fp1, fp2, af1)
          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       fp1           ,  !  7 
     &       fp2           ,  !  8 
     &       rtod*af1         !  9
        enddo
        write(lu,*)
        write(lu,*)
      enddo



      write(*,*) 'Writing to fort.8 :',
     &            '   _r0_   _r_   mp1   mp2   am1'
      lu = 8
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      mp1            mp2            am1      '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 9)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)
          call pval(deps(jvm11,k), deps(jvm22,k), deps(jvm12,k),
     &              mp1, mp2, am1)
          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       mp1           ,  !  7 
     &       mp2           ,  !  8 
     &       rtod*am1         !  9
        enddo
        write(lu,*)
        write(lu,*)
      enddo


      write(*,*) 'Writing to fort.9 :',
     &           '   _r0_   r-r0'
      lu = 9
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      ux             uy             uz       '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 6)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nk = 3
        else
c------- quad element
         nk = 4
        endif

        do i = 1, nk
          kel(i) = kelem(i,n) 
        enddo
        kel(nk+1) = kel(1)

        do i = 1, nk+1
          k = kel(i)

          ux = vars(ivrx,k) - pars(lvr0x,k)
          uy = vars(ivry,k) - pars(lvr0y,k)
          uz = vars(ivrz,k) - pars(lvr0z,k)
          write(lu,2000) 
     &       pars(lvr0x,k)   ,  !  1
     &       pars(lvr0y,k)   ,  !  2 
     &       pars(lvr0z,k)   ,  !  3 
     &       ux             ,  !  4 
     &       uy             ,  !  5 
     &       uz                !  6 
        enddo
        write(lu,*)
        write(lu,*)
      enddo


      write(*,*) 'Writing to fort.31 in ij format:',
     &       '   _r0_   _r_  _d_ _qxyz_ qn  psi'
      lu = 31
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      dx             dy             dz      ',
     & '      qx             qy             qz      ',
     & '      qn             psi '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 14)
      do j = 1, nj
        do i = 1, ni
          k = kij(i,j)

          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       vars(ivdx ,k)  ,  !  7
     &       vars(ivdy ,k)  ,  !  8 
     &       vars(ivdz ,k)  ,  !  9 
     &       pars(lvqx ,k)  ,  ! 10
     &       pars(lvqy ,k)  ,  ! 11
     &       pars(lvqz ,k)  ,  ! 12
     &       pars(lvqn ,k)  ,  ! 13
     &       vars(ivps ,k)     ! 14
        enddo
        write(lu,*)
      enddo

      write(*,*) 'Writing to fort.33 in ij format:',
     &           '   _r0_   _r_   eps   kap   ga_  psi'
      lu = 33
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '     eps11          eps22          eps12     ',
     & '     kap11          kap22          kap12     ',
     & '      gam1           gam2           psi '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 15)
      do j = 1, nj
        do i = 1, ni
          k = kij(i,j)
          d1 = vars(ivdx,k)*deps(jve1x,k)
     &       + vars(ivdy,k)*deps(jve1y,k)
     &       + vars(ivdz,k)*deps(jve1z,k)
          d2 = vars(ivdx,k)*deps(jve2x,k)
     &       + vars(ivdy,k)*deps(jve2y,k)
     &       + vars(ivdz,k)*deps(jve2z,k)

          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1   x0
     &       pars(lvr0y,k)  ,  !  2   y0
     &       pars(lvr0z,k)  ,  !  3   z0
     &       vars(ivrx ,k)  ,  !  4   x
     &       vars(ivry ,k)  ,  !  5   y
     &       vars(ivrz ,k)  ,  !  6   z
     &       deps(jve11,k)  ,  !  7   e11
     &       deps(jve22,k)  ,  !  8   e22
     &       deps(jve12,k)  ,  !  9   e12
     &       deps(jvk11,k)  ,  ! 10   k11
     &       deps(jvk22,k)  ,  ! 11   k22
     &       deps(jvk12,k)  ,  ! 12   k12
     &       deps(jvga1,k)  ,  ! 13   gam1
     &       deps(jvga2,k)  ,  ! 14   gam2
     &       vars(ivps ,k)     ! 15   psi
        enddo
        write(lu,*)
      enddo

      write(*,*) 'Writing to fort.34 in ij grid:',
     &           '   _r0_   _r_   _f_   _m_   fn_'
      lu = 34
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '      f11            f22            f12      ',
     & '      m11            m22            m12      ',
     & '      f1n            f2n     '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 14)
      do j = 1, nj
        do i = 1, ni
          k = kij(i,j)
          write(lu,2000) 
     &       pars(lvr0x,k)  ,  !  1
     &       pars(lvr0y,k)  ,  !  2 
     &       pars(lvr0z,k)  ,  !  3 
     &       vars(ivrx ,k)  ,  !  4
     &       vars(ivry ,k)  ,  !  5 
     &       vars(ivrz ,k)  ,  !  6 
     &       deps(jvf11,k)  ,  !  7 
     &       deps(jvf22,k)  ,  !  8 
     &       deps(jvf12,k)  ,  !  9
     &       deps(jvm11,k)  ,  ! 10 
     &       deps(jvm22,k)  ,  ! 11 
     &       deps(jvm12,k)  ,  ! 12
     &       deps(jvf1n,k)  ,  ! 13
     &       deps(jvf2n,k)     ! 14 
        enddo
        write(lu,*)
      enddo



c---- skip residual output
      return

cc=================================================================
cc     Residual output below (mainly useful for debugging)
cc     requires the residual data to be brought in 
cc     via the common block
c
c      write(*,*) 
c     & 'Writing to fort.16:  ',
c     & ' x0  y0   Rfx  Rfy  Rfz  Rmx  Rmy  Rmz  Rth'
c      lu = 16
c      write(lu,1100)
c     & '      r0x            r0y     ',
c     & '      Rfx            Rfy            Rfz      ',
c     & '      Rmx            Rmy            Rmz      ',
c     & '      Rth    '
cc       12345678901234 12345678901234 12345678901234 
c      write(lu,1200) (icol, icol=1, 9)
c      do n = 1, nelem
c        if(kelem(4,n).eq.0) then
cc------- tri element
c         nk = 3
c        else
cc------- quad element
c         nk = 4
c        endif
c
c        do i = 1, nk
c          kel(i) = kelem(i,n) 
c        enddo
c        kel(nk+1) = kel(1)
c
c        do i = 1, nk+1
c          k = kel(i)
c
c          aki = 1.0/ares(k)
c          write(lu,2000) 
c     &       pars(lvr0x,k)    ,  !  1
c     &       pars(lvr0y,k)    ,  !  2 
c     &       resc(1,k)*aki   ,  !  3 
c     &       resc(2,k)*aki   ,  !  4 
c     &       resc(3,k)*aki   ,  !  5 
c     &       resc(4,k)*aki   ,  !  6 
c     &       resc(5,k)*aki   ,  !  7 
c     &       resc(6,k)*aki   ,  !  8 
c     &       resc(7,k)*aki      !  9 
c        enddo
c        write(lu,*)
c        write(lu,*)
c      enddo
c
c      write(*,*) 'Writing to fort.17:  ',
c     &  ' x0  y0   Rr1  Rr2  Rph  Rn1  Rn2  Rrn '
c      lu = 17
c      write(lu,1100)
c     & '      r0x            r0y     ',
c     & '      Rr1            Rr2            Rph      ',
c     & '      Rn1            Rn2            Rrn      '
cc       12345678901234 12345678901234 12345678901234 
c      write(lu,1200) (icol, icol=1, 8)
c      do n = 1, nelem
c        if(kelem(4,n).eq.0) then
cc------- tri element
c         nk = 3
c        else
cc------- quad element
c         nk = 4
c        endif
c
c        do i = 1, nk
c          kel(i) = kelem(i,n) 
c        enddo
c        kel(nk+1) = kel(1)
c
c        do i = 1, nk+1
c          k = kel(i)
c
c          aki = 1.0/ares(k)
c          write(lu,2000) 
c     &       pars(lvr0x,k)       ,  !  1
c     &       pars(lvr0y,k)       ,  !  2 
c     &       resp0(irf1,k)*aki   ,  !  3 
c     &       resp0(irf2,k)*aki   ,  !  4 
c     &       resp0(irfn,k)*aki   ,  !  5 
c     &       resp0(irm1,k)*aki   ,  !  6 
c     &       resp0(irm2,k)*aki   ,  !  7 
c     &       resp0(irmn,k)*aki      !  8 
c        enddo
c        write(lu,*)
c        write(lu,*)
c      enddo
c
c      write(*,*) 'Writing to fort.18:  ',
c     &  ' x0  y0   Rr1  Rr2  Rph  Rn1  Rn2  Rrn   (+BC)'
c      lu = 18
c      write(lu,1100)
c     & '      r0x            r0y     ',
c     & '      Rr1            Rr2            Rph      ',
c     & '      Rn1            Rn2            Rrn      '
cc       12345678901234 12345678901234 12345678901234 
c      write(lu,1200) (icol, icol=1, 8)
c      do n = 1, nelem
c        if(kelem(4,n).eq.0) then
cc------- tri element
c         nk = 3
c        else
cc------- quad element
c         nk = 4
c        endif
c
c        do i = 1, nk
c          kel(i) = kelem(i,n) 
c        enddo
c        kel(nk+1) = kel(1)
c
c        do i = 1, nk+1
c          k = kel(i)
c
c          aki = 1.0/ares(k)
c          write(lu,2000) 
c     &       pars(lvr0x,k)       ,  !  1
c     &       pars(lvr0y,k)       ,  !  2 
c     &       resp(irf1,k)*aki   ,  !  3 
c     &       resp(irf2,k)*aki   ,  !  4 
c     &       resp(irfn,k)*aki   ,  !  5 
c     &       resp(irm1,k)*aki   ,  !  6 
c     &       resp(irm2,k)*aki   ,  !  7 
c     &       resp(irmn,k)*aki      !  8 
c        enddo
c        write(lu,*)
c        write(lu,*)
c      enddo
c
c      return
      end ! hsmout




      subroutine hsmoutv(lu,vfac, nnode,vars, ivdim,ivx,vec)
      include 'index.inc'

      real vars(ivtot,*)
      real vec(ivdim,*)

 3100 format(1x,30g15.7)

      ivy = ivx + 1
      ivz = ivx + 2

      do k = 1, nnode
        x0 = vars(ivrx,k)
        y0 = vars(ivry,k)
        z0 = vars(ivrz,k)

        write(lu,3100) 
     &       x0,
     &       y0,
     &       z0
        write(lu,3100) 
     &       x0 + vfac*vec(ivx,k), 
     &       y0 + vfac*vec(ivy,k),
     &       z0 + vfac*vec(ivz,k)
        write(lu,*)
        write(lu,*)
      enddo

      return
      end ! hsmoutv



      subroutine pval(fxx,fyy,fxy,f1,f2,a1)
c--------------------------------------------------
c     Calculates principal values (eigenvalues) 
c     of the input symmetric 2x2 tensor
c
c        fxx  fxy
c        fxy  fyy
c
c   Inputs:
c       fxx, fyy, fxy
c
c   Outputs:
c       f1  max eigenvalue
c       f2  min eigenvalue
c       a1  eigenvector angle
c
c   Note: the principal x,y directions are 
c        {  cos(a1) , sin(a1) }  for f1
c        { -sin(a1) , cos(a1) }  for f2
c--------------------------------------------------

      favg = 0.5*(fxx+fyy)
      fdev = sqrt( 0.25*(fxx-fyy)**2 + fxy**2 )

      if(favg .ge. 0.0) then
       f1 = favg + fdev
       f2 = favg - fdev
       a1 = 0.5*atan2( -2.0*fxy , fxx-fyy)
      else
       f1 = favg - fdev
       f2 = favg + fdev
       a1 = 0.5*atan2( -2.0*fxy , fxx-fyy) + 2.0*atan(1.0)
      endif

      return
      end ! pval



 
      subroutine hsmouts(lrcurv, ldrill,
     &                   nelem,kelem,vars,deps,pars,parg,
     &                   kdim,ldim,nedim,nddim,nmdim,
     &                   idim,jdim,ni,nj,kij, ns)
c---------------------------------------------------------------
c     Writes solution to files, on ns x ns grid on each element
c---------------------------------------------------------------
      include 'index.inc'

c---------------------------------------------------------------
c---- passed arrays
      logical lrcurv, ldrill

c---- node primary variables
      real vars(ivtot,kdim)

c---- node parameters
      real pars(lvtot,kdim)

c---- node secondary (dependent) variables
      real deps(jvtot,kdim)

c---- global parameters
      real parg(lgtot)

c---- element -> node  pointers
      integer kelem(4,nedim)

      integer kij(idim,jdim)

c---------------------------------------------------------------
c---- local work arrays
      parameter (ngdim=9+81)

      real r0(3,ngdim)
      real d0(3,ngdim)

      real a01(3,ngdim)
      real a02(3,ngdim)

      real a01i(3,ngdim)
      real a02i(3,ngdim)

      real g0(2,2,ngdim)
      real h0(2,2,ngdim)
      real l0(2,ngdim)
      real ginv0(2,2,ngdim)
      real gdet0(ngdim)
      real rxr0(3,ngdim)

      real r(3,ngdim),
     &     r_rj(3,3,4,ngdim), r_dj(3,3,4,ngdim), r_pj(3,4,ngdim)
      real d(3,ngdim),
     &     d_rj(3,3,4,ngdim), d_dj(3,3,4,ngdim), d_pj(3,4,ngdim)

      real a1(3,ngdim),
     &     a1_rj(3,3,4,ngdim), a1_dj(3,3,4,ngdim), a1_pj(3,4,ngdim)
      real a2(3,ngdim),
     &     a2_rj(3,3,4,ngdim), a2_dj(3,3,4,ngdim), a2_pj(3,4,ngdim)

      real a1i(3,ngdim),
     &     a1i_rj(3,3,4,ngdim), a1i_dj(3,3,4,ngdim), a1i_pj(3,4,ngdim)
      real a2i(3,ngdim),
     &     a2i_rj(3,3,4,ngdim), a2i_dj(3,3,4,ngdim), a2i_pj(3,4,ngdim)

      real g(2,2,ngdim),
     &     g_rj(2,2,3,4,ngdim), g_dj(2,2,3,4,ngdim), g_pj(2,2,4,ngdim)
      real h(2,2,ngdim),
     &     h_rj(2,2,3,4,ngdim), h_dj(2,2,3,4,ngdim), h_pj(2,2,4,ngdim)

      real l(2,ngdim)  ,
     &     l_rj(2,  3,4,ngdim), l_dj(2,  3,4,ngdim), l_pj(2,  4,ngdim)

      real ginv(2,2,ngdim),
     &     ginv_rj(2,2,3,4,ngdim), ginv_dj(2,2,3,4,ngdim),
     &     ginv_pj(2,2,4,ngdim)

      real gdet(ngdim),
     &     gdet_rj(3,4,ngdim), gdet_dj(3,4,ngdim), gdet_pj(4,ngdim)

      real rxr(3,ngdim), 
     &     rxr_rj(3,3,4,ngdim)

      real ga(2,ngdim) , eab(2,2,ngdim), kab(2,2,ngdim)
      real gam(3,ngdim), eps(3,3,ngdim), kap(3,3,ngdim)

      real xsi(ngdim), eta(ngdim), fwt(ngdim)
c---------------------------------------------------------------

      rtod = 45.0 / atan(1.0)
      zero = 0.0

 1100 format('#', a,a,a,a,a,a,a)
 1200 format('#', 20(6x,i2,7x))
 2000 format(1x,30g15.7)

      write(*,*) 'Writing to fort.51 :',
     &       '   _r0_   _r_  _gxyz_  _exy_  _kxy_'
      lu = 51
      write(lu,1100)
     & '      r0x            r0y            r0z      ',
     & '      rx             ry             rz       ',
     & '     gamx           gamy           gamz      '
c       12345678901234 12345678901234 12345678901234 
      write(lu,1200) (icol, icol=1, 6)
      do n = 1, nelem
        if(kelem(4,n).eq.0) then
c------- tri element
         nj = 3
         k1 = kelem(1,n)
         k2 = kelem(2,n)
         k3 = kelem(3,n)
         k4 = kelem(3,n)

         ng = 3

         xsi(1) = 0.1666666666666666666666667
         xsi(2) = 0.6666666666666666666666667
         xsi(3) = 0.1666666666666666666666667

         eta(1) = 0.1666666666666666666666667
         eta(2) = 0.1666666666666666666666667
         eta(3) = 0.6666666666666666666666667

         fwt(1) = 0.1666666666666666666666667
         fwt(2) = 0.1666666666666666666666667
         fwt(3) = 0.1666666666666666666666667

        else
c------- quad element
         nj = 4
         k1 = kelem(1,n)
         k2 = kelem(2,n)
         k3 = kelem(3,n)
         k4 = kelem(4,n)

         ng = 4

         xsi(1) = -0.577350269189626
         xsi(2) =  0.577350269189626
         xsi(3) =  0.577350269189626
         xsi(4) = -0.577350269189626

         eta(1) = -0.577350269189626
         eta(2) = -0.577350269189626
         eta(3) =  0.577350269189626
         eta(4) =  0.577350269189626

         fwt(1) = 1.0
         fwt(2) = 1.0
         fwt(3) = 1.0
         fwt(4) = 1.0
        endif

        ig0 = ng

        do j = 0, ns
          do i = 0, ns
            fi = float(i)/float(ns)
            fj = float(j)/float(ns)

            ng = ng+1
            ig = min( ng , ngdim )
            if(nj.eq.3) then
              xsi(ig) = min( fi , 1.0-fj )
              eta(ig) = min( fj , 1.0-fi )
              fwt(ig) = 0.
             else
              xsi(ig) = -1.0 + 2.0*fi
              eta(ig) = -1.0 + 2.0*fj
              fwt(ig) = 0.
            endif
          enddo
        enddo

        if(ng .gt. ngdim) then
         write(*,*) 'HSMOUTS: Array overflow.  Increase ngdim to', ngdim
         stop
        endif

        call hsmgeo(nj, ng, xsi,eta, fwt,
     &     lrcurv, ldrill,
     &     pars(lvr0x,k1),pars(lvr0x,k2),pars(lvr0x,k3),pars(lvr0x,k4),
     &     pars(lvn0x,k1),pars(lvn0x,k2),pars(lvn0x,k3),pars(lvn0x,k4),
     &     zero          , zero         , zero         , zero         ,
     &                    r0,  r_rj,  r_dj,  r_pj,
     &                    d0,  d_rj,  d_dj,  d_pj,
     &                   a01, a1_rj, a1_dj, a1_pj,
     &                   a02, a2_rj, a2_dj, a2_pj,
     &                   a01i, a1i_rj, a1i_dj, a1i_pj,
     &                   a02i, a2i_rj, a2i_dj, a2i_pj,
     &                   g0  ,   g_rj,   g_dj,   g_pj,
     &                   h0  ,   h_rj,   h_dj,   h_pj,
     &                   l0  ,   l_rj,   l_dj,   l_pj,
     &                   ginv0, ginv_rj, ginv_dj, ginv_pj,
     &                   gdet0, gdet_rj, gdet_dj, gdet_pj,
     &                   rxr0 ,  rxr_rj )

        call hsmgeo(nj, ng, xsi,eta, fwt,
     &     lrcurv, ldrill,
     &     vars(ivrx,k1), vars(ivrx,k2), vars(ivrx,k3), vars(ivrx,k4),
     &     vars(ivdx,k1), vars(ivdx,k2), vars(ivdx,k3), vars(ivdx,k4),
     &     vars(ivps,k1), vars(ivps,k2), vars(ivps,k3), vars(ivps,k4),
     &                    r,  r_rj,  r_dj,  r_pj,
     &                    d,  d_rj,  d_dj,  d_pj,
     &                   a1, a1_rj, a1_dj, a1_pj,
     &                   a2, a2_rj, a2_dj, a2_pj,
     &                   a1i, a1i_rj, a1i_dj, a1i_pj,
     &                   a2i, a2i_rj, a2i_dj, a2i_pj,
     &                   g  ,   g_rj,   g_dj,   g_pj,
     &                   h  ,   h_rj,   h_dj,   h_pj,
     &                   l  ,   l_rj,   l_dj,   l_pj,
     &                   ginv, ginv_rj, ginv_dj, ginv_pj,
     &                   gdet, gdet_rj, gdet_dj, gdet_pj,
     &                   rxr ,  rxr_rj )

        do ig = ig0+1, ng

          do ia = 1, 2
            ga(ia,ig) = l(ia,ig) - l0(ia,ig)
            do ib = 1, 2
              eab(ia,ib,ig) = (g(ia,ib,ig) - g0(ia,ib,ig))*0.5
              kab(ia,ib,ig) = (h(ia,ib,ig) - h0(ia,ib,ig))
            enddo
          enddo
        

          do kk = 1, 3
            gam(kk,ig) = ga(1,ig)*a01i(kk,ig)
     &                 + ga(2,ig)*a02i(kk,ig)
            do ll = 1, 3
              eps(kk,ll,ig) = eab(1,1,ig)*a01i(kk,ig)*a01i(ll,ig)
     &                      + eab(1,2,ig)*a01i(kk,ig)*a02i(ll,ig)
     &                      + eab(2,1,ig)*a02i(kk,ig)*a01i(ll,ig)
     &                      + eab(2,2,ig)*a02i(kk,ig)*a02i(ll,ig)
              kap(kk,ll,ig) = kab(1,1,ig)*a01i(kk,ig)*a01i(ll,ig)
     &                      + kab(1,2,ig)*a01i(kk,ig)*a02i(ll,ig)
     &                      + kab(2,1,ig)*a02i(kk,ig)*a01i(ll,ig)
     &                      + kab(2,2,ig)*a02i(kk,ig)*a02i(ll,ig)
            enddo
          enddo

          rxra = sqrt(rxr(1,ig)**2 + rxr(2,ig)**2 + rxr(3,ig)**2)
          rxr(1,ig) = rxr(1,ig)/rxra
          rxr(2,ig) = rxr(2,ig)/rxra
          rxr(3,ig) = rxr(3,ig)/rxra

          a1a  = sqrt(a1(1,ig)**2 + a1(2,ig)**2 + a1(3,ig)**2)
          a01a = sqrt(a01(1,ig)**2 + a01(2,ig)**2 + a01(3,ig)**2)

          write(lu,2000) 
     &       r0(1,ig)  ,  !  1
     &       r0(2,ig)  ,  !  2 
     &       r0(3,ig)  ,  !  3 
     &       r(1,ig)  ,  !  4
     &       r(2,ig)  ,  !  5 
     &       r(3,ig)  ,  !  6 
     &       d(1,ig)  ,  !  7
     &       d(2,ig)  ,  !  8 
     &       d(3,ig)  ,  !  9 
     &       rxr(1,ig)  ,  !  10
     &       rxr(2,ig)  ,  !  11 
     &       rxr(3,ig)  ,  !  12
     &       gam(1,ig),  !  13
     &       gam(2,ig),  !  14
     &       gam(3,ig),  !  15
     &       eps(1,1,ig),  ! 16
     &       eps(1,2,ig),  ! 17
     &       eps(2,1,ig),  ! 18
     &       eps(2,2,ig),  ! 19
     &       kap(1,1,ig),  ! 20
     &       kap(1,2,ig),  ! 21
     &       kap(2,1,ig),  ! 22
     &       kap(2,2,ig),  ! 23
     &       gdet(ig),     ! 24  |d1r|
     &       a01a,         ! 25
     &       a1a,          ! 26
     &       0.5*(a1a**2 - a01a**2)  ! 27

          if(mod(ig-ig0,ns+1) .eq. 0) then
           write(lu,*)
          endif

        enddo ! next ig

        write(lu,*)
        write(lu,*)
      enddo ! next element

      return
      end ! hsmouts
