
      subroutine hsmglr(igres, parg,
     &                  nnode, par, var,
     &                  nglv, glv,
     &                  nelem,kelem,
     &                  kdim,nedim,
     &                  gres,
     &                  nkgres, kgres, 
     &                  gres_var,
     &                  gres_glv,
     &                  aelem,
     &                  belem )
c--------------------------------------------------------------------
c     Version 0.93                                      30 Oct 2017
c
c     Returns global-equation residual and Jacobians.
c     All data arrays are referenced via pointers defined in index.inc
c
c  Inputs:
c      igres     index of global residual to be returned
c
c      parg(.)   global parameters (can be all zero)
c
c      nnode     number of nodes k = 1..nnode
c      par(.k)   node parameters
c      var(.k)   node primary variables
c
c      nglv      number of global variables
c      glv(.)    global variables
c
c      nelem      number of elements n = 1..nelem
c      kelem(.n)  node k indices of element n
c
c  Outputs:
c      gres  )   node primary variables
c      dep(.k)   node dependent variables
c
c      ddel    last max position Newton delta  (fraction of dref)
c      adel    last max normal angle Newton delta
c
c               ( converged if  ddel < dtol  and  adel < atol )
c
c--------------------------------------------------------------------
      include 'index.inc'

c--------------------------------------------------------------------
c---- passed variables and arrays
c-   (last dimensions currently set to physical sizes to enable OOB checking)
      real parg(lgtot)

      integer nnode
      real var(ivtot,kdim),
     &     par(lvtot,kdim)

      integer nglv
      real glv(*)

      integer nelem
      integer kelem(4,nedim)

      real gres
      real gres_var(ivtot,kdim)
      real gres_glv(*)

      real aelem(3), belem(3)

c--------------------------------------------------------------------
c---- local variables and arrays
      real rrj(3,4)
      real rc(3) , rc_rrj(3,3,4)
      real dr1(3), dr1_rrj(3,3,4)
      real dr2(3), dr2_rrj(3,3,4)
      real anc(3), anc_rrj(3,3,4)

      real rj(3,4)
      real ra(3), ra_rj(3,3,4)
      real rb(3), rb_rj(3,3,4)
      real aa(3), aa_rj(3,3,4), aa_rrj(3,3,4)
      real bb(3), bb_rj(3,3,4), bb_rrj(3,3,4)

      real vc1(3)
      real vc(3), vc_aa(3,3), vc_bb(3,3), vc_rj(3,3,4), vc_rrj(3,3,4)
      real vel(3), vrel(3)

      real sim(3)

      logical lbound
c--------------------------------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c--------------------------------------------------------------------

c---- frame velocity
      vel(1) = parg(lgvelx)
      vel(2) = parg(lgvely)
      vel(3) = parg(lgvelz)

c---- will set flow tangency residual on element index igres
      nres = igres

c---- nodes of this residual element
      kr1 = kelem(1,nres)
      kr2 = kelem(2,nres)
      kr3 = kelem(3,nres)
      kr4 = kelem(4,nres)

      do k = 1, 3
        iv = ivrx-1 + k
        rrj(k,1) = var(iv,kr1)
        rrj(k,2) = var(iv,kr2)
        rrj(k,3) = var(iv,kr3)
        rrj(k,4) = var(iv,kr4)
      enddo

      do k = 1, 3
        do jj = 1, 4
          do kk = 1, 3
            rc_rrj(k,kk,jj) = 0.
            dr1_rrj(k,kk,jj) = 0.
            dr2_rrj(k,kk,jj) = 0.
            anc_rrj(k,kk,jj) = 0.
          enddo
        enddo
      enddo

      do k = 1, 3
c------ cell bisectors
        dr1(k) = - rrj(k,4) + rrj(k,3)
     &           - rrj(k,1) + rrj(k,2)
        dr2(k) =   rrj(k,4) + rrj(k,3)
     &           - rrj(k,1) - rrj(k,2)

        dr1_rrj(k,k,1) = -1.0
        dr1_rrj(k,k,2) =  1.0
        dr1_rrj(k,k,3) =  1.0
        dr1_rrj(k,k,4) = -1.0

        dr2_rrj(k,k,1) = -1.0
        dr2_rrj(k,k,2) = -1.0
        dr2_rrj(k,k,3) =  1.0
        dr2_rrj(k,k,4) =  1.0

c------ control point
        rc(k) = 0.25*(rrj(k,1)+rrj(k,4))/2.0
     &        + 0.75*(rrj(k,2)+rrj(k,3))/2.0
        rc_rrj(k,k,1) = 0.25/2.0
        rc_rrj(k,k,2) = 0.75/2.0
        rc_rrj(k,k,3) = 0.75/2.0
        rc_rrj(k,k,4) = 0.25/2.0
      enddo

c---- normal vector at control point
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        anc(k) = dr1(ic)*dr2(jc) - dr1(jc)*dr2(ic)

        do jj = 1, 4
          do kk = 1, 3
            anc_rrj(k,kk,jj) = dr1_rrj(ic,kk,jj)*dr2(jc)
     &                       - dr1_rrj(jc,kk,jj)*dr2(ic)
     &                       + dr1(ic)*dr2_rrj(jc,kk,jj)
     &                       - dr1(jc)*dr2_rrj(ic,kk,jj)
c            anc_rrj(k,kk,jj) = 0.
          enddo
        enddo
      enddo

      aelem(1) = 0.25*anc(1)
      aelem(2) = 0.25*anc(2)
      aelem(3) = 0.25*anc(3)

c      anc(1) = 0.
c      anc(2) = 0.
c      anc(3) = 1.0

c------------------------------------------------------------------
c---- clear residual and Jacobians
      gres = 0.
      do k = 1, nnode
        do iv = 1, ivtot
          gres_var(iv,k) = 0.
        enddo
      enddo
      do iglv = 1, nglv
        gres_glv(iglv) = 0.
      enddo

c-------------------------------------
c---- go over all elements
      do 100 n = 1, nelem

c---- nodes of this element
      k1 = kelem(1,n)
      k2 = kelem(2,n)
      k3 = kelem(3,n)
      k4 = kelem(4,n)

      do k = 1, 3
        iv = ivrx-1 + k
        rj(k,1) = var(iv,k1)
        rj(k,2) = var(iv,k2)
        rj(k,3) = var(iv,k3)
        rj(k,4) = var(iv,k4)
      enddo

      do k = 1, 3
        vc(k) = 0.
        do jj = 1, 4
          do kk = 1, 3
            vc_rj(k,kk,jj) = 0.
            vc_rrj(k,kk,jj) = 0.
          enddo
        enddo
      enddo


      do image = -1, 1, 2
c      do image =  1, 1

      sim(1) = 1.0
      sim(2) = float(image)
      sim(3) = 1.0

      do k = 1, 3
        do jj = 1, 4
          do kk = 1, 3
            ra_rj(k,kk,jj) = 0.
            rb_rj(k,kk,jj) = 0.
            aa_rj(k,kk,jj) = 0.
            bb_rj(k,kk,jj) = 0.
            aa_rrj(k,kk,jj) = 0.
            bb_rrj(k,kk,jj) = 0.
          enddo
        enddo
      enddo

      do k = 1, 3
        ra(k) = (0.75*rj(k,1) + 0.25*rj(k,2))*sim(k)
        rb(k) = (0.75*rj(k,4) + 0.25*rj(k,3))*sim(k)

        ra_rj(k,k,1) = 0.75*sim(k)
        ra_rj(k,k,2) = 0.25*sim(k)

        rb_rj(k,k,3) = 0.25*sim(k)
        rb_rj(k,k,4) = 0.75*sim(k)

        aa(k) = ra(k) - rc(k)
        bb(k) = rb(k) - rc(k)
        do jj = 1, 4
          aa_rj(k,k,jj) = ra_rj(k,k,jj)
          bb_rj(k,k,jj) = rb_rj(k,k,jj)
          aa_rrj(k,k,jj) = -rc_rrj(k,k,jj)
          bb_rrj(k,k,jj) = -rc_rrj(k,k,jj)
        enddo
      enddo

      if(n.eq.nres .and. image.eq.1) then
        belem(1) = rb(1) - ra(1)
        belem(2) = rb(2) - ra(2)
        belem(3) = rb(3) - ra(3)
      endif

      beta = 0.
      lbound = .true.
      call vorvel(lbound, aa, bb, vc1, vc_aa, vc_bb)

      do k = 1, 3
        vc(k) = vc(k) + vc1(k)*float(image)
        do l = 1, 3
          vc_aa(k,l) = vc_aa(k,l)*float(image)
          vc_bb(k,l) = vc_bb(k,l)*float(image)
        enddo
      enddo

c      write(*,*) 
c      write(*,*) 'ra', ra
c      write(*,*) 'rb', rb
c      write(*,*) 'rc', rc
c      write(*,*) 'aa', aa
c      write(*,*) 'bb', bb
c      write(*,*) 'vc', vc
c      write(*,*) 'vc_aa', (vc_aa(3,l), l=1, 3)
c      write(*,*) 'vc_bb', (vc_bb(3,l), l=1, 3)
c      pause


      do k = 1, 3
        do jj = 1, 4
          do kk = 1, 3
            vc_rj(k,kk,jj) = vc_rj(k,kk,jj)
     &                     + vc_aa(k,1)*aa_rj(1,kk,jj)
     &                     + vc_aa(k,2)*aa_rj(2,kk,jj)
     &                     + vc_aa(k,3)*aa_rj(3,kk,jj)
     &                     + vc_bb(k,1)*bb_rj(1,kk,jj)
     &                     + vc_bb(k,2)*bb_rj(2,kk,jj)
     &                     + vc_bb(k,3)*bb_rj(3,kk,jj)
            vc_rrj(k,kk,jj) = vc_rrj(k,kk,jj)
     &                      + vc_aa(k,1)*aa_rrj(1,kk,jj)
     &                      + vc_aa(k,2)*aa_rrj(2,kk,jj)
     &                      + vc_aa(k,3)*aa_rrj(3,kk,jj)
     &                      + vc_bb(k,1)*bb_rrj(1,kk,jj)
     &                      + vc_bb(k,2)*bb_rrj(2,kk,jj)
     &                      + vc_bb(k,3)*bb_rrj(3,kk,jj)
          enddo
        enddo
      enddo

      enddo  ! next image

c---- index of global variable is element index
      iglv = n
      gamma = glv(iglv)

c---- fluid relative to control point
      vcdn =  vc(1)*anc(1) +  vc(2)*anc(2) +  vc(3)*anc(3)
      vedn = vel(1)*anc(1) + vel(2)*anc(2) + vel(3)*anc(3)

      gres = gres + gamma*vcdn - vedn
      gres_glv(iglv) = vcdn

c      if(n.eq.1) write(*,*)
c      write(*,*) 'vc', nres, n, vc, anc, gres_glv(iglv)

      do jj = 1, 4
        kj  = kelem(jj,n)
        krj = kelem(jj,nres)

        do kk = 1, 3
          iv = ivrx-1 + kk
          gres_var(iv,kj) = gres_var(iv,kj)
     &                   + gamma*vc_rj(1,kk,jj)*anc(1)
     &                   + gamma*vc_rj(2,kk,jj)*anc(2)
     &                   + gamma*vc_rj(3,kk,jj)*anc(3)
          gres_var(iv,krj) = gres_var(iv,krj)
     &                   + gamma*vc_rrj(1,kk,jj)*anc(1)
     &                   + gamma*vc_rrj(2,kk,jj)*anc(2)
     &                   + gamma*vc_rrj(3,kk,jj)*anc(3)
     &                  + (gamma*vc(1)-vel(1))*anc_rrj(1,kk,jj)
     &                  + (gamma*vc(2)-vel(2))*anc_rrj(2,kk,jj)
     &                  + (gamma*vc(3)-vel(3))*anc_rrj(3,kk,jj)
        enddo
      enddo

 100  continue
      return

      end


      subroutine vorvel(lbound, a, b, uvw, uvw_a, uvw_b )
c-------------------------------------------------------------------
c     calculates the velocity induced by a horseshoe vortex 
c     of unit strength, with no core radius.
c
c     the point where the velocity is calculated is at 0,0,0,
c     and the vortex legs trail in the 1,0,0 direction to infinity.
c
c     positive circulation is by righthand rule from a to b.
c
c  input:
c     lbound  if .true., include contribution of bound vortex leg
c     a(3)    coordinates of vertex #1 of the vortex
c     b(3)    coordinates of vertex #2 of the vortex
c
c  output: 
c     uvw(3)      induced velocity
c     uvw_a(3,3)  duvw/da  sensitivity
c     uvw_a(3,3)  duvw/db  sensitivity
c
c-------------------------------------------------------------------
      logical lbound
      real a(3), b(3), uvw(3), uvw_a(3,3), uvw_b(3,3)
c
c
      real axb(3), axb_a(3,3), axb_b(3,3)
c
      data pi4  / 12.5663706143591729538505735331180 /
c
      asq = a(1)**2 + a(2)**2 + a(3)**2
      bsq = b(1)**2 + b(2)**2 + b(3)**2
c
      amag = sqrt(asq)
      bmag = sqrt(bsq)
c
      do k = 1, 3
        uvw(k) = 0.
        do l = 1, 3
          uvw_a(k,l) = 0.
          uvw_b(k,l) = 0.
        enddo
      enddo
c
c---- contribution from the transverse bound leg
      if (lbound .and.  amag*bmag .ne. 0.0) then
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
c
        axb_a(1,1) = 0.
        axb_a(1,2) =  b(3)
        axb_a(1,3) = -b(2)
        axb_a(2,1) = -b(3)
        axb_a(2,2) = 0.
        axb_a(2,3) =  b(1)
        axb_a(3,1) =  b(2)
        axb_a(3,2) = -b(1)
        axb_a(3,3) = 0.
c
        axb_b(1,1) = 0.
        axb_b(1,2) = -a(3)
        axb_b(1,3) =  a(2)
        axb_b(2,1) =  a(3)
        axb_b(2,2) = 0.
        axb_b(2,3) = -a(1)
        axb_b(3,1) = -a(2)
        axb_b(3,2) =  a(1)
        axb_b(3,3) = 0.
c
        adb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
c
        den     =     amag*bmag + adb
        den_asq = 0.5*bmag/amag
        den_bsq = 0.5*amag/bmag
c
        t = (1.0/amag + 1.0/bmag) / den
c
        t_adb = -t/den
        t_asq = -t/den*den_asq - 0.5/(den*amag*asq)
        t_bsq = -t/den*den_bsq - 0.5/(den*bmag*bsq)
c
        do k = 1, 3
          uvw(k) = axb(k)*t
c
          do l = 1, 3
            uvw_a(k,l) = axb(k)*t_asq  *  a(l)*2.0
     &                 + axb(k)*t_adb  *  b(l)
     &                 + axb_a(k,l)*t
            uvw_b(k,l) = axb(k)*t_bsq  *  b(l)*2.0
     &                 + axb(k)*t_adb  *  a(l)
     &                 + axb_b(k,l)*t
          enddo
        enddo
      endif
c
c---- trailing leg attached to a
      if (amag .ne. 0.0) then
        axisq = a(3)**2 + a(2)**2
c
        adi = a(1)
        rsq = axisq
c
        t = - (1.0 - adi/amag) / rsq
c
        t_adi  =  1.0 / (rsq*amag)
        t_amag = -adi / (rsq*asq )
        t_rsq  =   -t /  rsq
c
        uvw(2) = uvw(2) + a(3)*t
        uvw(3) = uvw(3) - a(2)*t
c
        uvw_a(2,3) = uvw_a(2,3) + t
        uvw_a(3,2) = uvw_a(3,2) - t
c
        uvw_a(2,1) = uvw_a(2,1) + a(3)*t_adi
        uvw_a(3,1) = uvw_a(3,1) - a(2)*t_adi
c
        uvw_a(2,2) = uvw_a(2,2) + a(3)*t_rsq*a(2)*2.0
        uvw_a(2,3) = uvw_a(2,3) + a(3)*t_rsq*a(3)*2.0
        uvw_a(3,2) = uvw_a(3,2) - a(2)*t_rsq*a(2)*2.0
        uvw_a(3,3) = uvw_a(3,3) - a(2)*t_rsq*a(3)*2.0
c
        uvw_a(2,1) = uvw_a(2,1) + a(3)*t_amag*a(1)/amag
        uvw_a(2,2) = uvw_a(2,2) + a(3)*t_amag*a(2)/amag
        uvw_a(2,3) = uvw_a(2,3) + a(3)*t_amag*a(3)/amag
        uvw_a(3,1) = uvw_a(3,1) - a(2)*t_amag*a(1)/amag
        uvw_a(3,2) = uvw_a(3,2) - a(2)*t_amag*a(2)/amag
        uvw_a(3,3) = uvw_a(3,3) - a(2)*t_amag*a(3)/amag
c
      endif
c
c---- trailing leg attached to b
      if (bmag .ne. 0.0) then
        bxisq = b(3)**2 + b(2)**2
c
        bdi = b(1)
        rsq = bxisq
c
        t =   (1.0 - bdi/bmag) / rsq
c
        t_bdi  = -1.0 / (rsq*bmag)
        t_bmag =  bdi / (rsq*bsq )
        t_rsq  =   -t /  rsq
c
        uvw(2) = uvw(2) + b(3)*t
        uvw(3) = uvw(3) - b(2)*t
c
        uvw_b(2,3) = uvw_b(2,3) + t
        uvw_b(3,2) = uvw_b(3,2) - t
c
        uvw_b(2,1) = uvw_b(2,1) + b(3)*t_bdi
        uvw_b(3,1) = uvw_b(3,1) - b(2)*t_bdi
c
        uvw_b(2,2) = uvw_b(2,2) + b(3)*t_rsq*b(2)*2.0
        uvw_b(2,3) = uvw_b(2,3) + b(3)*t_rsq*b(3)*2.0
        uvw_b(3,2) = uvw_b(3,2) - b(2)*t_rsq*b(2)*2.0
        uvw_b(3,3) = uvw_b(3,3) - b(2)*t_rsq*b(3)*2.0
c
        uvw_b(2,1) = uvw_b(2,1) + b(3)*t_bmag*b(1)/bmag
        uvw_b(2,2) = uvw_b(2,2) + b(3)*t_bmag*b(2)/bmag
        uvw_b(2,3) = uvw_b(2,3) + b(3)*t_bmag*b(3)/bmag
        uvw_b(3,1) = uvw_b(3,1) - b(2)*t_bmag*b(1)/bmag
        uvw_b(3,2) = uvw_b(3,2) - b(2)*t_bmag*b(2)/bmag
        uvw_b(3,3) = uvw_b(3,3) - b(2)*t_bmag*b(3)/bmag
c
      endif
c
      do k = 1, 3
        uvw(k) = uvw(k)/pi4 
        do l = 1, 3
          uvw_a(k,l) = uvw_a(k,l)/pi4
          uvw_b(k,l) = uvw_b(k,l)/pi4
        enddo
      enddo
c
      return
      end ! vorvel
