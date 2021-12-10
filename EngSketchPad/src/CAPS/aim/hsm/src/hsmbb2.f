
      subroutine hsmbb(b, 
     &                 b1, b1_b, 
     &                 b2, b2_b )
c--------------------------------------------------------------
c     Calculates two unit vectors b1,b2 which are orthogonal
c     to the b vector and also orthognal to each other.
c
c  Inputs: 
c     b(.)    Input vector
c
c  Outputs:
c     b1(.)   Orthogonal vectors
c     b2(.)
c     b1_b(..)  d b1 / d b
c     b2_b(..)  d b2 / d b
c
c--------------------------------------------------------------

      real b(3),
     &     b1(3),
     &     b2(3)

      real b1_b(3,3), 
     &     b2_b(3,3)

      real ab(3)
      real bt(3), bt_b(3,3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c----------------------------------------------

      ab(1) = abs(b(1))
      ab(2) = abs(b(2))
      ab(3) = abs(b(3))
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      ab(1) = 0.0
c      ab(2) = 1.0
c      ab(3) = 1.0
cc     ab(1) = 1.0
cc     ab(2) = 0.0
cc     ab(3) = 1.0
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c---- set bt along smallest component of b
      do l = 1, 3
        bt_b(1,l) = 0.
        bt_b(2,l) = 0.
        bt_b(3,l) = 0.
      enddo

      if    (ab(1) .le. min( ab(2) , ab(3) ) ) then
        bt(1) = 1.0 - b(1)*b(1)
        bt(2) =     - b(2)*b(1)
        bt(3) =     - b(3)*b(1)
        bt_b(1,1) = -  2.0*b(1)
        bt_b(2,2) = -      b(1)
        bt_b(3,3) = -      b(1)
        bt_b(2,1) = -b(2)
        bt_b(3,1) = -b(3)

      elseif(ab(2) .le. min( ab(1) , ab(3) ) ) then
        bt(1) =     - b(1)*b(2)
        bt(2) = 1.0 - b(2)*b(2)
        bt(3) =     - b(3)*b(2)
        bt_b(1,1) = -      b(2)
        bt_b(2,2) = -  2.0*b(2)
        bt_b(3,3) = -      b(2)
        bt_b(1,2) = -b(1)
        bt_b(3,2) = -b(3)

      else
        bt(1) =     - b(1)*b(3)
        bt(2) =     - b(2)*b(3)
        bt(3) = 1.0 - b(3)*b(3)
        bt_b(1,1) = -      b(3)
        bt_b(2,2) = -      b(3)
        bt_b(3,3) = -  2.0*b(3)
        bt_b(1,3) = -b(1)
        bt_b(2,3) = -b(2)
      
      endif

c---- normalize bt to get b1
      bta = sqrt(bt(1)**2 + bt(2)**2 + bt(3)**2)
      do k = 1, 3
        b1(k)  =  bt(k)/bta
        b1_bta = -bt(k)/bta**2

        do l = 1, 3
          b1_b(k,l) = bt_b(k,l)/bta
     &              + b1_bta*( bt(1)*bt_b(1,l)
     &                       + bt(2)*bt_b(2,l)
     &                       + bt(3)*bt_b(3,l) ) / bta
        enddo
      enddo

c---- set bt = b x b1
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        bt(k)     = b(ic)*b1(jc)     - b(jc)*b1(ic)
        bt_b(k,1) = b(ic)*b1_b(jc,1) - b(jc)*b1_b(ic,1)
        bt_b(k,2) = b(ic)*b1_b(jc,2) - b(jc)*b1_b(ic,2)
        bt_b(k,3) = b(ic)*b1_b(jc,3) - b(jc)*b1_b(ic,3)

        bt_b(k,ic) = bt_b(k,ic) + b1(jc)
        bt_b(k,jc) = bt_b(k,jc) - b1(ic)
      enddo

c---- normalize bt to get b2
      bta = sqrt(bt(1)**2 + bt(2)**2 + bt(3)**2)
      do k = 1, 3
        b2(k)  =  bt(k)/bta
        b2_bta = -bt(k)/bta**2

        do l = 1, 3
          b2_b(k,l) = bt_b(k,l)/bta
     &              + b2_bta*( bt(1)*bt_b(1,l)
     &                       + bt(2)*bt_b(2,l)
     &                       + bt(3)*bt_b(3,l) ) / bta
        enddo
      enddo


      return
      end ! hsmbb



      subroutine hsmjtl(d1, d2,
     &                  tj, tj_d1, tj_d2, 
     &                  lj, lj_d1, lj_d2,
     &                  nj, nj_d1, nj_d2 )
c-------------------------------------------------------------------
c     Calculates two unit vectors t,l for formulating 
c     and projecting joint matching residuals.
c
c  Inputs: 
c     d1(.)   material-normal vector at node 1
c     d2(.)   material-normal vector at node 2
c
c  Outputs:
c     tj(.)   unit vector along d1-d2
c     lj(.)   unit vector along d1 x d2
c     nj(.)   unit vector along d1+d2
c
c-------------------------------------------------------------------
      real d1(3), d2(3)

      real lj(3), lj_d1(3,3), lj_d2(3,3)
      real tj(3), tj_d1(3,3), tj_d2(3,3)
      real nj(3), nj_d1(3,3), nj_d2(3,3)

c---- local work arrays
      real dxd(3), dxd_d1(3,3), dxd_d2(3,3)
      real dxda, dxda_d1(3), dxda_d2(3)

      real dmd(3), dmd_d1(3,3), dmd_d2(3,3)
      real dmda, dmda_d1(3), dmda_d2(3)

      real dpd(3), dpd_d1(3,3), dpd_d2(3,3)
      real dpda, dpda_d1(3), dpda_d2(3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /

      data ddmin / 1.0e-6 /

c---- d1 x d2
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        dxd(k) = d1(ic)*d2(jc)
     &         - d1(jc)*d2(ic)
        dxd_d1(k,ic) =  d2(jc)
        dxd_d1(k,jc) = -d2(ic)
        dxd_d1(k,k ) = 0.
        dxd_d2(k,ic) = -d1(jc)
        dxd_d2(k,jc) =  d1(ic)
        dxd_d2(k,k ) = 0.
      enddo

      dxda = sqrt(dxd(1)**2 + dxd(2)**2 + dxd(3)**2)
      dxdai = 1.0 / max( dxda , ddmin )
      do kk = 1, 3
        dxda_d1(kk) = ( dxd(1)*dxd_d1(1,kk)
     &                + dxd(2)*dxd_d1(2,kk)
     &                + dxd(3)*dxd_d1(3,kk) ) * dxdai
        dxda_d2(kk) = ( dxd(1)*dxd_d2(1,kk)
     &                + dxd(2)*dxd_d2(2,kk)
     &                + dxd(3)*dxd_d2(3,kk) ) * dxdai
      enddo

c-----------------------------------------------------------------
      if(dxda .lt. ddmin) then
c------ d1,d2 vectors are almost parallel, so joint is "flat" ...
c-       ... make tj,lj orthogonal to d2, but otherwise arbitrary
        call hsmbb(d2, 
     &             tj, tj_d2,
     &             lj, lj_d2 )
        do kk = 1, 3
          tj_d1(1,kk) = 0.
          tj_d1(2,kk) = 0.
          tj_d1(3,kk) = 0.
          lj_d1(1,kk) = 0.
          lj_d1(2,kk) = 0.
          lj_d1(3,kk) = 0.
        enddo

c------ nj = tj x lj
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          nj(k) = tj(ic)*lj(jc) - tj(jc)*lj(ic)
          do kk = 1, 3
            nj_d1(k,kk) = tj_d1(ic,kk)*lj(jc)
     &                  - tj_d1(jc,kk)*lj(ic)
     &                  + tj(ic)*lj_d1(jc,kk)
     &                  - tj(jc)*lj_d1(ic,kk)
            nj_d2(k,kk) = tj_d2(ic,kk)*lj(jc)
     &                  - tj_d2(jc,kk)*lj(ic)
     &                  + tj(ic)*lj_d2(jc,kk)
     &                  - tj(jc)*lj_d2(ic,kk)
          enddo
        enddo

c-----------------------------------------------------------------
      else
c------ lj vector is normalized d x d
        do k = 1, 3
          lj(k) = dxd(k)/dxda
          do kk = 1, 3
            lj_d1(k,kk) = (dxd_d1(k,kk) - lj(k)*dxda_d1(kk)) / dxda
            lj_d2(k,kk) = (dxd_d2(k,kk) - lj(k)*dxda_d2(kk)) / dxda
          enddo
        enddo

c------ tj vector is normalized d1-d2, or nj vector is normalized d1+d2,
c-       whichever is bigger
        do k = 1, 3
          dmd(k) = d1(k) - d2(k)
          dpd(k) = d1(k) + d2(k)
          do kk = 1, 3
            dmd_d1(k,kk) = 0.
            dmd_d2(k,kk) = 0.
            dpd_d1(k,kk) = 0.
            dpd_d2(k,kk) = 0.
          enddo
          dmd_d1(k,k) =  1.0
          dmd_d2(k,k) = -1.0
          dpd_d1(k,k) =  1.0
          dpd_d2(k,k) =  1.0
        enddo

        dmda = sqrt(dmd(1)**2 + dmd(2)**2 + dmd(3)**2)
        dpda = sqrt(dpd(1)**2 + dpd(2)**2 + dpd(3)**2)

c- - - - - - - - - - - - - - - - - 
        if(dpda .gt. dmda) then
c------- set nj parallel to dpd
         do kk = 1, 3
           dpda_d1(kk) = ( dpd(1)*dpd_d1(1,kk)
     &                   + dpd(2)*dpd_d1(2,kk)
     &                   + dpd(3)*dpd_d1(3,kk) ) / dpda
           dpda_d2(kk) = ( dpd(1)*dpd_d2(1,kk)
     &                   + dpd(2)*dpd_d2(2,kk)
     &                   + dpd(3)*dpd_d2(3,kk) ) / dpda
         enddo

         do k = 1, 3
           nj(k) = dpd(k)/dpda
           do kk = 1, 3
             nj_d1(k,kk) = (dpd_d1(k,kk) - nj(k)*dpda_d1(kk)) / dpda
             nj_d2(k,kk) = (dpd_d2(k,kk) - nj(k)*dpda_d2(kk)) / dpda
           enddo
         enddo

c- - - - - - - - - - - - - - - - - 
        else
c------- set nj parallel to dmd
         do kk = 1, 3
           dmda_d1(kk) = ( dmd(1)*dmd_d1(1,kk)
     &                   + dmd(2)*dmd_d1(2,kk)
     &                   + dmd(3)*dmd_d1(3,kk) ) / dmda
           dmda_d2(kk) = ( dmd(1)*dmd_d2(1,kk)
     &                   + dmd(2)*dmd_d2(2,kk)
     &                   + dmd(3)*dmd_d2(3,kk) ) / dmda
         enddo

         do k = 1, 3
           nj(k) = dmd(k)/dmda
           do kk = 1, 3
             nj_d1(k,kk) = (dmd_d1(k,kk) - nj(k)*dmda_d1(kk)) / dmda
             nj_d2(k,kk) = (dmd_d2(k,kk) - nj(k)*dmda_d2(kk)) / dmda
           enddo
         enddo
c- - - - - - - - - - - - - - - - - 
        endif

c------ however nj was defined, set tj = lj x nj
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          tj(k) = lj(ic)*nj(jc) - lj(jc)*nj(ic)
          do kk = 1, 3
            tj_d1(k,kk) = lj_d1(ic,kk)*nj(jc)
     &                  - lj_d1(jc,kk)*nj(ic)
     &                  + lj(ic)*nj_d1(jc,kk)
     &                  - lj(jc)*nj_d1(ic,kk)
            tj_d2(k,kk) = lj_d2(ic,kk)*nj(jc)
     &                  - lj_d2(jc,kk)*nj(ic)
     &                  + lj(ic)*nj_d2(jc,kk)
     &                  - lj(jc)*nj_d2(ic,kk)
          enddo
        enddo

c-----------------------------------------------------------------
      endif

      return
      end ! hsmjtl




      subroutine hsmdda(d1, d2,
     &                  lj, lj_d1, lj_d2,
     &                  aj, aj_d1, aj_d2 )
c-------------------------------------------------------------------
c     Calculates angle between two unit vectors d1, d2
c     Used for formulating joint matching residual.
c
c  Inputs: 
c     d1(.)   material-normal vector at node 1
c     d2(.)   material-normal vector at node 2
c     lj(.)   unit vector along d1 x d2
c
c  Outputs:
c     aj         angle between d1, d2
c     aj_d1(.)   d(aj)/d(d1)
c     aj_d2(.)   d(aj)/d(d2)
c
c-------------------------------------------------------------------

c---- passed arrays
      real d1(3), d2(3)
      real lj(3), lj_d1(3,3), lj_d2(3,3)
      real aj, aj_d1(3), aj_d2(3)

c---- local work arrays
      real dxd(3), dxd_d1(3,3), dxd_d2(3,3)

      real dcrs, dcrs_d1(3), dcrs_d2(3)
      real ddot, ddot_d1(3), ddot_d2(3)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c----------------------------------------------

c---- d1 x d2
      do kk = 1, 3
        ic = icrs(kk)
        jc = jcrs(kk)
        dxd(kk) = d1(ic)*d2(jc)
     &          - d1(jc)*d2(ic)
        dxd_d1(kk,ic) =  d2(jc)
        dxd_d1(kk,jc) = -d2(ic)
        dxd_d1(kk,kk) = 0.
        dxd_d2(kk,ic) = -d1(jc)
        dxd_d2(kk,jc) =  d1(ic)
        dxd_d2(kk,kk) = 0.
      enddo

c---- dxd . lj
      dcrs = dxd(1)*lj(1) + dxd(2)*lj(2) + dxd(3)*lj(3)
      do ll = 1, 3
        dcrs_d1(ll) = dxd_d1(1,ll)*lj(1)
     &              + dxd_d1(2,ll)*lj(2)
     &              + dxd_d1(3,ll)*lj(3)
     &              + dxd(1)      *lj_d1(1,ll)
     &              + dxd(2)      *lj_d1(2,ll)
     &              + dxd(3)      *lj_d1(3,ll)
        dcrs_d2(ll) = dxd_d2(1,ll)*lj(1)
     &              + dxd_d2(2,ll)*lj(2)
     &              + dxd_d2(3,ll)*lj(3)
     &              + dxd(1)      *lj_d2(1,ll)
     &              + dxd(2)      *lj_d2(2,ll)
     &              + dxd(3)      *lj_d2(3,ll)
      enddo

c---- d1 . d2
      ddot = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)
      do ll = 1, 3
        ddot_d1(ll) = d2(ll)
        ddot_d2(ll) = d1(ll)
      enddo

c---- angle
      aj = atan2( dcrs , ddot )

      aj_dcrs =  ddot / (dcrs**2 + ddot**2)
      aj_ddot = -dcrs / (dcrs**2 + ddot**2)
      do kk = 1, 3
        aj_d1(kk) = aj_dcrs*dcrs_d1(kk) + aj_ddot*ddot_d1(kk)
        aj_d2(kk) = aj_dcrs*dcrs_d2(kk) + aj_ddot*ddot_d2(kk)
      enddo

      return
      end ! hsmdda


