
       subroutine hsmgeo(neln, 
     &                   ng, xsig, etag, fwtg, 
     &                   lrcurv, ldrill,
     &                   r1, r2, r3, r4,
     &                   d1, d2, d3, d4,
     &                   p1, p2, p3, p4,
     &                    r,    r_rj,   r_dj,   r_pj,
     &                    d,    d_rj,   d_dj,   d_pj,
     &                   a1,   a1_rj,  a1_dj,  a1_pj,  
     &                   a2,   a2_rj,  a2_dj,  a2_pj,
     &                   a1i, a1i_rj, a1i_dj, a1i_pj,
     &                   a2i, a2i_rj, a2i_dj, a2i_pj,
     &                   g  ,   g_rj,   g_dj,   g_pj,
     &                   h  ,   h_rj,   h_dj,   h_pj,
     &                   l  ,   l_rj,   l_dj,   l_pj,
     &                   ginv, ginv_rj, ginv_dj, ginv_pj,
     &                   gdet, gdet_rj, gdet_dj, gdet_pj,
     &                   rxr ,  rxr_rj )
c---------------------------------------------------------------------------
c     Calculates basis vectors and metric tensors for a 3-node or 4-node
c     isoparametric/MITC element of a 2D manifold in 3D space,
c     at the xsig,etag element Gauss points
c
c     Nodes 1,2,3,4 are ordered counterclockwise around the element
c
c     All vectors are assumed to be in common xyz cartesian axes
c
c----------------------------------------------
c Inputs:
c   neln     number of element nodes (3 = tri, 4 = quad)
c
c   ng       number of Gauss points
c   xsig(.)  Gauss point xsi coordinate locations
c   etag(.)  Gauss point eta coordinate locations
c   fwtg(.)  Gauss integration weights 
c   
c   lrcurv     T  define vectors on curved surface
c              F  define vectors on bilinear surface
c   ldrill     T  include drilling DOF
c              F  exclude drilling DOF
c   r#(.)    node position vectors
c   d#(.)    node material normal vectors
c   p#       drilling rotations
c
c   Note: If neln=3, then r4,d4,p4 will not affect results
c
c----------------------------------------------
c Outputs:
c    r(.)    position vector
c    d(.)    unit director vector
c   a1(.)    covariant basis vector along xsi,  = dr/dxsi
c   a2(.)    covariant basis vector along eta,  = dr/deta
c   a1i(.)   contravariant basis vector of xsi
c   a2i(.)   contravariant basis vector of eta
c   g(..)    metric tensor     (1st fundamental form of manifold)
c   h(..)    curvature tensor  (2nd fundamental form of manifold)
c   l(.)     director lean vector  (d.a1,d.a2)
c   ginv(..) inverse metric tensor
c   gdet     metric tensor determinant
c   rxr(.)   d1r x d2r  bilinear area vector
c
c   []_rj    Jacobians wrt r1,r2,r3,r4 
c   []_dj    Jacobians wrt d1,d2,d3,d4 
c   []_pj    Jacobians wrt p1,p2,p3,p4 
c---------------------------------------------------------------------------
      implicit real (a-h,m,o-z)

c----------------------------------------------
c---- passed variables
      integer neln

      real xsig(ng), etag(ng), fwtg(ng)

      logical lrcurv, ldrill

      real r1(3),
     &     r2(3),
     &     r3(3),
     &     r4(3)
      real d1(3),
     &     d2(3),
     &     d3(3),
     &     d4(3)

      real r(3,ng), r_rj(3,3,4,ng), r_dj(3,3,4,ng), r_pj(3,4,ng)
      real d(3,ng), d_rj(3,3,4,ng), d_dj(3,3,4,ng), d_pj(3,4,ng)

      real a1(3,ng), a1_rj(3,3,4,ng), a1_dj(3,3,4,ng), a1_pj(3,4,ng)
      real a2(3,ng), a2_rj(3,3,4,ng), a2_dj(3,3,4,ng), a2_pj(3,4,ng)

      real a1i(3,ng), a1i_rj(3,3,4,ng), a1i_dj(3,3,4,ng), a1i_pj(3,4,ng)
      real a2i(3,ng), a2i_rj(3,3,4,ng), a2i_dj(3,3,4,ng), a2i_pj(3,4,ng)

      real g(2,2,ng), g_rj(2,2,3,4,ng), g_dj(2,2,3,4,ng), g_pj(2,2,4,ng)
      real h(2,2,ng), h_rj(2,2,3,4,ng), h_dj(2,2,3,4,ng), h_pj(2,2,4,ng)
      real l(2,ng)  , l_rj(2,  3,4,ng), l_dj(2,  3,4,ng), l_pj(2,  4,ng)
      real ginv(2,2,ng), ginv_rj(2,2,3,4,ng), ginv_dj(2,2,3,4,ng),
     &                   ginv_pj(2,2,  4,ng)
      real gdet(ng), gdet_rj(3,4,ng), gdet_dj(3,4,ng), gdet_pj(4,ng)
      real rxr(3,ng), rxr_rj(3,3,4,ng)

c----------------------------------------------
c---- local work arrays
c     parameter (ngdim = 4)
      parameter (ngdim = 85)

      real rj(3,4)
      real dj(3,4)
      real pj(4)

      real nnc(4)

      real nn(4)
      real d1nn(4), d2nn(4)

      real nnb(4)
      real d1nnb(4), d2nnb(4)

      real db(3), db_dj(3,3,4)
      real d1db(3), d1db_dj(3,3,4)
      real d2db(3), d2db_dj(3,3,4)

      real dba, dba_dj(3,4)
      real d1dba, d1dba_dj(3,4)
      real d2dba, d2dba_dj(3,4)

      real rb(3,ngdim), rb_rj(3,3,4,ngdim)
      real d1rb(3,ngdim), d1rb_rj(3,3,4,ngdim)
      real d2rb(3,ngdim), d2rb_rj(3,3,4,ngdim)

      real rn(3,ngdim), rn_rj(3,3,4,ngdim), rn_dj(3,3,4,ngdim)
      real d1rn(3,ngdim), d1rn_rj(3,3,4,ngdim), d1rn_dj(3,3,4,ngdim)
      real d2rn(3,ngdim), d2rn_rj(3,3,4,ngdim), d2rn_dj(3,3,4,ngdim)
 
      real rs(3,ngdim), rs_rj(3,3,4,ngdim), rs_dj(3,3,4,ngdim),
     &                  rs_pj(3,4,ngdim)
      real d1rs(3,ngdim), d1rs_rj(3,3,4,ngdim), d1rs_dj(3,3,4,ngdim),
     &                    d1rs_pj(3,4,ngdim)
      real d2rs(3,ngdim), d2rs_rj(3,3,4,ngdim), d2rs_dj(3,3,4,ngdim),
     &                    d2rs_pj(3,4,ngdim)

      real d1d(3,ngdim), d1d_rj(3,3,4,ngdim), d1d_dj(3,3,4,ngdim)
      real d2d(3,ngdim), d2d_rj(3,3,4,ngdim), d2d_dj(3,3,4,ngdim)

      real rxra_rj(3,4)

      real d1r(3,ngdim), d1r_rj(3,3,4,ngdim), d1r_dj(3,3,4,ngdim),
     &                   d1r_pj(3,4,ngdim)
      real d2r(3,ngdim), d2r_rj(3,3,4,ngdim), d2r_dj(3,3,4,ngdim),
     &                   d2r_pj(3,4,ngdim)

      real d1ra(ngdim), d1ra_rj(3,4,ngdim), d1ra_dj(3,4,ngdim),
     &                  d1ra_pj(4,ngdim)
      real d2ra(ngdim), d2ra_rj(3,4,ngdim), d2ra_dj(3,4,ngdim),
     &                  d2ra_pj(4,ngdim)

      real bnb(ngdim), bnb_rj(3,4,ngdim), bnb_dj(3,4,ngdim)
      real d1bnb(ngdim), d1bnb_rj(3,4,ngdim), d1bnb_dj(3,4,ngdim)
      real d2bnb(ngdim), d2bnb_rj(3,4,ngdim), d2bnb_dj(3,4,ngdim)

      real d1rbs(3)
      real d2rbs(3)

      real elen(4)
      real dlr(3)
      real de(3),               de_dj(3,3,4)
      real le(4), le_rj(4,3,4), le_dj(4,3,4)
      real be(4), be_rj(4,3,4), be_dj(4,3,4)


      real dxrj(3,4), dxrj_rj(3,4,3,4), dxrj_dj(3,4,3,4)
      real d1dxrj(3,4), d1dxrj_rj(3,4,3,4), d1dxrj_dj(3,4,3,4)
      real d2dxrj(3,4), d2dxrj_rj(3,4,3,4), d2dxrj_dj(3,4,3,4)

      real d1rbsq, d1rbsq_rj(3,4)
      real d2rbsq, d2rbsq_rj(3,4)

      real sb1, sb1_rj(3,4), sb1_dj(3,4), sb1_pj(4)
      real sb2, sb2_rj(3,4), sb2_dj(3,4), sb2_pj(4)

      logical lhavg

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c----------------------------------------------
c---- sqrt(2)
      data sq2 / 1.41421356237309514547462185873882845 /
c----------------------------------------------
c      include 'gausst.inc'
c      include 'gaussq.inc'

      lhavg = .true.
c     lhavg = .false.

      if(ng .gt. ngdim) then
       write(*,*) 'HSMGEO: Local work array overflow.',
     &  '  Increase  ngdim  to', ng
       stop
      endif

c---- put passed node quantities into more convenient node-indexed arrays
      do k = 1, 3
        rj(k,1) = r1(k)
        rj(k,2) = r2(k)
        rj(k,3) = r3(k)
        rj(k,4) = r4(k)

        dj(k,1) = d1(k)
        dj(k,2) = d2(k)
        dj(k,3) = d3(k)
        dj(k,4) = d4(k)
      enddo
      pj(1) = p1
      pj(2) = p2
      pj(3) = p3
      pj(4) = p4

c--------------------------------
c---- clear j=4 edge quantities to prevent undefined-variable referencing 
c-    for the tri element case
      j = 4
      le(j) = 0.
      be(j) = 0.
      do jj = 1, 4
        do kk = 1, 3
          le_rj(j,kk,jj) = 0.
          le_dj(j,kk,jj) = 0.
          be_rj(j,kk,jj) = 0.
          be_dj(j,kk,jj) = 0.
        enddo
      enddo

      do k = 1, 3
        dxrj(k,j) = 0.
        d1dxrj(k,j) = 0.
        d2dxrj(k,j) = 0.
        do jj = 1, 4
          do kk = 1, 3
            dxrj_rj(k,j,kk,jj) = 0.
            dxrj_dj(k,j,kk,jj) = 0.
            d1dxrj_rj(k,j,kk,jj) = 0.
            d1dxrj_dj(k,j,kk,jj) = 0.
            d2dxrj_rj(k,j,kk,jj) = 0.
            d2dxrj_dj(k,j,kk,jj) = 0.
          enddo
        enddo
      enddo

c--------------------------------
c---- clear j=4 Jacobians to prevent undefined-variable referencing 
c-    for the tri element case
      do ig = 1, ng
c      do jj = 1, 4  !!!
      jj = 4

      do k = 1, 3
        do kk = 1, 3
          r_rj(k,kk,jj,ig) = 0.
          r_dj(k,kk,jj,ig) = 0.
          d_rj(k,kk,jj,ig) = 0.
          d_dj(k,kk,jj,ig) = 0.
          a1_rj(k,kk,jj,ig) = 0.
          a1_dj(k,kk,jj,ig) = 0.
          a2_rj(k,kk,jj,ig) = 0.
          a2_dj(k,kk,jj,ig) = 0.
          a1i_rj(k,kk,jj,ig) = 0.
          a1i_dj(k,kk,jj,ig) = 0.
          a2i_rj(k,kk,jj,ig) = 0.
          a2i_dj(k,kk,jj,ig) = 0.
          rxr_rj(k,kk,jj,ig) = 0.
        enddo
        r_pj(k,jj,ig) = 0.
        d_pj(k,jj,ig) = 0.
        a1_pj(k,jj,ig) = 0.
        a2_pj(k,jj,ig) = 0.
        a1i_pj(k,jj,ig) = 0.
        a2i_pj(k,jj,ig) = 0.
      enddo

      do kk = 1, 3
        do ib = 1, 2
          do ia = 1, 2
            g_rj(ia,ib,kk,jj,ig) = 0.
            g_dj(ia,ib,kk,jj,ig) = 0.
            h_rj(ia,ib,kk,jj,ig) = 0.
            h_dj(ia,ib,kk,jj,ig) = 0.
            ginv_rj(ia,ib,kk,jj,ig) = 0.
            ginv_dj(ia,ib,kk,jj,ig) = 0.
          enddo
          l_rj(ib,kk,jj,ig) = 0.
          l_dj(ib,kk,jj,ig) = 0.
        enddo
        gdet_rj(kk,jj,ig) = 0.
        gdet_dj(kk,jj,ig) = 0.
      enddo

      do ib = 1, 2
        do ia = 1, 2
          g_pj(ia,ib,jj,ig) = 0.
          h_pj(ia,ib,jj,ig) = 0.
          ginv_pj(ia,ib,jj,ig) = 0.
        enddo
        l_pj(ib,jj,ig) = 0.
      enddo
      gdet_pj(jj,ig) = 0.

cc     enddo  ! next jj !!!

      enddo ! next ig

c---- set edge lengths in xi-eta plane
      if(neln.eq.3) then
       elen(1) = 1.0
       elen(2) = sq2
       elen(3) = 1.0
cc     elen(4) = 1.0

      else
       elen(1) = 2.0
       elen(2) = 2.0
       elen(3) = 2.0
       elen(4) = 2.0
      endif
c=========================================================================

c---- go over each edge
      do j = 1, neln
c------ nodes at edge endpoints
        jo = j
        jp = mod(j,neln) + 1

c------ coefficient of bubble contribution for this edge
c-      (director divergence along edge)
        be(j) = (dj(1,jp)-dj(1,jo))*(rj(1,jp)-rj(1,jo))/elen(j)
     &        + (dj(2,jp)-dj(2,jo))*(rj(2,jp)-rj(2,jo))/elen(j)
     &        + (dj(3,jp)-dj(3,jo))*(rj(3,jp)-rj(3,jo))/elen(j)
        do jj = 1, neln
          do kk = 1, 3
            be_rj(j,kk,jj) = 0.
            be_dj(j,kk,jj) = 0.
          enddo
        enddo
        do kk = 1, 3
          be_rj(j,kk,jp) =  (dj(kk,jp)-dj(kk,jo))/elen(j)
          be_rj(j,kk,jo) = -(dj(kk,jp)-dj(kk,jo))/elen(j)
          be_dj(j,kk,jp) =  (rj(kk,jp)-rj(kk,jo))/elen(j)
          be_dj(j,kk,jo) = -(rj(kk,jp)-rj(kk,jo))/elen(j)
        enddo
      enddo ! next edge

c-----------------------------------------------
c---- interpolate director leans (for transverse shears) using MITC

c---- first set director lean at midpoint of each edge
      do j = 1, neln
c------ nodes at edge endpoints
        jo = j
        jp = mod(j,neln) + 1

c------ d at edge midpoint
        desq = (dj(1,jp)+dj(1,jo))**2
     &       + (dj(2,jp)+dj(2,jo))**2
     &       + (dj(3,jp)+dj(3,jo))**2
        dea = sqrt(desq)
        do k = 1, 3
          de(k) = (dj(k,jp) + dj(k,jo)) / dea

          do jj = 1, neln
            de_dj(k,1,jj) = 0.
            de_dj(k,2,jj) = 0.
            de_dj(k,3,jj) = 0.
          enddo 

          do kk = 1, 3
            de_dj(k,kk,jp) = -de(k)/desq * (dj(kk,jp)+dj(kk,jo))
            de_dj(k,kk,jo) = -de(k)/desq * (dj(kk,jp)+dj(kk,jo))
          enddo
          de_dj(k,k,jp) = de_dj(k,k,jp) + 1.0/dea
          de_dj(k,k,jo) = de_dj(k,k,jo) + 1.0/dea
        enddo  ! next vector component k

c------ position directional derivative along edge
        dlr(1) = (rj(1,jp)-rj(1,jo))/elen(j)
        dlr(2) = (rj(2,jp)-rj(2,jo))/elen(j)
        dlr(3) = (rj(3,jp)-rj(3,jo))/elen(j)

c------ director lean along edge at midpoint
        le(j) = de(1)*dlr(1)
     &        + de(2)*dlr(2)
     &        + de(3)*dlr(3)
        do jj = 1, neln
          do kk = 1, 3
            le_rj(j,kk,jj) = 0.
            le_dj(j,kk,jj) = 
     &          de_dj(1,kk,jj)*dlr(1)
     &        + de_dj(2,kk,jj)*dlr(2)
     &        + de_dj(3,kk,jj)*dlr(3)
          enddo
        enddo
        le_rj(j,1,jp) = le_rj(j,1,jp) + de(1)/elen(j)
        le_rj(j,2,jp) = le_rj(j,2,jp) + de(2)/elen(j)
        le_rj(j,3,jp) = le_rj(j,3,jp) + de(3)/elen(j)
        le_rj(j,1,jo) = le_rj(j,1,jo) - de(1)/elen(j)
        le_rj(j,2,jo) = le_rj(j,2,jo) - de(2)/elen(j)
        le_rj(j,3,jo) = le_rj(j,3,jo) - de(3)/elen(j)
      enddo ! next edge

c====================================================================================
c---- go over Gauss points
      do 100 ig = 1, ng

      xsi = xsig(ig)
      eta = etag(ig)

c---- interpolate director lean edge values to Gauss points
      if(neln .eq. 3) then
c----- MITC3 interpolation
       l(1,ig) =  le(1)*(1.0-eta) - (le(3) + sq2*le(2))*eta
       l(2,ig) = -le(3)*(1.0-xsi) + (le(1) + sq2*le(2))*xsi
       do jj = 1, neln
         do kk = 1, 3
           l_rj(1,kk,jj,ig) =  le_rj(1,kk,jj)*(1.0-eta)
     &                       - le_rj(3,kk,jj)*eta
     &                   - sq2*le_rj(2,kk,jj)*eta
           l_rj(2,kk,jj,ig) = -le_rj(3,kk,jj)*(1.0-xsi)
     &                       + le_rj(1,kk,jj)*xsi
     &                   + sq2*le_rj(2,kk,jj)*xsi
 
           l_dj(1,kk,jj,ig) =  le_dj(1,kk,jj)*(1.0-eta)
     &                       - le_dj(3,kk,jj)*eta
     &                   - sq2*le_dj(2,kk,jj)*eta
           l_dj(2,kk,jj,ig) = -le_dj(3,kk,jj)*(1.0-xsi)
     &                       + le_dj(1,kk,jj)*xsi
     &                   + sq2*le_dj(2,kk,jj)*xsi
         enddo
         l_pj(1,jj,ig) = 0.
         l_pj(2,jj,ig) = 0.
       enddo

      else
c----- MITC4 interpolation
       nnc(1) =  0.5*(1.0-eta)
       nnc(3) = -0.5*(1.0+eta)

       nnc(2) =  0.5*(1.0+xsi)
       nnc(4) = -0.5*(1.0-xsi)

       l(1,ig) = le(1)*nnc(1)
     &         + le(3)*nnc(3)
       l(2,ig) = le(2)*nnc(2)
     &         + le(4)*nnc(4)
       do jj = 1, neln
         do kk = 1, 3
           l_rj(1,kk,jj,ig) = le_rj(1,kk,jj)*nnc(1)
     &                      + le_rj(3,kk,jj)*nnc(3)
           l_rj(2,kk,jj,ig) = le_rj(2,kk,jj)*nnc(2)
     &                      + le_rj(4,kk,jj)*nnc(4)

           l_dj(1,kk,jj,ig) = le_dj(1,kk,jj)*nnc(1)
     &                      + le_dj(3,kk,jj)*nnc(3)
           l_dj(2,kk,jj,ig) = le_dj(2,kk,jj)*nnc(2)
     &                      + le_dj(4,kk,jj)*nnc(4)
         enddo

         l_pj(1,jj,ig) = 0.
         l_pj(2,jj,ig) = 0.
       enddo

      endif

c---------------------------------------------------------------------
c---- remaining quantities interpolated from element nodes

      if(neln.eq.3) then
c----- bilinear node weights NN at gauss point
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

c----- quadratic edge weights NB at gauss point
       nnb(1) = 0.5*(1.0 - xsi - eta)*xsi
       nnb(2) = 0.5*xsi*eta * sq2
       nnb(3) = 0.5*(1.0 - xsi - eta)*eta
       nnb(4) = 0.

       d1nnb(1) =  0.5 - xsi
       d1nnb(2) =  0.5*eta * sq2
       d1nnb(3) = -0.5*eta
       d1nnb(4) = 0.
       
       d2nnb(1) = -0.5*xsi
       d2nnb(2) =  0.5*xsi * sq2
       d2nnb(3) =  0.5 - eta
       d2nnb(4) = 0.

      else
c----- bilinear node weights NN at gauss point
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

c----- quadratic edge weights NB at gauss point
       nnb(1) = 0.125*(1.0 - xsi**2)*(1.0 - eta)
       nnb(2) = 0.125*(1.0 - eta**2)*(1.0 + xsi)
       nnb(3) = 0.125*(1.0 - xsi**2)*(1.0 + eta)
       nnb(4) = 0.125*(1.0 - eta**2)*(1.0 - xsi)

       d1nnb(1) = -0.25*xsi*(1.0 - eta)
       d1nnb(2) =  0.125*(1.0 - eta**2)
       d1nnb(3) = -0.25*xsi*(1.0 + eta)
       d1nnb(4) = -0.125*(1.0 - eta**2)

       d2nnb(1) = -0.125*(1.0 - xsi**2)
       d2nnb(2) = -0.25*eta*(1.0 + xsi)
       d2nnb(3) =  0.125*(1.0 - xsi**2)
       d2nnb(4) = -0.25*eta*(1.0 - xsi)
      endif

c----------------------------------------------------------------
c---- bilinear rbar vector and derivatives drbar/dxsi, drbar/deta
      do k = 1, 3
        rb(k,ig)   = rj(k,1)*nn(1)
     &             + rj(k,2)*nn(2)
     &             + rj(k,3)*nn(3)
     &             + rj(k,4)*nn(4)
        d1rb(k,ig) = rj(k,1)*d1nn(1)
     &             + rj(k,2)*d1nn(2)
     &             + rj(k,3)*d1nn(3)
     &             + rj(k,4)*d1nn(4)
        d2rb(k,ig) = rj(k,1)*d2nn(1)
     &             + rj(k,2)*d2nn(2)
     &             + rj(k,3)*d2nn(3)
     &             + rj(k,4)*d2nn(4)
        do jj = 1, 4
          do kk = 1, 3
            rb_rj(k,kk,jj,ig) = 0.
            d1rb_rj(k,kk,jj,ig) = 0.
            d2rb_rj(k,kk,jj,ig) = 0.
          enddo
          rb_rj(k,k,jj,ig) = nn(jj)
          d1rb_rj(k,k,jj,ig) = d1nn(jj)
          d2rb_rj(k,k,jj,ig) = d2nn(jj)
        enddo
      enddo

c----------------------------------------------------------------
c---- bilinear d vector and derivatives
      do k = 1, 3
        db(k) = dj(k,1)*nn(1)
     &        + dj(k,2)*nn(2)
     &        + dj(k,3)*nn(3)
     &        + dj(k,4)*nn(4)
        d1db(k) = dj(k,1)*d1nn(1)
     &          + dj(k,2)*d1nn(2)
     &          + dj(k,3)*d1nn(3)
     &          + dj(k,4)*d1nn(4)
        d2db(k) = dj(k,1)*d2nn(1)
     &          + dj(k,2)*d2nn(2)
     &          + dj(k,3)*d2nn(3)
     &          + dj(k,4)*d2nn(4)
        do jj = 1, neln
          do kk = 1, 3
            db_dj(k,kk,jj) = 0.
            d1db_dj(k,kk,jj) = 0.
            d2db_dj(k,kk,jj) = 0.
          enddo
          db_dj(k,k,jj) = nn(jj)
          d1db_dj(k,k,jj) = d1nn(jj)
          d2db_dj(k,k,jj) = d2nn(jj)
        enddo
      enddo

c---- db magnitude
      dba = sqrt(db(1)**2 + db(2)**2 + db(3)**2)
      d1dba = ( db(1)*d1db(1)
     &       + db(2)*d1db(2)
     &       + db(3)*d1db(3) )/dba
      d2dba = ( db(1)*d2db(1)
     &       + db(2)*d2db(2)
     &       + db(3)*d2db(3) )/dba
      do jj = 1, neln
        do kk = 1, 3
          dba_dj(kk,jj) = nn(jj)*db(kk)/dba
          d1dba_dj(kk,jj) = (nn(jj)*d1db(kk) + db(kk)*d1nn(jj))/dba
     &                   - d1dba/dba*dba_dj(kk,jj)
          d2dba_dj(kk,jj) = (nn(jj)*d2db(kk) + db(kk)*d2nn(jj))/dba
     &                   - d2dba/dba*dba_dj(kk,jj)
        enddo
      enddo

c---- normalize to get unit d
      do k = 1, 3
        d(k,ig) = db(k)/dba
        d1d(k,ig) = (d1db(k) - d(k,ig)*d1dba)/dba
        d2d(k,ig) = (d2db(k) - d(k,ig)*d2dba)/dba
        do jj = 1, neln
          do kk = 1, 3
            d_rj(k,kk,jj,ig) = 0.
            d1d_rj(k,kk,jj,ig) = 0.
            d2d_rj(k,kk,jj,ig) = 0.

            d_dj(k,kk,jj,ig) = ( db_dj(k,kk,jj)
     &                         - d(k,ig)*dba_dj(kk,jj) )/dba
            d1d_dj(k,kk,jj,ig) = d1db_dj(k,kk,jj)/dba
     &                       - ( d_dj(k,kk,jj,ig)*d1dba
     &                          +d(k,ig)         *d1dba_dj(kk,jj) )/dba
     &                       - d1d(k,ig)*dba_dj(kk,jj)/dba
            d2d_dj(k,kk,jj,ig) = d2db_dj(k,kk,jj)/dba
     &                       - ( d_dj(k,kk,jj,ig)*d2dba
     &                          +d(k,ig)         *d2dba_dj(kk,jj) )/dba
     &                       - d2d(k,ig)*dba_dj(kk,jj)/dba
          enddo
          d_pj(k,jj,ig) = 0.
        enddo
      enddo

c----------------------------------------------------------------
c---- rb_xi x rb_eta
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        rxr(k,ig) = d1rb(ic,ig)*d2rb(jc,ig) - d1rb(jc,ig)*d2rb(ic,ig)
        do jj = 1, neln
          rxr_rj(k,k ,jj,ig) = 0.
          rxr_rj(k,ic,jj,ig) = 
     &           d1nn(jj)*d2rb(jc,ig) - d1rb(jc,ig)*d2nn(jj)
          rxr_rj(k,jc,jj,ig) = 
     &           d1rb(ic,ig)*d2nn(jj) - d1nn(jj)*d2rb(ic,ig)
        enddo
      enddo

c      rxra = sqrt(rxr(1)**2 + rxr(2)**2 + rxr(3)**2)
c      do jj = 1, neln
c        do kk = 1, 3
c          rxra_rj(kk,jj) = ( rxr(1,ig)*rxr_rj(1,kk,jj,ig)
c     &                     + rxr(2,ig)*rxr_rj(2,kk,jj,ig)
c     &                     + rxr(3,ig)*rxr_rj(3,kk,jj,ig) ) / rxra
c        enddo
c      enddo

c----------------------------------------------------------------
      if(lrcurv) then
c----- quadratic normal geometry offset
       bnb(ig) = be(1)*nnb(1)
     &         + be(2)*nnb(2)
     &         + be(3)*nnb(3)
     &         + be(4)*nnb(4)
       d1bnb(ig) = be(1)*d1nnb(1)
     &           + be(2)*d1nnb(2)
     &           + be(3)*d1nnb(3)
     &           + be(4)*d1nnb(4)
       d2bnb(ig) = be(1)*d2nnb(1)
     &           + be(2)*d2nnb(2)
     &           + be(3)*d2nnb(3)
     &           + be(4)*d2nnb(4)

c       if(ig .le. neln) then
c        write(*,*) ig, bnb(ig)
c       endif

       do jj = 1, neln
         do kk = 1, 3
           bnb_rj(kk,jj,ig) = be_rj(1,kk,jj)*nnb(1)
     &                      + be_rj(2,kk,jj)*nnb(2)
     &                      + be_rj(3,kk,jj)*nnb(3)
     &                      + be_rj(4,kk,jj)*nnb(4)
           bnb_dj(kk,jj,ig) = be_dj(1,kk,jj)*nnb(1)
     &                      + be_dj(2,kk,jj)*nnb(2)
     &                      + be_dj(3,kk,jj)*nnb(3)
     &                      + be_dj(4,kk,jj)*nnb(4)

           d1bnb_rj(kk,jj,ig) = be_rj(1,kk,jj)*d1nnb(1)
     &                        + be_rj(2,kk,jj)*d1nnb(2)
     &                        + be_rj(3,kk,jj)*d1nnb(3)
     &                        + be_rj(4,kk,jj)*d1nnb(4)
           d1bnb_dj(kk,jj,ig) = be_dj(1,kk,jj)*d1nnb(1)
     &                        + be_dj(2,kk,jj)*d1nnb(2)
     &                        + be_dj(3,kk,jj)*d1nnb(3)
     &                        + be_dj(4,kk,jj)*d1nnb(4)

           d2bnb_rj(kk,jj,ig) = be_rj(1,kk,jj)*d2nnb(1)
     &                        + be_rj(2,kk,jj)*d2nnb(2)
     &                        + be_rj(3,kk,jj)*d2nnb(3)
     &                        + be_rj(4,kk,jj)*d2nnb(4)
           d2bnb_dj(kk,jj,ig) = be_dj(1,kk,jj)*d2nnb(1)
     &                        + be_dj(2,kk,jj)*d2nnb(2)
     &                        + be_dj(3,kk,jj)*d2nnb(3)
     &                        + be_dj(4,kk,jj)*d2nnb(4)
         enddo
       enddo

c----- offset bubble using d vector
       do k = 1, 3
         rn(k,ig) = bnb(ig)*d(k,ig)
         d1rn(k,ig) = d1bnb(ig)*d(k,ig) + bnb(ig)*d1d(k,ig)
         d2rn(k,ig) = d2bnb(ig)*d(k,ig) + bnb(ig)*d2d(k,ig)
         do jj = 1, neln
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = bnb_rj(kk,jj,ig)*d(k,ig)
     &                         + bnb(ig)         *d_rj(k,kk,jj,ig)
             rn_dj(k,kk,jj,ig) = bnb_dj(kk,jj,ig)*d(k,ig)
     &                         + bnb(ig)         *d_dj(k,kk,jj,ig)

             d1rn_rj(k,kk,jj,ig) = d1bnb_rj(kk,jj,ig)*d(k,ig)
     &                           + bnb_rj(kk,jj,ig)*d1d(k,ig)
     &                           + bnb(ig)         *d1d_rj(k,kk,jj,ig)
             d1rn_dj(k,kk,jj,ig) = d1bnb_dj(kk,jj,ig)*d(k,ig)
     &                           + d1bnb(ig)         *d_dj(k,kk,jj,ig)
     &                           + bnb_dj(kk,jj,ig)*d1d(k,ig)
     &                           + bnb(ig)         *d1d_dj(k,kk,jj,ig)

             d2rn_rj(k,kk,jj,ig) = d2bnb_rj(kk,jj,ig)*d(k,ig)
     &                           + bnb_rj(kk,jj,ig)*d2d(k,ig)
     &                           + bnb(ig)         *d2d_rj(k,kk,jj,ig)
             d2rn_dj(k,kk,jj,ig) = d2bnb_dj(kk,jj,ig)*d(k,ig)
     &                           + d2bnb(ig)         *d_dj(k,kk,jj,ig)
     &                           + bnb_dj(kk,jj,ig)*d2d(k,ig)
     &                           + bnb(ig)         *d2d_dj(k,kk,jj,ig)
           enddo
         enddo
       enddo

       if(.false.) then
       do k = 1, 3
         rn(k,ig) = bnb(ig)*db(k)
         d1rn(k,ig) = d1bnb(ig)*db(k) + bnb(ig)*d1db(k)
         d2rn(k,ig) = d2bnb(ig)*db(k) + bnb(ig)*d2db(k)
         do jj = 1, neln
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = bnb_rj(kk,jj,ig)*db(k)
             rn_dj(k,kk,jj,ig) = bnb_dj(kk,jj,ig)*db(k)

             d1rn_rj(k,kk,jj,ig) = d1bnb_rj(kk,jj,ig)*db(k)
     &                           + bnb_rj(kk,jj,ig)*d1db(k)
             d1rn_dj(k,kk,jj,ig) = d1bnb_dj(kk,jj,ig)*db(k)
     &                           + d1bnb(ig)         *db_dj(k,kk,jj)
     &                           + bnb_dj(kk,jj,ig)*d1db(k)
     &                           + bnb(ig)         *d1db_dj(k,kk,jj)

             d2rn_rj(k,kk,jj,ig) = d2bnb_rj(kk,jj,ig)*db(k)
     &                           + bnb_rj(kk,jj,ig)*d2db(k)
             d2rn_dj(k,kk,jj,ig) = d2bnb_dj(kk,jj,ig)*db(k)
     &                           + d2bnb(ig)         *db_dj(k,kk,jj)
     &                           + bnb_dj(kk,jj,ig)*d2db(k)
     &                           + bnb(ig)         *d2db_dj(k,kk,jj)
           enddo
         enddo
       enddo
       endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
c----- bilinear: no bubble offset
       bnb(ig) = 0.
       d1bnb(ig) = 0.
       d2bnb(ig) = 0.
       do jj = 1, neln
         do kk = 1, 3
           bnb_rj(kk,jj,ig) = 0.
           bnb_dj(kk,jj,ig) = 0.
           d1bnb_rj(kk,jj,ig) = 0.
           d1bnb_dj(kk,jj,ig) = 0.
           d2bnb_rj(kk,jj,ig) = 0.
           d2bnb_dj(kk,jj,ig) = 0.
         enddo
       enddo

       do k = 1, 3
         rn(k,ig) = 0.
         d1rn(k,ig) = 0.
         d2rn(k,ig) = 0.

         do jj = 1, neln
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = 0.
             rn_dj(k,kk,jj,ig) = 0.

             d1rn_rj(k,kk,jj,ig) = 0.
             d1rn_dj(k,kk,jj,ig) = 0.

             d2rn_rj(k,kk,jj,ig) = 0.
             d2rn_dj(k,kk,jj,ig) = 0.
           enddo
         enddo
       enddo

      endif ! lrcurv
c----------------------------------------------------------------

      if(ldrill) then
c----- quadratic in-surface drilling geometry offset
       do j = 1, neln
         do k = 1, 3
           ic = icrs(k)
           jc = jcrs(k)
           dxrj(k,j)   = dj(ic,j)*(rb(jc,ig)-rj(jc,j))
     &                 - dj(jc,j)*(rb(ic,ig)-rj(ic,j))
           d1dxrj(k,j) = dj(ic,j)*d1rb(jc,ig)
     &                 - dj(jc,j)*d1rb(ic,ig)
           d2dxrj(k,j) = dj(ic,j)*d2rb(jc,ig)
     &                 - dj(jc,j)*d2rb(ic,ig)
           do jj = 1, neln
             do kk = 1, 3
               dxrj_rj(k,j,kk,jj) = dj(ic,j)*rb_rj(jc,kk,jj,ig)
     &                            - dj(jc,j)*rb_rj(ic,kk,jj,ig)
               d1dxrj_rj(k,j,kk,jj) = dj(ic,j)*d1rb_rj(jc,kk,jj,ig)
     &                              - dj(jc,j)*d1rb_rj(ic,kk,jj,ig)
               d2dxrj_rj(k,j,kk,jj) = dj(ic,j)*d2rb_rj(jc,kk,jj,ig)
     &                              - dj(jc,j)*d2rb_rj(ic,kk,jj,ig)

               dxrj_dj(k,j,kk,jj) = 0.
               d1dxrj_dj(k,j,kk,jj) = 0.
               d2dxrj_dj(k,j,kk,jj) = 0.
             enddo
           enddo
           dxrj_rj(k,j,ic,j) = dxrj_rj(k,j,ic,j) + dj(jc,j)
           dxrj_rj(k,j,jc,j) = dxrj_rj(k,j,jc,j) - dj(ic,j)
           dxrj_dj(k,j,ic,j) = dxrj_dj(k,j,ic,j) + (rb(jc,ig)-rj(jc,j))
           dxrj_dj(k,j,jc,j) = dxrj_dj(k,j,jc,j) - (rb(ic,ig)-rj(ic,j))

           d1dxrj_dj(k,j,ic,j) = d1dxrj_dj(k,j,ic,j) + d1rb(jc,ig)
           d1dxrj_dj(k,j,jc,j) = d1dxrj_dj(k,j,jc,j) - d1rb(ic,ig)

           d2dxrj_dj(k,j,ic,j) = d2dxrj_dj(k,j,ic,j) + d2rb(jc,ig)
           d2dxrj_dj(k,j,jc,j) = d2dxrj_dj(k,j,jc,j) - d2rb(ic,ig)
         enddo
       enddo

       do k = 1, 3
         rs(k,ig) = pj(1)*dxrj(k,1)*nn(1)
     &            + pj(2)*dxrj(k,2)*nn(2)
     &            + pj(3)*dxrj(k,3)*nn(3)
     &            + pj(4)*dxrj(k,4)*nn(4)
         d1rs(k,ig) = pj(1)*d1dxrj(k,1)*nn(1)
     &              + pj(2)*d1dxrj(k,2)*nn(2)
     &              + pj(3)*d1dxrj(k,3)*nn(3)
     &              + pj(4)*d1dxrj(k,4)*nn(4)
     &              + pj(1)*dxrj(k,1)*d1nn(1)
     &              + pj(2)*dxrj(k,2)*d1nn(2)
     &              + pj(3)*dxrj(k,3)*d1nn(3)
     &              + pj(4)*dxrj(k,4)*d1nn(4)
         d2rs(k,ig) = pj(1)*d2dxrj(k,1)*nn(1)
     &              + pj(2)*d2dxrj(k,2)*nn(2)
     &              + pj(3)*d2dxrj(k,3)*nn(3)
     &              + pj(4)*d2dxrj(k,4)*nn(4)
     &              + pj(1)*dxrj(k,1)*d2nn(1)
     &              + pj(2)*dxrj(k,2)*d2nn(2)
     &              + pj(3)*dxrj(k,3)*d2nn(3)
     &              + pj(4)*dxrj(k,4)*d2nn(4)
         do jj = 1, neln
           do kk = 1, 3
             rs_rj(k,kk,jj,ig) =
     &               pj(1)*dxrj_rj(k,1,kk,jj)*nn(1)
     &             + pj(2)*dxrj_rj(k,2,kk,jj)*nn(2)
     &             + pj(3)*dxrj_rj(k,3,kk,jj)*nn(3)
     &             + pj(4)*dxrj_rj(k,4,kk,jj)*nn(4)
             rs_dj(k,kk,jj,ig) =
     &               pj(1)*dxrj_dj(k,1,kk,jj)*nn(1)
     &             + pj(2)*dxrj_dj(k,2,kk,jj)*nn(2)
     &             + pj(3)*dxrj_dj(k,3,kk,jj)*nn(3)
     &             + pj(4)*dxrj_dj(k,4,kk,jj)*nn(4)
 
             d1rs_rj(k,kk,jj,ig) = 
     &                pj(1)*d1dxrj_rj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d1dxrj_rj(k,2,kk,jj)*nn(2)
     &              + pj(3)*d1dxrj_rj(k,3,kk,jj)*nn(3)
     &              + pj(4)*d1dxrj_rj(k,4,kk,jj)*nn(4)
     &              + pj(1)*dxrj_rj(k,1,kk,jj)*d1nn(1)
     &              + pj(2)*dxrj_rj(k,2,kk,jj)*d1nn(2)
     &              + pj(3)*dxrj_rj(k,3,kk,jj)*d1nn(3)
     &              + pj(4)*dxrj_rj(k,4,kk,jj)*d1nn(4)
             d2rs_rj(k,kk,jj,ig) = 
     &                pj(1)*d2dxrj_rj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d2dxrj_rj(k,2,kk,jj)*nn(2)
     &              + pj(3)*d2dxrj_rj(k,3,kk,jj)*nn(3)
     &              + pj(4)*d2dxrj_rj(k,4,kk,jj)*nn(4)
     &              + pj(1)*dxrj_rj(k,1,kk,jj)*d2nn(1)
     &              + pj(2)*dxrj_rj(k,2,kk,jj)*d2nn(2)
     &              + pj(3)*dxrj_rj(k,3,kk,jj)*d2nn(3)
     &              + pj(4)*dxrj_rj(k,4,kk,jj)*d2nn(4)

             d1rs_dj(k,kk,jj,ig) = 
     &                pj(1)*d1dxrj_dj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d1dxrj_dj(k,2,kk,jj)*nn(2)
     &              + pj(3)*d1dxrj_dj(k,3,kk,jj)*nn(3)
     &              + pj(4)*d1dxrj_dj(k,4,kk,jj)*nn(4)
     &              + pj(1)*dxrj_dj(k,1,kk,jj)*d1nn(1)
     &              + pj(2)*dxrj_dj(k,2,kk,jj)*d1nn(2)
     &              + pj(3)*dxrj_dj(k,3,kk,jj)*d1nn(3)
     &              + pj(4)*dxrj_dj(k,4,kk,jj)*d1nn(4)
             d2rs_dj(k,kk,jj,ig) = 
     &                pj(1)*d2dxrj_dj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d2dxrj_dj(k,2,kk,jj)*nn(2)
     &              + pj(3)*d2dxrj_dj(k,3,kk,jj)*nn(3)
     &              + pj(4)*d2dxrj_dj(k,4,kk,jj)*nn(4)
     &              + pj(1)*dxrj_dj(k,1,kk,jj)*d2nn(1)
     &              + pj(2)*dxrj_dj(k,2,kk,jj)*d2nn(2)
     &              + pj(3)*dxrj_dj(k,3,kk,jj)*d2nn(3)
     &              + pj(4)*dxrj_dj(k,4,kk,jj)*d2nn(4)
           enddo
           rs_pj(k,jj,ig) = dxrj(k,jj)*nn(jj)
           d1rs_pj(k,jj,ig) = d1dxrj(k,jj)*nn(jj)
     &                      + dxrj(k,jj)*d1nn(jj)
           d2rs_pj(k,jj,ig) = d2dxrj(k,jj)*nn(jj)
     &                      + dxrj(k,jj)*d2nn(jj)
         enddo
       enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
c----- no drill offset
       do k = 1, 3
         rs(k,ig) = 0.
         d1rs(k,ig) = 0.
         d2rs(k,ig) = 0.
         do jj = 1, neln
           do kk = 1, 3
             rs_rj(k,kk,jj,ig) = 0.
             rs_dj(k,kk,jj,ig) = 0.
             d1rs_rj(k,kk,jj,ig) = 0.
             d1rs_dj(k,kk,jj,ig) = 0.
             d2rs_rj(k,kk,jj,ig) = 0.
             d2rs_dj(k,kk,jj,ig) = 0.
           enddo
           rs_pj(k,jj,ig) = 0.
           d1rs_pj(k,jj,ig) = 0.
           d2rs_pj(k,jj,ig) = 0.
         enddo
       enddo

      endif ! ldrill
c--------------------------------------------------------------------

c---- geometry and derivatives
      do k = 1, 3
        r(k,ig) = rb(k,ig) + rn(k,ig) + rs(k,ig)
        d1r(k,ig) = d1rb(k,ig) + d1rn(k,ig) + d1rs(k,ig)
        d2r(k,ig) = d2rb(k,ig) + d2rn(k,ig) + d2rs(k,ig)

        do jj = 1, neln
          do kk = 1, 3
            r_rj(k,kk,jj,ig) = rb_rj(k,kk,jj,ig)
     &                       + rs_rj(k,kk,jj,ig)
     &                       + rn_rj(k,kk,jj,ig)
            d1r_rj(k,kk,jj,ig) = d1rb_rj(k,kk,jj,ig)
     &                         + d1rn_rj(k,kk,jj,ig)
     &                         + d1rs_rj(k,kk,jj,ig)
            d2r_rj(k,kk,jj,ig) = d2rb_rj(k,kk,jj,ig)
     &                         + d2rn_rj(k,kk,jj,ig)
     &                         + d2rs_rj(k,kk,jj,ig)

            r_dj(k,kk,jj,ig) = rs_dj(k,kk,jj,ig)
     &                       + rn_dj(k,kk,jj,ig)
            d1r_dj(k,kk,jj,ig) = d1rn_dj(k,kk,jj,ig)
     &                         + d1rs_dj(k,kk,jj,ig)
            d2r_dj(k,kk,jj,ig) = d2rn_dj(k,kk,jj,ig)
     &                         + d2rs_dj(k,kk,jj,ig)
          enddo
          r_pj(k,jj,ig) = rs_pj(k,jj,ig)
          d1r_pj(k,jj,ig) = d1rs_pj(k,jj,ig)
          d2r_pj(k,jj,ig) = d2rs_pj(k,jj,ig)
        enddo
      enddo

c---- d1r,d2r magnitudes
      d1ra(ig) = sqrt(d1r(1,ig)**2 + d1r(2,ig)**2 + d1r(3,ig)**2)
      d2ra(ig) = sqrt(d2r(1,ig)**2 + d2r(2,ig)**2 + d2r(3,ig)**2)
      do jj = 1, neln
        do kk = 1, 3
          d1ra_rj(kk,jj,ig) = ( d1r(1,ig)*d1r_rj(1,kk,jj,ig)
     &                        + d1r(2,ig)*d1r_rj(2,kk,jj,ig)
     &                        + d1r(3,ig)*d1r_rj(3,kk,jj,ig) )/d1ra(ig)
          d2ra_rj(kk,jj,ig) = ( d2r(1,ig)*d2r_rj(1,kk,jj,ig)
     &                        + d2r(2,ig)*d2r_rj(2,kk,jj,ig)
     &                        + d2r(3,ig)*d2r_rj(3,kk,jj,ig) )/d2ra(ig)
          d1ra_dj(kk,jj,ig) = ( d1r(1,ig)*d1r_dj(1,kk,jj,ig)
     &                        + d1r(2,ig)*d1r_dj(2,kk,jj,ig)
     &                        + d1r(3,ig)*d1r_dj(3,kk,jj,ig) )/d1ra(ig)
          d2ra_dj(kk,jj,ig) = ( d2r(1,ig)*d2r_dj(1,kk,jj,ig)
     &                        + d2r(2,ig)*d2r_dj(2,kk,jj,ig)
     &                        + d2r(3,ig)*d2r_dj(3,kk,jj,ig) )/d2ra(ig)
        enddo
        d1ra_pj(jj,ig) = ( d1r(1,ig)*d1r_pj(1,jj,ig)
     &                   + d1r(2,ig)*d1r_pj(2,jj,ig)
     &                   + d1r(3,ig)*d1r_pj(3,jj,ig) )/d1ra(ig)
        d2ra_pj(jj,ig) = ( d2r(1,ig)*d2r_pj(1,jj,ig)
     &                   + d2r(2,ig)*d2r_pj(2,jj,ig)
     &                   + d2r(3,ig)*d2r_pj(3,jj,ig) )/d2ra(ig)
      enddo


 100  continue ! next Gauss point

c==========================================================================

      asum = 0.
      do ig = 1, ng
        asum = asum + fwtg(ig)
      enddo


c---- integrate to get bubble basis vector magnitudes sb1,sb2
      sb1 = 0.
      sb2 = 0.
      do jj = 1, neln
        do kk = 1, 3
          sb1_rj(kk,jj) = 0.
          sb2_rj(kk,jj) = 0.
          sb1_dj(kk,jj) = 0.
          sb2_dj(kk,jj) = 0.
        enddo
        sb1_pj(jj) = 0.
        sb2_pj(jj) = 0.
      enddo

      do ig = 1, ng
c------ Gauss weight wk/A 
        fwt = fwtg(ig) / asum

        d1rbsq = d1rb(1,ig)**2
     &         + d1rb(2,ig)**2
     &         + d1rb(3,ig)**2
        d2rbsq = d2rb(1,ig)**2
     &         + d2rb(2,ig)**2
     &         + d2rb(3,ig)**2
        do jj = 1, neln
          do kk = 1, 3
            d1rbsq_rj(kk,jj) = 2.0*d1rb(1,ig)*d1rb_rj(1,kk,jj,ig)
     &                       + 2.0*d1rb(2,ig)*d1rb_rj(2,kk,jj,ig)
     &                       + 2.0*d1rb(3,ig)*d1rb_rj(3,kk,jj,ig)
            d2rbsq_rj(kk,jj) = 2.0*d2rb(1,ig)*d2rb_rj(1,kk,jj,ig)
     &                       + 2.0*d2rb(2,ig)*d2rb_rj(2,kk,jj,ig)
     &                       + 2.0*d2rb(3,ig)*d2rb_rj(3,kk,jj,ig)
          enddo
        enddo



c        fb1 = sqrt(d1rbsq + d1bnb(ig)**2)
c        fb2 = sqrt(d2rbsq + d2bnb(ig)**2)
c        fb1_d1rbsq = 0.5/fb1
c        fb2_d2rbsq = 0.5/fb2
c        fb1_d1bnb = d1bnb(ig)/fb1
c        fb2_d2bnb = d2bnb(ig)/fb2
c        sb1 = sb1 + fb1*fwt 
c        sb2 = sb2 + fb2*fwt 
c        do jj = 1, neln
c          do kk = 1, 3
c            sb1_rj(kk,jj) = 
c     &      sb1_rj(kk,jj) + (fb1_d1rbsq*d1rbsq_rj(kk,jj)
c     &                    +  fb1_d1bnb*d1bnb_rj(kk,jj,ig) ) * fwt
c            sb2_rj(kk,jj) = 
c     &      sb2_rj(kk,jj) + (fb2_d2rbsq*d2rbsq_rj(kk,jj)
c     &                    +  fb2_d2bnb*d2bnb_rj(kk,jj,ig) ) * fwt
c            sb1_dj(kk,jj) = 
c     &      sb1_dj(kk,jj) +  fb1_d1bnb*d1bnb_dj(kk,jj,ig)   * fwt
c            sb2_dj(kk,jj) = 
c     &      sb2_dj(kk,jj) +  fb2_d2bnb*d2bnb_dj(kk,jj,ig)   * fwt
c          enddo
c        enddo

        sb1 = sb1 + d1ra(ig)*fwt 
        sb2 = sb2 + d2ra(ig)*fwt 
        do jj = 1, neln
          do kk = 1, 3
            sb1_rj(kk,jj) = 
     &      sb1_rj(kk,jj) + d1ra_rj(kk,jj,ig)*fwt
            sb2_rj(kk,jj) = 
     &      sb2_rj(kk,jj) + d2ra_rj(kk,jj,ig)*fwt
            sb1_dj(kk,jj) = 
     &      sb1_dj(kk,jj) + d1ra_dj(kk,jj,ig)*fwt
            sb2_dj(kk,jj) = 
     &      sb2_dj(kk,jj) + d2ra_dj(kk,jj,ig)*fwt
          enddo
          sb1_pj(jj) = 
     &    sb1_pj(jj) + d1ra_pj(jj,ig)*fwt
          sb2_pj(jj) =
     &    sb2_pj(jj) + d2ra_pj(jj,ig)*fwt
        enddo
        
      enddo

c==================================================================================
c---- set final quantities at Gauss points
      do 200 ig = 1, ng

c--------------------------------------------------------------------
c----- basis vectors
       do k = 1, 3
         a1(k,ig) = d1r(k,ig)
         a2(k,ig) = d2r(k,ig)
         do jj = 1, neln
           do kk = 1, 3
             a1_rj(k,kk,jj,ig) = d1r_rj(k,kk,jj,ig)
             a2_rj(k,kk,jj,ig) = d2r_rj(k,kk,jj,ig)
             a1_dj(k,kk,jj,ig) = d1r_dj(k,kk,jj,ig)
             a2_dj(k,kk,jj,ig) = d2r_dj(k,kk,jj,ig)
           enddo
           a1_pj(k,jj,ig) = d1r_pj(k,jj,ig)
           a2_pj(k,jj,ig) = d2r_pj(k,jj,ig)
         enddo
       enddo

c--------------------------------------------------------------------
c---- covariant metric tensor g_ab
      g(1,1,ig) = a1(1,ig)*a1(1,ig)
     &          + a1(2,ig)*a1(2,ig)
     &          + a1(3,ig)*a1(3,ig)
      g(2,2,ig) = a2(1,ig)*a2(1,ig)
     &          + a2(2,ig)*a2(2,ig)
     &          + a2(3,ig)*a2(3,ig)
      g(1,2,ig) = a1(1,ig)*a2(1,ig)
     &          + a1(2,ig)*a2(2,ig)
     &          + a1(3,ig)*a2(3,ig)
      g(2,1,ig) = g(1,2,ig)

      do jj = 1, neln
        do kk = 1, 3
          g_rj(1,1,kk,jj,ig) = 2.0*a1(1,ig)*a1_rj(1,kk,jj,ig)
     &                       + 2.0*a1(2,ig)*a1_rj(2,kk,jj,ig)
     &                       + 2.0*a1(3,ig)*a1_rj(3,kk,jj,ig)
          g_rj(2,2,kk,jj,ig) = 2.0*a2(1,ig)*a2_rj(1,kk,jj,ig)
     &                       + 2.0*a2(2,ig)*a2_rj(2,kk,jj,ig)
     &                       + 2.0*a2(3,ig)*a2_rj(3,kk,jj,ig)
          g_rj(1,2,kk,jj,ig) = a1_rj(1,kk,jj,ig)*a2(1,ig)
     &                       + a1_rj(2,kk,jj,ig)*a2(2,ig)
     &                       + a1_rj(3,kk,jj,ig)*a2(3,ig)
     &                       + a1(1,ig)*a2_rj(1,kk,jj,ig)
     &                       + a1(2,ig)*a2_rj(2,kk,jj,ig)
     &                       + a1(3,ig)*a2_rj(3,kk,jj,ig)
          g_rj(2,1,kk,jj,ig) = g_rj(1,2,kk,jj,ig)

          g_dj(1,1,kk,jj,ig) = 2.0*a1(1,ig)*a1_dj(1,kk,jj,ig)
     &                       + 2.0*a1(2,ig)*a1_dj(2,kk,jj,ig)
     &                       + 2.0*a1(3,ig)*a1_dj(3,kk,jj,ig)
          g_dj(2,2,kk,jj,ig) = 2.0*a2(1,ig)*a2_dj(1,kk,jj,ig)
     &                       + 2.0*a2(2,ig)*a2_dj(2,kk,jj,ig)
     &                       + 2.0*a2(3,ig)*a2_dj(3,kk,jj,ig)
          g_dj(1,2,kk,jj,ig) = a1_dj(1,kk,jj,ig)*a2(1,ig)
     &                       + a1_dj(2,kk,jj,ig)*a2(2,ig)
     &                       + a1_dj(3,kk,jj,ig)*a2(3,ig)
     &                       + a1(1,ig)*a2_dj(1,kk,jj,ig)
     &                       + a1(2,ig)*a2_dj(2,kk,jj,ig)
     &                       + a1(3,ig)*a2_dj(3,kk,jj,ig)
          g_dj(2,1,kk,jj,ig) = g_dj(1,2,kk,jj,ig)
        enddo

        g_pj(1,1,jj,ig) = 2.0*a1(1,ig)*a1_pj(1,jj,ig)
     &                  + 2.0*a1(2,ig)*a1_pj(2,jj,ig)
     &                  + 2.0*a1(3,ig)*a1_pj(3,jj,ig)
        g_pj(2,2,jj,ig) = 2.0*a2(1,ig)*a2_pj(1,jj,ig)
     &                  + 2.0*a2(2,ig)*a2_pj(2,jj,ig)
     &                  + 2.0*a2(3,ig)*a2_pj(3,jj,ig)
        g_pj(1,2,jj,ig) = a1_pj(1,jj,ig)*a2(1,ig)
     &                  + a1_pj(2,jj,ig)*a2(2,ig)
     &                  + a1_pj(3,jj,ig)*a2(3,ig)
     &                  + a1(1,ig)*a2_pj(1,jj,ig)
     &                  + a1(2,ig)*a2_pj(2,jj,ig)
     &                  + a1(3,ig)*a2_pj(3,jj,ig)
        g_pj(2,1,jj,ig) = g_pj(1,2,jj,ig)
      enddo ! next jj

      gdet(ig)  =  g(1,1,ig)*g(2,2,ig) - g(1,2,ig)*g(2,1,ig)
      do jj = 1, neln
        do kk = 1, 3
          gdet_rj(kk,jj,ig) = 
     &       g_rj(1,1,kk,jj,ig)*g(2,2,ig) + g(1,1,ig)*g_rj(2,2,kk,jj,ig)
     &     - g_rj(1,2,kk,jj,ig)*g(2,1,ig) - g(1,2,ig)*g_rj(2,1,kk,jj,ig)
          gdet_dj(kk,jj,ig) = 
     &       g_dj(1,1,kk,jj,ig)*g(2,2,ig) + g(1,1,ig)*g_dj(2,2,kk,jj,ig)
     &     - g_dj(1,2,kk,jj,ig)*g(2,1,ig) - g(1,2,ig)*g_dj(2,1,kk,jj,ig)
        enddo
        gdet_pj(jj,ig) = 
     &       g_pj(1,1,jj,ig)*g(2,2,ig) + g(1,1,ig)*g_pj(2,2,jj,ig)
     &     - g_pj(1,2,jj,ig)*g(2,1,ig) - g(1,2,ig)*g_pj(2,1,jj,ig)
      enddo

c----------------------------------------------------------------
c---- contravariant metric tensor  g^ij, and Jacobians wrt r
      ginv(1,1,ig) =  g(2,2,ig) / gdet(ig)
      ginv(2,2,ig) =  g(1,1,ig) / gdet(ig)
      ginv(1,2,ig) = -g(1,2,ig) / gdet(ig)
      ginv(2,1,ig) = -g(2,1,ig) / gdet(ig)
      do jj = 1, neln
        do kk = 1, 3
          ginv_rj(1,1,kk,jj,ig) = ( g_rj(2,2,kk,jj,ig)
     &                     - ginv(1,1,ig)*gdet_rj(kk,jj,ig)) / gdet(ig)
          ginv_rj(2,2,kk,jj,ig) = ( g_rj(1,1,kk,jj,ig)
     &                     - ginv(2,2,ig)*gdet_rj(kk,jj,ig)) / gdet(ig)
          ginv_rj(1,2,kk,jj,ig) = (-g_rj(1,2,kk,jj,ig)
     &                     - ginv(1,2,ig)*gdet_rj(kk,jj,ig)) / gdet(ig)
          ginv_rj(2,1,kk,jj,ig) = (-g_rj(2,1,kk,jj,ig)
     &                     - ginv(2,1,ig)*gdet_rj(kk,jj,ig)) / gdet(ig)

          ginv_dj(1,1,kk,jj,ig) = ( g_dj(2,2,kk,jj,ig)
     &                     - ginv(1,1,ig)*gdet_dj(kk,jj,ig)) / gdet(ig)
          ginv_dj(2,2,kk,jj,ig) = ( g_dj(1,1,kk,jj,ig)
     &                     - ginv(2,2,ig)*gdet_dj(kk,jj,ig)) / gdet(ig)
          ginv_dj(1,2,kk,jj,ig) = (-g_dj(1,2,kk,jj,ig)
     &                     - ginv(1,2,ig)*gdet_dj(kk,jj,ig)) / gdet(ig)
          ginv_dj(2,1,kk,jj,ig) = (-g_dj(2,1,kk,jj,ig)
     &                     - ginv(2,1,ig)*gdet_dj(kk,jj,ig)) / gdet(ig)
        enddo
        ginv_pj(1,1,jj,ig) = ( g_pj(2,2,jj,ig)
     &                       - ginv(1,1,ig)*gdet_pj(jj,ig)) / gdet(ig)
        ginv_pj(2,2,jj,ig) = ( g_pj(1,1,jj,ig)
     &                       - ginv(2,2,ig)*gdet_pj(jj,ig)) / gdet(ig)
        ginv_pj(1,2,jj,ig) = (-g_pj(1,2,jj,ig)
     &                       - ginv(1,2,ig)*gdet_pj(jj,ig)) / gdet(ig)
        ginv_pj(2,1,jj,ig) = (-g_pj(2,1,jj,ig)
     &                       - ginv(2,1,ig)*gdet_pj(jj,ig)) / gdet(ig)
      enddo

c----------------------------------------------------------------
c---- contravariant basis vectors
      do k = 1, 3
        a1i(k,ig) = ginv(1,1,ig)*a1(k,ig) + ginv(1,2,ig)*a2(k,ig)
        a2i(k,ig) = ginv(2,1,ig)*a1(k,ig) + ginv(2,2,ig)*a2(k,ig)
        do jj = 1, neln
          do kk = 1, 3
            a1i_rj(k,kk,jj,ig) = ginv_rj(1,1,kk,jj,ig)*a1(k,ig)
     &                         + ginv_rj(1,2,kk,jj,ig)*a2(k,ig)
     &                         + ginv(1,1,ig)         *a1_rj(k,kk,jj,ig)
     &                         + ginv(1,2,ig)         *a2_rj(k,kk,jj,ig)
            a2i_rj(k,kk,jj,ig) = ginv_rj(2,1,kk,jj,ig)*a1(k,ig)
     &                         + ginv_rj(2,2,kk,jj,ig)*a2(k,ig)
     &                         + ginv(2,1,ig)         *a1_rj(k,kk,jj,ig)
     &                         + ginv(2,2,ig)         *a2_rj(k,kk,jj,ig)

            a1i_dj(k,kk,jj,ig) = ginv_dj(1,1,kk,jj,ig)*a1(k,ig)
     &                         + ginv_dj(1,2,kk,jj,ig)*a2(k,ig)
     &                         + ginv(1,1,ig)         *a1_dj(k,kk,jj,ig)
     &                         + ginv(1,2,ig)         *a2_dj(k,kk,jj,ig)
            a2i_dj(k,kk,jj,ig) = ginv_dj(2,1,kk,jj,ig)*a1(k,ig)
     &                         + ginv_dj(2,2,kk,jj,ig)*a2(k,ig)
     &                         + ginv(2,1,ig)         *a1_dj(k,kk,jj,ig)
     &                         + ginv(2,2,ig)         *a2_dj(k,kk,jj,ig)
          enddo

          a1i_pj(k,jj,ig) = ginv_pj(1,1,jj,ig)*a1(k,ig)
     &                    + ginv_pj(1,2,jj,ig)*a2(k,ig)
     &                    + ginv(1,1,ig)*a1_pj(k,jj,ig)
     &                    + ginv(1,2,ig)*a2_pj(k,jj,ig)
          a2i_pj(k,jj,ig) = ginv_pj(2,1,jj,ig)*a1(k,ig)
     &                    + ginv_pj(2,2,jj,ig)*a2(k,ig)
     &                    + ginv(2,1,ig)*a1_pj(k,jj,ig)
     &                    + ginv(2,2,ig)*a2_pj(k,jj,ig)
        enddo
      enddo

c===================================================================================
c---- set metric and curvature tensors for strains

      do k = 1, 3
        d1rbs(k) = d1rb(k,ig) + d1rs(k,ig)
        d2rbs(k) = d2rb(k,ig) + d2rs(k,ig)
      enddo


c---- covariant metric tensor g_ab
      g(1,1,ig) = d1rbs(1)*d1rbs(1)
     &          + d1rbs(2)*d1rbs(2)
     &          + d1rbs(3)*d1rbs(3)
      g(2,2,ig) = d2rbs(1)*d2rbs(1)
     &          + d2rbs(2)*d2rbs(2)
     &          + d2rbs(3)*d2rbs(3)
      g(1,2,ig) = d1rbs(1)*d2rbs(1)
     &          + d1rbs(2)*d2rbs(2)
     &          + d1rbs(3)*d2rbs(3)
      g(2,1,ig) = g(1,2,ig)

      do jj = 1, neln
        do kk = 1, 3
          g_rj(1,1,kk,jj,ig) =
     &         2.0*d1rbs(1)*(d1rb_rj(1,kk,jj,ig)+d1rs_rj(1,kk,jj,ig))
     &       + 2.0*d1rbs(2)*(d1rb_rj(2,kk,jj,ig)+d1rs_rj(2,kk,jj,ig))
     &       + 2.0*d1rbs(3)*(d1rb_rj(3,kk,jj,ig)+d1rs_rj(3,kk,jj,ig))
          g_rj(2,2,kk,jj,ig) =
     &         2.0*d2rbs(1)*(d2rb_rj(1,kk,jj,ig)+d2rs_rj(1,kk,jj,ig))
     &       + 2.0*d2rbs(2)*(d2rb_rj(2,kk,jj,ig)+d2rs_rj(2,kk,jj,ig))
     &       + 2.0*d2rbs(3)*(d2rb_rj(3,kk,jj,ig)+d2rs_rj(3,kk,jj,ig))
          g_rj(1,2,kk,jj,ig) = 
     &             (d1rb_rj(1,kk,jj,ig)+d1rs_rj(1,kk,jj,ig))*d2rbs(1)
     &           + (d1rb_rj(2,kk,jj,ig)+d1rs_rj(2,kk,jj,ig))*d2rbs(2)
     &           + (d1rb_rj(3,kk,jj,ig)+d1rs_rj(3,kk,jj,ig))*d2rbs(3)
     &            + d1rbs(1)*(d2rb_rj(1,kk,jj,ig)+d2rs_rj(1,kk,jj,ig))
     &            + d1rbs(2)*(d2rb_rj(2,kk,jj,ig)+d2rs_rj(2,kk,jj,ig))
     &            + d1rbs(3)*(d2rb_rj(3,kk,jj,ig)+d2rs_rj(3,kk,jj,ig))
          g_rj(2,1,kk,jj,ig) = g_rj(1,2,kk,jj,ig)

          g_dj(1,1,kk,jj,ig) =
     &         2.0*d1rbs(1)*d1rs_dj(1,kk,jj,ig)
     &       + 2.0*d1rbs(2)*d1rs_dj(2,kk,jj,ig)
     &       + 2.0*d1rbs(3)*d1rs_dj(3,kk,jj,ig)
          g_dj(2,2,kk,jj,ig) =
     &         2.0*d2rbs(1)*d2rs_dj(1,kk,jj,ig)
     &       + 2.0*d2rbs(2)*d2rs_dj(2,kk,jj,ig)
     &       + 2.0*d2rbs(3)*d2rs_dj(3,kk,jj,ig)
          g_dj(1,2,kk,jj,ig) =
     &             d1rs_dj(1,kk,jj,ig)*d2rbs(1)
     &           + d1rs_dj(2,kk,jj,ig)*d2rbs(2)
     &           + d1rs_dj(3,kk,jj,ig)*d2rbs(3)
     &           + d1rbs(1)*d2rs_dj(1,kk,jj,ig)
     &           + d1rbs(2)*d2rs_dj(2,kk,jj,ig)
     &           + d1rbs(3)*d2rs_dj(3,kk,jj,ig)
          g_dj(2,1,kk,jj,ig) = g_dj(1,2,kk,jj,ig)
        enddo
        g_pj(1,1,jj,ig) = 
     &         2.0*d1rbs(1)*d1rs_pj(1,jj,ig)
     &       + 2.0*d1rbs(2)*d1rs_pj(2,jj,ig)
     &       + 2.0*d1rbs(3)*d1rs_pj(3,jj,ig)
        g_pj(2,2,jj,ig) = 
     &         2.0*d2rbs(1)*d2rs_pj(1,jj,ig)
     &       + 2.0*d2rbs(2)*d2rs_pj(2,jj,ig)
     &       + 2.0*d2rbs(3)*d2rs_pj(3,jj,ig)
        g_pj(1,2,jj,ig) = 
     &             d1rs_pj(1,jj,ig)*d2rbs(1)
     &           + d1rs_pj(2,jj,ig)*d2rbs(2)
     &           + d1rs_pj(3,jj,ig)*d2rbs(3)
     &           + d1rbs(1)*d2rs_pj(1,jj,ig)
     &           + d1rbs(2)*d2rs_pj(2,jj,ig)
     &           + d1rbs(3)*d2rs_pj(3,jj,ig)
        g_pj(2,1,jj,ig) = g_pj(1,2,jj,ig)
      enddo ! next jj


c----------------------------------------------------------------
c---- curvature tensor h_ab
      h(1,1,ig) = d1rbs(1)*d1d(1,ig)
     &          + d1rbs(2)*d1d(2,ig)
     &          + d1rbs(3)*d1d(3,ig)
      h(2,2,ig) = d2rbs(1)*d2d(1,ig)
     &          + d2rbs(2)*d2d(2,ig)
     &          + d2rbs(3)*d2d(3,ig)

      if(lhavg) then
      h(1,2,ig) = 
     &  0.5*(  d1rbs(1)*d2d(1,ig)
     &       + d1rbs(2)*d2d(2,ig)
     &       + d1rbs(3)*d2d(3,ig)
     &       + d2rbs(1)*d1d(1,ig)
     &       + d2rbs(2)*d1d(2,ig)
     &       + d2rbs(3)*d1d(3,ig) )
      h(2,1,ig) = h(1,2,ig)
      else
      h(1,2,ig) = d1rbs(1)*d2d(1,ig)
     &          + d1rbs(2)*d2d(2,ig)
     &          + d1rbs(3)*d2d(3,ig)
      h(2,1,ig) = d2rbs(1)*d1d(1,ig)
     &          + d2rbs(2)*d1d(2,ig)
     &          + d2rbs(3)*d1d(3,ig)
      endif

      do jj = 1, neln
        do kk = 1, 3
          h_rj(1,1,kk,jj,ig) =
     &             (d1rb_rj(1,kk,jj,ig)+d1rs_rj(1,kk,jj,ig))*d1d(1,ig)
     &           + (d1rb_rj(2,kk,jj,ig)+d1rs_rj(2,kk,jj,ig))*d1d(2,ig)
     &           + (d1rb_rj(3,kk,jj,ig)+d1rs_rj(3,kk,jj,ig))*d1d(3,ig)
     &            + d1rbs(1)*d1d_rj(1,kk,jj,ig)
     &            + d1rbs(2)*d1d_rj(2,kk,jj,ig)
     &            + d1rbs(3)*d1d_rj(3,kk,jj,ig)
          h_rj(2,2,kk,jj,ig) =
     &             (d2rb_rj(1,kk,jj,ig)+d2rs_rj(1,kk,jj,ig))*d2d(1,ig)
     &           + (d2rb_rj(2,kk,jj,ig)+d2rs_rj(2,kk,jj,ig))*d2d(2,ig)
     &           + (d2rb_rj(3,kk,jj,ig)+d2rs_rj(3,kk,jj,ig))*d2d(3,ig)
     &            + d2rbs(1)*d2d_rj(1,kk,jj,ig)
     &            + d2rbs(2)*d2d_rj(2,kk,jj,ig)
     &            + d2rbs(3)*d2d_rj(3,kk,jj,ig)

          if(lhavg) then
          h_rj(1,2,kk,jj,ig) = 
     &      0.5*(  (d1rb_rj(1,kk,jj,ig)+d1rs_rj(1,kk,jj,ig))*d2d(1,ig)
     &           + (d1rb_rj(2,kk,jj,ig)+d1rs_rj(2,kk,jj,ig))*d2d(2,ig)
     &           + (d1rb_rj(3,kk,jj,ig)+d1rs_rj(3,kk,jj,ig))*d2d(3,ig)
     &            + d1rbs(1)*d2d_rj(1,kk,jj,ig)
     &            + d1rbs(2)*d2d_rj(2,kk,jj,ig)
     &            + d1rbs(3)*d2d_rj(3,kk,jj,ig)
     &           + (d2rb_rj(1,kk,jj,ig)+d2rs_rj(1,kk,jj,ig))*d1d(1,ig)
     &           + (d2rb_rj(2,kk,jj,ig)+d2rs_rj(2,kk,jj,ig))*d1d(2,ig)
     &           + (d2rb_rj(3,kk,jj,ig)+d2rs_rj(3,kk,jj,ig))*d1d(3,ig)
     &            + d2rbs(1)*d1d_rj(1,kk,jj,ig)
     &            + d2rbs(2)*d1d_rj(2,kk,jj,ig)
     &            + d2rbs(3)*d1d_rj(3,kk,jj,ig) )
          h_rj(2,1,kk,jj,ig) = h_rj(1,2,kk,jj,ig)
          else
          h_rj(1,2,kk,jj,ig) =
     &             (d1rb_rj(1,kk,jj,ig)+d1rs_rj(1,kk,jj,ig))*d2d(1,ig)
     &           + (d1rb_rj(2,kk,jj,ig)+d1rs_rj(2,kk,jj,ig))*d2d(2,ig)
     &           + (d1rb_rj(3,kk,jj,ig)+d1rs_rj(3,kk,jj,ig))*d2d(3,ig)
     &            + d1rbs(1)*d2d_rj(1,kk,jj,ig)
     &            + d1rbs(2)*d2d_rj(2,kk,jj,ig)
     &            + d1rbs(3)*d2d_rj(3,kk,jj,ig)
          h_rj(2,1,kk,jj,ig) =
     &             (d2rb_rj(1,kk,jj,ig)+d2rs_rj(1,kk,jj,ig))*d1d(1,ig)
     &           + (d2rb_rj(2,kk,jj,ig)+d2rs_rj(2,kk,jj,ig))*d1d(2,ig)
     &           + (d2rb_rj(3,kk,jj,ig)+d2rs_rj(3,kk,jj,ig))*d1d(3,ig)
     &            + d2rbs(1)*d1d_rj(1,kk,jj,ig)
     &            + d2rbs(2)*d1d_rj(2,kk,jj,ig)
     &            + d2rbs(3)*d1d_rj(3,kk,jj,ig)
          endif

          h_dj(1,1,kk,jj,ig) = 
     &              d1rs_dj(1,kk,jj,ig)*d1d(1,ig)
     &            + d1rs_dj(2,kk,jj,ig)*d1d(2,ig)
     &            + d1rs_dj(3,kk,jj,ig)*d1d(3,ig)
     &            + d1rbs(1)*d1d_dj(1,kk,jj,ig)
     &            + d1rbs(2)*d1d_dj(2,kk,jj,ig)
     &            + d1rbs(3)*d1d_dj(3,kk,jj,ig)
          h_dj(2,2,kk,jj,ig) =
     &              d2rs_dj(1,kk,jj,ig)*d2d(1,ig)
     &            + d2rs_dj(2,kk,jj,ig)*d2d(2,ig)
     &            + d2rs_dj(3,kk,jj,ig)*d2d(3,ig)
     &            + d2rbs(1)*d2d_dj(1,kk,jj,ig)
     &            + d2rbs(2)*d2d_dj(2,kk,jj,ig)
     &            + d2rbs(3)*d2d_dj(3,kk,jj,ig)
          if(lhavg) then
          h_dj(1,2,kk,jj,ig) = 
     &         0.5*(  d1rs_dj(1,kk,jj,ig)*d2d(1,ig)
     &              + d1rs_dj(2,kk,jj,ig)*d2d(2,ig)
     &              + d1rs_dj(3,kk,jj,ig)*d2d(3,ig)
     &              + d1rbs(1)*d2d_dj(1,kk,jj,ig)
     &              + d1rbs(2)*d2d_dj(2,kk,jj,ig)
     &              + d1rbs(3)*d2d_dj(3,kk,jj,ig)
     &              + d2rs_dj(1,kk,jj,ig)*d1d(1,ig)
     &              + d2rs_dj(2,kk,jj,ig)*d1d(2,ig)
     &              + d2rs_dj(3,kk,jj,ig)*d1d(3,ig)
     &              + d2rbs(1)*d1d_dj(1,kk,jj,ig)
     &              + d2rbs(2)*d1d_dj(2,kk,jj,ig)
     &              + d2rbs(3)*d1d_dj(3,kk,jj,ig) )
          h_dj(2,1,kk,jj,ig) = h_dj(1,2,kk,jj,ig)
          else
          h_dj(1,2,kk,jj,ig) =
     &              d1rs_dj(1,kk,jj,ig)*d2d(1,ig)
     &            + d1rs_dj(2,kk,jj,ig)*d2d(2,ig)
     &            + d1rs_dj(3,kk,jj,ig)*d2d(3,ig)
     &            + d1rbs(1)*d2d_dj(1,kk,jj,ig)
     &            + d1rbs(2)*d2d_dj(2,kk,jj,ig)
     &            + d1rbs(3)*d2d_dj(3,kk,jj,ig)
          h_dj(2,1,kk,jj,ig) =
     &              d2rs_dj(1,kk,jj,ig)*d1d(1,ig)
     &            + d2rs_dj(2,kk,jj,ig)*d1d(2,ig)
     &            + d2rs_dj(3,kk,jj,ig)*d1d(3,ig)
     &            + d2rbs(1)*d1d_dj(1,kk,jj,ig)
     &            + d2rbs(2)*d1d_dj(2,kk,jj,ig)
     &            + d2rbs(3)*d1d_dj(3,kk,jj,ig)
          endif
        enddo
        h_pj(1,1,jj,ig) =
     &             d1rs_pj(1,jj,ig)*d1d(1,ig)
     &           + d1rs_pj(2,jj,ig)*d1d(2,ig)
     &           + d1rs_pj(3,jj,ig)*d1d(3,ig)

        h_pj(2,2,jj,ig) =
     &             d2rs_pj(1,jj,ig)*d2d(1,ig)
     &           + d2rs_pj(2,jj,ig)*d2d(2,ig)
     &           + d2rs_pj(3,jj,ig)*d2d(3,ig)
        if(lhavg) then
        h_pj(1,2,jj,ig) =
     &         0.5*(  d1rs_pj(1,jj,ig)*d2d(1,ig)
     &              + d1rs_pj(2,jj,ig)*d2d(2,ig)
     &              + d1rs_pj(3,jj,ig)*d2d(3,ig)
     &              + d2rs_pj(1,jj,ig)*d1d(1,ig)
     &              + d2rs_pj(2,jj,ig)*d1d(2,ig)
     &              + d2rs_pj(3,jj,ig)*d1d(3,ig) )
        h_pj(2,1,jj,ig) = h_pj(1,2,jj,ig)
        else
        h_pj(1,2,jj,ig) =
     &           d1rs_pj(1,jj,ig)*d2d(1,ig)
     &         + d1rs_pj(2,jj,ig)*d2d(2,ig)
     &         + d1rs_pj(3,jj,ig)*d2d(3,ig)
        h_pj(2,1,jj,ig) =
     &           d2rs_pj(1,jj,ig)*d1d(1,ig)
     &         + d2rs_pj(2,jj,ig)*d1d(2,ig)
     &         + d2rs_pj(3,jj,ig)*d1d(3,ig)
        endif
      enddo


c@@@@@@@@@@@@@@@@@
c      gdet(ig) = d1ra(ig)
c@@@@@@@@@@@@@@@@@

 200  continue ! next ig

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c      ig = 1
c      gdet(ig) = sb1
c      do jj = 1, neln
c        do kk = 1, 3
c          gdet_rj(kk,jj,ig) = sb1_rj(kk,jj)
c          gdet_dj(kk,jj,ig) = sb1_dj(kk,jj)
c        enddo
c        gdet_pj(jj,ig) = sb1_pj(jj)
c      enddo
c      ig = 2
c      gdet(ig) = sb2
c      do jj = 1, neln
c        do kk = 1, 3
c          gdet_rj(kk,jj,ig) = sb2_rj(kk,jj)
c          gdet_dj(kk,jj,ig) = sb2_dj(kk,jj)
c        enddo
c        gdet_pj(jj,ig) = sb2_pj(jj)
c      enddo
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      return
      end ! hsmgeo



      subroutine hsmgeo2(ng, xsig, fwtg,
     &                   lrcurv, ldrill,
     &                   r1, r2,
     &                   d1, d2,
     &                   p1, p2,
     &                    r,  r_rj,  r_dj,  r_pj, 
     &                    d,  d_rj,  d_dj,  d_pj, 
     &                   a1, a1_rj, a1_dj, a1_pj, 
     &                    g,  g_rj,  g_dj,  g_pj, 
     &                    h,  h_rj,  h_dj,  h_pj, 
     &                    l,  l_rj,  l_dj,  l_pj, 
     &                   th, th_rj, th_dj, th_pj, 
     &                   lh, lh_rj, lh_dj, lh_pj )
c---------------------------------------------------------------------------
c     Calculates basis vectors and metrics for a line segment
c     in 3D space between r1 and r2, at the xsi coordinate location.
c
c     All vectors are assumed to be in common xyz cartesian axes
c
c----------------------------------------------
c Inputs:
c   
c   ng        number of Gauss points
c   xsig(.)   Gauss point locations  (-1..+1)
c
c   lrcurv    T  define vectors on curved surface
c             F  define vectors on bilinear surface
c   ldrill    T  include drilling DOF
c             F  exclude drilling DOF
c
c   r#(.)    node position vectors
c   d#(.)    node unit director vectors
c   p#       node drilling rotation angles
c
c----------------------------------------------
c Outputs:
c    r(.)    position vector
c    d(.)    director unit vector
c   a1(.)    covariant basis vector along xsi  = dr/dxsi
c    g       metric   = a1.a1
c    h       curvature
c    l       director lean
c    th(.)   edge-normal unit vector
c    lh(.)   edge-tangential unit vector
c
c   []_rj    Jacobians wrt r1,r2
c   []_dj    Jacobians wrt d1,d2
c   []_pj    Jacobians wrt p1,p2
c---------------------------------------------------------------------------
      implicit real (a-h,m,o-z)

c----------------------------------------------
c---- passed variables
      real xsig(ng), fwtg(ng)

      logical lrcurv, ldrill

      real r1(3), r2(3)
      real d1(3), d2(3)
      real p1, p2

      real r(3,ng), r_rj(3,3,2,ng), r_dj(3,3,2,ng), r_pj(3,2,ng)
      real d(3,ng), d_rj(3,3,2,ng), d_dj(3,3,2,ng), d_pj(3,2,ng)

      real a1(3,ng), a1_rj(3,3,2,ng), a1_dj(3,3,2,ng), a1_pj(3,2,ng)

      real g(ng), g_rj(3,2,ng), g_dj(3,2,ng), g_pj(2,ng)
      real h(ng), h_rj(3,2,ng), h_dj(3,2,ng), h_pj(2,ng)
      real l(ng), l_rj(3,2,ng), l_dj(3,2,ng), l_pj(2,ng)

      real th(3,ng), th_rj(3,3,2,ng), th_dj(3,3,2,ng), th_pj(3,2,ng)
      real lh(3,ng), lh_rj(3,3,2,ng), lh_dj(3,3,2,ng), lh_pj(3,2,ng)

c----------------------------------------------
c---- local work arrays
      parameter (ngdim = 3)
      real rj(3,2)
      real dj(3,2)
      real pj(2)

      real nn(2)
      real d1nn(2)

      real nnb
      real d1nnb

      real rb(3,ngdim), rb_rj(3,3,2,ngdim)
      real d1rb(3,ngdim), d1rb_rj(3,3,2,ngdim)

      real rn(3,ngdim), rn_rj(3,3,2,ngdim), rn_dj(3,3,2,ngdim)
      real d1rn(3,ngdim), d1rn_rj(3,3,2,ngdim), d1rn_dj(3,3,2,ngdim)
 
      real rs(3,ngdim), rs_rj(3,3,2,ngdim), rs_dj(3,3,2,ngdim),
     &                  rs_pj(3,2,ngdim)
      real d1rs(3,ngdim), d1rs_rj(3,3,2,ngdim), d1rs_dj(3,3,2,ngdim),
     &                    d1rs_pj(3,2,ngdim)

      real d1d(3,ngdim), d1d_rj(3,3,2,ngdim), d1d_dj(3,3,2,ngdim)

      real d1r(3,ngdim), d1r_rj(3,3,2,ngdim), d1r_dj(3,3,2,ngdim),
     &                   d1r_pj(3,2,ngdim)

      real bnb(ngdim), bnb_rj(3,2,ngdim), bnb_dj(3,2,ngdim)
      real d1bnb(ngdim), d1bnb_rj(3,2,ngdim), d1bnb_dj(3,2,ngdim)

      real be, be_rj(3,2), be_dj(3,2)

      real db(3), db_dj(3,3,2)
      real d1db(3), d1db_dj(3,3,2)

      real dba, dba_dj(3,2)
      real d1dba, d1dba_dj(3,2)

      real dxrj(3,2), dxrj_rj(3,2,3,2), dxrj_dj(3,2,3,2)
      real d1dxrj(3,2), d1dxrj_rj(3,2,3,2), d1dxrj_dj(3,2,3,2)

      real de(3), de_dj(3,3,2)
      real lh_a, lh_g

      real d1rbsq, d1rbsq_rj(3,4)

      real d1ra, d1ra_rj(3,4), d1ra_dj(3,4), d1ra_pj(4)

      real sb1, sb1_rj(3,4), sb1_dj(3,4)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c----------------------------------------------

      if(ng .gt. ngdim) then
       write(*,*) 'HSMGEO: Local work array overflow.',
     &  '  Increase  ngdim  to', ng
       stop
      endif

c---- put passed node quantities into more convenient node-indexed arrays
      do k = 1, 3
        rj(k,1) = r1(k)
        rj(k,2) = r2(k)

        dj(k,1) = d1(k)
        dj(k,2) = d2(k)
      enddo

      pj(1) = p1
      pj(2) = p2

c------------------------------------------------------
c---- coefficient of bubble contribution for this edge
c-    (director divergence along edge)
      be = (dj(1,2)-dj(1,1))*(rj(1,2)-rj(1,1))/2.0
     &   + (dj(2,2)-dj(2,1))*(rj(2,2)-rj(2,1))/2.0
     &   + (dj(3,2)-dj(3,1))*(rj(3,2)-rj(3,1))/2.0
      do jj = 1, 2
        do kk = 1, 3
          be_rj(kk,jj) = 0.
          be_dj(kk,jj) = 0.
        enddo
      enddo
      do kk = 1, 3
        be_rj(kk,2) =  (dj(kk,2)-dj(kk,1))/2.0
        be_rj(kk,1) = -(dj(kk,2)-dj(kk,1))/2.0
        be_dj(kk,2) =  (rj(kk,2)-rj(kk,1))/2.0
        be_dj(kk,1) = -(rj(kk,2)-rj(kk,1))/2.0
      enddo

c---- set director lean using MITC at edge midpoint
      desq = (dj(1,2)+dj(1,1))**2
     &     + (dj(2,2)+dj(2,1))**2
     &     + (dj(3,2)+dj(3,1))**2
      dea = sqrt(desq)
      do k = 1, 3
        de(k) = (dj(k,2) + dj(k,1)) / dea

        do jj = 1, 2
          de_dj(k,1,jj) = 0.
          de_dj(k,2,jj) = 0.
          de_dj(k,3,jj) = 0.
        enddo 

        do kk = 1, 3
          de_dj(k,kk,2) = -de(k)/desq * (dj(kk,2)+dj(kk,1))
          de_dj(k,kk,1) = -de(k)/desq * (dj(kk,2)+dj(kk,1))
        enddo
        de_dj(k,k,2) = de_dj(k,k,2) + 1.0/dea
        de_dj(k,k,1) = de_dj(k,k,1) + 1.0/dea
      enddo 

c---- director lean along edge is same as at midpoint
      do ig = 1, ng
        l(ig) = de(1)*(rj(1,2)-rj(1,1))/2.0
     &        + de(2)*(rj(2,2)-rj(2,1))/2.0
     &        + de(3)*(rj(3,2)-rj(3,1))/2.0
        do jj = 1, 2
          do kk = 1, 3
            l_rj(kk,jj,ig) = 0.
            l_dj(kk,jj,ig) = 
     &         de_dj(1,kk,jj)*(rj(1,2)-rj(1,1))/2.0
     &       + de_dj(2,kk,jj)*(rj(2,2)-rj(2,1))/2.0
     &       + de_dj(3,kk,jj)*(rj(3,2)-rj(3,1))/2.0
          enddo
          l_pj(jj,ig) = 0.
        enddo
        l_rj(1,2,ig) = l_rj(1,2,ig) + de(1)/2.0
        l_rj(2,2,ig) = l_rj(2,2,ig) + de(2)/2.0
        l_rj(3,2,ig) = l_rj(3,2,ig) + de(3)/2.0
        l_rj(1,1,ig) = l_rj(1,1,ig) - de(1)/2.0
        l_rj(2,1,ig) = l_rj(2,1,ig) - de(2)/2.0
        l_rj(3,1,ig) = l_rj(3,1,ig) - de(3)/2.0
      enddo

c======================================================================
c---- go over Gauss points
      do 100 ig = 1, ng

      xsi = xsig(ig)

c----------------------------------------------------------------
c---- linear interpolation weights "N" at gauss point
      nn(1) = 0.5*(1.0 - xsi)
      nn(2) = 0.5*(1.0 + xsi)

c---- xsi derivatives of N
      d1nn(1) = -0.5
      d1nn(2) =  0.5
c
c---- quadratic edge weight NB at Gauss point
      nnb = 0.25*(1.0 - xsi**2)
      d1nnb = -0.5*xsi

c----------------------------------------------------------------
c---- bilinear rb vector
      do k = 1, 3
        rb(k,ig) = rj(k,1)*nn(1)
     &           + rj(k,2)*nn(2)
        d1rb(k,ig) = rj(k,1)*d1nn(1)
     &             + rj(k,2)*d1nn(2)
        do jj = 1, 2
          do kk = 1, 3
            rb_rj(k,kk,jj,ig) = 0.
            d1rb_rj(k,kk,jj,ig) = 0.
          enddo
          rb_rj(k,k,jj,ig) = nn(jj)
          d1rb_rj(k,k,jj,ig) = d1nn(jj)
        enddo
      enddo

c---- blinear d vector
      do k = 1, 3
        db(k) = dj(k,1)*nn(1)
     &          + dj(k,2)*nn(2)
        d1db(k) = dj(k,1)*d1nn(1)
     &            + dj(k,2)*d1nn(2)
        do jj = 1, 2
          do kk = 1, 3
            db_dj(k,kk,jj) = 0.
            d1db_dj(k,kk,jj) = 0.
          enddo
          db_dj(k,k,jj) = nn(jj)
          d1db_dj(k,k,jj) = d1nn(jj)
        enddo
      enddo

c---- set db magnitude, and normalize to get d
      dba = sqrt(db(1)**2 + db(2)**2 + db(3)**2)
      d1dba = ( db(1)*d1db(1)
     &       + db(2)*d1db(2)
     &       + db(3)*d1db(3) )/dba
      do jj = 1, 2
        do kk = 1, 3
          dba_dj(kk,jj) = nn(jj)*db(kk)/dba
          d1dba_dj(kk,jj) = (nn(jj)*d1db(kk) + db(kk)*d1nn(jj))/dba
     &                   - d1dba/dba*dba_dj(kk,jj)
        enddo
      enddo

c---- interpolated unit d vector
      do k = 1, 3
        d(k,ig) = db(k)/dba
        d1d(k,ig) = (d1db(k) - d(k,ig)*d1dba)/dba
        do jj = 1, 2
          do kk = 1, 3
            d_rj(k,kk,jj,ig) = 0.
            d1d_rj(k,kk,jj,ig) = 0.

            d_dj(k,kk,jj,ig) = ( db_dj(k,kk,jj)
     &                         - d(k,ig)*dba_dj(kk,jj) )/dba
            d1d_dj(k,kk,jj,ig) = d1db_dj(k,kk,jj)/dba
     &                       - ( d_dj(k,kk,jj,ig)*d1dba
     &                          +d(k,ig)         *d1dba_dj(kk,jj) )/dba
     &                       - d1d(k,ig)*dba_dj(kk,jj)/dba
          enddo
          d_pj(k,jj,ig) = 0.
        enddo
      enddo

c----------------------------------------------------------------
      if(lrcurv) then
c----- add on quadratic bubble part to position and covariant vectors
       bnb(ig) = be*nnb
       d1bnb(ig) = be*d1nnb
       do jj = 1, 2
         do kk = 1, 3
           bnb_rj(kk,jj,ig) = be_rj(kk,jj)*nnb
           bnb_dj(kk,jj,ig) = be_dj(kk,jj)*nnb
           d1bnb_rj(kk,jj,ig) = be_rj(kk,jj)*d1nnb
           d1bnb_dj(kk,jj,ig) = be_dj(kk,jj)*d1nnb
         enddo
       enddo

       do k = 1, 3
         rn(k,ig) = bnb(ig)*d(k,ig)
         d1rn(k,ig) = d1bnb(ig)*d(k,ig) + bnb(ig)*d1d(k,ig)
         do jj = 1, 2
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = bnb_rj(kk,jj,ig)*d(k,ig)
     &                         + bnb(ig)         *d_rj(k,kk,jj,ig)
             rn_dj(k,kk,jj,ig) = bnb_dj(kk,jj,ig)*d(k,ig)
     &                         + bnb(ig)         *d_dj(k,kk,jj,ig)

             d1rn_rj(k,kk,jj,ig) = d1bnb_rj(kk,jj,ig)*d(k,ig)
     &                           + bnb_rj(kk,jj,ig)*d1d(k,ig)
     &                           + bnb(ig)         *d1d_rj(k,kk,jj,ig)
             d1rn_dj(k,kk,jj,ig) = d1bnb_dj(kk,jj,ig)*d(k,ig)
     &                           + d1bnb(ig)         *d_dj(k,kk,jj,ig)
     &                           + bnb_dj(kk,jj,ig)*d1d(k,ig)
     &                           + bnb(ig)         *d1d_dj(k,kk,jj,ig)
           enddo
         enddo
       enddo

       if(.false.) then
       do k = 1, 3
         rn(k,ig) = bnb(ig)*db(k)
         d1rn(k,ig) = d1bnb(ig)*db(k) + bnb(ig)*d1db(k)
         do jj = 1, 2
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = bnb_rj(kk,jj,ig)*db(k)
             rn_dj(k,kk,jj,ig) = bnb_dj(kk,jj,ig)*db(k)
     &                         + bnb(ig)         *db_dj(k,kk,jj)

             d1rn_rj(k,kk,jj,ig) = d1bnb_rj(kk,jj,ig)*db(k)
     &                           + bnb_rj(kk,jj,ig)*d1db(k)
             d1rn_dj(k,kk,jj,ig) = d1bnb_dj(kk,jj,ig)*db(k)
     &                           + d1bnb(ig)         *db_dj(k,kk,jj)
     &                           + bnb_dj(kk,jj,ig)*d1db(k)
     &                           + bnb(ig)         *d1db_dj(k,kk,jj)
           enddo
         enddo
       enddo
       endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
c----- bilinear: no bubble offset
       bnb(ig) = 0.
       d1bnb(ig) = 0.
       do jj = 1, 2
         do kk = 1, 3
           bnb_rj(kk,jj,ig) = 0.
           bnb_dj(kk,jj,ig) = 0.
           d1bnb_rj(kk,jj,ig) = 0.
           d1bnb_dj(kk,jj,ig) = 0.
         enddo
       enddo

       do k = 1, 3
         rn(k,ig) = 0.
         d1rn(k,ig) = 0.

         do jj = 1, 2
           do kk = 1, 3
             rn_rj(k,kk,jj,ig) = 0.
             rn_dj(k,kk,jj,ig) = 0.

             d1rn_rj(k,kk,jj,ig) = 0.
             d1rn_dj(k,kk,jj,ig) = 0.
           enddo
         enddo
       enddo

      endif ! lrcurv


c----------------------------------------------------------------

      if(ldrill) then
c----- quadratic in-surface drilling geometry offset
       do j = 1, 2
         do k = 1, 3
           ic = icrs(k)
           jc = jcrs(k)
           dxrj(k,j)   = dj(ic,j)*(rb(jc,ig)-rj(jc,j))
     &                 - dj(jc,j)*(rb(ic,ig)-rj(ic,j))
           d1dxrj(k,j) = dj(ic,j)*d1rb(jc,ig)
     &                 - dj(jc,j)*d1rb(ic,ig)
           do jj = 1, 2
             do kk = 1, 3
               dxrj_rj(k,j,kk,jj) = dj(ic,j)*rb_rj(jc,kk,jj,ig)
     &                            - dj(jc,j)*rb_rj(ic,kk,jj,ig)
               d1dxrj_rj(k,j,kk,jj) = dj(ic,j)*d1rb_rj(jc,kk,jj,ig)
     &                              - dj(jc,j)*d1rb_rj(ic,kk,jj,ig)
               dxrj_dj(k,j,kk,jj) = 0.
               d1dxrj_dj(k,j,kk,jj) = 0.
             enddo
           enddo
           dxrj_rj(k,j,ic,j) = dxrj_rj(k,j,ic,j) + dj(jc,j)
           dxrj_rj(k,j,jc,j) = dxrj_rj(k,j,jc,j) - dj(ic,j)
           dxrj_dj(k,j,ic,j) = dxrj_dj(k,j,ic,j) + (rb(jc,ig)-rj(jc,j))
           dxrj_dj(k,j,jc,j) = dxrj_dj(k,j,jc,j) - (rb(ic,ig)-rj(ic,j))

           d1dxrj_dj(k,j,ic,j) = d1dxrj_dj(k,j,ic,j) + d1rb(jc,ig)
           d1dxrj_dj(k,j,jc,j) = d1dxrj_dj(k,j,jc,j) - d1rb(ic,ig)
         enddo
       enddo

       do k = 1, 3
         rs(k,ig) = pj(1)*dxrj(k,1)*nn(1)
     &            + pj(2)*dxrj(k,2)*nn(2)
         d1rs(k,ig) = pj(1)*d1dxrj(k,1)*nn(1)
     &              + pj(2)*d1dxrj(k,2)*nn(2)
     &              + pj(1)*dxrj(k,1)*d1nn(1)
     &              + pj(2)*dxrj(k,2)*d1nn(2)
         do jj = 1, 2
           do kk = 1, 3
             rs_rj(k,kk,jj,ig) =
     &               pj(1)*dxrj_rj(k,1,kk,jj)*nn(1)
     &             + pj(2)*dxrj_rj(k,2,kk,jj)*nn(2)
             rs_dj(k,kk,jj,ig) =
     &               pj(1)*dxrj_dj(k,1,kk,jj)*nn(1)
     &             + pj(2)*dxrj_dj(k,2,kk,jj)*nn(2)

             d1rs_rj(k,kk,jj,ig) = 
     &                pj(1)*d1dxrj_rj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d1dxrj_rj(k,2,kk,jj)*nn(2)
     &              + pj(1)*dxrj_rj(k,1,kk,jj)*d1nn(1)
     &              + pj(2)*dxrj_rj(k,2,kk,jj)*d1nn(2)

             d1rs_dj(k,kk,jj,ig) = 
     &                pj(1)*d1dxrj_dj(k,1,kk,jj)*nn(1)
     &              + pj(2)*d1dxrj_dj(k,2,kk,jj)*nn(2)
     &              + pj(1)*dxrj_dj(k,1,kk,jj)*d1nn(1)
     &              + pj(2)*dxrj_dj(k,2,kk,jj)*d1nn(2)
           enddo
           rs_pj(k,jj,ig) = dxrj(k,jj)*nn(jj)
           d1rs_pj(k,jj,ig) = d1dxrj(k,jj)*nn(jj)
     &                      + dxrj(k,jj)*d1nn(jj)
         enddo
       enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
c----- no drill offset
       do k = 1, 3
         rs(k,ig) = 0.
         d1rs(k,ig) = 0.
         do jj = 1, 2
           do kk = 1, 3
             rs_rj(k,kk,jj,ig) = 0.
             rs_dj(k,kk,jj,ig) = 0.
             d1rs_rj(k,kk,jj,ig) = 0.
             d1rs_dj(k,kk,jj,ig) = 0.
           enddo
           rs_pj(k,jj,ig) = 0.
           d1rs_pj(k,jj,ig) = 0.
         enddo
       enddo

      endif ! ldrill
c--------------------------------------------------------------------

 100  continue ! next Gauss point

c==========================================================================

      asum = 0.
      do ig = 1, ng
        asum = asum + fwtg(ig)
      enddo

c---- integrate to get bubble basis vector magnitudes sb1,sb2
      sb1 = 0.
      do jj = 1, 2
        do kk = 1, 3
          sb1_rj(kk,jj) = 0.
          sb1_dj(kk,jj) = 0.
        enddo
      enddo

      do ig = 1, ng
c------ Gauss weight wk/A 
        fwt = fwtg(ig) / asum

        d1rbsq = d1rb(1,ig)**2
     &         + d1rb(2,ig)**2
     &         + d1rb(3,ig)**2
        do jj = 1, 2
          do kk = 1, 3
            d1rbsq_rj(kk,jj) = 2.0*d1rb(1,ig)*d1rb_rj(1,kk,jj,ig)
     &                       + 2.0*d1rb(2,ig)*d1rb_rj(2,kk,jj,ig)
     &                       + 2.0*d1rb(3,ig)*d1rb_rj(3,kk,jj,ig)
          enddo
        enddo

        fb1 = sqrt(d1rbsq + d1bnb(ig)**2)

        fb1_d1rbsq = 0.5/fb1
        fb1_d1bnb = d1bnb(ig)/fb1

        sb1 = sb1 + fb1*fwt 

        do jj = 1, 2
          do kk = 1, 3
            sb1_rj(kk,jj) = 
     &      sb1_rj(kk,jj) + (fb1_d1rbsq*d1rbsq_rj(kk,jj)
     &                    +  fb1_d1bnb*d1bnb_rj(kk,jj,ig) ) * fwt
            sb1_dj(kk,jj) = 
     &      sb1_dj(kk,jj) +  fb1_d1bnb*d1bnb_dj(kk,jj,ig)   * fwt
          enddo
        enddo
        
      enddo

c==================================================================================
c---- set final quantities at Gauss points
      do 200 ig = 1, ng

c--------------------------------------------------------------------
c----- geometry and derivatives
       do k = 1, 3
         r(k,ig) = rb(k,ig) + rn(k,ig) + rs(k,ig)
         d1r(k,ig) = d1rb(k,ig) + d1rn(k,ig) + d1rs(k,ig)

         do jj = 1, 2
           do kk = 1, 3
             r_rj(k,kk,jj,ig) = rb_rj(k,kk,jj,ig)
     &                        + rs_rj(k,kk,jj,ig)
     &                        + rn_rj(k,kk,jj,ig)
             d1r_rj(k,kk,jj,ig) = d1rb_rj(k,kk,jj,ig)
     &                          + d1rn_rj(k,kk,jj,ig)
     &                          + d1rs_rj(k,kk,jj,ig)

             r_dj(k,kk,jj,ig) = rs_dj(k,kk,jj,ig)
     &                        + rn_dj(k,kk,jj,ig)
             d1r_dj(k,kk,jj,ig) = d1rn_dj(k,kk,jj,ig)
     &                          + d1rs_dj(k,kk,jj,ig)
           enddo
           r_pj(k,jj,ig) = rs_pj(k,jj,ig)
           d1r_pj(k,jj,ig) = d1rs_pj(k,jj,ig)
         enddo
       enddo

c--------------------------------------------------------------------
c----- basis vectors
       d1ra = sqrt(d1r(1,ig)**2 + d1r(2,ig)**2 + d1r(3,ig)**2)
       do jj = 1, 2
         do kk = 1, 3
           d1ra_rj(kk,jj) = ( d1r(1,ig)*d1r_rj(1,kk,jj,ig)
     &                      + d1r(2,ig)*d1r_rj(2,kk,jj,ig)
     &                      + d1r(3,ig)*d1r_rj(3,kk,jj,ig) ) / d1ra
           d1ra_dj(kk,jj) = ( d1r(1,ig)*d1r_dj(1,kk,jj,ig)
     &                      + d1r(2,ig)*d1r_dj(2,kk,jj,ig)
     &                      + d1r(3,ig)*d1r_dj(3,kk,jj,ig) ) / d1ra
         enddo
         d1ra_pj(jj) = ( d1r(1,ig)*d1r_pj(1,jj,ig)
     &                 + d1r(2,ig)*d1r_pj(2,jj,ig)
     &                 + d1r(3,ig)*d1r_pj(3,jj,ig) ) / d1ra
       enddo

       do k = 1, 3
         a1(k,ig) = d1r(k,ig)
         do jj = 1, 2
           do kk = 1, 3
             a1_rj(k,kk,jj,ig) = d1r_rj(k,kk,jj,ig)
             a1_dj(k,kk,jj,ig) = d1r_dj(k,kk,jj,ig)
           enddo
           a1_pj(k,jj,ig) = d1r_pj(k,jj,ig)
         enddo
       enddo

c----------------------------------------------------------------
c---- covariant metric "tensor" g_ij, and Jacobians
      g(ig) = d1rb(1,ig)*d1rb(1,ig)
     &      + d1rb(2,ig)*d1rb(2,ig)
     &      + d1rb(3,ig)*d1rb(3,ig)
      do jj = 1, 2
        do kk = 1, 3
          g_rj(kk,jj,ig) = 2.0*d1rb(1,ig)*d1rb_rj(1,kk,jj,ig)
     &                   + 2.0*d1rb(2,ig)*d1rb_rj(2,kk,jj,ig)
     &                   + 2.0*d1rb(3,ig)*d1rb_rj(3,kk,jj,ig)
          g_dj(kk,jj,ig) = 0.
        enddo
        g_pj(jj,ig) = 0.
      enddo

c----------------------------------------------------------------
c---- curvature and Jacobians
      h(ig) = d1rb(1,ig)*d1d(1,ig)
     &      + d1rb(2,ig)*d1d(2,ig)
     &      + d1rb(3,ig)*d1d(3,ig)
      do jj = 1, 2
        do kk = 1, 3
          h_rj(kk,jj,ig) = d1rb_rj(1,kk,jj,ig)*d1d(1,ig)
     &                   + d1rb_rj(2,kk,jj,ig)*d1d(2,ig)
     &                   + d1rb_rj(3,kk,jj,ig)*d1d(3,ig)
     &                   + d1rb(1,ig)         *d1d_rj(1,kk,jj,ig)
     &                   + d1rb(2,ig)         *d1d_rj(2,kk,jj,ig)
     &                   + d1rb(3,ig)         *d1d_rj(3,kk,jj,ig)
          h_dj(kk,jj,ig) = d1rb(1,ig)         *d1d_dj(1,kk,jj,ig)
     &                   + d1rb(2,ig)         *d1d_dj(2,kk,jj,ig)
     &                   + d1rb(3,ig)         *d1d_dj(3,kk,jj,ig)
        enddo
        h_pj(jj,ig) = 0.
      enddo

c----------------------------------------------------------------
c---- unit edge-tangential vector l = a1 / sqrt(g)
      do k = 1, 3
        lh(k,ig) = a1(k,ig) / sqrt(g(ig))
        lh_a =   1.0 / sqrt(g(ig))
        lh_g = -0.5*lh(k,ig)/g(ig)
        do jj = 1, 2
          do kk = 1, 3
            lh_rj(k,kk,jj,ig) = lh_a*a1_rj(k,kk,jj,ig)
     &                        + lh_g* g_rj(kk,jj,ig)
            lh_dj(k,kk,jj,ig) = lh_a*a1_dj(k,kk,jj,ig)
     &                        + lh_g* g_dj(kk,jj,ig)
          enddo
          lh_pj(k,jj,ig) = lh_a*a1_pj(k,jj,ig) + lh_g*g_pj(jj,ig)
        enddo
      enddo

c----------------------------------------------------------------
c---- unit edge-normal vector,  t = l x d
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        th(k,ig) = lh(ic,ig)*d(jc,ig)
     &           - lh(jc,ig)*d(ic,ig)

        do jj = 1, 2
          do kk = 1, 3
            th_rj(k,kk,jj,ig) = lh_rj(ic,kk,jj,ig)*d(jc,ig)
     &                        - lh_rj(jc,kk,jj,ig)*d(ic,ig)
     &                        + lh(ic,ig)*d_rj(jc,kk,jj,ig)
     &                        - lh(jc,ig)*d_rj(ic,kk,jj,ig)
            th_dj(k,kk,jj,ig) = lh_dj(ic,kk,jj,ig)*d(jc,ig)
     &                        - lh_dj(jc,kk,jj,ig)*d(ic,ig)
     &                        + lh(ic,ig)*d_dj(jc,kk,jj,ig)
     &                        - lh(jc,ig)*d_dj(ic,kk,jj,ig)
          enddo
          th_pj(k,jj,ig) = lh_pj(ic,jj,ig)*d(jc,ig)
     &                   - lh_pj(jc,jj,ig)*d(ic,ig)
     &                   + lh(ic,ig)*d_pj(jc,jj,ig)
     &                   - lh(jc,ig)*d_pj(ic,jj,ig)
        enddo
      enddo
c----------------------------------------------------------------
 200  continue ! next Gauss point

      return
      end ! hsmgeo2

