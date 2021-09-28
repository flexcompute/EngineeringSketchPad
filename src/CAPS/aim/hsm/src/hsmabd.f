
       subroutine hsmabd(neln,
     &                   par1    , par2    , par3    , par4    , 
     &                   aaj, bbj, ddj, ssj )
c---------------------------------------------------------------------------
c     Calculates nodal stiffness and compliance matrices at the nodes
c     of a tri or quad element, in the element's local cartesian basis.
c
c Inputs:
c   neln        number of nodes  (3 = tri, 4 = quad)
c   par#(.)   node parameters  (par4 not used if neln = 3)
c
c Outputs:
c   aaj(....j)   A stiffness tensors at nodes
c   bbj(....j)   B stiffness tensors at nodes
c   ddj(....j)   D stiffness tensors at nodes
c   ssj(....j)   S transverse shear stiffness tensors at nodes
c---------------------------------------------------------------------------
      implicit real (a-h,m,o-z)
      include 'index.inc'

c----------------------------------------------
c---- passed variables
      real par1(lvtot),
     &     par2(lvtot),
     &     par3(lvtot),
     &     par4(lvtot)

      real aaj(2,2,2,2,4), bbj(2,2,2,2,4),  ddj(2,2,2,2,4)
      real ssj(2,2,4)

c----------------------------------------------
c---- local arrays
      real par(lvtot,4)

      real r0j(3,4)
      real n0j(3,4), e01j(3,4), e02j(3,4)

      real d1nn(4), d2nn(4)

      real d1r0(3),
     &     d2r0(3)

      real crs(3)
      real c1(3)
      real c2(3)
      real c3(3)

      real cdn, cxn(3)
      real v(3)
      real lsq, lam(3)
      real vxc1(3), vxc2(3), vxc3(3)
      real c1j(3), c2j(3), c3j(3)

      real ec(2,2)

c----------------------------------------------
c---- cross-product indices
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 /
      data jcrs / 3, 1, 2 /
c----------------------------------------------
c---- approximate machine zero
      data epsm / 1.0e-9 /
c     data epsm / 1.0e-18 /
c----------------------------------------------

c---- put passed node quantities into more convenient node-indexed arrays
      do lv = 1, lvtot
        par(lv,1) = par1(lv)
        par(lv,2) = par2(lv)
        par(lv,3) = par3(lv)
        par(lv,4) = par4(lv)
      enddo
c
c---- further unpack vectors
      do j = 1, 4
        r0j(1,j) = par(lvr0x,j)
        r0j(2,j) = par(lvr0y,j)
        r0j(3,j) = par(lvr0z,j)

        n0j(1,j) = par(lvn0x,j)
        n0j(2,j) = par(lvn0y,j)
        n0j(3,j) = par(lvn0z,j)

        e01j(1,j) = par(lve01x,j)
        e01j(2,j) = par(lve01y,j)
        e01j(3,j) = par(lve01z,j)

        e02j(1,j) = par(lve02x,j)
        e02j(2,j) = par(lve02y,j)
        e02j(3,j) = par(lve02z,j)
      enddo

c---- xsi, eta derivatives of interpolation weights at cell centroid
      if(neln.eq.3) then
       xsi = 0.3333333333333333333333333
       eta = 0.3333333333333333333333333
c      nn(1) = 1.0 - xsi - eta
c      nn(2) = xsi
c      nn(3) = eta
c      nn(4) = 0.

       d1nn(1) = -1.0
       d1nn(2) =  1.0
       d1nn(3) =  0.
       d1nn(4) = 0.

       d2nn(1) = -1.0
       d2nn(2) =  0.
       d2nn(3) =  1.0
       d2nn(4) = 0.

      else
       xsi = 0.0
       eta = 0.0
c      nn(1) = 0.25*(1.0 - xsi)*(1.0 - eta)
c      nn(2) = 0.25*(1.0 + xsi)*(1.0 - eta)
c      nn(3) = 0.25*(1.0 + xsi)*(1.0 + eta)
c      nn(4) = 0.25*(1.0 - xsi)*(1.0 + eta)

       d1nn(1) = -0.25*(1.0 - eta)
       d1nn(2) =  0.25*(1.0 - eta)
       d1nn(3) =  0.25*(1.0 + eta)
       d1nn(4) = -0.25*(1.0 + eta)

       d2nn(1) = -0.25*(1.0 - xsi)
       d2nn(2) = -0.25*(1.0 + xsi)
       d2nn(3) =  0.25*(1.0 + xsi)
       d2nn(4) =  0.25*(1.0 - xsi)

      endif

c---- centroid dr0/dxsi, dr0/deta
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

c---- element basis vectors
      d1ra = sqrt(d1r0(1)**2 + d1r0(2)**2 + d1r0(3)**2)
      c1(1) = d1r0(1)/d1ra
      c1(2) = d1r0(2)/d1ra
      c1(3) = d1r0(3)/d1ra

      crs(1) = d1r0(2)*d2r0(3) - d1r0(3)*d2r0(2)
      crs(2) = d1r0(3)*d2r0(1) - d1r0(1)*d2r0(3)
      crs(3) = d1r0(1)*d2r0(2) - d1r0(2)*d2r0(1)
      crsa = sqrt(crs(1)**2 + crs(2)**2 + crs(3)**2)
      c3(1) = crs(1)/crsa
      c3(2) = crs(2)/crsa
      c3(3) = crs(3)/crsa

      c2(1) = c3(2)*c1(3) - c3(3)*c1(2)
      c2(2) = c3(3)*c1(1) - c3(1)*c1(3)
      c2(3) = c3(1)*c1(2) - c3(2)*c1(1)

c---- go over the nodes
      do j = 1, neln

c-------------------------------------------------------------
c------ rotate c1,c2,c3 so c3 is aligned with n0j

c------ cxn = c3 x n0j
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          cxn(k) = c3(ic)*n0j(jc,j)
     &           - c3(jc)*n0j(ic,j)
        enddo

c------ |cxn| = |c||n| sin(t)
        cxna = sqrt(cxn(1)**2 + cxn(2)**2 + cxn(3)**2)

c------ c3.n0j = |c||n| cos(t)
        cdn = c3(1)*n0j(1,j)
     &      + c3(2)*n0j(2,j)
     &      + c3(3)*n0j(3,j)

c------ theta = atan[ |c||n| sin(t)  /  |c||n| cos(t) ]
        theta = atan2( cxna , cdn )

c----------------------------------------
c------ rotate c1,c2 by theta to get c1j,c2j vectors

        w = cos(0.5*theta)

        if(theta**2 .lt. epsm) then
c------- v for small theta case,   assumes  |c||n| = 1
         do k = 1, 3
           v(k) = 0.5*cxn(k)
         enddo
        else
c------- v for general case
         do k = 1, 3
           v(k) = sin(0.5*theta)*cxn(k)/cxna
         enddo

c- - - - - - - - - - - - - - - - - - - -
        endif

c------ v x c1, v x c2
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          vxc1(k) = v(ic)*c1(jc) - v(jc)*c1(ic)
          vxc2(k) = v(ic)*c2(jc) - v(jc)*c2(ic)
c         vxc3(k) = v(ic)*c3(jc) - v(jc)*c3(ic)
        enddo

c------ c1j = c1 +  2 v x ( v x c1  +  w c1 )
c-      c2j = c2 +  2 v x ( v x c2  +  w c2 )
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)

          c1j(k) = c1(k)
     &           + 2.0*v(ic)*(vxc1(jc) + w*c1(jc))
     &           - 2.0*v(jc)*(vxc1(ic) + w*c1(ic))
          c2j(k) = c2(k)
     &           + 2.0*v(ic)*(vxc2(jc) + w*c2(jc))
     &           - 2.0*v(jc)*(vxc2(ic) + w*c2(ic))
c         c3j(k) = c3(k)
c    &           + 2.0*v(ic)*(vxc3(jc) + w*c3(jc))
c    &           - 2.0*v(jc)*(vxc3(ic) + w*c3(ic))
          c3j(k) = n0j(k,j)
        enddo

c------ basis vector dot products
        ec(1,1) = e01j(1,j)*c1j(1)
     &          + e01j(2,j)*c1j(2)
     &          + e01j(3,j)*c1j(3)
        ec(1,2) = e01j(1,j)*c2j(1)
     &          + e01j(2,j)*c2j(2)
     &          + e01j(3,j)*c2j(3)
        ec(2,1) = e02j(1,j)*c1j(1)
     &          + e02j(2,j)*c1j(2)
     &          + e02j(3,j)*c1j(3)
        ec(2,2) = e02j(1,j)*c2j(1)
     &          + e02j(2,j)*c2j(2)
     &          + e02j(3,j)*c2j(3)


c------ Voigt contraction 
c       1 1 = 11 11
c       2 2 = 22 22
C
c       1 2 = 11 22
c       2 1 = 22 11
C
c       1 6 = 11 12
c       1 6 = 11 21
c       6 1 = 12 11
c       6 1 = 21 11
C
c       2 6 = 22 12
c       2 6 = 22 21
c       6 2 = 12 22
c       6 2 = 21 22
C
c       6 6 = 12 12
c       6 6 = 21 21
c       6 6 = 12 21
c       6 6 = 21 12
c
c       5 5 = 1n 1n
c       5 5 = n1 n1
c
c       4 4 = 2n 2n
c       4 4 = n2 n2
c
        do id = 1, 2
          do ic = 1, 2
            do ib = 1, 2
              do ia = 1, 2
                aaj(ia,ib,ic,id,j) =
     &                par(lva11,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lva22,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lva12,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lva12,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lva16,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
     &              + par(lva16,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
     &              + par(lva16,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lva16,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lva26,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lva26,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
     &              + par(lva26,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lva26,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lva66,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lva66,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lva66,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lva66,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
                bbj(ia,ib,ic,id,j) =                              
     &                par(lvb11,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvb22,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvb12,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvb12,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvb16,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvb16,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
     &              + par(lvb16,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvb16,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvb26,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvb26,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
     &              + par(lvb26,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvb26,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvb66,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvb66,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lvb66,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lvb66,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
                ddj(ia,ib,ic,id,j) =                              
     &                par(lvd11,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvd22,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvd12,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvd12,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvd16,j)*ec(1,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvd16,j)*ec(1,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
     &              + par(lvd16,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvd16,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(1,id)
     &              + par(lvd26,j)*ec(2,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvd26,j)*ec(2,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
     &              + par(lvd26,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvd26,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(2,id)
     &              + par(lvd66,j)*ec(1,ia)*ec(2,ib)*ec(1,ic)*ec(2,id)
     &              + par(lvd66,j)*ec(2,ia)*ec(1,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lvd66,j)*ec(1,ia)*ec(2,ib)*ec(2,ic)*ec(1,id)
c    &              + par(lvd66,j)*ec(2,ia)*ec(1,ib)*ec(1,ic)*ec(2,id)
              enddo
            enddo
          enddo
        enddo

        do ia = 1, 2
          do ib = 1, 2
            ssj(ia,ib,j) =
     &                par(lva55,j)*ec(1,ia)*ec(1,ib)
     &              + par(lva44,j)*ec(2,ia)*ec(2,ib)
          enddo
        enddo

      enddo ! next j

      if(neln.eq.3) then
c------ zero out dummy j=4 matrices
        j = 4
        do id = 1, 2
          do ic = 1, 2
            do ib = 1, 2
              do ia = 1, 2
                aaj(ia,ib,ic,id,j) = 0.
                bbj(ia,ib,ic,id,j) = 0.
                ddj(ia,ib,ic,id,j) = 0.
              enddo
            enddo
          enddo
        enddo

        do ia = 1, 2
          do ib = 1, 2
            ssj(ia,ib,j) = 0.
          enddo
        enddo
      endif

      return
      end ! hsmabd
