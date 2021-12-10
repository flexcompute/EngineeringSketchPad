
      subroutine sblud(nbdim,nblk,namat,ifrst,ilast,mfrst, a)
c-------------------------------------------------------------------------
c     Performs LU decomposition on a sparse block matrix in column storage.
c
c     This routine does not create new matrix elements due to fill-in.
c     Hence, the matrix array must be set up to contain all fill-in 
c     that will occur, with the expected fill-in blocks set to zero.
c
c     If no such fill-in storage is provided, the fill-in blocks will
c     be ignored.  In this case, this routine will in effect return 
c     an ILU decomposition.
c
c   Inputs:
c     nbdim     dimension of blocks
c     nblk      size of blocks
c     namat     size of matrix (number of blocks)
c     ifrst(j)  row index of topmost block in column j
c     ilast(j)  row index of bottommost block in column j
c     mfrst(j)  m index of block matrix element a(..m) at ifrst,j
c     a(..m)    matrix elements stored as blocks in columns
c
c   Outputs:
c     a(..m)    LU-decomposed matrix, to be used in subroutine SBBKS
c
c
c   Note:  The block indices in column j are  
c
c            m  =  mfrst(j) ... mfrst(j+1)-1
c
c          So the dummy-block index mfrst(namat+1) must be set accordingly.
c
c-------------------------------------------------------------------------
      integer ifrst(*),ilast(*),mfrst(*)
      real a(nbdim,nbdim,*)

 1200 format(1x, 3(i4,i4,3x) )

c      do j = 1, namat
c        ilast(j) = ifrst(j) + mfrst(j+1) - mfrst(j) - 1
c      enddo

c---- go over all pivot rows
      do 1000 ipiv = 1, namat
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks below diagonal
        do 100 je = 1, jpiv-1
c-------- if this block does not exist, don`t do anything
          if(ipiv .gt. ilast(je)) go to 100

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot block row
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify remaining blocks in pivot block row
          do 10 j = je+1, namat
c---------- if modifying block does not exist, don`t do anything
            if(imod .lt. ifrst(j)) go to 10

            m    = mfrst(j) + ipiv - ifrst(j)
            mmod = mfrst(j) + imod - ifrst(j)

            if(m .gt. mfrst(j+1)-1) then
             write(*,*) '? sblud: no fill space allocated, m > mlast(j)'
             write(*,*) '         ipiv, je, j =', ipiv, je, j
            endif

c---------- BLAS could be used for this matrix-matrix multiply,add
c           A  <--  A - Aeli*Amod
            do jj = 1, nblk
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*a(kk,jj,mmod)
                enddo
                a(ii,jj,m) = a(ii,jj,m) - sum
              enddo
            enddo

 10       continue ! next j
c- - - - - - - - - - - - - - - - - - - - - - - - - 
 100    continue ! next je (stops at jpiv-1)

c       mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
c       write(*,*) 'LU ...', ipiv, nblk, mpiv
c       do k = 1, nblk
c         write(*,'(1x,8g12.4)') (a(k,l,mpiv), l=1, nblk)
c       enddo

c------ LU-decompose diagonal pivot block Apiv
        mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
        call ludcmp(nbdim,nblk,a(1,1,mpiv))

c------------------------------------------------------------------------
c------ multiply remaining blocks in pivot row by Apiv^-1
        do 150 jmul = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(jmul)) go to 150

          mmul = mfrst(jmul) + ipiv - ifrst(jmul)
          call baksubb(nbdim,nblk,a(1,1,mpiv),a(1,1,mmul))
 150    continue
c------------------------------------------------------------------------
 1000 continue ! next pivot row

      return
      end ! sblud



      subroutine sbludb(nbdim,nblk,namat,ifrst,ilast,mfrst, a, b)
c-------------------------------------------------------------------------
c     Combines action of SBLUD and SBBKS
c
c   Inputs:
c     nbdim     dimension of blocks
c     nblk      size of blocks
c     namat     size of matrix (number of blocks)
c     ifrst(j)  row index of topmost block in column j
c     ilast(j)  row index of bottommost block in column j
c     mfrst(j)  m index of block matrix element a(..m) at ifrst,j
c     a(..m)    matrix elements stored as blocks in columns
c     b(.i)     r.h.s. vectors
c
c   Outputs:
c     a(..m)    LU-decomposed matrix, to be used in subroutine SBBKS
c     b(.i)     Solution vectors
c
c
c   Note:  The block indices in column j are  
c
c            m  =  mfrst(j) ... mfrst(j+1)-1
c
c          So the dummy-block index mfrst(namat+1) must be set accordingly.
c
c-------------------------------------------------------------------------
      integer ifrst(*),ilast(*),mfrst(*)
      real a(nbdim,nbdim,*)
      real b(nbdim,*)

 1200 format(1x, 3(i4,i4,3x) )

c      do j = 1, namat
c        ilast(j) = ifrst(j) + mfrst(j+1) - mfrst(j) - 1
c      enddo

c---- go over all pivot rows
      do 1000 ipiv = 1, namat
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks below diagonal
        do 100 je = 1, jpiv-1
c-------- if this block does not exist, don`t do anything
          if(ipiv .gt. ilast(je)) go to 100

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot block row
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify remaining blocks in pivot block row
          do 10 j = je+1, namat
c---------- if modifying block does not exist, don`t do anything
            if(imod .lt. ifrst(j)) go to 10

            m    = mfrst(j) + ipiv - ifrst(j)
            mmod = mfrst(j) + imod - ifrst(j)

            if(m .gt. mfrst(j+1)-1) then
             write(*,*) '? sblud: no fill space allocated, m > mlast(j)'
             write(*,*) '         ipiv, je, j =', ipiv, je, j
            endif

c---------- BLAS could be used for this matrix-matrix multiply,add
c           A  <--  A - Aeli*Amod
            do jj = 1, nblk
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*a(kk,jj,mmod)
                enddo
                a(ii,jj,m) = a(ii,jj,m) - sum
              enddo
            enddo

 10       continue ! next j
c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - 
 100    continue ! next je (stops at jpiv-1)

c------ LU-decompose diagonal pivot block Apiv
        mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
        call ludcmp(nbdim,nblk,a(1,1,mpiv))

c------------------------------------------------------------------------
c------ multiply remaining blocks in pivot row by Apiv^-1
        do 150 jmul = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(jmul)) go to 150

          mmul = mfrst(jmul) + ipiv - ifrst(jmul)
          call baksubb(nbdim,nblk,a(1,1,mpiv),a(1,1,mmul))
 150    continue

c------------------------------------------------------------------------
c------ multiply r.h.s. vector by Apiv^-1
        call baksub(nbdim,nblk,a(1,1,mpiv),b(1,ipiv))

c------------------------------------------------------------------------
 1000 continue ! next pivot row


c---- back-substitution 
      do 2000 ipiv = namat-1, 1, -1
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks above diagonal
        do 200 je = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(je)) go to 200

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot vector
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo
 200    continue ! next je

 2000 continue

      return
      end ! sbludb



      subroutine sbludbi(nbdim,nblk,namat,ifrst,ilast,mfrst, a, ipp, b)
c-------------------------------------------------------------------------
c     Combines action of SBLUD and SBBKS
c
c   Inputs:
c     nbdim     dimension of blocks
c     nblk      size of blocks
c     namat     size of matrix (number of blocks)
c     ifrst(j)  row index of topmost block in column j
c     ilast(j)  row index of bottommost block in column j
c     mfrst(j)  m index of block matrix element a(..m) at ifrst,j
c     a(..m)    matrix elements stored as blocks in columns
c     b(.i)     r.h.s. vectors
c
c   Outputs:
c     a(..m)    LU-decomposed matrix, to be used in subroutine SBBKS
c     ipp(.m)   pivoting permutation indices
c     b(.i)     Solution vectors
c
c
c   Note:  The block indices in column j are  
c
c            m  =  mfrst(j) ... mfrst(j+1)-1
c
c          So the dummy-block index mfrst(namat+1) must be set accordingly.
c
c-------------------------------------------------------------------------
      integer ifrst(*),ilast(*),mfrst(*)
      real a(nbdim,nbdim,*)
      real b(nbdim,*)
      integer ipp(nbdim,*)

 1200 format(1x, 3(i4,i4,3x) )

c      do j = 1, namat
c        ilast(j) = ifrst(j) + mfrst(j+1) - mfrst(j) - 1
c      enddo

c---- go over all pivot rows
      do 1000 ipiv = 1, namat
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks below diagonal
        do 100 je = 1, jpiv-1
c-------- if this block does not exist, don`t do anything
          if(ipiv .gt. ilast(je)) go to 100

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot block row
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify remaining blocks in pivot block row
          do 10 j = je+1, namat
c---------- if modifying block does not exist, don`t do anything
            if(imod .lt. ifrst(j)) go to 10

            m    = mfrst(j) + ipiv - ifrst(j)
            mmod = mfrst(j) + imod - ifrst(j)

            if(m .gt. mfrst(j+1)-1) then
             write(*,*) '? sblud: no fill space allocated, m > mlast(j)'
             write(*,*) '         ipiv, je, j =', ipiv, je, j
            endif

c---------- BLAS could be used for this matrix-matrix multiply,add
c           A  <--  A - Aeli*Amod
            do jj = 1, nblk
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*a(kk,jj,mmod)
                enddo
                a(ii,jj,m) = a(ii,jj,m) - sum
              enddo
            enddo

 10       continue ! next j
c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - 
 100    continue ! next je (stops at jpiv-1)

c        mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
c        do kk = 1, nblk
c          write(*,'(1x,8f10.5)') (a(kk,ll,mpiv),ll=1,nblk)
c        enddo

c------ LU-decompose diagonal pivot block Apiv
        mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
        call ludcmpi(nbdim,nblk,a(1,1,mpiv),ipp(1,mpiv))

c------------------------------------------------------------------------
c------ multiply remaining blocks in pivot row by Apiv^-1
        do 150 jmul = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(jmul)) go to 150

          mmul = mfrst(jmul) + ipiv - ifrst(jmul)
          call baksubbi(nbdim,nblk,a(1,1,mpiv),ipp(1,mpiv),a(1,1,mmul))
 150    continue

c------------------------------------------------------------------------
c------ multiply r.h.s. vector by Apiv^-1
        call baksubi(nbdim,nblk,a(1,1,mpiv),ipp(1,mpiv),b(1,ipiv))

c------------------------------------------------------------------------
 1000 continue ! next pivot row


c---- back-substitution 
      do 2000 ipiv = namat-1, 1, -1
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks above diagonal
        do 200 je = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(je)) go to 200

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot vector
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo
 200    continue ! next je

 2000 continue

      return
      end ! sbludbi




      subroutine sbbks(nbdim,nblk,namat,ifrst,ilast,mfrst, a, b)
c-------------------------------------------------------------------------
c     Performs fwd,back-substitution using LU-decomposed block matrix 
c     calculated by subroutine SBLUD
c
c   Inputs:
c     nbdim     dimension of blocks
c     nblk      size of blocks
c     namat     size of matrix (number of blocks)
c     ifrst(j)  row index of topmost block in column j
c     ilast(j)  row index of bottommost block in column j
c     mfrst(j)  m index of block matrix element a(..m) at ifrst,j
c     a(..m)    matrix blocks in columns, LU-factored
c     b(.i)     r.h.s. blocked vector
c
c   Outputs:
c     b(.i)     solution vector
c
c-------------------------------------------------------------------------
      integer ifrst(*),ilast(*),mfrst(*)
      real a(nbdim,nbdim,*)
      real b(nbdim,*)

 1200 format(1x, 3(i4,i4,3x) )

c      do j = 1, namat
c        ilast(j) = ifrst(j) + mfrst(j+1) - mfrst(j) - 1
c      enddo

c---- go over all pivot rows
      do 1000 ipiv = 1, namat
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks below diagonal
        do 100 je = 1, jpiv-1
c-------- if this block does not exist, don`t do anything
          if(ipiv .gt. ilast(je)) go to 100

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot vector
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - 
 100    continue ! next je (stops at jpiv-1)

c------------------------------------------------------------------------
c------ multiply r.h.s. vector by Apiv^-1
        mpiv = mfrst(jpiv) + ipiv - ifrst(jpiv)
        call baksub(nbdim,nblk,a(1,1,mpiv),b(1,ipiv))
c------------------------------------------------------------------------
 1000 continue  ! next pivot row


c---- back-substitution 
      do 2000 ipiv = namat-1, 1, -1
        jpiv = ipiv

c------------------------------------------------------------------------
c------ eliminate pivot-row blocks above diagonal
        do 200 je = jpiv+1, namat
c-------- if this block does not exist, don`t do anything
          if(ipiv .lt. ifrst(je)) go to 200

c-------- block ipiv,je to be eliminated
          meli = mfrst(je) + ipiv - ifrst(je)

c-------- block row to be used to modify pivot vector
          imod = je

c- - - - - - - - - - - - - - - - - - - - - - - - - 
c-------- modify r.h.s. vector in pivot block row
            i = ipiv

c---------- BLAS could be used for this matrix-vector multiply,add
c           b  <--  b - Aeli*bmod
              do ii = 1, nblk
                sum = 0.
                do kk = 1, nblk
                  sum = sum + a(ii,kk,meli)*b(kk,imod)
                enddo
                b(ii,i) = b(ii,i) - sum
              enddo
 200    continue ! next je

 2000 continue

      return
      end ! sbbks
