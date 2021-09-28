
      subroutine bmdump(lu,
     &                 irtot,kdim,nmdim,
     &                 neq,nnode,nmmat, ifrst,ilast,mfrst,amat)
c----------------------------------------------------------------
c     Writes matrix structure info to logical unit lu,
c     for use in matrix-plotting program bmplot
c----------------------------------------------------------------

c------------------------------------
c---- passed arrays
      integer ifrst(kdim),
     &        ilast(kdim),
     &        mfrst(kdim+1)
      real amat(irtot,irtot,nmdim)

      open(lu,status='unknown',form='unformatted')

      write(lu) neq, nnode, nmmat

      do j = 1, nnode
        ilast(j) = ifrst(j) + mfrst(j+1) - mfrst(j) - 1
        write(lu) ifrst(j), ilast(j)
      enddo

      do j = 1, nnode
        do i = ifrst(j), ilast(j)
          m = mfrst(j) - ifrst(j) + i

          do l = 1, neq
            do k = 1, neq
              if(amat(k,l,m) .ne. 0.0) write(lu) i, j, k, l
            enddo
          enddo

        enddo
      enddo
      write(lu) 0, 0, 0, 0
c
      close(lu)

      return
      end ! bmdump

