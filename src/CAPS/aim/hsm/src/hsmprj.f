      
      subroutine hsmfpr(bf1, bf2, bfn, bm1, bm2, bmn,
     &                  resc,
     &                  resp,
     &                  resp_bf1,
     &                  resp_bf2,
     &                  resp_bfn,
     &                  resp_bm1,
     &                  resp_bm2,
     &                  resp_bmn )
c--------------------------------------------------------------------
c     Projects HSM force residuals and their Jacobians 
c     from xyz basis into b1,b2,bn basis.
c
c     Assumes resc(.) contains the following residuals:
c
c      res( 1) = R^fx |     | resp(1) = R^f1
c      res( 2) = R^fy | --> | resp(2) = R^f2
c      res( 3) = R^fz |     | resp(3) = R^fn
c
c      res( 4) = R^mx |     | resp(4) = R^m1
c      res( 5) = R^my | --> | resp(5) = R^m2
c      res( 6) = R^mz |     | resp(6) = R^mn
c
c
c    Inputs:
c      bf1(.)  force projection vectors
c      bf2(.)
c      bfn(.)
c      bm1(.)  moment projection vectors
c      bm2(.)
c      bmn(.)
c      resc(.)  cartesian xyz residuals
c
c    Outputs:
c      resp(.)  projected b-basis residuals
c
c--------------------------------------------------------------------
      include 'index.inc'

      real bf1(3), bf2(3), bfn(3), bm1(3), bm2(3), bmn(3)
      real resc(irtot)
      real resp(irtot),
     &     resp_bf1(irtot,3),
     &     resp_bf2(irtot,3),
     &     resp_bfn(irtot,3),
     &     resp_bm1(irtot,3),
     &     resp_bm2(irtot,3),
     &     resp_bmn(irtot,3)

      do ir = 1, irtot
        do k = 1, 3
          resp_bf1(ir,k) = 0.
          resp_bf2(ir,k) = 0.
          resp_bfn(ir,k) = 0.
          resp_bm1(ir,k) = 0.
          resp_bm2(ir,k) = 0.
          resp_bmn(ir,k) = 0.
        enddo
      enddo

c---- R^f1
      ir = irf1
      resp(ir) = bf1(1)*resc(1)
     &         + bf1(2)*resc(2)
     &         + bf1(3)*resc(3)
      resp_bf1(ir,1) =  resc(1)
      resp_bf1(ir,2) =  resc(2)
      resp_bf1(ir,3) =  resc(3)

c---- R^f2
      ir = irf2
      resp(ir) = bf2(1)*resc(1)
     &         + bf2(2)*resc(2)
     &         + bf2(3)*resc(3)
      resp_bf2(ir,1) =  resc(1)
      resp_bf2(ir,2) =  resc(2)
      resp_bf2(ir,3) =  resc(3)

c---- R^fn
      ir = irfn
      resp(ir) = bfn(1)*resc(1)
     &         + bfn(2)*resc(2)
     &         + bfn(3)*resc(3)
      resp_bfn(ir,1) =  resc(1)
      resp_bfn(ir,2) =  resc(2)
      resp_bfn(ir,3) =  resc(3)

c---- R^m1
      ir = irm1
      resp(ir) = bm2(1)*resc(4)
     &         + bm2(2)*resc(5)
     &         + bm2(3)*resc(6)
      resp_bm2(ir,1) =  resc(4)
      resp_bm2(ir,2) =  resc(5)
      resp_bm2(ir,3) =  resc(6)

c---- R^m2
      ir = irm2
      resp(ir) = bm1(1)*resc(4)
     &         + bm1(2)*resc(5)
     &         + bm1(3)*resc(6)
      resp_bm1(ir,1) =  resc(4)
      resp_bm1(ir,2) =  resc(5)
      resp_bm1(ir,3) =  resc(6)

cc---- R^m1
c      ir = irm1
c      resp(ir) = bm1(1)*resc(4)
c     &         + bm1(2)*resc(5)
c     &         + bm1(3)*resc(6)
c      resp_bm1(ir,1) =  resc(4)
c      resp_bm1(ir,2) =  resc(5)
c      resp_bm1(ir,3) =  resc(6)

cc---- R^m2
c      ir = irm2
c      resp(ir) = bm2(1)*resc(4)
c     &         + bm2(2)*resc(5)
c     &         + bm2(3)*resc(6)
c      resp_bm2(ir,1) =  resc(4)
c      resp_bm2(ir,2) =  resc(5)
c      resp_bm2(ir,3) =  resc(6)

c---- R^mn
      ir = irmn
      resp(ir) = bmn(1)*resc(4)
     &         + bmn(2)*resc(5)
     &         + bmn(3)*resc(6)
      resp_bmn(ir,1) =  resc(4)
      resp_bmn(ir,2) =  resc(5)
      resp_bmn(ir,3) =  resc(6)

      return
      end ! hsmfpr

      
      subroutine hsmdpr(bf1, bf2, bfn, bm1, bm2,
     &                  res,
     &                  resp, 
     &                  resp_bf1,
     &                  resp_bf2,
     &                  resp_bfn,
     &                  resp_bm1,
     &                  resp_bm2 )
c--------------------------------------------------------------------
c     Projects HSM Dirichlet residuals and their Jacobians 
c     from xyz basis into bf,bm bases.
c
c     Assumes res(.) contains the following residuals:
c
c      res( 1) = R^rx |     | resp(1) = R^r1
c      res( 2) = R^ry | --> | resp(2) = R^r2
c      res( 3) = R^rz |     | resp(3) = R^rn
c
c      res( 4) = R^dx |     | resp(4) = R^d1
c      res( 5) = R^dy | --> | resp(5) = R^d2
c      res( 6) = R^dz |     | resp(6) = R^psi
c
c
c    Inputs:
c      bf1(.)   force projection vectors
c      bf1(.)
c      bfn(.)
c      bm1(.)   moment projection vectors
c      bm2(.)
c      res(.)  xyz residuals
c
c    Outputs:
c      resp(.)  projected residuals
c
c--------------------------------------------------------------------
      include 'index.inc'

      real bf1(3),
     &     bf2(3),
     &     bfn(3),
     &     bm1(3),
     &     bm2(3)
      real res(ivtot)
      real resp(irtot),
     &     resp_bf1(irtot,3),
     &     resp_bf2(irtot,3),
     &     resp_bfn(irtot,3),
     &     resp_bm1(irtot,3),
     &     resp_bm2(irtot,3)

      do k = 1, 3
        do ir = 1, irtot
          resp_bf1(ir,k) = 0.
          resp_bf2(ir,k) = 0.
          resp_bfn(ir,k) = 0.
          resp_bm1(ir,k) = 0.
          resp_bm2(ir,k) = 0.
        enddo
      enddo

c---- R^r1
      ir = irf1
      resp(ir) = bf1(1)*res(1)
     &         + bf1(2)*res(2)
     &         + bf1(3)*res(3)
      resp_bf1(ir,1) = res(1)
      resp_bf1(ir,2) = res(2)
      resp_bf1(ir,3) = res(3)

c---- R^r2
      ir = irf2
      resp(ir) = bf2(1)*res(1)
     &         + bf2(2)*res(2)
     &         + bf2(3)*res(3)
      resp_bf2(ir,1) = res(1)
      resp_bf2(ir,2) = res(2)
      resp_bf2(ir,3) = res(3)

c---- R^rn
      ir = irfn
      resp(ir) = bfn(1)*res(1)
     &         + bfn(2)*res(2)
     &         + bfn(3)*res(3)
      resp_bfn(ir,1) = res(1)
      resp_bfn(ir,2) = res(2)
      resp_bfn(ir,3) = res(3)

      
c---- R^d1
      ir = irm1
      resp(ir) = bm1(1)*res(4)
     &         + bm1(2)*res(5)
     &         + bm1(3)*res(6)
      resp_bm1(ir,1) = res(4)
      resp_bm1(ir,2) = res(5)
      resp_bm1(ir,3) = res(6)

c---- R^d2
      ir = irm2
      resp(ir) = bm2(1)*res(4)
     &         + bm2(2)*res(5)
     &         + bm2(3)*res(6)
      resp_bm2(ir,1) = res(4)
      resp_bm2(ir,2) = res(5)
      resp_bm2(ir,3) = res(6)

c---- R^psi
      ir = irmn
      resp(ir) = res(7)
      
      end ! hsmdpr



      subroutine hsmvpr(bf1, bf2, bfn, bm1, bm2,
     &                  resp_var, 
     &                  resp_dvp )
c--------------------------------------------------------------------
c     Changes HSM Jacobians from primary variables in xyz basis
c     into projected variables in b1,b2,bn basis.
c
c     Primary variables from projected variables:
c      dr   = bf1 dr1  +  bf2 dr2  +  bfn drn
c      dd   = bm1 dn1  +  bm2 dn2
c      dpsi =     dpsi
c
c     Projected variable array indexing:
c      dvp(1) = dr1
c      dvp(2) = dr2
c      dvp(3) = drn
c      dvp(4) = dd1
c      dvp(5) = dd2
c      dvp(6) = dpsi
c
c    Inputs:
c      bf1(.)   force projection vectors
c      bf1(.)
c      bfn(.)
c      bm1(.)   moment projection vectors
c      bm2(.)
c      resp_var(...)  dresp(.)/dvar(..)
c
c    Outputs:
c      resp_dvp(...)  dresp(.)/dvp(..)
c
c--------------------------------------------------------------------
      include 'index.inc'

      real bf1(3), bf2(3), bfn(3), bm1(3), bm2(3)
      real resp_var(irtot,ivtot)
      real resp_dvp(irtot,irtot)

c---- project variables for each equation
      do ir = 1, irtot

c         write(*,*) 'ir =', ir, bf1
c         write(*,*) (resp_var(ir,iv), iv=1, ivtot)
         

c------ dR / dr1
        resp_dvp(ir,irf1) = resp_var(ir,ivrx)*bf1(1)
     &                    + resp_var(ir,ivry)*bf1(2)
     &                    + resp_var(ir,ivrz)*bf1(3)

c------ dR / dr2
        resp_dvp(ir,irf2) = resp_var(ir,ivrx)*bf2(1)
     &                    + resp_var(ir,ivry)*bf2(2)
     &                    + resp_var(ir,ivrz)*bf2(3)

c------ dR / drn
        resp_dvp(ir,irfn) = resp_var(ir,ivrx)*bfn(1)
     &                    + resp_var(ir,ivry)*bfn(2)
     &                    + resp_var(ir,ivrz)*bfn(3)

c------ dR / dd1
        resp_dvp(ir,irm1) = resp_var(ir,ivdx)*bm1(1)
     &                    + resp_var(ir,ivdy)*bm1(2)
     &                    + resp_var(ir,ivdz)*bm1(3)

c------ dR / dd2
        resp_dvp(ir,irm2) = resp_var(ir,ivdx)*bm2(1)
     &                    + resp_var(ir,ivdy)*bm2(2)
     &                    + resp_var(ir,ivdz)*bm2(3)

c------ dR / dpsi
        resp_dvp(ir,irmn) = resp_var(ir,ivps)

      enddo

      return
      end ! hsmvpr


      subroutine hsmvpr1(ir,
     &                  bf1, bf2, bfn, bm1, bm2,
     &                  resp_var, 
     &                  resp_dvp )
c--------------------------------------------------------------------
c     Same as hsmvpr, but only processes one equation ir
c--------------------------------------------------------------------
      include 'index.inc'

      real bf1(3), bf2(3), bfn(3), bm1(3), bm2(3)
      real resp_var(irtot,ivtot)
      real resp_dvp(irtot,irtot)

c------ dR / dr1
        resp_dvp(ir,irf1) = resp_var(ir,ivrx)*bf1(1)
     &                    + resp_var(ir,ivry)*bf1(2)
     &                    + resp_var(ir,ivrz)*bf1(3)

c------ dR / dr2
        resp_dvp(ir,irf2) = resp_var(ir,ivrx)*bf2(1)
     &                    + resp_var(ir,ivry)*bf2(2)
     &                    + resp_var(ir,ivrz)*bf2(3)

c------ dR / drn
        resp_dvp(ir,irfn) = resp_var(ir,ivrx)*bfn(1)
     &                    + resp_var(ir,ivry)*bfn(2)
     &                    + resp_var(ir,ivrz)*bfn(3)

c------ dR / dd1
        resp_dvp(ir,irm1) = resp_var(ir,ivdx)*bm1(1)
     &                    + resp_var(ir,ivdy)*bm1(2)
     &                    + resp_var(ir,ivdz)*bm1(3)

c------ dR / dd2
        resp_dvp(ir,irm2) = resp_var(ir,ivdx)*bm2(1)
     &                    + resp_var(ir,ivdy)*bm2(2)
     &                    + resp_var(ir,ivdz)*bm2(3)

c------ dR / dps
        resp_dvp(ir,irmn) = resp_var(ir,ivps)

      return
      end ! hsmvpr1
