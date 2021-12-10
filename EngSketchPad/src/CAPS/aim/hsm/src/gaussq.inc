
C---- 4-point Gauss quadrature data for quad element
      parameter (ngaussq=4)
      real xsiq(ngaussq), etaq(ngaussq), fwtq(ngaussq)
      data xsiq / -0.577350269189626,  0.577350269189626,
     &             0.577350269189626, -0.577350269189626 /
      data etaq / -0.577350269189626, -0.577350269189626,
     &             0.577350269189626,  0.577350269189626 /
c      data xsiq / -0.8,  0.8,
c     &             0.8, -0.8 /
c      data etaq / -0.8, -0.8,
c     &             0.8,  0.8 /
      data fwtq / 1.0, 1.0,
     &            1.0, 1.0 /

