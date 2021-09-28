      subroutine cross(a,b,c)
c-------------------------------------
c     Returns cross product c = a x b
c-------------------------------------
      real a(3), b(3), c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end ! cross
