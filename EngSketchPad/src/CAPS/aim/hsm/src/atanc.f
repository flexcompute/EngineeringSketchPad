
      function atanc(y,x,thold)
c---------------------------------------------------------------
c     atan2 function with branch cut checking.
c
c     Increments position angle of point x,y from some previous
c     value thold due to a change in position, ensuring that the
c     position change does not cross the atan2 branch cut
c     (which is in the -x direction).  For example:
c
c       atanc( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
c       atan2( -1.0 , -1.0 )            returns  -.75*pi .
c
c     Typically, atanc is used to fill an array of angles:
c
c        theta(1) = atan2( y(1) , x(1) )
c        do i=2, n
c          theta(i) = atanc( y(i) , x(i) , theta(i-1) )
c        end do
c
c     This will prevent the angle array theta(i) from jumping by 
c     +/- 2 pi when the path x(i),y(i) crosses the negative x axis.
c
c     Input:
c       x,y     point position coordinates
c       thold   position angle of nearby point
c
c     Output:
c       atanc   position angle of x,y
c---------------------------------------------------------------
      data  pi / 3.14159265358979323846264338327950 /
      data tpi / 6.28318530717958647692528676655901 /
c
c---- set new position angle, ignoring branch cut in atan2 function for now
      thnew = atan2( y , x )
c
c---- set angle change and remove any multiples of +/-2 pi 
      dt = thnew - thold
      dtcorr = dt - sign(tpi,dt)*aint( abs(dt)/tpi + 0.5 )
c
c---- set correct new angle
      atanc = thold + dtcorr
c
      return
      end ! atanc

