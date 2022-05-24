	program main
c
c      EGADS: Electronic Geometry Aircraft Design System
c
c             Approximate a cylinderical surface
c
c      Copyright 2011-2022, Massachusetts Institute of Technology
c      Licensed under The GNU Lesser General Public License, version 2.1
c      See http://www.opensource.org/licenses/lgpl-2.1.php
c
c
	include 'egads.inc'
c
	pointer   (pxyzs, xyzs), (pivec, ivec), (prvec, rvec)
	integer   sizes(2), oclass, mtype, ivec(*)
	integer*8 context8, cyl8, tess8, bspl8
	real*8    info(13), xyzs(*), rvec(*)
c
c       initialize
c
	istat = IG_open(context8)
	write(*,*) 'IG_open = ', istat
c
c	make our own cylinderical surface
	info( 1) = 0.0
	info( 2) = 0.0
	info( 3) = 0.0
	info( 4) = 1.0
	info( 5) = 0.0
	info( 6) = 0.0
	info( 7) = 0.0
	info( 8) = 1.0
	info( 9) = 0.0
	info(10) = 0.0
	info(11) = 0.0
	info(12) = 1.0
	info(13) = 2.0
	istat = IG_makeGeometry(context8, SURFACE, CYLINDER, NULL,
     &                          0, info, cyl8)
	if(istat .EQ. 0) then
	  info(1)  = 0.0
	  info(2)  = 1.0
	  info(3)  = 0.0
	  info(4)  = 1.0
	  sizes(1) = 160
	  sizes(2) = 16
	  istat = IG_makeTessGeom(cyl8, info, sizes, tess8)
	  if (istat .EQ. 0) then
	    istat = IG_getTessGeom(tess8, sizes, pxyzs)
  	    write(*,*) 'IG_getTessGeom = ', istat, sizes(1), sizes(2)
            istat = IG_approximate(context8, 3, 1.d-6, sizes, xyzs,
     &                             bspl8)
  	    write(*,*) 'IG_approximate = ', istat
	    if (stat .EQ. 0) then
	      istat = IG_getGeometry(bspl8, oclass, mtype, NULL, 
     &                               pivec, prvec)
	      if (istat .EQ. 0) then
	        write(*,1000) oclass, mtype, ivec(1)
 1000	        format(' Oclass = ', i2, '  Mtype = ', i2,
     &                 ' flags = ', z2)
	        write(*,1001) 'U', ivec(2), ivec(3), ivec(4)
	        write(*,1001) 'V', ivec(5), ivec(6), ivec(7)
 1001		format(' ', a1, ' Deg = ', i1, ' #CP = ', i4,
     &                 ' #Knots = ', i3)
	        call IG_free(pivec)
	        call IG_free(prvec)
	      endif
	      istat = IG_deleteObject(bspl8)
  	      write(*,*) 'IG_deleteObject BSpline = ', istat
	    endif
	    istat = IG_deleteObject(tess8)
  	    write(*,*) 'IG_deleteObject tess = ', istat
	  else
	    write(*,*) 'IG_makeTessGeom = ', istat
	  endif
	  istat = IG_deleteObject(cyl8)
  	  write(*,*) 'IG_deleteObject cyl = ', istat
	else
	  write(*,*) 'IG_makeGeometry = ', istat
	endif
c
	istat = IG_close(context8)
	write(*,*) 'IG_close = ', istat
c
	call exit
	end
