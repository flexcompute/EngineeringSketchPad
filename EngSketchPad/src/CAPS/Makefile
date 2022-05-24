#
IDIR = $(ESP_ROOT)/include
LDIR = $(ESP_ROOT)/lib
BDIR = $(ESP_ROOT)/bin

ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libudunits2.dylib $(IDIR)/udunits2.h \
		$(LDIR)/libc3dio.a $(LDIR)/libCart3D.a \
		$(LDIR)/libgcc_s.1.dylib $(LDIR)/libquadmath.0.dylib \
		$(BDIR)/avl $(BDIR)/pplot $(BDIR)/pxplot $(BDIR)/xfoil \
		$(BDIR)/mastros.exe $(BDIR)/ASTRO.D01 $(BDIR)/ASTRO.IDX \
		$(BDIR)/mses $(BDIR)/mset

$(LDIR)/libudunits2.dylib:	udunits/DARWIN64/libudunits2.0.dylib
	/bin/cp udunits/DARWIN64/libudunits2.0.dylib $(LDIR)
	/bin/ln -s $(LDIR)/libudunits2.0.dylib $(LDIR)/libudunits2.dylib

$(IDIR)/udunits2.h:		udunits/include/udunits2.h
	/bin/cp udunits/include/*.h $(IDIR)

$(LDIR)/libc3dio.a:		Cart3D/OSX64_CLANG/libc3dio.a
	/bin/cp Cart3D/OSX64_CLANG/libc3dio.a $(LDIR)/libc3dio.a

$(LDIR)/libCart3D.a:		Cart3D/OSX64_CLANG/libCart3D.a
	/bin/cp Cart3D/OSX64_CLANG/libCart3D.a $(LDIR)/libCart3D.a

$(LDIR)/libgcc_s.1.dylib:	Libraries/DARWIN64/libgcc_s.1.dylib
	/bin/cp Libraries/DARWIN64/libgcc_s.1.dylib $(LDIR)/libgcc_s.1.dylib

$(LDIR)/libquadmath.0.dylib:	Libraries/DARWIN64/libquadmath.0.dylib
	/bin/cp Libraries/DARWIN64/libquadmath.0.dylib $(LDIR)/libquadmath.0.dylib

$(BDIR)/avl:			Executables/DARWIN64/avl
	/bin/cp Executables/DARWIN64/avl $(BDIR)/avl

$(BDIR)/pplot:			Executables/DARWIN64/pplot
	/bin/cp Executables/DARWIN64/pplot $(BDIR)/pplot

$(BDIR)/pxplot:			Executables/DARWIN64/pxplot
	/bin/cp Executables/DARWIN64/pxplot $(BDIR)/pxplot

$(BDIR)/xfoil:			Executables/DARWIN64/xfoil
	/bin/cp Executables/DARWIN64/xfoil $(BDIR)/xfoil

$(BDIR)/ASTRO.D01:		Executables/DARWIN64/ASTRO.D01
	/bin/cp Executables/DARWIN64/ASTRO.D01 $(BDIR)/ASTRO.D01

$(BDIR)/ASTRO.IDX:		Executables/DARWIN64/ASTRO.IDX
	/bin/cp Executables/DARWIN64/ASTRO.IDX $(BDIR)/ASTRO.IDX

$(BDIR)/mastros.exe:		Executables/DARWIN64/mastros.exe
	/bin/cp Executables/DARWIN64/mastros.exe $(BDIR)/mastros.exe
	install_name_tool -add_rpath $(LDIR) $(BDIR)/mastros.exe

$(BDIR)/mses:			Executables/DARWIN64/mses
	/bin/cp Executables/DARWIN64/mses $(BDIR)/mses

$(BDIR)/mset:			Executables/DARWIN64/mset
	/bin/cp Executables/DARWIN64/mset $(BDIR)/mset
endif

ifeq ($(ESP_ARCH),DARWINM1)
default:	$(LDIR)/libudunits2.dylib $(IDIR)/udunits2.h \
		$(LDIR)/libc3dio.a $(LDIR)/libCart3D.a \
		$(LDIR)/libgcc_s.1.dylib $(LDIR)/libquadmath.0.dylib \
		$(BDIR)/avl $(BDIR)/pplot $(BDIR)/pxplot $(BDIR)/xfoil \
		$(BDIR)/mastros.exe $(BDIR)/ASTRO.D01 $(BDIR)/ASTRO.IDX \
		$(BDIR)/mses $(BDIR)/mset

$(LDIR)/libudunits2.dylib:	udunits/DARWINM1/libudunits2.0.dylib
	/bin/cp udunits/DARWINM1/libudunits2.0.dylib $(LDIR)
	/bin/ln -s $(LDIR)/libudunits2.0.dylib $(LDIR)/libudunits2.dylib

$(IDIR)/udunits2.h:		udunits/include/udunits2.h
	/bin/cp udunits/include/*.h $(IDIR)

$(LDIR)/libc3dio.a:		Cart3D/M1_CLANG/libc3dio.a
	/bin/cp Cart3D/M1_CLANG/libc3dio.a $(LDIR)/libc3dio.a

$(LDIR)/libCart3D.a:		Cart3D/M1_CLANG/libCart3D.a
	/bin/cp Cart3D/M1_CLANG/libCart3D.a $(LDIR)/libCart3D.a

$(LDIR)/libgcc_s.1.dylib:	Libraries/DARWIN64/libgcc_s.1.dylib
	/bin/cp Libraries/DARWIN64/libgcc_s.1.dylib $(LDIR)/libgcc_s.1.dylib

$(LDIR)/libquadmath.0.dylib:	Libraries/DARWIN64/libquadmath.0.dylib
	/bin/cp Libraries/DARWIN64/libquadmath.0.dylib $(LDIR)/libquadmath.0.dylib

$(BDIR)/avl:			Executables/DARWIN64/avl
	/bin/cp Executables/DARWIN64/avl $(BDIR)/avl

$(BDIR)/pplot:			Executables/DARWIN64/pplot
	/bin/cp Executables/DARWIN64/pplot $(BDIR)/pplot

$(BDIR)/pxplot:			Executables/DARWIN64/pxplot
	/bin/cp Executables/DARWIN64/pxplot $(BDIR)/pxplot

$(BDIR)/xfoil:			Executables/DARWIN64/xfoil
	/bin/cp Executables/DARWIN64/xfoil $(BDIR)/xfoil

$(BDIR)/ASTRO.D01:		Executables/DARWIN64/ASTRO.D01
	/bin/cp Executables/DARWIN64/ASTRO.D01 $(BDIR)/ASTRO.D01

$(BDIR)/ASTRO.IDX:		Executables/DARWIN64/ASTRO.IDX
	/bin/cp Executables/DARWIN64/ASTRO.IDX $(BDIR)/ASTRO.IDX

$(BDIR)/mastros.exe:		Executables/DARWIN64/mastros.exe
	/bin/cp Executables/DARWIN64/mastros.exe $(BDIR)/mastros.exe
	install_name_tool -add_rpath $(LDIR) $(BDIR)/mastros.exe

$(BDIR)/mses:			Executables/DARWINM1/mses
	/bin/cp Executables/DARWINM1/mses $(BDIR)/mses

$(BDIR)/mset:			Executables/DARWINM1/mset
	/bin/cp Executables/DARWINM1/mset $(BDIR)/mset
endif

ifeq ($(ESP_ARCH),LINUX64)
UDUNITS=FOUND
  ifeq ("$(wildcard /usr/include/udunits*.h)","")
    ifeq ("$(wildcard /usr/include/udunits2/udunits*.h)","")
#     $(error UDUNITS2 not found -- Please load the package...)
      UDUNITS=NOTFOUND
    endif
  endif

ifeq ($(UDUNITS),FOUND)
default:	$(LDIR)/libc3dio.a $(LDIR)/libCart3D.a \
		$(BDIR)/avl $(BDIR)/pplot $(BDIR)/pxplot $(BDIR)/xfoil \
		$(BDIR)/mastros.exe $(BDIR)/ASTRO.D01 $(BDIR)/ASTRO.IDX \
		$(BDIR)/mses $(BDIR)/mset
else
default:	$(LDIR)/libudunits2.so $(IDIR)/udunits2.h \
		$(LDIR)/libc3dio.a $(LDIR)/libCart3D.a \
		$(BDIR)/avl $(BDIR)/pplot $(BDIR)/pxplot $(BDIR)/xfoil \
		$(BDIR)/mastros.exe $(BDIR)/ASTRO.D01 $(BDIR)/ASTRO.IDX \
		$(BDIR)/mses $(BDIR)/mset
endif

$(LDIR)/libudunits2.so:	udunits/LINUX64/libudunits2.so
	/bin/cp udunits/LINUX64/libudunits2.so $(LDIR)
	ln -s $(LDIR)/libudunits2.so $(LDIR)/libudunits2.so.0

$(IDIR)/udunits2.h:	udunits/include/udunits2.h
	/bin/cp udunits/include/*.h $(IDIR)

$(LDIR)/libc3dio.a:	Cart3D/X86_64_GCC/libc3dio.a
	/bin/cp Cart3D/X86_64_GCC/libc3dio.a $(LDIR)/libc3dio.a

$(LDIR)/libCart3D.a:	Cart3D/X86_64_GCC/libCart3D.a
	/bin/cp Cart3D/X86_64_GCC/libCart3D.a $(LDIR)/libCart3D.a

$(BDIR)/avl:		Executables/LINUX64/avl
	/bin/cp Executables/LINUX64/avl $(BDIR)/avl

$(BDIR)/pplot:		Executables/LINUX64/pplot
	/bin/cp Executables/LINUX64/pplot $(BDIR)/pplot

$(BDIR)/pxplot:		Executables/LINUX64/pxplot
	/bin/cp Executables/LINUX64/pxplot $(BDIR)/pxplot

$(BDIR)/xfoil:		Executables/LINUX64/xfoil
	/bin/cp Executables/LINUX64/xfoil $(BDIR)/xfoil

$(BDIR)/ASTRO.D01:	Executables/LINUX64/ASTRO.D01
	/bin/cp Executables/LINUX64/ASTRO.D01 $(BDIR)/ASTRO.D01

$(BDIR)/ASTRO.IDX:	Executables/LINUX64/ASTRO.IDX
	/bin/cp Executables/LINUX64/ASTRO.IDX $(BDIR)/ASTRO.IDX

$(BDIR)/mastros.exe:	Executables/LINUX64/mastros.exe
	/bin/cp Executables/LINUX64/mastros.exe $(BDIR)/mastros.exe

$(BDIR)/mses:		Executables/LINUX64/mses
	/bin/cp Executables/LINUX64/mses $(BDIR)/mses

$(BDIR)/mset:		Executables/LINUX64/mset
	/bin/cp Executables/LINUX64/mset $(BDIR)/mset
endif
