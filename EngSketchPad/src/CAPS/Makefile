#
IDIR = $(ESP_ROOT)/include
LDIR = $(ESP_ROOT)/lib
BDIR = $(ESP_ROOT)/bin

DEFAULTS = $(BDIR)/avl $(BDIR)/pplot $(BDIR)/pxplot $(BDIR)/xfoil \
	   $(BDIR)/mastros.exe $(BDIR)/ASTRO.D01 $(BDIR)/ASTRO.IDX \
	   $(BDIR)/mses $(BDIR)/mset

# Check for system udunits
UDUNITS=FOUND
ifeq ("$(wildcard /usr/include/udunits*.h)","")
  ifeq ("$(wildcard /usr/include/udunits2/udunits*.h)","")
#     $(error UDUNITS2 not found -- Please load the package...)
    UDUNITS=NOTFOUND
  endif
endif

ifeq ($(UDUNITS),NOTFOUND)
  DEFAULTS += $(IDIR)/udunits2.h
  ifneq (,$(findstring LINUX,$(ESP_ARCH)))
    DEFAULTS += $(LDIR)/libudunits2.so
  endif
  ifneq (,$(findstring DARWIN,$(ESP_ARCH)))
    DEFAULTS += $(LDIR)/libudunits2.dylib
  endif
endif

ifndef AFLR
DEFAULTS +=  $(LDIR)/aflr2AIM.so $(LDIR)/aflr3AIM.so $(LDIR)/aflr4AIM.so
endif

default:	$(DEFAULTS)

ifeq ($(ESP_ARCH),LINUX_AARCH64)
$(LDIR)/libc3dio.a:

$(LDIR)/libCart3D.a:
endif

ifeq ($(UDUNITS),NOTFOUND)
$(IDIR)/udunits2.h:	udunits/include/udunits2.h
	/bin/cp udunits/include/*.h $(IDIR)

# udunits DARWIN
$(LDIR)/libudunits2.dylib:	udunits/$(ESP_ARCH)/libudunits2.0.dylib
	/bin/cp udunits/$(ESP_ARCH)/libudunits2.0.dylib $(LDIR)
	/bin/ln -sf libudunits2.0.dylib $(LDIR)/libudunits2.dylib

# udunits LINUX
$(LDIR)/libudunits2.so:	udunits/$(ESP_ARCH)/libudunits2.so
	/bin/cp udunits/$(ESP_ARCH)/libudunits2.so $(LDIR)
	ln -sf libudunits2.so $(LDIR)/libudunits2.so.0
endif

$(BDIR)/avl:		Executables/$(ESP_ARCH)/avl
	/bin/cp Executables/$(ESP_ARCH)/avl $(BDIR)/avl

$(BDIR)/pplot:		Executables/$(ESP_ARCH)/pplot
	/bin/cp Executables/$(ESP_ARCH)/pplot $(BDIR)/pplot

$(BDIR)/pxplot:		Executables/$(ESP_ARCH)/pxplot
	/bin/cp Executables/$(ESP_ARCH)/pxplot $(BDIR)/pxplot

$(BDIR)/xfoil:		Executables/$(ESP_ARCH)/xfoil
	/bin/cp Executables/$(ESP_ARCH)/xfoil $(BDIR)/xfoil

$(BDIR)/ASTRO.D01:	Executables/$(ESP_ARCH)/ASTRO.D01
	/bin/cp Executables/$(ESP_ARCH)/ASTRO.D01 $(BDIR)/ASTRO.D01

$(BDIR)/ASTRO.IDX:	Executables/$(ESP_ARCH)/ASTRO.IDX
	/bin/cp Executables/$(ESP_ARCH)/ASTRO.IDX $(BDIR)/ASTRO.IDX

$(BDIR)/mastros.exe:	Executables/$(ESP_ARCH)/mastros.exe
	/bin/cp Executables/$(ESP_ARCH)/mastros.exe $(BDIR)/mastros.exe

$(BDIR)/mses:		Executables/$(ESP_ARCH)/mses
	/bin/cp Executables/$(ESP_ARCH)/mses $(BDIR)/mses

$(BDIR)/mset:		Executables/$(ESP_ARCH)/mset
	/bin/cp Executables/$(ESP_ARCH)/mset $(BDIR)/mset

ifndef AFLR
$(LDIR)/aflr2AIM.so:		Libraries/$(ESP_ARCH)/aflr2AIM.so
	/bin/cp Libraries/$(ESP_ARCH)/aflr2AIM.so $(LDIR)/aflr2AIM.so

$(LDIR)/aflr3AIM.so:		Libraries/$(ESP_ARCH)/aflr3AIM.so
	/bin/cp Libraries/$(ESP_ARCH)/aflr3AIM.so $(LDIR)/aflr3AIM.so

$(LDIR)/aflr4AIM.so:		Libraries/$(ESP_ARCH)/aflr4AIM.so
	/bin/cp Libraries/$(ESP_ARCH)/aflr4AIM.so $(LDIR)/aflr4AIM.so
endif
