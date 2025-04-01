#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

ifndef SEACAS

default: 
lint:
clean:

else

ifndef NETCDFINC
NETCDFINC = .
endif

ifneq ($(shell unset LD_PRELOAD && printf '\043include <netcdf.h>' | $(CXX) -I$(NETCDFINC) -E -M - 2>&1 > /dev/null && echo 0 || echo 1), 0)
$(error fatal error: 'netcdf.h' file not found. Please set NETCDFINC.)
endif

ifneq ($(wildcard $(SEACAS)/lib/libexodus.*),)
SEACASLIB=$(SEACAS)/lib
endif
ifneq ($(wildcard $(SEACAS)/lib64/libexodus.*),)
SEACASLIB=$(SEACAS)/lib64
endif
ifeq ($(SEACASLIB),)
$(error fatal error: could not find libexodus)
endif

ifdef EXODUS_SIDE_BLOCKS
EXODEF=-DEXODUS_SIDE_BLOCKS
endif

default: $(LDIR)/exodusWriter.so

VPATH = $(ODIR)

OBJS = exodusWriter.o

$(LDIR)/exodusWriter.so:	$(ODIR)/exodusWriter.o $(LDIR)/libaimUtil.a
	$(CXX) $(SOFLGS) $(CPDS) -o $(LDIR)/exodusWriter.so $(ODIR)/exodusWriter.o \
		-L$(SEACASLIB) -lexodus -Wl,-rpath $(SEACASLIB) $(NETCDFLIB) \
		-L$(LDIR) -laimUtil -legads -locsm -ludunits2 -ldl $(RPATH) -lm

$(ODIR)/exodusWriter.o:	exodusWriter.cpp $(IDIR)/aimUtil.h $(IDIR)/aimMesh.h
	$(CXX) -c $(COPTS) $(CPDS) $(DEFINE) $(EXODEF) -I../include -I$(IDIR) -I$(SEACAS)/include -I$(NETCDFINC) exodusWriter.cpp \
		-o $(ODIR)/exodusWriter.o

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:

clean:
	-(cd $(ODIR); rm -f $(OBJS) )
endif

cleanall:	clean
	-rm -f $(LDIR)/exodusWriter.so

