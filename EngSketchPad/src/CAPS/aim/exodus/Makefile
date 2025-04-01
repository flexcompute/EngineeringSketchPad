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
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
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

OBJS =	$(ODIR)/exodusAIM.o $(ODIR)/libmeshb7.o

$(LDIR)/exodusAIM.so:	$(OBJS)  $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/exodusAIM.so $(OBJS) -L$(LDIR) \
		-L$(SEACASLIB) -lexodus -Wl,-rpath $(SEACASLIB) $(NETCDFLIB) \
		-laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/exodusAIM.o:	exodusAIM.c $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) \
		-I../meshWriter/exodusWriter \
		-I../meshWriter/libmeshbWriter \
		-I../utils -I$(SEACAS)/include \
		-I$(NETCDFINC) \
		exodusAIM.c \
		-o $(ODIR)/exodusAIM.o

$(ODIR)/libmeshb7.o:	../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
			../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.h
	$(CC) -c $(COPTS) $(DEFINE) \
		-I../meshWriter/libmeshbWriter/libMeshb/sources \
		../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
		-o $(ODIR)/libmeshb7.o

lint:
	$(LINT) -I$(SEACAS)/include -I../../include -I$(IDIR) -I../utils -I../meshWriter/exodusWriter -I../meshWriter/libmeshbWriter $(wildcard *.c) -exportlocal -uniondef -nestcomment -boolops

clean:
	-rm -f $(OBJS)
endif

cleanall:	clean
	-rm -f $(LDIR)/exodusAIM.so

dox: doc/html/index.html
doc/html/index.html: \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen exodusAIM_dox.cfg; cd latex; make; mv refman.pdf exodusAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f exodusAIM.tag)
