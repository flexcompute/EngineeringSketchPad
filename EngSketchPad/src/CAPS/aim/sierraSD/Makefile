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
	@echo "Not compiling sierraSDAIM: SEACAS must be set -- Please fix the environment...."
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

OBJS =	$(ODIR)/sierraSDAIM.o

$(LDIR)/sierraSDAIM.so:	$(OBJS) $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/sierraSDAIM.so $(OBJS) -L$(LDIR) \
		-L$(SEACASLIB) -lexodus -Wl,-rpath $(SEACASLIB) $(NETCDFLIB) \
		-lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/sierraSDAIM.o:	sierraSDAIM.c $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../meshWriter/exodusWriter -I../utils  -I$(SEACAS)/include -I$(NETCDFINC) sierraSDAIM.c \
		-o $(ODIR)/sierraSDAIM.o

lint:
	$(LINT) -I$(SEACAS)/include -I../../include -I$(IDIR) -I../utils -I../meshWriter/exodusWriter $(wildcard *.c) -exportlocal -uniondef -nestcomment -boolops

clean:
	-rm -f $(OBJS)
endif

cleanall:	clean
	-rm -f $(LDIR)/sierraSDAIM.so

# Add dox when this AIM is completed
dox: 

todo_dox: doc/html/index.html
doc/html/index.html: \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen sierraSDAIM_dox.cfg; cd latex; make; mv refman.pdf sierraSDAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f sierraSDAIM.tag)
