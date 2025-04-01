#	Makefile for SplineDesPmtr
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
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

default:	$(LDIR)/splinedespmtr.so

$(LDIR)/splinedespmtr.so:	$(ODIR)/udpSplineDesPmtr.o
	$(CXX) $(SOFLGS) -o $(LDIR)/splinedespmtr.so $(ODIR)/udpSplineDesPmtr.o \
		-L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/udpSplineDesPmtr.o:	udpSplineDesPmtr.cpp \
		$(IDIR)/udpUtilities.c $(IDIR)/udpUtilities.h \
		$(IDIR)/egads.h $(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) udpSplineDesPmtr.cpp -o $@

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother

.PHONY: test
test : $(LDIR)/splinedespmtr.so
	$(MAKE) -C test

clean:
	-rm $(ODIR)/udpSplineDesPmtr.o

cleanall:	clean
	-rm $(LDIR)/splinedespmtr.so
