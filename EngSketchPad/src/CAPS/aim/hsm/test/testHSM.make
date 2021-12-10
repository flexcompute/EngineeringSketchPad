#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = ../src
TDIR  = $(ESP_ROOT)/bin
endif

all: # Do nothing for now...

default:    $(TDIR)/testHSM
#default:   $(TDIR)/testHSM $(TDIR)/hsmrun

#$(TDIR)/hsmrun:     $(ODIR)/hsmrun.o $(LDIR)/libhsm.a
#	$(FCOMP) -o $(TDIR)/hsmrun $(ODIR)/hsmrun.o -L$(LDIR) -lhsm

#$(ODIR)/hsmrun.o:   hsmrun.f
#	$(FCOMP) -c $(FOPTS) $(FFLAG) -fno-recursive hsmrun.f -o $(ODIR)/hsmrun.o

$(TDIR)/testHSM:	$(ODIR)/testHSM.o $(ODIR)/valence.o \
			$(LDIR)/librcm.a $(LDIR)/libhsm.a
#$(LDIR)/$(SHLIB)

#	$(CXX) -o $(TDIR)/testHSM $(ODIR)/testHSM.o \
#       $(ODIR)/valence.o -L$(LDIR) -lrcm -lhsm -legads \
#       -lgfortran $(RPATH) -Wl,-no_compact_unwind
ifeq ($(EFCOMP),ifort)
	$(FCOMP) -o $(TDIR)/testHSM $(ODIR)/testHSM.o \
		$(ODIR)/valence.o -L$(LDIR) \
		-lrcm -lhsm -legads $(FRPATH) $(CPPSLB) $(ESPFLIBS) -nofor-main
else
	$(FCOMP) -o $(TDIR)/testHSM $(ODIR)/testHSM.o \
		$(ODIR)/valence.o -L$(LDIR) \
		-lrcm -lhsm -legads $(FRPATH) $(CPPSLB) $(ESPFLIBS)
endif

$(LDIR)/libhsm.a:
	$(MAKE) -C ../src

$(LDIR)/librcm.a:
	$(MAKE) -C ../rcm

$(ODIR)/testHSM.o:	testHSM.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) testHSM.c -o $(ODIR)/testHSM.o

$(ODIR)/valence.o:	valence.c
	$(CC) -c $(COPTS) $(DEFINE) valence.c -o $(ODIR)/valence.o

clean:
	-rm $(ODIR)/testHSM.o $(ODIR)/valence.o

cleanall:   clean
	-rm $(TDIR)/testHSM
