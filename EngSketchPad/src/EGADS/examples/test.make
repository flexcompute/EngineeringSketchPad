#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

default:	$(TDIR)/testc $(TDIR)/testf $(TDIR)/parsec $(TDIR)/parsef

$(TDIR)/testc:		$(ODIR)/testc.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/testc $(ODIR)/testc.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/testc.o:	testc.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) testc.c -o $(ODIR)/testc.o

$(TDIR)/testf:		$(ODIR)/testf.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/testf $(ODIR)/testf.o -L$(LDIR) -lfgads -legads \
		$(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/testf.o:	testf.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) testf.f -o $(ODIR)/testf.o

$(TDIR)/parsec:		$(ODIR)/parsec.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/parsec $(ODIR)/parsec.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/parsec.o:	parsec.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) parsec.c -o $(ODIR)/parsec.o

$(TDIR)/parsef:		$(ODIR)/parsef.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/parsef $(ODIR)/parsef.o -L$(LDIR) -lfgads -legads \
		$(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/parsef.o:	parsef.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) parsef.f -o $(ODIR)/parsef.o

clean:
	-rm $(ODIR)/testc.o $(ODIR)/testf.o $(ODIR)/parsec.o $(ODIR)/parsef.o

cleanall:	clean
	-rm $(TDIR)/testc   $(TDIR)/testf   $(TDIR)/parsec   $(TDIR)/parsef
