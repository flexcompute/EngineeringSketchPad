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

$(TDIR)/fit:	$(ODIR)/fit.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/fit $(ODIR)/fit.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/fit.o:	fit.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) fit.c -o $(ODIR)/fit.o

clean:
	-rm $(ODIR)/fit.o 

cleanall:	clean
	-rm $(TDIR)/fit
