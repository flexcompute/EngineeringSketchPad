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

$(TDIR)/sweep:	$(ODIR)/sweep.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/sweep $(ODIR)/sweep.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/sweep.o:	sweep.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) sweep.c -o $(ODIR)/sweep.o

clean:
	-rm $(ODIR)/sweep.o

cleanall:	clean
	-rm $(TDIR)/sweep
