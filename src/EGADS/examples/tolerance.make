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

$(TDIR)/tolerance:	$(ODIR)/tolerance.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/tolerance $(ODIR)/tolerance.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/tolerance.o:	tolerance.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) tolerance.c -o $(ODIR)/tolerance.o

clean:
	-rm $(ODIR)/tolerance.o

cleanall:	clean
	-rm $(TDIR)/tolerance
