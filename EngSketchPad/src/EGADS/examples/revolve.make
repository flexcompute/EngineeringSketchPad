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

$(TDIR)/revolve:	$(ODIR)/revolve.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/revolve $(ODIR)/revolve.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/revolve.o:	revolve.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) revolve.c -o $(ODIR)/revolve.o

clean:
	-rm $(ODIR)/revolve.o 

cleanall:	clean
	-rm $(TDIR)/revolve
