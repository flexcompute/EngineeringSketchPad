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

$(TDIR)/globalTess:	$(ODIR)/globalTess.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/globalTess $(ODIR)/globalTess.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/globalTess.o:	globalTess.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) globalTess.c -o $(ODIR)/globalTess.o

clean:
	-rm $(ODIR)/globalTess.o

cleanall:	clean
	-rm $(TDIR)/globalTess
