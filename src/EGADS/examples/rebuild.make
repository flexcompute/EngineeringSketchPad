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

$(TDIR)/rebuild:	$(ODIR)/rebuild.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/rebuild $(ODIR)/rebuild.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/rebuild.o:	rebuild.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) rebuild.c -o $(ODIR)/rebuild.o

clean:
	-rm $(ODIR)/rebuild.o 

cleanall:	clean
	-rm $(TDIR)/rebuild
