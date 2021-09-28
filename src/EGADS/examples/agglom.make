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

$(TDIR)/agglom:	$(ODIR)/agglom.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/agglom $(ODIR)/agglom.o -L$(LDIR) $(RPATH) -legads -lm

$(ODIR)/agglom.o:	agglom.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) agglom.c -o $(ODIR)/agglom.o

clean:
	-rm $(ODIR)/agglom.o

cleanall:	clean
	-rm $(TDIR)/agglom
