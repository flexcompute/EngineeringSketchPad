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

$(TDIR)/edges:	$(ODIR)/edges.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/edges $(ODIR)/edges.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/edges.o:	edges.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) edges.c -o $(ODIR)/edges.o

clean:
	-rm $(ODIR)/edges.o 

cleanall:	clean
	-rm $(TDIR)/edges
