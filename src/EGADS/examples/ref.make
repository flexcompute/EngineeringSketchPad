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

$(TDIR)/ref:	$(ODIR)/ref.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/ref $(ODIR)/ref.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/ref.o:	ref.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) ref.c -o $(ODIR)/ref.o

clean:
	-rm $(ODIR)/ref.o

cleanall:	clean
	-rm $(TDIR)/ref
