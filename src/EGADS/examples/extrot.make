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

$(TDIR)/extrot:	$(ODIR)/extrot.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/extrot $(ODIR)/extrot.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/extrot.o:	extrot.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) extrot.c -o $(ODIR)/extrot.o

clean:
	-rm $(ODIR)/extrot.o

cleanall:	clean
	-rm $(TDIR)/extrot
