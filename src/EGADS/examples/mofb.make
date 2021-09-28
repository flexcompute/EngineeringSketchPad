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

$(TDIR)/mofb:	$(ODIR)/mofb.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/mofb $(ODIR)/mofb.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/mofb.o:	mofb.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) mofb.c -o $(ODIR)/mofb.o

clean:
	-rm $(ODIR)/mofb.o

cleanall:	clean
	-rm $(TDIR)/mofb
