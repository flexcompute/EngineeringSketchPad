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

$(TDIR)/fillet:	$(ODIR)/fillet.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/fillet $(ODIR)/fillet.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/fillet.o:	fillet.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) fillet.c -o $(ODIR)/fillet.o

clean:
	-rm $(ODIR)/fillet.o 

cleanall:	clean
	-rm $(TDIR)/fillet
