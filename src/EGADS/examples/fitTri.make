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

$(TDIR)/fitTri:	$(ODIR)/fitTri.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/fitTri $(ODIR)/fitTri.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/fitTri.o:	fitTri.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) fitTri.c -o $(ODIR)/fitTri.o

clean:
	-rm $(ODIR)/fitTri.o 

cleanall:	clean
	-rm $(TDIR)/fitTri
