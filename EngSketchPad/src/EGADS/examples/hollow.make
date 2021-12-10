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

$(TDIR)/hollow:	$(ODIR)/hollow.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/hollow $(ODIR)/hollow.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/hollow.o:	hollow.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) hollow.c -o $(ODIR)/hollow.o

clean:
	-rm $(ODIR)/hollow.o 

cleanall:	clean
	-rm $(TDIR)/hollow
