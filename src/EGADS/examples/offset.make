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

$(TDIR)/offset:	$(ODIR)/offset.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/offset $(ODIR)/offset.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/offset.o:	offset.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) offset.c -o $(ODIR)/offset.o

clean:
	-rm $(ODIR)/offset.o

cleanall:	clean
	-rm $(TDIR)/offset
