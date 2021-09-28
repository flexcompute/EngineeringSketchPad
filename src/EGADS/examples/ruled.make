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

$(TDIR)/ruled:	$(ODIR)/ruled.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/ruled $(ODIR)/ruled.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/ruled.o:	ruled.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) ruled.c -o $(ODIR)/ruled.o

clean:
	-rm $(ODIR)/ruled.o

cleanall:	clean
	-rm $(TDIR)/ruled
