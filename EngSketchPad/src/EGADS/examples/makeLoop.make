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

$(TDIR)/makeLoop:	$(ODIR)/makeLoop.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/makeLoop $(ODIR)/makeLoop.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/makeLoop.o:	makeLoop.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) makeLoop.c -o $(ODIR)/makeLoop.o

clean:
	-rm $(ODIR)/makeLoop.o

cleanall:	clean
	-rm $(TDIR)/makeLoop
