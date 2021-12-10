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

$(TDIR)/replace:	$(ODIR)/replace.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/replace $(ODIR)/replace.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/replace.o:	replace.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) replace.c -o $(ODIR)/replace.o

clean:
	-rm $(ODIR)/replace.o 

cleanall:	clean
	-rm $(TDIR)/replace
