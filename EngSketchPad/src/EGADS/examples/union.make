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

$(TDIR)/union:	$(ODIR)/union.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/union $(ODIR)/union.o -L$(LDIR) -lfgads -legads \
		$(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/union.o:	union.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) union.f -o $(ODIR)/union.o

clean:
	-rm $(ODIR)/union.o

cleanall:	clean
	-rm $(TDIR)/union
