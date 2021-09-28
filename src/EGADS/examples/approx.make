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

$(TDIR)/approx:	$(ODIR)/approx.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/approx $(ODIR)/approx.o -L$(LDIR) -lfgads \
		-legads $(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/approx.o:	approx.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) approx.f -o $(ODIR)/approx.o

clean:
	-rm $(ODIR)/approx.o

cleanall:	clean
	-rm $(TDIR)/approx
