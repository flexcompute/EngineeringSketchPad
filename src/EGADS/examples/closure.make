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

$(TDIR)/closure:	$(ODIR)/closure.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/closure $(ODIR)/closure.o -L$(LDIR) -lfgads \
		-legads $(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/closure.o:	closure.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) closure.f -o $(ODIR)/closure.o

clean:
	-rm $(ODIR)/closure.o

cleanall:	clean
	-rm $(TDIR)/closure
