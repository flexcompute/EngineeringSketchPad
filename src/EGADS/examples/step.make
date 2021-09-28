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

$(TDIR)/step:	$(ODIR)/step.o $(LDIR)/libfgads.a $(LDIR)/$(SHLIB)
	$(FCOMP) -o $(TDIR)/step $(ODIR)/step.o -L$(LDIR) -lfgads -legads \
		$(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/step.o:	step.f
	$(FCOMP) -c $(FOPTS) -I$(IDIR) step.f -o $(ODIR)/step.o

clean:
	-rm $(ODIR)/step.o

cleanall:	clean
	-rm $(TDIR)/step
