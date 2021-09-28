#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/examples/cCAPS
else
ODIR  = .
TDIR  = .
endif

$(TDIR)/pointwiseTest:	$(ODIR)/pointwiseTest.o $(LDIR)/$(CSHLIB) | $(TDIR) 
	$(CC) -o $(TDIR)/pointwiseTest $(ODIR)/pointwiseTest.o -L$(LDIR) -lcaps -locsm -legads \
		-ludunits2 $(RPATH) -lm -ldl

$(ODIR)/pointwiseTest.o:	pointwiseTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) pointwiseTest.c -o $(ODIR)/pointwiseTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/pointwiseTest.o

cleanall:	clean
	-rm -f $(TDIR)/pointwiseTest
