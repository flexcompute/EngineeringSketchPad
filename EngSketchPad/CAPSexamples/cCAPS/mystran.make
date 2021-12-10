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

$(TDIR)/mystranTest:    $(ODIR)/mystranTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/mystranTest $(ODIR)/mystranTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/mystranTest.o:  mystranTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) mystranTest.c -o $(ODIR)/mystranTest.o   

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/mystranTest.o 

cleanall:	clean
	-rm -f $(TDIR)/mystranTest
