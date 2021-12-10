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

$(TDIR)/interferenceTest:    $(ODIR)/interferenceTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/interferenceTest $(ODIR)/interferenceTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/interferenceTest.o:  interferenceTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) interferenceTest.c \
		-o $(ODIR)/interferenceTest.o   

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/interferenceTest.o 

cleanall:	clean
	-rm -f $(TDIR)/interferenceTest
