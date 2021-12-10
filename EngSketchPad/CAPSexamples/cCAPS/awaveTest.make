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

$(TDIR)/awaveTest:	$(ODIR)/awaveTest.o $(LDIR)/$(CSHLIB) | $(TDIR) 
	$(CC) -o $(TDIR)/awaveTest $(ODIR)/awaveTest.o -L$(LDIR) -lcaps -legads -locsm \
		-ludunits2 $(RPATH) -lm -ldl

$(ODIR)/awaveTest.o:	awaveTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) awaveTest.c -o $(ODIR)/awaveTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/awaveTest.o

cleanall:	clean
	-rm -f $(TDIR)/awaveTest
