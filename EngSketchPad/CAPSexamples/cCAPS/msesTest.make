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

$(TDIR)/msesTest:	$(ODIR)/msesTest.o $(LDIR)/$(CSHLIB) | $(TDIR) 
	$(CC) -o $(TDIR)/msesTest $(ODIR)/msesTest.o -L$(LDIR) -lcaps -legads -locsm \
		-ludunits2 $(RPATH) -lm -ldl

$(ODIR)/msesTest.o:	msesTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) msesTest.c -o $(ODIR)/msesTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/msesTest.o

cleanall:	clean
	-rm -f $(TDIR)/msesTest
