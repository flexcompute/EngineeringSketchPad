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

$(TDIR)/avlTest:	$(ODIR)/avlTest.o $(LDIR)/$(CSHLIB) | $(TDIR) 
	$(CC) -o $(TDIR)/avlTest $(ODIR)/avlTest.o -L$(LDIR) -lcaps -legads -locsm \
		-ludunits2 $(RPATH) -lm -ldl

$(ODIR)/avlTest.o:	avlTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) avlTest.c -o $(ODIR)/avlTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/avlTest.o

cleanall:	clean
	-rm -f $(TDIR)/avlTest
