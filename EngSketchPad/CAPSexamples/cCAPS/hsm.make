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

hsm:  $(TDIR)/hsmTest $(TDIR)/hsmSimplePlateTest $(TDIR)/hsmCantileverPlateTest \
	$(TDIR)/hsmJoinedPlateTest
	
$(TDIR)/hsmTest:    $(ODIR)/hsmTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/hsmTest $(ODIR)/hsmTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/hsmTest.o:  hsmTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) hsmTest.c -o $(ODIR)/hsmTest.o   


$(TDIR)/hsmSimplePlateTest:    $(ODIR)/hsmSimplePlateTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/hsmSimplePlateTest $(ODIR)/hsmSimplePlateTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/hsmSimplePlateTest.o:  hsmSimplePlateTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) hsmSimplePlateTest.c -o $(ODIR)/hsmSimplePlateTest.o  


$(TDIR)/hsmCantileverPlateTest:    $(ODIR)/hsmCantileverPlateTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/hsmCantileverPlateTest $(ODIR)/hsmCantileverPlateTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/hsmCantileverPlateTest.o:  hsmCantileverPlateTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) hsmCantileverPlateTest.c -o $(ODIR)/hsmCantileverPlateTest.o  


$(TDIR)/hsmJoinedPlateTest:    $(ODIR)/hsmJoinedPlateTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/hsmJoinedPlateTest $(ODIR)/hsmJoinedPlateTest.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(ODIR)/hsmJoinedPlateTest.o:  hsmJoinedPlateTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) hsmJoinedPlateTest.c -o $(ODIR)/hsmJoinedPlateTest.o  


ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/hsmTest.o $(ODIR)/hsmSimplePlateTest.o $(ODIR)/hsmCantileverPlateTest.o $(ODIR)/hsmJoinedPlateTest.o

cleanall:	clean
	-rm -f $(TDIR)/hsmTest $(TDIR)/hsmSimplePlateTest $(TDIR)/hsmCantileverPlateTest $(TDIR)/hsmJoinedPlateTest
