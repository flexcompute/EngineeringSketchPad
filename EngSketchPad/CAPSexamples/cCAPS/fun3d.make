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

fun3d: $(TDIR)/fun3dAFLR2Test $(TDIR)/fun3dTetgenTest \
	$(TDIR)/aeroelasticTest

$(TDIR)/fun3dAFLR2Test:  $(ODIR)/fun3dAFLR2Test.o $(LDIR)/$(CSHLIB) | $(TDIR)
	$(CC) -o $(TDIR)/fun3dAFLR2Test $(ODIR)/fun3dAFLR2Test.o \
	-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(TDIR)/fun3dTetgenTest:  $(ODIR)/fun3dTetgenTest.o $(LDIR)/$(CSHLIB) | $(TDIR)
	$(CC) -o $(TDIR)/fun3dTetgenTest $(ODIR)/fun3dTetgenTest.o \
	-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl

$(TDIR)/aeroelasticTest:  $(ODIR)/aeroelasticTest.o $(LDIR)/$(CSHLIB) | $(TDIR)
	$(CC) -o $(TDIR)/aeroelasticTest $(ODIR)/aeroelasticTest.o \
	-L$(LDIR) -lcaps -locsm -legads -ludunits2 $(RPATH) -lm -ldl


$(ODIR)/fun3dAFLR2Test.o:  fun3dAFLR2Test.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) fun3dAFLR2Test.c -o $(ODIR)/fun3dAFLR2Test.o

$(ODIR)/fun3dTetgenTest.o:  fun3dTetgenTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) fun3dTetgenTest.c -o $(ODIR)/fun3dTetgenTest.o
	
$(ODIR)/aeroelasticTest.o:  aeroelasticTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) aeroelasticTest.c -o $(ODIR)/aeroelasticTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/fun3dAFLR2Test.o $(ODIR)/aeroelasticTest.o $(ODIR)/fun3dTetgenTest.o

cleanall:	clean
	-rm -f $(TDIR)/fun3dAFLR2Test $(TDIR)/aeroelasticTest $(TDIR)/fun3dTetgenTest
