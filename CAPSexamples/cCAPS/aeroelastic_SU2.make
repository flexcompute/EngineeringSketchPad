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

$(TDIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest:	$(ODIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest.o $(LDIR)/$(CSHLIB) | $(TDIR) 
	$(CC) -o $(TDIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest $(ODIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest.o -L$(LDIR) -lcaps -locsm -legads \
		-ludunits2 $(RPATH) -lm -ldl

$(ODIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest.o:	aeroelasticSimple_Iterative_SU2_and_MystranTest.c $(IDIR)/caps.h | $(ODIR)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) aeroelasticSimple_Iterative_SU2_and_MystranTest.c -o $(ODIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest.o

ifdef ESP_BLOC
$(ODIR):
	mkdir -p $@

$(TDIR):
	mkdir -p $@
endif

clean:
	-rm -f $(ODIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest.o

cleanall:	clean
	-rm -f $(TDIR)/aeroelasticSimple_Iterative_SU2_and_MystranTest
