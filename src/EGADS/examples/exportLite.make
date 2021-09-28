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

$(TDIR)/exportLite:	$(ODIR)/exportLite.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/exportLite $(ODIR)/exportLite.o  \
		 -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/exportLite.o:	../src/egadsExport.c $(IDIR)/egads.h \
		$(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../util/uvmap \
		../src/egadsExport.c -o $(ODIR)/exportLite.o -DSTANDALONE

clean:
	-rm $(ODIR)/exportLite.o

cleanall:	clean
	-rm $(TDIR)/exportLite
