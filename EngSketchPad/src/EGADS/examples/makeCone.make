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

$(TDIR)/makeCone:	$(ODIR)/makeCone.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/makeCone $(ODIR)/makeCone.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/makeCone.o:	makeCone.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) makeCone.c -o $(ODIR)/makeCone.o

clean:
	-rm $(ODIR)/makeCone.o 

cleanall:	clean
	-rm $(TDIR)/makeCone
