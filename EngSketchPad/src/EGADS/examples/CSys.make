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

$(TDIR)/CSys:	$(ODIR)/CSys.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/CSys $(ODIR)/CSys.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/CSys.o:	CSys.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) CSys.c -o $(ODIR)/CSys.o

clean:
	-rm $(ODIR)/CSys.o

cleanall:	clean
	-rm $(TDIR)/CSys
