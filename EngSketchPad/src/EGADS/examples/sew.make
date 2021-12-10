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

$(TDIR)/sew:	$(ODIR)/sew.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/sew $(ODIR)/sew.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/sew.o:	sew.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) sew.c -o $(ODIR)/sew.o

clean:
	-rm $(ODIR)/sew.o

cleanall:	clean
	-rm $(TDIR)/sew
