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

$(TDIR)/egads2cart:	$(ODIR)/egads2cart.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/egads2cart $(ODIR)/egads2cart.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/egads2cart.o:	egads2cart.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egads2cart.c \
		-o $(ODIR)/egads2cart.o

clean:
	-rm $(ODIR)/egads2cart.o

cleanall:	clean
	-rm $(TDIR)/egads2cart
