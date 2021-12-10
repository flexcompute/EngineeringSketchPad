#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib

egads2tri:	egads2tri.o $(LDIR)/$(SHLIB)
	$(CXX) -o egads2tri egads2tri.o -L$(LDIR) -legads $(RPATH) -lm

egads2tri.o:	egads2tri.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egads2tri.c

clean:
	-rm egads2tri.o

cleanall:	clean
	-rm egads2tri
