#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib

tire:	tire.o $(LDIR)/$(SHLIB)
	$(CXX) -o tire tire.o -L$(LDIR) -legads $(RPATH) -lm

tire.o:	tire.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) tire.c

clean:
	-rm -f tire.o

cleanall:	clean
	-rm -f tire
