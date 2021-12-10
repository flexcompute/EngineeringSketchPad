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

$(TDIR)/egadsTopo_dot:	$(ODIR)/egadsTopo_dot.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/egadsTopo_dot $(ODIR)/egadsTopo_dot.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/egadsTopo_dot.o:	egadsTopo_dot.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egadsTopo_dot.c -o $(ODIR)/egadsTopo_dot.o

clean:
	-rm $(ODIR)/egadsTopo_dot.o

cleanall:	clean
	-rm $(TDIR)/egadsTopo_dot
