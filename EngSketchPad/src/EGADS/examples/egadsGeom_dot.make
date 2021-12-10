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

$(TDIR)/egadsGeom_dot:	$(ODIR)/egadsGeom_dot.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/egadsGeom_dot $(ODIR)/egadsGeom_dot.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/egadsGeom_dot.o:	egadsGeom_dot.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egadsGeom_dot.c -o $(ODIR)/egadsGeom_dot.o

clean:
	-rm -f $(ODIR)/egadsGeom_dot.o

cleanall:	clean
	-rm -f $(TDIR)/egadsGeom_dot
