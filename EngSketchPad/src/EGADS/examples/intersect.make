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

$(TDIR)/intersect:	$(ODIR)/intersect.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/intersect $(ODIR)/intersect.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/intersect.o:	intersect.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) intersect.c -o $(ODIR)/intersect.o

clean:
	-rm $(ODIR)/intersect.o 

cleanall:	clean
	-rm $(TDIR)/intersect
