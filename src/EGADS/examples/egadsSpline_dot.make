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

$(TDIR)/egadsSpline_dot:	$(ODIR)/egadsSpline_dot.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/egadsSpline_dot $(ODIR)/egadsSpline_dot.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/egadsSpline_dot.o:	egadsSpline_dot.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egadsSpline_dot.c -o $(ODIR)/egadsSpline_dot.o

clean:
	-rm -f $(ODIR)/egadsSpline_dot.o

cleanall:	clean
	-rm -f $(TDIR)/egadsSpline_dot
