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

$(TDIR)/egadsHLevel_dot:	$(ODIR)/egadsHLevel_dot.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/egadsHLevel_dot $(ODIR)/egadsHLevel_dot.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/egadsHLevel_dot.o:	egadsHLevel_dot.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egadsHLevel_dot.c -o $(ODIR)/egadsHLevel_dot.o

clean:
	-rm $(ODIR)/egadsHLevel_dot.o

cleanall:	clean
	-rm $(TDIR)/egadsHLevel_dot
