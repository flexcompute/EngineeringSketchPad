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

$(TDIR)/chamfer:	$(ODIR)/chamfer.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/chamfer $(ODIR)/chamfer.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/chamfer.o:	chamfer.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) chamfer.c -o $(ODIR)/chamfer.o

clean:
	-rm $(ODIR)/chamfer.o

cleanall:	clean
	-rm $(TDIR)/chamfer
