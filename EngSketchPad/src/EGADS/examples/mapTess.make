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

$(TDIR)/mapTess:	$(ODIR)/mapTess.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/mapTess $(ODIR)/mapTess.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/mapTess.o:	mapTess.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) mapTess.c -o $(ODIR)/mapTess.o

clean:
	-rm $(ODIR)/mapTess.o

cleanall:	clean
	-rm $(TDIR)/mapTess
