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

$(TDIR)/blend:	$(ODIR)/blend.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/blend $(ODIR)/blend.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/blend.o:	blend.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) blend.c -o $(ODIR)/blend.o

clean:
	-rm $(ODIR)/blend.o 

cleanall:	clean
	-rm $(TDIR)/blend
