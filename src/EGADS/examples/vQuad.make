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

$(TDIR)/vQuad:	$(ODIR)/vQuad.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vQuad $(ODIR)/vQuad.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(ODIR)/vQuad.o:	vQuad.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vQuad.c -o $(ODIR)/vQuad.o

clean:
	-rm $(ODIR)/vQuad.o

cleanall:	clean
	-rm $(TDIR)/vQuad
