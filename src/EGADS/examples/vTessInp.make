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

$(TDIR)/vTessInp:	$(ODIR)/vTessInp.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vTessInp $(ODIR)/vTessInp.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(ODIR)/vTessInp.o:	vTessInp.c $(IDIR)/egads.h $(IDIR)/wsserver.h \
			$(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vTessInp.c -o $(ODIR)/vTessInp.o

clean:
	-rm $(ODIR)/vTessInp.o

cleanall:	clean
	-rm $(TDIR)/vTessInp
