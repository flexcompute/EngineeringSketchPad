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

$(TDIR)/vGeom:	$(ODIR)/vGeom.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vGeom $(ODIR)/vGeom.o -L$(LDIR) -lwsserver -legads \
		-lpthread -lz $(RPATH) -lm

$(ODIR)/vGeom.o:	vGeom.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vGeom.c -o $(ODIR)/vGeom.o

clean:
	-rm $(ODIR)/vGeom.o

cleanall:	clean
	-rm $(TDIR)/vGeom
