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

$(TDIR)/vAttr:	$(ODIR)/vAttr.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vAttr $(ODIR)/vAttr.o \
		-L$(LDIR) -lwsserver -legads $(RPATH)

$(ODIR)/vAttr.o:	vAttr.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vAttr.c -o $(ODIR)/vAttr.o

clean:
	-rm $(ODIR)/vAttr.o

cleanall:	clean
	-rm $(TDIR)/vAttr
