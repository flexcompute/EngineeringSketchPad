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

$(TDIR)/vEffect:	$(ODIR)/vEffect.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vEffect $(ODIR)/vEffect.o \
		-L$(LDIR) -lwsserver -legads $(RPATH)

$(ODIR)/vEffect.o:	vEffect.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vEffect.c -o $(ODIR)/vEffect.o

clean:
	-rm $(ODIR)/vEffect.o

cleanall:	clean
	-rm $(TDIR)/vEffect
