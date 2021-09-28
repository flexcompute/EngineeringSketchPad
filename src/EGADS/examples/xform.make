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

$(TDIR)/xform:	$(ODIR)/xform.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/xform $(ODIR)/xform.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/xform.o:	xform.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) xform.c -o $(ODIR)/xform.o

clean:
	-rm $(ODIR)/xform.o

cleanall:	clean
	-rm $(TDIR)/xform
