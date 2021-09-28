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

$(TDIR)/makeFace3D:	$(ODIR)/makeFace3D.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/makeFace3D $(ODIR)/makeFace3D.o -L$(LDIR) -legads \
		$(RPATH) -lm

$(ODIR)/makeFace3D.o:	makeFace3D.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) makeFace3D.c \
		-o $(ODIR)/makeFace3D.o

clean:
	-rm $(ODIR)/makeFace3D.o

cleanall:	clean
	-rm $(TDIR)/makeFace3D
