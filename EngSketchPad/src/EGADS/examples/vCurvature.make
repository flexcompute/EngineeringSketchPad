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

$(TDIR)/vCurvature:	$(ODIR)/vCurvature.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vCurvature $(ODIR)/vCurvature.o \
		-L$(LDIR) -lwsserver -legads $(RPATH)

$(ODIR)/vCurvature.o:	vCurvature.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vCurvature.c -o $(ODIR)/vCurvature.o

clean:
	-rm $(ODIR)/vCurvature.o

cleanall:	clean
	-rm $(TDIR)/vCurvature
