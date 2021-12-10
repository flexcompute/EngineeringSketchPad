#
include ../include/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

default:	$(TDIR)/patchTest

$(TDIR)/patchTest:	$(ODIR)/patchTest.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/patchTest $(ODIR)/patchTest.o \
		-L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/patchTest.o:	egadsPatch.c ../include/egads.h \
			../include/egadsTypes.h ../include/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -DSTANDALONE egadsPatch.c \
		-o $(ODIR)/patchTest.o

clean:
	-rm $(ODIR)/patchTest.o

cleanall:	clean
	-rm $(TDIR)/patchTest
