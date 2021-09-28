#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
endif


$(TDIR)/egads2cgt:	$(ODIR)/egads2cgt.o $(ODIR)/aflr4egads.o
	$(CXX) -o $(TDIR)/egads2cgt $(ODIR)/egads2cgt.o $(ODIR)/aflr4egads.o \
		-L$(LDIR) -legads -lcaps -locsm $(RPATH) -lm -ldl

$(ODIR)/egads2cgt.o:	egads2cgt.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egads2cgt.c \
		-o $(ODIR)/egads2cgt.o

$(ODIR)/aflr4egads.o:	aflr4egads.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) aflr4egads.c \
		-o $(ODIR)/aflr4egads.o

lint:
	$(LINT) -I../../include -I$(IDIR) $(wildcard *.c) -uniondef -exportlocal -nullpass

clean:
	-rm $(ODIR)/egads2cgt.o $(ODIR)/aflr4egads.o

cleanall:	clean
	-rm $(TDIR)/egads2cgt
