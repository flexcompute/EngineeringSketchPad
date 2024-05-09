#
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

default: $(LDIR)/libmeshbWriter.so

VPATH = $(ODIR)

OBJS = libmeshbWriter.o

$(LDIR)/libmeshbWriter.so:	$(ODIR)/libmeshbWriter.o $(ODIR)/libmeshb7.o $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/libmeshbWriter.so $(ODIR)/libmeshbWriter.o $(ODIR)/libmeshb7.o \
		-L$(LDIR) -laimUtil -legads -locsm -ludunits2 -ldl $(RPATH) -lm

$(ODIR)/libmeshbWriter.o:	libmeshbWriter.c libmeshbWriter.h $(IDIR)/aimUtil.h $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I$(IDIR) libmeshbWriter.c \
		-o $(ODIR)/libmeshbWriter.o

$(ODIR)/libmeshb7.o:	libMeshb/sources/libmeshb7.c libMeshb/sources/libmeshb7.h
	$(CC) -c $(COPTS) $(DEFINE)  -IlibMeshb/sources -I../include -I$(IDIR) libMeshb/sources/libmeshb7.c \
		-o $(ODIR)/libmeshb7.o

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:
	$(LINT) -I../../../include -I$(IDIR) $(patsubst %.o,%.c,$(OBJS)) $(LINT_FLAGS)

clean:
	-(cd $(ODIR); rm -f $(OBJS) libmeshb7.o)

cleanall:	clean
	-rm -f $(LDIR)/libmeshbWriter.so
