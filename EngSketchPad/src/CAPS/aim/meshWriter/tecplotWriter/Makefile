#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
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

default: $(LDIR)/tecplotWriter.so

VPATH = $(ODIR)

OBJS = tecplotWriter.o

$(LDIR)/tecplotWriter.so:	$(ODIR)/tecplotWriter.o $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/tecplotWriter.so $(ODIR)/tecplotWriter.o \
		-L$(LDIR) -laimUtil -legads -locsm -ludunits2 -ldl $(RPATH) -lm

$(ODIR)/tecplotWriter.o:	tecplotWriter.c $(IDIR)/aimUtil.h $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I$(IDIR) tecplotWriter.c \
		-o $(ODIR)/tecplotWriter.o

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:
	$(LINT) -I../../../include -I$(IDIR) $(patsubst %.o,%.c,$(OBJS)) $(LINT_FLAGS)

clean:
	-(cd $(ODIR); rm -f $(OBJS) )

cleanall:	clean
	-rm -f $(LDIR)/tecplotWriter.so
