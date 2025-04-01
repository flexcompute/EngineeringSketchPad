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

default: $(LDIR)/bstlWriter.so

VPATH = $(ODIR)

OBJS = bstlWriter.o

$(LDIR)/bstlWriter.so:	$(ODIR)/bstlWriter.o $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/bstlWriter.so $(ODIR)/bstlWriter.o \
		-L$(LDIR) -laimUtil -legads -locsm -ludunits2 -ldl $(RPATH) -lm

$(ODIR)/bstlWriter.o:	bstlWriter.c $(IDIR)/aimUtil.h $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I$(IDIR) bstlWriter.c \
		-o $(ODIR)/bstlWriter.o

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:
	$(LINT) -I../../../include -I$(IDIR) $(patsubst %.o,%.c,$(OBJS)) $(LINT_FLAGS)

clean:
	-(cd $(ODIR); rm -f $(OBJS) )

cleanall:	clean
	-rm -f $(LDIR)/bstlWriter.so
