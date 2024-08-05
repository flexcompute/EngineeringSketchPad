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

# Directory to install header files
UTILDIR=$(IDIR)/aim/utils

VPATH = $(ODIR):cython

#OBJS  =	attrUtils.o meshUtils.o cfdUtils.o miscUtils.o feaUtils.o vlmUtils.o nastranUtils.o tecplotUtils.o arrayUtils.o deprecateUtils.o cardUtils.o nastranCards.o tempUtils.o jsonUtils.o pyscriptUtils.o
OBJS = attrUtils.o \
       meshUtils.o \
       cfdUtils.o \
       miscUtils.o \
       feaUtils.o \
       vlmUtils.o \
       nastranCards.o \
       nastranUtils.o \
       tecplotUtils.o \
       arrayUtils.o \
       deprecateUtils.o \
       cardUtils.o \
       jsonUtils.o

OBJSP =	vlmSpanSpace.o

ifdef PYTHONINC
OBJSPython = nastranOP2Reader.o
PYINC=-DHAVE_PYTHON -I$(PYTHONINC)

#Extract the python include directory
PYTHONIPATH0=$(word 1,$(filter -I%,$(PYTHONINC)))
ifeq ("$(PYTHONIPATH0)","")
PYTHONIPATH=$(PYTHONINC)
else
PYTHONIPATH=$(PYTHONIPATH0:-I%=%)
endif
PYINCDEP=$(PYTHONIPATH)/Python.h
endif

HFILES=$(patsubst %.h,$(UTILDIR)/%.h,$(wildcard *.h))

default: $(LDIR)/libutils.a $(HFILES)

$(LDIR)/libutils.a:	$(OBJS) $(OBJSPython) $(OBJSP)
	touch $(LDIR)/libutils.a
	rm $(LDIR)/libutils.a
	(cd $(ODIR); ar -rs $(LDIR)/libutils.a $(OBJS) $(OBJSPython) $(OBJSP) )

$(OBJS) $(OBJSPython): %.o: %.c %.h $(PYINCDEP)
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I$(IDIR) -I. $(PYINC) -Icython $< -o $(ODIR)/$@
	
$(OBJSP): %.o: %.cpp %.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I../include -I$(IDIR) -I. $< -o $(ODIR)/$@

$(UTILDIR):
	mkdir -p $(UTILDIR)

$(HFILES): $(UTILDIR)/%.h: %.h | $(UTILDIR)
	cp -p -f $< $@

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:
	$(LINT) -I../../include -I$(IDIR) $(patsubst %.o,%.c,$(OBJS)) $(LINT_FLAGS)

clean:
	-(cd $(ODIR); rm -f $(OBJS) $(OBJSPython) $(OBJSP) )

cleanall:	clean
	-rm -f $(LDIR)/libutils.a

.PHONY: cython
cython:
	make -C cython
