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

VPATH = $(ODIR)

#OBJS  =	attrUtils.o meshUtils.o cfdUtils.o miscUtils.o feaUtils.o vlmUtils.o nastranUtils.o tecplotUtils.o arrayUtils.o deprecateUtils.o cardUtils.o nastranCards.o tempUtils.o jsonUtils.o pyscriptUtils.o
OBJS  =	attrUtils.o meshUtils.o cfdUtils.o miscUtils.o feaUtils.o vlmUtils.o nastranUtils.o tecplotUtils.o arrayUtils.o deprecateUtils.o cardUtils.o nastranCards.o jsonUtils.o

OBJSP =	vlmSpanSpace.o

ifdef PYTHONINC
OBJSPython = nastranOP2Reader.o
PYINC=-DHAVE_PYTHON -I$(PYTHONINC)
endif

default: $(LDIR)/libutils.a $(IDIR)/meshTypes.h $(IDIR)/miscTypes.h


$(LDIR)/libutils.a:	$(OBJS) $(OBJSPython) $(OBJSP)
	touch $(LDIR)/libutils.a
	rm $(LDIR)/libutils.a
	(cd $(ODIR); ar -rs $(LDIR)/libutils.a $(OBJS) $(OBJSPython) $(OBJSP) )

$(OBJS) $(OBJSPython): %.o: %.c %.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I$(IDIR) -I. $(PYINC) $< -o $(ODIR)/$@
	
$(OBJSP): %.o: %.cpp %.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I../include -I$(IDIR) -I. $< -o $(ODIR)/$@

$(IDIR)/meshTypes.h: meshTypes.h
	touch $(IDIR)/meshTypes.h
	rm $(IDIR)/meshTypes.h
	cp -p meshTypes.h $(IDIR)/meshTypes.h

$(IDIR)/miscTypes.h: miscTypes.h
	touch $(IDIR)/miscTypes.h
	rm $(IDIR)/miscTypes.h
	cp -p miscTypes.h $(IDIR)/miscTypes.h

LINT_FLAGS=-uniondef -exportlocal -nullpass -retvalother -evalorder -nullderef \
           -initallelements -bufferoverflowhigh -boolops -dependenttrans -kepttrans -immediatetrans
lint:
	$(LINT) -I../../include -I$(IDIR) $(patsubst %.o,%.c,$(OBJS)) $(LINT_FLAGS)

clean:
	-(cd $(ODIR); rm -f $(OBJS) $(OBJSPython) $(OBJSP) )

cleanall:	clean
	-rm -f $(LDIR)/libutils.a

cython:
	(cython -2 nastranOP2Reader.pyx -I../../pyCAPS/src/ --capi-reexport-cincludes) 

