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

ifdef FCOMP
default:	$(LDIR)/frictionAIM.so $(TDIR)/friction
else
default:	$(LDIR)/frictionAIM.so
endif


$(LDIR)/frictionAIM.so:	$(ODIR)/frictionAIM.o $(LDIR)/libutils.a  $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/frictionAIM.so $(ODIR)/frictionAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/frictionAIM.o:	frictionAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -Wno-format -I$(IDIR) -I../utils frictionAIM.c \
		-o $(ODIR)/frictionAIM.o

$(TDIR)/friction:	$(ODIR)/friction.o
	$(FCOMP) -o $(TDIR)/friction $(ODIR)/friction.o $(ESPFLIBS)

$(ODIR)/friction.o:	friction_eja_mod.f
	$(FCOMP) -c $(FOPTS) friction_eja_mod.f -o $(ODIR)/friction.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -evalorder -nestcomment -boolops -retvalother

clean:
	-rm $(ODIR)/frictionAIM.o $(ODIR)/friction.o

cleanall:	clean
	-rm $(LDIR)/frictionAIM.so $(TDIR)/friction

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/friction_PyTest.py \
				../../examples/csmData/frictionWingTailFuselage.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen frictionAIM_dox.cfg; cd latex; make; mv refman.pdf frictionAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f frictionAIM.tag)	
