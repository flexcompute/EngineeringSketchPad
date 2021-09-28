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


$(LDIR)/nastranAIM.so:	$(ODIR)/nastranAIM.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/nastranAIM.so $(ODIR)/nastranAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm $(CPPSLB) $(PYTHONLIB)

$(ODIR)/nastranAIM.o:	nastranAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils nastranAIM.c \
		-o $(ODIR)/nastranAIM.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef

clean:
	-rm $(ODIR)/nastranAIM.o

cleanall:	clean
	-rm $(LDIR)/nastranAIM.so
	
dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/nastran_ThreeBar_PyTest.py \
				../../examples/pyCAPS/nastran_ThreeBarMultiLoad_PyTest.py \
				../../examples/pyCAPS/nastran_ThreeBarFreq_PyTest.py \
				../../examples/pyCAPS/nastran_ThreeBarDesign_PyTest.py \
				../../examples/pyCAPS/nastran_CompositeWingFreq_PyTest.py \
				../../examples/csmData/feaThreeBar.csm \
				../../examples/csmData/compositeWing.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen nastranAIM_dox.cfg; cd latex; make; mv refman.pdf nastranAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f nastranAIM.tag)
