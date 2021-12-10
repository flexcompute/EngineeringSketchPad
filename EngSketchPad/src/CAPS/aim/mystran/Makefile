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


$(LDIR)/mystranAIM.so:	$(ODIR)/mystranAIM.o $(ODIR)/mystranUtils.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/mystranAIM.so $(ODIR)/mystranAIM.o $(ODIR)/mystranUtils.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm $(PYTHONLIB)

$(ODIR)/mystranAIM.o:	mystranAIM.c mystranUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils mystranAIM.c \
		-o $(ODIR)/mystranAIM.o

$(ODIR)/mystranUtils.o:	mystranUtils.c mystranUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils mystranUtils.c \
		-o $(ODIR)/mystranUtils.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef

clean:
	-rm $(ODIR)/mystranAIM.o $(ODIR)/mystranUtils.o

cleanall:	clean
	-rm $(LDIR)/mystranAIM.so
	
dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/mystran_PyTest.py \
				../../examples/csmData/feaWingBEM.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen mystranAIM_dox.cfg; cd latex; make; mv refman.pdf mystranAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f mystranAIM.tag)
