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


$(LDIR)/xfoilAIM.so:	$(ODIR)/xfoilAIM.o $(LDIR)/libaimUtil.a \
			$(LDIR)/libutils.a
	$(CC) $(SOFLGS) -o $(LDIR)/xfoilAIM.so $(ODIR)/xfoilAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/xfoilAIM.o: xfoilAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils \
		xfoilAIM.c -o $(ODIR)/xfoilAIM.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother

clean:
	-rm $(ODIR)/xfoilAIM.o

cleanall:	clean
	-rm $(LDIR)/xfoilAIM.so

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/xfoil_PyTest.py \
				../../examples/csmData/airfoilSection.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen xfoilAIM_dox.cfg; cd latex; make; mv refman.pdf xfoilAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f xfoilAIM.tag)
