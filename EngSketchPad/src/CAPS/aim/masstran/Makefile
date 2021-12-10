#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
endif

default: $(LDIR)/masstranAIM.so #$(BDIR)/Masstran

$(LDIR)/masstranAIM.so:	$(ODIR)/masstranAIM.o $(LDIR)/libutils.a \
		$(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/masstranAIM.so $(ODIR)/masstranAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/masstranAIM.o:	masstranAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils masstranAIM.c \
		-o $(ODIR)/masstranAIM.o

#$(BDIR)/Masstran: Masstran.c
#	$(CC) $(COPTS) Masstran.c -o $(BDIR)/Masstran -lm

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother -boolops

clean:
	-rm -f $(ODIR)/masstranAIM.o 

cleanall:	clean
	-rm -f $(LDIR)/masstranAIM.so $(BDIR)/Masstran

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/masstran_PyTest.py \
				../../examples/csmData/feaWingBEM.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen masstranAIM_dox.cfg; cd latex; make; mv refman.pdf masstranAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f masstranAIM.tag)
