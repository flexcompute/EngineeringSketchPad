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


$(LDIR)/msesAIM.so:	$(ODIR)/msesAIM.o $(ODIR)/msesUtils.o $(LDIR)/libaimUtil.a \
			$(LDIR)/libutils.a
	$(CC) $(SOFLGS) -o $(LDIR)/msesAIM.so $(ODIR)/msesAIM.o $(ODIR)/msesUtils.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/%.o:	%.c msesUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils $< -o $@

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother

clean:
	-rm $(ODIR)/msesAIM.o $(ODIR)/msesUtils.o

cleanall:	clean
	-rm $(LDIR)/msesAIM.so

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/csmData/airfoilSection.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen msesAIM_dox.cfg; cd latex; make; mv refman.pdf msesAIM.pdf)
#				../../examples/pyCAPS/mses_PyTest.py \

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f msesAIM.tag)
