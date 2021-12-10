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
else
ODIR  = .
endif


$(LDIR)/interferenceAIM.so:	$(ODIR)/interferenceAIM.o $(ODIR)/cloudFns.o \
				$(LDIR)/libaimUtil.a
	touch $(LDIR)/interferenceAIM.so
	rm $(LDIR)/interferenceAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/interferenceAIM.so \
		$(ODIR)/interferenceAIM.o $(ODIR)/cloudFns.o -L$(LDIR) \
		-laimUtil -locsm -legads -ludunits2 -ldl $(CPPSLB) $(RPATH) -lm

$(ODIR)/interferenceAIM.o:	interferenceAIM.c cloud.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) interferenceAIM.c \
		-o $(ODIR)/interferenceAIM.o

$(ODIR)/cloudFns.o:	cloudFns.c cloud.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. cloudFns.c \
		-o $(ODIR)/cloudFns.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -bufferoverflowhigh -uniondef

clean:
	-rm $(ODIR)/interferenceAIM.o $(ODIR)/cloudFns.o

cleanall:	clean
	-rm $(LDIR)/interferenceAIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen interferenceAIM_dox.cfg; cd latex; make; mv refman.pdf interferenceAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f interferenceAIM.tag)
