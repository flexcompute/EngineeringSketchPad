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


$(LDIR)/tacsAIM.so:	$(ODIR)/tacsAIM.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/tacsAIM.so $(ODIR)/tacsAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm $(CPPSLB) $(PYTHONLIB)

$(ODIR)/tacsAIM.o:	tacsAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils tacsAIM.c \
		-o $(ODIR)/tacsAIM.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef

clean:
	-rm $(ODIR)/tacsAIM.o

cleanall:	clean
	-rm $(LDIR)/tacsAIM.so
	
dox: doc/html/index.html
doc/html/index.html: \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen tacsAIM_dox.cfg; cd latex; make; mv refman.pdf tacsAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f tacsAIM.tag)
