#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
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


default: $(LDIR)/flightstreamAIM.so

$(LDIR)/flightstreamAIM.so: $(ODIR)/flightstreamAIM.o
	$(CC) $(SOFLGS) $^ -o $@ -L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -lm -ldl

$(ODIR)/%.o : %.c
	$(CC) -c $(COPTS) $(DEFINE) -I. -I$(IDIR) -I../utils -I../meshWriter/wavefrontWriter $< -o $@

clean:
	-rm $(ODIR)/flightstreamAIM.o

cleanall:	clean
	-rm $(LDIR)/flightstreamAIM.so

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils -I../meshWriter/wavefrontWriter flightstreamAIM.c -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef

dox: doc/html/index.html
doc/html/index.html: \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc && doxygen flightstreamAIM_dox.cfg && cd latex && make && mv refman.pdf flightstreamAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f flightstreamAIM.tag)
