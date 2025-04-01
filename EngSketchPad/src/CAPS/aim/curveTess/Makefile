#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#
include $(ESP_ROOT)/include/$(ESP_ARCH)

IDIR=$(ESP_ROOT)/include
LDIR=$(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR=$(ESP_BLOC)/obj
TDIR=$(ESP_BLOC)/test
else
ODIR=.
TDIR=$(ESP_ROOT)/bin
endif

CFLAGS=
CPPFLAGS=$(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)/aim/utils
LDFLAGS=$(SOFLGS) -L$(LDIR) $(RPATH)
LIBS=-lutils -laimUtil -locsm -legads -ludunits2 -lm -ldl

all: $(LDIR)/curveTessAIM.so

clean:
	-rm $(ODIR)/curveTessAIM.o

cleanall:	clean
	-rm $(LDIR)/curveTessAIM.so

$(ODIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(LDIR)/curveTessAIM.so: $(ODIR)/curveTessAIM.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(LDFLAGS) $(ODIR)/curveTessAIM.o -o $@ $(LIBS)

	
lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -bufferoverflowhigh -uniondef


dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen curveTessAIM_dox.cfg; cd latex; make; mv refman.pdf curveTessAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f curveTessAIM.tag)
	