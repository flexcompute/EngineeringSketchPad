ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif

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

default: $(LDIR)/zaeroAIM.so

OBJECTS = $(ODIR)/zaeroAIM.o $(ODIR)/zaeroUtils.o \
          $(ODIR)/zaeroLinearFlutter.o \
          $(ODIR)/zaeroTrim.o $(ODIR)/zaeroTypes.o \
          $(ODIR)/zaeroCards.o $(ODIR)/zaeroGeneral.o \
          $(ODIR)/zaeroGraphical.o

HEADER  = zaeroCards.h zaeroDiscipline.h \
          zaeroGeneral.h zaeroGraphical.h zaeroTypes.h \
          zaeroUtils.h
LIBS    = -lutils -laimUtil -locsm -legads -ludunits2 -lm -ldl

.PHONY: test clean cleanall

$(LDIR)/zaeroAIM.so: $(OBJECTS)
	$(CC) $(SOFLGS) $(LDFLAGS) -o $@ $(OBJECTS) -L$(LDIR) $(LIBS)

$(OBJECTS): $(ODIR)/%.o: %.c $(HEADER)
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils $< -o $@

clean:
	-rm -f $(OBJECTS)

cleanall: clean
	-rm -f $(LDIR)/zaeroAIM.so

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils *.c -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen zaeroAIM_dox.cfg; cd latex; make; mv refman.pdf zaeroAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f zaeroAIM.tag)

