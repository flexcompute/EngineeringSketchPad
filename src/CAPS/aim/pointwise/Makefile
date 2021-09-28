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

ifdef GELITE
default:	$(LDIR)/pointwiseAIM.so $(TDIR)/egads2nmb
else
default:	$(LDIR)/pointwiseAIM.so
endif

$(LDIR)/pointwiseAIM.so:	$(ODIR)/pointwiseAIM.o $(ODIR)/hashElement.o \
			$(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	touch $(LDIR)/pointwiseAIM.so
	rm $(LDIR)/pointwiseAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/pointwiseAIM.so $(ODIR)/pointwiseAIM.o $(ODIR)/hashElement.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl \
		$(RPATH) -lm -lpthread

$(ODIR)/pointwiseAIM.o:	pointwiseAIM.c hashElement.h $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils pointwiseAIM.c \
		-o $(ODIR)/pointwiseAIM.o

$(ODIR)/hashElement.o:	hashElement.c hashElement.h $(IDIR)/capsTypes.h 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils hashElement.c \
		-o $(ODIR)/hashElement.o


$(TDIR)/egads2nmb:	$(ODIR)/egads2nmb.o
	$(CXX) -o $(TDIR)/egads2nmb $(ODIR)/egads2nmb.o \
		-L$(LDIR) -lnmb -legads -lpthread $(RPATH)

$(ODIR)/egads2nmb.o:	egads2nmb.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. egads2nmb.c \
		-o $(ODIR)/egads2nmb.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils pointwiseAIM.c hashElement.c -exportlocal -nestcomment -evalorder -retvalother -boolops -shiftimplementation -initallelements
	
clean:
	-rm -f $(ODIR)/egads2nmb.o $(ODIR)/pointwiseAIM.o $(ODIR)/hashElement.o

cleanall:	clean
	-rm -f $(TDIR)/egads2nmb $(LDIR)/pointwiseAIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen pointwiseAIM_dox.cfg; cd latex; make; mv refman.pdf pointwiseAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f pointwiseAIM.tag)
	
