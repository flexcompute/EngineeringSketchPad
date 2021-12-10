#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = src
TDIR  = $(ESP_ROOT)/bin
endif

ifdef FCOMP
default:	$(LDIR)/hsmAIM.so
else
default:
endif

$(LDIR)/hsmAIM.so:	$(ODIR)/hsmAIM.o $(ODIR)/hsmUtils.o $(ODIR)/hsmAdj.o \
			$(ODIR)/rcm.o $(ODIR)/genrcmi.o \
			$(LDIR)/libutils.a $(LDIR)/libhsm.a $(LDIR)/libaimUtil.a
			
	touch $(LDIR)/hsmAIM.so
	rm $(LDIR)/hsmAIM.so
	$(FCOMP) $(FSOFLG) -o $(LDIR)/hsmAIM.so \
		$(ODIR)/hsmAIM.o $(ODIR)/hsmUtils.o $(ODIR)/hsmAdj.o \
		$(ODIR)/rcm.o $(ODIR)/genrcmi.o \
		 -L$(LDIR) -lutils -laimUtil -locsm \
		-legads -ludunits2 -lhsm -ldl $(FRPATH) $(CPPSLB) $(ESPFLIBS)

$(ODIR)/hsmAIM.o:	hsmAIM.c $(IDIR)/capsTypes.h 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils hsmAIM.c \
		-o $(ODIR)/hsmAIM.o

$(ODIR)/hsmUtils.o:	hsmUtils.c $(IDIR)/capsTypes.h 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils hsmUtils.c \
		-o $(ODIR)/hsmUtils.o

$(ODIR)/hsmAdj.o:	hsmAdj.c $(IDIR)/capsTypes.h 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils hsmAdj.c \
		-o $(ODIR)/hsmAdj.o

$(ODIR)/rcm.o:	rcm/rcm.cpp
	$(CXX) -c $(CPPOPT) -Wno-comment rcm/rcm.cpp $(DEFINE) \
		-o $(ODIR)/rcm.o

$(ODIR)/genrcmi.o:	rcm/genrcmi.cpp
	$(CXX) -c $(CPPOPT) rcm/genrcmi.cpp $(DEFINE) \
		-o $(ODIR)/genrcmi.o

# Target for building HSM fortran files
HSMINC   = $(wildcard src/*.inc)
HSMSRC   = $(wildcard src/*.f)
HSMOBJS  = $(subst src,$(ODIR),$(HSMSRC:.f=.o))
$(HSMOBJS): $(HSMINC)
	$(MAKE) -C src $(subst $(ODIR)/,,$@)

$(LDIR)/libhsm.a: $(HSMOBJS)
	$(MAKE) -C src $@

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -bufferoverflowhigh

cleanAIM:
	-rm -f $(ODIR)/hsmAIM.o $(ODIR)/hsmUtils.o $(ODIR)/hsmAdj.o $(ODIR)/rcm.o $(ODIR)/genrcmi.o

cleanallAIM:	cleanAIM
	-rm -f $(LDIR)/hsmAIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen hsmAIM_dox.cfg; cd latex; make; mv refman.pdf hsmAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f hsmAIM.tag)

	
clean: cleanAIM
	make -C src clean
	
cleanall: cleanallAIM
	make -C src cleanall
