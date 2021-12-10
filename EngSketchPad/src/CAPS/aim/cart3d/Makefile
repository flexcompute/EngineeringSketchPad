#	Makefile for ESP/XDDM and Cart3D AIM
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
#ifndef CART3DLIB
#XDDM = $(CART3DLIB)
#endif
ifneq ($(wildcard xddm*),)
XDDM = .
else
XDDM = $(CART3DLIB)
endif
ifndef CART3DLIB
CART3DLIB = $(IDIR)
endif
XMLIB = $(shell xml2-config --libs)
# due to MAC problem on Xcode 8 & OSX 10.11
ifeq ($(ESP_ARCH),DARWIN64)
XMLIB = -lxml2
endif


ifeq ("$(XMLIB)","")
default:
	@echo "XML developer libarary missing. Not compiling cart3dAIM."
else
default:	$(LDIR)/cart3dAIM.so $(TDIR)/cart3dTest $(BDIR)/ESPxddm
endif


$(BDIR)/ESPxddm:	$(ODIR)/ESPxddm.o $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddm.o \
			$(ODIR)/bodyTess.o
	$(CC) -o $(BDIR)/ESPxddm $(ODIR)/ESPxddm.o $(ODIR)/writeTrix.o \
		$(ODIR)/surfTrix.o $(ODIR)/xddm.o $(ODIR)/bodyTess.o -L$(LDIR) -locsm -legads \
		-L$(CART3DLIB)/$(CART_ARCH) -lc3dio -lCart3D \
		-lpthread -ldl $(FRPATH) $(XMLIB) -lm

$(LDIR)/cart3dAIM.so:	$(ODIR)/cart3dAIM.o $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddm.o \
			$(ODIR)/bodyTess.o $(LDIR)/libaimUtil.a $(LDIR)/libutils.a
	$(CC) $(SOFLGS) -o $(LDIR)/cart3dAIM.so $(ODIR)/cart3dAIM.o \
		$(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddm.o $(ODIR)/bodyTess.o -L$(LDIR) -laimUtil \
		-locsm -legads -ludunits2 -lutils -L$(CART3DLIB)/$(CART_ARCH) \
		-lc3dio -lCart3D -ldl $(XMLIB) -lm

$(ODIR)/cart3dAIM.o:	cart3dAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -DCAPS_MESHUTILS -I$(IDIR) -I../utils -I$(CART3DLIB)/libCart3D -I$(XDDM)/xddm \
		cart3dAIM.c -o $(ODIR)/cart3dAIM.o

$(ODIR)/ESPxddm.o:	ESPxddm.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(XDDM)/xddm \
		-I. `xml2-config --cflags` ESPxddm.c -o $(ODIR)/ESPxddm.o

$(ODIR)/writeTrix.o:	writeTrix.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/c3dio \
		-I$(CART3DLIB)/libCart3D -I$(XDDM)/xddm \
		writeTrix.c -o $(ODIR)/writeTrix.o

$(ODIR)/surfTrix.o:	surfTrix.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/c3dio \
		-I$(CART3DLIB)/libCart3D \
		surfTrix.c -o $(ODIR)/surfTrix.o

$(ODIR)/xddm.o:	$(XDDM)/xddm/xddm.c $(XDDM)/xddm/xddm.h
	$(CC) -c $(COPTS) $(DEFINE) -pedantic -I$(XDDM)/xddm -I. \
		`xml2-config --cflags` $(XDDM)/xddm/xddm.c \
		 -o $(ODIR)/xddm.o

$(ODIR)/bodyTess.o:	bodyTess.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) bodyTess.c -o $(ODIR)/bodyTess.o

$(TDIR)/cart3dTest:	$(ODIR)/cart3dTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/cart3dTest $(ODIR)/cart3dTest.o \
		-L$(LDIR) -lcaps -locsm -legads $(RPATH) -lm -ldl

$(ODIR)/cart3dTest.o:	cart3dTest.c $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) cart3dTest.c \
		-o $(ODIR)/cart3dTest.o

lint:
	@echo "Cart3D AIM:"
	$(LINT) -I../../include -I$(IDIR) -I../utils -I$(XDDM)/xddm -I$(CART3DLIB)/libCart3D -I$(CART3DLIB)/c3dio cart3dTest.c cart3dAIM.c -incondefs -exportlocal -nestcomment -uniondef -redef
	@echo " "
	@echo "ESPxddm:"
	$(LINT) -I../../include -I$(IDIR) -I../utils -I$(XDDM)/xddm -I$(CART3DLIB)/libCart3D -I$(CART3DLIB)/c3dio -Dxddm_free=free ESPxddm.c bodyTess.c -incondefs -exportlocal -nestcomment -uniondef

clean:
	-rm -f $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddm.o $(ODIR)/bodyTess.o \
		$(ODIR)/ESPxddm.o $(ODIR)/cart3dAIM.o $(ODIR)/cart3dTest.o

cleanall:	clean
	-rm -f $(BDIR)/ESPxddm $(LDIR)/cart3dAIM.so $(TDIR)/cart3dTest
	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen cart3dAIM_dox.cfg; cd latex; make; mv refman.pdf cart3dAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f cart3dAIM.tag)	
