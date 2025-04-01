#	Makefile for ESP/XDDM and Cart3D AIM
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
ifndef CART3DLIB
CART3DLIB = c3dlibs
endif
XMLIB0 = $(shell xml2-config --libs)
XMLIB1 = $(subst -lm,,$(XMLIB0))
XMLIB  = $(subst -lpthread,,$(XMLIB1))
# due to MAC problem on Xcode 8 & OSX 10.11
ifeq ($(ESP_ARCH),DARWIN64)
XMLIB = -lxml2
endif

# Suppress ld: warning: could not create compact unwind
# Remove this when it's no longer necessary...
ifeq ($(ESP_ARCH),DARWIN_ARM64)
XMLIB += -Wl,-no_compact_unwind
endif

# Cart3D not availble on linux aarch64 yet
ifeq ($(ESP_ARCH),LINUX_AARCH64)
XMLIB= 
endif

# Cart3D not availble on cygwin
ifeq ("$(ESP_ARCH)","CYGWIN64")
XMLIB= 
endif

ifeq ("$(XMLIB)","")
default:
	@echo "XML developer libarary missing. Not compiling cart3dAIM."
else
default:	$(LDIR)/cart3dAIM.so $(TDIR)/cart3dTest $(BDIR)/ESPxddm

XMLINC = $(shell xml2-config --cflags)
endif

VPATH   += $(CART3DLIB)/c3dio
C3DIOSRC = c3dio_lib.c \
           c3d_vtk_trix.c \
           geomStructures.c \
           memory_util.c
C3DIOOBJ = $(addprefix $(ODIR)/,$(C3DIOSRC:.c=.o))

VPATH  += $(CART3DLIB)/xddm
XDDMSRC = xddm.c \
          xddm_read.c \
          xddm_write.c
XDDMOBJ = $(addprefix $(ODIR)/,$(XDDMSRC:.c=.o))

$(BDIR)/ESPxddm:	$(ODIR)/ESPxddm.o $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddm.o \
			$(ODIR)/bodyTess.o $(C3DIOOBJ) $(XDDMOBJ)
	$(CC) -o $(BDIR)/ESPxddm $(ODIR)/ESPxddm.o $(ODIR)/writeTrix.o \
		$(ODIR)/surfTrix.o $(ODIR)/bodyTess.o $(XDDMOBJ) $(C3DIOOBJ) \
		 -L$(LDIR) -locsm -legads \
		-lpthread -ldl $(FRPATH) $(XMLIB) -lm

$(LDIR)/cart3dAIM.so:	$(ODIR)/cart3dAIM.o $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddmUtils.o \
			$(ODIR)/bodyTess.o $(C3DIOOBJ) $(XDDMOBJ) $(LDIR)/libaimUtil.a $(LDIR)/libutils.a
	$(CC) $(SOFLGS) -o $(LDIR)/cart3dAIM.so $(ODIR)/cart3dAIM.o \
		$(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddmUtils.o $(ODIR)/bodyTess.o \
		$(C3DIOOBJ) $(XDDMOBJ) -L$(LDIR) -lutils -laimUtil \
		-locsm -legads -ludunits2 -ldl $(XMLIB) -lm

$(ODIR)/cart3dAIM.o:	cart3dAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -DCAPS_MESHUTILS -I$(IDIR) -I../utils -I$(CART3DLIB)/c3dio \
		-I$(CART3DLIB)/xddm \
		cart3dAIM.c -o $(ODIR)/cart3dAIM.o

$(ODIR)/ESPxddm.o:	ESPxddm.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/xddm $(XMLINC) \
		ESPxddm.c -o $(ODIR)/ESPxddm.o

$(ODIR)/xddmUtils.o:	xddmUtils.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/xddm \
		xddmUtils.c -o $(ODIR)/xddmUtils.o

$(ODIR)/writeTrix.o:	writeTrix.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/c3dio -I$(CART3DLIB)/xddm \
		writeTrix.c -o $(ODIR)/writeTrix.o

$(ODIR)/surfTrix.o:	surfTrix.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/c3dio \
		surfTrix.c -o $(ODIR)/surfTrix.o

$(ODIR)/bodyTess.o:	bodyTess.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) bodyTess.c -o $(ODIR)/bodyTess.o

$(TDIR)/cart3dTest:	$(ODIR)/cart3dTest.o $(LDIR)/$(CSHLIB)
	$(CC) -o $(TDIR)/cart3dTest $(ODIR)/cart3dTest.o \
		-L$(LDIR) -lcaps -locsm -legads $(RPATH) -lm -ldl

$(ODIR)/cart3dTest.o:	cart3dTest.c $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) cart3dTest.c \
		-o $(ODIR)/cart3dTest.o

$(ODIR)/%.o: %.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(CART3DLIB)/c3dio -I$(CART3DLIB)/xddm $(XMLINC) -w $< -o $@


lint:
	@echo "Cart3D AIM:"
	$(LINT) -I../../include -I$(IDIR) -I../utils -I$(CART3DLIB)/xddm -I$(CART3DLIB)/c3dio $(XMLINC) cart3dTest.c cart3dAIM.c -incondefs -exportlocal -nestcomment -uniondef -redef
	@echo " "
	@echo "ESPxddm:"
	$(LINT) -I../../include -I$(IDIR) -I../utils -I$(CART3DLIB)/xddm -I$(CART3DLIB)/c3dio $(XMLINC) -Dxddm_free=free ESPxddm.c bodyTess.c -incondefs -exportlocal -nestcomment -uniondef

clean:
	-rm -f $(ODIR)/writeTrix.o $(ODIR)/surfTrix.o $(ODIR)/xddmUtils.o $(ODIR)/bodyTess.o \
		$(ODIR)/ESPxddm.o $(ODIR)/cart3dAIM.o $(ODIR)/cart3dTest.o $(C3DIOOBJ) $(XDDMOBJ)

cleanall:	clean
	-rm -f $(BDIR)/ESPxddm $(LDIR)/cart3dAIM.so $(TDIR)/cart3dTest
	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen cart3dAIM_dox.cfg; cd latex; make; mv refman.pdf cart3dAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f cart3dAIM.tag)	
