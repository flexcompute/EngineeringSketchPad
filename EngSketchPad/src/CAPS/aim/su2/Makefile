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

OBJS =	$(ODIR)/su2AIM.o $(ODIR)/su2Utils.o $(ODIR)/su2_4_Cardinal.o \
	$(ODIR)/su2_5_Raven.o $(ODIR)/su2_6_Falcon.o $(ODIR)/su2_7_Blackbird.o


$(LDIR)/su2AIM.so:	$(OBJS) $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/su2AIM.so $(OBJS) \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/su2AIM.o:	su2AIM.c su2Utils.h $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../meshWriter/su2Writer -I../utils su2AIM.c \
		-o $(ODIR)/su2AIM.o

$(ODIR)/su2Utils.o:	su2Utils.c su2Utils.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils su2Utils.c \
		-o $(ODIR)/su2Utils.o

$(ODIR)/su2_4_Cardinal.o:	su2_4_Cardinal.c su2Utils.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils su2_4_Cardinal.c \
		-o $(ODIR)/su2_4_Cardinal.o

$(ODIR)/su2_5_Raven.o:	su2_5_Raven.c su2Utils.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils su2_5_Raven.c \
		-o $(ODIR)/su2_5_Raven.o

$(ODIR)/su2_6_Falcon.o:	su2_6_Falcon.c su2Utils.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils su2_6_Falcon.c \
		-o $(ODIR)/su2_6_Falcon.o

$(ODIR)/su2_7_Blackbird.o:	su2_7_Blackbird.c su2Utils.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils su2_7_Blackbird.c \
		-o $(ODIR)/su2_7_Blackbird.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils -I../meshWriter/su2Writer $(wildcard *.c) -exportlocal -uniondef -nestcomment -boolops

clean:
	-rm $(OBJS)

cleanall:	clean
	-rm $(LDIR)/su2AIM.so

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/su2_and_Tetgen_PyTest.py \
				../../examples/csmData/cfdMultiBody.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen su2AIM_dox.cfg; cd latex; make; mv refman.pdf su2AIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f su2AIM.tag)
