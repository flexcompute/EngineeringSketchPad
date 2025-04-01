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
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif


$(LDIR)/metrisAIM.so:	$(ODIR)/metrisAIM.o $(ODIR)/libmeshb7.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	touch $(LDIR)/metrisAIM.so
	rm $(LDIR)/metrisAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/metrisAIM.so $(ODIR)/metrisAIM.o $(ODIR)/libmeshb7.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl \
		$(RPATH) -lm

$(ODIR)/metrisAIM.o:	metrisAIM.c $(IDIR)/capsTypes.h \
			../utils/meshUtils.h ../meshWriter/libmeshbWriter/libmeshbWriter.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I../meshWriter/libmeshbWriter metrisAIM.c \
		-o $(ODIR)/metrisAIM.o

$(ODIR)/libmeshb7.o:	../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
						../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.h
	$(CC) -c $(COPTS) $(DEFINE) \
		-I../meshWriter/libmeshbWriter/libMeshb/sources \
		../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
		-o $(ODIR)/libmeshb7.o

lint:
	$(LINT) -I../../include -I../meshWriter/libmeshbWriter -I../meshWriter/libmeshbWriter/libMeshb/sources -I$(IDIR) -I../utils $(wildcard *.c) -uniondef -exportlocal

clean:
	-rm $(ODIR)/metrisAIM.o $(ODIR)/libmeshb7.o

cleanall:	clean
	-rm $(LDIR)/metrisAIM.so

	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen metrisAIM_dox.cfg; cd latex; make; mv refman.pdf metrisAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f metrisAIM.tag)
