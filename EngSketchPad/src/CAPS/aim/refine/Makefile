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
else
ODIR  = .
endif


$(LDIR)/refineAIM.so:	$(ODIR)/refineAIM.o $(ODIR)/libmeshb7.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	touch $(LDIR)/refineAIM.so
	rm $(LDIR)/refineAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/refineAIM.so $(ODIR)/refineAIM.o $(ODIR)/libmeshb7.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl \
		$(RPATH) -lm

$(ODIR)/refineAIM.o:	refineAIM.c $(IDIR)/capsTypes.h \
			../utils/meshUtils.h ../meshWriter/libmeshbWriter/libmeshbWriter.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I../meshWriter/libmeshbWriter refineAIM.c \
		-o $(ODIR)/refineAIM.o

$(ODIR)/libmeshb7.o:	../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
						../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.h
	$(CC) -c $(COPTS) $(DEFINE) \
		-I../meshWriter/libmeshbWriter/libMeshb/sources \
		../meshWriter/libmeshbWriter/libMeshb/sources/libmeshb7.c \
		-o $(ODIR)/libmeshb7.o

lint:
	$(LINT) -I../../include -I../meshWriter/libmeshbWriter -I../meshWriter/libmeshbWriter/libMeshb/sources -I$(IDIR) -I../utils $(wildcard *.c) -uniondef -exportlocal

clean:
	-rm $(ODIR)/refineAIM.o $(ODIR)/libmeshb7.o

cleanall:	clean
	-rm $(LDIR)/refineAIM.so

	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen refineAIM_dox.cfg; cd latex; make; mv refman.pdf refineAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f refineAIM.tag)
