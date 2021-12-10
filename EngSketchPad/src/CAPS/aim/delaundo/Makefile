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

 
$(LDIR)/delaundoAIM.so:	$(ODIR)/delaundoAIM.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	touch $(LDIR)/delaundoAIM.so
	rm $(LDIR)/delaundoAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/delaundoAIM.so $(ODIR)/delaundoAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl \
		$(RPATH) -lm

$(ODIR)/delaundoAIM.o:	delaundoAIM.c $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h  $(IDIR)/aimMesh.h ../utils/meshUtils.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils delaundoAIM.c \
		-o $(ODIR)/delaundoAIM.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -bufferoverflowhigh

clean:
	-rm $(ODIR)/delaundoAIM.o

cleanall:	clean
	-rm $(LDIR)/delaundoAIM.so

	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen delaundoAIM_dox.cfg; cd latex; make; mv refman.pdf delaundoAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f delaundoAIM.tag)
