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


$(LDIR)/astrosAIM.so:	$(ODIR)/astrosAIM.o $(ODIR)/astrosUtils.o $(ODIR)/astrosCards.o $(LDIR)/libutils.a $(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/astrosAIM.so $(ODIR)/astrosAIM.o  $(ODIR)/astrosUtils.o $(ODIR)/astrosCards.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm $(CPPSLB) $(PYTHONLIB)

$(ODIR)/astrosAIM.o:	astrosAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils astrosAIM.c \
		-o $(ODIR)/astrosAIM.o
		
$(ODIR)/astrosUtils.o:	astrosUtils.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils astrosUtils.c \
		-o $(ODIR)/astrosUtils.o
         
$(ODIR)/astrosCards.o:  astrosCards.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils astrosCards.c \
        -o $(ODIR)/astrosCards.o

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -boolops -bufferoverflowhigh -evalorder -uniondef -predboolothers

clean:
	-rm $(ODIR)/astrosAIM.o $(ODIR)/astrosUtils.o $(ODIR)/astrosCards.o

cleanall:	clean
	-rm $(LDIR)/astrosAIM.so
	
dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen astrosAIM_dox.cfg; cd latex; make; mv refman.pdf astrosAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f astrosAIM.tag)
