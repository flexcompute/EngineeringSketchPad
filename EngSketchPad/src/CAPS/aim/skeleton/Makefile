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


$(LDIR)/skeletonAIM.so:	$(ODIR)/skeletonAIM.o $(LDIR)/libutils.a \
		$(LDIR)/libaimUtil.a
	$(CC) $(SOFLGS) -o $(LDIR)/skeletonAIM.so $(ODIR)/skeletonAIM.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/skeletonAIM.o:	skeletonAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils skeletonAIM.c \
		-o $(ODIR)/skeletonAIM.o

test:
	$(VALGRIND_COMMAND) python -u case.py

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -uniondef -exportlocal

clean:
	-rm $(ODIR)/skeletonAIM.o 

cleanall:	clean
	-rm $(LDIR)/skeletonAIM.so
