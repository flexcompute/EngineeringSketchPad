#	Makefile for udfDumpByName
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
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
else
ODIR  = .
endif

$(LDIR)/dumpbyname.so:	$(ODIR)/udfDumpByName.o
	touch $(LDIR)/dumpbyname.so
	rm $(LDIR)/dumpbyname.so
	$(CXX) $(SOFLGS) -o $(LDIR)/dumpbyname.so $(ODIR)/udfDumpByName.o \
		-L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/udfDumpByName.o:	udfDumpByName.cpp $(IDIR)/udpUtilities.c $(IDIR)/udpUtilities.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) -I. \
		udfDumpByName.cpp -o $(ODIR)/udfDumpByName.o

lint:
#	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother

.PHONY: test
test : $(LDIR)/dumpbyname.so
	$(MAKE) -C test

clean:
	-rm $(ODIR)/udfDumpByName.o

cleanall:	clean
	-rm $(LDIR)/dumpbyname.so
