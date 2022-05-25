#       Makefile for serveESP
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

BINLIST =	$(BDIR)/serveESP
TIMLIST =	$(LDIR)/viewer.so

default:	$(BINLIST) $(TIMLIST)

#
#	binaries
#
$(BDIR)/serveESP:	$(ODIR)/serveESP.o $(LDIR)/$(OSHLIB)
	$(CXX) -o $(BDIR)/serveESP $(ODIR)/serveESP.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl $(RPATH) 

$(ODIR)/serveESP.o:	serveESP.c OpenCSM.h common.h tim.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. serveESP.c -o $(ODIR)/serveESP.o

$(LDIR)/viewer.so:	$(ODIR)/timViewer.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/viewer.so
	rm $(LDIR)/viewer.so
	$(CXX) $(SOFLGS) -o $(LDIR)/viewer.so $(ODIR)/timViewer.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timViewer.o:	timViewer.c tim.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. timViewer.c -o $(ODIR)/timViewer.o

#
#	other targets
#
# Static analysis
SCANDIR=/tmp/scanOpenCSM
SCANEXCLUDE=
include $(IDIR)/STANALYZER.make

lint:
	@echo "Checking serveESP..."
	$(LINT) -I$(IDIR) serveESP.c OpenCSM.c udp.c tim.c -allocmismatch -duplicatequals -macrovarprefixexclude -exportlocal -mustfreefresh -mayaliasunique -kepttrans -immediatetrans
	@echo " "
	@echo "Checking viewer.so..."
	$(LINT) -I$(IDIR) timViewer.c

clean:
	(cd $(ODIR); rm -f serveESP.o timViewer.o )

cleanall:	clean
	rm -f  $(BINLIST)
