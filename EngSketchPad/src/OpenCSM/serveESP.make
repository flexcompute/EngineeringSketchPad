#       Makefile for serveESP's CAPS-based TIMs
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

TIMLIST =	$(LDIR)/capsMode.so \
		$(LDIR)/flowchart.so \
		$(LDIR)/viewer.so

default:	$(TIMLIST)

#
#	binaries
#
$(LDIR)/capsMode.so:	$(ODIR)/timCapsMode.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/capsMode.so
	rm $(LDIR)/capsMode.so
	$(CXX) $(SOFLGS) -o $(LDIR)/capsMode.so $(ODIR)/timCapsMode.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timCapsMode.o:	timCapsMode.c tim.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. timCapsMode.c -o $(ODIR)/timCapsMode.o

$(LDIR)/flowchart.so:	$(ODIR)/timFlowchart.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/flowchart.so
	rm $(LDIR)/flowchart.so
	$(CXX) $(SOFLGS) -o $(LDIR)/flowchart.so $(ODIR)/timFlowchart.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timFlowchart.o:	timFlowchart.c tim.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. timFlowchart.c -o $(ODIR)/timFlowchart.o

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
	@echo "Checking capsMode.so..."
	$(LINT) -I$(IDIR) timCapsMode.c
	@echo " "
	@echo "Checking flowchart.so..."
	$(LINT) -I$(IDIR) timFlowchart.c
	@echo " "
	@echo "Checking viewer.so..."
	$(LINT) -I$(IDIR) timViewer.c

clean:
	(cd $(ODIR); rm -f timCapsMode.o timFlowchart.o timViewer.o )

cleanall:	clean
	rm -f  $(TIMLIST)
