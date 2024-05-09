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

#Adding the rpath is needed when compiling with Python
ifdef PYTHONINC
TIMLIST    += $(LDIR)/pyscript.so
PYTHONLPATH = $(filter -L%,$(PYTHONLIB))
PYTHONRPATH = $(PYTHONLPATH:-L%=-Wl,-rpath %)
endif

default:	$(TIMLIST)

#
#	binaries
#
$(LDIR)/capsMode.so:	$(ODIR)/timCapsMode.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/capsMode.so
	rm $(LDIR)/capsMode.so
	$(CXX) $(SOFLGS) -o $(LDIR)/capsMode.so $(ODIR)/timCapsMode.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timCapsMode.o:	timCapsMode.c OpenCSM.h tim.h $(IDIR)/emp.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. timCapsMode.c -o $(ODIR)/timCapsMode.o


$(LDIR)/flowchart.so:	$(ODIR)/timFlowchart.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/flowchart.so
	rm $(LDIR)/flowchart.so
	$(CXX) $(SOFLGS) -o $(LDIR)/flowchart.so $(ODIR)/timFlowchart.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timFlowchart.o:	timFlowchart.c OpenCSM.h tim.h $(IDIR)/emp.h $(IDIR)/wsserver.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. timFlowchart.c -o $(ODIR)/timFlowchart.o


$(LDIR)/pyscript.so:	$(ODIR)/timPyscript.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/pyscript.so
	rm $(LDIR)/pyscript.so
	$(CXX) $(SOFLGS) -o $(LDIR)/pyscript.so $(ODIR)/timPyscript.o $(PYTHONLIB) -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl $(RPATH) $(PYTHONRPATH)

$(ODIR)/timPyscript.o:	timPyscript.c OpenCSM.h tim.h $(IDIR)/emp.h $(IDIR)/caps.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. -I$(PYTHONINC) timPyscript.c -o $(ODIR)/timPyscript.o


$(LDIR)/viewer.so:	$(ODIR)/timViewer.o $(LDIR)/$(OSHLIB)
	touch $(LDIR)/viewer.so
	rm $(LDIR)/viewer.so
	$(CXX) $(SOFLGS) -o $(LDIR)/viewer.so $(ODIR)/timViewer.o -L$(LDIR) -lcaps -locsm -lwsserver -legads -ldl -lm

$(ODIR)/timViewer.o:	timViewer.c OpenCSM.h tim.h $(IDIR)/emp.h $(IDIR)/wsserver.h $(IDIR)/caps.h
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
	(cd $(ODIR); rm -f timCapsMode.o timFlowchart.o timPyscript.o timViewer.o )

cleanall:	clean
	rm -f  $(TIMLIST)
