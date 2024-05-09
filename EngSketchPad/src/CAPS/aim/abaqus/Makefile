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

ifdef PYTHONINC

#Adding the rpath is needed when compiling with ESP_Python
PYTHONLPATH=$(filter -L%,$(PYTHONLIB))
PYTHONRPATH=$(PYTHONLPATH:-L%=-Wl,-rpath %)

#Extract the python include directory
PYTHONIPATH0=$(word 1,$(filter -I%,$(PYTHONINC)))
ifeq ("$(PYTHONIPATH0)","")
PYTHONIPATH=$(PYTHONINC)
else
PYTHONIPATH=$(PYTHONIPATH0:-I%=%)
endif

$(LDIR)/abaqusAIM.so: $(ODIR)/abaqusAIM.o $(ODIR)/abaqusUtils.o $(LDIR)/libaimUtil.a \
			$(LDIR)/libutils.a $(ODIR)/abaqusFILReader.o
	echo "Building Abaqus with Python"
	$(CC) $(SOFLGS) -o $(LDIR)/abaqusAIM.so $(ODIR)/abaqusAIM.o $(ODIR)/abaqusUtils.o $(ODIR)/abaqusFILReader.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 $(PYTHONLIB) -lm -ldl $(PYTHONRPATH)
		
else 
$(LDIR)/abaqusAIM.so: $(ODIR)/abaqusAIM.o $(ODIR)/abaqusUtils.o $(LDIR)/libaimUtil.a \
			$(LDIR)/libutils.a
	echo "Building Abaqus without Python"
	$(CC) $(SOFLGS) -o $(LDIR)/abaqusAIM.so $(ODIR)/abaqusAIM.o $(ODIR)/abaqusUtils.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -lm -ldl

endif 

$(ODIR)/abaqusAIM.o: abaqusAIM.c abaqusUtils.h ../utils/meshUtils.h \
		../utils/meshUtils.h \
		../utils/miscUtils.h \
		../utils/feaUtils.h \
		$(IDIR)/capsTypes.h $(IDIR)/aimUtil.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils \
		abaqusAIM.c -o $(ODIR)/abaqusAIM.o
		
ifdef PYTHONINC
$(ODIR)/abaqusUtils.o: abaqusUtils.c abaqusUtils.h cython/abaqusFILReader.h $(IDIR)/capsTypes.h $(PYTHONIPATH)/Python.h
	$(CC) -c $(COPTS) $(DEFINE) -DHAVE_PYTHON -I$(PYTHONINC) -Icython -I$(IDIR) -I../utils \
		abaqusUtils.c -o $(ODIR)/abaqusUtils.o

$(ODIR)/abaqusFILReader.o: cython/abaqusFILReader.c cython/abaqusFILReader.h $(PYTHONIPATH)/Python.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils \
		-I$(PYTHONINC) cython/abaqusFILReader.c -o $(ODIR)/abaqusFILReader.o
else 
$(ODIR)/abaqusUtils.o: abaqusUtils.c abaqusUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils \
		abaqusUtils.c -o $(ODIR)/abaqusUtils.o
endif 

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils abaqusAIM.c abaqusUtils.c -exportlocal -nestcomment -evalorder -retvalother -boolops -bufferoverflowhigh -uniondef
		
clean:
	-rm $(ODIR)/abaqusAIM.o $(ODIR)/abaqusUtils.o
ifdef PYTHONINC
	-rm $(ODIR)/abaqusFILReader.o
endif


cleanall:	clean
	-rm $(LDIR)/abaqusAIM.so

.PHONY: cython
cython:
	make -C cython

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen abaqusAIM_dox.cfg; cd latex; make; mv refman.pdf abaqusAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f abaqusAIM.tag)
