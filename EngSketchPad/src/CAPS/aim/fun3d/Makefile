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
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
endif

ifdef PYTHONINC
include $(ESP_ROOT)/include/ESPPYTHON.make

$(LDIR)/fun3dAIM.so:	$(ODIR)/fun3dAIM.o $(ODIR)/fun3dUtils.o $(ODIR)/fun3dNamelist.o \
			$(LDIR)/libaimUtil.a $(LDIR)/libutils.a
	@echo "Building FUN3D with Python"
	$(CC) $(SOFLGS) -o $(LDIR)/fun3dAIM.so $(ODIR)/fun3dAIM.o $(ODIR)/fun3dUtils.o \
		$(ODIR)/fun3dNamelist.o -L$(LDIR) -lutils -laimUtil -locsm \
		-legads -ludunits2 $(PYTHONLIBFLAGS) -ldl -lm
		

$(ODIR)/fun3dAIM.o: fun3dAIM.c fun3dUtils.h cython/fun3dNamelist.h $(IDIR)/capsTypes.h $(PYTHONIPATH)/Python.h
	$(CC) -c $(COPTS) $(DEFINE) -DHAVE_PYTHON -DCYTHON_PEP489_MULTI_PHASE_INIT=0 \
		-I$(IDIR) -I../utils -I../meshWriter/ugridWriter \
		-I$(PYTHONINC) -Icython fun3dAIM.c -o $(ODIR)/fun3dAIM.o

$(ODIR)/fun3dNamelist.o: cython/fun3dNamelist.c cython/fun3dNamelist.h $(PYTHONIPATH)/Python.h
	$(CC) -c $(COPTS) $(DEFINE) -DCYTHON_PEP489_MULTI_PHASE_INIT=0 \
		-I$(IDIR) -I../utils \
		-I$(PYTHONINC) cython/fun3dNamelist.c -o $(ODIR)/fun3dNamelist.o

else 

$(LDIR)/fun3dAIM.so:	$(ODIR)/fun3dAIM.o $(ODIR)/fun3dUtils.o $(LDIR)/libaimUtil.a \
			$(LDIR)/libutils.a
	@echo "Building FUN3D without Python"
	$(CC) $(SOFLGS) -o $(LDIR)/fun3dAIM.so $(ODIR)/fun3dAIM.o $(ODIR)/fun3dUtils.o \
		-L$(LDIR) -lutils -laimUtil -locsm -legads -ludunits2 -ldl -lm

$(ODIR)/fun3dAIM.o: fun3dAIM.c fun3dUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I../meshWriter/ugridWriter \
		fun3dAIM.c -o $(ODIR)/fun3dAIM.o
endif 

$(ODIR)/fun3dUtils.o: fun3dUtils.c fun3dUtils.h $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I../meshWriter/ugridWriter \
		fun3dUtils.c -o $(ODIR)/fun3dUtils.o


lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils -I../meshWriter/ugridWriter fun3dAIM.c fun3dUtils.c -exportlocal -uniondef -nestcomment -bufferoverflowhigh -boolops

clean:
	-rm -f $(ODIR)/fun3dAIM.o $(ODIR)/fun3dUtils.o
ifdef PYTHONINC
	-rm -f $(ODIR)/fun3dNamelist.o
endif

cleanall:	clean
	-rm -f $(LDIR)/fun3dAIM.so

# Run Cython - development only	
.PHONY:  cython
cython:
	make -C cython

dox: doc/html/index.html
doc/html/index.html: \
				../../examples/pyCAPS/fun3d_and_Tetgen_PyTest.py \
				../../examples/csmData/cfdMultiBody.csm \
				$(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen fun3dAIM_dox.cfg; cd latex; make; mv refman.pdf fun3dAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f fun3dAIM.tag)
