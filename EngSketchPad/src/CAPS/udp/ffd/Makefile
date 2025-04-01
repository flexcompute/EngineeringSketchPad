#	Makefile for FFD
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

default:	$(LDIR)/FFDsplines.so $(LDIR)/FFDdeform.so

$(LDIR)/FFDsplines.so:	$(ODIR)/udpFFDsplines.o $(ODIR)/ffdUtil.o
	$(CXX) $(SOFLGS) -o $(LDIR)/FFDsplines.so $(ODIR)/udpFFDsplines.o $(ODIR)/ffdUtil.o \
		-L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/udpFFDsplines.o:	udpFFDsplines.cpp spline3d_eval.hpp \
		$(IDIR)/udpUtilities.c $(IDIR)/udpUtilities.h \
		$(IDIR)/egads.h $(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) udpFFDsplines.cpp -o $@

$(LDIR)/FFDdeform.so:	$(ODIR)/udpFFDdeform.o $(ODIR)/spline3d_eval.o $(ODIR)/ffdUtil.o
	$(CXX) $(SOFLGS) -o $(LDIR)/FFDdeform.so $(ODIR)/udpFFDdeform.o \
		$(ODIR)/spline3d_eval.o $(ODIR)/ffdUtil.o -L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/udpFFDdeform.o:	udpFFDdeform.cpp spline3d_eval.hpp \
		$(IDIR)/udpUtilities.c $(IDIR)/udpUtilities.h \
		$(IDIR)/egads.h $(IDIR)/egadsTypes.h $(IDIR)/egadsErrors.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) udpFFDdeform.cpp -o $@

$(ODIR)/ffdUtil.o:	ffdUtil.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) ffdUtil.c -o $@

$(ODIR)/spline3d_eval.o:	spline3d_eval.cpp $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) spline3d_eval.cpp -o $@

spline3d_eval_test:	spline3d_eval.cpp $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
		$(IDIR)/egadsErrors.h
	$(CXX) $(CPPOPT) $(DEFINE) -I$(IDIR) -DSTANDALONE spline3d_eval.cpp -o $@

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -evalorder -retvalother

.PHONY: test
test : $(LDIR)/FFDsplines.so $(LDIR)/FFDdeform.so
	$(MAKE) -C test

clean:
	-rm $(ODIR)/udpFFDsplines.o
	-rm $(ODIR)/udpFFDdeform.o
	-rm $(ODIR)/spline3d_eval.o
	-rm $(ODIR)/ffdUtil.o

cleanall:	clean
	-rm $(LDIR)/FFDsplines.so
	-rm $(LDIR)/FFDdeform.so
