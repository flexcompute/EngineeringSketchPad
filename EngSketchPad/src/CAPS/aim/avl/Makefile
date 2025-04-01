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
else
ODIR  = .
endif

# Allow implicit rule for object files in ODIR directory
VPATH = $(ODIR)

AVL_MRF_FILES=$(wildcard avlmrf/*.c)
AVL_MRF_OBJS=$(foreach src,$(AVL_MRF_FILES),$(ODIR)/$(basename $(notdir $(src))).o)

$(LDIR)/avlAIM.so:	$(ODIR)/avlAIM.o $(AVL_MRF_OBJS) $(LDIR)/libaimUtil.a $(LDIR)/libutils.a
	touch $(LDIR)/avlAIM.so
	rm $(LDIR)/avlAIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/avlAIM.so $(ODIR)/avlAIM.o $(AVL_MRF_OBJS) \
		-L$(LDIR) -laimUtil -lutils -locsm -legads -ludunits2 -ldl \
		$(CPPSLB) $(RPATH) -lm

$(ODIR)/avlAIM.o:	avlAIM.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils avlAIM.c -o $(ODIR)/avlAIM.o

$(ODIR)/%.o: avlmrf/%.c avlmrf/%.h avlmrf/read_util.h
	$(CC) -c $(COPTS) -Wno-unused-function -Iavlmrf $< -o $@

lint:
	$(LINT) -I../../include -I$(IDIR) -I../utils $(wildcard *.c) -exportlocal -nestcomment -bufferoverflowhigh -uniondef

clean:
	-rm $(ODIR)/avlAIM.o $(AVL_MRF_OBJS)

cleanall:	clean
	-rm $(LDIR)/avlAIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen avlAIM_dox.cfg; cd latex; make; mv refman.pdf avlAIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f avlAIM.tag)
