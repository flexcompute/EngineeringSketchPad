#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

$(TDIR)/multiContext:	$(ODIR)/multiContext.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/multiContext $(ODIR)/multiContext.o -L$(LDIR) \
		-legads $(RPATH) -lm

$(ODIR)/multiContext.o:	multiContext.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) multiContext.c \
		-o $(ODIR)/multiContext.o

clean:
	-rm $(ODIR)/multiContext.o 

cleanall:	clean
	-rm $(TDIR)/multiContext
