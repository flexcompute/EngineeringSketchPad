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

default:	$(TDIR)/vTesstatic $(TDIR)/vTess

$(TDIR)/vTess:	$(ODIR)/vTess.o $(ODIR)/retessFaces.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vTess $(ODIR)/vTess.o $(ODIR)/retessFaces.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(TDIR)/vTesstatic:	$(ODIR)/vTess.o $(ODIR)/retessFaces.o \
		$(LDIR)/libwsserver.a $(LDIR)/libegadstatic.a
	$(CXX) -o $(TDIR)/vTesstatic $(ODIR)/vTess.o $(ODIR)/retessFaces.o \
		-L$(LDIR) -lwsserver -legadstatic $(LIBPATH) $(LIBS) \
		-lpthread -lz $(RPATH) -lm

$(ODIR)/vTess.o:	vTess.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vTess.c -o $(ODIR)/vTess.o

$(ODIR)/retessFaces.o:	../util/retessFaces.c $(IDIR)/egads.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) ../util/retessFaces.c \
		-o $(ODIR)/retessFaces.o

clean:
	-rm $(ODIR)/vTess.o $(ODIR)/retessFaces.o

cleanall:	clean
	-rm $(TDIR)/vTess $(TDIR)/vTesstatic
