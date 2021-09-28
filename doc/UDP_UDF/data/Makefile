#	Makefile
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib

default:	tire tire.so

tire:	tire.o
	$(CC) -o tire tire.o -L$(LDIR) -legads -lm $(RPATH)

tire.o:	udpTire.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) udpTire.c -o tire.o

tire.so:	udpTire.o
	touch tire.so
	rm tire.so
	$(CC) $(SOFLGS) -o tire.so udpTire.o -L$(LDIR) -legads -lm $(RPATH)

udpTire.o:	udpTire.c
	$(CC) -c $(COPTS) $(DEFINE) -DUDP -I$(IDIR) udpTire.c

clean:
	-rm udpTire.o tire.o

cleanall:	clean
	-rm tire tire.so
