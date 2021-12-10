#	Makefile
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib

$(LDIR)/edgeAttr.so:	udfEdgeAttr.o
	touch $(LDIR)/edgeAttr.so
	rm $(LDIR)/edgeAttr.so
	$(CC) $(SOFLGS) -o $(LDIR)/edgeAttr.so udfEdgeAttr.o -L$(LDIR) \
		-legads -lm $(RPATH)

udfEdgeAttr.o:	udfEdgeAttr.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) udfEdgeAttr.c

clean:
	-rm udfEdgeAttr.o

cleanall:	clean
	-rm $(LDIR)/edgeAttr.so
