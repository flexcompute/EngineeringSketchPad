#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = .
endif


VPATH = $(ODIR)

$(TDIR)/ftest:	$(LDIR)/libwsserver.a $(ODIR)/test.o
	$(FCOMP) -o $(TDIR)/ftest $(ODIR)/test.o $(LDIR)/libwsserver.a \
		-lpthread -lz $(CPPSLB)

$(ODIR)/test.o:	test.f 
	$(FCOMP) -c $(FOPTS) -fno-range-check -fallow-argument-mismatch test.f \
		-I../include -o $(ODIR)/test.o

clean:
	-rm $(ODIR)/test.o
