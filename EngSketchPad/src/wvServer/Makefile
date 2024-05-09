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

ifneq (,$(findstring LINUX,$(ESP_ARCH)))
default:	$(TDIR)/server $(LDIR)/libwsserver.so
else
default:	$(TDIR)/server $(LDIR)/libwsserver.dylib
endif

VPATH = $(ODIR)

OBJS  = base64-decode.o handshake.o client-handshake.o libwebsockets.o \
        extension-deflate-stream.o md5.o extension-x-google-mux.o parsers.o \
        extension.o sha-1.o server.o wv.o browserMessage.o

$(TDIR)/server:	$(LDIR)/libwsserver.a $(ODIR)/servertest.o
	$(CXX) -o $(TDIR)/server $(ODIR)/servertest.o $(LDIR)/libwsserver.a  \
		-lpthread -lz -lm

$(ODIR)/servertest.o:	server.c 
	$(CC) -c $(COPTS) -DLWS_NO_FORK -DSTANDALONE server.c \
		-I$(IDIR) -I. -o $(ODIR)/servertest.o

$(LDIR)/libwsserver.a:	map.o $(OBJS) $(ODIR)/fwv.o
	touch $(LDIR)/libwsserver.a
	rm $(LDIR)/libwsserver.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libwsserver.a map.o $(OBJS) \
		$(ODIR)/fwv.o; $(RANLB) )

$(LDIR)/libwsserver.so:	$(OBJS) map.o
	touch $(LDIR)/libwsserver.so
	rm $(LDIR)/libwsserver.so
	(cd $(ODIR); $(CXX) -shared -Wl,-no-undefined \
		-o $(LDIR)/libwsserver.so $(OBJS) map.o -lz -lpthread)

$(LDIR)/libwsserver.dylib: $(OBJS) map.o
	touch $(LDIR)/libwsserver.dylib
	rm $(LDIR)/libwsserver.dylib
	(cd $(ODIR); $(CXX) -dynamiclib -o $(LDIR)/libwsserver.dylib \
		$(OBJS) map.o -lz \
		-install_name '@rpath/libwsserver.dylib' \
		-compatibility_version $(CASREV) \
		-current_version $(EGREV) )

$(OBJS):	extension-deflate-stream.h libwebsockets.h \
		extension-x-google-mux.h private-libwebsockets.h \
		$(IDIR)/wsserver.h $(IDIR)/wsss.h
.c.o:
	$(CC) -c $(COPTS) -DLWS_NO_FORK -I$(IDIR) $< -o $(ODIR)/$@

$(ODIR)/fwv.o:	fwv.c
	$(CC) -c $(COPTS) -DLWS_NO_FORK -I$(IDIR) fwv.c -o $(ODIR)/fwv.o

map.o:		map.cpp
	$(CXX) -c $(CPPOPT) map.cpp -o $(ODIR)/map.o

lint:
	@echo "Selective analysis:"
	$(LINT) server.c wv.c fwv.c -I../include -I. -DLWS_NO_FORK -DSTANDALONE -uniondef -evalorder -boolops -exitarg -nullpass -exportlocal -kepttrans -fullinitblock +charint -syntax
	@echo " "
	@echo "Full analysis:"
	$(LINT) server.c wv.c fwv.c libwebsockets.c sha-1.c base64-decode.c extension.c extension-deflate-stream.c extension-x-google-mux.c md5.c parsers.c handshake.c client-handshake.c browserMessage.c -I../include -I. -DLWS_NO_FORK -DSTANDALONE -boolops -exitarg -fullinitblock -nullpass -mayaliasunique +charint -exportlocal -evalorder -shiftimplementation -bufferoverflowhigh -casebreak -kepttrans -retvalother -uniondef -immediatetrans -observertrans -allocmismatch -syntax

clean:
	-(cd $(ODIR); rm map.o fwv.o $(OBJS) servertest.o )

cleanall:	clean
	-rm -f $(TDIR)/server $(LDIR)/libwsserver.a
ifneq (,$(findstring LINUX,$(ESP_ARCH)))
	-rm -f $(LDIR)/libwsserver.so
else
	-rm -f $(LDIR)/libwsserver.dylib
endif
