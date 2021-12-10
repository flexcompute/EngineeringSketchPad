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

# locations for openssl (MAC) & certificate directory
OSSINC = /usr/local/Cellar/openssl@1.1/1.1.1l_1/include
OSSLIB = /usr/local/Cellar/openssl@1.1/1.1.1l_1/lib
CRTDIR = -DLWS_OPENSSL_CLIENT_CERTS=\"@clientcertdir@\"

ifeq ("$(ESP_ARCH)","LINUX64")
default:	$(TDIR)/server $(LDIR)/libwsserver.so
else
default:	$(TDIR)/server $(LDIR)/libwsserver.dylib
endif

VPATH = $(ODIR)

OBJS  = base64-decode.o handshake.o client-handshake.o libwebsockets.o \
        extension-deflate-stream.o md5.o extension-x-google-mux.o parsers.o \
        extension.o sha-1.o browserMessage.o server.o wv.o

$(TDIR)/server:	$(LDIR)/libwsserver.a $(ODIR)/servertest.o
	$(CXX) -o $(TDIR)/server $(ODIR)/servertest.o $(LDIR)/libwsserver.a \
		-lpthread -lz -L$(OSSLIB) -lssl -lcrypto -lm

$(ODIR)/servertest.o:	server.c 
	$(CC) -c $(COPTS) -DLWS_OPENSSL_SUPPORT -DLWS_NO_FORK -DSTANDALONE \
		server.c -I$(IDIR) -I$(OSSINC) -I. -o $(ODIR)/servertest.o

$(LDIR)/libwsserver.a:	map.o $(OBJS) $(ODIR)/fwv.o
	touch $(LDIR)/libwsserver.a
	rm $(LDIR)/libwsserver.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libwsserver.a map.o $(OBJS) \
		$(ODIR)/fwv.o; $(RANLB) )

$(LDIR)/libwsserver.so:	$(OBJS) map.o
	touch $(LDIR)/libwsserver.so
	rm $(LDIR)/libwsserver.so
	(cd $(ODIR); $(CXX) -shared -Wl,-no-undefined \
		-o $(LDIR)/libwsserver.so $(OBJS) map.o -lz -lpthread \
		-L$(OSSLIB) -lssl -lcrypto -lm)

$(LDIR)/libwsserver.dylib: $(OBJS) map.o
	touch $(LDIR)/libwsserver.dylib
	rm $(LDIR)/libwsserver.dylib
	(cd $(ODIR); $(CXX) -dynamiclib -o $(LDIR)/libwsserver.dylib \
		$(OBJS) map.o -lz -L$(OSSLIB) -lssl -lcrypto -lm \
		-undefined error -install_name '@rpath/libwsserver.dylib' \
		-compatibility_version $(CASREV) \
		-current_version $(EGREV) )

$(OBJS):	extension-deflate-stream.h libwebsockets.h \
		extension-x-google-mux.h private-libwebsockets.h \
		$(IDIR)/wsserver.h $(IDIR)/wsss.h
.c.o:
	$(CC) -c $(COPTS) -DLWS_OPENSSL_SUPPORT -DLWS_NO_FORK $(CRTDIR) \
		-I$(IDIR) -I$(OSSINC) $< -o $(ODIR)/$@

$(ODIR)/fwv.o:	fwv.c
	$(CC) -c $(COPTS) -DLWS_NO_FORK -I$(IDIR) fwv.c -o $(ODIR)/fwv.o

map.o:		map.cpp
	$(CXX) -c $(CPPOPT) map.cpp -o $(ODIR)/map.o

lint:
	@echo "Selective analysis:"
	$(LINT) server.c wv.c fwv.c -I../include -I. -DLWS_OPENSSL_SUPPORT -DLWS_NO_FORK -DSTANDALONE -uniondef -evalorder -boolops -exitarg -nullpass -exportlocal -kepttrans -fullinitblock +charint -syntax
	@echo " "
	@echo "Full analysis:"
	$(LINT) server.c wv.c fwv.c libwebsockets.c sha-1.c base64-decode.c extension.c extension-deflate-stream.c extension-x-google-mux.c md5.c parsers.c handshake.c client-handshake.c -I../include -I. -DLWS_OPENSSL_SUPPORT -DLWS_NO_FORK -DSTANDALONE -boolops -exitarg -fullinitblock -nullpass -mayaliasunique +charint -exportlocal -evalorder -shiftimplementation -bufferoverflowhigh -casebreak -kepttrans -retvalother -uniondef -immediatetrans -observertrans -allocmismatch -syntax

clean:
	-(cd $(ODIR); rm map.o $(OBJS) servertest.o )

cleanall:	clean
	-rm $(TDIR)/server $(LDIR)/libwsserver.a
ifeq ("$(ESP_ARCH)","LINUX64")
	-rm $(LDIR)/libwsserver.so
else
	-rm $(LDIR)/libwsserver.dylib
endif
