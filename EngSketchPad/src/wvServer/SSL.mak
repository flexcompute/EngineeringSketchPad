#
IDIR = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
SDIR = $(MAKEDIR)
LDIR = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!ELSE
ODIR = .
TDIR = .
!ENDIF

# locations for openssl & default certificate directory
OSSINC = c:\Users\Bob\lib\OpenSSL\include
OSSLIB = c:\Users\Bob\lib\OpenSSL\lib\VC
DEFINE = $(DEFINE) /D_UNICODE /DUNICODE /DLWS_OPENSSL_CLIENT_CERTS=\"@clientcertdir@\"

OBJS = gettimeofday.obj websock-w32.obj base64-decode.obj handshake.obj \
       client-handshake.obj libwebsockets.obj extension-deflate-stream.obj \
       md5.obj extension-x-google-mux.obj parsers.obj extension.obj sha-1.obj \
       browserMessage.obj server.obj wv.obj

!IFDEF ESP_BLOC
default:	start $(TDIR)\server.exe $(LDIR)\wsserver.dll end
!ELSE
default:	$(TDIR)\server.exe $(LDIR)\wsserver.dll
!ENDIF

start:
	cd $(ODIR)
	xcopy $(SDIR)\*.c           /Q /Y
	xcopy $(SDIR)\*.cpp         /Q /Y
	xcopy $(SDIR)\*.h           /Q /Y
	xcopy $(SDIR)\*.def         /Q /Y

$(TDIR)\server.exe:	$(ODIR)\servertest.obj $(LDIR)\wsservers.lib \
			$(LDIR)\z.lib
	cl /Fe$(TDIR)\server.exe $(ODIR)\servertest.obj $(LDIR)\wsservers.lib \
		$(LDIR)\z.lib $(OSSLIB)\libssl64MD.lib \
		$(OSSLIB)\libcrypto64MD.lib ws2_32.lib
	$(MCOMP) /manifest $(TDIR)\server.exe.manifest \
		/outputresource:$(TDIR)\server.exe;1

$(ODIR)\servertest.obj:	server.c
	cl /c $(COPTS) $(DEFINE) /DSTANDALONE /I$(IDIR) /I$(IDIR)\winhelpers \
		/I$(OSSINC) /I. /I$(SDIR)\win32helpers server.c \
		/Fo$(ODIR)\servertest.obj

$(LDIR)\wsservers.lib:	$(ODIR)\map.obj $(OBJS) $(ODIR)\fwv.obj
	-del $(LDIR)\wsservers.lib
	lib /out:$(LDIR)\wsservers.lib $(ODIR)\map.obj $(OBJS) $(ODIR)\fwv.obj

$(LDIR)\wsserver.dll:	$(ODIR)\map.obj $(OBJS)
	-del $(LDIR)\wsserver.dll $(LDIR)\wsserver.lib $(LDIR)\wsserver.exp
	link /out:$(LDIR)\wsserver.dll /dll /def:wsserver.def $(ODIR)\map.obj \
		$(OBJS) $(LDIR)\z.lib $(OSSLIB)\libssl64MD.lib \
		$(OSSLIB)\libcrypto64MD.lib ws2_32.lib

$(OBJS):	extension-deflate-stream.h libwebsockets.h $(IDIR)\wsserver.h \
		extension-x-google-mux.h private-libwebsockets.h \
		$(IDIR)\wsss.h
.c.obj:
	cl /c $(COPTS) $(DEFINE) /I. /I$(IDIR) /I$(IDIR)\winhelpers \
		/I$(OSSINC) /I$(SDIR)\win32helpers /I$(SDIR)\zlib $<

$(ODIR)\map.obj:	map.cpp
	cl /c $(CPPOPT) $(DEFINE) map.cpp /Fo$(ODIR)\map.obj

$(ODIR)\fwv.obj:	fwv.c
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) /I. /I$(IDIR)\winhelpers \
		/I$(SDIR)\win32helpers fwv.c /Fo$(ODIR)\fwv.obj

$(LDIR)\z.lib:
	cd $(SDIR)\zlib
	nmake -f NMakefile
!IFDEF ESP_BLOC
	cd $(ODIR)
!ELSE
	cd $(SDIR)
!ENDIF

end:
	-del *.c *.h
	cd $(SDIR)

clean:
	cd $(SDIR)\zlib
	nmake -f NMakefile clean
!IFDEF ESP_BLOC
	cd $(ODIR)
!ELSE
	cd $(SDIR)
!ENDIF
	-del map.obj $(OBJS)
	cd $(SDIR)

cleanall:	clean
	-del $(LDIR)\wsserver.dll $(LDIR)\wsserver.lib $(LDIR)\wsserver.exp
	-del $(LDIR)\wsservers.lib $(TDIR)\server.exe $(ODIR)\servertest.obj
