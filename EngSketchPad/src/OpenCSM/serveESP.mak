#	Makefile for serveESP
#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib
BDIR  = $(ESP_ROOT)\bin
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
!ELSE
ODIR  = .
!ENDIF

BINLIST = 	$(BDIR)\serveESP.exe
TIMLIST =	$(LDIR)\viewer.dll

default:	$(BINLIST) $(TIMLIST)

#
#	binaries
#
$(BDIR)\serveESP.exe:	$(ODIR)\serveESP.obj $(LDIR)\ocsm.dll
	cl /Fe$(BDIR)\serveESP.exe $(ODIR)\serveESP.obj /link /LIBPATH:$(LDIR) caps.lib ocsm.lib wsserver.lib egads.lib
        $(MCOMP) /manifest $(BDIR)\serveESP.exe.manifest               /outputresource:$(BDIR)\serveESP.exe;1

$(ODIR)\serveESP.obj:	OpenCSM.h common.h serveESP.c tim.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) /I$(IDIR)\winhelpers serveESP.c /Fo$(ODIR)\serveESP.obj

$(LDIR)\viewer.dll:	$(ODIR)\timViewer.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\viewer.dll $(LDIR)\viewer.lib $(LDIR)\viewer.exp
	link /out:$(LDIR)\viewer.dll /dll /def:tim.def $(ODIR)\timViewer.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\viewer.dll.manifest /outputresource:$(LDIR)\viewer.dll;2

$(ODIR)\timViewer.obj:	timViewer.c tim.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timViewer.c /Fo$(ODIR)\timViewer.obj

#
#	other targets
#
clean:
	-del $(ODIR)\serveESP.obj
	-del $(ODIR)\timViewer.obj

cleanall:	clean
	-del $(BDIR)\serveESP.exe   $(BDIR)\sensCSM.exe    $(BDIR)\*.manifest
	-del $(LDIR)\viewer.dll     $(LDIR)\viewer.lib     $(LDIR)\viewer.exp
