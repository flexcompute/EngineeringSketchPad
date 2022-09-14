#	Makefile for serveESP's CAPS-based TIMs
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

TIMLIST =	$(LDIR)\capsMode.dll \
		$(LDIR)\flowchart.dll \
		$(LDIR)\viewer.dll

default:	$(TIMLIST)

#
#	binaries
#
$(LDIR)\capsMode.dll:	$(ODIR)\timCapsMode.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\capsMode.dll $(LDIR)\capsMode.lib $(LDIR)\capsMode.exp
	link /out:$(LDIR)\capsMode.dll /dll /def:tim.def $(ODIR)\timCapsMode.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\capsMode.dll.manifest /outputresource:$(LDIR)\capsMode.dll;2

$(ODIR)\timCapsMode.obj:	timCapsMode.c tim.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timCapsMode.c /Fo$(ODIR)\timCapsMode.obj

$(LDIR)\flowchart.dll:	$(ODIR)\timFlowchart.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\flowchart.dll $(LDIR)\flowchart.lib $(LDIR)\flowchart.exp
	link /out:$(LDIR)\flowchart.dll /dll /def:tim.def $(ODIR)\timFlowchart.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\flowchart.dll.manifest /outputresource:$(LDIR)\flowchart.dll;2

$(ODIR)\timFlowchart.obj:	timFlowchart.c tim.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timFlowchart.c /Fo$(ODIR)\timFlowchart.obj

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
	-del $(ODIR)\timCapsMode.obj
	-del $(ODIR)\timFlowchart.obj
	-del $(ODIR)\timViewer.obj

cleanall:	clean
	-del $(LDIR)\capsMode.dll   $(LDIR)\capsMode.lib   $(LDIR)\capsMode.exp
	-del $(LDIR)\flowchart.dll  $(LDIR)\flowchart.lib  $(LDIR)\flowchart.exp
	-del $(LDIR)\viewer.dll     $(LDIR)\viewer.lib     $(LDIR)\viewer.exp
