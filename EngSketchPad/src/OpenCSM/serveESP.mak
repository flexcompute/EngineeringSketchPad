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

!IFDEF PYTHONINC
TIMLIST =	$(LDIR)\capsMode.dll \
		$(LDIR)\flowchart.dll \
		$(LDIR)\pyscript.dll \
		$(LDIR)\viewer.dll
!ENDIF

default:	$(TIMLIST)

#
#	binaries
#
$(LDIR)\capsMode.dll:	$(ODIR)\timCapsMode.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\capsMode.dll $(LDIR)\capsMode.lib $(LDIR)\capsMode.exp
	link /out:$(LDIR)\capsMode.dll /dll /def:tim.def $(ODIR)\timCapsMode.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\capsMode.dll.manifest /outputresource:$(LDIR)\capsMode.dll;2

$(ODIR)\timCapsMode.obj:	timCapsMode.c OpenCSM.h tim.h $(IDIR)\emp.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timCapsMode.c /Fo$(ODIR)\timCapsMode.obj


$(LDIR)\flowchart.dll:	$(ODIR)\timFlowchart.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\flowchart.dll $(LDIR)\flowchart.lib $(LDIR)\flowchart.exp
	link /out:$(LDIR)\flowchart.dll /dll /def:tim.def $(ODIR)\timFlowchart.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\flowchart.dll.manifest /outputresource:$(LDIR)\flowchart.dll;2

$(ODIR)\timFlowchart.obj:	timFlowchart.c OpenCSM.h tim.h $(IDIR)\emp.h $(IDIR)\wsserver.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timFlowchart.c /Fo$(ODIR)\timFlowchart.obj


$(LDIR)\pyscript.dll:	$(ODIR)\timPyscript.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\pyscript.dll $(LDIR)\pyscript.lib $(LDIR)\pyscript.exp
	link /out:$(LDIR)\pyscript.dll /dll /def:timPyscript.def $(ODIR)\timPyscript.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib $(PYTHONLIB)
	$(MCOMP) /manifest $(LDIR)\pyscript.dll.manifest /outputresource:$(LDIR)\pyscript.dll;2

$(ODIR)\timPyscript.obj:	timPyscript.c OpenCSM.h tim.h $(IDIR)\emp.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) /I$(PYTHONINC) timPyscript.c /Fo$(ODIR)\timPyscript.obj


$(LDIR)\viewer.dll:	$(ODIR)\timViewer.obj $(LDIR)\ocsm.dll
	-del $(LDIR)\viewer.dll $(LDIR)\viewer.lib $(LDIR)\viewer.exp
	link /out:$(LDIR)\viewer.dll /dll /def:tim.def $(ODIR)\timViewer.obj $(LDIR)\caps.lib $(LDIR)\ocsm.lib $(LDIR)\wsserver.lib $(LDIR)\egads.lib
	$(MCOMP) /manifest $(LDIR)\viewer.dll.manifest /outputresource:$(LDIR)\viewer.dll;2

$(ODIR)\timViewer.obj:	timViewer.c OpenCSM.h tim.h $(IDIR)\emp.h $(IDIR)\wsserver.h $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) timViewer.c /Fo$(ODIR)\timViewer.obj

#
#	other targets
#
clean:
	-del $(ODIR)\timCapsMode.obj
	-del $(ODIR)\timFlowchart.obj
	-del $(ODIR)\timPyscript.obj
	-del $(ODIR)\timViewer.obj

cleanall:	clean
	-del $(LDIR)\capsMode.dll   $(LDIR)\capsMode.lib   $(LDIR)\capsMode.exp
	-del $(LDIR)\flowchart.dll  $(LDIR)\flowchart.lib  $(LDIR)\flowchart.exp
	-del $(LDIR)\pyscript.dll   $(LDIR)\pyscript.lib   $(LDIR)\pyscript.exp
	-del $(LDIR)\viewer.dll     $(LDIR)\viewer.lib     $(LDIR)\viewer.exp
