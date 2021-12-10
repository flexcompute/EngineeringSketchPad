#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
SDIR  = $(MAKEDIR)
UDIR  = $(SDIR)\..\util
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF

$(TDIR)\vTess.exe:	$(ODIR)\vTess.obj $(ODIR)\retessFaces.obj \
			$(LDIR)\egads.lib $(LDIR)\wsserver.lib
	cl /Fe$(TDIR)\vTess.exe $(ODIR)\vTess.obj $(ODIR)\retessFaces.obj \
		/link /LIBPATH:$(LDIR) wsserver.lib egads.lib
	$(MCOMP) /manifest $(TDIR)\vTess.exe.manifest \
		/outputresource:$(TDIR)\vTess.exe;1

$(ODIR)\vTess.obj:	vTess.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers vTess.c \
		/Fo$(ODIR)\vTess.obj

$(ODIR)\retessFaces.obj:	$(UDIR)\retessFaces.c $(IDIR)\egads.h \
				$(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) $(UDIR)\retessFaces.c \
		/Fo$(ODIR)\retessFaces.obj

clean:
	-del $(ODIR)\retessFaces.obj $(ODIR)\vTess.obj

cleanall:	clean
	-del $(TDIR)\vTess.exe $(TDIR)\vTess.exe.manifest
