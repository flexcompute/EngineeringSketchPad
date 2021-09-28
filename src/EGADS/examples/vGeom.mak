#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
SDIR  = $(MAKEDIR)
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF

$(TDIR)\vGeom.exe:	$(ODIR)\vGeom.obj $(LDIR)\egads.lib \
			$(LDIR)\wsserver.lib $(LDIR)\z.lib
	cl /Fe$(TDIR)\vGeom.exe $(ODIR)\vGeom.obj /link /LIBPATH:$(LDIR) \
		wsserver.lib z.lib egads.lib ws2_32.lib
	$(MCOMP) /manifest $(TDIR)\vGeom.exe.manifest \
		/outputresource:$(TDIR)\vGeom.exe;1

$(ODIR)\vGeom.obj:	vGeom.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers vGeom.c \
		/Fo$(ODIR)\vGeom.obj

clean:
	-del $(ODIR)\vGeom.obj

cleanall:	clean
	-del $(TDIR)\vGeom.exe $(TDIR)\vGeom.exe.manifest
