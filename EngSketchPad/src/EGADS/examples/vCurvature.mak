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

$(TDIR)\vCurvature.exe:	$(ODIR)\vCurvature.obj $(LDIR)\egads.lib \
			$(LDIR)\wsserver.lib
	cl /Fe$(TDIR)\vCurvature.exe $(ODIR)\vCurvature.obj \
		/link /LIBPATH:$(LDIR) wsserver.lib egads.lib
	$(MCOMP) /manifest $(TDIR)\vCurvature.exe.manifest \
		/outputresource:$(TDIR)\vCurvature.exe;1

$(ODIR)\vCurvature.obj:	vCurvature.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers vCurvature.c \
		/Fo$(ODIR)\vCurvature.obj

clean:
	-del $(ODIR)\vCurvature.obj

cleanall:	clean
	-del $(TDIR)\vCurvature.exe $(TDIR)\vCurvature.exe.manifest
