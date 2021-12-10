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

$(TDIR)\vEffect.exe:	$(ODIR)\vEffect.obj $(LDIR)\egads.lib \
			$(LDIR)\wsserver.lib
	cl /Fe$(TDIR)\vEffect.exe $(ODIR)\vEffect.obj /link /LIBPATH:$(LDIR) \
		wsserver.lib egads.lib
	$(MCOMP) /manifest $(TDIR)\vEffect.exe.manifest \
		/outputresource:$(TDIR)\vEffect.exe;1

$(ODIR)\vEffect.obj:	vEffect.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers vEffect.c \
		/Fo$(ODIR)\vEffect.obj

clean:
	-del $(ODIR)\vEffect.obj

cleanall:	clean
	-del $(TDIR)\vEffect.exe $(TDIR)\vEffect.exe.manifest
