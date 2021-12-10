#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF

$(TDIR)\blend.exe:	$(ODIR)\blend.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\blend.exe $(ODIR)\blend.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\blend.exe.manifest \
		/outputresource:$(TDIR)\blend.exe;1

$(ODIR)\blend.obj:	blend.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) blend.c \
		/Fo$(ODIR)\blend.obj

clean:
	-del $(ODIR)\blend.obj

cleanall:	clean
	-del $(TDIR)\blend.exe $(TDIR)\blend.exe.manifest
