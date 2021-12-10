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

$(TDIR)\multiContext.exe:	$(ODIR)\multiContext.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\multiContext.exe $(ODIR)\multiContext.obj \
		$(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\multiContext.exe.manifest \
		/outputresource:$(TDIR)\multiContext.exe;1

$(ODIR)\multiContext.obj:	multiContext.c $(IDIR)\egads.h \
		$(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) multiContext.c \
		/Fo$(ODIR)\multiContext.obj

clean:
	-del $(ODIR)\multiContext.obj

cleanall:	clean
	-del $(TDIR)\multiContext.exe $(TDIR)\multiContext.exe.manifest
