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

$(TDIR)\egadsGeom_dot.exe:	$(ODIR)\egadsGeom_dot.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\egadsGeom_dot.exe $(ODIR)\egadsGeom_dot.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\egadsGeom_dot.exe.manifest \
		/outputresource:$(TDIR)\egadsGeom_dot.exe;1

$(ODIR)\egadsGeom_dot.obj:	egadsGeom_dot.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egadsGeom_dot.c \
		/Fo$(ODIR)\egadsGeom_dot.obj

clean:
	-del $(ODIR)\egadsGeom_dot.obj

cleanall:	clean
	-del $(TDIR)\egadsGeom_dot.exe $(TDIR)\egadsGeom_dot.exe.manifest
