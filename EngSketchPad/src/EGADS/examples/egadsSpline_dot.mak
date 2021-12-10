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

$(TDIR)\egadsSpline_dot.exe:	$(ODIR)\egadsSpline_dot.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\egadsSpline_dot.exe $(ODIR)\egadsSpline_dot.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\egadsSpline_dot.exe.manifest \
		/outputresource:$(TDIR)\egadsSpline_dot.exe;1

$(ODIR)\egadsSpline_dot.obj:	egadsSpline_dot.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egadsSpline_dot.c \
		/Fo$(ODIR)\egadsSpline_dot.obj

clean:
	-del $(ODIR)\egadsSpline_dot.obj

cleanall:	clean
	-del $(TDIR)\egadsSpline_dot.exe $(TDIR)\egadsSpline_dot.exe.manifest
