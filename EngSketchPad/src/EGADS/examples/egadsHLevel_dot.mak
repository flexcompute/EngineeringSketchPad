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

$(TDIR)\egadsHLevel_dot.exe:	$(ODIR)\egadsHLevel_dot.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\egadsHLevel_dot.exe $(ODIR)\egadsHLevel_dot.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\egadsHLevel_dot.exe.manifest \
		/outputresource:$(TDIR)\egadsHLevel_dot.exe;1

$(ODIR)\egadsHLevel_dot.obj:	egadsHLevel_dot.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egadsHLevel_dot.c \
		/Fo$(ODIR)\egadsHLevel_dot.obj

clean:
	-del $(ODIR)\egadsHLevel_dot.obj

cleanall:	clean
	-del $(TDIR)\egadsHLevel_dot.exe $(TDIR)\egadsHLevel_dot.exe.manifest
