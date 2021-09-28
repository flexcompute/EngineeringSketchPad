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

$(TDIR)\egadsTopo_dot.exe:	$(ODIR)\egadsTopo_dot.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\egadsTopo_dot.exe $(ODIR)\egadsTopo_dot.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\egadsTopo_dot.exe.manifest \
		/outputresource:$(TDIR)\egadsTopo_dot.exe;1

$(ODIR)\egadsTopo_dot.obj:	egadsTopo_dot.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egadsTopo_dot.c \
		/Fo$(ODIR)\egadsTopo_dot.obj

clean:
	-del $(ODIR)\egadsTopo_dot.obj

cleanall:	clean
	-del $(TDIR)\egadsTopo_dot.exe $(TDIR)\egadsTopo_dot.exe.manifest
