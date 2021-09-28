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

$(TDIR)\egads2cart.exe:	$(ODIR)\egads2cart.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\egads2cart.exe $(ODIR)\egads2cart.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\egads2cart.exe.manifest \
		/outputresource:$(TDIR)\egads2cart.exe;1

$(ODIR)\egads2cart.obj:	egads2cart.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egads2cart.c \
		/Fo$(ODIR)\egads2cart.obj

clean:
	-del $(ODIR)\egads2cart.obj

cleanall:	clean
	-del $(TDIR)\egads2cart.exe $(TDIR)\egads2cart.exe.manifest
