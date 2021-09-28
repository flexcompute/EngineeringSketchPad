#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib

egads2tri.exe:	egads2tri.obj $(LDIR)\egads.lib
	cl egads2tri.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest egads2tri.exe.manifest /outputresource:egads2tri.exe;1

egads2tri.obj:	egads2tri.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) egads2tri.c 

clean:
	-del egads2tri.obj

cleanall:	clean
	-del egads2tri.exe egads2tri.exe.manifest
