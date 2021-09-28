#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib

tire.exe:	tire.obj $(LDIR)\egads.lib
	cl tire.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest tire.exe.manifest /outputresource:tire.exe;1

tire.obj:	egads2tri.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) tire.c 

clean:
	-del tire.obj

cleanall:	clean
	-del tire.exe tire.exe.manifest
