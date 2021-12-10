#
!include ..\include\$(ESP_ARCH).$(MSVC)
IDIR = $(ESP_ROOT)\include
LDIR = $(ESP_ROOT)\lib
!ifdef ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!else
ODIR = .
TDIR = $(ESP_ROOT)\bin
!endif


$(TDIR)\patchTest.exe:	$(ODIR)\patchTest.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\patchTest.exe $(ODIR)\patchTest.obj \
		$(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\patchTest.exe.manifest \
		/outputresource:$(TDIR)\patchTest.exe;1

$(ODIR)\patchTest.obj:	egadsPatch.c $(IDIR)\egads.h \
		$(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -DSTANDALONE egadsPatch.c \
		/Fo$(ODIR)\patchTest.obj

clean:
	-del $(ODIR)\patchTest.obj

cleanall:	clean
	-del $(TDIR)\patchTest.exe
