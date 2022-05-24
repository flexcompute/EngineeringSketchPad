#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\examples\cCAPS
!ELSE
ODIR  = .
TDIR  = .
!ENDIF

$(TDIR)\msesTest.exe:	$(ODIR)\msesTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\msesTest.exe $(ODIR)\msesTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\msesTest.obj:	msesTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) msesTest.c /Fo$(ODIR)\msesTest.obj

clean:
	-del $(ODIR)\msesTest.obj

cleanall:	clean
	-del $(TDIR)\msesTest.exe
