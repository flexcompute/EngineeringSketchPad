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

$(TDIR)\pointwiseTest.exe:	$(ODIR)\pointwiseTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\pointwiseTest.exe $(ODIR)\pointwiseTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\pointwiseTest.obj:	pointwiseTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) pointwiseTest.c /Fo$(ODIR)\pointwiseTest.obj

clean:
	-del $(ODIR)\pointwiseTest.obj

cleanall:	clean
	-del $(TDIR)\pointwiseTest.exe
