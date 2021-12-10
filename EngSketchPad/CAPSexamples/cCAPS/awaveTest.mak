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

$(TDIR)\awaveTest.exe:	$(ODIR)\awaveTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\awaveTest.exe $(ODIR)\awaveTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\awaveTest.obj:	awaveTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) awaveTest.c /Fo$(ODIR)\awaveTest.obj

clean:
	-del $(ODIR)\awaveTest.obj

cleanall:	clean
	-del $(TDIR)\awaveTest.exe
