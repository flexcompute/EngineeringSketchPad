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

$(TDIR)\frictionTest.exe:	$(ODIR)\frictionTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\frictionTest.exe $(ODIR)\frictionTest.obj \
		$(LIBPTH) caps.lib egads.lib

$(ODIR)\frictionTest.obj:	frictionTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) frictionTest.c /Fo$(ODIR)\frictionTest.obj

clean:
	-del $(ODIR)\frictionTest.obj

cleanall:	clean
	-del $(TDIR)\frictionTest.exe
