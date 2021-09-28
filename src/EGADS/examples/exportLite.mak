#
IDIR = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!ELSE
ODIR = .
TDIR = $(ESP_ROOT)\bin
!ENDIF

$(TDIR)\exportLite.exe:	$(ODIR)\exportLite.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\exportLite.exe $(ODIR)\exportLite.obj $(LIBPTH) egads.lib

$(ODIR)\exportLite.obj:	..\src\egadsExport.c $(IDIR)\egads.h \
		$(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h 
	cl /c $(COPTS) $(DEFINE) /I$(IDIR) /I..\util\uvmap -DSTANDALONE \
		..\src\egadsExport.c /Fo$(ODIR)\exportLite.obj

clean:
	-del $(ODIR)\exportLite.obj

cleanall:	clean
	-del $(TDIR)\exportLite.exe
