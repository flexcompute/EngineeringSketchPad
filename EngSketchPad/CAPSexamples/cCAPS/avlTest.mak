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

$(TDIR)\avlTest.exe:	$(ODIR)\avlTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\avlTest.exe $(ODIR)\avlTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\avlTest.obj:	avlTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) avlTest.c /Fo$(ODIR)\avlTest.obj

clean:
	-del $(ODIR)\avlTest.obj

cleanall:	clean
	-del $(TDIR)\avlTest.exe
