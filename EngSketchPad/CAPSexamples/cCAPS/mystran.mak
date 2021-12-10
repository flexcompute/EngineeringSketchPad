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

$(TDIR)\mystranTest.exe:    $(ODIR)\mystranTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)/mystranTest.exe $(ODIR)\mystranTest.obj \
		/link /LIBPATH:$(LDIR) caps.lib egads.lib

$(ODIR)\mystranTest.obj:    mystranTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) mystranTest.c /Fo$(ODIR)\mystranTest.obj 
    
clean:
	-del $(ODIR)\mystranTest.obj

cleanall:	clean
	-del $(TDIR)\mystranTest.exe
