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

$(TDIR)\interferenceTest.exe:    $(ODIR)\interferenceTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\interferenceTest.exe $(ODIR)\interferenceTest.obj \
		/link /LIBPATH:$(LDIR) caps.lib egads.lib

$(ODIR)\interferenceTest.obj:    interferenceTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) interferenceTest.c \
		/Fo$(ODIR)\interferenceTest.obj 
    
clean:
	-del $(ODIR)\interferenceTest.obj

cleanall:	clean
	-del $(TDIR)\interferenceTest.exe
