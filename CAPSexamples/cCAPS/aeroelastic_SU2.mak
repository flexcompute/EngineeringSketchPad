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

$(TDIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.exe:	$(ODIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.exe $(ODIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.obj:	aeroelasticSimple_Iterative_SU2_and_MystranTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) aeroelasticSimple_Iterative_SU2_and_MystranTest.c /Fo$(ODIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.obj

clean:
	-del $(ODIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.obj

cleanall:	clean
	-del $(TDIR)\aeroelasticSimple_Iterative_SU2_and_MystranTest.exe
