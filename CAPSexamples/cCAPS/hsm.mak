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

hsm:  $(TDIR)\hsmTest.exe $(TDIR)\hsmSimplePlateTest.exe $(TDIR)\hsmCantileverPlateTest.exe \
	$(TDIR)\hsmJoinedPlateTest.exe

$(TDIR)\hsmTest.exe:	$(ODIR)\hsmTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\hsmTest.exe $(ODIR)\hsmTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\hsmTest.obj:	hsmTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) hsmTest.c /Fo$(ODIR)\hsmTest.obj


$(TDIR)\hsmSimplePlateTest.exe:	$(ODIR)\hsmSimplePlateTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\hsmSimplePlateTest.exe $(ODIR)\hsmSimplePlateTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\hsmSimplePlateTest.obj:	hsmSimplePlateTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) hsmSimplePlateTest.c /Fo$(ODIR)\hsmSimplePlateTest.obj


$(TDIR)\hsmCantileverPlateTest.exe:	$(ODIR)\hsmCantileverPlateTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\hsmCantileverPlateTest.exe $(ODIR)\hsmCantileverPlateTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\hsmCantileverPlateTest.obj:	hsmCantileverPlateTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) hsmCantileverPlateTest.c /Fo$(ODIR)\hsmCantileverPlateTest.obj


$(TDIR)\hsmJoinedPlateTest.exe:	$(ODIR)\hsmJoinedPlateTest.obj $(LDIR)\caps.lib
	cl /Fe$(TDIR)\hsmJoinedPlateTest.exe $(ODIR)\hsmJoinedPlateTest.obj $(LIBPTH) caps.lib \
		egads.lib

$(ODIR)\hsmJoinedPlateTest.obj:	hsmJoinedPlateTest.c $(IDIR)\caps.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) hsmJoinedPlateTest.c /Fo$(ODIR)\hsmJoinedPlateTest.obj


clean:
	-del $(ODIR)\hsmTest.obj $(ODIR)\hsmSimplePlateTest.obj $(ODIR)\hsmCantileverPlateTest.obj \
		$(ODIR)\hsmJoinedPlateTest.obj

cleanall:	clean
	-del $(TDIR)\hsmTest.exe $(TDIR)\hsmSimplePlateTest.exe $(TDIR)\hsmCantileverPlateTest.exe \
		$(TDIR)\hsmJoinedPlateTest.exe
