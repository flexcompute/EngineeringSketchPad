#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF

!IF [ifort /help > ifort.txt]
default:	$(TDIR)\testc.exe $(TDIR)\parsec.exe
!ELSE
default:	$(TDIR)\testc.exe $(TDIR)\testf.exe $(TDIR)\parsec.exe \
		$(TDIR)\parsef.exe
!ENDIF

$(TDIR)\testc.exe:	$(ODIR)\testc.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\testc.exe $(ODIR)\testc.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\testc.exe.manifest \
		/outputresource:$(TDIR)\testc.exe;1

$(ODIR)\testc.obj:	testc.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) testc.c /Fo$(ODIR)\testc.obj

$(TDIR)\testf.exe:	$(ODIR)\testf.obj $(LDIR)\fgads.lib $(LDIR)\egads.lib
	ifort /traceback /Fe$(TDIR)\testf.exe $(ODIR)\testf.obj $(LIBPTH) \
		fgads.lib egads.lib $(FLIBS)
#	cl /MD /Fe$(TDIR)\testf.exe $(TDIR)\testf.obj $(LIBPTH) fgads.lib \
#		egads.lib libifcoremd.lib \
#		libifportmd.lib libmmd.lib /entry:mainCRTStartup
	$(MCOMP) /manifest $(TDIR)\testf.exe.manifest \
		/outputresource:$(TDIR)\testf.exe;1
	-del ifort.txt

$(ODIR)\testf.obj:	testf.f
	ifort /c $(FOPTS) -I$(IDIR) testf.f /Fo$(ODIR)\testf.obj

$(TDIR)\parsec.exe:	$(ODIR)\parsec.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\parsec.exe $(ODIR)\parsec.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\parsec.exe.manifest \
		/outputresource:$(TDIR)\parsec.exe;1

$(ODIR)\parsec.obj:	parsec.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) parsec.c /Fo$(ODIR)\parsec.obj

$(TDIR)\parsef.exe:	$(ODIR)\parsef.obj $(LDIR)\fgads.lib $(LDIR)\egads.lib
	ifort /traceback /Fe$(TDIR)\parsef.exe $(ODIR)\parsef.obj $(LIBPTH) \
		fgads.lib egads.lib $(FLIBS)
	$(MCOMP) /manifest $(TDIR)\parsef.exe.manifest \
		/outputresource:$(TDIR)\parsef.exe;1
	-del ifort.txt

$(ODIR)\parsef.obj:	parsef.f
	ifort /c $(FOPTS) -I$(IDIR) parsef.f /Fo$(ODIR)\parsef.obj

clean:
	-del $(ODIR)\testc.obj $(ODIR)\testf.obj
	-del $(ODIR)\parsec.obj $(ODIR)\parsef.obj

cleanall:	clean
	-del $(TDIR)\testc.exe $(TDIR)\testc.exe.manifest
	-del $(TDIR)\testf.exe $(TDIR)\testf.exe.manifest
	-del $(TDIR)\parsec.exe $(TDIR)\parsec.exe.manifest
	-del $(TDIR)\parsef.exe $(TDIR)\parsef.exe.manifest
