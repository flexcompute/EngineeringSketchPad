#
include ../include/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
IDIR  = $(ESP_ROOT)/include
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

# allow for clang static analyzer
ifneq (,$(findstring analyzer,$(CC)))
DEFINE += -DDEBUG
endif
ifneq (,$(findstring intercept,$(CC)))
DEFINE += -DDEBUG
endif

VPATH = $(ODIR)

OBJS  = liteBase.o liteMemory.o liteGeom.o liteTopo.o liteAttrs.o liteImport.o \
	liteString.o

ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libegadslite.dylib $(LDIR)/libegadslitestatic.a
else
ifeq ($(ESP_ARCH),DARWINM1)
default:	$(LDIR)/libegadslite.dylib $(LDIR)/libegadslitestatic.a
else
default:	$(LDIR)/libegadslite.so    $(LDIR)/libegadslitestatic.a
endif
endif

$(LDIR)/libegadslite.so:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o liteEffect.o \
				liteUVmap.o
	touch $(LDIR)/libegadslite.so
	rm $(LDIR)/libegadslite.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined \
		-o $(LDIR)/libegadslite.so $(OBJS) liteTess.o liteTris.o \
		liteQuads.o liteTessInp.o egadsRobust.o emp.o evaLite.o \
		ratLite.o liteRegQuads.o liteEffect.o liteUVmap.o \
		-lpthread -lm )

$(LDIR)/libegadslite.dylib:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o liteEffect.o \
				liteUVmap.o
	touch $(LDIR)/libegadslite.dylib
	rm $(LDIR)/libegadslite.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegadslite.dylib \
                $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobust.o emp.o evaLite.o ratLite.o liteRegQuads.o \
		liteEffect.o liteUVmap.o \
		-undefined error -install_name '@rpath/libegadslite.dylib' \
                -current_version $(EGREV) )

$(LDIR)/libegadslitestatic.a:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o liteEffect.o \
				liteUVmap.o
	touch $(LDIR)/libegadslitestatic.a
	rm $(LDIR)/libegadslitestatic.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libegadslitestatic.a $(OBJS) \
                liteTess.o liteTris.o liteQuads.o liteTessInp.o egadsRobust.o \
		emp.o evaLite.o ratLite.o liteRegQuads.o liteEffect.o \
		liteUVmap.o ; $(RANLB) )

$(OBJS): %.o:	%.c ../include/egadsErrors.h ../src/egadsInternals.h \
		../include/egadsTypes.h liteClasses.h liteDevice.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I. -I../src -I../util/uvmap \
		$< -o $(ODIR)/$@

liteTess.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsTess.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTess.c -o $(ODIR)/liteTess.o

liteTris.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsTris.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTris.c -o $(ODIR)/liteTris.o

liteQuads.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsQuads.c
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I. -I../src \
                ../src/egadsQuads.c -o $(ODIR)/liteQuads.o

liteTessInp.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsTessInp.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                -I../util ../src/egadsTessInp.c -o $(ODIR)/liteTessInp.o

liteEffect.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsEffect.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsEffect.c -o $(ODIR)/liteEffect.o

liteUVmap.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h liteUVmap.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                -I../util/uvmap liteUVmap.c -o $(ODIR)/liteUVmap.o

emp.o:		../util/emp.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../util/emp.c \
		-o $(ODIR)/emp.o

evaLite.o:	../util/evaluate.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
		../util/evaluate.c -o $(ODIR)/evaLite.o

ratLite.o:	../util/rational.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE ../util/rational.c \
		-o $(ODIR)/ratLite.o

liteRegQuads.o:	../util/regQuads.c ../util/regQuads.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I../util \
		../util/regQuads.c -o $(ODIR)/liteRegQuads.o

egadsRobust.o:	../src/egadsRobust.c
	$(CC) -c $(COPTS) $(DEFINE) ../src/egadsRobust.c \
		-o $(ODIR)/egadsRobust.o

# Static analysis
SCANDIR=/tmp/scanEGADSlite
SCANEXCLUDE=
include ../include/STANALYZER.make

lint:
	$(LINT) -I$(IDIR) -DLITE -I../src -I. -I../util/uvmap liteAttrs.c liteGeom.c liteMemory.c liteBase.c liteUVmap.c liteImport.c liteString.c liteTest.c liteTopo.c ../util/evaluate.c ../util/rational.c ../util/regQuads.c ../src/egadsEffect.c -evalorder -exportlocal -shiftimplementation -aliasunique -retvalother -uniondef -formatcode

clean:
	-(cd $(ODIR); rm $(OBJS) egadsRobust.o emp.o liteTess.o liteTris.o liteQuads.o liteTessInp.o evaLite.o ratLite.o liteRegQuads.o liteEffect.o liteUVmap.o )

cleanall:	clean
	touch $(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.so
	-rm $(LDIR)/libegadslite*
