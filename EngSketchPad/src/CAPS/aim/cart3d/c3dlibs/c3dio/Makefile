#  $Id: Makefile-C3DIO,v 1.2 2022/11/10 19:22:55 mnemec Exp $

# Makefile for C3DIO library
#
#         % make
#         % make test
#         % make install
#         % make new
#         % make all
#         % make clean
#         % make distclean
#         % make tgz
#
#----------------------------------------------------------------

SHELL      = /bin/sh
CMD        = libc3dio.a
TESTER     = test_c3dio
THEDATE    = `date "+%y.%m.%d"`
THEARCHIVE = $(CMD)--$(THEDATE)

AR      = ar
ARFLAGS = ruv
CC      = gcc
LINK    = $(CC)

# debug
CCFLAGS = -g -pedantic -Wall -fsanitize=address
LFLAGS  = -fsanitize=address

# default
CCFLAGS =
LFLAGS  =

INSTALL = install
INSTALLFLAGS = -v
INSTALLDIR = ..

LIBS = `xml2-config --libs`
INCLUDE =  -I./ `xml2-config --cflags`

CMAIN   =

CSOURCE = c3dio_lib.c \
          c3d_vtk_trix.c \
          geomStructures.c \
          memory_util.c

CHEADER = CaseInfo.h \
          ConvInfo.h \
          memory_util.h \
          limiters.h \
          c3d_global.h \
          basicTypes.h \
          cartCells.h \
          geomStructures.h \
          geomTypes.h \
          gridStructures.h \
          int64.h \
          dualtime.h \
          IOinfo.h \
          MGinfo.h \
          SolverInfo.h \
          PostProcInfo.h \
          GridInfo.h \
          sensors.h \
          c3dio_lib.h

CMAINOBJ = $(CMAIN:.c=.o)
COBJS    = $(CSOURCE:.c=.o)

TESTMAIN = test_c3dio.c
TESTMAINOBJ = $(TESTMAIN:.c=.o)

.c.o:
	$(CC)  $(CCFLAGS) $(INCLUDE) -c $<

default: $(CMD)

all: $(CMD) test

$(CMD): $(COBJS)
	$(AR) $(ARFLAGS) $@ $(COBJS)

$(TESTER): $(CMD)

test: $(TESTER)

$(TESTER): $(TESTMAINOBJ) $(CMD)
	$(LINK) $(LFLAGS) $(INCLUDE) -o $@ $(TESTMAINOBJ) $(CMD) $(LIBS)

install: all
	$(INSTALL) -d $(INSTALLDIR)/lib
	$(INSTALL) -d $(INSTALLDIR)/include
	$(INSTALL) -d $(INSTALLDIR)/bin
	$(INSTALL) $(INSTALLFLAGS) -m 644 $(CMD)     $(INSTALLDIR)/lib
	$(INSTALL) $(INSTALLFLAGS) -m 644 $(CHEADER) $(INSTALLDIR)/include
	$(INSTALL) $(INSTALLFLAGS) -m 744 $(TESTER)  $(INSTALLDIR)/bin

clean:
	-\rm -f *.o $(CMD) $(TESTER)

distclean: clean
	-\rm -rf $(INSTALLDIR)/lib $(INSTALLDIR)/bin $(INSTALLDIR)/include

new:
	make clean
	make all

tgz:
	-tar cvzf $(THEARCHIVE).tgz Makefile* *.c *.h
	-chmod 640 $(THEARCHIVE).tgz
