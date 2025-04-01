# $Id: Makefile-XDDM,v 1.3 2022/11/10 19:26:21 mnemec Exp $

# Makefile for XDDM library
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
CMD        = libxddm.a
TESTER     = test_xddm
THEDATE    = `date "+%y.%m.%d"`
THEARCHIVE =  $(CMD)--$(THEDATE)

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

CMAIN   = main.c

CSOURCE = xddm.c \
          xddm_read.c \
          xddm_write.c

CHEADER = xddm.h

CMAINOBJ = $(CMAIN:.c=.o)
COBJS    = $(CSOURCE:.c=.o) 

.c.o:
	$(CC)  $(CCFLAGS) $(INCLUDE) -c $<

default: $(CMD)

all: $(CMD) test

$(CMD): $(CMAINOBJ) $(COBJS)
	$(AR) $(ARFLAGS) $(CMD) $(COBJS)

$(TESTER): $(CMD)

test: $(TESTER)

$(TESTER): $(CMAINOBJ)
	$(LINK) $(LFLAGS) $(INCLUDE) -o $@ $(CMAINOBJ) $(CMD) $(LIBS)

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
	-tar cvzf $(THEARCHIVE).tgz Makefile *.c *.h
	-chmod 640 $(THEARCHIVE).tgz
