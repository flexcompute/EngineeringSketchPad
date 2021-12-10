# Makefile for XDDM library

# $Id: Makefile,v 1.1.1.1 2013/06/04 17:52:20 mnemec Exp $

# MN, STC, June 2013

SHELL  = /bin/sh
CMD    = xddm
THEDATE    = `date "+%y.%m.%d"`
THEARCHIVE =  $(CMD)Lib--$(THEDATE)

AR      = ar 
ARFLAGS = ruv
CC      = gcc
LINK    = $(CC)

CCFLAGS = -Wall -g -pedantic

INSTALL = install
INSTALLFLAGS = -v

LIBDIR  = 
LIBS    = `xml2-config --libs`
INCLUDE = -I./ `xml2-config --cflags`

CMAIN   = main.c

CSOURCE = xddm.c 

CHEADER = 

CMAINOBJ = $(CMAIN:.c=.o)
COBJS    = $(CSOURCE:.c=.o) 

.c.o:
	$(CC)  $(CCFLAGS) $(INCLUDE) -c $<

default: all

all: $(CMD)

$(CMD): $(CMAINOBJ) $(COBJS)
	$(AR) $(ARFLAGS) lib$(CMD).a $(COBJS) $(FOBJS)
	$(LINK) $(LFLAGS) $(INCLUDE) -o xddm_tester $(CMAINOBJ) libxddm.a $(LIBDIR) $(LIBS)

install: $(CMD)_$(ARCH)
	 $(INSTALL) $(INSTALLFLAGS) -m 750 $(CMD) $(HOME)/bin/$(CMD)

clean:
	/bin/rm -f *.o *.a xddm_tester

new:
	make clean
	make all	

tgz:
	-tar cvzf $(THEARCHIVE).tgz --exclude-vcs --exclude-backups \
            Makefile *.c *.h examples 
