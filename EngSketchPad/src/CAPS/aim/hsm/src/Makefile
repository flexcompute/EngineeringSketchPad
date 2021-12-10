#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
endif

VPATH = $(ODIR)

ifeq ($(EFCOMP),ifort)
FFLAG = -real-size 64 -check bounds 
NORECURS = -norecursive
else
FFLAG = -fdefault-real-8
FFLAG += -fbounds-check -ffpe-trap=invalid,zero
NORECURS = -fno-recursive
endif
FOBJS = hsmsol.o hsmdep.o hsmout.o hsmgeo.o hsmabd.o hsmeqn.o hsmbc.o \
	hsmprj.o hsmren.o hsmre1.o hsmrfm.o hsmbb2.o hsmglr.o \
	ludcmp.o sbsolve.o atanc.o bmdump.o cross.o

HSMINC = $(wildcard *.inc)

$(TDIR)/hsmrun:		$(LDIR)/libhsm.a hsmrun.o
	$(FCOMP) $(FOPTS) $(FFLAG) $(NORECURS) $(ODIR)/hsmrun.o -L$(LDIR) -lhsm \
		-o $(TDIR)/hsmrun 
	
hsmrun.o:	hsmrun.f
	$(FCOMP) -c $(FOPTS) $(FFLAG) $(NORECURS) hsmrun.f -o $(ODIR)/hsmrun.o
	
	
$(LDIR)/libhsm.a:	$(FOBJS)
	touch $(LDIR)/libhsm.a
	rm $(LDIR)/libhsm.a
	(cd $(ODIR); ar -rs $(LDIR)/libhsm.a $(FOBJS) )

$(FOBJS): %.o:	%.f $(HSMINC)
	$(FCOMP) -c $(FOPTS) $(FFLAG) $(ESPFFLAGS) $< -o $(ODIR)/$@

clean:
	(cd $(ODIR); rm -f $(FOBJS) hsmrun.o)

cleanall:	clean
	-rm -f $(TDIR)/hsmrun $(LDIR)/libhsm.a
