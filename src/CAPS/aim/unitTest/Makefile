#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

ifdef ESP_BLOC
TDIR  = $(ESP_BLOC)/aimUnitTest
COPY  = copy
else
TDIR  = .
COPY  =
endif

all: $(COPY) | $(TDIR) 
#	(cd cCAPS; make -f test.make)

clean:
#	(cd cCAPS; make -f test.make clean)

cleanall:
#	(cd cCAPS; make -f test.make cleanall)

ifdef PYTHONINC
test: cCAPStest pyCAPStest
else
test: cCAPStest
endif

cCAPStest: all $(COPY)
#	(cd $(TDIR)/cCAPS; make test)

ifdef PYTHONINC
pyCAPStest: all $(COPY)
	(cd $(TDIR)/pyCAPS; make test)
endif

ifdef ESP_BLOC

PYFILES=$(wildcard pyCAPS/*.* pyCAPS/Makefile)
PYTEST=$(patsubst pyCAPS/%,$(TDIR)/pyCAPS/%,$(PYFILES))

CSMFILES=$(wildcard csmData/*.csm)
CSMTEST=$(patsubst csmData/%,$(TDIR)/csmData/%,$(CSMFILES))

DATAFILES=$(wildcard datafiles/*.*)
DATATEST=$(patsubst datafiles/%,$(TDIR)/datafiles/%,$(DATAFILES))

copy: $(PYTEST) $(CSMTEST) $(DATATEST) | $(TDIR)/cCAPS

$(TDIR):
	mkdir -p $(TDIR)
$(TDIR)/cCAPS: | $(TDIR)
	mkdir -p $(TDIR)/cCAPS
$(TDIR)/pyCAPS: | $(TDIR)
	mkdir -p $(TDIR)/pyCAPS
$(TDIR)/csmData: | $(TDIR)
	mkdir -p $(TDIR)/csmData
$(TDIR)/datafiles: | $(TDIR)
	mkdir -p $(TDIR)/datafiles
	
$(PYTEST): $(TDIR)/pyCAPS/%: pyCAPS/% | $(TDIR)/pyCAPS
	cp -p -f $< $@

$(CSMTEST): $(TDIR)/csmData/%.csm: csmData/%.csm | $(TDIR)/csmData
	cp -p -f $< $@
	
$(DATATEST): $(TDIR)/datafiles/%: datafiles/% | $(TDIR)/datafiles
	cp -p -f $< $@

endif
