#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

ifdef ESP_BLOC
TDIR  = $(ESP_BLOC)/examples
COPY  = copy
else
TDIR  = .
COPY  =
endif

all: $(COPY) | $(TDIR) 
	(make -C cCAPS -f aeroelastic_SU2.make)
	(make -C cCAPS -f avlTest.make)
	(make -C cCAPS -f awaveTest.make)
	(make -C cCAPS -f frictionTest.make)
	(make -C cCAPS -f fun3d.make)
	(make -C cCAPS -f hsm.make)
	(make -C cCAPS -f interferenceTest.make)
	(make -C cCAPS -f msesTest.make)
	(make -C cCAPS -f mystran.make)
	(make -C cCAPS -f pointwiseTest.make)

clean:
	(make -C cCAPS -f aeroelastic_SU2.make clean)
	(make -C cCAPS -f avlTest.make clean)
	(make -C cCAPS -f awaveTest.make clean)
	(make -C cCAPS -f frictionTest.make clean)	
	(make -C cCAPS -f fun3d.make clean)
	(make -C cCAPS -f hsm.make clean)
	(make -C cCAPS -f interferenceTest.make clean)
	(make -C cCAPS -f msesTest.make clean)
	(make -C cCAPS -f mystran.make clean)
	(make -C cCAPS -f pointwiseTest.make clean)

cleanall:
	(make -C cCAPS -f aeroelastic_SU2.make cleanall)
	(make -C cCAPS -f avlTest.make cleanall)
	(make -C cCAPS -f awaveTest.make cleanall)
	(make -C cCAPS -f frictionTest.make cleanall)
	(make -C cCAPS -f fun3d.make cleanall)
	(make -C cCAPS -f hsm.make cleanall)
	(make -C cCAPS -f interferenceTest.make cleanall)
	(make -C cCAPS -f msesTest.make cleanall)
	(make -C cCAPS -f mystran.make cleanall)
	(make -C cCAPS -f pointwiseTest.make cleanall)

test: cCAPStest pyCAPStest

cCAPStest: all $(COPY)
	(cd $(TDIR)/regressionTest; ./execute_CTestRegression.sh $(TYPE))

pyCAPStest: all $(COPY)
	(cd $(TDIR)/regressionTest; ./execute_PyTestRegression.sh $(TYPE))

ifdef ESP_BLOC

PYFILES=$(wildcard pyCAPS/*.py)
PYTEST=$(patsubst pyCAPS/%.py,$(TDIR)/pyCAPS/%.py,$(PYFILES))

CSMFILES=$(wildcard csmData/*.csm)
CSMFILES+=$(wildcard csmData/*.udc)
CSMFILES+=$(wildcard csmData/x15)
CSMTEST=$(patsubst csmData/%,$(TDIR)/csmData/%,$(CSMFILES))

ANALYSISFILES=$(wildcard analysisData/nastran/*.*)
ANALYSISDATA=$(patsubst analysisData/%,$(TDIR)/analysisData/%,$(ANALYSISFILES))

CORSAIRFILES=$(wildcard corsairlite/*.p*)
CORSAIRFILES+=$(wildcard corsairlite/capsPhase/csm/*.csm)
CORSAIRFILES+=$(wildcard corsairlite/capsPhase/pyCAPS/*.py)
CORSAIRFILES+=$(wildcard corsairlite/multifidelity/csm/*.csm)
CORSAIRFILES+=$(wildcard corsairlite/multifidelity/pyCAPS/*.py)
CORSAIRTEST=$(patsubst corsairlite/%,$(TDIR)/corsairlite/%,$(CORSAIRFILES))

REGFILES=$(wildcard regressionTest/*.*)
REGTEST=$(patsubst regressionTest/%,$(TDIR)/regressionTest/%,$(REGFILES))

copy: $(TDIR)/analysisData $(TDIR)/csmData $(TDIR)/cCAPS $(TDIR)/pyCAPS $(TDIR)/regressionTest $(TDIR)/corsairlite \
	$(PYTEST) $(CSMTEST) $(ANALYSISDATA) $(CORSAIRTEST) $(REGTEST)

$(TDIR):
	mkdir -p $(TDIR)
$(TDIR)/cCAPS: | $(TDIR)
	mkdir -p $(TDIR)/cCAPS
$(TDIR)/pyCAPS: | $(TDIR)
	mkdir -p $(TDIR)/pyCAPS
$(TDIR)/analysisData: | $(TDIR)
	mkdir -p $(TDIR)/analysisData
$(TDIR)/analysisData/nastran: | $(TDIR)/analysisData
	mkdir -p $(TDIR)/analysisData/nastran
$(TDIR)/csmData: | $(TDIR)
	mkdir -p $(TDIR)/csmData
$(TDIR)/corsairlite: | $(TDIR)
	mkdir -p $(TDIR)/corsairlite
$(TDIR)/corsairlite/capsPhase/csm: | $(TDIR)
	mkdir -p $(TDIR)/corsairlite/capsPhase/csm
$(TDIR)/corsairlite/capsPhase/pyCAPS: | $(TDIR)
	mkdir -p $(TDIR)/corsairlite/capsPhase/pyCAPS
$(TDIR)/corsairlite/multifidelity/csm: | $(TDIR)
	mkdir -p $(TDIR)/corsairlite/multifidelity/csm
$(TDIR)/corsairlite/multifidelity/pyCAPS: | $(TDIR)
	mkdir -p $(TDIR)/corsairlite/multifidelity/pyCAPS
$(TDIR)/regressionTest: | $(TDIR)
	mkdir -p $(TDIR)/regressionTest
	
$(PYTEST): $(TDIR)/pyCAPS/%.py: pyCAPS/%.py | $(TDIR)/pyCAPS
	cp -p -f $< $@

$(ANALYSISDATA): $(TDIR)/analysisData/%: analysisData/% | $(TDIR)/analysisData $(TDIR)/analysisData/nastran
	cp -p -f $< $@

$(CSMTEST): $(TDIR)/csmData/%: csmData/% | $(TDIR)/csmData
	cp -r -p -f $< $@

$(CORSAIRTEST): $(TDIR)/corsairlite/%: corsairlite/% | $(TDIR)/corsairlite \
	$(TDIR)/corsairlite/capsPhase/csm \
	$(TDIR)/corsairlite/capsPhase/pyCAPS \
	$(TDIR)/corsairlite/multifidelity/csm \
	$(TDIR)/corsairlite/multifidelity/pyCAPS
	cp -p -f $< $@

$(REGTEST): $(TDIR)/regressionTest/%: regressionTest/% | $(TDIR)/regressionTest
	cp -p -f $< $@

endif
