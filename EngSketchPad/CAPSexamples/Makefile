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

ifdef PYTHONINC
test: cCAPStest pyCAPStest
else
test: cCAPStest
endif

cCAPStest: all $(COPY)
	(cd $(TDIR)/regressionTest; ./execute_CTestRegression.sh $(TYPE))

ifdef PYTHONINC
pyCAPStest: all $(COPY)
	(cd $(TDIR)/regressionTest; ./execute_PyTestRegression.sh $(TYPE))
endif

ifdef ESP_BLOC

PYFILES=$(wildcard pyCAPS/*.py)
PYTEST=$(patsubst pyCAPS/%.py,$(TDIR)/pyCAPS/%.py,$(PYFILES))

CSMFILES=$(wildcard csmData/*.csm)
CSMFILES+=$(wildcard csmData/*.udc)
CSMTEST=$(patsubst csmData/%,$(TDIR)/csmData/%,$(CSMFILES))

REGFILES=$(wildcard regressionTest/*.*)
REGTEST=$(patsubst regressionTest/%,$(TDIR)/regressionTest/%,$(REGFILES))

copy: $(TDIR)/csmData $(TDIR)/cCAPS $(TDIR)/pyCAPS $(TDIR)/regressionTest \
	$(PYTEST) $(CSMTEST) $(REGTEST)

$(TDIR):
	mkdir -p $(TDIR)
$(TDIR)/cCAPS: | $(TDIR)
	mkdir -p $(TDIR)/cCAPS
$(TDIR)/pyCAPS: $(TDIR)
	mkdir -p $(TDIR)/pyCAPS
$(TDIR)/csmData: | $(TDIR)
	mkdir -p $(TDIR)/csmData
$(TDIR)/regressionTest: | $(TDIR)
	mkdir -p $(TDIR)/regressionTest
	
$(PYTEST): $(TDIR)/pyCAPS/%.py: pyCAPS/%.py | $(TDIR)/pyCAPS
	cp -p -f $< $@

$(CSMTEST): $(TDIR)/csmData/%: csmData/% | $(TDIR)/csmData
	cp -p -f $< $@

$(REGTEST): $(TDIR)/regressionTest/%: regressionTest/% | $(TDIR)/regressionTest
	cp -p -f $< $@

endif
