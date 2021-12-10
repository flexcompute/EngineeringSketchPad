#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

all:
	(cd EGADS/src; make)
	(cd EGADS/lite; make)
	(cd wvServer; make)
	(cd EGADS/examples; make -f egads2cart.make)
	(cd EGADS/examples; make -f vCurvature.make)
	(cd EGADS/examples; make -f vGeom.make)
	(cd EGADS/examples; make -f vTess.make)
	(cd OpenCSM; make)
	(cd Slugs; make)
	(cd CAPS; make)
	(cd CAPS/src; make)
	(cd ../CAPSexamples; make)
	(cd CAPS/aim; make)
	@echo " "
	@echo " *** Build Completed! ***"

clean:
	(cd CAPS/aim; make clean)
	(cd ../CAPSexamples; make clean)
	(cd CAPS/src; make clean)
	(cd Slugs; make clean)
	(cd OpenCSM; make clean)
	(cd wvServer; make clean)
	(cd EGADS/examples; make -f egads2cart.make clean)
	(cd EGADS/examples; make -f vGeom.make clean)
	(cd EGADS/examples; make -f vTess.make clean)
	(cd EGADS/examples; make -f vCurvature.make clean)
	(cd EGADS/lite; make clean)
	(cd EGADS/src; make clean)

cleanall:
	(cd CAPS/aim; make cleanall)
	(cd ../CAPSexamples; make cleanall)
	(cd CAPS/src; make cleanall)
	(cd Slugs; make cleanall)
	(cd OpenCSM; make cleanall)
	(cd wvServer; make cleanall)
	(cd EGADS/examples; make -f egads2cart.make cleanall)
	(cd EGADS/examples; make -f vGeom.make cleanall)
	(cd EGADS/examples; make -f vTess.make cleanall)
	(cd EGADS/examples; make -f vCurvature.make cleanall)
	(cd EGADS/lite; make cleanall)
	(cd EGADS/src; make cleanall)
