ifneq (, $(shell which julia))
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
egads_types_ctts.jl: $(ESP_ROOT)/include/egadsErrors.h $(ESP_ROOT)/include/egadsTypes.h
	@echo " Creating file with egads constants"
	@julia --eval 'include("build.jl")'
endif

.PHONY: cleanall
cleanall:
	rm -f egads_types_ctts.jl
	@echo "removed constants file: egads_types_ctts.jl"

