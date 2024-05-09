ifneq (, $(shell which julia))
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
.PHONY: dev
dev: ../egadscommon/egads_types_ctts.jl
ifneq ($(shell julia --eval 'using egads' 2> /dev/null; echo $$?), 0)
	@echo "Dev egads : add package to Julia dependency"; 
	@julia --eval 'using Pkg; Pkg.develop(path=".")';
endif

.PHONY: test
test: dev
	@echo "Test egads"
	@$(VALGRIND_COMMAND) julia --eval 'using Pkg; Pkg.test("egads")'

.PHONY: cleanall
cleanall: clean
	@echo "rm egads"
	@julia --eval 'using Pkg; try Pkg.rm("egads") catch throw @info " egads not in Julia distribution. nothing to do here" end'

../egadscommon/egads_types_ctts.jl: $(ESP_ROOT)/include/egadsErrors.h $(ESP_ROOT)/include/egadsTypes.h
	@$(MAKE) -C ../egadscommon
else
.PHONY: dev
dev:
	@echo "julia executable not found!"

.PHONY: test
test:
	@echo "julia executable not found!"

.PHONY: cleanall
cleanall: clean
endif

.PHONY: clean
clean:
