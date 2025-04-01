#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

SHELL=/bin/bash

CSMFILES=$(wildcard ffd*.csm)
VFYFILES=$(patsubst %.csm,verify_$(CASREV).?/%.vfy,$(CSMFILES))
GSENFILES=$(patsubst %.csm,verify_$(CASREV).?/%.gsen,$(CSMFILES))
TSENFILES=$(patsubst %.csm,verify_$(CASREV).?/%.tsen,$(CSMFILES))

ifndef OUTLEVEL
OUTLEVEL=0
endif

test: testESP #testSensGeom testSensTess

testESP:
	@rm -f csmlog.txt
	@printf "============================================================\n"
	@printf "Running 'serveESP -batch -verify' on all .csm files\n"
	@set -o pipefail; \
	for f in $(CSMFILES); do \
	  printf "============================================================\n"; \
	  printf "%s;\n" $$f | tee -a csmlog.txt; \
	  $(ESP_ROOT)/bin/serveESP $$f -batch -verify -outLevel $(OUTLEVEL) $(DUMPEGADS) | tee -a csmlog.txt \
	    && printf "%-35s passed (as expected)\n" $$f \
	    || printf "%-35s failed (as NOT expected)\n" $$f; \
	  printf "\n"; \
	done
	@printf "================================         Found errors             ================================\n"
	@awk '/\.csm;$$/ { csm = $$0 }; /build not completed|ERROR:: problem in/ { printf("%-35s %s\n", csm, $$0) }' csmlog.txt
	@printf "==================================================================================================\n"
	@rm -f csmlog.txt


testSensGeom:
	@rm -f csmlog.txt
	@printf "============================================================\n"
	@printf "Running 'sensCSM -geom' on all .csm files\n"
	@set -o pipefail; \
	for f in $(CSMFILES); do \
	  printf "============================================================\n"; \
	  printf "%s;\n" $$f | tee -a csmlog.txt; \
	  $(ESP_ROOT)/bin/sensCSM $$f -geom -outLevel $(OUTLEVEL) | tee -a csmlog.txt \
	    && printf "%-35s passed (as expected)\n" $$f \
	    || printf "%-35s failed (as NOT expected)\n" $$f; \
	  printf "\n"; \
	done
	@printf "================================         Found errors             ================================\n"
	@awk '/\.csm;$$/ { csm = $$0 }; /build not completed|ERROR:: problem in/ { printf("%-35s %s\n", csm, $$0) }' csmlog.txt
	@printf "==================================================================================================\n"
	@rm -f csmlog.txt


testSensTess:
	@rm -f csmlog.txt
	@printf "============================================================\n"
	@printf "Running 'sensCSM -tess' on all .csm files\n"
	@set -o pipefail; \
	for f in $(CSMFILES); do \
	  printf "============================================================\n"; \
	  printf "%s;\n" $$f | tee -a csmlog.txt; \
	  $(ESP_ROOT)/bin/sensCSM $$f -tess -outLevel $(OUTLEVEL) | tee -a csmlog.txt \
	    && printf "%-35s passed (as expected)\n" $$f \
	    || printf "%-35s failed (as NOT expected)\n" $$f; \
	  printf "\n"; \
	done
	@printf "================================         Found errors             ================================\n"
	@awk '/\.csm;$$/ { csm = $$0 }; /build not completed|ERROR:: problem in/ { printf("%-35s %s\n", csm, $$0) }' csmlog.txt
	@printf "==================================================================================================\n"
	@rm -f csmlog.txt


addVerify: $(VFYFILES) #$(GSENFILES) $(TSENFILES)

verify_$(CASREV).?/%.vfy: %.csm
	@printf "============================================================\n"
	@printf "Updating $< verification data\n"
	@$(ESP_ROOT)/bin/serveESP -batch -addVerify -outLevel 0 $<
	@printf "\n"

verify_$(CASREV).?/%.gsen: %.csm
	@printf "============================================================\n"
	@printf "Updating $< geom sens data\n"
	@$(ESP_ROOT)/bin/sensCSM $< -geom -addVerify -outLevel 0
	@printf "\n"

verify_$(CASREV).?/%.tsen: %.csm
	@printf "============================================================\n"
	@printf "Updating $< tess sens data\n"
	@$(ESP_ROOT)/bin/sensCSM $< -tess -addVerify -outLevel 0
	@printf "\n"
	
