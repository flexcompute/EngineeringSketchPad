
ifneq (,$(findstring DARWIN,$(ESP_ARCH)))
# Python version of scan-build
#
# Install with: pip install scan-build
#
SCANBUILD=intercept-build --override-compiler --cdb=$(SCANDIR)compile_commands.json $(MAKE) -f $(word 1,$(MAKEFILE_LIST)) CC=intercept-cc CXX=intercept-c++ && analyze-build --cdb=$(SCANDIR)compile_commands.json -o $(SCANDIR) $(SCANEXCLUDE)

.PHONY: scan-build
scan-build:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	$(SCANBUILD) --status-bugs

.PHONY: scan-view
scan-view:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	$(SCANBUILD)
	-@open $(SCANDIR)/*/index.html || true
endif

ifeq ("$(ESP_ARCH)","LINUX64")
# Perl version of scan-build
#
# Install on ubuntu: sudo apt-get install clang-tools
#
.PHONY: scan-build
scan-build:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	scan-build -o $(SCANDIR) $(SCANEXCLUDE) --status-bugs $(MAKE) -f $(word 1,$(MAKEFILE_LIST))

.PHONY: scan-view
scan-view:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	scan-build -o $(SCANDIR) $(SCANEXCLUDE) $(MAKE) -f $(word 1,$(MAKEFILE_LIST))
	-@xdg-open $(SCANDIR)/*/index.html || true
endif
