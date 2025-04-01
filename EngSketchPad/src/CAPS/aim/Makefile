#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
#
LDIR  = $(ESP_ROOT)/lib

export AIM_UTILSDIR = $(realpath ./utils)

# Some versions of make leave the '/' others remove it...
DIRS = $(filter-out utils/, $(filter-out utils, $(sort $(dir $(wildcard */Makefile)))))

.PHONY:		utils $(DIRS)

default:	$(DIRS)

lint:		$(DIRS)
	$(MAKE) -C utils $@

dox:		$(DIRS)

doxclean:	$(DIRS)

clean:		$(DIRS)
	$(MAKE) -C utils $@

cleanall:	$(DIRS)
	$(MAKE) -C utils $@

utils:
	$(MAKE) -C utils


ifeq ("$(MAKECMDGOALS)","")
$(DIRS): utils
else
$(DIRS):
endif
	$(MAKE) -C $@ $(MAKECMDGOALS)
