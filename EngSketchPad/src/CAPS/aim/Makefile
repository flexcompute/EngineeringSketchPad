#
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#
#
LDIR  = $(ESP_ROOT)/lib

export AIM_UTILSDIR = $(realpath ./utils)

DIRS = $(filter-out utils, $(filter-out NMakefile, $(filter-out Makefile, $(wildcard *))))

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
