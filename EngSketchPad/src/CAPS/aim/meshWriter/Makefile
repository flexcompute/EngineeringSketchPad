#
DIRS = $(filter-out NMakefile, $(filter-out Makefile, $(wildcard *)))

.PHONY:		$(DIRS)

default:	$(DIRS)

$(DIRS):
	make -C $@ $(MAKECMDGOALS)

lint:		$(DIRS)

clean:		$(DIRS)

cleanall:	$(DIRS)

