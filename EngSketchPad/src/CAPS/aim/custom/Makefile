#
DIRS = $(filter-out NMakefile, $(filter-out Makefile, $(wildcard *)))

.PHONY:		$(DIRS)

default:	$(DIRS)

$(DIRS):
	-(cd $@; make $(MAKECMDGOALS))

clean:		$(DIRS)

cleanall:	$(DIRS)

