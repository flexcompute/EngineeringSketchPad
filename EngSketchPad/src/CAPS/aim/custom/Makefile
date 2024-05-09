#
DIRS = $(filter-out NMakefile, $(filter-out Makefile, $(wildcard *)))

.PHONY:		$(DIRS)

default:	$(DIRS)

$(DIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

lint:		$(DIRS)

clean:		$(DIRS)

cleanall:	$(DIRS)

dox:		$(DIRS)

doxclean:	$(DIRS)