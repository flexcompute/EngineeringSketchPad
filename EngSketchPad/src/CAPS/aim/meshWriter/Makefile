#
DIRS = $(sort $(dir $(wildcard *Writer/Makefile)))

.PHONY:		$(DIRS)

default:	$(DIRS)

$(DIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

lint:		$(DIRS)

clean:		$(DIRS)

cleanall:	$(DIRS)

dox:
doxclean:

