
DIRS = $(sort $(dir $(wildcard */Makefile)))

.PHONY:		$(DIRS)

default:	$(DIRS)

lint:		$(DIRS)

test:		$(DIRS)

dox:		$(DIRS)

doxclean:	$(DIRS)

clean:		$(DIRS)

cleanall:	$(DIRS)

$(DIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)
