#	Makefile for EGADS_dot
#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

VPATH=$(ODIR)

all: circle1_dot circle2_dot naca_dot ruled_naca_dot

circle1_dot:	$(ODIR)/circle1_dot.o $(ODIR)/pingBodies.o
	$(CC) -o circle1_dot $(ODIR)/circle1_dot.o $(ODIR)/pingBodies.o \
		-L$(LDIR) -legads -lm $(RPATH)

circle2_dot:	$(ODIR)/circle2_dot.o $(ODIR)/pingBodies.o
	$(CC) -o circle2_dot $(ODIR)/circle2_dot.o $(ODIR)/pingBodies.o \
		-L$(LDIR) -legads -lm $(RPATH)

naca_dot:	$(ODIR)/naca_dot.o $(ODIR)/pingBodies.o
	$(CC) -o naca_dot $(ODIR)/naca_dot.o $(ODIR)/pingBodies.o \
		-L$(LDIR) -legads -lm $(RPATH)

ruled_naca_dot:	$(ODIR)/ruled_naca_dot.o $(ODIR)/pingBodies.o
	$(CC) -o ruled_naca_dot $(ODIR)/ruled_naca_dot.o $(ODIR)/pingBodies.o \
		-L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/%.o:	%.c 
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. $< -o $@

test: circle1_dot circle2_dot naca_dot ruled_naca_dot
	$(VALGRIND_COMMAND) ./circle1_dot
	$(VALGRIND_COMMAND) ./circle2_dot
	$(VALGRIND_COMMAND) ./naca_dot
	$(VALGRIND_COMMAND) ./ruled_naca_dot


clean:
	rm -f $(ODIR)/circle1_dot.o \
		$(ODIR)/circle2_dot.o \
		$(ODIR)/naca_dot.o \
		$(ODIR)/pingBodies.o \
		$(ODIR)/ruled_naca_dot.o

cleanall:	clean
	rm -f circle1_dot circle2_dot naca_dot ruled_naca_dot
