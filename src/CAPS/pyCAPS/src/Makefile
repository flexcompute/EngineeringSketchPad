#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
DDIR  = $(ESP_BLOC)/dist
BDIR  = $(ESP_BLOC)/obj/build
else
ODIR  = .
TDIR  = $(ESP_ROOT)/bin
DDIR  = $(ESP_ROOT)/dist
BDIR  = ./build
endif

ifeq ($(OS),Windows_NT)
  EXT=dll
else
  EXT=so
endif

ifdef PYTHONINC

default: $(LDIR)/pyCAPS.$(EXT)

$(LDIR)/pyCAPS.$(EXT):	$(ODIR)/pyCAPS.o $(LDIR)/$(SHLIB) $(LDIR)/$(CSHLIB)
	touch $(LDIR)/pyCAPS.$(EXT)
	rm $(LDIR)/pyCAPS.$(EXT)
	$(CC) $(SOFLGS) -o $(LDIR)/pyCAPS.$(EXT) $(ODIR)/pyCAPS.o \
		-L$(LDIR) -lcaps -locsm -legads -ludunits2 -lwsserver $(PYTHONLIB) -lpthread -lm -lz $(RPATH)

$(ODIR)/pyCAPS.o:	pyCAPS.c $(IDIR)/capsTypes.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../../include -I/usr/include/udunits2\
		-I$(PYTHONINC) -w pyCAPS.c -o $(ODIR)/pyCAPS.o
else
default:
endif
	 
clean:
	-rm -f $(ODIR)/pyCAPS.o

cleanall:	clean
	-rm -f $(LDIR)/pyCAPS.$(EXT)

clean-wheel:
	(python setup.py clean --build-temp=$(BDIR) )

cleanall-wheel: clean-wheel # Remove all wheels - keep none
	(python setup.py egg_info --egg-base=$(DDIR) rotate --keep 0 --match ".whl, .egg-info" --dist-dir=$(DDIR) )
	
# Make a wheel instead of compiling directly - using default python
wheel:
	(python setup.py build --build-base=$(BDIR) egg_info --egg-base=$(BDIR) bdist_wheel "--CC=$(CC)" "--COPTS=$(COPTS) -w" "--SOFLGS=$(SOFLGS)" --dist-dir=$(DDIR) )


# Use Python 2
wheel2: 
	(python2 setup.py build --build-base=$(BDIR) egg_info --egg-base=$(BDIR) bdist_wheel "--CC=$(CC)" "--COPTS=$(COPTS) -w" "--SOFLGS=$(SOFLGS)" --dist-dir=$(DDIR) )

# Use Python 3
wheel3: 
	(python3 setup.py build --build-base=$(BDIR) egg_info --egg-base=$(BDIR) bdist_wheel "--CC=$(CC)" "--COPTS=$(COPTS) -w" "--SOFLGS=$(SOFLGS)" --dist-dir=$(DDIR) )
	
# Make both 2 and 3 version
wheels: wheel2 wheel3

# Install pyCAPS with default python
install: wheel
	(pip install pyCAPS --force-reinstall --upgrade --no-index --find-links=$(DDIR) $(INSTALL)) 

# Use Python 2
install2: wheel2
	(pip2 install pyCAPS --force-reinstall --upgrade --no-index --find-links=$(DDIR) $(INSTALL))

# Use Python 3
install3: wheel3
	(pip3 install pyCAPS --force-reinstall --upgrade --no-index --find-links=$(DDIR) $(INSTALL))

# Make both 2 and 3 version
installs: install2 install3


# Run Cython - development only	
cython:
	(cython -Wextra pyCAPS.pyx)
	-rm pyCAPS.h

# Run Cython and compile pyCAPS.c - development only
all: cython default

# Create doxygen docomentation
dox:
	(cd ../doc; doxygen pyCAPS_dox.cfg; cd latex; make; mv refman.pdf pyCAPS.pdf)

# Cleanup documentation
doxclean:
	(cd ../doc; rm -f INPUT; rm -rf html latex; rm -f pyCAPS.tag)
