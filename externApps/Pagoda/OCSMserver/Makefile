#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib

default:	Demo1_server Demo1_client

Demo1_server:	Demo1_server.c
	mpicc -o Demo1_server -I$(IDIR) Demo1_server.c \
		-L$(LDIR) $(LIBPATH) -locsm -legads $(RPATH) -lm

Demo1_client:	Demo1_client.c
	mpicc -o Demo1_client -I$(IDIR) Demo1_client.c \
		-L$(LDIR) -legadslite $(RPATH) -lm

clean:
	-rm Demo1_server Demo1_client

test:
	mpirun -np 1 Demo1_server tutorial1_new.csm : -np 3 Demo1_client
