#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

default:	simple_client simple_server

simple_client:	simple_client.o
	mpicc -o simple_client simple_client.o -L$(LDIR) -legadslite $(RPATH) -lm

simple_client.o:	simple_client.c 
	mpicc -c $(COPTS) $(DEFINE) -I$(IDIR) simple_client.c

simple_server:	simple_server.o
	mpicc -o simple_server simple_server.o $(LIBPATH) \
		-L$(LDIR) -legads $(RPATH) -lm

simple_server.o:	simple_server.c 
	mpicc -c $(COPTS) $(DEFINE) -I$(IDIR) simple_server.c

run1:	simple_client simple_server
	mpirun -np 1 simple_server stream.egads : -np 1 simple_client

run:	simple_client simple_server
	mpirun -np 1 simple_server stream.egads : -np 2 simple_client

clean:
	-rm simple_client.o simple_server.o

cleanall:	clean
	-rm simple_client simple_server
