
ifdef PYTHONLIB
#Add rpath to PYTHONLIB
PYTHONLPATH=$(filter -L%,$(PYTHONLIB))
PYTHONRPATH=$(PYTHONLPATH:-L%=-Wl,-rpath %)
PYTHONLIBFLAGS=$(PYTHONLIB) $(PYTHONRPATH)
endif

ifdef PYTHONINC
#Extract the python include directory for Python.h dependency
PYTHONIPATH0=$(word 1,$(filter -I%,$(PYTHONINC)))
ifeq ("$(PYTHONIPATH0)","")
PYTHONIPATH=$(PYTHONINC)
else
PYTHONIPATH=$(PYTHONIPATH0:-I%=%)
endif
endif