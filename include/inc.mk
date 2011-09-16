
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def.mk

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE := $(shell cat nu.mk )
NOTUSEDSOURCE+=$(wildcard *.inc.*)  \
	$(wildcard *.i.*) \
	$(wildcard *.ref.*) \
	$(wildcard *.old.*) \
	$(wildcard *.o.*) \
	$(wildcard *.other.*) \
	$(wildcard *.save.*) 

#SOURCE=$(filter-out,$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE)),$(GENFFILES))
SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

LDFLAGS = -L.
