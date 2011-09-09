
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE=$(wildcard *.inc.*)  \
	$(wildcard *.i.*) \
	$(wildcard *.ref.*) \
	$(wildcard *.old.*) \
	$(wildcard *.other.*)

#SOURCE=$(filter-out,$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE)),$(GENFFILES))
SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

LDFLAGS = -L.