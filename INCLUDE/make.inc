
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE=$(wildcard *.inc.*) $(wildcard *.i.*)

SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

OBJS+= dv.o rca.o

LDFLAGS = -L.
