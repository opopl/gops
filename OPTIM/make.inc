
ROOTPATH=$(PPATH)/../
INCPATH = $(ROOTPATH)/INCLUDE/

include $(INCPATH)/def

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE=

SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

FFLAGS+= -I $(PPATH)/NEB
FFLAGS+= -I $(PPATH)/CONNECT 
FFLAGS+= -I $(PPATH)/AMH  

LDFLAGS = -L.
