
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def.mk

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE := $(shell test -f nu.mk && cat nu.mk )
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

#libs {{{
BLAS_SOURCE=$(shell find $(BLAS_PATH) -name \*.f )
LAPACK_SOURCE=$(shell find $(LAPACK_PATH) -name \*.f )
AMH_SOURCE=$(shell find $(AMH_PATH) -name \*.f -o -name \*.f90 )
BLAS_OBS := 
#}}}

LDFLAGS = -L.
