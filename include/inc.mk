
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def.mk

#SOURCE and OBJS {{{

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
#NOTUSEDSOURCE := $(filter-out $(shell cat nu.mk ),$(wildcard #*) ) 
#NOTUSEDSOURCE := $(shell echo `cat nu.mk | sed '/^#/d'`  ) 
NOTUSEDSOURCE := $(shell cat nu.mk ) 
NOTUSEDSOURCE += $(shell find old -name \*.f -o -name \*.f90 -o -name \*.F ) 
NOTUSEDSOURCE += $(wildcard *.inc.*)  \
	$(wildcard *.i.*) \
	$(wildcard *.ref.*) \
	$(wildcard *.old.*) \
	$(wildcard *.o.*) \
	$(wildcard *.other.*) \
	$(wildcard *.save.*) 

#SOURCE=$(filter-out,$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE)),$(GENFFILES))
SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(sort $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE)))))
#OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))
#}}}

#libs {{{
BLAS_SOURCE=$(shell find $(BLAS_PATH) -name \*.f )
LAPACK_SOURCE=$(shell find $(LAPACK_PATH) -name \*.f )
AMH_SOURCE=$(shell find $(AMH_PATH) -name \*.f -o -name \*.f90 )

BLAS_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BLAS_SOURCE))))
LAPACK_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(LAPACK_SOURCE))))
AMH_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(AMH_SOURCE))))
#}}}

LDFLAGS = -L.
