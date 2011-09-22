
ROOTPATH=$(PPATH)/../

include $(INCPATH)/def.mk

#SOURCE and OBJS {{{

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE := $(shell cat $(F_NU))
NOTUSEDSOURCE += $(shell find old -name \*.f -o -name \*.f90 -o -name \*.F ) 
NOTUSEDSOURCE += $(wildcard *.inc.*)  \
	$(wildcard *.i.*) \
	$(wildcard *.ref.*) \
	$(wildcard *.old.*) \
	$(wildcard *.o.*) \
	$(wildcard *.other.*) \
	$(wildcard *.save.*) 

SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(sort $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE)))))
#}}}

#libs {{{

LPUSED := $(strip $(shell cat $(F_LPU)))
LPBASE_SOURCE := $(filter-out $(LPUSED),$(SOURCE))
LPBASE_OBJS := $(sort $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(LPBASE_SOURCE)))))
LPU_OBJS := $(sort $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(LPUSED)))))

BLAS_SOURCE=$(shell find $(BLAS_PATH) -name \*.f )
LAPACK_SOURCE=$(shell find $(LAPACK_PATH) -name \*.f )
AMH_SOURCE=$(shell find $(AMH_PATH) -name \*.f -o -name \*.f90 )

BLAS_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BLAS_SOURCE))))
LAPACK_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(LAPACK_SOURCE))))
AMH_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(AMH_SOURCE))))
#}}}

LDFLAGS = -L.
