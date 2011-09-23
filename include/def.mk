
# definitions {{{

# compilers 

FC=$(FC90)
FCN=$(FC)

# directories  {{{
BINPATH=$(ROOTPATH)/bin

# general library path
LIBPATH=$(ROOTPATH)/lib/
# path for the objects
OBJSPATH=$(ROOTPATH)/obj/$(PROGNAME)/fc_$(FC)/
# general modules path
MODGENPATH=$(ROOTPATH)/mod/$(PROGNAME)/
# modules path for specific compiler
MODPATH=$(MODGENPATH)/fc_$(FC)/
# library path for *.a compiled with $(FC)
LIBAPATH=$(LIBPATH)/fc_$(FC)/

LAPACKPATH=$(ROOTPATH)/lapack
BLASPATH=$(ROOTPATH)/blas
SPATH=$(ROOTPATH)/scripts/
SAPATH=$(SPATH)/all/

# }}}

PROG=$(BINPATH)/$(PROGNAME)
LDFLAGS = -L.
DEFS=
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P
L=

F_LPU := $(INCPATH)/lpu_$(PROGNAME).mk
F_NU := $(INCPATH)/nu_$(PROGNAME).mk

# Libraries and related {{{

ARCH = ar
ARCHFLAGS = cr
ARCHDELFLAGS = d
RANLIB = ranlib

LPBASE := $(LIBAPATH)/base_$(PROGNAME).a
# LAPACK/BLAS Linear Algebra Libs
LBLAS := $(LIBAPATH)/libmyblas.a
LLAPACK := $(LIBAPATH)/libmylapack.a
LLIBS := $(LBLAS) $(LLAPACK) $(LBLAS)
LIBS  := $(LLIBS) $(LBASE)
#}}}

DEPS := deps.mk
SEARCH_PATH =  -I.. -I$(MODPATH)

#DEFS+=MPI

PRF=$(SAPATH)/porfuncs.sh
RCA=$(SAPATH)/rca.sh
DV=$(SAPATH)/dv.sh
DVOPTS=fflags "$(FFLAGS)" prog $(PROGNAME) fc_exec "$(FC)" make_opts "$(MAKE_OPTS)" 
HEADER=$(SAPATH)/header.sh
MKDEP=$(SAPATH)/mkdep.pl

AUXF = dv.f90 rca.f90

#}}}

export DL := "debug"

# pgi {{{
ifeq ($(FC),pgf90)

MODFLAG:= -module $(MODPATH)
DFFLAGS:= $(F0)

FFLAGS := -Mextend -O0 -Mnoframe -g -traceback
DFFLAGS += -Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -Mpgicoff -traceback

#not working yet {{{
ifeq ($(DL),debug)

FFLAGS += -Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -Mpgicoff -traceback

endif

ifeq ($(DL),noopt)
FFLAGS += -Mextend -O0 -Mnoframe 
endif
#}}}

NOOPT = -O0 -Mextend
LDFLAGS= -L.
FREEFORMAT_FLAG= -Mfree
SWITCH=pgi

endif
# }}}
# nagfor {{{

ifeq ($(FC),nagfor)

## {{{

#FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all

#ifeq ($(DL),noopt)
#FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O0
#endif

#ifeq ($(DL),opt)
#FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O3 
#endif

#ifeq ($(DL),debug)
 #FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte
#endif# }}}

MODFLAG= -mdir $(MODPATH)
FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte 
FFLAGS = -132 -g90 

NOOPT= -O0 -132  -kind=byte
LDFLAGS= -L.
SWITCH=nag

endif
# }}}
# ifort {{{

ifeq ($(FC),ifort) 

# FC = mpif77 
# FC = mpif90  
# DEFS+= -DMPI 
#
#FFLAGS= -132 -C -g -traceback -debug full
 #FFLAGS= -132 -O0 -g -traceback -fpe:0 -check all
 FFLAGS= -132 -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
#
#
ifeq ($(DL),debug) 
FFLAGS= -132 -C -g -traceback -debug full
endif

ifeq ($(DL),opt) 
FFLAGS= -132 -O4
endif

NOOPT= -132 -O0
SWITCH=ifort
LDFLAGS= -L.
FREEFORMAT_FLAG= -free
EXTRA_FLAGS=-I

endif

# }}}
# gfortran {{{

ifeq ($(FC),gfortran)

#F0=-ffree-form
F0=-I.

FFLAGS= -ffixed-line-length-265 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
FFLAGS= -ffixed-line-length-none -O0 

ifeq ($(DL),noopt)
FFLAGS= -ffixed-line-length-265 -O0 
endif

ifeq ($(DL),opt)
FFLAGS= -ffixed-line-length-132 -O3 -ftree-vectorize
endif

ifeq ($(DL),debug)
FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
endif

NOOPT= -O0 -ffixed-line-length-132
SWITCH=gfortran
LDFLAGS = -lblas -llapack
FREEFORMAT_FLAG= -ffree-form
EXTRA_FLAGS=-I
FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic 
#-Wall
FFLAGS+=$(F0)


endif
# }}}
#
FFLAGS += $(MODFLAG) 

