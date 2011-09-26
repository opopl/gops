
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
# pgi {{{
ifeq ($(FC),pgf90)

MODFLAG:= -I. -module $(MODPATH)
FFLAGS_o := -fast -Mipa=fast,inline -Msmartalloc
FFLAGS_g := -Mextend -O0 -Mnoframe -g -traceback
FFLAGS = $(FFLAGS_g) 
LDFLAGS= -L.
SWITCH=pgi

endif
# }}}
# nagfor {{{

ifeq ($(FC),nagfor)

MODFLAG= -mdir $(MODPATH)
FFLAGS_g = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte 
FFLAGS_g = -132 -g90 
FFLAGS = $(FFLAGS_g) 
LDFLAGS= -L.
SWITCH=nag

endif
# }}}
# ifort {{{

ifeq ($(FC),ifort) 

FFLAGS_g= -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
FFLAGS_o= -132 -O4
FFLAGS := $(FFLAGS_o) 
SWITCH=ifort
LDFLAGS= -L.
MODFLAG=-I. -module $(MODPATH)

endif

# }}}
# gfortran {{{

ifeq ($(FC),gfortran)

#F0=-ffree-form
MODFLAG=-I. -J$(MODPATH)
SWITCH=gfortran
LDFLAGS =-L.
FFLAGS_g= -ffixed-line-length-none -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic 
FFLAGS_o= -ffixed-line-length-none -O3 -ftree-vectorize
FFLAGS := $(FFLAGS_g) 

endif
# }}}
#
FFLAGS += $(MODFLAG) 
FFLAGS_o += $(MODFLAG) 
FFLAGS_g += $(MODFLAG) 

