

# definitions {{{

FC=$(FC90)

BINPATH=$(ROOTPATH)/bin
LIBPATH=$(ROOTPATH)/lib/
LIBAPATH=$(LIBPATH)/fc_$(FC)/
LAPACKPATH=$(ROOTPATH)/lapack
BLASPATH=$(ROOTPATH)/blas
SPATH=$(ROOTPATH)/scripts/
SAPATH=$(SPATH)/all/
PROG=$(BINPATH)/$(PROGNAME)

LDFLAGS = -L.
DEFS =
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P
L=

LBASE :=$(LIBAPATH)/baseG.a
# LAPACK/BLAS Linear Algebra Libs
LBLAS := $(LIBAPATH)/libmyblas.a
LLAPACK := $(LIBAPATH)/libmylapack.a
LLIBS := $(LBLAS) $(LLAPACK) $(LBLAS)
LIBS  := $(LLIBS) $(LBASE)
DEPS := deps.mk
SEARCH_PATH =  -I..

DL := debug 
#DEFS+=MPI

PRF=$(SAPATH)/porfuncs.sh
RCA=$(SAPATH)/rca.sh
DV=$(SAPATH)/dv.sh
HEADER=$(SAPATH)/header.sh
MKDEP=$(SAPATH)/mkdep.pl
AUXF = dv.f90 rca.f90

DVOPTS=fflags "$(FFLAGS)" prog $(PROGNAME) fc_exec "$(FC)" make_opts "$(MAKE_OPTS)" 
#}}}

# pgi {{{
ifeq ($(FC),pgf90)

#FFLAGS= -Mextend -O0 -Munroll -Mnoframe 
FFLAGS= -Mextend -O0 -Mnoframe 
NOOPT = -O0 -Mextend
LDFLAGS= -L.
FREEFORMAT_FLAG= -Mfree
EXTRA_FLAGS=-module
SWITCH=pgi

endif
# }}}
# nagfor {{{

ifeq ($(FC),nagfor)

FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all

ifeq ($(DL),noopt)
FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O0
endif

ifeq ($(DL),opt)
FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O3 
endif

ifeq ($(DL),debug)
 FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte
endif

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

