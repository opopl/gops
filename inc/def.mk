
# definitions {{{
# directories  {{{
BINPATH=$(ROOTPATH)/bin
# general library path
LIBPATH=$(ROOTPATH)/lib/
# general modules path
MODGENPATH=$(ROOTPATH)/mod/$(PROGNAME)/

LAPACKPATH=$(ROOTPATH)/lapack
BLASPATH=$(ROOTPATH)/blas
SPATH=$(ROOTPATH)/scripts/
SAPATH=$(SPATH)/all/

# }}}
# compiler include file {{{
COMP := comp.mk
include $(COMP)
#}}}
#PROG LDFLAGS SEARCH_PATH CPP L  {{{
PROG=$(BINPATH)/$(PROGNAME)
LDFLAGS = -L.
SEARCH_PATH =  -I.. -I$(MODPATH)
DEFS=
CPP = /lib/cpp
CPFLAGS = -traditional -P
L=
#}}}
# files {{{
F_LPU := $(INCPATH)/lpu_$(PROGNAME).mk
F_NU := $(INCPATH)/nu_$(PROGNAME).mk
DEPS := deps.mk
INITFILES := $(DEPS) $(F_LPU) $(F_NU) t.mk inc.mk $(COMP)
# }}}
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
#LIBS  := $(LLIBS) $(LPBASE)
LIBS  := $(LLIBS)
#}}}
# scripts {{{
PRF=$(SAPATH)/porfuncs.sh
RCA=$(SAPATH)/rca.sh
DV=$(SAPATH)/dv.sh
DVOPTS=fflags "$(FFLAGS)" prog $(PROGNAME) fc_exec "$(FC)" make_opts "$(MAKE_OPTS)" 
HEADER=$(SAPATH)/header.sh
MKDEP=$(SAPATH)/mkdep.pl
MKCOMP=$(SAPATH)/mkcomp.sh
#}}}

#AUXF = dv.f90 rca.f90
AUXF = dv.f90

#}}}

