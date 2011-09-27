#
# comp.mk. Compiler include file
# Generated: Tue Sep 27 16:44:28 BST 2011 
# Compiler: pgi
#
export FC=pgf90
# path for the objects
export OBJSPATH=$(ROOTPATH)/obj/$(PROGNAME)/fc_$(FC)/
# modules path for specific compiler
export MODPATH=$(MODGENPATH)/fc_$(FC)/
# library path for *.a compiled with $(FC)
export LIBAPATH=$(LIBPATH)/fc_$(FC)/
MODFLAG:= -I. -module $(MODPATH)
FFLAGS_o := -fast -Mipa=fast,inline -Msmartalloc
FFLAGS_g := -Mextend -O0 -Mnoframe -g -traceback
LDFLAGS= -L.
SWITCH=pgi
FFLAGS = $(FFLAGS_g) 
FFLAGS+=$(MODFLAG)
