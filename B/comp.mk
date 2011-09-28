#
# comp.mk. Compiler include file
# Generated: Wed Sep 28 00:59:45 BST 2011 
# Compiler: ifort
#
export FC=ifort
# path for the objects
export OBJSPATH=$(ROOTPATH)/obj/$(PROGNAME)/fc_$(FC)/
# modules path for specific compiler
export MODPATH=$(MODGENPATH)/fc_$(FC)/
# library path for *.a compiled with $(FC)
export LIBAPATH=$(LIBPATH)/fc_$(FC)/
FFLAGS_g= -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
FFLAGS_o= -132 -O4
SWITCH=ifort
LDFLAGS= -L.
MODFLAG=-I. -module $(MODPATH)
