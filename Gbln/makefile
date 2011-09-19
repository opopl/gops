#vim:ft=make

PROGNAME=$(shell basename $(PWD)  )
GENFFILES = porfuncs.f90 dv.f90 
#
#
PPATH=$(PWD)
INCPATH = $(PPATH)/../include/

#include file
include $(INCPATH)/inc.mk
include inc.mk
#
# targets
include $(INCPATH)/t.mk
include t.mk

$(PROG): $(LIBS)

$(DEPS): rca.f90 dv.f90 porfuncs.f90 nu.mk

include $(DEPS)
